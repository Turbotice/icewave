"""
chipiron_files_debug.py -- scan the whole Chipiron / Phonefleet expedition tree and tell you WHAT TO CHECK.

It reuses the loader of CASCADES26_expedition_v2.py (same test descriptions, same references), then runs checks on
every file of every test folder, every day and every .edc, and writes:
    audit_report_<stamp>.txt      this console report (grouped by problem, with what to do)
    audit_issues_<stamp>.csv      one row per problem found (severity, check, date, time, detail, action)
    audit_tests_<stamp>.csv       one row per test (usb duration, edc, rest reading, plateau ratio, n errors...)
    audit_edc_<stamp>.csv         one row per .edc used (AO functions, scalings, speed range, ramps, #tests)

RUN (from the folder holding CASCADES26_expedition_v2.py and the other chipiron_*.py files)
    py chipiron_audit.py
    py chipiron_audit.py --root "D:/Phonefleet_data" --days 2026-11-13 2026-11-14
    py chipiron_audit.py --no-decode          # file-level checks only (fast: no usb decoding)
    py chipiron_audit.py --only-errors        # hide INFO / WARN in the console (they stay in the csv)

SEVERITY
    ERROR  data unusable or certainly wrong (empty usb, edc not found, medium unknown, usb much shorter than planned)
    WARN   probably wrong / needs a decision (no edc_used.txt, reads > 10 % above the command, no temperature...)
    INFO   worth knowing (3-10 % calibration offset, near-match reference, no matching reference...)

WHAT IS CHECKED AUTOMATICALLY (see CHECKS below for the exact list) AND WHAT CAN ONLY BE CHECKED BY YOU
    The last block of the report prints the MANUAL CHECKLIST (ESCON Studio really had that .edc? phone files complete?...).
"""
import argparse
import os
import re
import sys
import types
from datetime import datetime
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")          # no windows: this script only scans

import numpy as np
import pandas as pd

EXPEDITION_FILE = "CASCADES26_expedition.py"     # looked for next to this script
# ---- thresholds (edit here)
REST_TOL_RPM = 25.0              # |median speed| before the first command above this = offset at rest
RATIO_INFO = 0.03                # plateau / peak reading differs from the command by more than this -> INFO
RATIO_WARN = 0.10                # ... more than this -> WARN (wrong .edc scaling or ADC_VREF)
SHORT_USB = 0.90                 # usb duration < this x planned duration -> ERROR
COVERAGE_MIN = 0.85              # decoded telemetry cycles / expected cycles below this -> WARN
GAP_MIN_S = 2.0                  # time gap between decoded cycles above max(this, 20 periods) -> WARN
SAT_FRAC = 0.01                  # fraction of ADC readings >= 1020 counts above this -> WARN
BAD_CYCLE_FRAC = 0.02            # fraction of cycles failing the byte consistency check above this -> WARN
CLOCK_OFFSET_MAX_S = 120.0       # |input-log /start - phone folder time| above this -> WARN
SPIKE_FRAC = 0.01                # fraction of plateau samples > 10 % away from the plateau median -> INFO
AO_TOP_V = 4.0                   # ESCON analog outputs span 0..4 V once shifted (chipiron_patch_edc)
RANGE_MARGIN = 0.02              # commanded value within this fraction of the AO range edge = clipping risk


# Exceptions.xlsx sits next to this script by default. Each column is a date; cells contain
# "HH:MM:SS -- instruction". Instructions are deliberately narrow: they suppress only the
# checks explicitly named, or exclude a test/day when the note says "skip it/day".
EXCEPTIONS_FILE = Path(__file__).resolve().with_name("Exceptions.xlsx")


def _exception_date(value):
    if isinstance(value, (datetime, pd.Timestamp)):
        return pd.Timestamp(value).strftime("%Y-%m-%d")
    try:
        return pd.to_datetime(str(value), dayfirst=True).strftime("%Y-%m-%d")
    except Exception:
        return None


def load_exceptions(path=EXCEPTIONS_FILE):
    """Read the Exceptions workbook into day-wide and test-specific directives."""
    path = Path(path)
    rules = {"skip_days": set(), "tests": {}}
    if not path.is_file():
        return rules
    try:
        raw = pd.read_excel(path, header=None, dtype=object)
    except Exception as exc:
        print(f"[audit] WARNING: cannot read {path.name}: {exc}")
        return rules

    for col in range(raw.shape[1]):
        day = _exception_date(raw.iat[0, col]) if raw.shape[0] else None
        if not day:
            continue
        for row in range(1, raw.shape[0]):
            value = raw.iat[row, col]
            if pd.isna(value) or not str(value).strip():
                continue
            note = str(value).strip()
            tm = re.search(r"(?<!\d)(\d{1,2}:\d{2}:\d{2})(?!\d)", note)
            if not tm:
                if "skip this day" in note.lower():
                    rules["skip_days"].add(day)
                continue
            hhmmss = tm.group(1)
            # The workbook has occasionally had a month/day vs day/month transposition.
            # Keep the exact date as primary; a reverse-date alias is handled at lookup time.
            rules["tests"].setdefault((day, hhmmss), []).append(note)
    return rules


def exception_notes(rules, date, time):
    """Exact date/time match, then a unique transposed-date match for legacy Excel headers."""
    t = str(time or "").replace("-", ":")
    hits = rules["tests"].get((str(date), t), [])
    if hits:
        return hits
    try:
        d = pd.Timestamp(date)
        swapped = f"{d.year:04d}-{d.day:02d}-{d.month:02d}"
        return rules["tests"].get((swapped, t), [])
    except Exception:
        return []


def exception_policy(notes):
    """Return (skip_test, skip_checks, skip_warnings) from human-readable notes."""
    skip_test, skip_checks, skip_warnings = False, set(), False
    for note in notes:
        low = note.lower()
        if "skip it" in low and "don't skip it" not in low and "do not skip it" not in low and "skip warnings" not in low:
            skip_test = True
        if "skip warnings" in low:
            skip_warnings = True
        if "test is correct" in low and ("don't skip it or warn" in low or "do not skip it or warn" in low):
            skip_checks.add("cmd_failed")
        if "acquisition is known to be incomplete" in low:
            skip_checks.add("usb_short")
        if "skip the multi-files" in low or "skip multi-files" in low:
            skip_checks.add("multi_files")
        if "usb-short" in low or "usb short" in low:
            skip_checks.add("usb_short")
        m = re.search(r"skip(?: the)? ([a-z0-9_, -]+?) warning", low)
        if m:
            names = m.group(1).replace(" and ", ",")
            for token in re.split(r"[, ]+", names.strip()):
                if token:
                    skip_checks.add(token.replace("-", "_"))
    return skip_test, skip_checks, skip_warnings


def usb_files_for_test(e):
    """All USB CSV parts belonging to a test folder (split recordings are expected)."""
    if not e.folder.is_dir():
        return []
    return sorted((p for p in e.folder.iterdir()
                   if p.is_file() and p.suffix.lower() == ".csv" and "usb" in p.name.lower()),
                  key=lambda p: p.name.lower())


def total_usb_duration(files):
    """Sum the recorded duration of each USB part; gaps between split files are not counted."""
    total = 0.0
    found = False
    for file in files:
        try:
            s = pd.to_numeric(pd.read_csv(file, header=None, usecols=[0]).iloc[:, 0], errors="coerce")
            s = s[np.isfinite(s) & (s > 0)]
            if len(s) >= 2:
                total += float((s.max() - s.min()) / 1e6)
                found = True
        except Exception:
            continue
    return total if found else None


def failed_motor_commands(e):
    """Ignore failed terminal x/rest commands; they do not mean a motor test command failed."""
    log = getattr(e, "inputs", None)
    if log is None or log.empty:
        return None
    if "ok" not in log.columns or "command" not in log.columns:
        return None
    failed = log[pd.to_numeric(log["ok"], errors="coerce").fillna(1).eq(0)]
    if failed.empty:
        return False
    for _, r in failed.iterrows():
        command = str(r.get("command", "")).strip().lower()
        note = str(r.get("note", "")).strip().lower()
        if command == "x" or note.startswith("rest"):
            continue
        if command[:1] in ("c", "o", "r"):
            return True
    return False

# ============================================================================ the checks (id -> severity, meaning, what to do)
CHECKS = {
    # ---- test folder: files
    "empty_folder":   ("ERROR", "test folder is empty", "Delete the folder or pull the files again from the phone (mobile_v3: ls / pull)."),
    "usb_missing":    ("ERROR", "no Usb csv in the folder", "Look for it on the phone (mobile_v3 `ls`, `pull`). Without it there is no speed / current."),
    "usb_empty":      ("ERROR", "Usb csv has no (or < 4) valid rows: only placeholder rows", "Pull the file again. If it is still empty the telemetry was never logged (USB cable / Arduino reset / phone app): the data is lost, use a repeat of this test."),
    "usb_undecodable": ("ERROR", "Usb csv has rows but no full telemetry cycle can be decoded", "Check the firmware version / sampling_period_ms.txt: a wrong period or a different byte layout breaks the decoder."),
    "usb_short":      ("ERROR", "Usb data stops long before the planned end of the test", "Is the file truncated on the phone (buffer full?) or did you stop the acquisition early? Compare with `ls` on the phone and with the input log; if it is the phone buffer, long tests need a longer buffer or several recordings."),
    "usb_gaps":       ("WARN", "holes in the decoded telemetry", "Lost bytes / Wi-Fi or USB dropouts. Look at the plot around that time before using the run."),
    "usb_coverage":   ("WARN", "far fewer decoded cycles than expected from the sampling period", "Is sampling_period_ms.txt the period of the firmware that was flashed? Else many cycles were lost."),
    "adc_saturated":  ("WARN", "ADC readings saturate (>= 1020 counts)", "The ESCON analog output exceeds what the Arduino reads: check AO scaling in the .edc / ESCON Studio."),
    "adc_flat":       ("WARN", "A1 never changes although the motor was commanded", "AO1 not connected / wrong wiring / motor disabled."),
    "bad_cycles":     ("WARN", "many cycles fail the byte consistency check", "Corrupted bytes: the affected samples are NaN. Check the USB cable / baud rate."),
    "no_sampling_period": ("WARN", "no sampling_period_ms.txt", "The decoder assumes the default period. Write the period of the flashed firmware in this file (250 or 10)."),
    "no_input_log":   ("ERROR", "no input_*.csv: the commands sent are unknown", "Test type, onset time and reference matching are impossible. If you remember the test, write its input log by hand; else exclude it."),
    "multi_files":    ("WARN", "several input csv in the same folder", "Check that only one input log belongs to this test."),
    "no_sensors":     ("INFO", "accelerometer / gyroscope / magnetic_field file missing", "Only matters if you use that sensor."),
    # ---- edc
    "no_edc_txt":     ("WARN", "no edc_used.txt: speed cannot be converted to rpm (plots fall back to VOLTS)", "Find which .edc was loaded in ESCON Studio during this test (lab notes, neighbouring tests of the same series, ESCON Studio logs) and write its file name on the first line of edc_used.txt."),
    "edc_not_found":  ("ERROR", "edc_used.txt names an .edc that cannot be found", "Copy that .edc next to the data or add its folder to EDC_SEARCH_FOLDERS in the expedition script."),
    "edc_wiring":     ("ERROR", "AO1 / AO2 functions of the .edc are not (speed, current): rpm columns would be wrong", "Check which analog output feeds A1 and A2 on this day; the script assumes A1 <- AO1 speed, A2 <- AO2 current."),
    "edc_range":      ("WARN", "commanded speeds go beyond the analog-output range of the .edc", "The measured speed is clipped at the edge of the AO range (e.g. a negative speed on a 0..4 V unipolar scaling). Use an .edc with a wider range for such tests."),
    "rest_offset":    ("WARN", "speed reads far from 0 before the first command", "Either the motor really turns at the rest command (set value offset, see 'offset adjust' in the .edc) or the conversion (edc / ADC_VREF) is wrong. Compare with ESCON Studio at rest."),
    "ratio_warn":     ("WARN", "measured plateau / peak is > 10 % away from the commanded speed", "Almost certainly the wrong .edc scaling for this test (ratio of two AO scalings) or a wrong ADC_VREF. Check which .edc was REALLY loaded in ESCON Studio."),
    "ratio_info":     ("INFO", "measured plateau / peak is 3-10 % away from the commanded speed", "Systematic calibration offset (ADC_VREF, ESCON AO accuracy). Same offset on the reference tests = fine when normalising by them."),
    "sine_attenuated": ("INFO", "sine amplitude much smaller than commanded", "Low-pass of the speed loop / ESCON ramp limit at this frequency, load, or a wrong scale. Compare with the reference."),
    "plateau_spikes": ("INFO", "many samples far from the plateau value", "Noise bursts or stick-slip in the medium; or telemetry glitches. Look at the raw plot."),
    "no_physical":    ("INFO", "speed could not be computed even with an .edc", "See the detail (usually the decoder or the AO functions)."),
    # ---- test description / logs
    "test_unknown":   ("WARN", "test type unknown or no motor command found", "Open the input log: was it a rest recording, or a command format the loader does not know?"),
    "cmd_failed":     ("WARN", "some commands were not delivered (ok = 0)", "The motor may not have done what the log says: compare with the Usb speed."),
    "clock_offset":   ("WARN", "laptop clock and phone folder time differ a lot", "Temperature / photo matching uses the laptop clock: check the clocks (or fix TEMP_CLOCK_OFFSET_S)."),
    "loader_issue":   ("WARN", "other issue reported by the loader", "See the detail."),
    # ---- medium / station / temperature
    "medium_unknown": ("ERROR", "medium unknown", "Fill the metadata txt (Milieu) and/or LOG_days_media_Chipiron.xlsx, or use MEDIUM_OVERRIDES."),
    "station_missing": ("WARN", "no station for a frazil / saltwater test", "Check LOG_UMA_Chipiron.xlsx, the metadata txt and the Notes time ranges, or use STATION_OVERRIDES."),
    "no_temperature": ("WARN", "no temperature file matched", "Is there an Exelog csv for that hour? Check the Exelog clock vs the laptop (TEMP_CLOCK_OFFSET_S) and the EXELOG folder of the day."),
    "temp_overlap":   ("INFO", "temperature matched by overlap, not by start time", "Fine for a long Exelog file; check the clock if unexpected."),
    "temp_no_sample": ("WARN", "temperature file matched but no sample inside the test window", "Exelog clock shifted, or the logger was off during the test."),
    # ---- references
    "ref_unusable":   ("ERROR", "REFERENCE test whose Usb cannot be used", "Repeat that reference test, or exclude it (REFERENCE_EXCLUDE) and rely on another reference type."),
    "ref_all_unusable": ("WARN", "all the references with this input are unusable", "Same as above: this test has no usable reference of that type."),
    "no_reference":   ("INFO", "no reference test with the same input", "Compare against the nearest-input reference by hand, or repeat the test on the reference day."),
    "ref_near":       ("INFO", "reference is a near match (inputs differ slightly)", "Check the ~NEAR fields before normalising."),
    # ---- day level
    "no_metadata":    ("WARN", "no metadata txt for the day", "Write metadataDDMMYYYY.txt (medium, blade, temperatures, heights)."),
    "no_exelog":      ("WARN", "no Exelog temperature files for the day", "Copy the Exelog csv into <day>/EXELOG*/."),
    "temp_unmatched": ("INFO", "Exelog files not matched to any test", "Normal for files recorded outside the tests; else check the clocks."),
    "no_images":      ("INFO", "grain-size / sampling-height images missing", "Copy the pictures of the day into grain_size/ and sampling_height/."),
    "no_image_csv":   ("INFO", "images not processed (no vertical_and_grain_size csv)", "Run the image processing for that station."),
    "crosscheck":     ("WARN", "metadata txt / media xlsx / UMA xlsx disagree", "Fix the typo in whichever file is wrong (the report line tells which)."),
}

MANUAL_CHECKLIST = [
    ("phone files (Usb / accelerometer / gyroscope / magnetic_field csv)",
     ["Is the file complete? Compare the file list on the phone (mobile_v3 `ls`) with what is saved; large files can be cut during `pull`.",
      "Usb duration far shorter than the test (see usb_short): phone buffer full vs acquisition stopped early (only you know what you did).",
      "Phone clock vs laptop clock (tsync): the report only compares the folder name with the input log."]),
    ("input_*.csv (commands sent)",
     ["Was the run what the log says (e.g. Ctrl-C during a program)? The log only records what was SENT.",
      "Commands with ok = 0 were not delivered; check the motor really moved in the Usb speed."]),
    ("edc_used.txt and the .edc files",
     ["The audit can read the .edc but cannot know what ESCON Studio REALLY had loaded: open ESCON Studio notes / photos for each series.",
      "Compare the per-day .edc timeline of the report with your lab notebook: the AO range is changed from series to series.",
      "Open each .edc in ESCON Studio once: AO scaling, ramp (accel / decel), mode of operation, 'offset adjust' enabled or not.",
      "ADC_VREF (5.06 in the expedition script, 4.85 in the plotting script): measure the 5 V pin of the Arduino on the phone supply and put ONE value everywhere."]),
    ("metadata<DDMMYYYY>.txt",
     ["Station(s), medium, blade, temperatures, heights: typos are not detectable beyond the xlsx cross-check.",
      "Multiple media / stations in a day need the Notes time ranges in LOG_days_media_Chipiron.xlsx."]),
    ("LOG_days_media_Chipiron.xlsx / LOG_UMA_Chipiron.xlsx",
     ["Every test day present, one row per day, Notes time ranges in the 'from 10am to 2:40pm' format.",
      "Station <-> day assignments (the report lists every disagreement as a cross-check)."]),
    ("Exelog temperature csv",
     ["Logger clock vs laptop clock (set TEMP_CLOCK_OFFSET_S if shifted); channel -> probe mapping (ch2..ch5); 'High' = open input."]),
    ("images (grain_size/, sampling_height/) and vertical_and_grain_size_day*_station*.csv",
     ["Right station and day in the file names; values plausible; pictures actually taken for every station sampled."]),
    ("reference days (13/11 empty, 14/11 ninja, 16/11 salt water)",
     ["Every reference test really carries the medium of its day (the audit warns if the logs say otherwise).",
      "Repeat the reference tests that are unusable (see ref_unusable)."]),
]


# ============================================================================ load the expedition module without running its cells
def load_expedition_module(path):
    path = Path(path)
    if not path.is_file():
        raise SystemExit(f"cannot find {path}: put chipiron_audit.py next to {EXPEDITION_FILE} or use --expedition FILE")
    sys.path.insert(0, str(path.resolve().parent))
    src = path.read_text(encoding="utf-8-sig")
    head = src.split("# %% LOAD the expedition")[0]      # definitions only: the LOAD / PLOT / DEBUG cells are not run
    mod = types.ModuleType("cascades_defs")
    mod.__file__ = str(path)
    sys.modules["cascades_defs"] = mod
    exec(compile(head, str(path), "exec"), mod.__dict__)
    return mod


# ============================================================================ issue collector
class Audit:
    def __init__(self, exceptions=None):
        self.issues = []
        self.exceptions = exceptions or load_exceptions()
        self.current_skip_checks = set()
        self.current_skip_warnings = False

    def add(self, check, e=None, detail="", date=None, time=None):
        sev, what, action = CHECKS[check]
        if check in self.current_skip_checks or (self.current_skip_warnings and sev == "WARN"):
            return
        self.issues.append(dict(severity=sev, check=check, date=date or (e.date if e else ""), time=time or (e.time if e else ""),
                                medium=getattr(e, "medium", ""), test=(getattr(e, "test_key", None) or getattr(e, "test_type", "")) if e else "",
                                detail=detail, meaning=what, action=action))

    def count(self, sev):
        return sum(1 for i in self.issues if i["severity"] == sev)


# ============================================================================ .edc inspection (once per file)
def inspect_edc(m, path, cache):
    key = str(path)
    if key in cache:
        return cache[key]
    info = dict(path=key, name=Path(key).name, ok=False, error="", ao1="", ao2="", ao1_slope=np.nan, ao1_offset=np.nan,
                speed_lo=np.nan, speed_hi=np.nan, ao2_slope=np.nan, accel=np.nan, decel=np.nan, op_mode="", polarity="",
                offset_adjust="", wiring_ok=False)
    try:
        from escon_edc_v3 import load_edc
        edc = load_edc(path)
        a1, a2 = edc.ao["AO1"], edc.ao["AO2"]
        info.update(ok=True, ao1=f"{a1.function} ({a1.quantity})", ao2=f"{a2.function} ({a2.quantity})",
                    ao1_slope=a1.slope_V_per_unit, ao1_offset=a1.offset_V, ao2_slope=a2.slope_V_per_unit,
                    op_mode=edc.enum_label("PARAM_ModeOfOperation") or "", polarity=edc.pwm_polarity() or "",
                    accel=edc.try_value("PARAM_SpeedRampAcceleration") or np.nan,
                    decel=edc.try_value("PARAM_SpeedRampDeceleration") or np.nan)
        if edc.offset_adjust:
            info["offset_adjust"] = f"{edc.offset_adjust[0]:g} rpm ({'ENABLED' if edc.offset_adjust[1] else 'not enabled'})"
        info["wiring_ok"] = a1.quantity == "speed" and a2.quantity == "current"
        if a1.usable:
            lo, hi = sorted(float(x) for x in a1.volts_to_phys([0.0, AO_TOP_V]))
            info.update(speed_lo=lo, speed_hi=hi)
    except Exception as exc:
        info["error"] = f"{type(exc).__name__}: {exc}"
    cache[key] = info
    return info


# ============================================================================ calibration checks (pure function: testable)
def commanded_range(e):
    """(lo, hi) of the commanded speeds [rpm] from the test description, or None."""
    t = e.test_type
    f = lambda v: v if v is not None else np.nan
    if t == "sine" and e.mean_rpm is not None and e.amp_rpm is not None:
        return e.mean_rpm - e.amp_rpm, e.mean_rpm + e.amp_rpm
    if t in ("sweep", "sweep_stepped") and e.mean_rpm is not None and e.amp_start_rpm is not None:
        a = max(abs(e.amp_start_rpm), abs(e.amp_end_rpm or 0))
        return e.mean_rpm - a, e.mean_rpm + a
    if t == "hold" and e.hold_rpm is not None:
        return float(e.hold_rpm), float(e.hold_rpm)
    if t in ("ramp", "steps") and e.speed_min_rpm is not None:
        return min(0.0, float(e.speed_min_rpm)), float(e.speed_max_rpm)
    return None


def calibration_checks(e, t_rel, sp):
    """Pure function on arrays (t_rel = s since the first motor command, sp = speed [rpm]).
    Returns (list of (check, detail), dict of measures)."""
    out, meas = [], {}
    t_rel, sp = np.asarray(t_rel, float), np.asarray(sp, float)
    ok = np.isfinite(t_rel) & np.isfinite(sp)
    t_rel, sp = t_rel[ok], sp[ok]
    if len(sp) < 5:
        return out, meas
    # ---- rest reading before the first command
    pre = (t_rel < -1.0)
    if pre.sum() >= 3:
        rest = float(np.median(sp[pre]))
        meas["rest_rpm"] = rest
        if abs(rest) > REST_TOL_RPM:
            out.append(("rest_offset", f"median {rest:+.0f} rpm over {int(pre.sum())} samples before the first command (limit {REST_TOL_RPM:g})"))
    # ---- plateau / peak against the command
    t, cr = e.test_type, commanded_range(e)
    if t in ("hold", "steps") and e.levels:
        val, dur = e.levels[0]
        dur = dur if dur is not None and dur > 0 else 60.0
        w = (t_rel >= min(15.0, 0.5 * dur)) & (t_rel <= dur)
        if w.sum() >= 5 and abs(val) > 1:
            plateau = float(np.median(sp[w]))
            ratio = plateau / val
            meas["plateau_ratio"] = ratio
            meas["plateau_rpm"] = plateau
            _ratio_issue(out, ratio, f"plateau reads {plateau:.0f} rpm for a command of {val:.0f} rpm (x{ratio:.3f})")
            dev = np.abs(sp[w] - plateau) > 0.10 * max(abs(plateau), 50.0)
            if dev.mean() > SPIKE_FRAC:
                out.append(("plateau_spikes", f"{100 * dev.mean():.1f} % of the plateau samples are > 10 % from {plateau:.0f} rpm"))
    elif t == "ramp" and cr is not None and abs(cr[1]) > 1:
        peak = float(np.percentile(sp[t_rel >= 0], 99)) if (t_rel >= 0).sum() > 5 else np.nan
        if np.isfinite(peak):
            ratio = peak / cr[1]
            meas["plateau_ratio"] = ratio
            meas["plateau_rpm"] = peak
            _ratio_issue(out, ratio, f"top of the ramp reads {peak:.0f} rpm for a command of {cr[1]:.0f} rpm (x{ratio:.3f})")
    elif t in ("sine", "sweep", "sweep_stepped") and cr is not None and (cr[1] - cr[0]) > 2:
        amp_cmd = 0.5 * (cr[1] - cr[0])
        w = t_rel >= 0
        if w.sum() >= 10:
            hi, lo = float(np.percentile(sp[w], 99)), float(np.percentile(sp[w], 1))
            amp_meas = 0.5 * (hi - lo)
            ratio = amp_meas / amp_cmd
            meas["plateau_ratio"] = ratio
            meas["plateau_rpm"] = amp_meas
            if ratio > 1 + RATIO_INFO:                                  # a speed loop cannot exceed its command: over-reading = scale
                _ratio_issue(out, ratio, f"measured amplitude {amp_meas:.0f} rpm for a command of {amp_cmd:.0f} rpm (x{ratio:.3f})")
            elif ratio < 0.5 and t == "sine":
                out.append(("sine_attenuated", f"amplitude {amp_meas:.0f} rpm = {100 * ratio:.0f} % of the command ({amp_cmd:.0f} rpm)"))
    return out, meas


def _ratio_issue(out, ratio, msg):
    dev = abs(ratio - 1.0)
    if dev > RATIO_WARN:
        out.append(("ratio_warn", msg))
    elif dev > RATIO_INFO:
        out.append(("ratio_info", msg))


# ============================================================================ checks of one test
def audit_test(m, e, A, edc_cache, decode=True):
    notes = exception_notes(A.exceptions, e.date, e.time)
    skip_test, skip_checks, skip_warnings = exception_policy(notes)
    if e.date in A.exceptions.get("skip_days", set()):
        skip_test = True
    if skip_test:
        return dict(date=e.date, time=e.time, medium=e.medium, station=e.station,
                    test_type=e.test_type, test_key=e.test_key, skipped_by_exception=True,
                    n_issues=0, n_errors=0), 0

    previous_skip, previous_warning = A.current_skip_checks, A.current_skip_warnings
    A.current_skip_checks = set(skip_checks)
    A.current_skip_warnings = bool(skip_warnings)
    row = dict(date=e.date, time=e.time, medium=e.medium, station=e.station, test_type=e.test_type, test_key=e.test_key,
               is_reference=e.is_reference, ref_type=e.ref_type, edc=e.edc_name, usb_duration_s=e.usb_duration_s,
               planned_duration_s=e.duration_s, usb_ok=None, n_cycles=None, rest_rpm=np.nan, plateau_ratio=np.nan,
               plateau_rpm=np.nan)
    n0 = len(A.issues)

    def add(check, detail=""):
        A.add(check, e, detail)

    files = [f for f in e.folder.iterdir() if f.is_file()] if e.folder.is_dir() else []
    if not files:
        add("empty_folder")
        return row, 0
    # ---- files
    csvs = [f for f in files if f.suffix.lower() == ".csv"]
    # Multiple USB CSVs are normal when a long acquisition is split into parts.
    if len(list(e.folder.glob("input_*.csv"))) > 1:
        add("multi_files", "several input_*.csv")
    if e.usb_file is None:
        add("usb_missing")
    if e.input_log_file is None:
        add("no_input_log")
    if not (e.folder / "sampling_period_ms.txt").is_file():
        add("no_sampling_period")
    miss = [k for k in ("accelerometer", "gyroscope", "magnetic_field") if k not in e.sensor_files]
    if miss:
        add("no_sensors", ", ".join(miss))
    # ---- edc
    edc_info = None
    if not (e.folder / "edc_used.txt").is_file():
        add("no_edc_txt")
    elif e.cfg is None or e.edc_path is None:
        add("edc_not_found", f"edc_used.txt says '{e.edc_name}'")
    else:
        edc_info = inspect_edc(m, e.edc_path, edc_cache)
        edc_info["n_tests"] = edc_info.get("n_tests", 0) + 1
        if not edc_info["ok"]:
            add("edc_not_found", f"{e.edc_name} cannot be read: {edc_info['error']}")
        elif not edc_info["wiring_ok"]:
            add("edc_wiring", f"{e.edc_name}: AO1 = {edc_info['ao1']}, AO2 = {edc_info['ao2']}")
        else:
            cr = commanded_range(e)
            lo, hi = edc_info["speed_lo"], edc_info["speed_hi"]
            if cr is not None and np.isfinite(lo) and np.isfinite(hi):
                span, tol = hi - lo, RANGE_MARGIN * (hi - lo)
                if cr[1] > hi - tol or cr[0] < lo + tol:
                    # a command of exactly 0 at the edge of a unipolar range is fine
                    if not (abs(cr[0] - lo) < 1e-6 and cr[0] >= 0):
                        add("edc_range", f"commanded {cr[0]:.0f}..{cr[1]:.0f} rpm, AO range of {e.edc_name} = {lo:.0f}..{hi:.0f} rpm")
    # ---- test description / log
    if e.test_type in ("unknown", "rest"):
        add("test_unknown", f"type '{e.test_type}'")
    for iss in e.issues:
        if iss.startswith("cross-check:") or iss in ("no edc_used.txt", "no usb file", "no input log", "no temperature match") \
                or iss.startswith("usb not decodable") or iss.startswith("medium unknown") or "no motor command" in iss \
                or iss.startswith("no input log: test type unknown"):
            continue
        if iss.startswith("input log unreadable"):
            add("no_input_log", iss)
        elif "not delivered" in iss:
            failed = failed_motor_commands(e)
            if failed is True or failed is None:
                add("cmd_failed", iss)
        else:
            add("loader_issue", iss)
    if e.clock_offset_s is not None and abs(e.clock_offset_s) > CLOCK_OFFSET_MAX_S:
        add("clock_offset", f"input-log /start is {e.clock_offset_s:+.0f} s from the phone folder time")
    # ---- medium / station / temperature
    if e.medium == "unknown":
        add("medium_unknown")
    if e.station is None and e.medium in m.STATION_MEDIA:
        add("station_missing")
    if e.temp_df is None:
        add("no_temperature")
    else:
        if e.temp_match == "overlap":
            add("temp_overlap")
        if not e.temp_mean_c:
            add("temp_no_sample")
    # ---- USB duration across every split part (the expedition loader keeps the first path)
    usb_parts = usb_files_for_test(e)
    combined_usb_duration = total_usb_duration(usb_parts)
    if combined_usb_duration is not None:
        row["usb_duration_s"] = combined_usb_duration
    # ---- usb content
    if e.usb_file is not None and decode:
        ok, why = m.usb_usable(e)
        row["usb_ok"] = ok
        if not ok:
            add("usb_empty" if "valid row" in why else "usb_undecodable", why)
        else:
            _usb_checks(m, e, A, add, row, edc_info, combined_usb_duration)
    elif e.usb_file is not None:
        try:
            n = len(e.usb_raw())
            row["usb_ok"] = n >= 4
            if n < 4:
                add("usb_empty", f"only {n} valid row(s) in the usb csv")
        except Exception as exc:
            add("usb_undecodable", f"{type(exc).__name__}: {exc}")
    row["n_issues"] = len(A.issues) - n0
    row["n_errors"] = sum(1 for i in A.issues[n0:] if i["severity"] == "ERROR")
    A.current_skip_checks, A.current_skip_warnings = previous_skip, previous_warning
    return row, len(A.issues) - n0


def _usb_checks(m, e, A, add, row, edc_info, combined_usb_duration=None):
    d = e.usb(physical=False)
    full = d[d["kind"] == "full"] if "kind" in d.columns else d
    row["n_cycles"] = len(full)
    period = e.sampling_period_ms or 250.0
    onset = m.test_onset_s(e) or 0.0
    # duration vs plan
    duration = combined_usb_duration if combined_usb_duration is not None else e.usb_duration_s
    if e.duration_s and duration is not None and duration < SHORT_USB * (onset + e.duration_s):
        add("usb_short", f"Usb lasts {duration:.0f} s, the test needs about {onset + e.duration_s:.0f} s "
                         f"(first command at +{onset:.0f} s, planned {e.duration_s:.0f} s)")
    # coverage, gaps, saturation, checksum
    if e.usb_duration_s:
        expected = e.usb_duration_s * 1000.0 / period
        if expected > 20 and len(full) / expected < COVERAGE_MIN:
            add("usb_coverage", f"{len(full)} full cycles decoded, ~{expected:.0f} expected at {period:g} ms")
    t = full["t_s"].to_numpy(float)
    if len(t) > 3:
        dt = np.diff(t)
        thr = max(GAP_MIN_S, 20 * period / 1000.0)
        big = dt > thr
        if big.any():
            k = int(np.argmax(dt))
            add("usb_gaps", f"{int(big.sum())} gap(s) > {thr:.1f} s, largest {dt[k]:.1f} s at t = {t[k]:.0f} s")
    for col in ("A1_adc", "A2_adc"):
        if col in full.columns:
            v = full[col].to_numpy(float)
            v = v[np.isfinite(v)]
            if len(v) and (v >= 1020).mean() > SAT_FRAC:
                add("adc_saturated", f"{col}: {100 * (v >= 1020).mean():.1f} % of the readings >= 1020")
    if "A1_adc" in full.columns and len(full):
        a1 = full["A1_adc"].to_numpy(float)
        if np.isnan(a1).mean() > BAD_CYCLE_FRAC:
            add("bad_cycles", f"{100 * np.isnan(a1).mean():.1f} % of the cycles fail the consistency check")
        fin = a1[np.isfinite(a1)]
        if e.test_type not in ("rest", "unknown") and len(fin) > 10 and np.ptp(fin) < 3:
            add("adc_flat", f"A1 varies by {np.ptp(fin):.0f} counts only")
    # calibration against the command (needs the .edc)
    if edc_info is not None and edc_info.get("ok") and edc_info.get("wiring_ok") and e.test_type not in ("rest", "unknown"):
        try:
            dp = e.usb(physical=True)
            col = next(c for c in dp.columns if c.startswith("A1_") and c.endswith("_rpm"))
            dp = dp[dp["kind"] == "full"] if "kind" in dp.columns else dp
            dp = dp.dropna(subset=[col])
            found, meas = calibration_checks(e, dp["t_s"].to_numpy(float) - onset, dp[col].to_numpy(float))
            for check, detail in found:
                add(check, detail)
            row.update({k: v for k, v in meas.items() if k in row})
        except Exception as exc:
            add("no_physical", f"{type(exc).__name__}: {exc}")


# ============================================================================ day / global / reference checks
def audit_days(m, exps, A):
    for date, di in (exps.days or {}).items():
        if "metadata" in di["missing"] or "metadata (unreadable)" in di["missing"]:
            A.add("no_metadata", date=date, detail=f"{di['n_tests']} test(s) that day")
        if di["n_temp_files"] == 0 and di["n_tests"]:
            A.add("no_exelog", date=date)
        if di["temp_unmatched"]:
            A.add("temp_unmatched", date=date, detail=", ".join(di["temp_unmatched"][:10]))
        if any("images" in x for x in di["missing"]):
            A.add("no_images", date=date, detail=", ".join(x for x in di["missing"] if "images" in x))
        if any("processed image csv" in x for x in di["missing"]):
            A.add("no_image_csv", date=date, detail=", ".join(x for x in di["missing"] if "processed" in x))
    block = False
    for ln in (exps.summary_text or "").splitlines():             # the cross-check block of the loader summary
        if ln.startswith("CROSS-CHECKS"):
            block = True
            continue
        if block and ln.startswith("  ?? "):
            mm = re.match(r"  \?\? (\d{4}-\d{2}-\d{2}): (.*)", ln)
            A.add("crosscheck", date=mm[1] if mm else "", detail=mm[2] if mm else ln[5:])
        elif block and ln.strip() and not ln.startswith("  "):
            block = False


def audit_references(m, exps, A, decode):
    pool = exps.ref_pool if exps.ref_pool is not None else exps
    if decode:
        for r in pool:
            if r.is_reference and not m.usb_usable(r)[0]:
                A.add("ref_unusable", r, f"{r.ref_type}: {m.usb_usable(r)[1]}")
    for e in exps:
        if e.is_reference or e.test_type in ("unknown", "rest"):
            continue
        if not any(e.ref_available.values()):
            A.add("no_reference", e, "no empty / ninja / saltwater reference with this input")
            continue
        if decode:
            for rt, lst in e.ref_available.items():
                if lst and not any(m.usb_usable(r)[0] for r in lst):
                    A.add("ref_all_unusable", e, f"{rt}: {len(lst)} candidate(s), none decodable")
                elif lst:
                    nf = m.near_fields(e, lst[0])
                    if nf:
                        A.add("ref_near", e, f"{rt}: {'; '.join(nf)}")


def edc_timeline(exps):
    """{date: [(first time, last time, edc name or None, n tests)]} : consecutive tests with the same .edc are merged."""
    out = {}
    for date in sorted({e.date for e in exps}):
        seq = []
        for e in sorted((x for x in exps if x.date == date), key=lambda x: x.time):
            name = e.edc_name if (e.folder / "edc_used.txt").is_file() else None
            if seq and seq[-1][2] == name:
                seq[-1][1], seq[-1][3] = e.time, seq[-1][3] + 1
            else:
                seq.append([e.time, e.time, name, 1])
        out[date] = seq
    return out


# ============================================================================ report
def build_report(A, exps, edc_cache, timeline, rows, only_errors=False):
    L = []
    P = L.append
    n_e, n_w, n_i = A.count("ERROR"), A.count("WARN"), A.count("INFO")
    bad_tests = {(r["date"], r["time"]) for r in rows if r.get("n_errors")}
    P("=" * 100)
    P(f"CHIPIRON DATA AUDIT   {datetime.now():%Y-%m-%d %H:%M:%S}")
    P(f"{len(exps)} test(s) in {len({e.date for e in exps})} day(s): {n_e} ERROR, {n_w} WARN, {n_i} INFO "
      f"({len(bad_tests)} test(s) with at least one ERROR)")
    P("=" * 100)
    order = {"ERROR": 0, "WARN": 1, "INFO": 2}
    by_check = {}
    for i in A.issues:
        by_check.setdefault(i["check"], []).append(i)
    for check in sorted(by_check, key=lambda c: (order[CHECKS[c][0]], -len(by_check[c]))):
        sev, what, action = CHECKS[check]
        if only_errors and sev != "ERROR":
            continue
        items = by_check[check]
        P("")
        P(f"[{sev}] {check}: {what}   ({len(items)})")
        P(f"    -> what to do: {action}")
        days = {}
        for i in items:
            days.setdefault(i["date"], []).append(i)
        show_detail = len(items) <= 15
        for date in sorted(days):
            its = days[date]
            if show_detail:
                for i in its:
                    P(f"    {date} {i['time']:<9} {i['medium']:<15} {i['detail']}")
            else:
                P(f"    {date}: " + ", ".join(i["time"] for i in its if i["time"]) + (f"  [{its[0]['detail']}]" if not its[0]["time"] else ""))
    P("")
    P("=" * 100)
    P(".EDC FILES USED (check against your notes; the audit cannot know what ESCON Studio really had)")
    P("=" * 100)
    for info in sorted(edc_cache.values(), key=lambda x: x["name"]):
        if not info["ok"]:
            P(f"  {info['name']:<40} UNREADABLE: {info['error']}")
            continue
        P(f"  {info['name']:<40} {info.get('n_tests', 0):>4} test(s) | AO1 {info['ao1']}, AO2 {info['ao2']}")
        P(f"      speed range {info['speed_lo']:.0f}..{info['speed_hi']:.0f} rpm (AO 0..{AO_TOP_V:g} V; {info['ao1_slope']:.4g} V/rpm, "
          f"offset {info['ao1_offset']:g} V) | {info['polarity']} | {info['op_mode']} | ramp {info['accel']:g}/{info['decel']:g} rpm/s"
          + (f" | offset adjust {info['offset_adjust']}" if info["offset_adjust"] else ""))
    P("")
    P(".EDC TIMELINE PER DAY (compare with your lab notebook: a change of .edc between series should match what you did)")
    for date, seq in timeline.items():
        P(f"  {date}")
        for a, b, name, n in seq:
            P(f"      {a}" + (f"..{b}" if b != a else "") + f"  ({n} test(s))  " + (name or "NO edc_used.txt"))
    P("")
    P("=" * 100)
    P("MANUAL CHECKLIST (cannot be checked automatically)")
    P("=" * 100)
    for what, items in MANUAL_CHECKLIST:
        P(f"  {what}")
        for it in items:
            P(f"      [ ] {it}")
    P("")
    P("Files written next to this report: audit_issues_*.csv (all problems), audit_tests_*.csv (one row per test), audit_edc_*.csv")
    return "\n".join(L)


# ============================================================================ main
def main(argv=None):
    ap = argparse.ArgumentParser(description="Scan the Chipiron data tree and list what to check")
    ap.add_argument("--root", default=None, help="Phonefleet_data folder (default: ROOT_FOLDER of the expedition script)")
    ap.add_argument("--days", nargs="*", default=None, help="only these days, e.g. 2026-11-13 2026-11-14")
    ap.add_argument("--expedition", default=None, help=f"path of {EXPEDITION_FILE}")
    ap.add_argument("--out", default=None, help="folder for the report and csv (default: <root>/audit)")
    ap.add_argument("--no-decode", action="store_true", help="file-level checks only (no usb decoding)")
    ap.add_argument("--only-errors", action="store_true", help="console: show ERROR groups only")
    args = ap.parse_args(argv)

    here = Path(__file__).resolve().parent
    m = load_expedition_module(args.expedition or here / EXPEDITION_FILE)
    root = Path(args.root) if args.root else Path(m.ROOT_FOLDER)
    decode = not args.no_decode
    print(f"[audit] loading the expedition from {root} ...")
    exps = m.load_expedition(root, args.days, None, verbose=False)
    print(f"[audit] {len(exps)} test(s); running the checks{' (usb decoding can take a minute)' if decode else ''} ...")

    exceptions = load_exceptions()
    A, edc_cache, rows = Audit(exceptions), {}, []
    print(f"[audit] Exceptions: {len(exceptions['skip_days'])} skipped day(s), "
          f"{len(exceptions['tests'])} test-specific note(s) from {EXCEPTIONS_FILE.name}")
    for i, e in enumerate(exps, 1):
        row, _ = audit_test(m, e, A, edc_cache, decode)
        rows.append(row)
        if i % 50 == 0:
            print(f"[audit]   {i}/{len(exps)}")
    audit_days(m, exps, A)
    audit_references(m, exps, A, decode)
    timeline = edc_timeline(exps)

    report = build_report(A, exps, edc_cache, timeline, rows, args.only_errors)
    print("\n" + report)

    out = Path(args.out) if args.out else root / "audit"
    out.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    (out / f"audit_report_{stamp}.txt").write_text(build_report(A, exps, edc_cache, timeline, rows, False), encoding="utf-8")
    pd.DataFrame(A.issues).to_csv(out / f"audit_issues_{stamp}.csv", index=False)
    pd.DataFrame(rows).to_csv(out / f"audit_tests_{stamp}.csv", index=False)
    pd.DataFrame(list(edc_cache.values())).to_csv(out / f"audit_edc_{stamp}.csv", index=False)
    print(f"\n[audit] written in {out.resolve()}")
    return 1 if A.count("ERROR") else 0


if __name__ == "__main__":
    sys.exit(main())