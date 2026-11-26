"""
chipiron_expedition_v1.py -- STEP 1: load every Chipiron / Phonefleet test of the expedition into Experiment objects.

ONE EXPERIMENT = ONE <date>/<HH-MM-SS> folder = one Experiment object holding
    * identity            date, time, station
    * day metadata        metadataDDMMYYYY.txt of the day (medium, blades, temperatures, heights...)
    * the test itself     read from the input_*.csv log of the folder (commands sent to the Arduino):
                          test_type, direction, freq_hz, amp_rpm, mean_rpm, duration_s, levels...
    * raw data (LAZY)     usb (raw / decoded / physical units), accelerometer, gyroscope, magnetic_field, gps
    * temperature         the Exelog file matched in time (uses match_temperature_usb.py: read_temperature, match)
    * day logs (xlsx)     LOG_days_media_Chipiron.xlsx (day -> medium) and LOG_UMA_Chipiron.xlsx (station -> days processed)
                          are CROSS-CHECKED with the metadata txt: every difference raises a "CROSS-CHECK" warning
                          (typo?) and is listed in a dedicated block of the summary and in exp.issues.
                          A day may hold SEVERAL media and/or stations: the Notes column of the media log gives the
                          time ranges ("frazil 101: from 10am to 11:30am. Empty; from 1:30pm to 2:50pm. ...") and every
                          test gets the medium and station of the range containing its start time.
    * images              grain-size and sampling-height pictures of the day (any image extension) and, if already
                          processed, the values of vertical_and_grain_size_day<xxxxxx>_station<xxx>.csv

REFERENCE DAYS (hardcoded in REFERENCE_DAYS below)
    13/11 = empty motor (empty_ref), 14/11 = labgrown slush (ninja_ref), 16/11 = salt water (saltwater_ref).
    EVERY test of those days is marked is_reference=True, ref_type=<empty_ref|ninja_ref|saltwater_ref>, and gets a
    test_key built from its INPUT LOG (all test types): e.g. "sine mean_rpm=0, amp_rpm=50, freq_hz=1".
    Every other test gets ref_available = {ref_type: [matching reference tests]} (same test type + same parameters,
    within a tolerance; the number of cycles / duration is ignored for sine and hold).
    Plot modes: no_ref | empty_ref | ninja_ref | saltwater_ref  ->  exps.pairs("empty_ref") gives [(test, ref or None)].

EXPECTED TREE (anything missing / empty is skipped with a warning, never a crash)
    ROOT/YYYY-MM-DD/metadata*.txt            (else ROOT/metadata<DDMMYYYY>.txt)
    ROOT/YYYY-MM-DD/HH-MM-SS/                usb / sensor csv, input_*.csv, edc_used.txt, sampling_period_ms.txt
    ROOT/YYYY-MM-DD/EXELOG*/DA-HH-MM.csv    temperature files of the day
    ROOT/YYYY-MM-DD/grain_size/              images   (+ vertical_and_grain_size_day*_station*.csv anywhere in the day)
    ROOT/YYYY-MM-DD/sampling_height/         images

USE
    from chipiron_expedition_v1 import load_expedition, between
    exps = load_expedition()                                  # prints + returns the summary, nothing heavy is read
    exps = load_expedition(days=["2026-11-15"], summary_folder="D:/summaries")
    df   = exps.table()                                       # one row per test (pandas), for quick looks
    sel  = exps.filter(medium="frazil", test_type="sine", freq_hz=between(0.1, 1))
    sel  = exps.filter(station=["400", "401"], date="2026-11-15", direction="down")
    e    = sel[0];  d = e.usb(physical=True);  T = e.temperature();  e.inputs  # commands sent
    for e in sel: ...                                         # plots are then easy to build
    # references (normalisation of the plots)
    exps.references("empty_ref")                              # all tests of the empty-motor day
    exps.find_reference(e, "ninja_ref")                       # reference test with the same input as e (or None)
    exps.find_references(e, "saltwater_ref")                  # list (several repeats of the same reference test)
    for e, r in exps.filter(medium="frazil", test_type="sine").pairs("empty_ref"):   # r is None if no match
        ...                                                   # r.usb(physical=True) etc. as for any Experiment
    exps.filter(is_reference=True, ref_type="ninja_ref")

SPYDER USE (everything is hardcoded in the cells at the bottom of the file)
    1. run the file once (F5) or the cells 'Helpers' and 'LOAD': loads E, saves the console in console_log_*.txt,
       prints the reference report (which test found which reference, per medium)
    2. run any '# %% PLOT n' cell (filters: reference type, medium, station, test type, ...); copy a cell for a new plot

USABLE USB DATA (v2)
    usb_usable(e) -> (ok, reason): a test is usable if its usb csv decodes into at least one full telemetry cycle.
    * references are chosen USABLE FIRST (PREFER_USABLE_REFS): among the reference tests with the same input, an
      undecodable one is only used if no usable one exists (then ref_available / find_reference / pairs follow that order);
    * decodability_report(E) lists every test whose usb cannot be decoded (+ which references that affects);
    * debug_match_plots closes each figure after saving it (no 40 Qt windows, no blocking plt.show()).

DEPENDENCIES  numpy, pandas, match_temperature_usb.py (temperature), chipiron_helpers_v4.py + escon_edc_v3.py (optional:
              rpm from the .edc), chipiron_decode_v4.py or _v3.py (optional, only when you call usb()), Pillow (optional,
              image time stamps).
"""
import argparse
import contextlib
import inspect
import io
import re
import sys
import unicodedata
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

# ======================= USER CONFIGURATION =======================
ROOT_FOLDER = Path(r"C:\Users\Gulliver\Desktop\Chipiron\git_clone_turbotice_icewave\icewave\icewave\chipiron\Phonefleet_data")
SUMMARY_FOLDER = None            # e.g. Path("D:/CASCADES26/summaries"): the summary .txt is saved there (None = not saved)
EXTRA_TEMP_FOLDERS = []          # more folders with Exelog csv (searched for every day, matched by real time stamps)
EDC_SEARCH_FOLDERS = []          # extra folders where the .edc named in edc_used.txt may be found
START_THRESHOLD_S = 300          # temperature file accepted as "same experiment" if its start is this close [s]
TEMP_CLOCK_OFFSET_S = 0.0        # added to the Exelog times if the logger clock is shifted [s]
TEMP_WINDOW_SOURCE = "log"       # time window used to match temperature: "log" (laptop clock, /start../stop of the input
                                 # log) or "usb" (phone folder name + usb duration)
DAYS_MEDIA_LOG = Path(r"C:\Users\Gulliver\Desktop\Chipiron\git_clone_turbotice_icewave\icewave\icewave\chipiron\Phonefleet_data\LOG_days_media_Chipiron.xlsx")            # None = ROOT/LOG_days_media_Chipiron.xlsx  (.xlsx or .csv: Day | Medium)
UMA_LOG = Path(r"C:\Users\Gulliver\Desktop\Chipiron\git_clone_turbotice_icewave\icewave\icewave\chipiron\Sorties_UMA\LOG_UMA_Chipiron.xlsx")                  # None = ROOT/../Sorties_UMA/LOG_UMA_Chipiron.xlsx  (Station | Sampled by | Day processed by Chipiron)
CROSSCHECK_PRIORITY = "txt"      # medium: which source wins when metadata txt and media xlsx disagree ("txt" or "xlsx");
                                 # (the Notes time ranges, when present, always decide the medium of each test;
                                 #  stations come from the UMA log, else from the txt)
NOTES_TOL_MIN = 2                # a test starting up to this many minutes outside a Notes range is still assigned to it
STATION_MEDIA = ("frazil", "saltwater")   # media for which a station number makes sense (not air-empty / labgrown)
STATION_OVERRIDES = {}           # {("2026-11-04", "10-30-00"): "101"} or {"2026-11-04": "101"}; wins over logs
ADC_VREF = 5.06 #4.85                  # Arduino ADC full scale [V] (as in Chipiron_plotting_v3)
WIRING = {"A1": "AO1", "A2": "AO2"}
IMAGE_EXTS = {".jpg", ".jpeg", ".png", ".tif", ".tiff", ".bmp", ".heic", ".webp"}
GRAIN_COL_RE = r"grain|diam|size"                 # column of the processed csv taken as grain size (first match)
HEIGHT_COL_RE = r"height|vertical|hauteur|thick"  # column taken as sampling height (first match)
REST_TOL_DUTY = 0.25             # |duty - rest duty| below this = motor at rest [%]
# ---- reference days: EVERY test of these days is a reference test, matched to the other tests by its input
REFERENCE_DAYS = {"2026-11-13": "empty_ref",       # empty motor (blade in air)
                  "2026-11-14": "ninja_ref",       # labgrown slush (ninja)
                  "2026-11-16": "saltwater_ref"}   # salt water
REFERENCE_EXCLUDE = set()        # {("2026-11-14", "12-26-04")}: a test of a reference day that must NOT be used as reference
REF_EXPECTED_MEDIUM = {"empty_ref": "air-empty", "ninja_ref": "labgrown_frazil", "saltwater_ref": "saltwater"}
                                 # a warning is raised if a reference-day test has another medium (typo in the logs?)
REF_MODES = ("no_ref", "empty_ref", "ninja_ref", "saltwater_ref")   # the normalisation modes of the plots
REF_RTOL = 0.002                 # two tests have "the same input" if every parameter agrees within this (tight) relative tolerance...
REF_ATOL = {"freq_hz": 1e-3, "duration_s": 1.0, "cycles": 0.5, "plateau_s": 1.0, "n_levels": 0.0}   # ...or this absolute one
REF_ATOL_DEFAULT = 0.5           # [rpm] floor for speeds / amplitudes
REF_SLOPE_FALLBACK = 66.0        # [rpm per % duty] used when a test has no .edc (bipolar 2640 rpm file)
REF_DUTY_STEP_C = 0.1            # [% duty] one step of "c<duty>" (speeds of hold / ramp / steps): tolerance = 1 step
REF_DUTY_STEP_O = 0.01           # [% duty] one step of "o<mean>,<amp>,<freq>" (amplitudes of sine / sweep)
REF_NEAR_IGNORE = ("plateau_s", "duration_s", "cycles")   # fields not reported as a "near match" difference
PREFER_USABLE_REFS = True        # among references with the same input, a reference whose usb decodes comes BEFORE one that does not
                                 # (costs one usb decode per reference test at load time; False = old behaviour, closest input only)
MEDIUM_OVERRIDES = {}            # {("2026-11-15", "13-04-16"): "air-empty"} or {"2026-11-15": "frazil"}; wins over metadata
# ==================================================================

MEDIA = ("air-empty", "saltwater", "frazil", "labgrown_frazil")
TIME_RE = re.compile(r"^\d{2}-\d{2}-\d{2}$")
DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")
N = r"(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"
SETUP_NOTES = ("arduino mode", "0 rpm before start", "rest after experiment")
SENSOR_PATTERNS = {"accelerometer": "accelero", "gyroscope": "gyro", "magnetic_field": "magneti", "gps": "gps"}
SENSOR_COLS = {"accelerometer": ["t_us", "ax", "ay", "az"], "gyroscope": ["t_us", "gx", "gy", "gz"],
               "magnetic_field": ["t_us", "mx", "my", "mz"], "gps": ["t_us", "gps_1", "gps_2", "gps_3"]}
SINE_RE = re.compile(rf"sine mean {N} rpm, amp {N} rpm, {N} Hz, {N} cycles")
SWEEP_RE = re.compile(rf"continuous amp sweep {N}->{N} rpm \((lin|log)\), mean {N}, {N} Hz, {N} s")
STEPPED_RE = re.compile(rf"amp sweep (\d+)/(\d+): mean {N}, amp {N} rpm, {N} Hz, {N} cycles")
BLADE_RE = re.compile(rf"'([^']+)'\s*radius\s*{N}\s*cm\s*,?\s*height\s*{N}\s*cm", re.I)
IMG_CSV_RE = re.compile(r"vertical_and_grain_size_day(?P<day>\d+)_station(?P<st>[A-Za-z0-9]+)\.csv$", re.I)


# ============================================================================ small helpers
def _norm(s):
    """lower case, no accents."""
    return "".join(c for c in unicodedata.normalize("NFD", str(s)) if unicodedata.category(c) != "Mn").lower()


def _nums(s):
    return [float(x.replace(",", ".")) for x in re.findall(r"-?\d+(?:[.,]\d+)?", str(s))]


def _num_or_range(s):
    """'28.5cm' -> 28.5 ; '14-20cm' -> 17.0 (mean of the range) ; None if no number."""
    m = re.match(r"^\s*(\d+(?:[.,]\d+)?)\s*-\s*(\d+(?:[.,]\d+)?)", str(s))
    if m:
        return 0.5 * (float(m[1].replace(",", ".")) + float(m[2].replace(",", ".")))
    v = _nums(s)
    return v[0] if v else None


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def load_edc_quiet(path):
    """escon_edc_v3.load_edc without its 'offset adjust ... stored but NOT enabled / ENABLED' warning (not wanted here).
    Any other line printed by the loader is kept."""
    from escon_edc_v3 import load_edc
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        edc = load_edc(path)
    for ln in buf.getvalue().splitlines():
        if "offset adjust" not in ln.lower():
            print(ln)
    return edc


def usb_part_key(p):
    """Order the parts of a recording split by the phone: ...Usb-1-..., ...Usb-2-... (same rule as Chipiron_plotting_v5)."""
    m = re.search(r"usb-(\d+)", p.name, re.I)
    return (int(m.group(1)) if m else 10 ** 9, p.name)


class Report:
    """Collects warnings (printed as they come) and per-day information for the summary."""

    def __init__(self, verbose=True):
        self.verbose = verbose
        self.warnings = []
        self.crosschecks = []
        self.days = {}

    def warn(self, where, msg):
        self.warnings.append(f"{where}: {msg}")
        if self.verbose:
            print(f"  ! {where}: {msg}")


# ============================================================================ xlsx logs (days/media, station/days)
DATE_FORMATS = ("%d-%m-%Y", "%d/%m/%Y", "%d-%m-%y", "%d/%m/%y", "%Y-%m-%d", "%d.%m.%Y")


def _read_table(path):
    p = Path(path)
    if p.suffix.lower() in (".xlsx", ".xlsm", ".xls"):
        return pd.read_excel(p, header=None, dtype=object)
    return pd.read_csv(p, header=None, sep=None, engine="python", dtype=str, encoding="utf-8-sig")


def _station_str(v):
    if v is None or (isinstance(v, float) and np.isnan(v)) or str(v).strip().lower() in ("", "nan"):
        return None
    try:
        f = float(v)
        return str(int(f)) if f == int(f) else str(v).strip()
    except (TypeError, ValueError):
        return str(v).strip()


def _stations_in(value):
    """Parse one or several station numbers from metadata."""
    if value is None:
        return []
    return list(dict.fromkeys(
        s for s in (
            _station_str(x)
            for x in re.split(r"[,;/+&]", str(value))
        )
        if s is not None
    ))


def _dates_in(cell, rep, where):
    """Cell with one or several dates (datetime, or text 'dd-mm-yyyy, dd-mm-yyyy') -> list of 'YYYY-MM-DD'."""
    if cell is None or (isinstance(cell, float) and np.isnan(cell)):
        return []
    if isinstance(cell, (datetime, pd.Timestamp)):
        return [pd.Timestamp(cell).strftime("%Y-%m-%d")]
    out = []
    for tok in re.split(r"[,;]", str(cell)):
        tok = tok.strip()
        if not tok:
            continue
        for fmt in DATE_FORMATS:
            try:
                out.append(datetime.strptime(tok.split(" ")[0], fmt).strftime("%Y-%m-%d"))
                break
            except ValueError:
                continue
        else:
            rep.warn(where, f"CROSS-CHECK: unreadable day '{tok}' (typo?)")
    return out


def _split_media(s):
    """'empty, labgrown frazil' -> ['air-empty', 'labgrown_frazil'] (unknown words are kept as typed, to be flagged)."""
    parts = re.split(r"[,;/+&]|\bet\b|\band\b", str(s), flags=re.I)
    return [norm_medium(x) for x in parts if x.strip()]


TIME_TOKEN = r"\d{1,2}(?:[:h]\d{2})?\s*(?:[ap]\.?m\.?)?"
SEG_RE = re.compile(rf"([A-Za-z][A-Za-z_\- ]*?)\s*(?:(\d{{2,4}})\s*)?[:;]\s*from\s+({TIME_TOKEN})\s+(?:to|until)\s+({TIME_TOKEN})", re.I)


def _minutes(tok):
    """'10am' / '2:40pm' / '14:30' / '14h30' -> minutes since midnight (None if unreadable)."""
    m = re.match(r"(\d{1,2})(?:[:h](\d{2}))?\s*([ap])?", str(tok).strip().lower().replace(".", ""))
    if not m:
        return None
    h, mi, ap = int(m[1]), int(m[2] or 0), m[3]
    if ap == "p" and h < 12:
        h += 12
    if ap == "a" and h == 12:
        h = 0
    return h * 60 + mi


def _hhmm(m):
    return f"{int(m) // 60:02d}:{int(m) % 60:02d}"


def parse_notes(text, flag):
    """'frazil 101: from 10am to 11:30am. Empty; from 1:30pm to 2:50pm.' -> [dict(medium, station, t0, t1, text)]"""
    segs = []
    for m in SEG_RE.finditer(str(text)):
        label, st, a, b = m[1].strip(), m[2], _minutes(m[3]), _minutes(m[4])
        toks = _split_media(label)
        if len(toks) != 1:
            flag(f"notes: cannot read the medium in '{m[0].strip()}' (typo?)")
            continue
        if toks[0] not in MEDIA:
            flag(f"notes: unrecognised medium '{toks[0]}' in '{m[0].strip()}' (typo?)")
        if a is None or b is None or b <= a:
            flag(f"notes: time range not valid in '{m[0].strip()}' (typo?)")
            continue
        segs.append(dict(medium=toks[0], station=st, t0=a, t1=b, text=m[0].strip()))
    return sorted(segs, key=lambda x: x["t0"])


def read_day_logs(root, rep, media_path=None, uma_path=None):
    """-> dict(media={date: [raw texts]}, uma={date: {stations}}, media_file, uma_file). Missing file = warning, empty dict."""
    root = Path(root)
    logs = dict(media={}, notes={}, uma={}, media_file=None, uma_file=None)
    mp = Path(media_path or DAYS_MEDIA_LOG or root / "LOG_days_media_Chipiron.xlsx")
    up = Path(uma_path or UMA_LOG or root.parent / "Sorties_UMA" / "LOG_UMA_Chipiron.xlsx")
    for key, path in (("media", mp), ("uma", up)):
        if not path.is_file():
            rep.warn("logs", f"{path.name} not found at {path} (cross-check with it skipped)")
            continue
        try:
            df = _read_table(path)
        except Exception as exc:
            rep.warn("logs", f"{path.name} unreadable ({exc}): cross-check with it skipped")
            continue
        logs[f"{key}_file"] = path
        for i, r in df.iterrows():
            c = list(r.values) + [None, None, None]
            if i == 0 and str(c[0]).strip().lower() in ("day", "station", "date", "jour"):
                continue
            where = f"{path.name} row {i + 1}"
            if key == "media":
                days = _dates_in(c[0], rep, where)
                txt = "" if c[1] is None or (isinstance(c[1], float) and np.isnan(c[1])) else str(c[1]).strip()
                for d in days:
                    if d in logs["media"] and txt != logs["media"][d][0]:
                        rep.warn(where, f"CROSS-CHECK: day {d} appears twice with different media "
                                        f"('{logs['media'][d][0]}' / '{txt}')")
                    logs["media"].setdefault(d, []).append(txt)
                    note = "" if c[2] is None or (isinstance(c[2], float) and np.isnan(c[2])) else str(c[2]).strip()
                    if note:
                        logs["notes"].setdefault(d, []).append(note)
                if not txt:
                    rep.warn(where, "CROSS-CHECK: empty medium")
            else:
                st = _station_str(c[0])
                days = _dates_in(c[2], rep, where)
                if st is None:
                    rep.warn(where, "CROSS-CHECK: row without station")
                if not days:
                    rep.warn(where, f"CROSS-CHECK: station {st} has no (readable) day")
                for d in days:
                    logs["uma"].setdefault(d, set()).add(st)
    return logs


# ============================================================================ day metadata
META_FIELDS = [
    ("station", r"^station", "str"),
    ("medium_raw", r"^(milieu|medium)", "str"),
    ("blade", r"^pales?", "str"),
    ("lab_temp_c", r"temperature.*labo", "num"),
    ("frazil_temps_c", r"temperature.*frazil", "ch"),
    ("tank_height_cm", r"hauteur totale", "num"),
    ("blade_emergence_cm", r"emergence", "num"),
    ("water_height_start_cm", r"hauteur d'eau.*debut", "num"),
    ("water_height_end_cm", r"hauteur d'eau.*fin", "num"),
    ("frazil_thickness_end_cm", r"epaisseur frazil", "num"),
]


def norm_medium(s):
    t = _norm(s or "")
    if not t.strip():
        return "unknown"
    if re.search(r"\bair\b|vide|empty", t):
        return "air-empty"
    if "frazil" in t and re.search(r"lab|grown|cultiv|artific|made|fait", t):
        return "labgrown_frazil"
    if "frazil" in t:
        return "frazil"
    if re.search(r"sal|sea|eau|water|mer\b", t):
        return "saltwater"
    return t.strip()


def parse_metadata(path):
    """-> (meta dict with canonical fields + 'raw' {original key: text}, blade dimensions {name: {...}})."""
    text = Path(path).read_text(encoding="utf-8-sig", errors="replace")
    raw, blades = {}, {}
    for ln in text.splitlines():
        ln = ln.strip()
        m = BLADE_RE.search(ln)
        if m:
            blades[_norm(m[1]).strip()] = dict(radius_cm=float(m[2]), height_cm=float(m[3]))
            continue
        if not ln or ln.startswith("#") or ":" not in ln:
            continue
        k, v = ln.split(":", 1)
        if v.strip():                                          # empty fields are skipped
            raw[k.strip()] = v.strip()
    meta = {"raw": raw}
    for fld, rx, kind in META_FIELDS:
        for k, v in raw.items():
            if re.search(rx, _norm(k)):
                if kind == "str":
                    meta[fld] = v
                elif kind == "num":
                    meta[fld] = _num_or_range(v)
                else:
                    meta[fld] = {f"CH{c}": float(t) for c, t in re.findall(r"CH(\d)\s*(-?\d+(?:\.\d+)?)", v, re.I)}
                break
    b = blades.get(_norm(meta.get("blade", "")).strip())
    if b:
        meta["blade_radius_cm"], meta["blade_height_cm"] = b["radius_cm"], b["height_cm"]
    return meta, blades


def crosscheck(date, meta, logs, n_tests, rep):
    """Compare metadata txt / media xlsx (+ Notes) / UMA xlsx for one day.
    Returns dict(stations, media, segs, medium_raw, flags). Several media / stations in a day are NORMAL (the Notes
    give the time ranges); they only raise a flag when the Notes are missing or inconsistent."""
    flags = []

    def flag(msg):
        flags.append(msg)
        rep.crosschecks.append(f"{date}: {msg}")
        rep.warn(date, f"CROSS-CHECK: {msg}")

    txt_stations = _stations_in(meta.get("station")) if meta else []
    txt_med = meta.get("medium_raw") if meta else None
    xl_st = sorted(logs["uma"].get(date, []))
    xl_med = logs["media"].get(date)
    notes = " ".join(logs["notes"].get(date, []))
    # ---- station: metadata txt vs UMA log
    if logs["uma_file"] is not None:
        if txt_stations and xl_st:
            missing = sorted(set(txt_stations) - set(xl_st))
            if missing:
                flag(f"STATION MISMATCH: metadata txt says {', '.join(txt_stations)}, "
                     f"{logs['uma_file'].name} says {', '.join(xl_st)}")
        if txt_stations and not xl_st and n_tests:
            flag(f"metadata txt gives station(s) {', '.join(txt_stations)} but the UMA log has none for this day")
    # ---- medium: metadata txt vs media column
    txt_set = _split_media(txt_med) if txt_med else []
    xl_set = _split_media(xl_med[0]) if xl_med else []
    for src, toks in (("metadata txt", txt_set), ("media xlsx", xl_set)):
        for t in toks:
            if t not in MEDIA:
                flag(f"unrecognised medium '{t}' in {src} (typo?)")
    if logs["media_file"] is not None and xl_med is None and n_tests:
        flag(f"day not in {logs['media_file'].name}")
    if txt_set and xl_set and set(txt_set) != set(xl_set):
        part = " (the txt lists only part of the media of the day)" if set(txt_set) < set(xl_set) else ""
        flag(f"MEDIUM MISMATCH: metadata txt says '{txt_med}' ({', '.join(txt_set)}), "
             f"{logs['media_file'].name} says '{xl_med[0]}' ({', '.join(xl_set)}){part}")
    if meta is None and n_tests and not xl_med:
        flag("no medium information at all (no metadata txt, not in the media log)")
    # ---- Notes: time ranges
    segs = parse_notes(notes, flag) if notes else []
    if notes and not segs:
        flag(f"Notes could not be read (expected 'medium [station]: from 10am to 2:40pm'): '{notes[:80]}'")
    seg_media = {x["medium"] for x in segs}
    if segs and xl_set and seg_media != set(xl_set):
        flag(f"NOTES vs MEDIUM column: notes give {sorted(seg_media)}, the Medium column gives {sorted(set(xl_set))}")
    seg_st = {x["station"] for x in segs if x["station"]}
    if seg_st and xl_st:
        for st in sorted(seg_st - set(xl_st)):
            flag(f"STATION {st} is in the notes but not in {logs['uma_file'].name} for this day")
        for st in sorted(set(xl_st) - seg_st):
            flag(f"STATION {st} is in {logs['uma_file'].name} but not in the notes of this day")
    for a, b in zip(segs, segs[1:]):
        if b["t0"] < a["t1"]:
            flag(f"notes: time ranges overlap ('{a['text']}' / '{b['text']}')")
    # ---- what the day contains
    first, second = ((txt_set, xl_set) if CROSSCHECK_PRIORITY == "txt" else (xl_set, txt_set))
    ok = lambda toks: bool(toks) and all(t in MEDIA for t in toks)
    media = list(dict.fromkeys(x["medium"] for x in segs)) or (first if ok(first) else second if ok(second) else (first or second))
    stations = xl_st or txt_stations or sorted(seg_st)
    if logs["uma_file"] is not None and not xl_st and n_tests and any(m in STATION_MEDIA for m in media):
        flag(f"day not in {logs['uma_file'].name} although it has {'/'.join(m for m in media if m in STATION_MEDIA)} tests (no station)")
    if (len(media) > 1 or len(stations) > 1) and not segs and n_tests:
        flag(f"several media {media} / stations {stations} that day but no readable time ranges in the Notes: "
             "tests cannot be assigned (use MEDIUM_OVERRIDES / STATION_OVERRIDES)")
    return dict(stations=stations, media=media, segs=segs, medium_raw=txt_med or (xl_med[0] if xl_med else None), flags=flags)


def _minutes_of(ts):
    return ts.hour * 60 + ts.minute + ts.second / 60


def resolve_test(cc, e, rep):
    """Medium and station of one test from the day's time ranges (Notes), else from the day's single medium/station."""
    where = f"{e.date} {e.time}"
    t0 = pd.Timestamp(e.log_start) if e.log_start is not None and pd.notna(e.log_start) else e.start
    m0, m1 = _minutes_of(t0), _minutes_of(e.window[1])
    seg = next((x for x in cc["segs"] if x["t0"] <= m0 < x["t1"]), None)
    if seg is None:
        seg = next((x for x in cc["segs"] if x["t0"] - NOTES_TOL_MIN <= m0 <= x["t1"] + NOTES_TOL_MIN), None)
    medium = station = None
    if seg:
        medium, station = seg["medium"], seg["station"]
        if m1 > seg["t1"] + NOTES_TOL_MIN:
            e.issues.append(f"test ends at {_hhmm(m1)}, after the end of its notes range ({seg['text']}): check")
            rep.warn(where, f"test runs past the end of its notes range ({seg['text']})")
    else:
        if cc["segs"]:
            rngs = "; ".join(f"{x['medium']} {x['station'] or ''} {_hhmm(x['t0'])}-{_hhmm(x['t1'])}" for x in cc["segs"])
            e.issues.append(f"start {_hhmm(m0)} is outside every time range of the notes ({rngs})")
            rep.warn(where, f"start {_hhmm(m0)} is not covered by the notes ({rngs})")
        elif len(cc["media"]) == 1:
            medium = cc["media"][0]
    ov = MEDIUM_OVERRIDES.get((e.date, e.time)) or MEDIUM_OVERRIDES.get(e.date)
    if ov:
        medium = norm_medium(ov)
    if medium is None and len(cc["media"]) == 1 and not cc["segs"]:
        medium = cc["media"][0]
    if station is None and medium in STATION_MEDIA:
        if len(cc["stations"]) == 1:
            station = cc["stations"][0]
        elif len(cc["stations"]) > 1:
            e.issues.append(f"several stations that day {cc['stations']} and none in the notes for this time")
    ovs = STATION_OVERRIDES.get((e.date, e.time)) or STATION_OVERRIDES.get(e.date)
    if ovs:
        station = str(ovs)
    e.medium, e.station, e.medium_raw = medium or "unknown", station, cc["medium_raw"]
    if medium is None:
        e.issues.append("medium unknown")


def find_meta(day_dir, root, date):
    ddmmyyyy = datetime.strptime(date, "%Y-%m-%d").strftime("%d%m%Y")
    c = sorted(day_dir.glob("metadata*.txt")) + sorted(root.glob(f"metadata{ddmmyyyy}*.txt"))
    return next((p for p in c if ddmmyyyy in p.name), c[0] if c else None)


# ============================================================================ input log -> test description
DESC_DEFAULT = dict(test_type="unknown", direction=None, rotation=None, freq_hz=None, mean_rpm=None, amp_rpm=None,
                    amp_start_rpm=None, amp_end_rpm=None, sweep_shape=None, hold_rpm=None, speed_min_rpm=None,
                    speed_max_rpm=None, n_levels=None, plateau_s=None, ramp_shape=None, cycles=None, duration_s=None,
                    log_start=None, log_stop=None, n_commands=0, levels=None)


def read_input_log(path):
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    for c in ("t_unix_s", "datetime", "command", "ok", "source", "note"):
        if c not in df.columns:
            df[c] = ""
    df["t"] = pd.to_numeric(df["t_unix_s"], errors="coerce")
    df["dt"] = pd.to_datetime(df["datetime"], errors="coerce")
    df["ok"] = df["ok"].astype(str).str.strip().isin(["1", "True", "true"])
    return df


def _args(cmd):
    try:
        return [float(x) for x in cmd[1:].split(",") if x.strip()]
    except ValueError:
        return []


def _infer_slope(ex):
    """rpm per % duty from a note that gives both (sine / sweeps), when the .edc is not available."""
    for x in ex:
        a = _args(x["c"])
        m = SINE_RE.search(x["n"])
        if m and x["c"][0] == "o" and len(a) >= 2 and a[1] > 0:
            return abs(float(m[2])) / a[1]
        m = SWEEP_RE.search(x["n"])
        if m and x["c"][0] == "r" and len(a) >= 3:
            for rpm, duty in ((float(m[1]), a[1]), (float(m[2]), a[2])):
                if duty > 0 and rpm:
                    return abs(rpm) / duty
        m = STEPPED_RE.search(x["n"])
        if m and x["c"][0] == "o" and len(a) >= 2 and a[1] > 0:
            return abs(float(m[4])) / a[1]
    return None


def _spacing(vals):
    """'linear' / 'log' / 'irregular' for a sorted list of positive plateau values."""
    v = np.asarray(sorted(set(np.round(vals, 3))), float)
    if len(v) < 3:
        return None
    if np.allclose(np.diff(v), np.diff(v)[0], rtol=0.05, atol=1e-6):
        return "linear"
    if v[0] > 0 and np.allclose(v[1:] / v[:-1], (v[1] / v[0]), rtol=0.05):
        return "log"
    return "irregular"


def _direction(seq):
    s = np.diff(np.asarray(seq, float))
    if not len(s) or np.allclose(s, 0):
        return "flat"
    return "up" if np.all(s >= -1e-9) else "down" if np.all(s <= 1e-9) else "updown"


def describe_test(log, cfg):
    """Input log DataFrame (or None) -> (dict of test fields, list of issue strings)."""
    d, issues = dict(DESC_DEFAULT), []
    if log is None or log.empty:
        return d, ["no input log: test type unknown"]
    cmds, start, stop, mode, failed = [], None, None, None, 0
    for r in log.itertuples(index=False):
        c, n, src = str(r.command).strip(), str(r.note).strip(), str(r.source).strip()
        if c.startswith("["):
            if "start" in c:
                start = r.dt
            elif "stop" in c:
                stop = r.dt
            continue
        if src == "session" and n in SETUP_NOTES:
            if n == "arduino mode":
                mode = c
            continue
        if not r.ok:
            failed += 1
            continue
        cmds.append(dict(dt=r.dt, c=c, n=n))
    d["log_start"], d["log_stop"] = start, stop
    if failed:
        issues.append(f"{failed} command(s) not delivered (ok = 0)")
    rest_duty = cfg.rest_duty if cfg is not None else (10.0 if mode == "mc" else 50.0)
    is_rest_cmd = lambda x: x["c"] == "x" or x["n"].startswith("rest")
    ex = [x for x in cmds if not is_rest_cmd(x) and x["c"][:1] in ("c", "o", "r")]
    d["n_commands"] = len(ex)
    if not ex:
        d["test_type"] = "rest"
        return d, ["no motor command: motor at rest during the acquisition"]
    t_first = ex[0]["dt"]
    ends = [x["dt"] for x in cmds if is_rest_cmd(x) and pd.notna(x["dt"]) and pd.notna(t_first) and x["dt"] > t_first]
    t_end = min(ends) if ends else stop
    if pd.notna(t_first) and t_end is not None and pd.notna(t_end):
        d["duration_s"] = (t_end - t_first).total_seconds()
    slope = cfg.slope if cfg is not None else _infer_slope(ex)
    if cfg is None and slope is None:
        issues.append("no .edc and no rpm in the notes: speeds not available (duty only)")

    def rpm(duty):
        if cfg is not None:
            return cfg.value_at(duty)
        return (duty - rest_duty) * slope if slope else None

    fam = {x["c"][0] for x in ex}
    if len(fam) > 1:
        issues.append(f"mixed command types {sorted(fam)}: classified as '{'sweep' if 'r' in fam else 'sine' if 'o' in fam else 'program'}'")

    if "r" in fam:                                                         # ---- continuous amplitude sweep
        x = next(x for x in ex if x["c"][0] == "r")
        m, a = SWEEP_RE.search(x["n"]), _args(x["c"])
        if m:
            lo, hi, shape, mean, f, T = float(m[1]), float(m[2]), m[3], float(m[4]), float(m[5]), float(m[6])
        elif len(a) >= 5:
            lo = a[1] * slope if slope else None
            hi = a[2] * slope if slope else None
            mean = rpm(a[0])
            f, T, shape = a[3], a[4], "log" if len(a) > 5 and a[5] == 1 else "lin"
        else:
            issues.append("sweep command not understood")
            d["test_type"] = "sweep"
            return d, issues
        d.update(test_type="sweep", amp_start_rpm=lo, amp_end_rpm=hi, sweep_shape={"lin": "linear", "log": "log"}[shape],
                 mean_rpm=mean, freq_hz=f, duration_s=T, cycles=T * f,
                 direction=None if lo is None else "up" if hi > lo else "down" if hi < lo else "flat")
        return d, issues

    if "o" in fam and any(STEPPED_RE.search(x["n"]) for x in ex):          # ---- stepped amplitude sweep
        ms = [STEPPED_RE.search(x["n"]) for x in ex if STEPPED_RE.search(x["n"])]
        amps = [float(m[4]) for m in ms]
        d.update(test_type="sweep_stepped", amp_start_rpm=amps[0], amp_end_rpm=amps[-1], mean_rpm=float(ms[0][3]),
                 freq_hz=float(ms[0][5]), cycles=float(ms[0][6]), direction=_direction(amps) if len(amps) > 1 else "flat",
                 sweep_shape=_spacing(amps), n_levels=len(amps))
        d["duration_s"] = d["cycles"] / d["freq_hz"] * len(amps) if d["freq_hz"] else d["duration_s"]
        return d, issues

    if "o" in fam:                                                         # ---- sine (firmware)
        x = next(x for x in ex if x["c"][0] == "o")
        m, a = SINE_RE.search(x["n"]), _args(x["c"])
        d["test_type"] = "sine"
        if m:
            d.update(mean_rpm=float(m[1]), amp_rpm=abs(float(m[2])), freq_hz=float(m[3]), cycles=float(m[4]))
            d["duration_s"] = d["cycles"] / d["freq_hz"] if d["freq_hz"] else d["duration_s"]
        elif len(a) >= 3:
            d.update(mean_rpm=rpm(a[0]), amp_rpm=a[1] * slope if slope else None, freq_hz=a[2])
        else:
            issues.append("sine command not understood")
        return d, issues

    levels = []                                                            # ---- constant duties: hold / steps / ramp
    for i, x in enumerate(cmds):
        if x["c"][0] != "c":
            continue
        duty = (_args(x["c"]) or [None])[0]
        if duty is None:
            continue
        nxt = cmds[i + 1]["dt"] if i + 1 < len(cmds) else (t_end if t_end is not None else stop)
        dur = (nxt - x["dt"]).total_seconds() if pd.notna(nxt) and pd.notna(x["dt"]) else None
        levels.append(dict(duty=duty, dur=dur, note=x["n"], rest=abs(duty - rest_duty) < REST_TOL_DUTY or x["n"].startswith("rest")))
    work = [l for l in levels if not l["rest"]]
    if any(x["n"].startswith("stream sine") for x in ex):
        duties = [l["duty"] for l in work]
        d["test_type"] = "sine"
        if duties and slope:
            d["amp_rpm"] = 0.5 * (max(duties) - min(duties)) * abs(slope)
            d["mean_rpm"] = rpm(0.5 * (max(duties) + min(duties)))
        issues.append("stream-mode sine: frequency not written in the log")
        return d, issues
    if not work:
        d["test_type"] = "rest"
        return d, issues
    vals = [rpm(l["duty"]) if slope or cfg is not None else l["duty"] - rest_duty for l in work]
    for l, v in zip(work, vals):
        l["val"] = v
    uniq = sorted(set(round(v, 1) for v in vals))
    ramp_kind = any(l["note"].startswith("ramp") for l in work)
    d["test_type"] = "hold" if len(uniq) == 1 else "ramp" if (ramp_kind or len(uniq) >= 3) else "steps"
    mag = [abs(v) for v in vals]
    d.update(n_levels=len(uniq), speed_min_rpm=float(min(vals)), speed_max_rpm=float(max(vals)),
             rotation="positive" if min(vals) > 0 else "negative" if max(vals) < 0 else "both")
    durs = [l["dur"] for l in work if l["dur"] is not None]
    if durs:
        d["plateau_s"] = float(np.median(durs))
        d["duration_s"] = float(np.sum(durs))
    if d["test_type"] == "hold":
        d["hold_rpm"] = float(vals[0])
    else:
        d["direction"] = _direction(mag)
        d["ramp_shape"] = "continuous" if ramp_kind else _spacing(mag)
    d["levels"] = [(l["val"], l["dur"]) for l in work]
    return d, issues


# ============================================================================ motor configuration (.edc)
def load_cfg(edc_txt, root, cache, rep, where):
    """edc_used.txt -> (MotorConfig or None, edc name, edc path or None). Cached by file name."""
    try:
        lines = [l.strip() for l in Path(edc_txt).read_text(encoding="utf-8", errors="replace").splitlines() if l.strip()]
    except OSError as exc:
        rep.warn(where, f"edc_used.txt unreadable ({exc})")
        return None, None, None
    name = next((l for l in lines if not l.startswith("#")), None)
    full = next((l.split(":", 1)[1].strip() for l in lines if l.lower().startswith("# full path:")), None)
    if not name:
        return None, None, None
    if name in cache:
        return cache[name][0], name, cache[name][1]
    cands = [Path(full)] if full else []
    for base in [root.parent, root] + [Path(p) for p in EDC_SEARCH_FOLDERS]:
        for pat in (name, f"*/{name}"):
            cands += list(base.glob(pat))
    path = next((c for c in cands if c.is_file()), None)
    cfg = None
    if path is None:
        rep.warn(where, f".edc '{name}' not found (rpm values come from the log notes when possible)")
    else:
        try:
            import chipiron_helpers_v4 as cm
            cfg = cm.MotorConfig.from_edc_object(load_edc_quiet(path), path)
        except Exception as exc:
            rep.warn(where, f".edc '{name}' could not be read ({exc})")
    cache[name] = (cfg, path)
    return cfg, name, path


# ============================================================================ the Experiment object
FLAT_FIELDS = ["date", "time", "station", "medium", "test_type", "direction", "rotation", "freq_hz", "mean_rpm", "amp_rpm",
               "amp_start_rpm", "amp_end_rpm", "sweep_shape", "hold_rpm", "speed_min_rpm", "speed_max_rpm", "n_levels",
               "plateau_s", "ramp_shape", "cycles", "duration_s", "blade", "blade_radius_cm", "blade_height_cm",
               "grain_size", "sampling_height", "is_reference", "ref_type", "test_key", "folder"]


@dataclass
class Experiment:
    date: str
    time: str
    folder: Path
    start: pd.Timestamp                                  # from the folder name (phone clock)
    station: Any = None
    medium: str = "unknown"
    medium_raw: Any = None
    meta: dict = field(default_factory=dict, repr=False)
    blade: Any = None
    blade_radius_cm: Any = None
    blade_height_cm: Any = None
    # files
    usb_file: Any = None                                 # first usb part (kept for compatibility)
    usb_files: list = field(default_factory=list, repr=False)   # ALL usb parts, in order (Usb-1, Usb-2, ...)
    sensor_files: dict = field(default_factory=dict, repr=False)
    input_log_file: Any = None
    inputs: Any = field(default=None, repr=False)        # DataFrame of the commands sent
    edc_name: Any = None
    edc_path: Any = None
    cfg: Any = field(default=None, repr=False)
    sampling_period_ms: Any = None
    usb_duration_s: Any = None
    # test description
    test_type: str = "unknown"
    direction: Any = None
    rotation: Any = None
    freq_hz: Any = None
    mean_rpm: Any = None
    amp_rpm: Any = None
    amp_start_rpm: Any = None
    amp_end_rpm: Any = None
    sweep_shape: Any = None
    hold_rpm: Any = None
    speed_min_rpm: Any = None
    speed_max_rpm: Any = None
    n_levels: Any = None
    plateau_s: Any = None
    ramp_shape: Any = None
    cycles: Any = None
    duration_s: Any = None
    levels: Any = field(default=None, repr=False)
    n_commands: int = 0
    log_start: Any = None
    log_stop: Any = None
    clock_offset_s: Any = None                           # log /start minus phone folder time [s]
    # temperature
    temp_file: Any = None
    temp_match: Any = None
    temp_delta_start_s: Any = None
    temp_columns: list = field(default_factory=list)
    temp_mean_c: dict = field(default_factory=dict)
    temp_df: Any = field(default=None, repr=False)
    # images
    grain_images: list = field(default_factory=list, repr=False)
    height_images: list = field(default_factory=list, repr=False)
    image_results_file: Any = None
    image_results: Any = field(default=None, repr=False)
    grain_size: Any = None
    sampling_height: Any = None
    # reference logic
    is_reference: bool = False                           # test of a REFERENCE_DAYS day
    ref_type: Any = None                                 # empty_ref / ninja_ref / saltwater_ref for reference tests
    test_key: Any = None                                 # text signature of the INPUT (same key = same protocol)
    ref_available: dict = field(default_factory=dict, repr=False)   # {ref_type: [matching reference Experiments]}
    issues: list = field(default_factory=list)
    _cache: dict = field(default_factory=dict, repr=False)

    # ---- windows / flat info
    @property
    def usb_window(self):
        return self.start, self.start + pd.Timedelta(seconds=self.usb_duration_s or 0)

    @property
    def log_window(self):
        if self.log_start is not None and self.log_stop is not None and pd.notna(self.log_start) and pd.notna(self.log_stop):
            return pd.Timestamp(self.log_start), pd.Timestamp(self.log_stop)
        return None

    @property
    def window(self):
        return (self.log_window or self.usb_window) if TEMP_WINDOW_SOURCE == "log" else (self.usb_window or self.log_window)

    def info(self):
        d = {k: getattr(self, k) for k in FLAT_FIELDS}
        d.update(has_usb=self.usb_file is not None, n_usb_parts=len(self.usb_files), has_input_log=self.inputs is not None,
                 has_temperature=self.temp_df is not None, temp_match=self.temp_match,
                 n_grain_images=len(self.grain_images), n_height_images=len(self.height_images),
                 has_image_results=self.image_results is not None, n_issues=len(self.issues),
                 refs_available=",".join(k for k, v in self.ref_available.items() if v))
        return d

    def params_str(self):
        f = lambda v, u="": "?" if v is None else f"{v:g}{u}"
        t = self.test_type
        if t == "sine":
            return f"mean {f(self.mean_rpm)} amp {f(self.amp_rpm)} rpm, {f(self.freq_hz)} Hz, {f(self.duration_s)} s"
        if t == "sweep":
            return f"amp {f(self.amp_start_rpm)}->{f(self.amp_end_rpm)} rpm ({self.sweep_shape}, {self.direction}), {f(self.freq_hz)} Hz, {f(self.duration_s)} s"
        if t == "sweep_stepped":
            return f"{self.n_levels} amps {f(self.amp_start_rpm)}->{f(self.amp_end_rpm)} rpm ({self.direction}), {f(self.freq_hz)} Hz"
        if t == "hold":
            return f"{f(self.hold_rpm)} rpm for {f(self.duration_s)} s"
        if t in ("ramp", "steps"):
            return (f"{f(self.speed_min_rpm)}..{f(self.speed_max_rpm)} rpm, {self.n_levels} levels, {self.direction}, "
                    f"{self.ramp_shape}, {f(self.plateau_s)} s/level")
        return ""

    # ---- lazy raw data
    def usb_raw(self):
        """Raw usb rows of ALL the parts of the recording (the phone splits long recordings: Usb-1, Usb-2, ...),
        concatenated in order with ONE continuous t_s (as Chipiron_plotting_v5.read_usb_rows_all)."""
        files = self.usb_files or ([self.usb_file] if self.usb_file is not None else [])
        if not files:
            raise FileNotFoundError(f"{self.date} {self.time}: no usb file")
        dec = _decoder()
        dfs = []
        for f in files:
            try:
                d = dec.read_usb_rows(f, verbose=False)
            except TypeError:
                d = dec.read_usb_rows(f)
            if len(d):
                dfs.append(d)
        if not dfs:
            return dec.read_usb_rows(files[0], verbose=False)       # empty frame, same as before
        if len(dfs) == 1:
            return dfs[0]
        df = pd.concat(dfs, ignore_index=True)
        df["t_s"] = (df["t_us"] - df["t_us"].iloc[0]) / 1e6
        return df

    def usb(self, physical=False):
        """Decoded usb frame (one row per telemetry cycle): t_s, A1_adc, A2_adc, A1_V, A2_V (+ speed/current if physical)."""
        key = ("usb", physical)
        if key in self._cache:
            return self._cache[key]
        dec = _decoder()
        kw = {"vref": ADC_VREF}
        params = inspect.signature(dec.decode_usb).parameters
        pname = next((p for p in params if "period" in p), None)
        if pname and self.sampling_period_ms:
            kw[pname] = self.sampling_period_ms
        if "verbose" in params:
            kw["verbose"] = False
        df = dec.decode_usb(self.usb_raw(), **kw)
        if physical:
            if self.edc_path is None:
                raise FileNotFoundError("no .edc available for this test: cannot convert to physical units")
            load_edc_quiet(self.edc_path).apply_to_frame(df, wiring=WIRING)
        self._cache[key] = df
        return df

    def sensor(self, kind):
        """accelerometer / gyroscope / magnetic_field / gps -> DataFrame with t_s (placeholders removed)."""
        if kind not in self.sensor_files:
            raise FileNotFoundError(f"{self.date} {self.time}: no {kind} file")
        if ("s", kind) not in self._cache:
            df = pd.read_csv(self.sensor_files[kind], header=None, skipinitialspace=True)
            if df.shape[1] > 1 and df.iloc[:, -1].isna().all():
                df = df.iloc[:, :-1]
            names = SENSOR_COLS[kind]
            df.columns = names if df.shape[1] == len(names) else [f"col{i}" for i in range(df.shape[1])]
            df = df[df.iloc[:, 0] > 0].reset_index(drop=True)
            if len(df):
                df["t_s"] = (df.iloc[:, 0] - df.iloc[0, 0]) / 1e6
            self._cache[("s", kind)] = df
        return self._cache[("s", kind)]

    def temperature(self, crop=True, pad_s=0.0):
        """Matched Exelog data; crop=True keeps the test window (+/- pad_s). Adds t_s = s since the window start."""
        if self.temp_df is None:
            raise FileNotFoundError(f"{self.date} {self.time}: no temperature file matched")
        w0, w1 = self.window
        df = self.temp_df
        if crop:
            df = df[(df["time"] >= w0 - pd.Timedelta(seconds=pad_s)) & (df["time"] <= w1 + pd.Timedelta(seconds=pad_s))]
        df = df.copy()
        df["t_s"] = (df["time"] - w0).dt.total_seconds()
        return df

    def nearest_image(self, kind="grain"):
        """(path, seconds from the test start) of the image of that kind closest in time (needs EXIF time + Pillow)."""
        best = None
        for p in (self.grain_images if kind == "grain" else self.height_images):
            t = _image_time(p)
            if t is not None:
                dt = abs((t - self.start).total_seconds())
                if best is None or dt < best[1]:
                    best = (p, dt)
        return best

    def __repr__(self):
        return (f"Experiment({self.date} {self.time}, st {self.station}, {self.medium}, {self.test_type}"
                f"{': ' + self.params_str() if self.params_str() else ''})")


def _decoder():
    for name in ("chipiron_decode_v4", "chipiron_decode_v3"):
        try:
            return __import__(name)
        except ImportError:
            continue
    raise ImportError("chipiron_decode_v4.py / chipiron_decode_v3.py not found next to this script")


def _image_time(p):
    try:
        from PIL import Image
        ex = Image.open(p).getexif()
        s = ex.get_ifd(0x8769).get(36867) or ex.get(306)
        return pd.Timestamp(datetime.strptime(s, "%Y:%m:%d %H:%M:%S")) if s else None
    except Exception:
        return None


# ============================================================================ the set of experiments
def between(lo, hi):
    """Use in filter(): freq_hz=between(0.1, 1)  (inclusive)."""
    return lambda v: v is not None and not pd.isna(v) and lo <= v <= hi


def _match(v, w):
    """filter() comparison: callable | numbers (isclose, '301' == 301) | text."""
    if callable(w):
        return bool(w(v))
    try:
        return bool(np.isclose(float(v), float(w), rtol=1e-6, atol=1e-9))
    except (TypeError, ValueError):
        return str(v) == str(w)


# ============================================================================ reference tests
# fields of the INPUT that define "the same test", per test type (duration / cycles ignored for sine and hold)
REF_FIELDS = {
    "sine": ("mean_rpm", "amp_rpm", "freq_hz"),
    "sweep": ("mean_rpm", "amp_start_rpm", "amp_end_rpm", "sweep_shape", "freq_hz", "duration_s"),
    "sweep_stepped": ("mean_rpm", "amp_start_rpm", "amp_end_rpm", "sweep_shape", "n_levels", "freq_hz", "cycles"),
    "hold": ("hold_rpm",),
    "ramp": ("speed_min_rpm", "speed_max_rpm", "n_levels", "ramp_shape", "direction", "plateau_s"),
    "steps": ("speed_min_rpm", "speed_max_rpm", "n_levels", "direction"),   # plateau time ignored (a reference run can be stopped early)
}


def make_test_key(e):
    """Text signature of the input of a test, built from the input log: 'sine mean_rpm=0, amp_rpm=50, freq_hz=1'.
    None if the type is unknown; 'rest' for a recording without motor command."""
    if e.test_type == "rest":
        return "rest"
    fields = REF_FIELDS.get(e.test_type)
    if fields is None:
        return None
    bits = []
    for k in fields:
        v = getattr(e, k)
        bits.append(f"{k}={v:g}" if isinstance(v, (int, float, np.number)) and not isinstance(v, bool) else f"{k}={v}")
    return f"{e.test_type} " + ", ".join(bits)


AMP_FIELDS = ("mean_rpm", "amp_rpm", "amp_start_rpm", "amp_end_rpm")      # set by "o" / "r" commands (duty, 2 decimals)


def ref_tol(k, a, b):
    """Absolute tolerance of field k between tests a and b: 1 duty step of the command that sets it (never looser)."""
    slopes = [abs(x.cfg.slope) for x in (a, b) if x.cfg is not None and x.cfg.slope]
    slope = max(slopes) if slopes else REF_SLOPE_FALLBACK
    if k in AMP_FIELDS:
        base = max(REF_ATOL_DEFAULT, 1.05 * REF_DUTY_STEP_O * slope)
    elif k in ("hold_rpm", "speed_min_rpm", "speed_max_rpm"):
        base = max(REF_ATOL_DEFAULT, 1.05 * REF_DUTY_STEP_C * slope)
    else:
        base = REF_ATOL.get(k, REF_ATOL_DEFAULT)
    return base


def _field_diff(k, a, b):
    """None if the field k agrees, else |difference| / tolerance (> 1) ; 'str' for a different text / None."""
    va, vb = getattr(a, k), getattr(b, k)
    if va is None or vb is None:
        return None if va is vb else "str"
    if isinstance(va, str) or isinstance(vb, str):
        return None if va == vb else "str"
    tol = ref_tol(k, a, b) + REF_RTOL * abs(vb)
    d = abs(va - vb)
    return None if d <= tol else d / tol


def mismatch_fields(a, b):
    """Fields of the input that differ between tests a and b (empty list = same test)."""
    if a.test_type != b.test_type or a.test_type not in REF_FIELDS:
        return ["test_type"]
    return [k for k in REF_FIELDS[a.test_type] if _field_diff(k, a, b) is not None]


def ref_distance(a, b):
    """How far b is from a (0 = identical input): sum of |difference| / tolerance over the fields. Used to pick the best reference."""
    if a.test_type != b.test_type or a.test_type not in REF_FIELDS:
        return float("inf")
    s = 0.0
    for k in REF_FIELDS[a.test_type]:
        va, vb = getattr(a, k), getattr(b, k)
        if isinstance(va, (int, float, np.number)) and isinstance(vb, (int, float, np.number)):
            s += abs(va - vb) / (ref_tol(k, a, b) + REF_RTOL * abs(vb))
    return s


def near_fields(a, b):
    """Fields that differ a little (inside the tolerance but not identical): the reference is a NEAR match, not an exact one."""
    out = []
    for k in REF_FIELDS.get(a.test_type, ()):
        va, vb = getattr(a, k), getattr(b, k)
        if k in REF_NEAR_IGNORE or not isinstance(va, (int, float, np.number)) or not isinstance(vb, (int, float, np.number)):
            continue
        if abs(va - vb) > max(0.05, 1e-3 * abs(vb)):
            out.append(f"{k} {va:g} vs {vb:g}")
    return out


def usb_usable(e):
    """(ok, reason): does the usb csv of this test decode into at least one full telemetry cycle? Cached on the test.
    ImportError of the decoder is NOT counted as a failure (returns True, 'not checked')."""
    key = ("usb_usable",)
    if key in e._cache:
        return e._cache[key]
    if e.usb_file is None:
        res = (False, "no usb file")
    else:
        try:
            n_raw = len(e.usb_raw())
            if n_raw < 4:
                res = (False, f"only {n_raw} valid row(s) in the usb csv (placeholder rows removed)")
            else:
                d = e.usb(physical=False)
                n_full = int((d["kind"] == "full").sum()) if "kind" in d.columns else len(d)
                if len(d) == 0 or n_full == 0:
                    res = (False, f"{n_raw} rows but no full telemetry cycle decoded ({len(d)} cycle(s))")
                else:
                    res = (True, f"{n_raw} rows, {n_full} full cycle(s)")
        except ImportError as exc:
            res = (True, f"not checked ({exc})")
        except Exception as exc:
            res = (False, f"{type(exc).__name__}: {exc}")
    e._cache[key] = res
    return res


def ref_sort_key(e, r):
    """Order of the candidate references r for test e: usable usb first (PREFER_USABLE_REFS), then closest input, then chronological."""
    bad = PREFER_USABLE_REFS and not usb_usable(r)[0]
    return (bad, ref_distance(e, r), r.date, r.time)


def decodability_report(exps, only_refs=False, csv_path=None):
    """Try to decode the usb of every loaded test (incl. the reference days) and list the ones that cannot be decoded,
    split into reference / other tests, with the reason, and what it means for the matching."""
    pool = exps.ref_pool if getattr(exps, "ref_pool", None) is not None else exps
    items = [e for e in pool if e.is_reference or not only_refs]
    print("\n" + "=" * 100)
    print(f"USB DECODABILITY CHECK: {len(items)} test(s)")
    print("=" * 100)
    bad = []
    for i, e in enumerate(items, 1):
        ok, why = usb_usable(e)
        if not ok:
            bad.append((e, why))
            msg = f"usb not decodable: {why}"
            if msg not in e.issues:
                e.issues.append(msg)
        if i % 50 == 0:
            print(f"  ... {i}/{len(items)} checked")
    n_ref = sum(1 for e, _ in bad if e.is_reference)
    print(f"{len(items) - len(bad)} decodable, {len(bad)} NOT decodable ({n_ref} reference test(s), {len(bad) - n_ref} other test(s))")
    for title, is_ref in (("REFERENCE tests that cannot be decoded", True), ("OTHER tests that cannot be decoded", False)):
        sub = [(e, w) for e, w in bad if e.is_reference == is_ref]
        if not sub:
            continue
        print(f"\n--- {title} ({len(sub)}) ---")
        for e, why in sub:
            print(f"  {e.date} {e.time} {e.medium:<15} st:{e.station or '-':<4} {e.test_key or e.test_type}")
            print(f"        -> {why}")
            if is_ref:
                affected = [t for t in exps if not t.is_reference
                            and any(r is e for r in t.ref_available.get(e.ref_type, []))]
                rescued = sum(1 for t in affected if any(usb_usable(r)[0] for r in t.ref_available[e.ref_type]))
                print(f"        {len(affected)} test(s) list it as {e.ref_type}: {rescued} have another USABLE reference with the same input, "
                      f"{len(affected) - rescued} have none (their {e.ref_type} comparison is impossible)")
    if csv_path and bad:
        pd.DataFrame([dict(date=e.date, time=e.time, medium=e.medium, station=e.station, is_reference=e.is_reference,
                           ref_type=e.ref_type, test_key=e.test_key, usb_file=str(e.usb_file), reason=w) for e, w in bad]
                     ).to_csv(csv_path, index=False)
        print(f"\nList written: {Path(csv_path).resolve()}")
    print("=" * 100)
    return [e for e, _ in bad]


def same_test(a, b):
    """True if a and b received the same input (same test type, same parameters within the tolerances above)."""
    return not mismatch_fields(a, b)


def mark_references(exps, rep):
    """Set is_reference / ref_type / test_key on the tests of one day (medium already resolved)."""
    for e in exps:
        e.test_key = make_test_key(e)
        rt = REFERENCE_DAYS.get(e.date)
        if rt is None or (e.date, e.time) in REFERENCE_EXCLUDE:
            continue
        e.is_reference, e.ref_type = True, rt
        want = REF_EXPECTED_MEDIUM.get(rt)
        if want and e.medium not in (want, "unknown"):
            e.issues.append(f"reference day {rt}: medium is '{e.medium}', expected '{want}' (typo in the logs?)")
            rep.warn(f"{e.date} {e.time}", f"CROSS-CHECK: {rt} day but medium '{e.medium}' (expected '{want}')")
        if e.test_key is None:
            e.issues.append("reference test whose input could not be identified (cannot be matched)")


def assign_references(exps, pool):
    """For every non-reference test: ref_available = {ref_type: [reference tests of the pool with the same input]}."""
    refs = [r for r in pool if r.is_reference]
    for e in exps:
        e.ref_available = {}
        if e.is_reference:
            continue
        for rt in REF_MODES[1:]:
            e.ref_available[rt] = sorted((r for r in refs if r.ref_type == rt and same_test(e, r)),
                                         key=lambda r: ref_sort_key(e, r))   # usable usb first, then closest input


def ref_tag(e):
    """Short tag for the summary line."""
    if e.is_reference:
        return f"[REF {e.ref_type}]"
    ok = [k.replace("_ref", "") for k, v in e.ref_available.items() if v]
    return "refs:" + (",".join(ok) if ok else "-")


def ref_summary(exps):
    """Summary block: what the reference days contain and how many other tests have a matching reference."""
    pool = exps.ref_pool if getattr(exps, "ref_pool", None) is not None else exps
    L = ["", "REFERENCE TESTS (REFERENCE_DAYS at the top of the script):"]
    for day, rt in REFERENCE_DAYS.items():
        r = [e for e in pool if e.date == day and e.is_reference]
        if not r:
            L.append(f"  {rt:<14} {day}: NO TEST LOADED (day folder missing, or not in --days)")
            continue
        keys = list(dict.fromkeys(e.test_key for e in r))
        L.append(f"  {rt:<14} {day}: {len(r)} test(s), {len(keys)} distinct input(s)")
        L += [f"      {k}   x{sum(1 for e in r if e.test_key == k)}" for k in keys]
    tests = [e for e in exps if not e.is_reference and e.test_type not in ("unknown", "rest")]
    L.append(f"  Non-reference tests with a matching reference (out of {len(tests)}):")
    for rt in REF_MODES[1:]:
        L.append(f"      {rt:<14} {sum(1 for e in tests if e.ref_available.get(rt))}")
    none = [e for e in tests if not any(e.ref_available.values())]
    L.append(f"  Tests with NO reference at all: {', '.join(f'{e.date} {e.time}' for e in none[:25]) or 'none'}"
             + (" ..." if len(none) > 25 else ""))
    return L


class ExperimentSet(list):
    """List of Experiment with .table(), .filter(), .references(), .find_reference(), .pairs() and the summary text."""
    summary_text = ""
    days = None
    ref_pool = None          # every test loaded (incl. the reference days, even when --days excludes them)

    def references(self, ref=None):
        """Reference tests (all, or of one ref_type: empty_ref / ninja_ref / saltwater_ref)."""
        if ref is not None and ref not in REF_MODES[1:]:
            raise ValueError(f"ref must be one of {REF_MODES[1:]}")
        pool = self.ref_pool if self.ref_pool is not None else self
        out = ExperimentSet(r for r in pool if r.is_reference and (ref is None or r.ref_type == ref))
        out.summary_text, out.days, out.ref_pool = self.summary_text, self.days, self.ref_pool
        return out

    def find_references(self, e, ref):
        """All reference tests of that ref_type with the same input as e: usable usb first (PREFER_USABLE_REFS), then closest input, then chronological."""
        return sorted((r for r in self.references(ref) if same_test(e, r)), key=lambda r: ref_sort_key(e, r))

    def find_reference(self, e, ref):
        """The first matching reference test (closest input, then chronological), or None. ref='no_ref' -> None."""
        if ref == "no_ref":
            return None
        m = self.find_references(e, ref)
        if len(m) > 1:
            print(f"  note: {len(m)} {ref} tests match {e.date} {e.time} ({e.test_key}): using {m[0].date} {m[0].time}; "
                  "use find_references() to get them all")
        return m[0] if m else None

    def pairs(self, ref="no_ref", **crit):
        """[(test, reference or None)] for the NON-reference tests of the set (after filter(**crit) if given).
        ref = no_ref | empty_ref | ninja_ref | saltwater_ref."""
        if ref not in REF_MODES:
            raise ValueError(f"ref must be one of {REF_MODES}")
        sel = self.filter(**crit) if crit else self
        out = [(e, self.find_reference(e, ref)) for e in sel if not e.is_reference]
        if ref != "no_ref":
            print(f"[pairs] {ref}: {sum(r is not None for _, r in out)} / {len(out)} test(s) have a reference")
        return out

    def table(self):
        return pd.DataFrame([e.info() for e in self])

    def filter(self, predicate=None, **crit):
        """crit: field=value | field=[values] | field=between(lo, hi) | field=callable. predicate(exp) -> bool optional."""
        out = ExperimentSet()
        for e in self:
            ok = predicate(e) if predicate else True
            info = e.info()
            for k, want in crit.items():
                v = info.get(k, getattr(e, k, None))
                if isinstance(want, (list, tuple, set)):
                    ok = ok and any(_match(v, w) for w in want)
                else:
                    ok = ok and _match(v, want)
            if ok:
                out.append(e)
        out.summary_text, out.days, out.ref_pool = self.summary_text, self.days, self.ref_pool
        return out


# ============================================================================ loading one day
def scan_temperature(day_dir, date, rep, mt):
    recs, folders = [], [p for p in day_dir.iterdir() if p.is_dir() and p.name.lower().startswith("exelog")]
    for f in EXTRA_TEMP_FOLDERS:
        if Path(f).is_dir():
            folders.append(Path(f))
    if not folders:
        return recs                                      # no temperature recorded that day: not an issue
    mt.TEMP_CLOCK_OFFSET_S = TEMP_CLOCK_OFFSET_S
    hint = datetime.strptime(date, "%Y-%m-%d")
    n_files = 0
    for fo in folders:
        for p in sorted(fo.rglob("*.csv")):
            n_files += 1
            try:
                df = mt.read_temperature(p, date_hint=hint, name_start=None)
                recs.append(dict(path=p, df=df, columns=df.attrs["columns"], column="|".join(df.attrs["columns"]),
                                 start=df["time"].min(), end=df["time"].max()))
            except Exception as exc:
                rep.warn(date, f"temperature file {p.name} unreadable ({exc})")
    return recs


def collect_images(day_dir):
    out = {"grain": [], "height": []}

    def kind_of(name):
        n = _norm(name)
        return "grain" if "grain" in n else "height" if re.search(r"height|vertical|sampling|hauteur", n) else None

    for p in sorted(day_dir.iterdir()):
        if p.is_dir() and not TIME_RE.match(p.name) and not p.name.lower().startswith("exelog"):
            for f in sorted(p.rglob("*")):
                if f.is_file() and f.suffix.lower() in IMAGE_EXTS:
                    k = kind_of(p.name) or kind_of(f.name)
                    if k:
                        out[k].append(f)
        elif p.is_file() and p.suffix.lower() in IMAGE_EXTS and kind_of(p.name):
            out[kind_of(p.name)].append(p)
    return out


def read_image_results(day_dir, date, stations, rep):
    """vertical_and_grain_size_day<xxxxxx>_station<xxx>.csv files of the day -> {station: result dict}."""
    files = [f for f in sorted(day_dir.rglob("vertical_and_grain_size_*.csv")) if IMG_CSV_RE.search(f.name)]
    out = {}
    d = datetime.strptime(date, "%Y-%m-%d")
    tokens = {d.strftime(f) for f in ("%d%m%y", "%y%m%d", "%d%m%Y", "%Y%m%d")}
    for f in files:
        mm = IMG_CSV_RE.search(f.name)
        st = mm["st"]
        if mm["day"] not in tokens:
            rep.warn(date, f"CROSS-CHECK: {f.name}: day token '{mm['day']}' does not match {date} (typo?)")
        if stations and st not in stations:
            rep.warn(date, f"CROSS-CHECK: {f.name}: station {st} is not among the stations of this day {stations} (typo?)")
        if st in out:
            rep.warn(date, f"two processed image csv for station {st}: keeping {out[st]['path'].name}, ignoring {f.name}")
            continue
        try:
            df = pd.read_csv(f, sep=None, engine="python")
        except Exception as exc:
            rep.warn(date, f"{f.name} unreadable ({exc})")
            continue
        res = dict(path=f, df=df, station=st, grain_col=None, height_col=None, grain_size=None, sampling_height=None)
        for key, rx in (("grain", GRAIN_COL_RE), ("height", HEIGHT_COL_RE)):
            for c in df.columns:
                v = pd.to_numeric(df[c], errors="coerce")
                if re.search(rx, _norm(c)) and v.notna().any():
                    res[f"{key}_col"] = c
                    res["grain_size" if key == "grain" else "sampling_height"] = float(v.mean())
                    break
        if res["grain_col"] is None or res["height_col"] is None:
            rep.warn(date, f"{f.name}: column(s) not recognised for "
                           f"{[k for k in ('grain', 'height') if res[k + '_col'] is None]} (columns: {list(df.columns)}); "
                           "the DataFrame is still attached as exp.image_results")
        out[st] = res
    return out


def load_test(tdir, date, root, meta, blades_ok, edc_cache, rep):
    t = tdir.name
    where = f"{date} {t}"
    files = [f for f in tdir.iterdir() if f.is_file()]
    csvs = [f for f in files if f.suffix.lower() == ".csv"]
    usb_parts = sorted((f for f in csvs if "usb" in f.name.lower() and not f.name.lower().startswith("input_")), key=usb_part_key)
    usb = usb_parts[0] if usb_parts else None
    sensors = {k: next((f for f in csvs if pat in f.name.lower()), None) for k, pat in SENSOR_PATTERNS.items()}
    sensors = {k: v for k, v in sensors.items() if v is not None}
    start = pd.Timestamp(f"{date} {t.replace('-', ':')}")
    e = Experiment(date=date, time=t, folder=tdir, start=start, usb_file=usb, usb_files=usb_parts, sensor_files=sensors)
    if usb is None:
        e.issues.append("no usb file")
        rep.warn(where, "no usb file")
    else:
        try:
            lo = hi = None
            for f in usb_parts:                          # duration over ALL the parts
                s_ = pd.to_numeric(pd.read_csv(f, header=None, usecols=[0]).iloc[:, 0], errors="coerce")
                s_ = s_[s_ > 0]
                if len(s_):
                    lo = s_.min() if lo is None else min(lo, s_.min())
                    hi = s_.max() if hi is None else max(hi, s_.max())
            e.usb_duration_s = float((hi - lo) / 1e6) if lo is not None else None
            if len(usb_parts) > 1:
                print(f"  {where}: usb split in {len(usb_parts)} parts, concatenated ({e.usb_duration_s:.0f} s)")
        except Exception as exc:
            e.issues.append(f"usb unreadable: {exc}")
            rep.warn(where, f"usb file unreadable ({exc})")
    p = tdir / "sampling_period_ms.txt"
    if p.is_file():
        e.sampling_period_ms = _f(p.read_text().strip())
    # edc
    edc_txt = tdir / "edc_used.txt"
    if edc_txt.is_file():
        e.cfg, e.edc_name, e.edc_path = load_cfg(edc_txt, root, edc_cache, rep, where)
    else:
        e.issues.append("no edc_used.txt")
    # input log -> test description
    logf = next(iter(sorted(tdir.glob("input_*.csv"))), None)
    log = None
    if logf is None:
        e.issues.append("no input log")
        rep.warn(where, "no input log: test type unknown")
    else:
        e.input_log_file = logf
        try:
            log = read_input_log(logf)
            e.inputs = log
        except Exception as exc:
            e.issues.append(f"input log unreadable: {exc}")
            rep.warn(where, f"input log unreadable ({exc})")
    desc, issues = describe_test(log, e.cfg)
    for k, v in desc.items():
        setattr(e, k, v)
    e.issues += issues
    if e.log_start is not None and pd.notna(e.log_start):
        e.clock_offset_s = (pd.Timestamp(e.log_start) - start).total_seconds()
    # day metadata
    if meta:
        for k in ("station", "blade", "blade_radius_cm", "blade_height_cm", "medium_raw"):
            setattr(e, k, meta.get(k))
        e.meta = meta
    return e                                             # station / medium are resolved in load_day (cross-checked)


def load_day(day_dir, root, rep, edc_cache, mt, logs):
    date = day_dir.name
    dinfo = dict(date=date, meta_file=None, station=None, medium=None, blade=None, n_tests=0, n_temp_files=0, checks=[],
                 temp_unmatched=[], n_grain=0, n_height=0, image_csv=None, missing=[])
    mp = find_meta(day_dir, root, date)
    meta = None
    if mp is None:
        rep.warn(date, "no metadata txt for this day")
        dinfo["missing"].append("metadata")
    else:
        try:
            meta, _ = parse_metadata(mp)
            dinfo.update(meta_file=mp.name, blade=meta.get("blade"))
        except Exception as exc:
            rep.warn(date, f"metadata unreadable ({exc})")
            dinfo["missing"].append("metadata (unreadable)")
    exps = []
    for tdir in sorted(p for p in day_dir.iterdir() if p.is_dir() and TIME_RE.match(p.name)):
        if not any(tdir.iterdir()):
            rep.warn(f"{date} {tdir.name}", "empty folder, ignored")
            continue
        try:
            exps.append(load_test(tdir, date, root, meta, True, edc_cache, rep))
        except Exception as exc:
            rep.warn(f"{date} {tdir.name}", f"could not be loaded ({type(exc).__name__}: {exc})")
    dinfo["n_tests"] = len(exps)
    if not exps:
        rep.warn(date, "no test folder found")
        dinfo["missing"].append("tests")
    # cross-check txt / xlsx logs, then medium and station of every test (from the Notes time ranges)
    cc = crosscheck(date, meta, logs, len(exps), rep)
    dinfo["checks"] = cc["flags"]
    dinfo["station"] = "/".join(cc["stations"]) or None
    dinfo["medium"] = (" | ".join(f"{x['medium']}{' ' + x['station'] if x['station'] else ''} "
                                  f"{_hhmm(x['t0'])}-{_hhmm(x['t1'])}" for x in cc["segs"])
                       or " + ".join(cc["media"]) or None)
    for e in exps:
        resolve_test(cc, e, rep)
        e.issues += [f"cross-check: {f}" for f in cc["flags"]]
    mark_references(exps, rep)
    # temperature
    if exps and mt is not None:
        try:
            temps = scan_temperature(day_dir, date, rep, mt)
            dinfo["n_temp_files"] = len(temps)
            if temps:
                mt.START_THRESHOLD_S = START_THRESHOLD_S
                wins = [e.window for e in exps]
                res = mt.match([dict(start=w[0], end=w[1]) for w in wins], temps)
                for i, (j, kind) in res.items():
                    e, tr = exps[i], temps[j]
                    e.temp_file, e.temp_match, e.temp_df, e.temp_columns = tr["path"], kind, tr["df"], tr["columns"]
                    e.temp_delta_start_s = (tr["start"] - wins[i][0]).total_seconds()
                    sub = tr["df"][(tr["df"]["time"] >= wins[i][0]) & (tr["df"]["time"] <= wins[i][1])]
                    e.temp_mean_c = {c: float(sub[c].mean()) for c in tr["columns"] if sub[c].notna().any()}
                    if not e.temp_mean_c:
                        e.issues.append("temperature file matched but has no sample inside the test window")
                used = {j for j, _ in res.values()}
                dinfo["temp_unmatched"] = [t["path"].name for j, t in enumerate(temps) if j not in used]
        except Exception as exc:
            rep.warn(date, f"temperature matching failed ({type(exc).__name__}: {exc})")
    elif exps:
        rep.warn(date, "match_temperature_usb.py not importable: temperature not linked")
    # images
    imgs = collect_images(day_dir)
    dinfo["n_grain"], dinfo["n_height"] = len(imgs["grain"]), len(imgs["height"])
    irs = read_image_results(day_dir, date, cc["stations"], rep)
    dinfo["image_csv"] = ", ".join(r["path"].name for r in irs.values()) or None
    for k, name in (("grain", "grain-size images"), ("height", "sampling-height images")):
        if not imgs[k]:
            dinfo["missing"].append(name)
    if imgs["grain"] or imgs["height"]:
        if not irs:
            dinfo["missing"].append("processed image csv (images not yet processed)")
        elif [x for x in cc["stations"] if x not in irs]:
            dinfo["missing"].append(f"processed image csv for station(s) {[x for x in cc['stations'] if x not in irs]}")
    for e in exps:
        e.grain_images, e.height_images = imgs["grain"], imgs["height"]
        ir = irs.get(e.station) if e.station else None
        if ir:
            e.image_results_file, e.image_results = ir["path"], ir["df"]
            e.grain_size, e.sampling_height = ir["grain_size"], ir["sampling_height"]
    return exps, dinfo


# ============================================================================ summary
def build_summary(root, exps, rep):
    L = ["=" * 100, f"CHIPIRON EXPEDITION LOAD SUMMARY  ({datetime.now():%Y-%m-%d %H:%M:%S})", f"root: {root}", "=" * 100]
    for date, di in rep.days.items():
        de = [e for e in exps if e.date == date]
        L += ["", f"--- {date}: {di['n_tests']} test(s) | station {di['station'] or '?'} | medium {di['medium'] or '?'} "
                  f"| blade {di['blade'] or '?'} | metadata {di['meta_file'] or 'MISSING'}",
              f"    images: grain {di['n_grain']}, sampling height {di['n_height']}, processed csv {di['image_csv'] or '-'}"
              f" | Exelog files {di['n_temp_files']}"
              + (f" (not matched: {', '.join(di['temp_unmatched'][:8])}{' ...' if len(di['temp_unmatched']) > 8 else ''})"
                 if di["temp_unmatched"] else "")]
        if di["missing"]:
            L.append(f"    missing: {', '.join(di['missing'])}")
        for f in di.get("checks", []):
            L.append(f"    CROSS-CHECK: {f}")
        for e in de:
            L.append(f"    {e.time}  {e.medium:<15} st:{e.station or '-':<4} {e.test_type:<13} {e.params_str():<62} usb:{'Y' if e.usb_file else '-'} "
                     f"log:{'Y' if e.inputs is not None else '-'} edc:{'Y' if e.cfg is not None else '-'} "
                     f"T:{e.temp_match or '-':<7} {ref_tag(e):<24}"
                     + (f" Tmean:{', '.join(f'{k.replace('Channel ', 'ch')}={v:.1f}' for k, v in e.temp_mean_c.items())}"
                        if e.temp_mean_c else ""))
            for iss in e.issues:
                L.append(f"          . {iss}")
    n = len(exps)
    cnt = lambda f: sum(1 for e in exps if f(e))
    L += ["", "=" * 100, f"TOTAL: {n} test(s) in {len(rep.days)} day(s)",
          f"  with usb {cnt(lambda e: e.usb_file)} | input log {cnt(lambda e: e.inputs is not None)} | .edc {cnt(lambda e: e.cfg is not None)}"
          f" | temperature {cnt(lambda e: e.temp_df is not None)} | grain images {cnt(lambda e: e.grain_images)}"
          f" | height images {cnt(lambda e: e.height_images)} | processed image csv {cnt(lambda e: e.image_results is not None)}"]
    for title, key in (("test types", "test_type"), ("media", "medium")):
        c = pd.Series([getattr(e, key) for e in exps]).value_counts()
        L.append(f"  {title}: " + ", ".join(f"{k} {v}" for k, v in c.items()))
    L.append(f"  WITHOUT input log (type unknown): {', '.join(f'{e.date} {e.time}' for e in exps if e.inputs is None) or 'none'}")
    L.append(f"  medium unknown: {', '.join(f'{e.date} {e.time}' for e in exps if e.medium == 'unknown') or 'none'}")
    L += ref_summary(exps)
    L += ["", f"CROSS-CHECKS metadata txt / media xlsx / UMA xlsx to verify ({len(rep.crosschecks)}):"]
    L += [f"  ?? {c}" for c in rep.crosschecks] or ["  all consistent"]
    if rep.warnings:
        L += ["", f"WARNINGS ({len(rep.warnings)}):"] + [f"  ! {w}" for w in rep.warnings]
    L.append("=" * 100)
    return "\n".join(L)


# ============================================================================ main entry point
def load_expedition(root=None, days=None, summary_folder=None, verbose=True, media_log=None, uma_log=None):
    """Scan root, return an ExperimentSet. Prints the summary (verbose) and saves it if summary_folder (or SUMMARY_FOLDER)."""
    root = Path(root or ROOT_FOLDER)
    if not root.is_dir():
        raise FileNotFoundError(f"root folder not found: {root}")
    try:
        import match_temperature_usb as mt
    except ImportError:
        mt = None
        print("! match_temperature_usb.py not found next to this script: temperature will not be linked")
    rep, exps, edc_cache = Report(verbose), ExperimentSet(), {}
    logs = read_day_logs(root, rep, media_log, uma_log)
    pool_extra = []                                        # reference days not in `days`: loaded for the pool only
    for day_dir in sorted(p for p in root.iterdir() if p.is_dir() and DATE_RE.match(p.name)):
        if days and day_dir.name not in days:
            if day_dir.name in REFERENCE_DAYS:
                if verbose:
                    print(f"[{day_dir.name}] reference day ({REFERENCE_DAYS[day_dir.name]}): loaded for the reference pool only")
                try:
                    pool_extra.extend(load_day(day_dir, root, rep, edc_cache, mt, logs)[0])
                except Exception as exc:
                    rep.warn(day_dir.name, f"reference day skipped ({type(exc).__name__}: {exc})")
            continue
        if verbose:
            print(f"[{day_dir.name}]")
        try:
            de, di = load_day(day_dir, root, rep, edc_cache, mt, logs)
        except Exception as exc:
            rep.warn(day_dir.name, f"day skipped ({type(exc).__name__}: {exc})")
            continue
        rep.days[day_dir.name] = di
        exps.extend(de)
    if not days:                                           # days of the xlsx logs that have no folder
        folders = {p.name for p in root.iterdir() if p.is_dir() and DATE_RE.match(p.name)}
        for key, name in (("media", logs["media_file"]), ("uma", logs["uma_file"])):
            for d in sorted(set(logs[key]) - folders):
                msg = f"day {d} is in {name.name} but has no folder in {root.name}"
                rep.crosschecks.append(f"{d}: {msg}")
                rep.warn(d, f"CROSS-CHECK: {msg}")
    folders_now = {p.name for p in root.iterdir() if p.is_dir() and DATE_RE.match(p.name)}
    for d, rt in REFERENCE_DAYS.items():
        if d not in folders_now:
            rep.warn(d, f"{rt} day has no folder in {root.name}: no reference of that type")
    exps.ref_pool = ExperimentSet(list(exps) + pool_extra)
    assign_references(exps, exps.ref_pool)
    exps.days = rep.days
    exps.summary_text = build_summary(root, exps, rep)
    if verbose:
        print("\n" + exps.summary_text)
    folder = summary_folder or SUMMARY_FOLDER
    if folder:
        folder = Path(folder).expanduser()
        folder.mkdir(parents=True, exist_ok=True)
        out = folder / f"expedition_summary_{datetime.now():%Y%m%d-%H%M%S}.txt"
        out.write_text(exps.summary_text, encoding="utf-8")
        print(f"Summary saved: {out.resolve()}")
    return exps


# %% Helpers (run once): console log, reference report, plotting engine ==========
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- console log (everything printed -> txt)
_CONSOLE = globals().get("_CONSOLE", {})          # survives a re-run of the file in the same Spyder console


class _Tee:
    """Writes to the real stream AND to the log file."""

    def __init__(self, stream, fh):
        self._stream, self._fh = stream, fh

    def write(self, s):
        self._stream.write(s)
        try:
            self._fh.write(s)
        except ValueError:                         # file already closed
            pass
        return len(s)

    def flush(self):
        self._stream.flush()
        try:
            self._fh.flush()
        except ValueError:
            pass

    def __getattr__(self, name):
        return getattr(self._stream, name)


def start_console_log(folder=None):
    """Everything printed from now on is also written to <folder>/console_log_<date-time>.txt (one file per session)."""
    if _CONSOLE.get("fh") is not None and not _CONSOLE["fh"].closed:
        print(f"[console log] already active: {_CONSOLE['path']}")
        return _CONSOLE["path"]
    folder = Path(folder or ROOT_FOLDER)
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f"console_log_{datetime.now():%Y%m%d-%H%M%S}.txt"
    fh = open(path, "w", encoding="utf-8", buffering=1)
    _CONSOLE.update(fh=fh, path=path, out=sys.stdout, err=sys.stderr)
    sys.stdout, sys.stderr = _Tee(sys.stdout, fh), _Tee(sys.stderr, fh)
    print(f"[console log] writing to {path}")
    return path


def stop_console_log():
    """Restore the normal console and close the log file."""
    if _CONSOLE.get("fh") is None or _CONSOLE["fh"].closed:
        return
    print(f"[console log] closed: {_CONSOLE['path']}")
    sys.stdout, sys.stderr = _CONSOLE["out"], _CONSOLE["err"]
    _CONSOLE["fh"].close()


# ---------------------------------------------------------------- reference report (debug)
def _fmt_diff(e, r, fields):
    return ", ".join(f"{k}: {getattr(e, k)} vs {getattr(r, k)}" for k in fields)


def closest_reference(e, pool, rt):
    """Explain why e has no reference of type rt: closest reference test of the same test_type and what differs."""
    if e.test_type in ("unknown", "rest"):
        return f"test type '{e.test_type}': cannot be matched"
    cands = [r for r in pool if r.is_reference and r.ref_type == rt and r.test_type == e.test_type]
    if not cands:
        return f"the {rt} day has no '{e.test_type}' test at all"
    best = min(cands, key=lambda r: ref_distance(e, r))
    return f"closest {best.date} {best.time} differs in [{_fmt_diff(e, best, mismatch_fields(e, best))}]"


def reference_report(exps, per_test=True, media=None):
    """Console report: for each medium, which tests found their reference (empty / ninja / saltwater) and which did not."""
    pool = exps.ref_pool if exps.ref_pool is not None else exps
    rtypes = REF_MODES[1:]
    print("\n" + "=" * 100)
    print("REFERENCE REPORT  (which test found which reference test)")
    print("=" * 100)
    for day, rt in REFERENCE_DAYS.items():
        n = sum(1 for r in pool if r.date == day and r.is_reference)
        print(f"  {rt:<14} = day {day}: {n} reference test(s) loaded")
    for medium in sorted({e.medium for e in exps}):
        if media and medium not in media:
            continue
        sub = [e for e in exps if e.medium == medium]
        tests = [e for e in sub if not e.is_reference]
        nref = len(sub) - len(tests)
        print(f"\n--- medium '{medium}': {len(sub)} test(s)"
              + (f", {nref} of them ARE reference tests (not matched against anything)" if nref else ""))
        if not tests:
            continue
        found = {rt: sum(1 for e in tests if e.ref_available.get(rt)) for rt in rtypes}
        print(f"    {len(tests)} test(s) to match -> found: " + ", ".join(f"{rt} {n}/{len(tests)}" for rt, n in found.items()))
        if per_test:
            for e in tests:
                cells = []
                for rt in rtypes:
                    m = e.ref_available.get(rt) or []
                    if m:
                        nf = near_fields(e, m[0])
                        cells.append(f"{rt[:-4]}: {m[0].date[5:]} {m[0].time}" + (f" (+{len(m) - 1})" if len(m) > 1 else "")
                                     + (f" ~NEAR[{'; '.join(nf)}]" if nf else ""))
                    else:
                        cells.append(f"{rt[:-4]}: -")
                print(f"    {e.date} {e.time} st:{e.station or '-':<4} {e.test_key or e.test_type}")
                print("          " + " | ".join(cells))
        none = [e for e in tests if not any(e.ref_available.values())]
        if none:
            print(f"    {len(none)} test(s) with NO reference of any type:")
            for e in none:
                print(f"      {e.date} {e.time} {e.test_key or e.test_type}")
                for rt in rtypes:
                    print(f"          {rt}: {closest_reference(e, pool, rt)}")
    print("=" * 100)


# ---------------------------------------------------------------- plotting engine
def test_onset_s(e):
    """Seconds between /start and the first motor command of the test (laptop clock), from the input log. None if unknown."""
    log = e.inputs
    if log is None or log.empty:
        return None
    t0, rest = None, (e.cfg.rest_duty if e.cfg is not None else None)
    for r in log.itertuples(index=False):
        c, n = str(r.command).strip(), str(r.note).strip()
        if c.startswith("["):
            if "start" in c:
                t0 = r.t
            continue
        if t0 is None or not r.ok or c[:1] not in ("c", "o", "r") or n.startswith("rest") or "rest" in n.split(" ")[:2]:
            continue
        if c[0] == "c" and rest is not None:
            a = _args(c)
            if a and abs(a[0] - rest) < REST_TOL_DUTY:
                continue
        return float(r.t - t0)
    return None


def load_run(e):
    """Decoded run of one test -> dict(t, speed, current, ylab, ilab) with t = 0 at the motor command. None if unusable."""
    try:
        d = e.usb(physical=True)
        sc = next(c for c in d.columns if c.startswith("A1_") and c.endswith("_rpm"))
        cc = next(c for c in d.columns if c.startswith("A2_") and c.endswith("_A"))
        lab = (r"$\omega$ (rpm)", r"$I$ (A)")
    except (StopIteration, FileNotFoundError, KeyError):
        try:
            d = e.usb(physical=False)
        except Exception as exc:
            print(f"  ! {e.date} {e.time}: cannot decode usb ({exc})")
            return None
        sc, cc, lab = "A1_V", "A2_V", ("$V_{A1}$ (V)", "$V_{A2}$ (V)")
        print(f"  ! {e.date} {e.time}: no .edc -> volts plotted (no Delta I in A)")
    except Exception as exc:
        print(f"  ! {e.date} {e.time}: cannot decode usb ({exc})")
        return None
    d = d[d["kind"] == "full"].dropna(subset=[sc, cc]) if "kind" in d.columns else d.dropna(subset=[sc, cc])
    on = test_onset_s(e)
    if on is None:
        print(f"  ! {e.date} {e.time}: motor command time unknown, t = 0 at the start of the acquisition")
        on = 0.0
    return dict(t=d["t_s"].to_numpy(float) - on, speed=d[sc].to_numpy(float), current=d[cc].to_numpy(float),
                ylab=lab[0], ilab=lab[1])


def plot_pair(e, ref, rtest, pre_s=3.0, post_s=2.0):
    """One figure: test (+ its reference). Returns the figure or None."""
    a = load_run(e)
    if a is None:
        return None
    b = load_run(rtest) if rtest is not None else None
    t_end = (e.duration_s or a["t"].max()) + post_s
    wa = (a["t"] >= -pre_s) & (a["t"] <= t_end)
    wb = None if b is None else (b["t"] >= -pre_s) & (b["t"] <= t_end)
    fig, ax = plt.subplots(2, 2, figsize=(13, 8), layout="constrained")
    for k, key, lab in ((0, "speed", a["ylab"]), (1, "current", a["ilab"])):
        if b is not None:
            ax[k, 0].plot(b["t"][wb], b[key][wb], color="0.6", lw=1.5, label=f"{ref} ({rtest.date} {rtest.time})")
        ax[k, 0].plot(a["t"][wa], a[key][wa], color="tab:blue", lw=1.5, label="test")
        ax[k, 0].set(xlabel=r"$t$ (s)", ylabel=lab)
        ax[k, 0].grid(alpha=0.3)
    ax[0, 0].legend(fontsize=8)
    if b is not None:
        ax[0, 1].scatter(b["speed"][wb], b["current"][wb], s=6, color="0.7", label=ref, rasterized=True)
    ax[0, 1].scatter(a["speed"][wa], a["current"][wa], s=6, color="tab:blue", label="test", rasterized=True)
    ax[0, 1].set(xlabel=a["ylab"], ylabel=a["ilab"])
    ax[0, 1].grid(alpha=0.3)
    if b is not None:
        o = np.argsort(b["t"])
        i_ref = np.interp(a["t"], b["t"][o], b["current"][o], left=np.nan, right=np.nan)    # same elapsed time
        dI = a["current"] - i_ref
        m = wa & np.isfinite(dI)
        sc = ax[1, 1].scatter(a["speed"][m], dI[m], c=a["t"][m], s=6, cmap="viridis", rasterized=True)
        ax[1, 1].axhline(0, color="k", ls="--", lw=1)
        ax[1, 1].set(xlabel=a["ylab"], ylabel=r"$\Delta I$ = test $-$ " + ref)
        ax[1, 1].grid(alpha=0.3)
        fig.colorbar(sc, ax=ax[1, 1], label=r"$t$ (s)")
        print(f"  {e.date} {e.time}: Delta I mean {np.nanmean(dI[m]):+.3g}, range {np.nanmin(dI[m]):+.3g}..{np.nanmax(dI[m]):+.3g} "
              f"({int(m.sum())} points; {int((wa & ~np.isfinite(dI)).sum())} outside the reference time range)")
    else:
        ax[1, 1].set_visible(False)
    fig.suptitle(f"{e.date} {e.time}  st {e.station or '-'}  {e.medium}  |  {e.test_key}  |  {ref}", fontsize=9)
    return fig


def select_tests(exps, filters):
    """Apply the filters one by one and print how many tests survive each step (and what exists when it reaches 0)."""
    crit = {k: v for k, v in (filters or {}).items() if v is not None}
    cur = exps
    print(f"[select] {len(cur)} test(s) loaded")
    for k, v in crit.items():
        nxt = cur.filter(**{k: v})
        print(f"[select]   {k} = {v!r}: {len(cur)} -> {len(nxt)}")
        if len(nxt) == 0:
            vals = pd.Series([str(e.info().get(k, getattr(e, k, None))) for e in cur]).value_counts()
            print(f"[select]   NOTHING LEFT. Values of '{k}' among the {len(cur)} remaining test(s): "
                  + ", ".join(f"{a} ({n})" for a, n in vals.items()))
        cur = nxt
        if len(cur) == 0:
            break
    return cur


def plot_section(exps, ref="no_ref", filters=None, require_ref=True, include_reference_tests=False,
                 max_figs=12, save=True, out_folder=None, pre_s=3.0, post_s=2.0):
    """Select tests with `filters`, find their reference (`ref`), draw one figure per test.
    ref                      no_ref | empty_ref | ninja_ref | saltwater_ref
    filters                  dict field -> value | [values] | between(lo, hi) | None (= not used). Any field of exps.table().
    require_ref              True: with ref != no_ref, only tests that HAVE that reference are drawn
    include_reference_tests  True (with ref='no_ref' only): also draw the tests of the reference days themselves"""
    if ref not in REF_MODES:
        raise ValueError(f"ref must be one of {REF_MODES}")
    sel = select_tests(exps, filters)
    if ref == "no_ref" and include_reference_tests:
        items = [(e, None) for e in sel]
    else:
        n_refday = sum(1 for e in sel if e.is_reference)
        if n_refday:
            print(f"[plot] {n_refday} selected test(s) belong to a reference day and are not drawn "
                  "(use ref='no_ref', include_reference_tests=True to draw them)")
        items = [(e, None if ref == "no_ref" else exps.find_reference(e, ref)) for e in sel if not e.is_reference]
    if ref != "no_ref":
        pool = exps.ref_pool if exps.ref_pool is not None else exps
        print(f"[plot] {ref}: {sum(r is not None for _, r in items)} / {len(items)} selected test(s) have a reference")
        for e, r in items:
            if r is None:
                print(f"   no match: {e.date} {e.time} {e.test_key} -> {closest_reference(e, pool, ref)}")
        if require_ref:
            items = [(e, r) for e, r in items if r is not None]
    if len(items) > max_figs:
        print(f"[plot] {len(items)} figure(s) requested, only the first {max_figs} are drawn (max_figs)")
        items = items[:max_figs]
    print(f"[plot] drawing {len(items)} figure(s)")
    out = Path(out_folder) if out_folder else ROOT_FOLDER.parent / "figs_ref"
    for e, r in items:
        fig = plot_pair(e, ref, r, pre_s, post_s)
        if fig is not None and save:
            try:
                out.mkdir(parents=True, exist_ok=True)
                path = out / f"{ref}_{e.date}_{e.time}.png"
                k = 2
                while path.exists():                              # never overwrite
                    path = out / f"{ref}_{e.date}_{e.time}_{k}.png"
                    k += 1
                fig.savefig(path, dpi=200)
                print(f"  saved {path}")
            except OSError as exc:
                print(f"  could not save the figure ({exc})")
    plt.show()


# %% LOAD the expedition (+ console log + reference report) ======================
# Run this cell first. It defines E (the ExperimentSet) used by all the plotting sections below.
LOAD_DAYS = None                  # None = every day, or e.g. ["2026-11-15"] (the reference days are always loaded for the pool)
TABLE_CSV = None                  # e.g. "tests.csv" to also write the one-row-per-test table
SAVE_CONSOLE_LOG = True           # everything printed is also saved as console_log_<date-time>.txt
CONSOLE_LOG_FOLDER = ROOT_FOLDER  # the Phonefleet_data folder
REPORT_PER_TEST = True            # True: 2 lines per test (which reference was found); False: only counts + tests without reference
REPORT_MEDIA = None               # None = all media, or e.g. ["frazil", "labgrown_frazil"]
CHECK_DECODE = True               # try to decode the usb of EVERY test and list the ones that cannot be decoded
DECODE_REPORT_CSV = None          # e.g. "usb_not_decodable.csv" to also save that list

if SAVE_CONSOLE_LOG:
    start_console_log(CONSOLE_LOG_FOLDER)
E = load_expedition(ROOT_FOLDER, LOAD_DAYS, SUMMARY_FOLDER)
if TABLE_CSV:
    E.table().to_csv(TABLE_CSV, index=False)
    print(f"Table written: {Path(TABLE_CSV).resolve()}")
reference_report(E, per_test=REPORT_PER_TEST, media=REPORT_MEDIA)
if CHECK_DECODE:
    NOT_DECODABLE = decodability_report(E, csv_path=DECODE_REPORT_CSV)


# %% PLOT 1: frazil sines against the empty-motor reference ======================
# Copy this whole cell to make another plotting section. None = filter not used.
# Fields you can filter on: medium, station, test_type, date, time, direction, freq_hz, mean_rpm, amp_rpm,
#   amp_start_rpm, amp_end_rpm, hold_rpm, speed_min_rpm, speed_max_rpm, n_levels, plateau_s, ramp_shape, cycles, duration_s,
#   blade, is_reference, ref_type, test_key ... (see E.table().columns).
# Values: one value, a list [a, b], or between(lo, hi) for numbers (e.g. freq_hz=between(0.1, 1)).
# test_type seen so far: sine, sweep, sweep_stepped, ramp, steps (long constant steps), hold (only 1, not frazil), rest.
# --- temporarily disabled (debug section below). Remove the leading '# ' of the lines below to re-enable. ---
# plot_section(
#     E,
#     ref="empty_ref",                  # no_ref | empty_ref | ninja_ref | saltwater_ref
#     filters=dict(
#         medium=["frazil"],            # air-empty | saltwater | frazil | labgrown_frazil | unknown
#         station=None,                 # e.g. ["301", "400"]
#         test_type=["sine"],
#         date=None,                    # e.g. ["2026-11-15"]
#         direction=None,               # up | down | updown | flat (sweeps / ramps / steps)
#         freq_hz=None,                 # e.g. between(0.1, 1)
#         amp_rpm=None,                 # e.g. [50, 500]
#         mean_rpm=None,
#     ),
#     require_ref=True,                 # only draw the tests that have the chosen reference
#     include_reference_tests=False,    # only used with ref="no_ref"
#     max_figs=12,
#     save=True,
#     out_folder=None,                  # None = <chipiron folder>/figs_ref
#     pre_s=3.0,                        # seconds shown before the motor command
#     post_s=2.0,                       # seconds shown after the end of the test
# )



# %% DEBUG: visual check of the test <-> reference matching ======================
# For each reference medium and each reference TEST, one figure with two subplots:
#   left  = all the frazil tests (every station / day) that were matched to this reference test
#   right = the reference curve itself (same axes limits, so the two sides can be compared directly)
# The matching is the one of the report above (closest input first). A "~" in the legend = near match (see ~NEAR in the report).
# Needs the Helpers and LOAD cells to have been run. Everything printed is also in the console log.
DBG_MEDIUM = "frazil"                                      # medium of the tests drawn on the left (None = all non-reference tests)
DBG_REF_TYPES = ["empty_ref", "ninja_ref", "saltwater_ref"]
DBG_TEST_TYPES = None                                      # None = all, or e.g. ["sine"], ["sine", "sweep"]
DBG_QUANTITY = "speed"                                     # "speed" or "current"
DBG_MAX_FIGS = 200                                         # safety: at most this many figures (all groups are still listed in the console)
DBG_SAVE = True                                            # True: save the figures (never overwrites) and CLOSE each one after saving
DBG_SHOW = False                                           # True: keep the figures open and block in plt.show() (use with few figures only)
DBG_OUT_FOLDER = None                                      # None = <chipiron folder>/figs_match_debug
DBG_PRE_S, DBG_POST_S = 3.0, 2.0                           # seconds shown before the motor command / after the end of the test


def debug_match_plots(exps, medium="frazil", ref_types=("empty_ref", "ninja_ref", "saltwater_ref"), test_types=None,
                      quantity="speed", max_figs=200, save=True, out_folder=None, pre_s=3.0, post_s=2.0, show=False):
    if quantity not in ("speed", "current"):
        raise ValueError("quantity must be 'speed' or 'current'")
    plt.rcParams["figure.max_open_warning"] = 0
    groups = {}                                            # (ref_type, ref date, ref time) -> (ref, [tests])
    for e in exps:
        if e.is_reference or (medium and e.medium != medium) or (test_types and e.test_type not in test_types):
            continue
        for rt in ref_types:
            m = e.ref_available.get(rt) or []
            if m:
                groups.setdefault((rt, m[0].date, m[0].time), (m[0], []))[1].append(e)
    order = {rt: i for i, rt in enumerate(ref_types)}
    items = sorted(groups.items(), key=lambda kv: (order[kv[0][0]], kv[0][1], kv[0][2]))
    print("\n" + "=" * 100)
    print(f"DEBUG MATCH PLOTS: {sum(len(v[1]) for _, v in items)} (test, reference) pair(s) in {len(items)} group(s) "
          f"[medium={medium}, test_types={test_types or 'all'}]")
    print("=" * 100)
    for (rt, d, t), (ref, tests) in items:
        print(f"{rt:<14} ref {d} {t}  {ref.test_key}   <-  {len(tests)} test(s)"
              + ("" if usb_usable(ref)[0] else "   [REFERENCE USB NOT DECODABLE: " + usb_usable(ref)[1] + "]"))
        for e in tests:
            nf = near_fields(e, ref)
            print(f"      {e.date} {e.time} st:{e.station or '-':<4} {e.test_key}" + (f"   ~NEAR[{'; '.join(nf)}]" if nf else "")
                  + ("" if usb_usable(e)[0] else "   [USB NOT DECODABLE: not drawn]"))
    if len(items) > max_figs:
        print(f"[debug] {len(items)} groups, only the first {max_figs} are drawn (DBG_MAX_FIGS)")
    out = Path(out_folder) if out_folder else ROOT_FOLDER.parent / "figs_match_debug"
    ls_cycle = ["-", "--", ":", "-."]
    keep = show or not save                                # nothing saved and nothing shown would be pointless
    for (rt, d, t), (ref, tests) in items[:max_figs]:
        r = load_run(ref)
        runs = [(e, load_run(e)) for e in tests]
        runs = [(e, a) for e, a in runs if a is not None]
        fig, (axl, axr) = plt.subplots(1, 2, figsize=(14, 5), sharex=True, sharey=True, layout="constrained")
        durs = [x.duration_s for x in tests + [ref] if x.duration_s]
        t_end = (max(durs) if durs else max([a["t"].max() for _, a in runs] + [1.0])) + post_s
        stations = sorted({str(e.station or "-") for e, _ in runs})
        cmap = plt.get_cmap("tab10")
        seen = {}
        ylab = None
        for e, a in runs:
            st = str(e.station or "-")
            k = seen[st] = seen.get(st, -1) + 1
            w = (a["t"] >= -pre_s) & (a["t"] <= t_end)
            near = "~" if near_fields(e, ref) else ""
            axl.plot(a["t"][w], a[quantity][w], lw=1, color=cmap(stations.index(st) % 10), ls=ls_cycle[k % 4],
                     label=f"{near}{e.date[5:]} {e.time} st{st}")
            ylab = a["ylab"] if quantity == "speed" else a["ilab"]
        if r is not None:
            w = (r["t"] >= -pre_s) & (r["t"] <= t_end)
            axr.plot(r["t"][w], r[quantity][w], lw=1, color="k", label=f"{ref.date[5:]} {ref.time}")
            ylab = ylab or (r["ylab"] if quantity == "speed" else r["ilab"])
            axr.legend(fontsize=8)
        else:
            axr.text(0.5, 0.5, "reference could not be decoded", transform=axr.transAxes, ha="center")
        if runs:
            axl.legend(fontsize=7, ncol=2 if len(runs) > 12 else 1)
        axl.set_title(f"{medium or 'tests'}: {len(runs)} test(s) matched", fontsize=9)
        axr.set_title(f"{rt}: {ref.date} {ref.time} (st {ref.station or '-'}, {ref.medium})", fontsize=9)
        for ax in (axl, axr):
            ax.set(xlabel=r"$t$ (s)")
            ax.grid(alpha=0.3)
        axl.set_ylabel(ylab or quantity)
        fig.suptitle(ref.test_key, fontsize=10)
        if save:
            try:
                out.mkdir(parents=True, exist_ok=True)
                path = out / f"match_{rt}_{ref.date}_{ref.time}_{quantity}.png"
                k = 2
                while path.exists():
                    path = out / f"match_{rt}_{ref.date}_{ref.time}_{quantity}_{k}.png"
                    k += 1
                fig.savefig(path, dpi=150)
                print(f"  saved {path}")
            except OSError as exc:
                print(f"  could not save the figure ({exc})")
        if not keep:
            plt.close(fig)                                 # free the window / memory once it is saved
    if keep:
        plt.show()


debug_match_plots(E, medium=DBG_MEDIUM, ref_types=DBG_REF_TYPES, test_types=DBG_TEST_TYPES, quantity=DBG_QUANTITY,
                  max_figs=DBG_MAX_FIGS, save=DBG_SAVE, out_folder=DBG_OUT_FOLDER, pre_s=DBG_PRE_S, post_s=DBG_POST_S,
                  show=DBG_SHOW)


# %% Stop the console log (optional) =============================================
# stop_console_log()