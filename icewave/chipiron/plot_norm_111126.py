#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_norm.py -- sine-speed runs OR ramp tests in frazil, alone or compared with the SAME kind of run in air.

TWO TESTS  (choose with  --test sine | ramp_test ,  or TEST below ,  or answer the prompt)
TWO MODES  (--mode no_ref | ref ; asked if not given) -- valid for BOTH tests now

  sine       frazil run (+ air run in ref mode). Figs 1-3, and in ref mode fig 4 (normalisation by the air).

  ramp_test  ONE frazil run (RAMP_RUN) made with chipiron_run_v4.py (holds / ramps), + in ref mode ONE air run
             (AIR_RAMP_RUN) made the same way. The commands of each run (input_*.csv, notes "ramp ..." /
             "hold ...") tell when the ramps happen.
             * ONLY THE FIRST ramp and the first ramp in the opposite direction are kept, the plateau between
               them is cut at its middle (same rule for both runs). Up / down in two colours.
                Fig 1  speed   vs time   (+ commanded speed, dashed, if SHOW_COMMAND)
                Fig 2  current vs time
                Fig 3  current vs speed  (up / down: hysteresis loop; ref mode: air in grey behind)
                Fig 4  (ref mode) frazil current MINUS air current, vs speed, up / down in two colours

FIG 4 OF THE RAMP TEST: RAW TIME-MATCHED AIR SUBTRACTION
  For every frazil measurement, the air current is evaluated at the SAME elapsed time since the
  start of the first ramp, then Delta I = I_frazil - I_air.
  No spike filtering, no speed matching, no binning, and no median/averaging are applied.
  The air curve is linearly interpolated in time only if its sample timestamps do not fall exactly
  on the frazil timestamps. Points outside the air time range are omitted.
  For a physically meaningful comparison, the two runs should use the same ramp program and timing.

AXIS LABELS   quantity, then its unit just after it, then what it is proportional to
                 speed   : omega (rpm) ~ gamma_dot     current : I (A) ~ tau

FIG 4 NORMALISATION METHODS (sine, ref mode only)
  ratio        I / I_air
  diff         I - I_air                 (same phase)
  diff_scaled  (I - I_air) / half-range of I_air
  diff_speed   I - I_air at the same measured speed and same branch (accelerating / decelerating)
  regression   I - (a I_air + b)
  PHASE_AVERAGE_FIG4 = True : phase-bin both runs (median), N_PHASE_BINS points per cycle.
  PHASE_AVERAGE_FIG4 = False: keep every frazil measurement; the air current is interpolated at the same
                              phase inside the same cycle (diff_speed: at the same speed and branch).
  Fig 4 always uses the physical current (A), so it needs the .edc.

DEPENDENCIES   numpy, pandas, matplotlib, chipiron_decode_v4.py, escon_edc_v3.py (in CODE_FOLDER)
RUN            py plot_norm.py                         (asks the test, then the mode)
               py plot_norm.py --test sine --mode ref
               py plot_norm.py --test ramp_test --mode ref [--force-time]
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# =============================== USER CONFIGURATION ===============================
try:
    CODE_FOLDER = Path(__file__).resolve().parent
except NameError:
    CODE_FOLDER = Path(".").resolve()

DATA_ROOT = Path(r"D:/CASCADES26/chipiron/Phonefleet_data")            # contains <date>/<time>/
EDC_FOLDER = Path(r"D:/CASCADES26/chipiron/edc_ESCON_EC_50_5")
# searched recursively for the .edc names
SAVE_FOLDER = Path(r"D:/CASCADES26/chipiron/figs_frazil_vs_air_ramptest_1111_1311")

TEST = None                # "sine" | "ramp_test" | None (= ask; the command line --test overrides this)

# ---- sine test
FRAZIL = ("2026-11-11", "16-46-03")      # (date folder, time folder)  -> always used (sine)
AIR = ("2026-11-13", "11-15-27")         # reference (empty motor + vane) -> used ONLY in ref mode (sine)
MODE = "no_ref"            # "no_ref" (frazil alone, NO air correction) | "ref" (air overlay + subtraction) | None (= ask).
                       # The command line --mode overrides this. Works for BOTH sine AND ramp_test.

# ---- ramp test
RAMP_RUN = ("2026-11-11", "14-19-28")    # (date folder, time folder) of the frazil ramp run
AIR_RAMP_RUN = ("2026-11-13", "10-08-28")  # air reference ramp run (used ONLY in ref mode)
# Ramp subtraction is raw and time-matched only: no filtering, binning, or averaging.
RAMP_T_TOL_S = 1.0         # time method / input check: commands sent at times differing by more than this are reported
RAMP_FORCE_TIME = False    # True, or pass --force-time: draw raw time subtraction even if command timing differs
RAMP_PRE_S = 0.0           # grey part shown before the first ramp (s); 0 = none
RAMP_POST_S = 0.0          # grey part shown after the end of the second ramp (s); 0 = none
SHOW_COMMAND = True        # fig 1: dashed line = commanded speed (needs the .edc / REAL_UNITS)
RAMP_DT_MAX_S = 1.5        # only used if the input log has NO notes (old runs): commands closer than this = ramp
COLOR_UP = "tab:red"       # first / increasing ramp (+ first half of the plateau)
COLOR_DOWN = "tab:blue"    # decreasing ramp (+ second half of the plateau)
COLOR_OUT = "0.65"         # before / after, and air reference in fig 3

REAL_UNITS = True          # True: rpm and A (from the .edc). False: Arduino volts A1 / A2
SAVE_PLOTS = True
FIG_FORMATS = ("png",)     # e.g. ("png", "pdf")
FIG_DPI = 200
OVERWRITE = False          # False: never overwrite, add _2, _3...

ADC_VREF = 5.06
WIRING = {"A1": "AO1", "A2": "AO2"}   # Arduino pin -> ESCON analog output (AO1 = speed, AO2 = current)
STRICT_INPUT_CHECK = True  # sine ref mode: True = stop if the two runs did not get the same input; False = warn

SKIP_CYCLES = 0            # first cycles ignored in figs 3-4 (start-up transient)
TRIM_TO_SAME_DURATION = True   # ref mode: trim the longer run to the same number of cycles
SINE_ONSET_SHIFT_S = 0.0   # fine shift of t = 0 (s) if the command / data clocks look offset (sine AND ramp_test)
TIME_PLOT_PRE_S = 3.0      # time shown before the sine command (s) in figs 1-2
TIME_PLOT_POST_S = 2.0     # time shown after the end of the sine (s)
SHARE_Y = True             # ref mode: same y range for frazil and air subplots
SHOW_TITLES = False        # figure titles (the slide already has one). Panel names Frazil / Air are always shown in ref mode

# Fig 3 (sine)
COLOR_BY = "phase"         # "phase" (cyclic, follows the loop) or "time" (absolute time, shows drift)
CMAP_PHASE = "hsv"
CMAP_TIME = "viridis"
POINT_SIZE = 18

# Fig 4 (sine, ref mode)
PHASE_AVERAGE_FIG4 = False # True = phase average; False = raw points
N_PHASE_BINS = 50
N_SPEED_BINS = 25
RATIO_MIN_FRAC = 0.10      # ratio: mask points where |I_air| < this fraction of max |I_air|
# None = ask | "diff" | ["diff", "diff_speed"] | "all"
NORM_METHOD = None

# ---- fonts: BIG so that figures are readable on a slide
FONT_LABEL = 34            # axis labels
FONT_TICK = 30             # tick labels
FONT_TITLE = 32            # panel / figure titles
FONT_LEGEND = 28
LINE_WIDTH = 2.2           # curves in the time plots
TEXT_FONT = ["Arial", "Helvetica", "DejaVu Sans"]   # first one found is used (text AND math symbols)

FIGSIZE_2ROWS = (13, 11)
FIGSIZE_2COLS = (18, 8.5)
FIGSIZE_1 = (11, 8)

COLOR_FRAZIL = "tab:blue"
COLOR_AIR = "tab:orange"
# ==================================================================================


sys.path.insert(0, str(CODE_FOLDER))
from chipiron_decode_v4 import read_usb_rows, decode_usb  # noqa: E402

OSC_RE = re.compile(r"o([\d.]+),([\d.]+),([\d.]+)")
CYCLES_RE = re.compile(r"([\d.]+) cycles")

LABEL = {"frazil": "Frazil", "air": "Air"}      # panel names (ref mode)
SOURCES = {"frazil": FRAZIL, "air": AIR}


# ------------------------------------------------------------------ test / mode
def get_test():
    """--test on the command line > TEST above > prompt (Enter = sine)."""
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("--test", choices=["sine", "ramp_test"], default=None)
    args, _ = ap.parse_known_args()
    test = args.test or TEST
    if test is None:
        print("Test:\n  1) sine      - sine-speed run (alone or against an air reference)"
              "\n  2) ramp_test - first ramp of a ramp run, up / down in two colours (alone or against an air ramp)")
        test = {"2": "ramp_test", "ramp_test": "ramp_test"}.get(input("Choice [1]: ").strip().lower(), "sine")
    if test not in ("sine", "ramp_test"):
        raise SystemExit(f"ERROR: unknown test '{test}' (use sine or ramp_test).")
    print(f"\n>>> Test: {test}")
    return test


def get_mode():
    """--mode on the command line > MODE above > prompt (Enter = no_ref)."""
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("--mode", choices=["no_ref", "ref"], default=None)
    args, _ = ap.parse_known_args()       # tolerant: also works inside Spyder / VS Code interactive
    mode = args.mode or MODE
    if mode is None:
        print("Mode:\n  1) no_ref - frazil run only"
              "\n  2) ref    - frazil + air reference (subplots / subtraction)")
        mode = {"2": "ref", "ref": "ref"}.get(input("Choice [1]: ").strip().lower(), "no_ref")
    if mode not in ("no_ref", "ref"):
        raise SystemExit(f"ERROR: unknown mode '{mode}' (use no_ref or ref).")
    print(f"\n>>> Mode: {mode}")
    return mode


def get_ramp_methods():
    """Ramp reference subtraction is deliberately time-matched only."""
    return ["time"]


def force_time_requested():
    """Whether to plot time-matched subtraction despite a command/timing mismatch."""
    return bool(RAMP_FORCE_TIME or "--force-time" in sys.argv)


# ------------------------------------------------------------------ loading
def part_key(p):
    m = re.search(r"Usb-(\d+)", p.name, re.I)
    return (int(m.group(1)) if m else 10 ** 9, p.name)


def _load_common(date, time_, kind):
    """Folder content shared by both tests: usb parts, input log, edc name, sampling period."""
    folder = DATA_ROOT / date / time_
    if not folder.is_dir():
        raise FileNotFoundError(f"[{kind}] folder not found: {folder}")
    usb = sorted((p for p in folder.glob("*.csv")
                  if "usb" in p.name.lower() and not p.name.startswith("input_")), key=part_key)
    if not usb:
        raise FileNotFoundError(f"[{kind}] no Usb csv in {folder}")
    logs = sorted(folder.glob("input_*.csv"))
    if not logs:
        raise FileNotFoundError(f"[{kind}] no input_*.csv in {folder}")
    log = pd.read_csv(logs[0], dtype=str).fillna("")
    edc_txt = (folder / "edc_used.txt").read_text(encoding="utf-8").splitlines()
    full = next((l.split(":", 1)[1].strip() for l in edc_txt if l.lower().startswith("# full path")), "")
    period = float((folder / "sampling_period_ms.txt").read_text(encoding="utf-8").split()[0])
    return dict(kind=kind, folder=folder, usb=usb, log_file=logs[0], log=log, edc_name=edc_txt[0].strip(),
                edc_full=full, period=period)


def load_meta(date, time_, kind):
    """Sine run: common content + sine parameters taken from the input log."""
    run = _load_common(date, time_, kind)
    log = run["log"]
    cmd = log["command"]
    i_start = log.index[cmd == "[acquisition /start]"]
    i_osc = log.index[cmd.str.match(r"^o[\d.]+,") & (log["source"] == "http")]
    i_end = log.index[(cmd == "x") & (log["source"] == "http")]
    if len(i_start) != 1 or len(i_osc) != 1 or len(i_end) != 1:
        raise ValueError(f"[{kind}] {run['log_file'].name}: expected one /start, one 'o...' and one 'x' command "
                         f"(found {len(i_start)}, {len(i_osc)}, {len(i_end)}).")
    t = log["t_unix_s"].astype(float)
    mean_d, amp_d, freq = (float(x) for x in OSC_RE.fullmatch(cmd[i_osc[0]]).groups())
    note = log.loc[i_osc[0], "note"]
    m = CYCLES_RE.search(note)
    t_on = t[i_osc[0]] - t[i_start[0]]                       # sine command, seconds after /start
    t_x = t[i_end[0]] - t[i_start[0]]
    n_cyc = float(m.group(1)) if m else (t_x - t_on) * freq
    print(f"[{kind}] {date} {time_}: {len(run['usb'])} Usb file(s), sine {freq:g} Hz, {n_cyc:g} cycles, "
          f"edc {run['edc_name']}, period {run['period']:g} ms")
    run.update(mean_duty=mean_d, amp_duty=amp_d, freq=freq, note=note, n_cycles=n_cyc,
               t_on=t_on + SINE_ONSET_SHIFT_S)
    return run


# ---- ramp test: commands -> ramps
def ramp_commands(log, kind):
    """DataFrame of the duty commands c<duty> sent by the run script after /start:
    t (s after /start), duty (%), kind ('hold' / 'ramp' / 'rest'), seg (label)."""
    cmd, src, note = log["command"], log["source"], log["note"]
    t = log["t_unix_s"].astype(float)
    i_start = log.index[cmd == "[acquisition /start]"]
    if len(i_start) != 1:
        raise ValueError(f"[{kind}] expected exactly one '[acquisition /start]' in the input log "
                         f"(found {len(i_start)}).")
    t0 = t[i_start[0]]
    rows = []
    for i in log.index[log.index > i_start[0]]:
        if src[i] != "http" or not re.fullmatch(r"c[\d.]+", cmd[i]):
            continue
        parts = note[i].split(" ", 1)
        rows.append(dict(t=t[i] - t0, duty=float(cmd[i][1:]), kind=parts[0],
                         seg=parts[1] if len(parts) > 1 else ""))
    c = pd.DataFrame(rows)
    if c.empty:
        raise ValueError(f"[{kind}] no c<duty> command found after /start in the input log.")
    if not c["kind"].isin(["hold", "ramp"]).any():          # old run: no notes -> guess from the spacing
        dt_next = np.r_[np.diff(c["t"].to_numpy()), np.inf]
        c["kind"] = np.where(dt_next <= RAMP_DT_MAX_S, "ramp", "hold")
        print(f"[{kind}] WARNING: the input log has no 'ramp' / 'hold' notes: commands closer than "
              f"{RAMP_DT_MAX_S:g} s are taken as ramp steps.")
    return c


def find_ramp_groups(c):
    """Ramps = consecutive 'ramp' commands whose duty keeps going in the same direction.
    Returns a list of dicts(first, last, sign) (indices in c; sign +1 = duty increasing, -1 = decreasing)."""
    kinds, duty = c["kind"].to_numpy(), c["duty"].to_numpy()
    groups, cur = [], None
    for i in range(len(c)):
        if kinds[i] != "ramp":
            cur = None
            continue
        prev = duty[i - 1] if i > 0 else duty[i]
        s = int(np.sign(duty[i] - prev))
        if cur is None or (s != 0 and cur["sign"] not in (0, s)):
            cur = dict(first=i, last=i, sign=s)
            groups.append(cur)
        else:
            cur["last"] = i
            if cur["sign"] == 0:
                cur["sign"] = s
    return [g for g in groups if g["sign"] != 0]


def _t_end(c, g):
    """End of a ramp = time of the command that follows it."""
    k = g["last"] + 1
    return float(c["t"].iloc[k]) if k < len(c) else float(c["t"].iloc[g["last"]]) + 1.0


def pick_first_ramp_pair(c, groups, kind):
    """First ramp (+ following ramps of the same direction) and the first ramp of the opposite direction.
    Times in s after /start. The plateau between them is cut at its middle."""
    g0 = groups[0]
    k = next((i for i, g in enumerate(groups) if g["sign"] == -g0["sign"]), None)
    name = lambda s: "up" if s > 0 else "down"
    print(f"[{kind}] {len(groups)} ramp(s) detected:")
    for i, g in enumerate(groups):
        print(f"    #{i + 1}: {name(g['sign']):<4} t = {c['t'].iloc[g['first']]:.1f} -> {_t_end(c, g):.1f} s, "
              f"duty {c['duty'].iloc[max(g['first'] - 1, 0)]:g} -> {c['duty'].iloc[g['last']]:g} %")
    if k is None:
        print(f"[{kind}] WARNING: no ramp in the opposite direction found: only the first ramp is plotted.")
        a1, b0, b1, dir2 = _t_end(c, g0), None, None, None
        t_split = a1
        t_end = a1
    else:
        a1 = _t_end(c, groups[k - 1])
        b0, b1, dir2 = float(c["t"].iloc[groups[k]["first"]]), _t_end(c, groups[k]), name(groups[k]["sign"])
        t_split, t_end = 0.5 * (a1 + b0), b1
        print(f"[{kind}] kept: first ramp ({name(g0['sign'])}) then the first {dir2} ramp; plateau between "
              f"them = {b0 - a1:.1f} s, cut at its middle (t = {t_split:.1f} s).")
    return dict(t_start=float(c["t"].iloc[g0["first"]]), t_split=t_split, t_end=t_end,
                dir1=name(g0["sign"]), dir2=dir2, n_ramps=len(groups))


def load_meta_ramp(date, time_, kind):
    run = _load_common(date, time_, kind)
    c = ramp_commands(run["log"], kind)
    groups = find_ramp_groups(c)
    if not groups:
        raise ValueError(f"[{kind}] no ramp found among the {len(c)} commands of {run['log_file'].name} "
                         f"(kinds: {sorted(set(c['kind']))}). Is this a ramp run (not a staircase / sine)?")
    r = pick_first_ramp_pair(c, groups, kind)
    t_on = r["t_start"] + SINE_ONSET_SHIFT_S                 # t = 0 : start of the first ramp
    run["ramp"] = dict(dir1=r["dir1"], dir2=r["dir2"], n_ramps=r["n_ramps"],
                       t_split=r["t_split"] - t_on, t_end=r["t_end"] - t_on)
    run["cmds"] = c.assign(t_rel=c["t"] - t_on)
    run["t_on"] = t_on
    print(f"[{kind}] {date} {time_}: {len(run['usb'])} Usb file(s), edc {run['edc_name']}, "
          f"period {run['period']:g} ms")
    return run


def check_same_input(a, b):
    """Sine ref mode: same Arduino commands, same .edc, same sampling period. Only the number of cycles may differ."""
    problems = []
    ca = [c for c in a["log"]["command"] if not c.startswith("[")]
    cb = [c for c in b["log"]["command"] if not c.startswith("[")]
    if ca != cb:
        problems.append(f"Arduino commands differ:\n    {a['kind']}: {ca}\n    {b['kind']}: {cb}")
    if a["edc_name"] != b["edc_name"]:
        problems.append(f".edc differs: {a['edc_name']} vs {b['edc_name']}")
    if a["period"] != b["period"]:
        problems.append(f"sampling period differs: {a['period']:g} ms vs {b['period']:g} ms")
    strip = lambda s: CYCLES_RE.sub("N cycles", s)
    if strip(a["note"]) != strip(b["note"]):
        problems.append(f"command notes differ: '{a['note']}' vs '{b['note']}'")
    if problems:
        msg = "The two runs did NOT get the same input:\n  - " + "\n  - ".join(problems)
        if STRICT_INPUT_CHECK:
            raise SystemExit("ERROR: " + msg + "\n  (wrong AIR folder? or run with --mode no_ref)")
        print("WARNING:", msg)
    else:
        print(f"Input check OK: same commands {ca}, same .edc ({a['edc_name']}), same period ({a['period']:g} ms).")
    if a["n_cycles"] != b["n_cycles"]:
        print(f"NOTE: different duration: {a['kind']} {a['n_cycles']:g} cycles, {b['kind']} {b['n_cycles']:g} cycles"
              + (" -> trimmed to the shorter one." if TRIM_TO_SAME_DURATION else "."))


def check_same_input_ramp(a, b):
    """Report whether the two ramp runs used the same commands and timing.
    A mismatch is only a warning: raw time-matched subtraction is still plotted."""
    if a["edc_name"] != b["edc_name"]:
        raise SystemExit(f"ERROR: .edc differs: {a['edc_name']} vs {b['edc_name']} (mapping / ESCON ramp not the same).")
    if a["period"] != b["period"]:
        raise SystemExit(f"ERROR: sampling period differs: {a['period']:g} ms vs {b['period']:g} ms.")
    print(f"Input check: same .edc ({a['edc_name']}), same period ({a['period']:g} ms).")
    ca, cb = a["cmds"], b["cmds"]
    same = len(ca) == len(cb) and np.allclose(ca["duty"].to_numpy(), cb["duty"].to_numpy())
    if same:
        dt = np.abs(ca["t_rel"].to_numpy() - cb["t_rel"].to_numpy())
        if dt.max() > RAMP_T_TOL_S:
            print(f"WARNING: same duty commands but their timing differs by up to {dt.max():.1f} s "
                  f"(> {RAMP_T_TOL_S:g} s). Raw time-matched subtraction will still be plotted.")
            same = False
        else:
            print(f"Same duty commands in both runs ({len(ca)} commands, timing within {dt.max():.2f} s).")
    else:
        print("\n!! The two ramp runs did NOT receive the same program:")
        for run in (a, b):
            c = run["cmds"]
            rp = c[c["kind"] == "ramp"]
            ramps = sorted(set(rp["seg"]))
            print(f"   {run['kind']}: {len(c)} commands, duty {c['duty'].min():g}..{c['duty'].max():g} %, "
                  f"ramps: {ramps}")
        print("   -> raw subtraction at identical elapsed times will still be plotted, but different programs "
              "can make the comparison less physically meaningful.")
    if a["ramp"]["dir1"] != b["ramp"]["dir1"] or a["ramp"]["dir2"] != b["ramp"]["dir2"]:
        raise SystemExit(f"ERROR: ramp directions differ ({a['ramp']['dir1']}/{a['ramp']['dir2']} vs "
                         f"{b['ramp']['dir1']}/{b['ramp']['dir2']}): up and down cannot be matched.")
    return same


def find_edc(run):
    """Look for the .edc named in edc_used.txt in EDC_FOLDER (recursive), else at its recorded full path."""
    hits = list(EDC_FOLDER.rglob(run["edc_name"])) if EDC_FOLDER.is_dir() else []
    if hits:
        return hits[0]
    p = Path(run["edc_full"]) if run["edc_full"] else None
    if p is not None and p.is_file():
        return p
    raise FileNotFoundError(f"[{run['kind']}] .edc '{run['edc_name']}' not found in {EDC_FOLDER} "
                            f"(nor at {run['edc_full'] or 'no recorded path'})")


def load_run_edcs(runs, ref=False):
    """.edc of each run (needed for rpm / A, for the command overlay, and for fig 4)."""
    edcs = {}
    for k, r in runs.items():
        try:
            path = find_edc(r)
            from escon_edc_v3 import load_edc
            edcs[k] = load_edc(path)
            print(f"[{k}] .edc: {path}")
        except Exception as exc:
            if REAL_UNITS:
                raise SystemExit(f"ERROR: REAL_UNITS is True but the .edc cannot be used: {exc}")
            print(f"WARNING: {exc}  (volts plotted anyway{'; fig 4 will be skipped' if ref else ''})")
            edcs[k] = None
    return edcs


def decode_run(run, edc):
    """Decode the Usb file(s), convert with the .edc, add speed / current columns (+ *_phys), t,
    and phase / cycle (sine) or 'phase' 0 / 1 / 2 (ramp_test)."""
    dfs = [read_usb_rows(p, verbose=False) for p in run["usb"]]
    raw = pd.concat([d for d in dfs if len(d)], ignore_index=True)
    if raw.empty:
        raise ValueError(f"[{run['kind']}] no valid USB rows found.")
    raw["t_s"] = (raw["t_us"] - raw["t_us"].iloc[0]) / 1e6
    d = decode_usb(raw, vref=ADC_VREF, period_ms=run["period"], verbose=False)
    d = d[d["kind"] == "full"].dropna(subset=["A1_adc", "A2_adc"]).reset_index(drop=True)
    print(f"[{run['kind']}] {len(d)} frames, framing '{d.attrs.get('framing')}', "
          f"time axis '{d.attrs.get('time_mode')}', {d['t_s'].iloc[-1]:.1f} s")
    d["speed_V"], d["current_V"] = d["A1_V"], d["A2_V"]
    if edc is not None:
        edc.apply_to_frame(d, wiring=WIRING)
        sc = [c for c in d.columns if c.startswith("A1_") and c.endswith("_rpm")]
        cc = [c for c in d.columns if c.startswith("A2_") and c.endswith("_A")]
        if not sc or not cc:
            raise ValueError(f"[{run['kind']}] with WIRING {WIRING} the .edc gives no speed on A1 / current on A2. "
                             "Check WIRING.")
        d["speed_phys"], d["current_phys"] = d[sc[0]], d[cc[0]]
    d["t"] = d["t_s"] - run["t_on"]                          # 0 = sine command / start of the first ramp
    if run.get("freq"):                                      # sine
        d["ph"] = (d["t"] * run["freq"]) % 1.0
        d["cycle"] = np.floor(d["t"] * run["freq"]).astype(int)  # cycle number since the sine onset
    if "ramp" in run:                                        # ramp test: 0 = outside, 1 = first ramp + 1/2 plateau, 2 = rest
        r = run["ramp"]
        t = d["t"]
        d["phase"] = np.select([(t >= 0) & (t <= r["t_split"]), (t > r["t_split"]) & (t <= r["t_end"])], [1, 2], 0)
    return d


# ------------------------------------------------------------------ plotting helpers
def setup_style():
    """Sans-serif text. Math text ($...$: Greek letters, dot, units) uses the SAME sans font -> no LaTeX-looking text."""
    plt.rcParams.update({
        "text.usetex": False,
        "font.family": "sans-serif", "font.sans-serif": TEXT_FONT,
        "mathtext.fontset": "custom", "mathtext.rm": "sans", "mathtext.it": "sans:italic",
        "mathtext.bf": "sans:bold", "mathtext.default": "it",
        "axes.labelsize": FONT_LABEL, "axes.titlesize": FONT_TITLE, "font.size": FONT_TICK,
        "xtick.labelsize": FONT_TICK, "ytick.labelsize": FONT_TICK, "legend.fontsize": FONT_LEGEND,
        "axes.linewidth": 1.8, "xtick.major.width": 1.8, "ytick.major.width": 1.8,
        "xtick.major.size": 9, "ytick.major.size": 9, "xtick.major.pad": 7, "ytick.major.pad": 7,
        "axes.labelpad": 10, "grid.linewidth": 1.0,
    })


def labels():
    """Axis labels: quantity, its unit right after it, then what it is proportional to.
    Speed ~ shear rate, current ~ shear stress."""
    if REAL_UNITS:
        return dict(speed=r"$\omega$ (rpm) $\propto \dot{\gamma}$", current=r"$I$ (A) $\propto \tau$")
    return dict(speed=r"$V_{A1}$ (V) $\propto \dot{\gamma}$", current=r"$V_{A2}$ (V) $\propto \tau$")


LABEL_TIME = r"$t$ (s)"
LABEL_PHASE = r"$\varphi = (t\,\mathrm{mod}\,T)/T$"


def cols():
    return ("speed_phys", "current_phys") if REAL_UNITS else ("speed_V", "current_V")


def save_fig(fig, name):
    if not SAVE_PLOTS:
        return
    SAVE_FOLDER.mkdir(parents=True, exist_ok=True)
    for ext in FIG_FORMATS:
        path = SAVE_FOLDER / f"{name}.{ext}"
        k = 2
        while path.exists() and not OVERWRITE:
            path = SAVE_FOLDER / f"{name}_{k}.{ext}"
            k += 1
        fig.savefig(path, dpi=FIG_DPI, bbox_inches="tight")
        print(f"  saved {path}")


def window(d, t_lo, t_hi):
    return d[(d["t"] >= t_lo) & (d["t"] <= t_hi)]


def tag(keys):
    return "_with_air" if len(keys) == 2 else "_frazil_only"


# ------------------------------------------------------------------ figs 1 and 2 (sine)
def fig_time(dd, keys, ycol, ylabel, title, name, dur):
    if len(keys) == 2:
        fig, axes = plt.subplots(2, 1, figsize=FIGSIZE_2ROWS, sharex=True, sharey=SHARE_Y, layout="constrained")
    else:
        fig, axes = plt.subplots(1, 1, figsize=FIGSIZE_1, layout="constrained")
    axes = np.atleast_1d(axes)
    colors = {"frazil": COLOR_FRAZIL, "air": COLOR_AIR}
    for ax, key in zip(axes, keys):
        d = window(dd[key], -TIME_PLOT_PRE_S, dur + TIME_PLOT_POST_S)
        ax.plot(d["t"], d[ycol], color=colors[key], lw=LINE_WIDTH)
        ax.set_ylabel(ylabel)
        if len(keys) == 2:
            ax.set_title(LABEL[key], loc="left")
        ax.grid(True, alpha=0.3)
    axes[-1].set_xlabel(LABEL_TIME)
    if SHOW_TITLES:
        fig.suptitle(title, fontsize=FONT_TITLE)
    save_fig(fig, name + tag(keys))
    return fig


# ------------------------------------------------------------------ fig 3 (sine)
def fig_loop(dd, keys, dur, runs):
    sx, sy = cols()
    lab = labels()
    if len(keys) == 2:
        fig, axes = plt.subplots(1, 2, figsize=FIGSIZE_2COLS, sharex=True, sharey=SHARE_Y, layout="constrained")
    else:
        fig, axes = plt.subplots(1, 1, figsize=FIGSIZE_1, layout="constrained")
    axes = np.atleast_1d(axes)
    if COLOR_BY == "time":
        cmap, clabel, vmin, vmax = CMAP_TIME, LABEL_TIME, SKIP_CYCLES / runs["frazil"]["freq"], dur
    else:
        cmap, clabel, vmin, vmax = CMAP_PHASE, LABEL_PHASE, 0.0, 1.0
    sc = None
    for ax, key in zip(axes, keys):
        d = window(dd[key], SKIP_CYCLES / runs[key]["freq"], dur)
        c = d["t"] if COLOR_BY == "time" else d["ph"]
        sc = ax.scatter(d[sx], d[sy], c=c, cmap=cmap, vmin=vmin, vmax=vmax, s=POINT_SIZE, rasterized=True)
        ax.set_xlabel(lab["speed"])
        if len(keys) == 2:
            ax.set_title(LABEL[key])
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel(lab["current"])
    cb = fig.colorbar(sc, ax=axes.ravel().tolist(), pad=0.02)
    cb.set_label(clabel, fontsize=FONT_LABEL)
    cb.ax.tick_params(labelsize=FONT_TICK)
    save_fig(fig, "fig3_current_vs_speed" + tag(keys))
    return fig


# ------------------------------------------------------------------ ramp test figures
def ramp_colors(r):
    col = {"up": COLOR_UP, "down": COLOR_DOWN}
    return col[r["dir1"]], (col[r["dir2"]] if r["dir2"] else None)


def ramp_window(d, r):
    return d[(d["t"] >= -RAMP_PRE_S) & (d["t"] <= r["t_end"] + RAMP_POST_S)]


def fig_ramp_time(d, run, ycol, ylabel, name, command=None):
    """y vs time. Grey = before / after, colour 1 = first ramp + first half of the plateau,
    colour 2 = second half of the plateau + opposite ramp. Dotted line = where the plateau is cut."""
    r = run["ramp"]
    c1, c2 = ramp_colors(r)
    w = ramp_window(d, r)
    fig, ax = plt.subplots(figsize=FIGSIZE_1, layout="constrained")
    for part in (w[(w["phase"] == 0) & (w["t"] < 0)], w[(w["phase"] == 0) & (w["t"] >= 0)]):
        if len(part):
            ax.plot(part["t"], part[ycol], color=COLOR_OUT, lw=LINE_WIDTH)
    p1, p2 = w[w["phase"] == 1], w[w["phase"] == 2]
    if len(p2):
        p1 = pd.concat([p1, p2.iloc[:1]])                    # no visual gap between the two colours
    ax.plot(p1["t"], p1[ycol], color=c1, lw=LINE_WIDTH, label=f"{r['dir1']} ramp")
    if len(p2):
        ax.plot(p2["t"], p2[ycol], color=c2, lw=LINE_WIDTH, label=f"{r['dir2']} ramp")
        ax.axvline(r["t_split"], color="k", ls=":", lw=1.5)
    #if command is not None:
        #tc, yc = command
        #m = (tc >= -RAMP_PRE_S) & (tc <= r["t_end"] + RAMP_POST_S)
        #ax.step(tc[m], yc[m], where="post", color="k", ls="--", lw=1.5, label="command")
    ax.set_xlabel(LABEL_TIME)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    save_fig(fig, name)
    return fig


def fig_ramp_loop(d, run, d_air=None):
    """current vs speed, up / down in two colours (only the two coloured phases).
    ref mode: the air ramps (same two phases) are drawn in grey behind."""
    sx, sy = cols()
    lab = labels()
    r = run["ramp"]
    c1, c2 = ramp_colors(r)
    fig, ax = plt.subplots(figsize=FIGSIZE_1, layout="constrained")
    if d_air is not None:
        a = d_air[d_air["phase"].isin([1, 2])]
        ax.scatter(a[sx], a[sy], color=COLOR_OUT, s=POINT_SIZE * 0.5, label="air", rasterized=True, zorder=1)
    p1, p2 = d[d["phase"] == 1], d[d["phase"] == 2]
    ax.scatter(p1[sx], p1[sy], color=c1, s=POINT_SIZE, label=f"{r['dir1']} ramp", rasterized=True, zorder=2)
    if len(p2):
        ax.scatter(p2[sx], p2[sy], color=c2, s=POINT_SIZE, label=f"{r['dir2']} ramp", rasterized=True, zorder=2)
    ax.set_xlabel(lab["speed"])
    ax.set_ylabel(lab["current"])
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", markerscale=2)
    save_fig(fig, "fig3_ramp_current_vs_speed" + ("_with_air" if d_air is not None else ""))
    return fig


# ------------------------------------------------------------------ ramp test: raw time-matched air subtraction (fig 4)
RAMP_NORM_INFO = {
    "time": "raw frazil current minus air current at the same elapsed time",
}


def air_current_at_time_ramp(d_air, t):
    """Air current at the same elapsed time. Interpolation only aligns sample timestamps."""
    a = d_air.dropna(subset=["current_phys", "t"]).sort_values("t")
    if len(a) < 2:
        return np.full(len(t), np.nan)

    ta = a["t"].to_numpy(dtype=float)
    ia = a["current_phys"].to_numpy(dtype=float)

    # Keep raw measurements; remove duplicate timestamps only (no averaging).
    ta, unique_idx = np.unique(ta, return_index=True)
    ia = ia[unique_idx]

    return np.interp(np.asarray(t, dtype=float), ta, ia, left=np.nan, right=np.nan)


def normalise_ramp(d_f, d_a, method="time"):
    """Keep every frazil sample and subtract air current at the same elapsed time."""
    if method != "time":
        raise ValueError("Ramp subtraction supports only method='time'.")

    w = d_f[d_f["phase"].isin([1, 2])].dropna(subset=["t", "speed_phys", "current_phys"]).copy()
    w["I_air"] = air_current_at_time_ramp(d_a, w["t"].to_numpy())
    w["value"] = w["current_phys"] - w["I_air"]

    n_bad = int(w["value"].isna().sum())
    print(f"  [time/raw] {len(w) - n_bad} / {len(w)} frazil samples subtracted, "
          f"{n_bad} outside the air time range.")
    for ph, nm in ((1, "first ramp"), (2, "second ramp")):
        values = w.loc[(w["phase"] == ph) & w["value"].notna(), "value"]
        if len(values):
            print(f"  [time/raw] phase {ph} ({nm}): Delta I mean {values.mean():+.3g} A, "
                  f"range {values.min():+.3g}..{values.max():+.3g} A")

    return w[w["value"].notna()].copy()


def fig_norm_ramp(tables, run):
    """Raw Delta I vs speed; no ramp samples are binned or averaged."""
    lab = labels()
    r = run["ramp"]
    c1, c2 = ramp_colors(r)
    n = len(tables)
    fig, axes = plt.subplots(1, n, figsize=(FIGSIZE_1[0] * n if n > 1 else FIGSIZE_1[0], FIGSIZE_1[1]),
                             layout="constrained", squeeze=False, sharey=True)
    for ax, (method, tab) in zip(axes[0], tables.items()):
        for ph, col, nm in ((1, c1, r["dir1"]), (2, c2, r["dir2"])):
            p = tab[tab["phase"] == ph]
            if len(p):
                p = p.sort_values("t")
                ax.plot(p["speed_phys"], p["value"], color=col, lw=LINE_WIDTH,
                        label=f"{nm} ramp", marker=".", markersize=4, rasterized=True)
        ax.axhline(0.0, color="k", ls="--", lw=1.5)
        ax.set_xlabel(lab["speed"])
        if n > 1:
            ax.set_title("raw time-matched subtraction")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", markerscale=2)
    axes[0][0].set_ylabel(r"$\Delta I$ (A) $\propto \Delta\tau$")
    save_fig(fig, "fig4_ramp_dI_" + "_".join(tables) + "_frazil_minus_air")
    return fig


# ------------------------------------------------------------------ phase averaging / branches (sine)
def phase_average(d, n_bins):
    b = np.minimum((d["ph"].to_numpy() * n_bins).astype(int), n_bins - 1)
    g = (d.groupby(b)[["speed_phys", "current_phys", "speed_V", "current_V"]].median().reindex(range(n_bins)))
    g["ph"] = (np.arange(n_bins) + 0.5) / n_bins
    return g


def branch_of_phase(ph):
    """+1: commanded speed increasing. -1: commanded speed decreasing."""
    return np.where(np.cos(2 * np.pi * np.asarray(ph)) >= 0, 1, -1)


def air_current_at_phase_raw(air, frazil):
    """For every frazil point: air current at the same phase inside the same cycle (periodic interpolation).
    NaN where the air run has no usable data in that cycle."""
    out = np.full(len(frazil), np.nan)
    air_groups = {int(c): g for c, g in air.groupby("cycle")}
    for cycle in np.unique(frazil["cycle"]):
        if cycle not in air_groups:
            continue
        ia = air_groups[cycle]
        mask = frazil["cycle"].to_numpy() == cycle
        phase_f = frazil.loc[mask, "ph"].to_numpy()
        phase_a, current_a = ia["ph"].to_numpy(), ia["current_phys"].to_numpy()
        valid = np.isfinite(phase_a) & np.isfinite(current_a)
        phase_a, current_a = phase_a[valid], current_a[valid]
        if len(phase_a) < 2:
            continue
        order = np.argsort(phase_a)
        phase_a, current_a = phase_a[order], current_a[order]
        phase_a, unique_idx = np.unique(phase_a, return_index=True)
        current_a = current_a[unique_idx]
        if len(phase_a) < 2:
            continue
        phase_ext = np.concatenate([phase_a - 1.0, phase_a, phase_a + 1.0])   # periodic: phase 0 / 1 handled
        out[mask] = np.interp(phase_f, phase_ext, np.tile(current_a, 3), left=np.nan, right=np.nan)
    return out


def air_current_at_speed(air_w, speed, branch):
    """Air current at the requested speed, accelerating / decelerating branches separately. NaN outside the air range."""
    out = np.full(len(speed), np.nan)
    for b in (1, -1):
        a = air_w[branch_of_phase(air_w["ph"]) == b]
        if len(a) < 10:
            continue
        edges = np.linspace(a["speed_phys"].min(), a["speed_phys"].max(), N_SPEED_BINS + 1)
        idx = np.clip(np.digitize(a["speed_phys"], edges) - 1, 0, N_SPEED_BINS - 1)
        g = a.groupby(idx)["current_phys"].median()
        centers = 0.5 * (edges[:-1] + edges[1:])
        m = branch == b
        out[m] = np.interp(speed[m], centers[g.index.to_numpy()], g.to_numpy(), left=np.nan, right=np.nan)
    return out


# ------------------------------------------------------------------ normalisation methods (sine fig 4)
NORM_INFO = {
    "ratio": "I / I_air",
    "diff": "I - I_air (same phase)",
    "diff_scaled": "(I - I_air) / half-range of I_air",
    "diff_speed": "I - I_air (same speed and same branch)",
    "regression": "I - (a * I_air + b)",
}
NORM_TITLE = {"ratio": "ratio", "diff": "diff", "diff_scaled": "diff scaled",
              "diff_speed": "diff at same speed", "regression": "regression"}


def normalise(dd, dur, runs, method):
    """-> (table with columns ph, speed, value ; y label ; reference level)."""
    wf = window(dd["frazil"], SKIP_CYCLES / runs["frazil"]["freq"], dur).copy()
    wa = window(dd["air"], SKIP_CYCLES / runs["air"]["freq"], dur).copy()

    if PHASE_AVERAGE_FIG4:
        gf, ga = phase_average(wf, N_PHASE_BINS), phase_average(wa, N_PHASE_BINS)
        i_f, i_a = gf["current_phys"].to_numpy(), ga["current_phys"].to_numpy()
        out = pd.DataFrame({"ph": gf["ph"].to_numpy(), "speed": gf["speed_phys"].to_numpy()})
    else:
        i_f = wf["current_phys"].to_numpy()
        i_a = air_current_at_phase_raw(wa, wf)
        out = pd.DataFrame({"ph": wf["ph"].to_numpy(), "speed": wf["speed_phys"].to_numpy()})

    if method == "ratio":
        max_air = np.nanmax(np.abs(i_a))
        ok = np.isfinite(i_a) & (np.abs(i_a) >= RATIO_MIN_FRAC * max_air)
        values = np.full(len(i_f), np.nan)
        valid = ok & np.isfinite(i_f) & np.isfinite(i_a)
        values[valid] = i_f[valid] / i_a[valid]
        out["value"] = values
        print(f"  ratio: {int((~valid).sum())} point(s) masked (missing or small air current).")
        return out, r"$I / I_{\mathrm{air}}$ (-)", 1.0

    if method == "diff":
        out["value"] = i_f - i_a
        return out, r"$\Delta I$ (A) $\propto \Delta\tau$", 0.0

    if method == "diff_scaled":
        valid_air = i_a[np.isfinite(i_a)]
        scale = (np.nanquantile(valid_air, 0.95) - np.nanquantile(valid_air, 0.05)) / 2
        if not np.isfinite(scale) or scale == 0:
            raise ValueError("diff_scaled: invalid air-current half-range.")
        out["value"] = (i_f - i_a) / scale
        print(f"  diff_scaled: air-current half-range = {scale:.3g} A")
        return out, r"$\Delta I / \delta I_{\mathrm{air}}$ (-)", 0.0

    if method == "regression":
        ok = np.isfinite(i_f) & np.isfinite(i_a)
        if ok.sum() < 2:
            raise ValueError("regression: not enough valid measurements.")
        a, b = np.polyfit(i_a[ok], i_f[ok], 1)
        out["value"] = i_f - (a * i_a + b)
        print(f"  regression: I_frazil = {a:.3f} * I_air + {b:.3f} A (fitted on {int(ok.sum())} points)")
        return out, r"$I - (a\,I_{\mathrm{air}} + b)$ (A)", 0.0

    if method == "diff_speed":
        i_ref = air_current_at_speed(wa, wf["speed_phys"].to_numpy(), branch_of_phase(wf["ph"]))
        resid = wf["current_phys"].to_numpy() - i_ref
        print(f"  diff_speed: {int(np.isnan(resid).sum())} of {len(resid)} frazil points outside the air speed range (dropped).")
        if PHASE_AVERAGE_FIG4:
            b = np.minimum((wf["ph"].to_numpy() * N_PHASE_BINS).astype(int), N_PHASE_BINS - 1)
            t = (pd.DataFrame({"speed": wf["speed_phys"].to_numpy(), "value": resid})
                 .groupby(b).median().reindex(range(N_PHASE_BINS)))
            out["value"] = t["value"].to_numpy()
        else:
            out["value"] = resid
        return out, r"$\Delta I(\omega)$ (A) $\propto \Delta\tau$", 0.0

    raise ValueError(f"Unknown normalisation method '{method}'")


def choose_methods():
    names = list(NORM_INFO)
    sel = NORM_METHOD
    if sel is None:
        print("\n" + "=" * 70)
        print("FIGURE 4 - NORMALISATION BY THE AIR REFERENCE")
        print("=" * 70)
        for i, k in enumerate(names, 1):
            print(f"\n  {i}) {k}\n     {NORM_INFO[k]}")
        print("\n  a) All five methods\n\nSelect 1 to 4 methods.  Examples: 2 | 2 4 | 1 2 3 4 | a")
        raw = input("\nYour choice [2 4]: ").strip() or "2 4"
        if raw.lower() in ("a", "all"):
            sel = list(names)
        else:
            sel = []
            for x in raw.replace(",", " ").split():
                if x.isdigit():
                    i = int(x)
                    if not 1 <= i <= len(names):
                        raise SystemExit(f"ERROR: method number {i} is invalid. Choose 1 to {len(names)}.")
                    sel.append(names[i - 1])
                else:
                    sel.append(x)
    elif isinstance(sel, str):
        sel = list(names) if sel.lower() == "all" else [sel]
    else:
        sel = list(sel)
    bad = [m for m in sel if m not in names]
    if bad:
        raise SystemExit(f"ERROR: unknown methods {bad}. Choose among {names}.")
    if not sel:
        raise SystemExit("ERROR: select at least one method.")
    if len(set(sel)) != len(sel):
        raise SystemExit("ERROR: duplicate methods selected.")
    if len(sel) > 4 and sel != names:
        raise SystemExit("ERROR: select a maximum of four methods, or explicitly choose 'all'.")
    print("\nFigure 4 methods: " + ", ".join(sel))
    return sel


# ------------------------------------------------------------------ fig 4 (sine)
def fig_norm(dd, dur, runs, methods):
    if "current_phys" not in dd["frazil"] or "current_phys" not in dd["air"]:
        print("Fig 4 skipped: needs the .edc (physical current) for both runs.")
        return None
    n = len(methods)
    if n == 1:
        nrows, ncols, figsize = 1, 1, FIGSIZE_1
    elif n == 2:
        nrows, ncols, figsize = 1, 2, FIGSIZE_2COLS
    else:
        nrows, ncols = int(np.ceil(n / 2)), 2
        figsize = (FIGSIZE_2COLS[0], FIGSIZE_2COLS[1] * 1.7)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, layout="constrained", squeeze=False)
    axes = axes.ravel()
    sc = None
    for ax, method in zip(axes, methods):
        print(f"\nFig 4 [{method}]: {NORM_INFO[method]}")
        tab, ylabel, ref = normalise(dd, dur, runs, method)
        v = tab["value"].to_numpy()
        valid = np.isfinite(v) & np.isfinite(tab["speed"].to_numpy())
        if int(valid.sum()) == 0:
            print(f"  WARNING: no valid points for {method}.")
            ax.set_title(NORM_TITLE[method])
            ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
            continue
        print(f"  Points plotted: {int(valid.sum())} / {len(tab)}  "
              + (f"(phase averaging, {N_PHASE_BINS} bins)" if PHASE_AVERAGE_FIG4 else "(raw, no phase averaging)"))
        if PHASE_AVERAGE_FIG4:
            ax.plot(tab["speed"], v, color="0.6", lw=1.5, zorder=1)
        sc = ax.scatter(tab.loc[valid, "speed"], v[valid], c=tab.loc[valid, "ph"], cmap=CMAP_PHASE, vmin=0, vmax=1,
                        s=(6 * POINT_SIZE if PHASE_AVERAGE_FIG4 else POINT_SIZE), zorder=2,
                        edgecolor=("k" if PHASE_AVERAGE_FIG4 else "none"), linewidth=0.6)
        ax.axhline(ref, color="k", ls="--", lw=1.5)
        ax.set_xlabel(labels()["speed"])
        ax.set_ylabel(ylabel)
        if n > 1:
            ax.set_title(NORM_TITLE[method])
        ax.grid(True, alpha=0.3)
    for ax in axes[n:]:
        ax.set_visible(False)
    if sc is not None:
        cb = fig.colorbar(sc, ax=list(axes[:n]), pad=0.02)
        cb.set_label(LABEL_PHASE, fontsize=FONT_LABEL)
        cb.ax.tick_params(labelsize=FONT_TICK)
    mode = "phaseavg" if PHASE_AVERAGE_FIG4 else "raw"
    save_fig(fig, "fig4_current_norm_" + "_".join(methods) + f"_{mode}")
    return fig


# ------------------------------------------------------------------ main: ramp test
def main_ramp(ref=False):
    run = load_meta_ramp(*RAMP_RUN, "ramp")
    runs = {"ramp": run}
    run_air = None
    same_cmds = True
    force_time = force_time_requested()
    if ref:
        run_air = load_meta_ramp(*AIR_RAMP_RUN, "air")
        runs["air"] = run_air
        same_cmds = check_same_input_ramp(run, run_air)
    edcs = load_run_edcs(runs, ref)

    # Keep every decoded USB sample: no spike/glitch filter is applied.
    d = decode_run(run, edcs["ramp"])
    d_air = decode_run(run_air, edcs["air"]) if ref else None
    r = run["ramp"]
    sx, sy = cols()
    lab = labels()

    for ph, name in ((1, f"{r['dir1']} ramp + 1st half of plateau"), (2, f"{r['dir2']} ramp + 2nd half of plateau")):
        w = d[d["phase"] == ph]
        if len(w):
            print(f"[ramp] phase {ph} ({name}): {w['t'].iloc[0]:.1f} -> {w['t'].iloc[-1]:.1f} s, "
                  f"speed {w[sx].min():.3g}..{w[sx].max():.3g}, current {w[sy].min():.3g}..{w[sy].max():.3g}  "
                  f"({len(w)} points)")
    if ref:
        ra = run_air["ramp"]
        for ph, name in ((1, ra["dir1"]), (2, ra["dir2"])):
            w = d_air[d_air["phase"] == ph]
            if len(w):
                print(f"[air]  phase {ph} ({name}): {w['t'].iloc[0]:.1f} -> {w['t'].iloc[-1]:.1f} s, "
                      f"speed {w[sx].min():.3g}..{w[sx].max():.3g}, current {w[sy].min():.3g}..{w[sy].max():.3g}  "
                      f"({len(w)} points)")

    command = None
    if SHOW_COMMAND and REAL_UNITS and edcs["ramp"] is not None:
        c = run["cmds"]
        command = (c["t_rel"].to_numpy(), np.asarray(edcs["ramp"].setpoint_from_duty(c["duty"].to_numpy()), float))

    figs = [fig_ramp_time(d, run, sx, lab["speed"], "fig1_ramp_speed_vs_time", command),
            fig_ramp_time(d, run, sy, lab["current"], "fig2_ramp_current_vs_time"),
            fig_ramp_loop(d, run, d_air)]

    # ---- fig 4: raw frazil current minus air current at the same elapsed time
    if ref:
        if "current_phys" not in d or "current_phys" not in d_air:
            print("Fig 4 skipped: needs the .edc (physical current) for both runs.")
        elif not same_cmds and not force_time:
            print("\nFig 4 NOT drawn: the two runs' command timing differs by more than "
                  f"{RAMP_T_TOL_S:g} s. Re-run with --force-time to draw the raw time subtraction anyway.")
        else:
            print("\n" + "=" * 65 + "\nFIGURE 4 (ramp): RAW FRAZIL MINUS AIR\n" + "=" * 65)
            if not same_cmds and force_time:
                print("  WARNING: --force-time requested. Subtracting the air current at the same elapsed time "
                      "despite the command/timing mismatch; interpret the result cautiously.")
            table = normalise_ramp(d, d_air, method="time")
            if len(table):
                print("Fig 4: drawing raw time-matched subtraction (no filtering, binning, or averaging).")
                figs.append(fig_norm_ramp({"time": table}, run))
            else:
                print("Fig 4: no valid points in the overlapping time range.")

    if not SAVE_PLOTS:
        print("\nSAVE_PLOTS is False: nothing saved.")
    plt.show()


# ------------------------------------------------------------------ main
def main():
    setup_style()
    test = get_test()
    mode = get_mode()
    ref = mode == "ref"
    if test == "ramp_test":
        return main_ramp(ref)

    keys = ("frazil", "air") if ref else ("frazil",)
    if not ref:
        print("no_ref: the AIR folder is ignored (not even opened).")

    # 1. metadata of the runs that are needed
    runs = {k: load_meta(*SOURCES[k], k) for k in keys}
    if ref:
        check_same_input(runs["frazil"], runs["air"])

    # 2. .edc of each run (needed for rpm / A, and for fig 4)
    edcs = load_run_edcs(runs, ref)

    # 3. decode
    dd = {k: decode_run(runs[k], edcs[k]) for k in keys}

    # 4. duration of the analysed part
    T = 1.0 / runs["frazil"]["freq"]
    n_cyc = min(r["n_cycles"] for r in runs.values()) if (ref and TRIM_TO_SAME_DURATION) else None
    dur = {k: (n_cyc if n_cyc is not None else runs[k]["n_cycles"]) * T for k in keys}
    dur_all = min(dur.values()) if (ref and TRIM_TO_SAME_DURATION) else max(dur.values())
    print(f"Sine: {runs['frazil']['freq']:g} Hz ({runs['frazil']['note']}), "
          f"analysis from cycle {SKIP_CYCLES} to {dur_all / T:g}.")

    # 5. sanity numbers
    sx, sy = cols()
    for k in keys:
        w = window(dd[k], SKIP_CYCLES * T, dur_all)
        amp = lambda s: (s.quantile(0.95) - s.quantile(0.05)) / 2
        print(f"[{k}] speed: mean {w[sx].mean():.3g}, half-range {amp(w[sx]):.3g} | "
              f"current: mean {w[sy].mean():.3g}, half-range {amp(w[sy]):.3g}   ({len(w)} points)")

    # 6. figures 1-3
    lab = labels()
    figs = [fig_time(dd, keys, sx, lab["speed"], "Speed vs time", "fig1_speed_vs_time", dur_all),
            fig_time(dd, keys, sy, lab["current"], "Current vs time", "fig2_current_vs_time", dur_all),
            fig_loop(dd, keys, dur_all, runs)]

    # 7. figure 4 (ref mode only)
    if ref:
        print("\n" + "=" * 65 + "\nFIGURE 4 DATA PROCESSING\n" + "=" * 65)
        if PHASE_AVERAGE_FIG4:
            print(f"Phase averaging ENABLED: {N_PHASE_BINS} bins per cycle.")
        else:
            print("Phase averaging DISABLED: raw frazil points; air current interpolated at the same phase, same cycle.")
        figs.append(fig_norm(dd, dur_all, runs, choose_methods()))

    if not SAVE_PLOTS:
        print("\nSAVE_PLOTS is False: nothing saved.")
    plt.show()


if __name__ == "__main__":
    main()