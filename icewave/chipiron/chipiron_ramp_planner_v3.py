#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
chipiron_ramp_planner_v3.py -- plan a staircase ramp (linear or log) in rpm, rad/s or in % (motor command)
and save it as CSV for chipiron_run_v3.py.  

WHAT IT DOES
  * loads the ESCON .edc (or reuses the config exported by the session) to convert rpm <-> duty %, and prints
    the main parameters;
  * prints the raw ramp: requested value, duty requested, duty actually sent, speed actually obtained
    (including the Arduino resolution of 0.39 % duty) and the error;
  * warns if a plateau cannot be reached with the ESCON speed ramp, if points fall outside 10..90 %, or if
    rounding the duty merges points;
  * asks whether to save it as CSV, and where. The file name contains the key parameters of ramp and mapping.
    The CSV is read by chipiron_run_v3.py (mode load_ramp_planner, or --csv ... --yes).

DEPENDENCIES   chipiron_helpers_v3.py (-> escon_edc_v3.py), numpy.
USED BY        chipiron_run_v3.py (reads the CSV). The session script can pass such a CSV too ("reuse a CSV").

OPTIONS  (anything not given on the command line is asked)
  --edc FILE      ESCON .edc (ignored when called from the session)
  --unit rpm|rad/s|%      unit of the start / end values
  --type linear|log       spacing of the plateaus (log needs values > 0 and of the same sign)
  --start V --end V       first / last plateau, in --unit
  --points N              number of plateaus (>= 2)
  --time T                time per plateau [s]
  --folder DIR            folder of the CSV (default ramps)
  --name FILE             file name (default: built from the parameters)
  --no-save               only print the ramp, do not save

EXAMPLES  (use `py` on Windows, `python3` on Linux)
  # fully interactive
  py chipiron_ramp_planner_v3.py
  # logarithmic ramp 10 -> 500 rpm, 12 plateaus of 60 s, saved in ramps/ without any question
  py chipiron_ramp_planner_v3.py --edc parameters_Sep26.edc --unit rpm --type log --start 10 --end 500 \\
        --points 12 --time 60 --folder ramps
  # linear ramp in rad/s
  py chipiron_ramp_planner_v3.py --unit rad/s --type linear --start 1 --end 40 --points 10 --time 30
  # ramp directly in % duty (motor command), only look at it
  py chipiron_ramp_planner_v3.py --unit % --start 52 --end 80 --points 8 --time 20 --no-save
  # then run it
  py chipiron_run_v3.py --csv ramps/<file>.csv --yes
"""
import argparse
import math
import sys
from pathlib import Path

import numpy as np

import chipiron_helpers_v3 as cm
from chipiron_helpers_v3 import ask, ask_float, ask_int, ask_yn

UNITS = {"rpm": "rpm", "%": "pct", "rad/s": "rads"}     # unit -> token used in the file name


def plan_values(ramp_type, v0, v1, n):
    """n plateau values between v0 and v1, linearly or logarithmically spaced."""
    if ramp_type == "linear":
        return np.linspace(v0, v1, n)
    return np.geomspace(v0, v1, n)


def main():
    ap = argparse.ArgumentParser(description="Ramp planner (rpm / % / rad/s); see the module docstring")
    ap.add_argument("--edc", default=None, help="ESCON .edc file (ignored when called from the session)")
    ap.add_argument("--unit", choices=list(UNITS), default=None)
    ap.add_argument("--type", dest="ramp_type", choices=["linear", "log"], default=None)
    ap.add_argument("--start", type=float, default=None)
    ap.add_argument("--end", type=float, default=None)
    ap.add_argument("--points", type=int, default=None)
    ap.add_argument("--time", type=float, default=None, help="time per plateau [s]")
    ap.add_argument("--folder", default=None)
    ap.add_argument("--name", default=None)
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args()

    cfg = cm.get_config(args.edc)

    # ---------------------------------------------------------------- user input (asked only if not given)
    units = ["rpm", "rad/s", "%"] if cfg.is_speed else ["%"]
    unit = args.unit or ask(f"Plan the ramp in ({' / '.join(units)})", units[0])
    if unit not in units:
        sys.exit(f"Unit must be one of {units}.")
    ramp_type = args.ramp_type or ask("Ramp type (linear/log)", "linear").lower()
    if ramp_type not in ("linear", "log"):
        sys.exit("Ramp type must be 'linear' or 'log'.")
    v0 = args.start if args.start is not None else ask_float(f"Starting value [{unit}]")
    v1 = args.end if args.end is not None else ask_float(f"Ending value [{unit}]")
    n = args.points if args.points is not None else ask_int("Number of points", 10)
    T = args.time if args.time is not None else ask_float("Time per point [s]", 60)

    if n < 2:
        sys.exit("Number of points must be at least 2.")
    if T <= 0:
        sys.exit("Time per point must be > 0.")
    if ramp_type == "log" and (v0 <= 0 or v1 <= 0):
        sys.exit("Logarithmic ramps require values > 0 (for a bipolar mapping plan the ramp in rpm, not in %).")

    # ---------------------------------------------------------------- compute (in the chosen unit -> duty)
    req = plan_values(ramp_type, v0, v1, n)
    try:
        if unit == "%":
            duty_req = req
        else:
            rpm_req = req * 60 / (2 * math.pi) if unit == "rad/s" else req
            duty_req = np.array([cfg.duty_for(x) for x in rpm_req])
    except ValueError as exc:
        sys.exit(f"ERROR: {exc}")
    if np.any(duty_req < cm.DUTY_MIN - 1e-6) or np.any(duty_req > cm.DUTY_MAX + 1e-6):
        sys.exit(f"ERROR: some points are outside the valid duty range {cm.DUTY_MIN:g}..{cm.DUTY_MAX:g} %.")

    duty_sent = np.array([cm.q(d) for d in duty_req])
    rows = [{"dur": T, "duty": float(d), "kind": "hold", "seg": f"pt{i}"} for i, d in enumerate(duty_sent, 1)]
    dup = int(n - len(set(duty_sent.tolist())))

    # ---------------------------------------------------------------- output
    total_s = n * T
    print("\n" + "=" * 78)
    print("RAMP SUMMARY")
    print("=" * 78)
    print(f"Ramp type         : {ramp_type}")
    print(f"Planned in        : {unit}   ({v0:g} -> {v1:g})")
    print(f"Number of points  : {n}")
    print(f"Time per point    : {T:g} s")
    print(f"Total time        : {total_s:g} s = {total_s / 60:.2f} min = {total_s / 3600:.2f} h")
    if ramp_type == "log":
        dec = abs(math.log10(v1) - math.log10(v0))
        print(f"Decades           : {dec:.3f}   ({(n - 1) / dec:.2f} points/decade, "
              f"factor {req[1] / req[0]:.6g} between points)")
    else:
        print(f"Step              : {req[1] - req[0]:.6g} {unit}")

    print("\n  #   requested        duty req  duty sent" + ("   speed obtained   error" if cfg.is_speed else ""))
    for i in range(n):
        line = f"{i + 1:3d}  {req[i]:12.6g} {unit:<5} {duty_req[i]:7.2f}  {duty_sent[i]:8g}"
        if cfg.is_speed:
            got = cfg.value_at(cm.hw_duty(duty_sent[i]))     # what the Arduino really outputs (0.39 % steps)
            want = cfg.value_at(duty_req[i])
            err = 100 * (got - want) / want if want else 0.0
            line += f"   {got:10.2f} rpm  {err:+6.2f} %"
        print(line)
    print("=" * 78)

    if dup:
        print(f"WARNING: rounding the duty to {cm.DUTY_DECIMALS} decimal(s) merges {dup} point(s) "
              f"(1 % duty = {cfg.slope:g} {cfg.unit}). Use fewer points, a wider range, or a steeper ESCON scaling.")
    for w in cm.check_timeline(cfg, rows):
        print("WARNING:", w)

    # ---------------------------------------------------------------- save?
    if args.no_save or (args.folder is None and not ask_yn("\nSave this ramp as CSV?", "y")):
        print("Not saved.")
        return
    interactive = args.folder is None             # --folder given = non-interactive: default file name, no question
    folder = Path(args.folder or ask("Folder", "ramps"))
    stem = Path(cfg.edc_name).stem
    name = (f"ramp_{ramp_type}_{v0:g}to{v1:g}{UNITS[unit]}_{n}pts_{T:g}s_{cfg.kind}"
            f"_zero{cfg.zero_duty:g}pct_{cfg.slope:.4g}{cfg.unit}per-pct_{stem}.csv")
    name = cm.safe_name(args.name or (ask("File name", cm.safe_name(name)) if interactive else name))
    if not name.lower().endswith(".csv"):
        name += ".csv"
    path = folder / name
    k = 1
    while path.exists():                       # never overwrite
        path = folder / f"{Path(name).stem}_{k}.csv"
        k += 1
    meta = {"type": "ramp_planner", "ramp_type": ramp_type, "planned_in": unit, "start": f"{v0:g}",
            "end": f"{v1:g}", "n_points": n, "time_per_point_s": f"{T:g}"}
    cm.write_plan_csv(path, rows, cfg, meta)
    print(f"Saved: {path.resolve()}")
    print("Run it with:  py chipiron_run_v3.py --csv", f'"{path}"', "--yes")


if __name__ == "__main__":
    main()