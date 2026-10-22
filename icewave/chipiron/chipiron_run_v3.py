"""
chipiron_run_v3.py -- tunable motor program for Chipiron: holds, ramps, staircases, up/down cycles, repeats,
or a CSV made by chipiron_ramp_planner_v3.py.  

HOW IT WORKS
  1. loads the motor configuration (PWM <-> rpm mapping, rest duty, ESCON speed ramp...):
       - alone      : asks for the .edc (or --edc file) and loads it;
       - from the session (chipiron_session_v3.py): receives it through the env var CHIPIRON_MOTOR_CFG,
                      the .edc is NOT loaded again; only the main parameters are printed.
  2. builds a program = list of (duration, duty) commands, either by hand (mode manual, block by block) or from
     a CSV (mode load_ramp_planner);
  3. checks every plateau / ramp against the ESCON speed-ramp limiter of the .edc (warning + minimum time);
  4. optionally saves the program as CSV, then runs it: one "c<duty>" command per row, at the right time.
  The motor is always sent back to rest (the duty giving 0 rpm in the .edc: 50 % bipolar, 10 % unipolar) and
  the acquisition is stopped (try/finally), also on Ctrl-C or error. Commands are clamped to 10..90 %.

DEPENDENCIES
  chipiron_helpers_v3.py (config, links, CSV, ramp-limiter model)  ->  escon_edc_v3.py (only when reading an .edc)
  Arduino firmware: PWM_full_v3.ino (commands 'c<duty>' and 'x'; mode mc or mo, matching the .edc).
  Needs the phone HTTP server (http link) or the LoRa bridge (lora link, pyserial).
CALLED BY
  chipiron_session_v3.py (choice 3 and "reuse a CSV"), in two steps so no question is asked during the recording:
     py chipiron_run_v3.py --build-only programs/p.csv          (before the acquisition: build + save only)
     py chipiron_run_v3.py --csv programs/p.csv --yes --no-record   (during the acquisition: run, no question)

COMMAND LINE OPTIONS
  --mode manual|load_ramp_planner   how the program is made (asked if omitted; load_ramp_planner if --csv)
  --csv FILE            CSV from the ramp planner or a program saved by this script
  --yes                 with --csv: run immediately, no question at all
  --build-only OUT.csv  build the program, save it to OUT.csv and exit (no motor command)
  --edc FILE            ESCON .edc (ignored when called from the session)
  --dt SECONDS          time between two commands inside a ramp (default 1.0)
  --rest-before S / --rest-after S   rest (0 rpm) added before / after the program (asked if omitted)
  --link http|lora      how commands reach the Arduino (default http)
  --phone-ip IP         phone IP (default: env CHIPIRON_PHONE_IP, else chipiron_helpers_v3.DEFAULT_PHONE_IP)
  --serial-port PORT    LoRa bridge port, e.g. COM5 or /dev/ttyUSB0 (with --link lora)
  --no-record           do NOT call /start and /stop on the phone (acquisition handled elsewhere: the session)
  --no-reset            do not send the rest command at the end (leave the motor as is)
  -s / --song           play SONG (test.mp3) with `play` during the run
  --verbose             print every HTTP request

EXAMPLES  (use `py` on Windows, `python3` on Linux)
  # 1. everything interactive: asks the .edc, then the blocks, then "Run it now?"
  py chipiron_run_v3.py
  # 2. standalone with the .edc given, phone on another IP, no start/stop of the phone (session records)
  py chipiron_run_v3.py --edc parameters_Sep26.edc --phone-ip 10.153.247.36 --no-record
  # 3. run a staircase made by the ramp planner, asking nothing
  py chipiron_run_v3.py --csv ramps/ramp_log_10to500rpm_10pts_60s_bipolar.csv --yes
  # 4. load the planner CSV, but still add blocks / change rest times interactively
  py chipiron_run_v3.py --mode load_ramp_planner --csv ramps/myramp.csv
  # 5. prepare a program without moving the motor, run it later
  py chipiron_run_v3.py --build-only programs/test1.csv --rest-before 30 --rest-after 10
  py chipiron_run_v3.py --csv programs/test1.csv --yes
  # 6. over LoRa
  py chipiron_run_v3.py --csv programs/test1.csv --yes --link lora --serial-port COM5

BLOCKS OF THE MANUAL MODE (menu; values in rpm or %, 'u' switches the unit; every block can be repeated)
  1) hold       value 400 rpm for 60 s
                    > 1 | Value 400 | Duration 60 | Repeat 1
  2) ramp       from 0 to 600 rpm, by duration (e.g. 120 s) or by rate (e.g. 5 rpm/s), step = --dt
                    > 2 | From 0 | To 600 | [d]uration d | 120 | Repeat 1
  3) staircase  typed list "100 200 400" / linear 5 plateaus from 100 to 500 / logarithmic 8 from 10 to 800,
                time per plateau (same or one by one), optional mirror staircase (comes back down),
                optional rest at 0 between plateaus
                    > 3 | 3 (log) | First 10 | Last 800 | Plateaus 8 | mirror n | same time y | 60 | gap 0
  4) up/down cycles   A -> B -> A with ramp times and holds, N cycles
                    > 4 | A -200 | B 200 | up 20 s | hold 5 | down 20 s | hold 5 | cycles 3
  5) repeat the WHOLE program so far N times (optional pause at rest between runs)
  6) load a CSV from the ramp planner (or a saved program) and append it
  7) remove the last block      u) switch rpm <-> %      9) done
  After 9: rest before / after, FINAL PROGRAM summary with ESCON warnings, "Save as CSV?", "Run it now?".
"""
import argparse
import math
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import chipiron_helpers_v3 as cm
from chipiron_helpers_v3 import ask, ask_float, ask_int, ask_yn, q

SONG = "test.mp3"
RAMP_DT_DEFAULT = 1.0        # s between two commands during a ramp
REST_BEFORE_S = 30           # default rest (0 rpm) at the beginning of a program (baseline in the recording)
REST_AFTER_S = 0             # default rest at the end


# ================================================================== run the program
def execute(rows, cfg, link, song=False, reset=True):
    """Start the link, send each row's command at its scheduled time, then ALWAYS rest + stop."""
    events, t = [], 0.0
    for r in rows:
        events.append((t, r))
        t += r["dur"]
    total = t
    link.start()
    if song:
        subprocess.Popen(["play", SONG], text=True)
    try:
        t0 = time.monotonic()
        for t_ev, r in events:
            wait = t0 + t_ev - time.monotonic()
            if wait > 0:
                time.sleep(wait)
            link.send(f"c{cm.fmt_duty(r['duty'])}")
            rpm = f"~{cfg.value_at(r['duty']):.0f} {cfg.unit}" if cfg.is_speed else ""
            print(f"[{time.monotonic() - t0:7.1f} / {total:.0f} s] c{cm.fmt_duty(r['duty'])} {rpm:>12}  "
                  f"{r['kind']:<4} {r['seg']}")
        wait = t0 + total - time.monotonic()
        if wait > 0:
            time.sleep(wait)
    except KeyboardInterrupt:
        print("\nInterrupted (Ctrl-C): going to rest.")
    finally:
        try:
            if reset:
                link.send(f"c{cm.fmt_duty(cfg.rest_duty)}")
                print(f"Motor sent back to rest (c{cm.fmt_duty(cfg.rest_duty)}).")
        except Exception as exc:
            print(f"!! could not send the rest command: {exc}  -> STOP THE MOTOR MANUALLY")
        finally:
            link.stop()


# ================================================================== block builders
class Builder:
    """Interactive program builder. self.blocks = list of (label, rows)."""

    def __init__(self, cfg, dt):
        self.cfg = cfg
        self.dt = dt
        self.unit = "rpm" if cfg.is_speed else "%"
        self.blocks = []

    # ---- unit handling
    def to_duty(self, v):
        """Typed value (rpm or %) -> duty %, always checked against 10..90 % (ValueError otherwise)."""
        if self.unit == "%":
            if not cm.DUTY_MIN - 1e-6 <= v <= cm.DUTY_MAX + 1e-6:
                raise ValueError(f"{v:g} % is outside {cm.DUTY_MIN:g}..{cm.DUTY_MAX:g} %")
            return v
        return self.cfg.duty_for(v)

    def show(self, duty):
        s = f"{duty:.2f} %"
        return f"{self.cfg.value_at(duty):.1f} {self.cfg.unit} ({s})" if self.cfg.is_speed else s

    # ---- row factories
    def hold(self, duty, dur, seg):
        return [{"dur": dur, "duty": q(duty), "kind": "hold", "seg": seg}]

    def ramp(self, d0, d1, T, seg):
        """Staircase approximation of a linear ramp: one command every self.dt seconds."""
        n = max(1, math.ceil(T / self.dt - 1e-9))
        rows = []
        for k in range(n):
            t_end = min((k + 1) * self.dt, T)
            rows.append({"dur": t_end - k * self.dt, "duty": q(d0 + (d1 - d0) * t_end / T), "kind": "ramp", "seg": seg})
        return cm.merge_rows(rows)

    def ask_ramp_time(self, d0, d1, label="Ramp"):
        """Ask duration (or rate), warn if the ESCON limiter cannot follow."""
        cfg = self.cfg
        tmin = cm.min_time(cfg, cfg.value_at(d0), cfg.value_at(d1)) if cfg.is_speed else 0.0
        if tmin:
            print(f"  ESCON speed ramp: this ramp cannot be faster than {tmin:.1f} s "
                  f"(accel {cfg.accel:g}, decel {cfg.decel:g} rpm/s).")
        by_rate = ask("  Define the ramp by [d]uration or by [r]ate", "d").lower().startswith("r")
        while True:
            if by_rate:
                rate = ask_float(f"  Rate [{self.unit}/s]", 5)
                span = abs(d1 - d0) * (1 if self.unit == "%" else abs(cfg.slope))
                T = span / abs(rate) if rate else 0
                print(f"  -> duration {T:.1f} s")
            else:
                T = ask_float(f"  {label} duration [s]", max(round(tmin * 1.2, 1), 10.0))
            if T <= 0:
                print("  duration must be > 0.")
            elif tmin and T < tmin - 1e-6:
                print(f"  WARNING: {T:.1f} s < {tmin:.1f} s: the ESCON ramp limiter will flatten this ramp.")
                if not ask_yn("  Keep it anyway?", "n"):
                    continue
                return T
            else:
                return T

    def repeat(self, rows, label):
        """Ask how many times a block is repeated (optional pause at rest between repeats)."""
        n = ask_int("  Repeat this block how many times (1 = once)", 1)
        if n <= 1:
            return rows, label
        pause = ask_float("  Pause at rest between repeats [s]", 0)
        out = []
        for k in range(1, n + 1):
            out += [dict(r, seg=f"{r['seg']} [rep {k}/{n}]") for r in rows]
            if pause > 0 and k < n:
                out += self.hold(self.cfg.rest_duty, pause, "pause")
        return out, f"{label} x{n}"

    # ---- block types
    def add_hold(self):
        d = self.to_duty(ask_float(f"  Value [{self.unit}]"))
        T = ask_float("  Duration [s]", 60)
        label = f"hold {self.show(d)} for {T:g} s"
        return self.hold(d, T, label), label

    def add_ramp(self):
        d0 = self.to_duty(ask_float(f"  From [{self.unit}]", 0 if self.unit == "rpm" else self.cfg.rest_duty))
        d1 = self.to_duty(ask_float(f"  To [{self.unit}]"))
        T = self.ask_ramp_time(d0, d1)
        label = f"ramp {self.show(d0)} -> {self.show(d1)} in {T:.1f} s"
        return self.ramp(d0, d1, T, label), label

    def add_staircase(self):
        print("  Plateau values:  1) typed list   2) linear   3) logarithmic")
        how = ask("  Choice", "1")
        if how == "1":
            vals = [float(x.replace(",", ".")) for x in ask(f"  Values [{self.unit}], separated by spaces").split()]
        else:
            v0, v1 = ask_float(f"  First [{self.unit}]"), ask_float(f"  Last [{self.unit}]")
            n = ask_int("  Number of plateaus", 5)
            if n < 2:
                raise ValueError("a staircase needs at least 2 plateaus.")
            if how == "3":
                if v0 <= 0 or v1 <= 0:
                    raise ValueError("log spacing needs values > 0.")
                vals = [v0 * (v1 / v0) ** (i / (n - 1)) for i in range(n)]
            else:
                vals = [v0 + (v1 - v0) * i / (n - 1) for i in range(n)]
        if not vals:
            raise ValueError("no values.")
        duties = [self.to_duty(v) for v in vals]
        if ask_yn("  Also come back down (mirror staircase)?", "n"):
            duties += duties[-2::-1]
        if ask_yn("  Same time on every plateau?", "y"):
            Ts = [ask_float("  Time per plateau [s]", 60)] * len(duties)
        else:
            Ts = [ask_float(f"  Time at {self.show(d)} [s]", 60) for d in duties]
        gap = ask_float("  Rest at 0 between plateaus [s]", 0)
        rows = []
        for i, (d, T) in enumerate(zip(duties, Ts), 1):
            rows += self.hold(d, T, f"step {i}: {self.show(d)} {T:g} s")
            if gap > 0 and i < len(duties):
                rows += self.hold(self.cfg.rest_duty, gap, "gap")
        return rows, f"staircase of {len(duties)} plateaus"

    def add_cycle(self):
        a = self.to_duty(ask_float(f"  Low value A [{self.unit}]", 0 if self.unit == "rpm" else self.cfg.rest_duty))
        b = self.to_duty(ask_float(f"  High value B [{self.unit}]"))
        print("  A -> B ramp:")
        Tup = self.ask_ramp_time(a, b, "Up")
        hb = ask_float("  Hold at B [s]", 0)
        print("  B -> A ramp:")
        Tdn = self.ask_ramp_time(b, a, "Down")
        ha = ask_float("  Hold at A [s]", 0)
        n = ask_int("  Number of cycles", 1)
        rows = []
        for k in range(1, n + 1):
            rows += self.ramp(a, b, Tup, f"cycle {k}/{n} up")
            if hb > 0:
                rows += self.hold(b, hb, f"cycle {k}/{n} hold B")
            rows += self.ramp(b, a, Tdn, f"cycle {k}/{n} down")
            if ha > 0:
                rows += self.hold(a, ha, f"cycle {k}/{n} hold A")
        return rows, f"{n} cycle(s) {self.show(a)} <-> {self.show(b)}"

    def add_csv(self, path=None):
        """Append a CSV from the ramp planner (or a program saved by this script); checks the mapping."""
        path = path or ask("CSV from chipiron_ramp_planner_v3.py (or a saved program)", self.latest_csv())
        rows, meta = cm.read_plan_csv(path.strip('"'), self.cfg)
        print(f"  Loaded {len(rows)} rows, {cm.total_time(rows):g} s  ({meta.get('type', 'csv')}, "
              f"planned in {meta.get('planned_in', '?')}, mapping {meta.get('mapping', '?')}, edc {meta.get('edc', '?')})")
        bad = cm.mapping_mismatch(self.cfg, meta)
        if bad:
            print("  WARNING:", bad)
            if not ask_yn("  Use it anyway?", "n"):
                raise ValueError("CSV not used.")
        return rows, f"csv {Path(path).name}"

    @staticmethod
    def latest_csv():
        c = list(Path("ramps").glob("*.csv")) + list(Path("programs").glob("*.csv"))
        return str(max(c, key=lambda p: p.stat().st_mtime)) if c else ""

    # ---- menu
    def menu(self):
        actions = {"1": self.add_hold, "2": self.add_ramp, "3": self.add_staircase,
                   "4": self.add_cycle, "6": self.add_csv}
        while True:
            self.print_program()
            print(f"\nAdd a block  (values in {self.unit}; 'u' switches rpm <-> %):")
            print("  1) hold   2) ramp   3) staircase   4) up/down cycle(s)   6) load CSV (ramp planner)")
            print("  5) repeat the WHOLE program so far   7) remove last block   9) done")
            ch = ask("Choice", "9" if self.blocks else "1").lower()
            if ch == "9":
                if self.blocks:
                    return
                print("The program is empty.")
            elif ch == "u":
                if self.cfg.is_speed:
                    self.unit = "%" if self.unit == "rpm" else "rpm"
                    print(f"Now typing values in {self.unit}.")
                else:
                    print("Only % is available (set value is not a speed).")
            elif ch == "7":
                if self.blocks:
                    print("Removed:", self.blocks.pop()[0])
            elif ch == "5":
                rows = [r for _, rs in self.blocks for r in rs]
                if not rows:
                    continue
                n = ask_int("  Total number of runs of the whole program (including the existing one)", 2)
                pause = ask_float("  Pause at rest between runs [s]", 0)
                out = []
                for k in range(1, n + 1):
                    out += [dict(r, seg=f"{r['seg']} [run {k}/{n}]") for r in rows]
                    if pause > 0 and k < n:
                        out += self.hold(self.cfg.rest_duty, pause, "pause")
                self.blocks = [(f"whole program x{n}", out)]
            elif ch in actions:
                try:
                    rows, label = actions[ch]()
                    if ch != "6":
                        rows, label = self.repeat(rows, label)
                    self.blocks.append((label, rows))
                except (ValueError, OSError) as exc:
                    print(f"  ERROR: {exc}  (block not added)")
            else:
                print("Unknown choice.")

    def print_program(self):
        if not self.blocks:
            print("\nProgram: (empty)")
            return
        print("\nProgram so far:")
        t = 0.0
        for i, (label, rows) in enumerate(self.blocks, 1):
            d = cm.total_time(rows)
            print(f"  {i:2d}. {label:<60} {d:8.1f} s  (starts at {t:.0f} s)")
            t += d
        print(f"  Total: {t:.1f} s = {t / 60:.1f} min")


# ================================================================== finalise / save
def finalise(builder, rest_before=None, rest_after=None):
    """Add rest before/after, print summary and ESCON checks. Returns the final list of rows."""
    cfg = builder.cfg
    body = [r for _, rs in builder.blocks for r in rs]
    if rest_before is None:
        rest_before = ask_float("\nRest at 0 rpm BEFORE the program [s]", REST_BEFORE_S)
    if rest_after is None:
        rest_after = ask_float("Rest at 0 rpm AFTER the program [s]", REST_AFTER_S)
    rows = (builder.hold(cfg.rest_duty, rest_before, "rest before") if rest_before > 0 else []) + body + \
           (builder.hold(cfg.rest_duty, rest_after, "rest after") if rest_after > 0 else [])
    rows = cm.merge_rows(rows)
    duties = [r["duty"] for r in rows]
    print(f"\nFINAL PROGRAM: {len(rows)} commands, {cm.total_time(rows):.1f} s "
          f"({cm.total_time(rows) / 60:.1f} min), duty range {min(duties):g}..{max(duties):g} %"
          + (f", speed range {cfg.value_at(min(duties)):.0f}..{cfg.value_at(max(duties)):.0f} {cfg.unit}" if cfg.is_speed else ""))
    for w in cm.check_timeline(cfg, rows):
        print("WARNING:", w)
    return rows


def maybe_save(rows, cfg, path=None, how=None):
    """Save the program as CSV (asks unless a path is given)."""
    if path is None:
        if not ask_yn("Save this program as CSV (reloadable with --csv)?", "y"):
            return
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        default = f"programs/program_{stamp}_{cm.total_time(rows):.0f}s_{cfg.kind}_{Path(cfg.edc_name).stem}.csv"
        path = ask("File", default)
    p = cm.write_plan_csv(path, rows, cfg, {"type": "run_v3_program", "planned_in": how or "mixed"})
    print(f"Saved: {Path(p).resolve()}")


# ================================================================== main
def gen_parser():
    p = argparse.ArgumentParser(description="Run a tunable program on Chipiron (see the module docstring for examples)")
    p.add_argument("--mode", choices=["manual", "load_ramp_planner"], default=None)
    p.add_argument("--csv", default=None, help="CSV from the ramp planner / saved program")
    p.add_argument("--yes", action="store_true", help="with --csv: run without any question")
    p.add_argument("--build-only", dest="build_only", default=None, metavar="OUT.csv",
                   help="build the program, save it to OUT.csv and exit (no motor command)")
    p.add_argument("--edc", default=None, help="ESCON .edc (ignored when called from the session)")
    p.add_argument("--dt", type=float, default=RAMP_DT_DEFAULT, help="seconds between commands in a ramp")
    p.add_argument("--rest-before", type=float, default=None, help="rest at 0 rpm before the program [s]")
    p.add_argument("--rest-after", type=float, default=None, help="rest at 0 rpm after the program [s]")
    p.add_argument("--link", choices=["http", "lora"], default="http")
    p.add_argument("--phone-ip", default=None, help="phone IP (default: env CHIPIRON_PHONE_IP or the helpers default)")
    p.add_argument("--serial-port", default=None, help="LoRa bridge serial port (with --link lora)")
    p.add_argument("-s", "--song", dest="song", action="store_true", help="play SONG during the run")
    p.add_argument("--no-reset", dest="reset", action="store_false", help="don't send the rest command at the end")
    p.add_argument("--no-record", dest="record", action="store_false", help="don't call /start and /stop on the phone")
    p.add_argument("--verbose", action="store_true", help="print every HTTP request")
    return p.parse_args()


def main(args):
    cfg = cm.get_config(args.edc)
    builder = Builder(cfg, args.dt)

    def link():
        return cm.make_link(args.link, args.phone_ip, args.serial_port, args.record, args.verbose)

    # ---- non-interactive run of a saved program (used by the session script)
    if args.csv and args.yes:
        rows, meta = cm.read_plan_csv(args.csv, cfg)
        bad = cm.mapping_mismatch(cfg, meta)
        if bad:
            sys.exit(f"ERROR: {bad}")
        print(f"Program {args.csv}: {len(rows)} commands, {cm.total_time(rows):.1f} s")
        for w in cm.check_timeline(cfg, rows):
            print("WARNING:", w)
        execute(rows, cfg, link(), args.song, args.reset)
        return

    mode = args.mode or ("load_ramp_planner" if args.csv else None)
    if mode is None:
        print("Mode:\n  1) manual            - build the program block by block"
              "\n  2) load_ramp_planner - load a CSV made by chipiron_ramp_planner_v3.py")
        mode = {"1": "manual", "2": "load_ramp_planner"}.get(ask("Choice", "1"), "manual")
    if mode == "load_ramp_planner":
        try:
            rows, label = builder.add_csv(args.csv)
            builder.blocks.append((label, rows))
        except (ValueError, OSError) as exc:
            sys.exit(f"ERROR: {exc}")
    builder.menu()
    rows = finalise(builder, args.rest_before, args.rest_after)
    if args.build_only:
        maybe_save(rows, cfg, args.build_only, builder.unit)
        return
    maybe_save(rows, cfg, None, builder.unit)
    if ask_yn("Run it now?", "n"):
        execute(rows, cfg, link(), args.song, args.reset)


if __name__ == "__main__":
    main(gen_parser())