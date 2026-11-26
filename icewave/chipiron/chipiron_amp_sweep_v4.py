"""
chipiron_amp_sweep_v4.py -- sine at FIXED frequency, amplitude swept (rpm). Firmware mode only
(Arduino mode mo + bipolar .edc, or any mean that keeps mean +/- amp inside 10..90 % in mo).

TWO KINDS OF SWEEP
  STEPPED     (--amps / --amp-range)   amplitude constant during each segment of N whole cycles
  CONTINUOUS  (--sweep LO HI)          ONE command "r<mean>,<amp0>,<amp1>,<freq>,<T>,<shape>": the Arduino
                                       (PWM_full_v4.ino) changes the amplitude smoothly from LO to HI in T seconds,
                                       linearly or geometrically (--log). No HTTP jitter, no steps. After T the sine
                                       keeps running at HI; this script then sends x (rest).
                                       Needs PWM_full_v4.ino (the r command does not exist in v3).

STEPPED: HOW
  For each amplitude A_k: send ONE "o<mean_duty>,<amp_duty>,<freq>" and wait one SEGMENT, then the next.
  Segment k starts at t0 + k*(T + rest), T = cycles / freq, so the timing does not drift. Because T is an integer
  number of periods, the old sine ends at the mean going up and the new one starts there too: no jump in the
  command (up to the HTTP jitter, tens of ms). Optional rest (x) between amplitudes. Always ends with x (rest).

INPUT LOG
  Every command is written in the input log by the helpers (env CHIPIRON_INPUT_LOG, set by the session).
  Stepped: cut the recording into one window per amplitude, discard the first cycle(s) of each segment.
  Continuous: ONE command; the amplitude at time t after it is
      linear  amp(t) = A0 + (A1 - A0) * t / T        log  amp(t) = A0 * (A1 / A0) ** (t / T)      (t <= T)
  (t counted from the logged time of the r command + the HTTP latency, a few tens of ms).

DURATION
  Stepped (choose AT MOST ONE):
    --cycles N       N periods per amplitude (default 5; integer). Time per amplitude = N / freq.
    --seg-time S     about S s per amplitude: N = round(S * freq), at least 1 cycle.
    --total-time S   about S s for the WHOLE sweep (rests included), rounded to whole cycles as above.
  Continuous (choose ONE):
    --sweep-time S   duration of the sweep [s]
    --cycles N       N periods: T = N / freq

OPTIONS
  --freq F               frequency [Hz] (required)
  --mean-rpm M           mean speed [rpm] (default 0)
  --amps A1 A2 ...       stepped: amplitudes [rpm], in this order
  --amp-range LO HI N    stepped: N amplitudes from LO to HI, linear (or logarithmic with --log)
  --sweep LO HI          continuous: amplitude from LO to HI [rpm] (linear, or geometric with --log)
  --sweep-time S         continuous: duration [s]
  --log                  logarithmic spacing / geometric envelope (LO, HI > 0)
  --rest-between S       stepped: rest at 0 rpm between amplitudes [s] (default 0)
  --check-only           all the checks, print the plan and exit: NO command sent, no phone contact
  --yes                  no question (ramp-limiter warning accepted, no "Run it?")
  --no-record            do not call /start and /stop on the phone (the session records)
  --link http|lora, --phone-ip IP, --serial-port PORT, --edc FILE, --verbose   as in chipiron_oscillate_v4.py

CHECKS (before anything moves)
  * mean +/- amplitude inside 10..90 % duty for EVERY amplitude (continuous: both ends);
  * peak acceleration amp*2*pi*f against the ESCON speed ramp of the .edc (continuous: at the largest amplitude):
    above it the limiter flattens the sine -> warning and confirmation (skipped with --yes);
  * amplitude below 4 PWM counts (0.39 % duty each): warning, the sine is very coarse.

EXAMPLES  (use `py` on Windows, `python3` on Linux)
  # stepped: amplitudes typed, 0.1 Hz, 10 cycles each (100 s per amplitude)
  py chipiron_amp_sweep_v4.py --freq 0.1 --amps 100 200 300 400 --cycles 10 --no-record --yes
  # stepped: 8 log-spaced amplitudes 50..600 rpm, whole sweep about 20 min, 10 s rest between
  py chipiron_amp_sweep_v4.py --freq 0.2 --amp-range 50 600 8 --log --total-time 1200 --rest-between 10 \\
        --no-record --yes
  # CONTINUOUS: amplitude 50 -> 600 rpm, geometric, over 10 min, 0.1 Hz (60 cycles)
  py chipiron_amp_sweep_v4.py --freq 0.1 --sweep 50 600 --log --sweep-time 600 --no-record --yes
  # CONTINUOUS, linear 0 -> 300 rpm in 40 cycles of 0.2 Hz (200 s), only look at the plan
  py chipiron_amp_sweep_v4.py --freq 0.2 --sweep 0 300 --cycles 40 --check-only
"""
import argparse
import math
import time

import chipiron_helpers_v4 as cm


def gen_parser():
    p = argparse.ArgumentParser(description="Fixed-frequency sine, amplitude sweep (firmware mode)")
    p.add_argument("--freq", type=float, required=True, help="frequency [Hz]")
    p.add_argument("--mean-rpm", type=float, default=0.0, help="mean speed [rpm]")
    p.add_argument("--amps", type=float, nargs="+", help="stepped: amplitudes [rpm], in the order they are run")
    p.add_argument("--amp-range", type=float, nargs=3, metavar=("LO", "HI", "N"), help="stepped: N amplitudes from LO to HI")
    p.add_argument("--sweep", type=float, nargs=2, metavar=("LO", "HI"),
                   help="CONTINUOUS sweep: amplitude from LO to HI [rpm] (needs PWM_full_v4.ino)")
    p.add_argument("--sweep-time", type=float, default=None, help="continuous: duration of the sweep [s]")
    p.add_argument("--log", action="store_true", help="log spacing (stepped) / geometric envelope (continuous)")
    p.add_argument("--cycles", type=float, default=None, help="cycles per amplitude (stepped) or in the sweep (continuous)")
    p.add_argument("--seg-time", type=float, default=None, help="stepped: about this many seconds per amplitude")
    p.add_argument("--total-time", type=float, default=None, help="stepped: about this many seconds for the whole sweep")
    p.add_argument("--rest-between", type=float, default=0.0, help="stepped: rest at 0 rpm between amplitudes [s]")
    p.add_argument("--edc", default=None)
    p.add_argument("--link", choices=["http", "lora"], default="http")
    p.add_argument("--phone-ip", default=None)
    p.add_argument("--serial-port", default=None)
    p.add_argument("--no-record", action="store_true", help="don't call /start and /stop (session records)")
    p.add_argument("--check-only", action="store_true", help="checks + plan, then exit: nothing is sent")
    p.add_argument("--yes", action="store_true", help="no confirmation questions")
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


def amplitudes(a):
    if a.amps:
        return list(a.amps)
    if a.amp_range:
        lo, hi, n = a.amp_range[0], a.amp_range[1], int(a.amp_range[2])
        if n < 2:
            raise SystemExit("--amp-range needs N >= 2.")
        if a.log:
            if lo <= 0 or hi <= 0:
                raise SystemExit("log spacing needs LO and HI > 0.")
            return [lo * (hi / lo) ** (i / (n - 1)) for i in range(n)]
        return [lo + (hi - lo) * i / (n - 1) for i in range(n)]
    raise SystemExit("give --amps, --amp-range or --sweep.")


def pick_cycles(a, n):
    """Stepped: cycles per amplitude from --cycles / --seg-time / --total-time (at most one of them)."""
    if sum(x is not None for x in (a.cycles, a.seg_time, a.total_time)) > 1:
        raise SystemExit("give only ONE of --cycles, --seg-time, --total-time.")
    if a.seg_time is not None or a.total_time is not None:
        seg = a.seg_time if a.seg_time is not None else (a.total_time - (n - 1) * a.rest_between) / n
        if seg <= 0:
            raise SystemExit("--total-time is too short for this number of amplitudes and rests.")
        cycles = max(1, round(seg * a.freq))
        print(f"Requested ~{seg:.1f} s per amplitude -> {cycles} cycle(s) = {cycles / a.freq:.1f} s "
              "(whole number of periods).")
        return float(cycles)
    cycles = 5.0 if a.cycles is None else a.cycles
    if abs(cycles - round(cycles)) > 1e-6:
        print(f"WARNING: {cycles:g} cycles is not an integer: the sine will jump at each amplitude change.")
    return cycles


def sweep_duration(a):
    """Continuous: duration T [s] from --sweep-time or --cycles (exactly one)."""
    if (a.sweep_time is None) == (a.cycles is None):
        raise SystemExit("continuous sweep: give ONE of --sweep-time or --cycles.")
    if a.seg_time is not None or a.total_time is not None or a.rest_between:
        raise SystemExit("--seg-time / --total-time / --rest-between belong to the stepped sweep.")
    T = a.sweep_time if a.sweep_time is not None else a.cycles / a.freq
    if T <= 0:
        raise SystemExit("the sweep duration must be > 0.")
    return T


def check_amplitudes(a, cfg, mean_duty, amps):
    """Checks on every amplitude in `amps` (continuous: the two ends). Exits on error, asks on ramp warning."""
    ramp = min((r for r in (cfg.accel, cfg.decel) if r), default=None)
    hw_step_rpm = abs(cfg.slope) * 100 / cm.PWM_COUNTS      # one Arduino PWM count in rpm
    error, ramp_warn = False, False
    for A in amps:
        d = abs(A / cfg.slope)
        if mean_duty - d < cm.DUTY_MIN or mean_duty + d > cm.DUTY_MAX:
            print(f"ERROR: amplitude {A:g} rpm -> duty {mean_duty - d:.1f}..{mean_duty + d:.1f} %, outside "
                  f"{cm.DUTY_MIN:g}..{cm.DUTY_MAX:g} %.")
            error = True
        acc = abs(A) * 2 * math.pi * a.freq
        if ramp and acc > ramp:
            print(f"WARNING: {A:g} rpm at {a.freq:g} Hz -> peak acceleration {acc:.0f} rpm/s > ESCON ramp {ramp:g} "
                  f"rpm/s (max amplitude {ramp / (2 * math.pi * a.freq):.0f} rpm): the sine will be flattened.")
            ramp_warn = True
        if abs(A) < 4 * hw_step_rpm:
            print(f"WARNING: {A:g} rpm is only {abs(A) / hw_step_rpm:.1f} PWM counts ({hw_step_rpm:.1f} rpm each): "
                  "the sine will be very coarse.")
    if error:
        raise SystemExit("aborted (see above).")
    if ramp_warn and not a.yes and not cm.ask_yn("Continue anyway?", "n"):
        raise SystemExit("aborted.")


def run_continuous(a, cfg, mean_duty):
    lo, hi = a.sweep
    if a.amps or a.amp_range:
        raise SystemExit("--sweep is continuous: do not combine it with --amps / --amp-range.")
    if lo < 0 or hi < 0:
        raise SystemExit("--sweep amplitudes must be >= 0.")
    if a.log and (lo <= 0 or hi <= 0):
        raise SystemExit("--log needs LO and HI > 0.")
    T = sweep_duration(a)
    check_amplitudes(a, cfg, mean_duty, [lo, hi])      # envelope is monotonic: the ends are the extremes
    n_cyc = T * a.freq
    print(f"CONTINUOUS sweep ({'geometric' if a.log else 'linear'}): amplitude {lo:g} -> {hi:g} rpm "
          f"({abs(lo / cfg.slope):.2f} -> {abs(hi / cfg.slope):.2f} % duty)")
    print(f"mean {a.mean_rpm:g} rpm, {a.freq:g} Hz, {T:.1f} s = {n_cyc:.1f} cycles ({T / 60:.1f} min); "
          "the sine then stays at the final amplitude until the rest command.")
    if a.check_only:
        print("Check OK (--check-only: nothing sent).")
        return
    if not a.yes and not cm.ask_yn("Run it?", "n"):
        raise SystemExit("aborted.")
    cmd = (f"r{mean_duty:.2f},{abs(lo / cfg.slope):.2f},{abs(hi / cfg.slope):.2f},"
           f"{a.freq:.4f},{T:.1f},{1 if a.log else 0}")
    link = cm.make_link(a.link, a.phone_ip, a.serial_port, record=not a.no_record, verbose=a.verbose)
    link.start()
    time.sleep(2)
    try:
        t0 = time.monotonic()
        link.send(cmd, note=f"continuous amp sweep {lo:g}->{hi:g} rpm ({'log' if a.log else 'lin'}), mean "
                            f"{a.mean_rpm:g}, {a.freq:g} Hz, {T:.1f} s")
        print(f"sent {cmd}")
        while True:
            left = t0 + T - time.monotonic()
            if left <= 0:
                break
            print(f"[{time.monotonic() - t0:7.1f} / {T:.0f} s]")
            time.sleep(min(10.0, left))
    except KeyboardInterrupt:
        print("\nInterrupted (Ctrl-C): going to rest.")
    finally:
        try:
            link.send("x", note="rest (end of continuous sweep)")
        except Exception as exc:
            print("Rest command failed - STOP THE MOTOR MANUALLY!", exc)
        link.stop()


def main(a):
    cfg = cm.get_config(a.edc)
    if not cfg.is_speed:
        raise SystemExit("the PWM set value is not a speed.")
    if a.freq <= 0 or (a.cycles is not None and a.cycles <= 0):
        raise SystemExit("--freq and --cycles must be > 0.")
    try:
        mean_duty = cfg.duty_for(a.mean_rpm)
    except ValueError as exc:
        raise SystemExit(f"mean speed: {exc}")

    if a.sweep:
        run_continuous(a, cfg, mean_duty)
        return
    if a.sweep_time is not None:
        raise SystemExit("--sweep-time is for the continuous sweep (--sweep LO HI).")

    amps = amplitudes(a)
    cycles = pick_cycles(a, len(amps))
    check_amplitudes(a, cfg, mean_duty, amps)

    T = cycles / a.freq
    T_seg = T + a.rest_between
    total = len(amps) * T + (len(amps) - 1) * a.rest_between
    print(f"{len(amps)} amplitudes [rpm]: " + ", ".join(f"{A:.4g}" for A in amps))
    print(f"mean {a.mean_rpm:g} rpm, {a.freq:g} Hz, {cycles:g} cycles = {T:.1f} s per amplitude, total {total:.0f} s "
          f"({total / 60:.1f} min)")
    if a.check_only:
        print("Check OK (--check-only: nothing sent).")
        return
    if not a.yes and not cm.ask_yn("Run it?", "n"):
        raise SystemExit("aborted.")

    link = cm.make_link(a.link, a.phone_ip, a.serial_port, record=not a.no_record, verbose=a.verbose)
    link.start()
    time.sleep(2)
    try:
        t0 = time.monotonic()
        for k, A in enumerate(amps):
            wait = t0 + k * T_seg - time.monotonic()
            if wait > 0:
                time.sleep(wait)
            link.send(f"o{mean_duty:.2f},{abs(A / cfg.slope):.2f},{a.freq:.4f}",
                      note=f"amp sweep {k + 1}/{len(amps)}: mean {a.mean_rpm:g}, amp {A:.4g} rpm, {a.freq:g} Hz, "
                           f"{cycles:g} cycles")
            print(f"[{time.monotonic() - t0:7.1f} / {total:.0f} s] amplitude {k + 1}/{len(amps)}: {A:.4g} rpm")
            if a.rest_between > 0 and k < len(amps) - 1:
                time.sleep(max(0.0, t0 + k * T_seg + T - time.monotonic()))
                link.send("x", note="rest between amplitudes")
        wait = t0 + (len(amps) - 1) * T_seg + T - time.monotonic()
        if wait > 0:
            time.sleep(wait)
    except KeyboardInterrupt:
        print("\nInterrupted (Ctrl-C): going to rest.")
    finally:
        try:
            link.send("x", note="rest (end of sweep)")
        except Exception as exc:
            print("Rest command failed - STOP THE MOTOR MANUALLY!", exc)
        link.stop()


if __name__ == "__main__":
    main(gen_parser())