"""
chipiron_oscillate_v3.py -- send a sinusoidal speed command to Chipiron. /!\ LORA NOT TESTED.
    speed(t) = mean + amp * sin(2*pi*f*t)        realised as a slowly varying PWM duty cycle

TWO MODES (--mode)
  firmware (default)  ONE message "o<mean_duty>,<amp_duty>,<freq>": the Arduino builds the sine itself at 100 Hz.
                      Needs PWM_full_v3.ino in Arduino mode mo (oscillation, rest duty 50 % = bipolar .edc).
  stream              works with any firmware / mode (also unipolar .edc, Arduino mode mc): a new "c<duty>"
                      command every --dt seconds. Coarse (duty steps of 1 / DUTY_DECIMALS), only for periods
                      >= ~20 s.

TWO LINKS (--link)
  http (default)  laptop -> wifi -> phone HTTP server -> USB -> Arduino
  lora            laptop -> USB -> lora_tx_bridge Arduino ~~radio~~ motor Arduino (USE_LORA 1 in the firmware);
                  needs pyserial and --serial-port. The phone acquisition is NOT started by this link.

MAPPING rpm <-> duty
  Comes from the ESCON .edc (chipiron_helpers_v3.get_config): read from --edc, or received from the session
  (env CHIPIRON_MOTOR_CFG). Typical values of parameters_Sep26.edc:
      bipolar : 10 % -> -800 rpm, 50 % -> 0, 90 % -> +800 rpm
      unipolar: 10 % ->    0 rpm,            90 % -> +800 rpm
  Manual override (no .edc): give BOTH --rpm-at-10 and --rpm-at-90.

SAFETY CHECKS before moving the motor
  * mean +/- amp must stay inside 10..90 % duty (else error);
  * peak acceleration of the sine = amp * 2*pi*f  is compared with the ESCON speed ramp of the .edc: above it
    the limiter flattens the sine -> warning and confirmation (skipped with --yes, which the session adds);
  * on exit (also Ctrl-C / error) the motor goes back to rest: "x" (firmware mode) or c<rest duty> (stream).

DEPENDENCIES
  chipiron_helpers_v3.py (-> escon_edc_v3.py when reading an .edc), PWM_full_v3.ino on the Arduino.
CALLED BY
  chipiron_session_v3.py (choice 2) with --no-record --yes; can also be run alone.

EXAMPLES  (use `py` on Windows, `python3` on Linux)
  # sine around 0 rpm (bipolar .edc, Arduino in mode mo), +/-300 rpm, 0.5 Hz, 10 cycles, phone records itself
  py chipiron_oscillate_v3.py --edc parameters_Sep26.edc --mean-rpm 0 --amp-rpm 300 --freq 0.5 --cycles 10
  # same, acquisition already running (session) and no question
  py chipiron_oscillate_v3.py --mean-rpm 0 --amp-rpm 300 --freq 0.5 --cycles 10 --no-record --yes
  # sine around 400 rpm, +/-100 rpm, 0.2 Hz, 5 cycles
  py chipiron_oscillate_v3.py --mean-rpm 400 --amp-rpm 100 --freq 0.2 --cycles 5
  # old firmware / unipolar .edc: stream mode, a command every 1 s, period 50 s
  py chipiron_oscillate_v3.py --mode stream --dt 1 --mean-rpm 400 --amp-rpm 100 --freq 0.02 --cycles 3
  # over LoRa
  py chipiron_oscillate_v3.py --link lora --serial-port COM5 --mean-rpm 400 --amp-rpm 100 --freq 0.1
  # mapping typed by hand (no .edc), phone on another IP, with music
  py chipiron_oscillate_v3.py --rpm-at-10 -800 --rpm-at-90 800 --phone-ip 192.168.43.1 --song \\
        --mean-rpm 0 --amp-rpm 200 --freq 0.2
"""
import argparse
import math
import subprocess
import sys
import time

import chipiron_helpers_v3 as cm

SONG = "test.mp3"


def gen_parser():
    p = argparse.ArgumentParser(description="Sinusoidal speed command for Chipiron (see the module docstring)")
    p.add_argument("--link", choices=["http", "lora"], default="http")
    p.add_argument("--serial-port", help="LoRa bridge serial port, e.g. /dev/ttyUSB0 or COM5")
    p.add_argument("--phone-ip", default=None, help="phone IP (default: env CHIPIRON_PHONE_IP or helpers default)")
    p.add_argument("--mode", choices=["firmware", "stream"], default="firmware")
    p.add_argument("--mean-rpm", type=float, required=True, help="mean speed [rpm]")
    p.add_argument("--amp-rpm", type=float, required=True, help="amplitude [rpm]")
    p.add_argument("--freq", type=float, required=True, help="frequency [Hz]")
    p.add_argument("--cycles", type=float, default=5, help="number of periods (default 5)")
    p.add_argument("--dt", type=float, default=1.0, help="stream mode: seconds between commands")
    p.add_argument("--edc", default=None, help="ESCON .edc (ignored when called from the session)")
    p.add_argument("--rpm-at-10", type=float, default=None, help="manual mapping: rpm at 10 %% duty (with --rpm-at-90)")
    p.add_argument("--rpm-at-90", type=float, default=None, help="manual mapping: rpm at 90 %% duty (with --rpm-at-10)")
    p.add_argument("--song", action="store_true", help="play SONG with `play` during the run")
    p.add_argument("--yes", action="store_true", help="do not ask for confirmation when the ESCON ramp is too slow")
    p.add_argument("--no-record", action="store_true",
                   help="don't call /start and /stop on the phone (acquisition handled by the session)")
    p.add_argument("--verbose", action="store_true", help="print every HTTP request")
    return p.parse_args()


def get_cfg(args):
    """Motor configuration: manual mapping if both --rpm-at-* are given, else .edc / session."""
    if (args.rpm_at_10 is None) != (args.rpm_at_90 is None):
        raise SystemExit("give BOTH --rpm-at-10 and --rpm-at-90 (or neither, to use the .edc).")
    if args.rpm_at_10 is not None:
        if args.rpm_at_10 == args.rpm_at_90:
            raise SystemExit("--rpm-at-10 and --rpm-at-90 must differ.")
        slope = (args.rpm_at_90 - args.rpm_at_10) / 80.0
        cfg = cm.MotorConfig("manual", "(--rpm-at-10/90)", "speed", "rpm", args.rpm_at_10 - slope * 10, slope)
        print(cfg.summary())
        return cfg
    return cm.get_config(args.edc)


def main(args):
    cfg = get_cfg(args)
    if not cfg.is_speed:
        raise SystemExit("the PWM set value of this .edc is not a speed: nothing to oscillate in rpm.")
    if args.freq <= 0 or args.cycles <= 0:
        raise SystemExit("--freq and --cycles must be > 0.")
    try:
        mean_duty = cfg.duty_for(args.mean_rpm)
    except ValueError as exc:
        raise SystemExit(f"mean speed: {exc}")
    amp_duty = args.amp_rpm / cfg.slope

    # ---- sanity checks before moving a motor
    lo, hi = mean_duty - abs(amp_duty), mean_duty + abs(amp_duty)
    if lo < cm.DUTY_MIN or hi > cm.DUTY_MAX:
        raise SystemExit(f"Duty would leave [{cm.DUTY_MIN:g},{cm.DUTY_MAX:g}] %: {lo:.1f}..{hi:.1f} %. "
                         f"Speeds reachable with this .edc: {cfg.value_at(cm.DUTY_MIN):g}..{cfg.value_at(cm.DUTY_MAX):g} rpm.")
    period = 1.0 / args.freq
    duration = args.cycles * period
    peak_acc = abs(args.amp_rpm) * 2 * math.pi * args.freq
    print(f"duty = {mean_duty:.2f} +/- {abs(amp_duty):.2f} %, period {period:.2f} s, total {duration:.1f} s, "
          f"peak acceleration {peak_acc:.0f} rpm/s")
    if args.mode == "stream" and period < 20:
        print("WARNING: stream mode is too coarse/slow for periods below ~20 s.")
    ramp = min((r for r in (cfg.accel, cfg.decel) if r), default=None)
    if ramp is None:
        print("(ESCON speed ramp not found: cannot check the sine against it)")
    elif peak_acc > ramp:
        max_amp = ramp / (2 * math.pi * args.freq)
        print(f"WARNING: peak acceleration {peak_acc:.0f} rpm/s > ESCON speed ramp {ramp:g} rpm/s: the limiter will "
              f"flatten the sine. At {args.freq:g} Hz keep the amplitude below {max_amp:.0f} rpm, or raise the ramp "
              "in ESCON Studio (chipiron_patch_edc_v3.py --ramp).")
        if not args.yes and not cm.ask_yn("Continue anyway?", "n"):
            raise SystemExit("aborted.")
    else:
        print(f"Ramp check OK: peak acceleration {peak_acc:.0f} rpm/s < ESCON speed ramp {ramp:g} rpm/s.")

    link = cm.make_link(args.link, args.phone_ip, args.serial_port, record=not args.no_record, verbose=args.verbose)
    link.start()
    if args.song:
        subprocess.Popen(["play", SONG], text=True)
    time.sleep(2)
    try:
        if args.mode == "firmware":
            link.send(f"o{mean_duty:.2f},{abs(amp_duty):.2f},{args.freq:.4f}")
            time.sleep(duration)                      # the Arduino does the rest
        else:
            t0 = time.monotonic()
            k = 0
            while time.monotonic() - t0 < duration:
                t = time.monotonic() - t0
                duty = mean_duty + amp_duty * math.sin(2 * math.pi * args.freq * t)
                link.send(f"c{cm.fmt_duty(cm.q(duty))}")
                k += 1
                time.sleep(max(0.0, t0 + k * args.dt - time.monotonic()))
    except KeyboardInterrupt:
        print("\nInterrupted (Ctrl-C): going to rest.")
    finally:
        # Always return to ~0 rpm, even on Ctrl-C or an error.
        try:
            link.send("x" if args.mode == "firmware" else f"c{cm.fmt_duty(cfg.rest_duty)}")
        except Exception as exc:
            print("Rest command failed - check the motor manually!", exc)
        link.stop()


if __name__ == "__main__":
    main(gen_parser())