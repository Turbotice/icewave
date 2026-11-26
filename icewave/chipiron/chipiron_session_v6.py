"""
chipiron_session_v6.py -- one Phonefleet measurement session, two modes.  

  save_mode     no experiment. Starts at the ls step, then saves data:  ls -> save
  full_session  status -> Arduino mode -> (questions) -> 0 rpm + wait -> start -> experiment -> rest -> stop -> ls -> save

RUN
    py chipiron_session_v6.py                              (asks the mode)
    py chipiron_session_v6.py save_mode
    py chipiron_session_v6.py full_session                 (asks the Arduino mode, default from the .edc)
    py chipiron_session_v6.py full_session mo              (Arduino in oscillation mode, bipolar .edc)
    py chipiron_session_v6.py full_session mc              (Arduino in constant mode, unipolar .edc)
    py chipiron_session_v6.py full_session mo parameters_Sep26.edc   (.edc may be anywhere on the line)
    py chipiron_session_v6.py full_session skip            (do not send any mode: you know it is already set)
  Without an .edc argument the script asks, default = the .edc used last time.

SAMPLING PERIOD (both modes)
  Right after the phone IP the script asks the telemetry sampling period in ms = TELEM_PERIOD_MS of the
  PWM_full_v4.ino flashed on the Arduino (250 = 4 Hz old rate, 10 = 100 Hz, default = the one used last time).
  It is written as sampling_period_ms.txt in every folder saved (this session's folders always; older folders only if
  they have none), so Chipiron_plotting_v4.py / chipiron_decode_v4.py decode each acquisition with the right period.
  The name of the .edc used is saved as edc_used.txt in the same folders as the input log (first line = file name),
  so you know which file to give at the .edc prompt when plotting (full_session only).

WHAT full_session DOES
  1. loads the .edc ONCE into a MotorConfig (chipiron_helpers_v4) and hands it to every script it launches through
     the env var CHIPIRON_MOTOR_CFG (+ CHIPIRON_PHONE_IP): they never ask for the .edc again;
  2. checks the phone, then sends the Arduino mode "mc" / "mo" (PWM_full_v4.ino). The mode's rest duty (mc 10 %,
     mo 50 %) must give 0 rpm with the .edc, else the session stops with an ERROR. The Arduino resets (and boots
     DISABLED) every time the phone USB port is re-opened, so the mode is sent at each session; it also enables
     the motor. The Arduino does not answer: only the HTTP reply of the phone is shown;
  3. EXPERIMENT choices, all questions are asked BEFORE the acquisition starts:
       1) manual   : you run the experiment yourself; 'e' / 'd' enable / disable the motor, Enter stops
       2) oscillate: chipiron_oscillate_v4.py (sine; firmware mode needs Arduino mo + bipolar .edc, stream mode
                     works with anything). The sine is checked against the ESCON ramp and the duty range.
       3) run      : chipiron_run_v4.py. Step A (before the acquisition) builds the program with --build-only
                     (or you give an existing CSV, e.g. from chipiron_ramp_planner_v3.py); step B (during the
                     acquisition) replays it with --csv --yes --no-record.
       4) another script (--no-record is added automatically for the scripts of NO_RECORD_SCRIPTS)
       5) amp sweep: chipiron_amp_sweep_v4.py (sine at fixed frequency, amplitude STEPPED or CONTINUOUS; needs
                     Arduino mo + bipolar .edc; continuous needs PWM_full_v4.ino). The session first asks the
                     sweep type, then:
                       STEPPED    mean, frequency, rest between amplitudes, amplitudes (typed list or
                                  'range LO HI N [log]'), DURATION (c<cycles> per amplitude, s<seconds> per
                                  amplitude, or t<total seconds>).
                                  Answers e.g.:  s / 0 / 0.1 / 0 / 100 200 300 400 / c10
                                  or:            s / 0 / 0.2 / 10 / range 50 600 8 log / t1200
                       CONTINUOUS mean, frequency, amplitude LO HI [rpm], geometric (log) or linear envelope,
                                  duration (s<seconds> of the sweep, or c<cycles>). ONE 'r' command is sent.
                                  Answers e.g.:  c / 0 / 0.1 / 50 600 / y / s600
                     Then the script is run with --check-only BEFORE the acquisition (duty range, ESCON ramp,
                     total time printed; nothing is sent), and with --yes --no-record during the acquisition.
  4. sends the 0 rpm command (c<rest duty>) and waits ZERO_WAIT_S s, THEN calls /start;
  5. after the experiment (also Ctrl-C) the motor is sent back to rest, the acquisition is stopped, the files are
     listed and saved as ROOT/YYYY-MM-DD/HH-MM-SS/<original gobannos filename> (gps, accelerometer, gyroscope,
     magnetic_field and Usb files together).

DEPENDENCIES  (all in the same folder)
    connect_v3.py            phone ID -> IP (IP_PREFIX / NETWORK / DEFAULT_PHONELIST)
    run_gobannos_v3.py       HTTP calls to the phone (status, start, stop, ls, get-file, usb-cmd)
    escon_edc_v3.py          .edc reader
    chipiron_helpers_v4.py   MotorConfig, ENV names, rest duty, duty formatting
    chipiron_run_v4.py, chipiron_oscillate_v4.py, chipiron_amp_sweep_v4.py   launched as sub-processes (choices 2, 3, 5)
    Arduino firmware PWM_full_v4.ino plugged into the phone.
"""
import atexit
import math
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import date
from pathlib import Path

import connect_v3 as connect
import run_gobannos_v3 as gob
import chipiron_helpers_v4 as cm
from chipiron_helpers_v4 import ask
from escon_edc_v3 import load_edc
from chipiron_decode_v4 import TELEM_PERIOD_MS, MIN_PERIOD_MS, PERIOD_FILE, write_period_file

TS_RE = re.compile(r"(\d{4}-\d{2}-\d{2})T(\d{2})_(\d{2})_(\d{2})")
MODES = ("save_mode", "full_session")
ARDUINO_MODES = ("mc", "mo", "skip")
SCRIPT_DIR = Path(__file__).resolve().parent
OSCILLATE_SCRIPT = "chipiron_oscillate_v4.py"
RUN_SCRIPT = "chipiron_run_v4.py"
AMP_SWEEP_SCRIPT = "chipiron_amp_sweep_v4.py"
# Scripts that accept --no-record (the phone is already recording, they must not call /start and /stop):
NO_RECORD_SCRIPTS = {OSCILLATE_SCRIPT, RUN_SCRIPT, AMP_SWEEP_SCRIPT}
PROGRAM_FOLDER = SCRIPT_DIR / "programs"          # programs built by choice 3 (--build-only)
LAST_IP_FILE = SCRIPT_DIR / ".phonefleet_last_ip.txt"   # remembers the phone IP used last
LAST_PERIOD_FILE = SCRIPT_DIR / ".phonefleet_last_period_ms.txt"   # remembers the telemetry period used last
PERIOD_MS = TELEM_PERIOD_MS                      # telemetry period [ms] of the flashed firmware, set by ask_sampling_period()
ZERO_WAIT_S = 30                                  # seconds to wait at 0 rpm before /start (tunable)
FIRMWARE_REST_DUTY = {"mc": 10.0, "mo": 50.0}     # rest duty of each Arduino mode (PWM_full_v4.ino)

INPUT_LOG_FOLDER = SCRIPT_DIR / "input_logs"      # commands sent are logged here during the session, then copied
                                                  # next to the phone data (input_<date>_<time>.csv) at the save step
INPUT_LOG = None                                  # Path of this session's log, set by setup_input_log()
FILES_BEFORE = None                               # files already on the phone before /start (to find this session's files)

CFG = None                                        # MotorConfig, filled by load_motor_config()
EDC_PATH = None                                   # Path of the .edc loaded for this session (its name is saved with the data)
TEMP_FILES = []                                   # temp JSON files to delete at exit


# ---------------------------------------------------------------- helpers
def fatal(msg):
    print(f"\nERROR: {msg}")
    sys.exit(1)


def _cleanup():
    for f in TEMP_FILES:
        try:
            os.remove(f)
        except OSError:
            pass


atexit.register(_cleanup)


def parse_name(filename):
    """Return (YYYY-MM-DD, HH-MM-SS) from a gobannos filename, or (None, None)."""
    m = TS_RE.search(os.path.basename(filename))
    if not m:
        return None, None
    d, h, mi, s = m.groups()
    return d, f"{h}-{mi}-{s}"


def norm_time(raw):
    """'15:11:58', '15_11_58' or '15-11-58' -> '15-11-58'."""
    return re.sub(r"[:_]", "-", raw.strip())


def status_str(ip):
    st = gob.get_status(ip)
    return st.decode() if isinstance(st, bytes) else str(st)


def wait_for_status(ip, target, timeout=15):
    t0 = time.time()
    while time.time() - t0 < timeout:
        if status_str(ip) == target:
            return True
        time.sleep(0.5)
    return False


def send_usb_cmd(ip, cmd, note=""):
    """Send <cmd> to the Arduino through the phone (same call as 'usb <command>' in mobile_v3.py).
    Raises RuntimeError if the phone cannot be reached."""
    print(f"  gob.usb_cmd({ip!r}, {cmd!r})")
    try:
        out = gob.usb_cmd(ip, cmd)
    except Exception:
        cm.log_command(cmd, ok=False, source="session", note=note)
        raise
    cm.log_command(cmd, ok=True, source="session", note=note)
    if out:
        print(f"  phone replied: {out[:100]}")
    return out


def ask_phone_ip():
    """Ask for the full phone IP, e.g. 10.48.96.2 (default = last one used, remembered in LAST_IP_FILE).
    Nothing is built from a phone ID / network number. Returns (ip, last octet for display)."""
    try:
        default = LAST_IP_FILE.read_text(encoding="utf-8").strip()
    except OSError:
        default = ""
    while True:
        raw = ask("Phone IP address (e.g. 10.48.96.2)", default).strip()
        raw = re.sub(r"^https?://", "", raw).split(":")[0]
        if re.fullmatch(r"\d{1,3}(\.\d{1,3}){3}", raw) and all(int(x) <= 255 for x in raw.split(".")):
            break
        print("  -> please type the full IPv4 address, like 10.48.96.2")
    try:
        LAST_IP_FILE.write_text(raw, encoding="utf-8")
    except OSError:
        pass
    return raw, raw.split(".")[-1]


def ask_sampling_period():
    """Ask the telemetry sampling period [ms] = TELEM_PERIOD_MS of the flashed PWM_full_v4.ino (250 = old 4 Hz,
    10 = 100 Hz). Default = last one used (LAST_PERIOD_FILE), else chipiron_decode_v4.TELEM_PERIOD_MS."""
    global PERIOD_MS
    try:
        default = float(LAST_PERIOD_FILE.read_text(encoding="utf-8").strip())
    except (OSError, ValueError):
        default = TELEM_PERIOD_MS
    while True:
        raw = ask("Telemetry sampling period [ms] (TELEM_PERIOD_MS in the flashed PWM_full_v4.ino; 250 = old firmware)",
                  f"{default:g}")
        try:
            p = float(raw.replace(",", "."))
            if p >= MIN_PERIOD_MS:
                break
        except ValueError:
            pass
        print(f"  -> type a period in ms, at least {MIN_PERIOD_MS:g} (e.g. 10 for 100 Hz, 250 for 4 Hz)")
    PERIOD_MS = p
    try:
        LAST_PERIOD_FILE.write_text(f"{p:g}", encoding="utf-8")
    except OSError:
        pass
    print(f"Sampling period {p:g} ms = {1000 / p:.3g} Hz (written next to the saved data as {PERIOD_FILE}).")
    if p < 250:
        print("  (reminder: the Arduino must run the v4 firmware: older firmware sends 250 ms frames without start flag)")
    return p


def choose_mode():
    """Mode from argv if given, otherwise ask."""
    if len(sys.argv) > 1:
        arg = sys.argv[1].strip().lower()
        if arg in MODES:
            return arg
        print(f"Unknown mode '{sys.argv[1]}'. Choose one of: {', '.join(MODES)}")
    print("Mode:")
    print("  1) save_mode     - no experiment, just list and save files from the phone")
    print("  2) full_session  - status, start, experiment, stop, then save files")
    return {"1": "save_mode", "2": "full_session"}.get(ask("Choice", "1"), "save_mode")


# ---------------------------------------------------------------- ESCON .edc -> MotorConfig
def pop_edc_arg():
    """Remove and return the first command-line argument ending in .edc (None if absent)."""
    for i, a in enumerate(sys.argv[1:], start=1):
        if a.lower().endswith(".edc"):
            return sys.argv.pop(i)
    return None


def load_motor_config(arg=None):
    """Load the .edc once into the global CFG and export it (env) for the sub-processes."""
    global CFG, EDC_PATH
    path = (arg or ask("ESCON .edc parameter file loaded in ESCON Studio", cm._default_edc())).strip('"').strip("'")
    if not path:
        fatal("no .edc file given (needed to know the rpm <-> duty mapping).")
    path = Path(path).expanduser()
    EDC_PATH = path.resolve()
    try:
        edc = load_edc(path)
        CFG = cm.MotorConfig.from_edc_object(edc, path)
    except Exception as exc:
        fatal(f"cannot use {path}: {exc}")
    try:
        cm.LAST_EDC_FILE.write_text(str(path.resolve()), encoding="utf-8")
    except OSError:
        pass
    print(edc.summary())
    print(CFG.summary())
    if not CFG.is_speed:
        fatal("the .edc does not define a PWM speed set value (needs a speed operation mode and a non-zero rpm/% "
              "scaling on the PWM input).")
    if not cm.DUTY_MIN - 0.01 <= CFG.zero_duty <= cm.DUTY_MAX + 0.01:
        fatal(f"0 rpm would be at {CFG.zero_duty:.1f} % duty, outside the valid PWM range "
              f"{cm.DUTY_MIN:g}..{cm.DUTY_MAX:g} %.")
    if CFG.kind == "other":
        print("WARNING: neither the standard unipolar (0 rpm at 10 %) nor bipolar mapping: Arduino modes mc/mo "
              "will not match (use 'skip' and check the rest duty yourself).")
    TEMP_FILES.append(CFG.dump_temp())
    os.environ[cm.ENV_CFG] = TEMP_FILES[-1]


def zero_command():
    """Arduino command giving 0 rpm with the loaded .edc, e.g. 'c50.0' (bipolar) or 'c10.0' (unipolar)."""
    return f"c{cm.fmt_duty(CFG.rest_duty)}"


def mode_problem(arduino_mode):
    """Error text if the Arduino mode's rest duty does not give 0 rpm with the .edc, else None."""
    rest = FIRMWARE_REST_DUTY.get(arduino_mode)
    if rest is None or abs(CFG.zero_duty - rest) <= 0.5:
        return None
    return (f"Arduino mode '{arduino_mode}' rests at {rest:g} % duty, but the .edc ({CFG.kind} mapping) "
            f"gives 0 rpm at {CFG.zero_duty:.1f} %. Use mode "
            f"'{'mo' if CFG.kind == 'bipolar' else 'mc'}' with this .edc, or load the matching .edc in ESCON Studio.")


# ---------------------------------------------------------------- steps
def step_status(ip, phone):
    print(f"\n[1] Status of phone {phone} ({ip})")
    try:
        st = status_str(ip)
    except Exception as exc:
        print(f"Cannot reach phone: {exc}")
        print("Check Wi-Fi / network number / phone ID in connect_v3.py.")
        sys.exit(1)
    print(f"Status: {st}")
    return st


def step_arduino_mode(ip, phone, arduino_mode=None):
    """Send mc / mo (also enables the motor) after checking it matches the .edc. Returns the mode used."""
    print("\n[1b] Arduino mode (PWM_full_v4.ino)")
    print("  mc   = constant mode     (rest duty 10 %, unipolar ESCON config)")
    print("  mo   = oscillation mode  (rest duty 50 %, bipolar ESCON config)")
    print("  skip = do not send anything (the Arduino boots with the motor DISABLED: send 'e' yourself)")
    print("  Make sure the ESCON configuration matches the mode, and the motor is stopped.")
    if arduino_mode not in ARDUINO_MODES:
        arduino_mode = ask("Arduino mode (mc / mo / skip)", "mo" if CFG.kind == "bipolar" else "mc").lower()
    if arduino_mode not in ARDUINO_MODES:
        fatal(f"unknown Arduino mode '{arduino_mode}'.")
    if arduino_mode == "skip":
        print("Arduino mode left as is. WARNING: if the Arduino was reset, the motor stays disabled until 'e'.")
        return "skip"
    problem = mode_problem(arduino_mode)
    if problem:
        fatal(problem)
    try:
        send_usb_cmd(ip, arduino_mode, note="arduino mode")
        time.sleep(0.5)
        print(f"Mode command '{arduino_mode}' sent (motor armed at the rest duty).")
    except Exception as exc:
        print(f"Could not send '{arduino_mode}': {exc}")
        print("Check the Arduino is plugged into the phone and the phone server is running.")
        if ask("Continue anyway? (y/n)", "n").lower() != "y":
            sys.exit(1)
    return arduino_mode


def find_script(name):
    """Look for a script next to this file, then in the current folder."""
    p = Path(name)
    if p.is_file():
        return p.resolve()
    for base in (SCRIPT_DIR, Path.cwd()):
        cand = base / name
        if cand.is_file():
            return cand.resolve()
    return None


def _need_script(name):
    script = find_script(name)
    if script is None:
        print(f"Cannot find {name} (looked next to {Path(__file__).name} and in {Path.cwd()}).")
        return None
    print(f"Using {script}")
    return script


def build_command(choice, arduino_mode):
    """Ask EVERYTHING needed (before the acquisition) and return the command list to run later,
    None for 'I run it myself', or 'abort'."""
    if choice == "1":
        return None

    if choice == "2":                                   # ---- oscillation
        script = _need_script(OSCILLATE_SCRIPT)
        if script is None:
            return "abort"
        default_mode = "firmware" if (arduino_mode == "mo" and CFG.kind == "bipolar") else "stream"
        osc_mode = ask("Oscillation mode: firmware (Arduino builds the sine, needs mo) / stream (c every dt s)",
                       default_mode).lower()
        if osc_mode not in ("firmware", "stream"):
            print("Unknown oscillation mode.")
            return "abort"
        if osc_mode == "firmware" and arduino_mode == "mc":
            print("\nERROR: firmware mode sends the 'o' command, which the Arduino only accepts in mode mo "
                  "(it is in mc). Restart the session with mo, or use stream mode.")
            return "abort"
        if osc_mode == "firmware" and CFG.kind != "bipolar":
            print(f"\nERROR: firmware mode (mo) needs the bipolar mapping but the .edc is {CFG.kind}. "
                  "Use stream mode, or load a bipolar .edc.")
            return "abort"
        if osc_mode == "firmware" and arduino_mode == "skip":
            print("  WARNING: Arduino mode was skipped: make sure it is in mo, otherwise the sine is ignored.")
        mean = ask("Mean speed (rpm)", "0" if CFG.kind == "bipolar" else "400")
        amp = ask("Amplitude (rpm)", "100")
        freq = ask("Frequency (Hz)", "0.2")
        cycles = ask("Number of cycles", "5")
        try:
            m, a, f, n = float(mean), float(amp), float(freq), float(cycles)
            if f <= 0 or n <= 0:
                raise ValueError("frequency and number of cycles must be > 0")
            lo, hi = CFG.duty_for(m) - abs(a / CFG.slope), CFG.duty_for(m) + abs(a / CFG.slope)
            if lo < cm.DUTY_MIN or hi > cm.DUTY_MAX:
                raise ValueError(f"duty would be {lo:.1f}..{hi:.1f} %, outside {cm.DUTY_MIN:g}..{cm.DUTY_MAX:g} %")
            print(f"  Expected duration of the oscillation: {n / f:.1f} s")
            peak_acc = abs(a) * 2 * math.pi * f
            ramp = min((r for r in (CFG.accel, CFG.decel) if r), default=None)
            if ramp is None:
                print("  (ESCON speed ramp not found in the .edc: cannot check the sine against it)")
            elif peak_acc > ramp:
                print(f"  WARNING: peak acceleration of the sine = {peak_acc:.0f} rpm/s > ESCON speed ramp "
                      f"{ramp:g} rpm/s: the ramp will flatten the sine. For {f:g} Hz keep the amplitude "
                      f"below {ramp / (2 * math.pi * f):.0f} rpm, or raise the ramp in ESCON Studio.")
                if ask("  Continue anyway? (y/n)", "n").lower() != "y":
                    return "abort"
            else:
                print(f"  Ramp check OK: peak acceleration {peak_acc:.0f} rpm/s < ESCON speed ramp {ramp:g} rpm/s.")
        except ValueError as exc:
            print(f"\nERROR: {exc}")
            return "abort"
        extra = ask("Extra arguments (optional, e.g. --song --dt 0.5)")
        args = ["--mode", osc_mode, "--mean-rpm", mean, "--amp-rpm", amp, "--freq", freq, "--cycles", cycles,
                "--yes"] + extra.split()

    elif choice == "3":                                 # ---- run program (two steps)
        script = _need_script(RUN_SCRIPT)
        if script is None:
            return "abort"
        existing = ask("Existing program / ramp-planner CSV to run (empty = build one now)").strip('"').strip("'")
        if existing:
            if not Path(existing).is_file():
                print(f"\nERROR: {existing} not found.")
                return "abort"
            program = Path(existing).resolve()
        else:
            PROGRAM_FOLDER.mkdir(exist_ok=True)
            program = PROGRAM_FOLDER / f"session_{time.strftime('%Y%m%d-%H%M%S')}.csv"
            print("  Step A: building the program now (nothing is sent to the motor, the acquisition is not started).")
            res = subprocess.run([sys.executable, str(script), "--build-only", str(program)], check=False)
            if res.returncode != 0 or not program.is_file():
                print(f"\nERROR: program not built (return code {res.returncode}).")
                return "abort"
        cmd = [sys.executable, str(script), "--csv", str(program), "--yes"]
        print(f"  Step B (during the acquisition) will run: {' '.join(cmd[1:])}")
        return cmd + ["--no-record"]

    elif choice == "5":                                 # ---- amplitude sweep (stepped or continuous)
        script = _need_script(AMP_SWEEP_SCRIPT)
        if script is None:
            return "abort"
        if arduino_mode == "mc":
            print("\nERROR: the amplitude sweep sends the 'o' / 'r' command, which the Arduino only accepts in mode mo "
                  "(it is in mc). Restart the session with mo.")
            return "abort"
        if CFG.kind != "bipolar":
            print(f"\nERROR: the amplitude sweep (mo) needs the bipolar mapping but the .edc is {CFG.kind}.")
            return "abort"
        if arduino_mode == "skip":
            print("  WARNING: Arduino mode was skipped: make sure it is in mo, otherwise the sine is ignored.")
        kind = ask("Sweep type: [s]tepped (plateaus of constant amplitude) / [c]ontinuous "
                   "(amplitude glides smoothly, ONE 'r' command, needs PWM_full_v4.ino)", "c").strip().lower()
        continuous = kind.startswith("c")
        mean = ask("Mean speed (rpm)", "0")
        freq = ask("Frequency (Hz)", "0.1")
        try:
            float(mean), float(freq)
            if continuous:
                lohi = ask("Amplitude from LO to HI [rpm] (two numbers)", "50 600").replace(",", " ").split()
                if len(lohi) != 2:
                    raise ValueError("give exactly two numbers: LO HI")
                float(lohi[0]), float(lohi[1])
                geo = ask("Geometric (log) envelope? (y/n; n = linear; log needs LO, HI > 0)", "y").lower().startswith("y")
                dur = ask("Duration: s<seconds> of the sweep, or c<cycles>", "s600").strip().lower()
                flag = {"s": "--sweep-time", "c": "--cycles"}.get(dur[:1])
                if flag is None:
                    raise ValueError("duration must start with s or c (e.g. s600, c60)")
                float(dur[1:])
                amp_args = ["--sweep", lohi[0], lohi[1]] + (["--log"] if geo else [])
                rest = None
            else:
                rest = ask("Rest at 0 rpm between amplitudes [s]", "0")
                float(rest)
                spec = ask("Amplitudes: typed list in rpm (100 200 300) or 'range LO HI N [log]'",
                           "100 200 300 400").split()
                dur = ask("Duration: c<cycles> per amplitude, s<seconds> per amplitude, or t<total seconds>",
                          "c5").strip().lower()
                if spec and spec[0].lower() == "range":
                    nums = [float(x) for x in spec[1:4]]
                    if len(nums) != 3:
                        raise ValueError("range needs LO HI N")
                    amp_args = ["--amp-range"] + spec[1:4] + (["--log"] if "log" in [x.lower() for x in spec[4:]] else [])
                else:
                    [float(x) for x in spec]
                    if not spec:
                        raise ValueError("no amplitude given")
                    amp_args = ["--amps"] + spec
                flag = {"c": "--cycles", "s": "--seg-time", "t": "--total-time"}.get(dur[:1])
                if flag is None:
                    raise ValueError("duration must start with c, s or t (e.g. c10, s60, t1200)")
                float(dur[1:])
        except ValueError as exc:
            print(f"\nERROR: {exc}")
            return "abort"
        extra = ask("Extra arguments (optional, e.g. --verbose)")
        args = ["--mean-rpm", mean, "--freq", freq] + amp_args + [flag, dur[1:]]
        if rest is not None:
            args += ["--rest-between", rest]
        args += extra.split()
        if continuous:
            print("  NOTE: the continuous sweep sends the 'r' command: the Arduino must run PWM_full_v4.ino "
                  "(older firmware ignores it and the motor stays at rest).")
        print("  Checking the sweep (nothing is sent, the acquisition is not started)...")
        res = subprocess.run([sys.executable, str(script)] + args + ["--check-only"], check=False)
        if res.returncode != 0:
            print(f"\nERROR: sweep check failed or refused (return code {res.returncode}).")
            return "abort"
        args.append("--yes")

    else:                                               # ---- another script
        name = ask("Script to run")
        if not name:
            return None
        script = _need_script(name)
        if script is None:
            return "abort"
        args = ask("Arguments (as you would type them after the script name)").split()

    cmd = [sys.executable, str(script)] + args
    if script.name in NO_RECORD_SCRIPTS and "--no-record" not in cmd:
        cmd.append("--no-record")
    return cmd


def step_zero_rpm(ip):
    """Send the 0 rpm command of the loaded .edc, then wait ZERO_WAIT_S s before /start."""
    cmd = zero_command()
    print(f"\n[2a] Motor to 0 rpm: sending '{cmd}' (0 rpm at {CFG.zero_duty:.2f} % duty, {CFG.kind} mapping)")
    try:
        send_usb_cmd(ip, cmd, note="0 rpm before start")
    except Exception as exc:
        print(f"Could not send '{cmd}': {exc}")
        if ask("Continue anyway? (y/n)", "n").lower() != "y":
            sys.exit(1)
    print(f"Waiting {ZERO_WAIT_S} s at 0 rpm before starting the acquisition...")
    t0 = time.time()
    while time.time() - t0 < ZERO_WAIT_S:
        time.sleep(min(10, max(0.0, ZERO_WAIT_S - (time.time() - t0))))
        left = ZERO_WAIT_S - (time.time() - t0)
        if left > 1:
            print(f"  {left:.0f} s left")


def manual_loop(ip):
    """Manual experiment: 'e' enables the motor, 'd' disables it, Enter stops. Other Arduino commands too."""
    print("Experiment running... commands: 'e' = enable motor, 'd' = disable motor, 'x' = rest, "
          "'c<duty>' (e.g. c60.0), 'o<mean>,<amp>,<freq>' and 'r<mean>,<amp0>,<amp1>,<freq>,<T>[,<shape>]' (mode mo), Enter = STOP the acquisition.")
    while True:
        key = input("> ").strip().lower()
        if key == "":
            return
        if key[0] in ("e", "d", "x", "c", "o", "r"):
            try:
                send_usb_cmd(ip, key, note="manual")
                print(f"command '{key}' sent.")
            except Exception as exc:
                print(f"Could not send '{key}': {exc}")
        else:
            print("Unknown key: use e, d, x, c<duty>, o<mean>,<amp>,<freq>, r<mean>,<amp0>,<amp1>,<freq>,<T>[,<shape>] or Enter.")


def step_rest(ip):
    """Best effort: motor back to rest after the experiment, whatever happened."""
    try:
        send_usb_cmd(ip, zero_command(), note="rest after experiment")
        print(f"Motor sent back to rest ({zero_command()}).")
    except Exception as exc:
        print(f"!! could not send the rest command: {exc}  -> CHECK / STOP THE MOTOR MANUALLY")


def step_start_experiment(ip, phone, arduino_mode=None):
    print("\n[2] Start experiment")
    print("  1) Start acquisition, then I run the experiment myself (press Enter to stop)")
    print(f"  2) Start acquisition, then run {OSCILLATE_SCRIPT} (sine command)")
    print(f"  3) Start acquisition, then run {RUN_SCRIPT} (holds / ramps / staircases / CSV)")
    print("  4) Start acquisition, then run another script")
    print(f"  5) Start acquisition, then run {AMP_SWEEP_SCRIPT} (sine, amplitude stepped OR continuous sweep; needs mo + bipolar)")
    choice = ask("Choice", "2" if CFG.kind == "bipolar" else "3")
    if choice not in ("1", "2", "3", "4", "5"):
        print("Unknown choice.")
        sys.exit(1)

    cmd = build_command(choice, arduino_mode)     # all questions BEFORE the acquisition starts
    if cmd == "abort":
        print("Nothing started.")
        sys.exit(1)

    step_zero_rpm(ip)                             # 0 rpm, wait, THEN start the acquisition
    gob.individual_start(ip)
    cm.log_command("[acquisition /start]", source="session")
    time.sleep(1)
    print(f"Status after start: {status_str(ip)}")
    try:
        if cmd is None:
            manual_loop(ip)
        else:
            print(f"Running: {' '.join(cmd)}")
            res = subprocess.run(cmd, check=False)
            print(f"Script finished with return code {res.returncode}")
    except KeyboardInterrupt:
        print("\nInterrupted (Ctrl-C).")
    finally:
        step_rest(ip)


def step_stop(ip, phone):
    print("\n[3] Stop acquisition")
    if status_str(ip) != "STOPPED":
        gob.individual_stop(ip)
    cm.log_command("[acquisition /stop]", source="session")
    if wait_for_status(ip, "STOPPED"):
        print("Phone STOPPED.")
    else:
        print(f"Warning: phone status is still {status_str(ip)}")
    time.sleep(1)


def step_ls(ip):
    print("\n[4] Files on phone")
    try:
        files = gob.get_file_list(ip)
    except Exception as exc:
        print(f"Cannot list files: {exc}")
        print("Check Wi-Fi / network number / phone ID in connect_v3.py.")
        sys.exit(1)
    for i, f in enumerate(files):
        print(i, f)
    print(f"({len(files)} files)")
    return files


def setup_input_log():
    """New log file for this session; its path is handed to run / oscillate through the environment."""
    global INPUT_LOG
    INPUT_LOG_FOLDER.mkdir(exist_ok=True)
    INPUT_LOG = INPUT_LOG_FOLDER / f"input_session_{time.strftime('%Y%m%d-%H%M%S')}.csv"
    os.environ[cm.ENV_LOG] = str(INPUT_LOG)
    print(f"Commands sent will be logged in {INPUT_LOG}")


def copy_input_log(targets, saved_folders):
    """Copy the session's input log as input_<date>_<time>.csv next to the phone files of THIS session.
    targets = {folder: (date, time)}. If this session's files could not be identified (phone unreachable
    before /start) and exactly one folder was saved, use that one."""
    if INPUT_LOG is None or not INPUT_LOG.is_file():
        return
    if not targets and FILES_BEFORE is None and len(saved_folders) == 1:
        folder = next(iter(saved_folders))
        targets = {folder: (folder.parent.name, folder.name)}
    if not targets:
        print(f"Input log not copied (no file of this session was saved). It stays in {INPUT_LOG}")
    for folder, (d, t) in targets.items():
        dest = folder / f"input_{d or 'unknown_date'}_{t or 'unknown_time'}.csv"
        shutil.copy2(INPUT_LOG, dest)               # same session = same content: overwriting is harmless
        print(f"  saved {dest}  (commands sent)")


EDC_FILE = "edc_used.txt"


def write_edc_name(targets, saved_folders):
    """Write the name of the .edc used in this session as edc_used.txt next to the phone files of THIS session
    (same folders as the input log). First line = file name, second line = full path (comment)."""
    if EDC_PATH is None:
        return
    if not targets and FILES_BEFORE is None and len(saved_folders) == 1:
        folder = next(iter(saved_folders))
        targets = {folder: (folder.parent.name, folder.name)}
    for folder in targets:
        dest = folder / EDC_FILE
        try:
            dest.write_text(f"{EDC_PATH.name}\n# full path: {EDC_PATH}\n", encoding="utf-8")
            print(f"  saved {dest}  ({EDC_PATH.name})")
        except OSError as exc:
            print(f"  could not write {EDC_FILE} in {folder}: {exc}")


def save_files(ip, phone, files, indexes, root):
    """Download files[i] for i in indexes into root/date/time/."""
    root = Path(root)
    saved_folders = set()
    targets = {}                                 # folders holding files recorded during THIS session
    for i in indexes:
        name = files[i]
        d, t = parse_name(name)
        folder = root / (d or "unknown_date") / (t or "unknown_time")
        folder.mkdir(parents=True, exist_ok=True)
        try:
            raw = gob.get_file(ip, name)
            data = raw if isinstance(raw, bytes) else str(raw).encode("utf-8")
            dest = folder / os.path.basename(name).lstrip("-")
            dest.write_bytes(data)
            saved_folders.add(folder)
            if FILES_BEFORE is not None and name not in FILES_BEFORE:
                targets[folder] = (d, t)
            print(f"  saved {dest}")
        except Exception as exc:
            print(f"  could not retrieve {name}: {exc}")
    for folder in sorted(saved_folders):             # sampling period next to the data (always for this session's folders)
        if folder in targets or not (folder / PERIOD_FILE).is_file():
            try:
                write_period_file(folder, PERIOD_MS)
                print(f"  saved {folder / PERIOD_FILE}  ({PERIOD_MS:g} ms)")
            except (OSError, ValueError) as exc:
                print(f"  could not write {PERIOD_FILE} in {folder}: {exc}")
    copy_input_log(targets, saved_folders)
    write_edc_name(targets, saved_folders)
    print(f"\nDone. {len(saved_folders)} folder(s) written under {root.resolve()}")


def step_save(ip, phone, files):
    if not files:
        print("No files on phone, nothing to save.")
        return
    print("\n[5] Save files")
    print("  1) one file (by index)")
    print("  2) all files on the phone")
    print("  3) all files of a given date")
    print("  4) all files of a given date and time")
    choice = ask("Choice", "4")
    root = ask("Root folder on laptop", "Phonefleet_data")

    today = date.today().isoformat()
    indexes = []

    if choice == "1":
        raw = ask("File index (from the list above)")
        try:
            idx = int(raw)
        except ValueError:
            print("Not an integer.")
            return
        if not 0 <= idx < len(files):
            print(f"Index must be between 0 and {len(files) - 1}.")
            return
        indexes = [idx]
    elif choice == "2":
        indexes = list(range(len(files)))
    elif choice == "3":
        d = ask("Date YYYY-MM-DD", today)
        indexes = [i for i, f in enumerate(files) if parse_name(f)[0] == d]
    elif choice == "4":
        d = ask("Date YYYY-MM-DD", today)
        t = norm_time(ask("Time HH:MM:SS (e.g. 15:11:58)"))
        indexes = [i for i, f in enumerate(files) if parse_name(f) == (d, t)]
    else:
        print("Unknown choice.")
        return

    if not indexes:
        print("No matching files.")
        return
    print(f"{len(indexes)} file(s) selected.")
    save_files(ip, phone, files, indexes, root)


# ---------------------------------------------------------------- main
def main():
    edc_arg = pop_edc_arg()          # before choose_mode(), which reads argv[1] / argv[2]
    mode = choose_mode()
    print(f"\n>>> Mode: {mode}")

    ip, phone = ask_phone_ip()
    print(f"Phone {ip}")
    ask_sampling_period()

    if mode == "full_session":
        load_motor_config(edc_arg)
        os.environ[cm.ENV_IP] = ip                   # every sub-script talks to THIS phone
        setup_input_log()
        step_status(ip, phone)
        arduino_mode = sys.argv[2].strip().lower() if len(sys.argv) > 2 else None
        arduino_mode = step_arduino_mode(ip, phone, arduino_mode)
        global FILES_BEFORE
        try:
            FILES_BEFORE = set(gob.get_file_list(ip))   # to recognise the files of this session later
        except Exception as exc:
            print(f"Could not list the files before the start ({exc}): the input log will only be copied if one folder is saved.")
        step_start_experiment(ip, phone, arduino_mode)
        step_stop(ip, phone)

    files = step_ls(ip)
    step_save(ip, phone, files)


if __name__ == "__main__":
    main()