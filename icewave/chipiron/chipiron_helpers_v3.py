"""
chipiron_helpers_v3.py -- shared helpers of the Chipiron motor suite (v3). 

WHAT IT PROVIDES
  * MotorConfig        PWM duty (%) <-> speed (rpm) mapping, rest duty, ESCON speed-ramp limits, motor type...
                       all read from the ESCON .edc file (the one loaded in ESCON Studio).
  * get_config()       standalone  -> asks for the .edc and loads it (through escon_edc_v3.load_edc)
                       from session -> chipiron_session_v3 exports the config it already loaded in a small
                                       JSON file whose path is in the env var CHIPIRON_MOTOR_CFG; the .edc is
                                       then NOT read a second time (no question asked during an acquisition).
  * HttpLink/LoraLink  how a command text ("c50", "o50,20,0.2", "x"...) reaches the Arduino:
                         http : laptop -> wifi -> phone HTTP server (port 8080) -> USB -> Arduino
                         lora : laptop -> USB -> LoRa bridge Arduino ~~radio~~ motor Arduino (needs pyserial)
  * check_timeline()   simulates the ESCON speed-ramp limiter on a list of rows and warns when a plateau is
                       not reached or a ramp is steeper than the limiter allows
  * write_plan_csv() / read_plan_csv()   the CSV exchanged between the ramp planner and the run script
  * hw_duty()          the duty cycle really produced by the Arduino (8-bit analogWrite, 100/255 % steps)

A "row" is a dict: {"dur": seconds, "duty": PWM duty cycle in %, "kind": "hold"|"ramp", "seg": label}.

DEPENDENCIES
  stdlib only.  Lazy imports: escon_edc_v3 (when an .edc must be read), pyserial (LoRa link only).

USED BY
  chipiron_run_v3.py, chipiron_oscillate_v3.py, chipiron_ramp_planner_v3.py, chipiron_session_v3.py

ENVIRONMENT VARIABLES (set by chipiron_session_v3.py when it launches another script)
  CHIPIRON_MOTOR_CFG   path of the JSON motor configuration to reuse (skips the .edc question)
  CHIPIRON_PHONE_IP    phone IP to use (so every script talks to the same phone as the session)

EXAMPLES (from Python)
  import chipiron_helpers_v3 as cm
  cfg = cm.get_config("parameters_Sep26.edc")      # prints the summary, returns a MotorConfig
  cfg.duty_for(400)                                # duty (%) that gives 400 rpm
  cfg.value_at(75)                                 # rpm obtained at 75 % duty
  cfg.rest_duty                                    # duty that gives 0 rpm (50.0 for the bipolar .edc)
  link = cm.make_link("http", ip="10.153.247.36", record=False)
  link.send("c60.0")                               # 60 % duty
  link.send("o50.00,20.00,0.2000")                 # sine 50 +/- 20 % at 0.2 Hz (firmware mode mo only)
  link.send("x")                                   # back to the rest duty
"""
import csv
import json
import math
import os
import re
import tempfile
import time
import urllib.request
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import quote

# ------------------------------------------------------------------ constants
ENV_CFG = "CHIPIRON_MOTOR_CFG"              # JSON motor config handed over by the session
ENV_IP = "CHIPIRON_PHONE_IP"                # phone IP handed over by the session
PHONE_PORT = 8080                           # Gobannos HTTP port (same as connect_v3.PORT)
DEFAULT_PHONE_IP = "10.153.247.36"          # used when neither --phone-ip nor CHIPIRON_PHONE_IP is given
#                                             (= connect_v3 with IP_PREFIX "10", network 247, phone 36)
SCRIPT_DIR = Path(__file__).resolve().parent
LAST_EDC_FILE = SCRIPT_DIR / ".phonefleet_last_edc.txt"   # remembers the .edc used last (shared with the session)
DUTY_MIN, DUTY_MAX = 10.0, 90.0             # valid ESCON PWM input range (%); outside = "PWM out of range" error
DUTY_DECIMALS = 1                           # decimals sent in "c<duty>". PWM_full_v3.ino parses decimals (atof),
#                                             but its real resolution is 100/255 = 0.39 % (see hw_duty()).
#                                             Set 0 only if an old firmware using atoi() (PWM_control_v3) is loaded.
PWM_COUNTS = 255                            # Arduino analogWrite resolution (8 bit)
MOTOR_TYPES = {0: "maxon EC motor (brushless)", 1: "maxon DC motor (brushed)"}   # enum EMotorType of the .edc
OP_MODES = {1: "current controller", 2: "speed controller (closed loop)", 3: "speed controller (open loop)"}


# ------------------------------------------------------------------ small input helpers (command-line prompts)
def ask(prompt, default=None):
    """input() with a [default]; Enter returns the default."""
    suffix = f" [{default}]" if default not in (None, "") else ""
    raw = input(f"{prompt}{suffix}: ").strip()
    return raw if raw else (default if default is not None else "")


def ask_float(prompt, default=None):
    """Ask until a number is typed (comma or dot accepted)."""
    while True:
        raw = ask(prompt, None if default is None else f"{default:g}")
        try:
            return float(raw.replace(",", "."))
        except ValueError:
            print("  -> please type a number.")


def ask_int(prompt, default=None):
    """Ask until an integer is typed."""
    while True:
        raw = ask(prompt, None if default is None else str(default))
        try:
            return int(raw)
        except ValueError:
            print("  -> please type an integer.")


def ask_yn(prompt, default="y"):
    """Yes/no question; returns True for anything starting with y."""
    return ask(f"{prompt} (y/n)", default).lower().startswith("y")


def q(duty):
    """Quantise a duty cycle to the number of decimals that is really sent to the Arduino."""
    return round(float(duty), DUTY_DECIMALS)


def fmt_duty(duty):
    """Text of a duty for the 'c' command: '50' (DUTY_DECIMALS = 0) or '50.0' (DUTY_DECIMALS = 1)."""
    return str(int(round(duty))) if DUTY_DECIMALS == 0 else f"{duty:.{DUTY_DECIMALS}f}"


def hw_duty(duty):
    """Duty cycle really produced by the Arduino: analogWrite(PWM, round(duty*255/100)) -> steps of 0.39 %."""
    return int(float(duty) * PWM_COUNTS / 100 + 0.5) * 100.0 / PWM_COUNTS


# ------------------------------------------------------------------ motor configuration
@dataclass
class MotorConfig:
    """Everything the scripts need to know about the ESCON set-up (read from the .edc)."""
    source: str            # "edc" | "session" | "manual"
    edc_name: str
    quantity: str          # "speed" | "current"
    unit: str              # "rpm" | "A" | "?"
    offset: float          # set value at 0 % duty
    slope: float           # set value per 1 % duty
    accel: float = None    # ESCON speed ramp acceleration (rpm/s)
    decel: float = None    # ESCON speed ramp deceleration (rpm/s)
    motor_type: str = "?"
    op_mode: str = "?"
    max_speed: float = None      # rpm
    cont_current: float = None   # A

    # ---- mapping
    @property
    def is_speed(self):
        return self.quantity == "speed"

    def value_at(self, duty):
        """Set value (rpm) obtained at this duty (%)."""
        return self.offset + self.slope * duty

    def duty_for(self, value):
        """Duty (%) that gives this set value. Raises ValueError if it falls outside 10..90 %."""
        d = (value - self.offset) / self.slope
        if not DUTY_MIN - 1e-6 <= d <= DUTY_MAX + 1e-6:
            raise ValueError(f"{value:g} {self.unit} needs duty {d:.2f} %, outside {DUTY_MIN:g}..{DUTY_MAX:g} %")
        return d

    @property
    def zero_duty(self):
        """Exact duty giving a zero set value (50 % bipolar, 10 % unipolar)."""
        return (0.0 - self.offset) / self.slope

    @property
    def rest_duty(self):
        """Duty to send to stop the motor, quantised like every other command."""
        return q(self.zero_duty)

    @property
    def kind(self):
        """'unipolar' (0 rpm at 10 %), 'bipolar' (0 rpm in the middle) or 'other'."""
        z = self.zero_duty
        if self.slope > 0 and abs(z - DUTY_MIN) <= 0.5:
            return "unipolar"
        if self.slope > 0 and DUTY_MIN + 0.5 < z < DUTY_MAX - 0.5:
            return "bipolar"
        return "other"

    # ---- display
    def summary(self):
        """Multi-line text with the main parameters (printed by get_config)."""
        u = self.unit
        f = lambda v, s="": "n/a" if v is None else f"{v:g}{s}"
        L = ["-" * 66,
             f"MOTOR / ESCON CONFIGURATION   (source: {self.source}, file: {self.edc_name})",
             "-" * 66,
             f"Motor type            : {self.motor_type}",
             f"Operation mode        : {self.op_mode}",
             f"PWM set value         : {self.quantity} [{u}]  -> {self.kind.upper()} mapping",
             f"  mapping             : {u} = {self.offset:g} + {self.slope:g} * duty[%]",
             f"  offset (at 0 % duty): {self.offset:g} {u}",
             f"  value at 10 % duty  : {self.value_at(10):g} {u}",
             f"  value at 90 % duty  : {self.value_at(90):g} {u}",
             f"  zero ({u}) at duty    : {self.zero_duty:.2f} %  (rest command c{fmt_duty(self.rest_duty)})",
             f"  1 % duty step       : {self.slope:g} {u}  (duty decimals sent: {DUTY_DECIMALS}, "
             f"hardware step 0.39 % = {abs(self.slope) * 100 / PWM_COUNTS:.2f} {u})",
             f"ESCON speed ramp      : acceleration {f(self.accel, ' rpm/s')}, deceleration {f(self.decel, ' rpm/s')}",
             f"Max motor speed       : {f(self.max_speed, ' rpm')}",
             f"Continuous current    : {f(self.cont_current, ' A')}",
             "-" * 66]
        if self.kind == "other":
            L.append("WARNING: neither the standard unipolar (0 at 10 %) nor bipolar mapping.")
        return "\n".join(L)

    def meta(self):
        """Key/value pairs written in the header of the CSV files (None values are skipped by the writer)."""
        return {"edc": self.edc_name, "mapping": self.kind, "quantity": self.quantity, "unit": self.unit,
                "offset": f"{self.offset:g}", "slope_per_pct": f"{self.slope:g}",
                "zero_duty_pct": f"{self.zero_duty:.3f}", "ramp_accel_rpm_s": self.accel,
                "ramp_decel_rpm_s": self.decel}

    # ---- (de)serialisation, used to hand the config from the session to its sub-processes
    def dump_temp(self):
        """Write the config as JSON in a temp file; returns its path (to put in ENV_CFG)."""
        fd, path = tempfile.mkstemp(prefix="chipiron_motor_", suffix=".json")
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump(asdict(self), fh)
        return path

    @classmethod
    def from_json(cls, path):
        """Rebuild a config written by dump_temp()."""
        with open(path, encoding="utf-8") as fh:
            d = json.load(fh)
        d["source"] = "session"
        return cls(**d)

    @classmethod
    def from_edc_object(cls, edc, path, source="edc"):
        """Build from an escon_edc_v3.EsconConfig (already loaded: no second read of the file)."""
        sp = edc.setpoint
        if sp is None or not sp.slope_per_pct:
            raise ValueError("the .edc defines no usable PWM set value (non-zero scaling needed).")

        def num(uid):
            v = edc.try_value(uid)
            return None if v is None else float(v)

        qty = sp.quantity
        return cls(source=source, edc_name=Path(str(path)).name, quantity=qty,
                   unit={"speed": "rpm", "current": "A"}.get(qty, "?"),
                   offset=float(sp.offset), slope=float(sp.slope_per_pct),
                   accel=num("PARAM_SpeedRampAcceleration"), decel=num("PARAM_SpeedRampDeceleration"),
                   motor_type=edc.enum_label("PARAM_MotorType") or "?",
                   op_mode=edc.enum_label("PARAM_ModeOfOperation") or "?",
                   max_speed=num("PARAM_MotorMaximumSpeed"), cont_current=num("PARAM_MotorContinousCurrentLimit"))


def _default_edc():
    """Last used .edc, else the most recent .edc next to this file / in the current folder."""
    try:
        p = LAST_EDC_FILE.read_text(encoding="utf-8").strip()
        if p and Path(p).is_file():
            return p
    except OSError:
        pass
    cands = list(SCRIPT_DIR.glob("*.edc")) + list(Path.cwd().glob("*.edc"))
    return str(max(cands, key=lambda p: p.stat().st_mtime)) if cands else ""


def manual_config():
    """Fallback when no .edc can be read: type the mapping by hand."""
    print("Manual mapping (no .edc available).")
    r10 = ask_float("Speed at 10 % duty [rpm]", 0)
    r90 = ask_float("Speed at 90 % duty [rpm]", 1000)
    if r90 == r10:
        raise SystemExit("ERROR: speeds at 10 % and 90 % must differ.")
    slope = (r90 - r10) / 80.0
    acc = ask_float("ESCON speed ramp acceleration [rpm/s]", 10)
    dec = ask_float("ESCON speed ramp deceleration [rpm/s]", acc)
    return MotorConfig("manual", "(manual entry)", "speed", "rpm", r10 - slope * 10, slope, acc, dec)


def get_config(edc_arg=None):
    """Session -> config received through the environment (NO .edc reload). Standalone -> load the .edc.
    Always prints the main parameters. Falls back to a manual mapping if the .edc cannot be used."""
    env = os.environ.get(ENV_CFG)
    if env and Path(env).is_file():
        cfg = MotorConfig.from_json(env)
        print(f"Motor configuration received from chipiron_session_v3 ({cfg.edc_name}): "
              "the .edc is NOT loaded again.")
    else:
        path = (edc_arg or ask("ESCON .edc parameter file (the one loaded in ESCON Studio)", _default_edc()))
        path = Path(path.strip('"').strip("'")).expanduser()
        try:
            from escon_edc_v3 import load_edc
            cfg = MotorConfig.from_edc_object(load_edc(path), path)
            try:
                LAST_EDC_FILE.write_text(str(path.resolve()), encoding="utf-8")
            except OSError:
                pass
        except Exception as exc:
            print(f"Cannot use {path}: {exc}")
            cfg = manual_config()
    print(cfg.summary())
    return cfg


# ------------------------------------------------------------------ links: how commands reach the Arduino
def resolve_phone_ip(cli_ip=None):
    """Priority: --phone-ip > env CHIPIRON_PHONE_IP > DEFAULT_PHONE_IP. 'http://1.2.3.4:8080' is accepted."""
    ip = cli_ip or os.environ.get(ENV_IP) or DEFAULT_PHONE_IP
    return re.sub(r"^https?://", "", str(ip)).split(":")[0]


class HttpLink:
    """Commands through the phone's HTTP server: GET http://<ip>:8080/usb-cmd/<command>.

    record=True  -> start() calls /start and stop() calls /stop (+ /status): the script drives the acquisition.
    record=False -> start()/stop() do nothing (the acquisition is handled elsewhere, e.g. by the session).
    verbose=True -> print every HTTP request and the phone's answer (useful for diagnostics)."""
    name = "http"

    def __init__(self, ip=None, port=PHONE_PORT, record=True, verbose=False, timeout=5):
        self.ip = resolve_phone_ip(ip)
        self.base = f"http://{self.ip}:{port}"
        self.record = record
        self.verbose = verbose
        self.timeout = timeout

    def _get(self, path):
        url = f"{self.base}{path}"
        try:
            with urllib.request.urlopen(url, timeout=self.timeout) as resp:
                body = resp.read().decode("utf-8", errors="replace")
                if self.verbose:
                    print(f"HTTP {resp.status}: GET {path}" + (f"\n  response: {body[:500]}" if body.strip() else ""))
                return body
        except HTTPError as exc:
            body = exc.read().decode("utf-8", errors="replace")
            raise RuntimeError(f"phone server rejected GET {path}: HTTP {exc.code}; response={body[:500]!r}") from exc
        except OSError as exc:          # URLError, timeout, connection refused...
            raise RuntimeError(f"cannot reach the phone server at {url}: {exc}") from exc

    def start(self):
        if not self.record:
            print("--no-record: not calling /start (acquisition handled elsewhere).")
            return
        self._get("/start")

    def send(self, raw):
        """Send one Arduino command. ',' and '.' are kept literal: the phone server does not decode '%2C',
        so the Arduino would never see the commas of 'o<mean>,<amp>,<freq>'."""
        encoded = quote(raw, safe=",.")
        if self.verbose:
            print(f"Sending Arduino command: {raw!r}  (URL path: /usb-cmd/{encoded})")
        return self._get(f"/usb-cmd/{encoded}")

    def stop(self):
        if not self.record:
            return
        try:
            self._get("/stop")
            time.sleep(0.1)
            print("phone status:", self._get("/status"))
        except Exception as exc:
            print(f"Phone stop/status request failed: {exc}")


class LoraLink:
    """Commands through the LoRa bridge Arduino (lora_tx_bridge) plugged into the laptop.
    The phone acquisition is NOT started by this link (the phone may be out of wifi range): start it yourself."""
    name = "lora"

    def __init__(self, port, baud=115200):
        import serial                       # pip install pyserial
        self.ser = serial.Serial(port, baud, timeout=4)
        time.sleep(2)                       # the Arduino resets when the serial port is opened

    def start(self):
        print("LoRa link: start the acquisition on the phone yourself.")

    def send(self, raw):
        self.ser.write((raw + "\n").encode())
        reply = self.ser.readline().decode(errors="replace").strip()
        print("bridge:", reply)
        return reply

    def stop(self):
        self.ser.close()


def make_link(kind="http", ip=None, serial_port=None, record=True, verbose=False):
    """Factory used by run / oscillate scripts (--link http|lora, --phone-ip, --serial-port, --no-record)."""
    if kind == "lora":
        if not serial_port:
            raise SystemExit("--serial-port is required with --link lora (e.g. COM5 or /dev/ttyUSB0)")
        return LoraLink(serial_port)
    return HttpLink(ip=ip, record=record, verbose=verbose)


# ------------------------------------------------------------------ ESCON speed-ramp limiter model
def _legs(s, r, acc, dec):
    """Split a speed change s -> r into legs (start, end, rate): braking to 0 then accelerating if it changes sign."""
    if s * r < 0:
        return [(s, 0.0, dec), (0.0, r, acc)]
    return [(s, r, acc)] if abs(r) >= abs(s) else [(s, r, dec)]


def min_time(cfg, s, r):
    """Shortest time the ESCON limiter needs to go from speed s to r (rpm). 0 if no limit known."""
    if not (cfg.is_speed and cfg.accel and cfg.decel):
        return 0.0
    return sum(abs(b - a) / rate for a, b, rate in _legs(s, r, cfg.accel, cfg.decel))


def _slew(cfg, s, r, T):
    """Speed reached after T s when the target jumps from s to r, and the time that would be needed."""
    need = min_time(cfg, s, r)
    if need <= T or need == 0:
        return r, need
    left, cur = T, s
    for a, b, rate in _legs(s, r, cfg.accel, cfg.decel):
        t_leg = abs(b - a) / rate
        if t_leg <= left:
            left -= t_leg
            cur = b
        else:
            cur = a + (b - a) * left / t_leg
            break
    return cur, need


def check_timeline(cfg, rows, max_warnings=8):
    """Simulate the ESCON limiter on the rows (start at rest). Returns a list of warning strings."""
    if not (cfg.is_speed and cfg.accel and cfg.decel):
        return ["(no speed-ramp limit available: ramp / plateau check skipped)"]
    warns, s = [], 0.0
    lag_by_seg = {}
    for r in rows:
        tgt = cfg.value_at(r["duty"])
        s_end, need = _slew(cfg, s, tgt, r["dur"])
        if r["kind"] == "hold" and need > r["dur"] + 0.05:
            warns.append(f"plateau '{r['seg']}': {tgt:g} rpm for {r['dur']:g} s, but the ESCON ramp needs "
                         f"{need:.1f} s from {s:g} rpm -> plateau not reached (stops at {s_end:.0f} rpm).")
        elif r["kind"] == "ramp":
            lag = abs(s_end - tgt)
            lag_by_seg[r["seg"]] = max(lag_by_seg.get(r["seg"], 0.0), lag)
        s = s_end
    tol = 3 * abs(cfg.slope)
    for seg, lag in lag_by_seg.items():
        if lag > tol:
            warns.append(f"ramp '{seg}': commanded slope exceeds the ESCON limit "
                         f"(accel {cfg.accel:g} / decel {cfg.decel:g} rpm/s); real speed lags by up to "
                         f"{lag:.0f} rpm. Lengthen the ramp or raise the ramp in ESCON Studio.")
    if len(warns) > max_warnings:
        warns = warns[:max_warnings] + [f"... and {len(warns) - max_warnings} more warning(s)."]
    return warns


# ------------------------------------------------------------------ rows helpers
def merge_rows(rows):
    """Merge consecutive rows with the same duty / kind / segment (no need to resend the same command)."""
    out = []
    for r in rows:
        if out and out[-1]["duty"] == r["duty"] and out[-1]["kind"] == r["kind"] and out[-1]["seg"] == r["seg"]:
            out[-1] = dict(out[-1], dur=out[-1]["dur"] + r["dur"])
        else:
            out.append(dict(r))
    return out


def total_time(rows):
    """Total duration (s) of a list of rows."""
    return sum(r["dur"] for r in rows)


# ------------------------------------------------------------------ CSV plan I/O
COLUMNS = ["t_start_s", "duration_s", "duty_pct", "rpm", "omega_rad_s", "kind", "segment"]


def write_plan_csv(path, rows, cfg, extra_meta=None):
    """Write rows as CSV: '# key: value' header lines (mapping of the .edc...) then one line per command."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    meta = {"created": datetime.now().isoformat(timespec="seconds"), **cfg.meta(), **(extra_meta or {})}
    with open(path, "w", newline="", encoding="utf-8") as fh:
        for k, v in meta.items():
            if v is not None:
                fh.write(f"# {k}: {v}\n")
        w = csv.writer(fh)
        w.writerow(COLUMNS)
        t = 0.0
        for r in rows:
            rpm = cfg.value_at(r["duty"]) if cfg.is_speed else ""
            om = rpm * 2 * math.pi / 60 if cfg.is_speed else ""
            w.writerow([f"{t:.3f}", f"{r['dur']:.3f}", f"{r['duty']:g}",
                        f"{rpm:.4f}" if rpm != "" else "", f"{om:.5f}" if om != "" else "", r["kind"], r["seg"]])
            t += r["dur"]
    return path


def read_plan_csv(path, cfg=None):
    """Return (rows, meta). duty_pct is what is sent; if absent it is rebuilt from the rpm column (needs cfg)."""
    meta, lines = {}, []
    with open(path, encoding="utf-8-sig") as fh:
        for line in fh:
            if line.startswith("#"):
                if ":" in line:
                    k, v = line[1:].split(":", 1)
                    meta[k.strip()] = v.strip()
            elif line.strip():
                lines.append(line)
    rows = []
    for i, d in enumerate(csv.DictReader(lines), start=1):
        dur = float(d.get("duration_s") or d.get("time_per_point_s") or 0)
        if d.get("duty_pct"):
            duty = float(d["duty_pct"])
        elif d.get("rpm") and cfg is not None:
            duty = cfg.duty_for(float(d["rpm"]))
        else:
            raise ValueError(f"row {i}: no duty_pct (and no rpm + mapping to convert it).")
        if dur <= 0:
            raise ValueError(f"row {i}: duration must be > 0.")
        if not DUTY_MIN - 1e-6 <= duty <= DUTY_MAX + 1e-6:
            raise ValueError(f"row {i}: duty {duty:g} % outside {DUTY_MIN:g}..{DUTY_MAX:g} %.")
        rows.append({"dur": dur, "duty": q(duty), "kind": d.get("kind") or "hold", "seg": d.get("segment") or f"row{i}"})
    if not rows:
        raise ValueError("no data rows in the CSV.")
    return rows, meta


def mapping_mismatch(cfg, meta):
    """Text if the CSV was planned with another mapping than the current one, else None."""
    try:
        off, slope = float(meta["offset"]), float(meta["slope_per_pct"])
    except (KeyError, ValueError):
        return "CSV has no mapping header: cannot check it was planned with the current .edc."
    d10 = abs((off + slope * 10) - cfg.value_at(10))
    d90 = abs((off + slope * 90) - cfg.value_at(90))
    if max(d10, d90) > 1:
        return (f"CSV was planned with another mapping ({meta.get('edc', '?')}: 10 % -> {off + slope * 10:g}, "
                f"90 % -> {off + slope * 90:g} {cfg.unit}) than the current one "
                f"(10 % -> {cfg.value_at(10):g}, 90 % -> {cfg.value_at(90):g}). The duty values are sent as they are, "
                "so the speeds will NOT be the planned ones.")
    return None


def safe_name(s):
    """File-name-safe version of a string."""
    return re.sub(r"[^A-Za-z0-9._+\-=@]", "_", s)