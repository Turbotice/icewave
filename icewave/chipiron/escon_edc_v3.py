"""
escon_edc_v3.py -- read an ESCON Studio parameter file (.edc), print a readable summary,
optionally save it as .txt, and give the conversions
"Arduino ADC counts -> volts -> real units" and "PWM duty % -> set value".

DEPENDENCIES   numpy, stdlib (xml).  No other Chipiron script.
USED BY        chipiron_helpers_v3.py (get_config), chipiron_session_v3.py, chipiron_plotting_v3.py

COMMAND LINE
    py escon_edc_v3.py                          # asks for the file, then asks if/where to save the summary
    py escon_edc_v3.py parameters_Sep26.edc     # asks only if/where to save
    py escon_edc_v3.py parameters_Sep26.edc --folder "D:/ESCON_summaries"   # saves the .txt, no question
    py escon_edc_v3.py parameters_Sep26.edc --no-save                       # console only

AS A MODULE
    from escon_edc_v3 import load_edc
    edc = load_edc("parameters_Sep26.edc")
    print(edc.report())                         # full key-parameter summary (incl. speed controller P gains)
    edc.save_summary("D:/ESCON_summaries")      # -> .txt, returns the path
    edc.setpoint_from_duty(50)                  # -> 0 rpm (bipolar mapping)   [duty in %]
    edc.setpoint_from_duty([10, 50, 90])        # -> array([-800., 0., 800.])
    rpm = edc.ao["AO2"].volts_to_phys(volts)    # AO2 = actual speed in the saved file
    edc.apply_to_frame(df)                      # adds speed / current columns to a decoded Usb DataFrame
    edc.enum_label("PARAM_MotorType")           # 'maxon DC motor'
    edc.try_value("PARAM_SpeedRampAcceleration")  # 2000.0 (rpm/s) or None if missing

Assumptions (the .edc does not contain the <conditionalSupport> definitions, so the unit that applies to a
parameter is chosen from the analog-output *function*):
  - function name contains "Current" -> A, "Speed" -> rpm, "Temperature" -> degC
  - AO "Scaling" has unit V/A, V/rpm or V/degC accordingly (multiplier read from the file)
  - the two-point definition (FirstSource/Target, SecondSource/Target) is used as a cross-check of
    Scaling/Offset; a warning is printed if they disagree.

OFFSET ADJUST 
  The file may hold "Set Value Offset Adjust Fixed Value" (PARAM_SetValueOffsetAdjustFixedValue, written by
  chipiron_patch_edc_v3.py --offset). It only acts if the functionality is enabled in the analog-input
  execution mask (bit 0x0010). If enabled it is ADDED to the PWM mapping offset (assumption: additive, verify
  in ESCON Studio); if not enabled it is ignored and a warning says so.
"""
from __future__ import annotations

import argparse
import datetime
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# Default folder proposed when the user is asked where to save the .txt summary
DEFAULT_SUMMARY_FOLDER = Path("./edc_summaries")
OFFSET_ADJUST_MASK_BIT = 0x0010     # bit of PARAM_AnalogInputsFunctionalitiesExecutionMask ("Offset Adjust Set Value")

_QUANTITIES = [  # (text in function name, quantity, unit, condition suffix used in the .edc)
    ("Current", "current", "A", "Current"),
    ("Speed", "speed", "rpm", "Speed"),
    ("Temperature", "temperature", "degC", "Temp"),
]


@dataclass
class AOChannel:
    """One ESCON analog output: volts = offset_V + slope_V_per_unit * quantity."""
    name: str                  # "AO1" / "AO2"
    function_code: int
    function: str              # e.g. "CurrentActualValue"
    quantity: str              # "current" / "speed" / "temperature" / "none" / "fixed"
    unit: str
    slope_V_per_unit: float = np.nan
    offset_V: float = np.nan

    @property
    def usable(self):
        return (self.quantity in ("current", "speed", "temperature")
                and bool(np.isfinite(self.slope_V_per_unit)) and self.slope_V_per_unit != 0)

    def volts_to_phys(self, volts):
        """Measured volts -> physical value (rpm, A, degC)."""
        return (np.asarray(volts, dtype=float) - self.offset_V) / self.slope_V_per_unit

    def phys_to_volts(self, value):
        """Physical value -> expected volts at the analog output."""
        return self.offset_V + self.slope_V_per_unit * np.asarray(value, dtype=float)


@dataclass
class Setpoint:
    """PWM set value: setpoint = offset + slope * duty(%)."""
    quantity: str
    unit: str
    slope_per_pct: float
    offset: float

    def __call__(self, duty_pct):
        return self.offset + self.slope_per_pct * np.asarray(duty_pct, dtype=float)


def _g(v, nd=6):
    """Compact number formatting (None -> 'n/a')."""
    return "n/a" if v is None else f"{v:.{nd}g}"


def _tok(v, nd=0):
    """Number -> filename-safe token (-800 -> m800)."""
    if v is None:
        return "NA"
    return f"{v:.{nd}f}".replace("-", "m").replace(".", "p")


class EsconConfig:
    """Parsed .edc: raw parameters, enum labels, analog outputs (self.ao) and PWM set value (self.setpoint)."""

    def __init__(self, path):
        self.path = str(path)
        root = ET.parse(path).getroot()
        self.params = {p.get("uniqueID"): p for p in root.iter("parameter")}
        if not self.params:
            raise ValueError(f"no <parameter> found in {path}: is it an ESCON .edc file?")
        self.ao_functions = {}
        self.enums = {}        # enum uniqueID -> {value: english label}
        for e in root.iter("enum"):
            labels = {}
            for ee in e.findall("enumEntry"):
                val = int(ee.findtext("value"))
                en = [l.text.strip() for l in ee.findall("label") if l.get("lang") == "en" and l.text]
                labels[val] = en[0] if en else ee.get("uniqueID")
            self.enums[e.get("uniqueID")] = labels
            if e.get("uniqueID") == "EAnalogOutputFunctionality":
                for ee in e.findall("enumEntry"):
                    self.ao_functions[int(ee.findtext("value"))] = ee.get("uniqueID").split("_", 1)[1]
        hdr = root.find(".//ProfileBody")
        self.file_name = hdr.get("fileName") if hdr is not None else None
        self.file_date = hdr.get("fileCreationDate") if hdr is not None else None
        self.product = (root.findtext(".//productName") or "").strip() or None
        self.warnings = []
        self.offset_adjust = None            # (value rpm, enabled?) when an offset adjust is stored in the file
        self.ao = {"AO1": self._build_ao(1), "AO2": self._build_ao(2)}
        self.setpoint = self._build_setpoint()

    # ---- low-level access -------------------------------------------------
    def raw(self, uid):
        """Integer stored in <actualValue> (hex '0x..' accepted)."""
        return int(self.params[uid].findtext("actualValue").strip(), 0)

    def mult(self, uid, cond=None):
        """Multiplier of the unit of parameter uid (the one tied to a condition ending with `cond`)."""
        units = self.params[uid].findall("unit")
        pick = None
        for u in units:
            ref = u.find("conditionalSupportIDRef")
            if cond is None and ref is None:
                pick = u
                break
            if cond is not None and ref is not None and ref.get("uniqueIDRef").endswith(cond):
                pick = u
                break
        if pick is None and len(units) == 1:
            pick = units[0]
        if pick is None:
            raise KeyError(f"{uid}: no unit matching condition '{cond}'")
        m = pick.findtext("multiplier")
        return float(m) if m else 1.0

    def value(self, uid, cond=None):
        """Physical value = raw integer x unit multiplier."""
        return self.raw(uid) * self.mult(uid, cond)

    def try_value(self, uid, cond=None):
        """Like value() but returns None if the parameter / unit is missing from the file."""
        if uid not in self.params:
            return None
        try:
            return self.value(uid, cond)
        except (KeyError, ValueError, AttributeError):
            return None

    def enum_label(self, uid):
        """English label of an enum-typed parameter (e.g. motor type, mode of operation)."""
        p = self.params.get(uid)
        if p is None:
            return None
        ref = p.find("dataTypeIDRef")
        try:
            raw = self.raw(uid)
        except (ValueError, AttributeError):
            return None
        if ref is None:
            return None
        return self.enums.get(ref.get("uniqueIDRef"), {}).get(raw, f"unknown code {raw}")

    def _warn(self, msg):
        self.warnings.append(msg)
        print(f"[edc] warning: {msg}")

    # ---- analog outputs (ESCON -> Arduino A1/A2) ---------------------------
    def _build_ao(self, n):
        code = (self.raw("PARAM_AnalogOutputsConfiguration") >> (4 * (n - 1))) & 0xF   # AO1 = low nibble
        fname = self.ao_functions.get(code, f"code{code}")
        ch = AOChannel(f"AO{n}", code, fname, "none", "")
        if fname == "FixValue":
            ch.quantity = "fixed"
            return ch
        for key, quantity, unit, cond in _QUANTITIES:
            if key in fname:
                break
        else:
            return ch                                    # function "None" or unknown
        ch.quantity, ch.unit = quantity, unit
        p = f"PARAM_AnalogOutput{n}"
        c = f"COND_AnalogOutput0{n}Configuration_{cond}"
        ch.slope_V_per_unit = self.value(f"{p}Scaling", c)
        ch.offset_V = self.value(f"{p}Offset")
        # cross-check with the two-point definition
        t1, t2 = self.value(f"{p}FirstTarget"), self.value(f"{p}SecondTarget")
        s1, s2 = self.value(f"{p}FirstSource", c), self.value(f"{p}SecondSource", c)
        if s2 != s1:
            slope2 = (t2 - t1) / (s2 - s1)
            off2 = t1 - slope2 * s1
            if not (np.isclose(slope2, ch.slope_V_per_unit, rtol=1e-2) and np.isclose(off2, ch.offset_V, atol=0.01)):
                self._warn(f"AO{n}: Scaling/Offset ({ch.slope_V_per_unit:g} V/{unit}, {ch.offset_V:g} V) "
                           f"disagree with the two-point definition ({slope2:g} V/{unit}, {off2:g} V); using Scaling/Offset")
        return ch

    # ---- PWM set value (Arduino -> ESCON) ------------------------------------
    def _build_setpoint(self):
        try:
            mode = self.raw("PARAM_ModeOfOperation")      # 1 current, 2 speed closed loop, 3 speed open loop
            if self.raw("PARAM_AnalogInputsConfiguration") >> 8 & 0xF != 1:
                self._warn("PWM input is not configured as 'Set Value' in the analog-inputs configuration")
            if mode == 1:
                slope = self.value("PARAM_SetValueCurrentScaling", "AnalogInputPwm_CurrentSetvalue")
                return Setpoint("current", "A", slope, self.value("PARAM_SetValueCurrentOffset"))
            slope = self.value("PARAM_SetValueSpeedScaling", "AnalogInputPwm_SpeedSetpoint")
            sp = Setpoint("speed", "rpm", slope, self.value("PARAM_SetValueSpeedOffset"))
            self._apply_offset_adjust(sp)
            return sp
        except KeyError as exc:
            self._warn(f"could not build the PWM set value mapping ({exc})")
            return None

    def _apply_offset_adjust(self, sp):
        """Include 'Set Value Offset Adjust Fixed Value' in the mapping only if it is enabled in the file."""
        adj = self.try_value("PARAM_SetValueOffsetAdjustFixedValue", "SpeedRegulator")
        if not adj:
            return
        try:
            mask = self.raw("PARAM_AnalogInputsFunctionalitiesExecutionMask")    # bit field, no unit multiplier
        except (KeyError, ValueError, AttributeError):
            mask = 0
        enabled = bool(mask & OFFSET_ADJUST_MASK_BIT)
        self.offset_adjust = (adj, enabled)
        if enabled:
            sp.offset += adj
            self._warn(f"offset adjust of {adj:g} rpm is ENABLED: added to the PWM mapping offset "
                       "(assumed additive, check in ESCON Studio)")
        else:
            self._warn(f"offset adjust fixed value {adj:g} rpm is stored in the file but NOT enabled "
                       f"(analog-input execution mask {mask:#06x}, bit {OFFSET_ADJUST_MASK_BIT:#06x} is 0): "
                       "ignored in the PWM mapping. Check in ESCON Studio.")

    def setpoint_from_duty(self, duty_pct):
        """Set value (rpm or A) for a duty cycle in %, scalar or array."""
        if self.setpoint is None:
            raise RuntimeError("no PWM set value mapping available")
        return self.setpoint(duty_pct)

    def pwm_polarity(self):
        """'bipolar', 'unipolar (positive)', 'unipolar (negative)' or None, from the 10 % and 90 % set values."""
        sp = self.setpoint
        if sp is None:
            return None
        lo, hi = float(sp(10)), float(sp(90))
        if min(lo, hi) < 0 < max(lo, hi):
            return "bipolar"
        return "unipolar (positive)" if min(lo, hi) >= 0 else "unipolar (negative)"

    # ---- apply to a decoded Usb DataFrame ----------------------------------
    def apply_to_frame(self, df, wiring=None, vref=5.0):
        """Add physical columns to a decoded Usb DataFrame (needs A1_V / A2_V, see chipiron_decode_v3).
        wiring: which ESCON output feeds which Arduino pin. Default = wiring verified on 2026-10-20
        (ESCON Studio, motor held at 600 rpm, compared with the Arduino readings):
        A1 <- AO1 (actual current), A2 <- AO2 (actual speed). If you rewire, pass wiring explicitly."""
        wiring = wiring or {"A1": "AO1", "A2": "AO2"}
        labels = {}
        for pin, ao in wiring.items():
            c = self.ao[ao]
            if not c.usable:
                self._warn(f"{ao} is '{c.function}' in the .edc: nothing to convert for {pin}")
                continue
            col = f"{pin}_{c.quantity}_{c.unit}"
            df[col] = c.volts_to_phys(df[f"{pin}_V"])
            labels[col] = f"{c.quantity} ({c.unit})"
            print(f"[edc] {pin} <- {ao} ({c.function}): {col} = ({pin}_V - {c.offset_V:g}) / {c.slope_V_per_unit:g}")
        df.attrs["labels"] = labels
        return df

    # ---- summaries ---------------------------------------------------------
    def summary(self, vref=5.0, adc_bits=10):
        """Short block: analog-output and PWM conversions (used inside report())."""
        lsb = vref / (2 ** adc_bits - 1)
        lines = [f"ESCON file {self.file_name or Path(self.path).name} ({self.file_date or 'date n/a'})"]
        for name, c in self.ao.items():
            if c.usable:
                lines.append(f"  {name}: {c.function}: V = {c.offset_V:g} + {c.slope_V_per_unit:g} * {c.quantity}[{c.unit}]"
                             f"   (1 ADC count = {lsb / c.slope_V_per_unit:.3g} {c.unit})")
            else:
                lines.append(f"  {name}: {c.function}")
        sp = self.setpoint
        if sp:
            lines.append(f"  PWM set value ({sp.quantity}): {sp.offset:g} + {sp.slope_per_pct:g} * duty[%] {sp.unit}"
                         f"   (10 % -> {sp(10):g}, 50 % -> {sp(50):g}, 90 % -> {sp(90):g} {sp.unit})")
        return "\n".join(lines)

    def report(self, vref=5.0, adc_bits=10):
        """Full key-parameter summary (string), same text for console and .txt."""
        bar = "=" * 72
        sub = lambda t: f"\n{t}\n" + "-" * 72
        L = [bar, "ESCON PARAMETER SUMMARY", bar,
             f"Source file    : {Path(self.path).resolve()}",
             f"Profile        : {self.file_name or Path(self.path).name} (created {self.file_date or 'n/a'})",
             f"Device         : {self.product or 'n/a'}",
             f"Report written : {datetime.datetime.now():%Y-%m-%d %H:%M:%S}"]

        # motor / mode
        enc = self.try_value("PARAM_EncoderPulseNumber")
        sensor = self.enum_label("PARAM_SpeedFeedbackSource") or "n/a"
        if enc is not None and "Encoder" in sensor:
            sensor += f" ({enc:g} counts/turn)"
        L += [sub("MOTOR / MODE"),
              f"Motor type            : {self.enum_label('PARAM_MotorType') or 'n/a'}",
              f"Mode of operation     : {self.enum_label('PARAM_ModeOfOperation') or 'n/a'}",
              f"Speed sensor          : {sensor}",
              f"Nominal current       : {_g(self.try_value('PARAM_MotorContinousCurrentLimit'))} A",
              f"Max. output current   : {_g(self.try_value('PARAM_CurrentLimit'))} A",
              f"Max. permissible speed: {_g(self.try_value('PARAM_MotorMaximumSpeed'))} rpm"]

        # PWM set value
        L.append(sub("PWM SET VALUE  (Arduino -> ESCON)"))
        sp = self.setpoint
        if sp is None:
            L.append("PWM mapping not available (see warnings)")
        else:
            zero = -sp.offset / sp.slope_per_pct if sp.slope_per_pct else None
            L += [f"Mapping               : set value = {sp.offset:g} + {sp.slope_per_pct:g} * duty[%]  ({sp.unit})",
                  f"Offset                : {sp.offset:g} {sp.unit}",
                  f"Scaling               : {sp.slope_per_pct:g} {sp.unit}/%",
                  f"Set value at 10 % PWM : {float(sp(10)):g} {sp.unit}",
                  f"Set value at 50 % PWM : {float(sp(50)):g} {sp.unit}",
                  f"Set value at 90 % PWM : {float(sp(90)):g} {sp.unit}",
                  f"Duty giving zero      : {_g(zero)} %",
                  f"Polarity              : {self.pwm_polarity()}"]
            if self.offset_adjust:
                adj, on = self.offset_adjust
                L.append(f"Offset adjust (fixed) : {adj:g} {sp.unit}, {'ENABLED (included above)' if on else 'NOT enabled (ignored)'}")

        # ramps
        L += [sub("SPEED RAMP  (speed modes only)"),
              f"Acceleration          : {_g(self.try_value('PARAM_SpeedRampAcceleration'))} rpm/s",
              f"Deceleration          : {_g(self.try_value('PARAM_SpeedRampDeceleration'))} rpm/s",
              f"Deceleration (stop)   : {_g(self.try_value('PARAM_SpeedRampDecelerationStop'))} rpm/s"]

        # speed controller gains (raw integers, no unit; Speed 1 is used with "Current Limiter", Speed 2 with
        # "Inner Current Control Loop"; chipiron_patch_edc_v3.py --speed-p-gain writes the one in use)
        def _raw(uid):
            try:
                return self.raw(uid)
            except (KeyError, ValueError, AttributeError):
                return None
        ctrl = _raw("PARAM_ModeOfRegulator")
        in_use = "2" if ctrl == 2 else "1"
        L += [sub("SPEED CONTROLLER GAINS"),
              f"Controller mode       : {self.enum_label('PARAM_ModeOfRegulator') or 'n/a'}"
              f"  -> speed {in_use} regulator in use",
              f"Speed 1 P gain        : {_g(_raw('PARAM_Speed1RegulatorP_Gain'))}"
              + ("   <- in use" if in_use == "1" else ""),
              f"Speed 2 P gain        : {_g(_raw('PARAM_Speed2RegulatorP_Gain'))}"
              + ("   <- in use" if in_use == "2" else "")]

        # analog outputs
        L += [sub("ANALOG OUTPUTS -> ADC conversion  (Vref %.3g V, %d bit)" % (vref, adc_bits))]
        L += self.summary(vref, adc_bits).splitlines()[1:]

        if self.warnings:
            L += [sub("WARNINGS")] + [f"  - {w}" for w in self.warnings]
        L += ["", "Values are read from the .edc file (not measured on the controller).", bar]
        return "\n".join(L)

    def summary_filename(self):
        """File name that encodes the main values, e.g. parameters_Sep26_DC_bipolar_m800to800rpm_ramp2000.txt"""
        mtype = (self.enum_label("PARAM_MotorType") or "").upper()
        motor = "DC" if "DC" in mtype else "EC" if "EC" in mtype else "MotorNA"
        pol = (self.pwm_polarity() or "polNA").split(" ")[0]
        sp = self.setpoint
        lo, hi, unit = (float(sp(10)), float(sp(90)), sp.unit) if sp else (None, None, "")
        ramp = self.try_value("PARAM_SpeedRampAcceleration")
        return f"{Path(self.path).stem}_{motor}_{pol}_{_tok(lo)}to{_tok(hi)}{unit}_ramp{_tok(ramp)}.txt"

    def save_summary(self, folder, **kw):
        """Write report() as .txt into `folder` (created if needed). Returns the Path."""
        folder = Path(folder).expanduser()
        folder.mkdir(parents=True, exist_ok=True)
        out = folder / self.summary_filename()
        out.write_text(self.report(**kw), encoding="utf-8")
        return out


def load_edc(path):
    """Parse an .edc file and return an EsconConfig."""
    return EsconConfig(path)


def _clean_path(s):
    """Strip quotes added by drag-and-drop / copy-as-path on Windows."""
    return s.strip().strip('"').strip("'")


def main(argv=None):
    ap = argparse.ArgumentParser(description="Read an ESCON .edc file, print a summary, optionally save it as .txt")
    ap.add_argument("edc", nargs="?", help="path to the .edc file (asked if omitted)")
    ap.add_argument("--folder", help="save the .txt summary in this folder without asking")
    ap.add_argument("--no-save", action="store_true", help="console summary only, do not ask about saving")
    args = ap.parse_args(argv)

    path = _clean_path(args.edc or input("Path to the .edc file: "))
    if not Path(path).is_file():
        print(f"[edc] error: file not found: {path}")
        return 1

    print(f"[edc] loading {path}")
    try:
        cfg = load_edc(path)
    except (ET.ParseError, ValueError, KeyError) as exc:
        print(f"[edc] error: could not read the file ({exc})")
        return 1
    print(f"[edc] loaded {len(cfg.params)} parameters, {len(cfg.enums)} enum definitions")
    print(f"[edc] analog outputs: " + ", ".join(f"{k}={v.function}" for k, v in cfg.ao.items()))
    print(f"[edc] PWM mapping: {'ok' if cfg.setpoint else 'not available'}, {len(cfg.warnings)} warning(s)\n")

    print(cfg.report())
    print()

    if args.no_save:
        print("[edc] not saved (--no-save)")
        return 0
    folder = args.folder
    if folder is None:
        if input("Save this summary as a .txt file? [y/N]: ").strip().lower() not in ("y", "yes", "o", "oui"):
            print("[edc] not saved")
            return 0
        folder = _clean_path(input(f"Folder to save in [{DEFAULT_SUMMARY_FOLDER}]: ")) or DEFAULT_SUMMARY_FOLDER
    try:
        out = cfg.save_summary(folder)
    except OSError as exc:
        print(f"[edc] error: could not write to '{folder}' ({exc})")
        return 1
    print(f"[edc] summary saved to: {out.resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())