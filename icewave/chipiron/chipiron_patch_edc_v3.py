"""
chipiron_patch_edc_v3.py -- make a modified COPY of an ESCON .edc parameter file.
Chek in ESCON Studio: open the result, CHECK the values in the GUI, then write them to the controller.

WHY
  The Arduino ADC only reads 0..5 V but a default ESCON config sets the analog outputs to -4..+4 V, so negative
  speed / current read as 0, and the speed output only uses a small part of the range. The "AO shift" re-maps
  both outputs to 0..4 V (0 rpm = 2 V, 0 A = 2 V), readable by the Arduino. The other options set the PWM
  set-value mapping (bipolar), the ESCON speed ramp, the motor type, a set-value offset and the speed P gain.

DEPENDENCIES   stdlib only.  Works on the text of the file (keeps everything else untouched).
USED WITH      escon_edc_v3.py (re-read the patched file to check it), chipiron_helpers_v3.py (consumes the mapping).

OPTIONS  (--in and --out required; the output file must not exist unless --force)
  --ao-shift / --no-ao-shift   remap analog outputs (default ON)
        speed  : -speed_range rpm -> 0 V, +speed_range -> 4 V       (--speed-range, default 1000 rpm)
        current: -current_range A -> 0 V, +current_range -> 4 V     (--current-range, default 5 A)
  --bipolar-setpoint R   PWM set value becomes bipolar: 10 % duty = -R rpm, 50 % = 0, 90 % = +R rpm
  --ramp RPM_PER_S       speed ramp acceleration AND deceleration
  --motor-type EC|DC     write the motor type (only when given!). Enum of the .edc: 0 = EC (brushless),
                         1 = DC (brushed). A wrong value can make the controller drive the motor wrongly.
  --offset RPM           "Set value offset adjust" fixed value (can be negative). It only acts if that
                         functionality is enabled in ESCON Studio; the script prints the current mask state.
  --speed-p-gain N       speed controller P gain (integer, no unit, as shown in ESCON Studio). Written to the
                         regulator actually in use: Speed 1 (mode of controller "Current Limiter") or Speed 2
                         ("Inner Current Control Loop"); the script says which one. Too low a gain flattens /
                         skews fast sines; too high gives noise or oscillation: raise it step by step.
  --force                overwrite --out if it exists

Each parameter is stored in the file as an integer = value / unit multiplier (e.g. rpm x 0.1). The .edc holds BOTH a
"two-point" and a "scaling + offset" description of each mapping; this script changes both.

EXAMPLES  
  # only the analog-output shift (default), brushless motor left as in the file
  py chipiron_patch_edc_v3.py --in parameters_Sep26.edc --out parameters_Sep26_AOshift.edc
  # bipolar +/-800 rpm, ramp 2000 rpm/s, DC motor, no AO change
  py chipiron_patch_edc_v3.py --in parameters_Sep26.edc --out parameters_bip.edc --no-ao-shift \
        --bipolar-setpoint 800 --ramp 2000 --motor-type DC
  # everything, plus a -11 rpm offset
  py chipiron_patch_edc_v3.py --in parameters_Sep26.edc --out parameters_off.edc --bipolar-setpoint 800 \
        --ramp 2000 --offset -11
  # wider analog ranges (speed +/-1500 rpm, current +/-10 A)
  py chipiron_patch_edc_v3.py --in a.edc --out b.edc --speed-range 1500 --current-range 10
  # raise the speed P gain only (file already has the AO shift, so --no-ao-shift)
  py chipiron_patch_edc_v3.py --in a.edc --out a_gain20.edc --no-ao-shift --speed-p-gain 20
"""
import argparse
import os
import re
import sys


def set_actual(text, uid, raw):
    """Replace <actualValue> of the parameter with this uniqueID. Returns (new_text, old_value)."""
    pat = re.compile(r'(<parameter uniqueID="%s"(?:(?!</parameter>).)*?<actualValue>)([^<]*)(</actualValue>)'
                     % re.escape(uid), re.S)
    m = pat.search(text)
    if not m:
        raise SystemExit(f"parameter {uid} not found in the file")
    old = m.group(2)
    return text[:m.start(2)] + str(int(round(raw))) + text[m.end(2):], old


def get_actual(text, uid):
    """Current <actualValue> of a parameter as text (None if absent)."""
    m = re.search(r'<parameter uniqueID="%s"(?:(?!</parameter>).)*?<actualValue>([^<]*)</actualValue>'
                  % re.escape(uid), text, re.S)
    return m.group(1).strip() if m else None


def apply(text, changes):
    """Apply [(uniqueID, raw integer)] and print a before -> after line for each."""
    for uid, raw in changes:
        text, old = set_actual(text, uid, raw)
        print(f"  {uid:48s} {old:>10s} -> {int(round(raw))}")
    return text


def ao_changes(speed_range, current_range):
    # unit multipliers: scaling V/rpm x1e-6, V/A x1e-3, offsets and targets V x1e-3,
    # sources rpm x0.1 and A x1e-4
    v0, v4 = 0.0, 4.0                                    # output voltage at the two points (V)
    sp_scale = (v4 - v0) / (2 * speed_range)             # V per rpm
    cu_scale = (v4 - v0) / (2 * current_range)           # V per A
    mid = (v0 + v4) / 2                                  # voltage at zero (2 V)
    return [
        ("PARAM_AnalogOutput1Scaling", cu_scale / 1e-3),
        ("PARAM_AnalogOutput1Offset", mid / 1e-3),
        ("PARAM_AnalogOutput1FirstTarget", v0 / 1e-3),
        ("PARAM_AnalogOutput1FirstSource", -current_range / 1e-4),
        ("PARAM_AnalogOutput1SecondTarget", v4 / 1e-3),
        ("PARAM_AnalogOutput1SecondSource", current_range / 1e-4),
        ("PARAM_AnalogOutput2Scaling", sp_scale / 1e-6),
        ("PARAM_AnalogOutput2Offset", mid / 1e-3),
        ("PARAM_AnalogOutput2FirstTarget", v0 / 1e-3),
        ("PARAM_AnalogOutput2FirstSource", -speed_range / 0.1),
        ("PARAM_AnalogOutput2SecondTarget", v4 / 1e-3),
        ("PARAM_AnalogOutput2SecondSource", speed_range / 0.1),
    ]


def bipolar_changes(R):
    # duty 10 % -> -R rpm, 90 % -> +R rpm : slope R/40 rpm per %, offset -1.25 R rpm
    return [
        ("PARAM_SetValueSpeedScaling", (R / 40.0) / 1e-3),        # rpm/% x0.001
        ("PARAM_SetValueSpeedOffset", (-1.25 * R) / 0.1),         # rpm x0.1
        ("PARAM_SetValueSpeedFirstTarget", -R / 0.1),
        ("PARAM_SetValueSpeedSecondTarget", R / 0.1),
    ]


def ramp_changes(rpm_s):
    return [("PARAM_SpeedRampAcceleration", rpm_s / 0.1),
            ("PARAM_SpeedRampDeceleration", rpm_s / 0.1)]


def motor_type_changes(mtype):
    # Enum EMotorType of the ESCON .edc (checked in the file): 0 = maxon EC motor (brushless), 1 = maxon DC motor.
    return [("PARAM_MotorType", 1 if mtype.upper() == "DC" else 0)]


def offset_changes(rpm):
    # Set value offset adjust, fixed value (0x4144:03), unit 0.1 rpm in speed mode
    return [("PARAM_SetValueOffsetAdjustFixedValue", rpm / 0.1)]


def speed_gain_changes(text, gain):
    # Speed controller P gain: 0x4201:01 = Speed 1 regulator, 0x4202:01 = Speed 2 regulator. Raw integer, no unit.
    # Speed 1 is used with "Current Limiter" (mode of controller = 3), Speed 2 with "Inner Current Control Loop" (= 2).
    if not 1 <= gain <= 65535:
        sys.exit("--speed-p-gain must be an integer between 1 and 65535.")
    op = get_actual(text, "PARAM_ModeOfOperation")
    if op is not None and int(op, 0) == 1:
        print("WARNING: the file is in current-controller mode: the speed P gain has no effect there.")
    ctrl = get_actual(text, "PARAM_ModeOfRegulator")
    use2 = ctrl is not None and int(ctrl, 0) == 2
    uid = "PARAM_Speed2RegulatorP_Gain" if use2 else "PARAM_Speed1RegulatorP_Gain"
    print(f"NOTE: mode of controller = {ctrl} -> speed {'2' if use2 else '1'} regulator is the one in use: writing {uid}.")
    return [(uid, gain)]


def main():
    p = argparse.ArgumentParser(description="Patch a copy of an ESCON .edc (see the module docstring)")
    p.add_argument("--in", dest="inp", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--ao-shift", dest="ao", action="store_true", default=True)
    p.add_argument("--no-ao-shift", dest="ao", action="store_false")
    p.add_argument("--speed-range", type=float, default=1000.0, help="rpm at +/- full scale of the speed output")
    p.add_argument("--current-range", type=float, default=5.0, help="A at +/- full scale of the current output")
    p.add_argument("--bipolar-setpoint", type=float, metavar="R_RPM")
    p.add_argument("--ramp", type=float, metavar="RPM_PER_S")
    p.add_argument("--motor-type", choices=["EC", "DC"], default=None,
                   help="written ONLY when given (EC = brushless, DC = brushed)")
    p.add_argument("--offset", type=float, metavar="RPM", help="fixed set-value offset in rpm (can be negative)")
    p.add_argument("--speed-p-gain", type=int, metavar="N",
                   help="speed controller P gain (integer, as shown in ESCON Studio)")
    p.add_argument("--force", action="store_true", help="overwrite --out if it exists")
    a = p.parse_args()

    if os.path.abspath(a.inp) == os.path.abspath(a.out):
        sys.exit("--out must differ from --in (the original is never modified).")
    if os.path.exists(a.out) and not a.force:
        sys.exit(f"{a.out} already exists (use --force to overwrite).")

    with open(a.inp, encoding="utf-8") as f:
        text = f.read()
    changes = []
    if a.ao:
        changes += ao_changes(a.speed_range, a.current_range)
    if a.bipolar_setpoint:
        changes += bipolar_changes(a.bipolar_setpoint)
    if a.ramp:
        changes += ramp_changes(a.ramp)
    if a.motor_type:
        print(f"WARNING: writing motor type {a.motor_type}: make sure it is the motor really connected.")
        changes += motor_type_changes(a.motor_type)
    if a.offset is not None:
        changes += offset_changes(a.offset)
    if a.speed_p_gain is not None:
        changes += speed_gain_changes(text, a.speed_p_gain)
    if not changes:
        sys.exit("nothing to change (all options off).")

    print("changes:")
    text = apply(text, changes)
    if a.offset is not None:
        mask = get_actual(text, "PARAM_AnalogInputsFunctionalitiesExecutionMask")
        on = mask is not None and int(mask, 0) & 0x0010
        print(f"NOTE: offset adjust is {'ENABLED' if on else 'NOT enabled'} in the file (analog-input execution mask "
              f"{mask}, bit 0x0010). {'' if on else 'The value is written but has no effect until you enable it in ESCON Studio.'}")
    with open(a.out, "w", encoding="utf-8") as f:
        f.write(text)
    print("written", a.out, "-> open it in ESCON Studio, check the values, then write to the controller.")


if __name__ == "__main__":
    main()