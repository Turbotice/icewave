# Run & save data  -- Chipiron (23/10/2026 -- Matilde)

Set of python/arduino scripts to run an Escon motor (on Chipiron or else) form laptop/fairphone, and collect data. /!\ everything  was tested on an EC Escon motor at PMMH lab, not on Chipiron yet. A1/A2 - speed/current, ranges and wiring must be checked on real motor. 

Python packages: `numpy`, `pandas`, `matplotlib` (plotting only), `scipy` (plotting fits only), `pyserial` (LoRa only).

## Files and who depends on whom
| File | Role | Imports |
|---|---|---|
| `chipiron_helpers_v3.py` | motor config (.edc → mapping, rest duty, ESCON ramp), links http/LoRa, plan CSV, ramp-limiter check | `escon_edc_v3` (lazy) |
| `escon_edc_v3.py` | reads the ESCON `.edc`: mapping, ramp, **speed P gains**, analog-output conversion (summary + .txt) | – |
| `chipiron_run_v3.py` | holds / ramps / staircases / cycles / CSV, runs the motor | helpers |
| `chipiron_oscillate_v3.py` | sine command (firmware or stream mode) | helpers |
| `chipiron_ramp_planner_v3.py` | plans a staircase ramp, saves CSV for run | helpers |
| `chipiron_session_v3.py` | full measurement: Arduino mode → 0 rpm → record → experiment → save | connect, run_gobannos, escon_edc, helpers; launches run / oscillate |
| `chipiron_patch_edc_v3.py` | makes a patched copy of an `.edc` (AO shift, bipolar mapping, ramp, motor type, offset, **speed P gain**) | – |
| `chipiron_decode_v3.py` | decodes the Usb csv (ADC counts / volts; A1 = current, A2 = speed, Vref 4.88 V) | – |
| `chipiron_plotting_v3.py` | plots one acquisition (cells `# %%`), + step / sine response fits (tau vs speed P gain) | decode (+ escon_edc optional) |
| `connect_v3.py` | phone ID → IP | – |
| `run_gobannos_v3.py` | HTTP calls to the phone | – |
| `mobile_v3.py` | manual menu (status, start/stop, ls, pull, usb, time) | connect, run_gobannos, tsync, rw_data |
| `tsync_v3.py` (+ `rw_data_v3.py`) | clock sync  | connect, rw_data |
| `PWM_full_v3/PWM_full_v3.ino` | Arduino firmware  | – |


## Typical workflow
1. (once) ESCON Studio: load your `.edc` configuration for motor (example: bipolar: 10 % = −800, 50 % = 0, 90 % = +800 rpm). Optionally `chipiron_patch_edc_v3.py` to create a new set of parameters (including the speed P gain, see below).
2. Upload `PWM_full_v3.ino` to the Arduino Mega (motor stays **disabled** until a mode command arrives).
3. Plan a ramp (optional): ex. `py chipiron_ramp_planner_v3.py --edc parameters_Sep26.edc`
4. Measure:ex. `py chipiron_session_v3.py full_session mo parameters_Sep26.edc`
5. Plot: `py chipiron_plotting_v3.py` .

## Arduino commands (text sent as /usb-cmd/<cmd>)
`mc` / `mo` mode (+enable) · `c60.0` duty % · `o50,20,0.2` sine (mo only) · `x` rest · `e` enable · `d` disable

## Telemetry: A1/A2 wiring and ADC reference 
- **A1 ← ESCON AO1 = actual current, A2 ← AO2 = actual speed** (lab test, motor held at 600 rpm in ESCON Studio, compared with the Arduino readings). In `PWM_full_v3.ino` the variables `telemSpeed` / `telemCurrent` were not renamed and are therefore swapped (`telemSpeed` holds A1 = current). The wiring is not detected: it is set in `WIRING` (plotting) and in the default `wiring` of `apply_to_frame` (escon_edc_v3).
- **ADC reference Vref = 4.88 V**, measured with a multimeter between GND and 5V of the lab Arduino (using 5 V gave a ~2.4 % error on every voltage / speed / current). Set in `ADC_VREF` (`chipiron_decode_v3.py`, `chipiron_plotting_v3.py`); 10-bit ADC assumed.
- With the AO-shift `.edc` (analog outputs 0–4 V): current = (V − 2)/0.4 A, speed = (V − 2)/0.002 rpm.
- /!\ **To be re-checked on Chipiron**: the A1/A2 wiring and Vref.

## Speed controller P gain (.edc)
- In the lab `.edc` the mode of controller is "Current Limiter", so ESCON uses the **Speed 1 regulator** (`PARAM_Speed1RegulatorP_Gain`, 0x4201:01); with "Inner Current Control Loop" it would be Speed 2 (0x4202:01). The P gain is a plain integer (no unit). It is **1** in `parameters_bipolar800_EC_ramp2000_offset-11_Oct26_ao_shift.edc`, which is low and can flatten / skew the sine.
- Read it: `py escon_edc_v3.py file.edc --no-save` → block "SPEED CONTROLLER GAINS" (both gains, the one in use is marked).
- Change it (on a copy; the script writes the regulator actually in use and says which):
  `py chipiron_patch_edc_v3.py --in a.edc --out a_gain20.edc --no-ao-shift --speed-p-gain 20`
  (`--no-ao-shift` if the file already has the AO shift; `--force` to overwrite an existing `--out`). 
- Then open the result in ESCON Studio, check the value, write it to the controller. Raise it step by step (e.g. 1, 5, 20, 50): too high gives noise / oscillation.
- To pick a value: `chipiron_plotting_v3.py` step-fit and sine-fit cells give tau (and amplitude ratio / lag) for each gain label.
- /!\ **To be checked on Chipiron**: the gain and mode of controller of its own `.edc` (can differ from the lab one) and the largest gain without noise.

## Warnings
/!\ LORA not tested yet.

## .edc files:
Those are saved configurations to be loaded in Escon Studio Software for a given motor. 

parameters.edc: ESCON motor/controller configuration for Chipiron (DC motor). Saved on 13 december 26. (Stéphane? Bicwin 26?)
parameters_SEP26.edc: ESCON motor/controller configuration for EC motor (PMMH lab). Saved on september 2026 (PMMH lab).
parameters_bipolar800_EC_ramp2000_offset-11_Oct26.edc : ESCON motor/controller configuration for EC motor (PMMH lab): bipolar 10% > -800 rpm, 90% > +800 rpm, offset -11 rpm, ramp 2000 rpm/s. Saved Oct. 2026 (PMMH lab).
parameters_bipolar800_EC_ramp2000_offset-11_Oct26_ao_shift.edc : same, with the analog outputs shifted to 0–4 V (AO1 = current, AO2 = speed, see above); speed P gain 1 (Speed 1 regulator). The -11 rpm offset adjust is stored but NOT enabled in this file (ignored by the scripts, enable it in ESCON Studio if wanted). Patched copies with another gain (`--speed-p-gain`) are named `..._gainN.edc`.
parameters_unipolar800_EC_ramp10_offset-11_Oct26.edc: ESCON motor/controller configuration for EC motor (PMMH lab): unipolar 10% > 0 rpm, 90% > +800 rpm, offset -11 rpm, ramp 10 rpm/s. Saved Oct. 2026 (PMMH lab).
parameters_unipolar800_EC_ramp2000_offset-11_Oct26.edc; ESCON motor/controller configuration for EC motor (PMMH lab): unipolar 10% > 0 rpm, 90% > +800 rpm, offset -11 rpm, ramp 2000 rpm/s. Saved Oct. 2026 (PMMH lab).


## Parallel frazil characterization
solid_fraction_frazil_MausDeLaRosa.py : Solid fraction and salinity of grease ice: Methods 1 to 4 of
Maus & De la Rosa (2012), 'Salinity and solid fraction of frazil and grease ice', J. Glaciol. 58(209), 594-612. Use frazil_synthetic_dataset_test.csv as a test synthetic input (
yet to be tested on real data). 


chipiron_io.py -- shared I/O for the general Chipiron data structure.   NOT TESTED ON REAL h5py / REAL DATA YET.
