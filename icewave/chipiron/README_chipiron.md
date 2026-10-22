# Run & save data  -- Chipiron (19/10/2026 -- Matilde)

Set of python/arduino scripts to run an Escon motor (on Chipiron or else) form laptop/fairphone, and collect data.

Python packages: `numpy`, `pandas`, `matplotlib` (plotting only), `pyserial` (LoRa only).

## Files and who depends on whom
| File | Role | Imports |
|---|---|---|
| `chipiron_helpers_v3.py` | motor config (.edc → mapping, rest duty, ESCON ramp), links http/LoRa, plan CSV, ramp-limiter check | `escon_edc_v3` (lazy) |
| `escon_edc_v3.py` | reads the ESCON `.edc` | – |
| `chipiron_run_v3.py` | holds / ramps / staircases / cycles / CSV, runs the motor | helpers |
| `chipiron_oscillate_v3.py` | sine command (firmware or stream mode) | helpers |
| `chipiron_ramp_planner_v3.py` | plans a staircase ramp, saves CSV for run | helpers |
| `chipiron_session_v3.py` | full measurement: Arduino mode → 0 rpm → record → experiment → save | connect, run_gobannos, escon_edc, helpers; launches run / oscillate |
| `chipiron_patch_edc_v3.py` | makes a patched copy of an `.edc` | – |
| `chipiron_decode_v3.py` | decodes the Usb csv (ADC A1/A2) | – |
| `chipiron_plotting_v3.py` | plots one acquisition (cells `# %%`) | decode (+ escon_edc optional) |
| `connect_v3.py` | phone ID → IP | – |
| `run_gobannos_v3.py` | HTTP calls to the phone | – |
| `mobile_v3.py` | manual menu (status, start/stop, ls, pull, usb, time) | connect, run_gobannos, tsync, rw_data |
| `tsync_v3.py` (+ `rw_data_v3.py`) | clock sync  | connect, rw_data |
| `PWM_full_v3/PWM_full_v3.ino` | Arduino firmware  | – |


## Typical workflow
1. (once) ESCON Studio: load your `.edc` configuration for motor (example: bipolar: 10 % = −800, 50 % = 0, 90 % = +800 rpm). Optionally `chipiron_patch_edc_v3.py` to create a new set of parameters.
2. Upload `PWM_full_v3.ino` to the Arduino Mega (motor stays **disabled** until a mode command arrives).
3. Plan a ramp (optional): ex. `py chipiron_ramp_planner_v3.py --edc parameters_Sep26.edc`
4. Measure:ex. `py chipiron_session_v3.py full_session mo parameters_Sep26.edc`
5. Plot: `py chipiron_plotting_v3.py` .

## Arduino commands (text sent as /usb-cmd/<cmd>)
`mc` / `mo` mode (+enable) · `c60.0` duty % · `o50,20,0.2` sine (mo only) · `x` rest · `e` enable · `d` disable

## Warnings
/!\ LORA not tested yet.

## .edc files:
Those are saved configurations to be loaded in Escon Studio Software for a given motor. 

parameters.edc: ESCON motor/controller configuration for Chipiron (DC motor). Saved on 13 december 26. (Stéphane? Bicwin 26?)
parameters_SEP26.edc: ESCON motor/controller configuration for EC motor (PMMH lab). Saved on september 2026 (PMMH lab).
parameters_bipolar800_EC_ramp2000_offset-11_Oct26.edc : ESCON motor/controller configuration for EC motor (PMMH lab): bipolar 10% > -800 rpm, 90% > +800 rpm, offset -11 rpm, ramp 2000 rpm/s. Saved Oct. 2026 (PMMH lab).
parameters_unipolar800_EC_ramp10_offset-11_Oct26.edc: ESCON motor/controller configuration for EC motor (PMMH lab): unipolar 10% > 0 rpm, 90% > +800 rpm, offset -11 rpm, ramp 10 rpm/s. Saved Oct. 2026 (PMMH lab).
parameters_unipolar800_EC_ramp2000_offset-11_Oct26.edc; ESCON motor/controller configuration for EC motor (PMMH lab): unipolar 10% > 0 rpm, 90% > +800 rpm, offset -11 rpm, ramp 2000 rpm/s. Saved Oct. 2026 (PMMH lab).


## Parallel frazil characterization
solid_fraction_frazil_MausDeLaRosa.py : Solid fraction and salinity of grease ice: Methods 1 to 4 of
Maus & De la Rosa (2012), 'Salinity and solid fraction of frazil and grease ice', J. Glaciol. 58(209), 594-612. Use frazil_synthetic_dataset_test.csv as a test synthetic input (
yet to be tested on real data). 


chipiron_io.py -- shared I/O for the general Chipiron data structure.   NOT TESTED ON REAL h5py / REAL DATA YET.
