# Chipiron: field rheometer (ESCON motor + Arduino + Fairphone/Gobannos)

* updated on 26/11/2026 (Matilde).
*Full explanations of the scripts (pipeline, commands, decoding, sine-skewing debugging): `Running_Chipiron.tex` / `Running_Chipiron.pdf`. Up-to-date scripts: <https://github.com/Turbotice/icewave/tree/main/icewave/chipiron>*

> **Status.** Everything was first tested on an EC ESCON motor at the PMMH lab. The v4 scripts and firmware were then used on Chipiron during the CASCADES26 expedition (Nov. 2026). **LoRa has never been tested.** Items marked `/!\` or listed in [section 12](#12-warnings-limitations-and-things-to-check) must be re-checked on the real set-up.

## Contents
1. [Overview](#1-overview)
2. [Hardware](#2-hardware)
3. [Start-up order and quick manual control (browser / URL)](#3-start-up-order-and-quick-manual-control-browser--url)
4. [Scripts: who does what and who depends on whom](#4-scripts-who-does-what-and-who-depends-on-whom)
5. [Typical workflow](#5-typical-workflow)
6. [Arduino firmware, commands and telemetry](#6-arduino-firmware-commands-and-telemetry)
7. [Saved data and file formats](#7-saved-data-and-file-formats)
8. [ESCON `.edc` files and controller settings](#8-escon-edc-files-and-controller-settings)
9. [Sine skewing and speed-controller gain (lab results)](#9-sine-skewing-and-speed-controller-gain-lab-results)
10. [Post-processing for CASCADES26](#10-post-processing-for-cascades26)
11. [Parallel frazil characterization](#11-parallel-frazil-characterization)
12. [Warnings, limitations and things to check](#12-warnings-limitations-and-things-to-check)
13. [Version history](#13-version-history)

---

## 1. Overview

```
laptop (Python, VS Code) --Wi-Fi/hotspot, HTTP :8080--> Fairphone (Gobannos app = HTTP server)
        --USB--> Arduino Mega 2560 --PWM duty 10-90 % (pin 13)--> ESCON 70/10 --> motor (Chipiron: DC, lab: EC)
                      ^---- A1 / A2 (0-5 V ADC) <-- ESCON analog outputs AO1 / AO2 (speed, current)
```

- The laptop sends HTTP requests to the phone (`/status`, `/start`, `/usb-cmd/<command>`, `/stop`, `/list-files`, `/get-file/<name>`).
- The **Gobannos** app on the phone forwards `<command>` (`mo`, `c60.0`, `o50,20,0.2`, `x`, `e`, `d`, ...) to the Arduino over USB, records the phone sensors (accelerometer, gyroscope, magnetic field, GPS) and logs every byte sent back by the Arduino (the "Usb" file).
- The Arduino (`PWM_full_v4.ino`) outputs a PWM duty cycle (ESCON set value), drives the enable line, and sends back the two analog readings (A1, A2) as 4-byte frames.
- The ESCON turns the duty cycle into a speed set value using the mapping stored in its `.edc`. The **same `.edc` is read by the Python scripts**, so rpm <-> duty conversions, the rest duty (duty giving 0 rpm), the ESCON speed-ramp limits and the AO -> rpm / A conversions are always those of the controller.
- Termux is **not** used for the experiment (Gobannos only); it is kept for diagnosis (SSH, file transfer, see [3.5](#35-optional-ssh-and-termux-diagnosis-only)).

---

## 2. Hardware

### 2.1 Architecture for controlling Chipiron
- **Smartphone** (FP4 used for the lab tests, **now FP5**) running Gobannos.
  - Accessible from a local network knowing its IP (the laptop is connected to the phone's Wi-Fi hotspot; broadcast wifi ?).
  - Accessible from anywhere using ZeroTier / Qrvpn (needs the MAC address for the correspondence).
  - Can be reached by SSH: requires Termux running + phone identifier.
- **USB hub** to
  - connect to an Arduino Mega 2560,
  - power the phone externally (external battery ?).
- **Arduino Mega 2560**
  - pin 13 (PWM): duty cycle from 10 % (0 speed, unipolar) to 90 % (max speed). With the bipolar `.edc`: 10 % = -R, 50 % = 0, 90 % = +R.
  - pin 52 (DO1): ESCON enable (digital input 2); pin 51 (DO2): direction (ignored by the ESCON as configured).
  - A1: speed data from the controller; A2: current data from the controller (**see 2.4: A2 is read on pin A3**).
- **Controller ESCON 70/10**
  - Power 24 V (not very sensitive to the value: should be good without tension control).
  - Wired to the motor.
  - Hall sensor connected (current measurement -> power ? depends on the input power ? **check**).
  - Optical encoder for speed (motor speed without the gear !).
  - Digital IN/OUT: can be configured from the ESCON software (configured on Babasse).
  - Analog IN/OUT: can be configured from the ESCON software (configured on Babasse). *Old note: AO2 = enable (High), AO1 = PWM is now in speed control mode. For the wiring used for the data (AO1 = speed, AO2 = current) see 2.4.*
  - Controlling options, in particular speed range and output ranges, from the ESCON software -> needed some tests in a granular medium.

### 2.2 Motor
- Specifications in the datasheet. **Max current 9 A**. `/!\` We will blow up mechanics and other things well before 9 A: set a **4-5 A current limit** in ESCON Studio.
- Gear ratio: 1184/18 ~ 1/66.
- Chipiron speed range: the first note says max speed is set to **2640 rpm before reduction, i.e. ~40 rpm after the gear** (2640/66 = 40); the `.edc` files are named `bipolar2640`. The 2025/26 configuration of Stéphane is quoted in `Running_Chipiron` as "PWM ratio 2460 rpm <-> 90 %, speed ratio 2460/4 rpm/V, current ratio 30/4 A/V". `/!\` 2640 vs 2460: check the number in the `.edc` actually used (`py escon_edc_v3.py <file>.edc --no-save`).
- Lab motor: EC (brushless), bipolar 10 % = -800 rpm, 50 % = 0, 90 % = +800 rpm (1 % duty = 20 rpm).

### 2.3 Temperature measurement in the Chipiron bucket (Exelog)
- **CH1**: air temperature (type K thermocouple) (doubtful). Distance from probe end/bottom: 40 cm, in the air.
- **CH2**: 30 cm from the bottom, frazil / air limit, surface.
- **CH3**: frazil temperature in the tank, near the surface (type K), frazil, 10 cm from the bottom (between water and frazil).
- **CH4**: bottom of the tank (type K), water, 0 cm.
- **CH5**: Pt100 in the air (never functional).
- An open input reads `High` in the Exelog csv (not numeric: dropped by the loader).

### 2.4 A1 / A2 wiring, ADC reference
- **Chipiron (expedition data): A1 = speed (AO1), A2 = current (AO2).** AO2 is wired on **pin A3** of the Arduino because of dead pins (`PWM_full_v4.ino`: `ai1 = A1`, `ai2 = A3`). Stéphane: speed pair is sent first, then the current pair ~50 ms later (old 250 ms firmware layout, see 6.3).
- **Every `.edc` used for the expedition data has AO1 = actual speed, AO2 = actual current**; checked on 2026-11-08 on the staircase run, the 250 ms sine run and the 10 ms sine run of 2026-11-06 (a 500 rpm commanded sine comes out as 498 rpm on A1, and as 6000+ rpm if A1/A2 are swapped).
- **Lab motor, early tests (20/10/2026, `..._ao_shift.edc`)**: motor held at 600 rpm in ESCON Studio and compared with the Arduino readings gave **A1 <- AO1 = current, A2 <- AO2 = speed**. The v3 firmware variables `telemSpeed` / `telemCurrent` were not renamed and are therefore swapped. This only applies to files recorded with that `.edc` / wiring; the wiring is **not detected**, it depends on the function set on AO1/AO2 in the `.edc` and on the cables. It is set in `WIRING` (plotting, `plot_norm`, expedition loader) and in the `wiring` argument of `apply_to_frame` (`escon_edc_v3`).
- **ADC reference Vref** is the real voltage of the Arduino 5 V pin (depends on how it is powered). Measure it with a multimeter between GND and 5V, **powered as in the experiment**. Lab Arduino: **4.88 V** (using 5 V gave a ~2.4 % error on every voltage, speed and current). Values currently in the scripts: 5.06 V (decoder, plotting, `plot_norm`; "measured on the Amundsen 8/11", real Chipiron), 4.85 V in `CASCADES26_expedition.py` (see warnings). 10-bit ADC assumed.
- With the AO-shift `.edc` (analog outputs 0-4 V, 2 V = zero): current = (V - 2)/0.4 A, speed = (V - 2)/0.002 rpm (for +/-1000 rpm and +/-5 A ranges); V = adc * Vref/1023. One ADC count ~ 2.4 rpm for a 1000 rpm range (~6 rpm for the 2460 rpm Chipiron range).

---

## 3. Start-up order and quick manual control (browser / URL)

*(Replaces "How to use Gobannos to control the motor"; the old `192.168.<network>.<phoneip>` addresses are now `http://<phone-ip>:8080/...`, see `connect_v3.py`.)*

### 3.1 Beware of the order in which you connect the elements
1. Launch Gobannos (phone hotspot on, Gobannos open **in the foreground**; it shows the phone IP).
2. Wire the USB hub.
3. Wire the Arduino (accept the Android USB permission pop-up for Gobannos; do not mix with a Serial USB terminal app).
4. Add the power supply (**plug the battery last**).

### 3.2 Laptop connection check (PowerShell)
```
ping <phone-ip>
Test-NetConnection -ComputerName <phone-ip> -Port 8080     # expect TcpTestSucceeded : True
```
Browser: `http://<phone-ip>:8080/status` should show `STOPPED` (or a running state). Find the phone IP: Gobannos displays it; or in Termux `termux-wifi-connectioninfo` (field `ip`, needs `pkg install termux-api` + Termux:API app); or on the laptop `ipconfig` -> default gateway of the Wi-Fi section ("Passerelle par défaut"). If the USB status line in Gobannos shows "demande de permission USB non envoyée", the phone cannot talk to the Arduino: replug and accept the permission.

### 3.3 Run an acquisition
Use the scripts (`chipiron_session_v6.py`, or `mobile_v3.py`; in the legacy phonefleet repo `mobile.py`) or a browser:
```
http://<phone-ip>:8080/status      (check status)
http://<phone-ip>:8080/start       (start acquisition)
http://<phone-ip>:8080/stop        (stop acquisition)
```

### 3.4 Control the motor speed from a URL
```
http://<phone-ip>:8080/usb-cmd/c60.0
```
`60.0` is the duty cycle in % of the PWM, from **10 (0 speed, unipolar) to 90 (max speed)**. Max speed was set to 2640 rpm before reduction, ~40 rpm after (see 2.2). With a bipolar `.edc` rest is `c50.0`. Full list of commands in [section 6.2](#62-commands-text-sent-as-usb-cmdcmd). Typical manual check (bipolar `.edc`, 1 % = 20 rpm on the lab motor): `/usb-cmd/mo` -> `/start` -> `/usb-cmd/c55.0` (~ +100 rpm) -> `c45.0` (~ -100 rpm) -> `o50,5,0.1` (sine 0 +/- 100 rpm, 10 s period) -> `x` (0 rpm) -> `d` (disabled) -> `e` -> `/stop`. Check on ESCON Studio the PWM value, speed and enable state after each line. Unipolar: `mc`, `c20` (~100 rpm), `c10` (0 rpm). Commands sent from the browser, `curl` or `mobile_v3.py` are **not** written in the input log.

### 3.5 Optional: SSH and Termux (diagnosis only)
Tested: SSH from a Windows laptop to Termux on the Fairphone, browse phone files, locate previously recorded USB serial measurements, transfer them over the local Wi-Fi/hotspot.
- On the phone (Termux): `pkg update`, `pkg install openssh`, `passwd` (PMMH default password), `whoami` (e.g. `u0_a235`), `sshd`; Termux's sshd listens on **port 8022**; check with `pgrep -a sshd`.
- From Windows PowerShell: `ssh -p 8022 u0_a235@<phone-ip>`; home is `/data/data/com.termux/files/home`.
- Shared storage: `termux-setup-storage` (grant permission), then `~/storage/shared` (Download, Documents, DCIM, Pictures, Recordings, RotorlaboData); the Download directory is `~/storage/shared/Download`. This gives access to shared storage only, not to the private files of other apps.
- Search files: `find ~ -type f \( -iname '*serial*' -o -iname '*rotor*' -o -iname '*.csv' -o -iname '*.txt' -o -iname '*.log' \) 2>/dev/null | head -200`
- Copy to the laptop: `scp -P 8022 -r u0_a235@<phone-ip>:storage/shared/Download/Matilde "C:\Users\matil\Downloads\"` (one file: same without `-r`).
- `phonefleet` (<https://github.com/Turbotice/phonefleet>, `pip install -e .` in the folder with `setup.py`) is **no longer needed for the v3/v4 scripts**, only for Termux/SSH debugging.

---

## 4. Scripts: who does what and who depends on whom

Python packages: `numpy`, `pandas`, `matplotlib`, `scipy` (step-fit cell of the plotting script), `openpyxl` (xlsx logs of the expedition loader), `h5py` (`chipiron_io` only), `Pillow` (optional, image time stamps), `uncertainties` (solid fraction), `pyserial` (LoRa only), `nmap` program (only for `connect_v3.py` scan). Everything is run from VS Code (`py script.py`, or cell by cell for the `# %%` scripts).

### 4.1 Running the motor
| File | Role | Imports |
|---|---|---|
| `chipiron_helpers_v4.py` | motor config (`.edc` -> mapping, rest duty, ESCON ramp), http/LoRa links, **input log** of the commands sent, plan CSV, ramp-limiter check | `escon_edc_v3` (lazy) |
| `escon_edc_v3.py` | reads an ESCON `.edc`: motor type, mode, PWM mapping, speed ramp, **speed P gains**, analog-output conversions; prints a summary / saves a `.txt` | - |
| `chipiron_patch_edc_v3.py` | makes a patched **copy** of an `.edc` (AO shift to 0-4 V, bipolar set value, ramp, motor type, offset) | - |
| `chipiron_run_v4.py` | holds / ramps / staircases / up-down cycles / repeats / CSV programs (`--build-only`, `--csv ... --yes`) | helpers |
| `chipiron_ramp_planner_v4.py` | plans a linear / log staircase in rpm, rad/s or %, saves the CSV for `run` | helpers |
| `chipiron_oscillate_v4.py` | sine command (firmware or stream mode) | helpers |
| `chipiron_amp_sweep_v4.py` | sine at fixed frequency, amplitude **stepped** or **continuous** (`r` command), `--check-only` | helpers |
| `chipiron_session_v6.py` | full measurement: phone IP + sampling period -> `.edc` -> Arduino mode -> 0 rpm + wait -> record -> experiment -> rest -> stop -> ls -> save (+ input log, `edc_used.txt`, period file) | connect, run_gobannos, escon_edc, helpers, decode; launches run / oscillate / amp_sweep |
| `PWM_full_v4/PWM_full_v4.ino` | Arduino firmware (modes mc/mo, `c`, `o`, `r`, `x`, `e`, `d`, telemetry with start-of-frame flag, tunable rate) | - |
| `mobile_v3.py` | manual menu (network, phones, status, start/stop, ls, pull, usb, time) | connect, run_gobannos, tsync, rw_data |
| `connect_v3.py` | phone ID -> IP (`IP_PREFIX`, `NETWORK`, `DEFAULT_PHONELIST`); `scan()` with nmap. Used by `mobile_v3`, `tsync_v3` and the session script for its helper functions | - |
| `run_gobannos_v3.py` | HTTP calls to the phone (status, start, stop, ls, get-file, usb-cmd, clean) | - |
| `tsync_v3.py` (+ `rw_data_v3.py`) | UDP clock synchronisation laptop <-> phone (ports 5000/5001); csv writer | connect, rw_data |

### 4.2 Decoding and plotting
| File | Role | Imports |
|---|---|---|
| `chipiron_decode_v4.py` | decodes the Usb csv (A1/A2 ADC) for **any** telemetry period; framings `flag` (v4 firmware), `nibble` (flagless fast data), `gap` (old 250 ms firmware); `TELEM_PERIOD_MS` must equal the flashed firmware's | - |
| `Chipiron_plotting_v5.py` | plots one or several acquisitions (cells `# %%`), **split recordings concatenated**, decoded ADC / volts / rpm and A, first-order step fit | decode (+ escon_edc optional) |

### 4.3 Post-processing (CASCADES26)
| File | Role | Imports |
|---|---|---|
| `scan_chipiron_experiments.py` | catalogue of all tests from the input logs (`data_summary/chipiron_exp_summary.csv`) + check that every time folder has an input csv | - |
| `CASCADES26_expedition.py` | loads every test into `Experiment` objects (metadata, test type, usb, temperature, images, xlsx cross-checks, **reference matching**) + plotting cells | match_temperature_usb, helpers, escon_edc, decode (lazy) |
| `chipiron_files_debug.py` | **data audit**: scans the whole tree with the expedition loader and lists what to check (report + csv, `Exceptions.xlsx` support) | CASCADES26_expedition, escon_edc |
| `match_temperature_usb.py` | matches Exelog temperature files to Usb acquisitions, writes `usb_temperature_match.csv` | - |
| `plot_norm_111126.py` | sine or ramp tests in frazil, alone or against air (`--test sine\|ramp_test --mode no_ref\|ref`) | decode, escon_edc |
| `quick_data_visu.py` | Beamer overview slides (A1, A2, A2 vs A1 per acquisition) + zip for Overleaf | - |
| `solid_fraction_frazil_MausDeLaRosa_v1.py` | solid fraction and salinity of grease ice (Methods 0-4) | chipiron_io (optional) |
| `plot_phi_cascades.py` | plots phi_s from the results files | - |
| `chipiron_io.py` | shared I/O (metadata txt, HDF5 rheology, sample csv). **Not tested on real h5py / real data yet** | - |

---

## 5. Typical workflow

1. **(once) ESCON Studio**: load your `.edc` (bipolar example: 10 % = -R, 50 % = 0, 90 % = +R rpm). Optionally `chipiron_patch_edc_v3.py` to create a new set of parameters. Check **which AO gives speed and which gives current**. Keep the controller disabled while loading; double-check every parameter in the GUI. Read an `.edc` first with `py escon_edc_v3.py file.edc --no-save`.
2. **Arduino**: upload `PWM_full_v4.ino` to the Mega (board Arduino Mega 2560, from the Arduino IDE on a laptop, once) and **note `TELEM_PERIOD_MS`** (the motor stays **disabled** until a mode command arrives). Measure the ADC reference `Vref` (GND to 5V, powered as in the experiment). Then unplug from the laptop and plug into the phone.
3. **Phone**: hotspot on, Gobannos in the foreground, plug the Arduino, accept USB permission (see [3.1](#31-beware-of-the-order-in-which-you-connect-the-elements)). Laptop on the hotspot, connection check (3.2).
4. **Plan a ramp (optional)**: `py chipiron_ramp_planner_v4.py --edc <file>.edc`
5. **Measure**: `py chipiron_session_v6.py full_session mo <file>.edc`
   It asks the **full phone IP** (default: last one) and the **telemetry sampling period** (must equal `TELEM_PERIOD_MS` of the flashed firmware), checks the Arduino mode against the `.edc` (rest duty must give 0 rpm), then the experiment: **1 manual, 2 sine, 3 program (holds/ramps/steps/CSV), 4 other script, 5 amplitude sweep**. All questions are asked **before** the acquisition; then 0 rpm command + `ZERO_WAIT_S` = 30 s wait, `/start`, experiment (scripts are called with `--no-record`), motor back to rest in any case (also on Ctrl-C), `/stop`, `ls`, save.
   Variants: `py chipiron_session_v6.py` (asks the mode), `save_mode` (no experiment: ls + save), `full_session mc <file>.edc` (unipolar), `full_session skip` (send no mode: after a reset the motor stays disabled until `e`, `mo` or `mc`).
6. **Plot**: `py Chipiron_plotting_v5.py` (answer 250 as period only for old 250 ms files without `sampling_period_ms.txt`).
7. **Post-process**: `py scan_chipiron_experiments.py <Phonefleet_data>`, `py match_temperature_usb.py`, `load_expedition()`, `py chipiron_files_debug.py`, `py plot_norm_111126.py`, ... (section 10).
8. **End of the test**: motor at rest on the ESCON display, no error, phone status `STOPPED`; if the motor is left for a long time send `d` or switch the supply off; keep the program CSV (`programs/session_<date>.csv`) next to the data; push modified scripts to git; save and label the ESCON parameters.

### 5.1 Experiment types and which Arduino mode / `.edc` they need
| Test | Arduino mode | ESCON mapping | Script run by the session | Rest duty |
|---|---|---|---|---|
| Oscillation, firmware mode (smooth sine, built by the Arduino at 100 Hz) | `mo` | bipolar | `chipiron_oscillate_v4.py --mode firmware` | 50 % |
| Oscillation, stream mode (coarse, not recommended; only for periods >~ 20 s) | `mo` or `mc` | bipolar or unipolar | `chipiron_oscillate_v4.py --mode stream` | from `.edc` |
| Constant speed, ramps, steps, up/down cycles, ramp-planner CSV | `mo` or `mc` | bipolar or unipolar | `chipiron_run_v4.py` | from `.edc` |
| Amplitude sweep at fixed frequency (stepped / continuous) | `mo` | bipolar | `chipiron_amp_sweep_v4.py` (continuous needs v4 firmware) | 50 % |
| Manual (URL / keyboard commands, any other script) | any | any | none or your script | from `.edc` |

The session stops with an ERROR if the rest duty of the Arduino mode does not match the `.edc`, and refuses firmware-mode oscillation with `mc` or with a unipolar `.edc`.

### 5.2 Useful command examples
```
py escon_edc_v3.py parameters_Sep26.edc --folder edc_summaries          # summary saved as .txt
py chipiron_patch_edc_v3.py --in a.edc --out b.edc --bipolar-setpoint 800 --ramp 2000 --motor-type DC --offset -11
py chipiron_oscillate_v4.py --edc a.edc --mean-rpm 0 --amp-rpm 100 --freq 0.1 --cycles 5
py chipiron_run_v4.py --build-only programs/test1.csv --rest-before 30 --rest-after 10
py chipiron_run_v4.py --csv programs/test1.csv --yes
py chipiron_ramp_planner_v4.py --unit rpm --type log --start 10 --end 500 --points 12 --time 60 --folder ramps
py chipiron_amp_sweep_v4.py --freq 0.1 --sweep 50 600 --log --sweep-time 600 --no-record --yes   # continuous
py chipiron_amp_sweep_v4.py --freq 0.1 --amps 100 200 300 400 --cycles 10 --no-record --yes      # stepped
py tsync_v3.py --phone 36 --n 100 --iter 3 --folder Phonefleet/Gobannos_Tsync
```
Standalone `run` / `oscillate` / `amp_sweep` ask for the `.edc` (or `--edc`), use the phone IP `--phone-ip` > `CHIPIRON_PHONE_IP` > default `10.153.247.36`, and call `/start` and `/stop` themselves unless `--no-record`. Standalone they do **not** send the rest command at the beginning: send `c50.0` (or `mo`) yourself first. The session hands the `.edc` configuration and phone IP to its sub-processes through `CHIPIRON_MOTOR_CFG`, `CHIPIRON_PHONE_IP` and the command log through `CHIPIRON_INPUT_LOG`.

### 5.3 Safety and pre-checks (before any motor test)
Motor unloaded or free to turn, **a hand on the supply switch**, small speeds first (a few hundred rpm), current limit set in ESCON Studio (4-5 A suggested), **battery plugged last**. `d` (disable) makes the motor coast and is **not an emergency stop**. The `.edc` loaded in ESCON Studio must be the one given to the scripts (nothing verifies it). Fault tests (see `Running_Chipiron`): `mc` + bipolar `.edc` must stop with an error; firmware oscillation with `mc` or amplitude beyond +/-R must error before anything starts; unplugging the Arduino must make commands fail visibly (final rest prints "STOP THE MOTOR MANUALLY"); after an Arduino reset (replug) the motor must be disabled until `mo`, `mc` or `e`.

### 5.4 Program builder blocks (`chipiron_run_v4.py`, manual mode)
Values in rpm or %, `u` switches the unit; every block can be repeated N times with an optional pause at rest; every value is checked against 10-90 % duty: 1 hold, 2 ramp (by duration or rate, warns if shorter than the ESCON ramp allows), 3 staircase (typed list / linear / log, optional mirror, same time or one by one, optional rest between plateaus), 4 up/down cycle(s), 5 repeat the whole program, 6 load CSV (ramp planner or saved program, mapping checked), 7 remove last block, 9 done (rest before/after, FINAL PROGRAM summary with ESCON warnings, "Save as CSV?", "Run it now?"). Ramps are staircases (`--dt >= 1 s` recommended). The ramp-planner CSV header stores the mapping so `run` refuses a CSV planned with another `.edc`; columns `t_start_s, duration_s, duty_pct, rpm, omega_rad_s, kind, segment`.

### 5.5 Sine checks (before anything moves)
mean +/- amplitude must fit in 10-90 % duty; the peak acceleration `A * 2*pi*f` is compared with the ESCON speed ramp (above it the limiter flattens the sine: warning with the largest usable amplitude `a_ramp/(2*pi*f)` and a confirmation). Rule of thumb: rad/s -> rpm = x 60/2pi ~ x 10. The ESCON speed-controller gain is a second, independent limit (section 9).

---

## 6. Arduino firmware, commands and telemetry

### 6.1 Firmware `PWM_full_v4.ino` (Arduino Mega 2560)
- Two running modes chosen at run time by a text command (no re-upload): `mc` constant (rest duty 10 %, unipolar ESCON config) and `mo` oscillation (rest duty 50 %, bipolar config). **The ESCON configuration must match the mode.**
- **Safe start**: at power-up / reset the motor is **disabled** (`START_ENABLED 0`), default mode is oscillation (bipolar). Sending `mc` or `mo` selects the mode, writes its rest duty and *then* enables the motor ("arming"); `e` enables without changing mode. (Before, the Arduino booted in constant mode with the motor enabled: with a bipolar `.edc` the 10 % rest duty meant full negative speed.)
- The Arduino resets every time it is re-plugged or its USB port is re-opened: send the mode again at each session (the session script does it).
- Switches at the top of the file: `DEFAULT_MODE`, `START_ENABLED`, `REST_DUTY_CONST` (10), `REST_DUTY_OSC` (50), `USE_LORA` (0, untested), `TELEM_PERIOD_MS`.
- v4 novelties: command `r` (continuous amplitude sweep), command buffers enlarged from 32 to 64 characters, tunable telemetry rate, start-of-frame flag. Duty updated at 100 Hz; `analogWrite` is 8 bit -> duty steps of 100/255 = 0.39 % (~7.8 rpm with 20 rpm/%).
- Firmware silently clamps duty, mean and amplitude to 10-90 %; the sine starts at the mean going up (no phase control).

### 6.2 Commands (text sent as `/usb-cmd/<cmd>`)
| Command | Mode | Effect |
|---|---|---|
| `mc` / `mo` | any | select mode, go to rest duty (10 / 50 %), **enable** the motor |
| `c<duty>` | both | constant duty in %, decimals accepted (`c50`, `c60.5`); with the bipolar lab `.edc`: `c50.0` = 0 rpm, `c60.0` = +200 rpm, `c40.0` = -200 rpm |
| `o<mean>,<amp>,<freq>` | `mo` only | sine on the duty: `o50,20,0.5` = 50 % +/- 20 % at 0.5 Hz |
| `r<mean>,<amp0>,<amp1>,<freq>,<T>[,<shape>]` | `mo` only | continuous amplitude sweep: amplitude from amp0 to amp1 (% duty) in T s, shape 0 linear / 1 geometric (needs amp0, amp1 > 0); after T the sine stays at amp1 until `x` / `c` / `o`; e.g. `r50.00,2.50,30.00,0.1000,600.0,1` |
| `x` | both | back to the rest duty of the current mode |
| `e` | both | enable the ESCON (enable line high) |
| `d` | both | rest duty, then disable the ESCON: motor coasts, **no braking** (send `x` and wait first for a controlled stop) |

Commas and dots stay literal in the URL (the phone server does not decode `%2C`; the firmware also tolerates `%2C`). Two commands arriving within ~5 ms can be merged. Duty values are clamped to 10-90 % by the Arduino (and refused by the scripts outside this range); hardware resolution 0.39 %. Channels to send a command: browser URL, `curl.exe -v --max-time 5 http://<ip>:8080/usb-cmd/mo`, `mobile_v3.py` (`usb mo`, sent to all active phones), the session manual prompt (`e`, `d`, `x`, `c60.0`, `o...`, `r...`), the scripts, or your own Python (`chipiron_helpers_v4.make_link(...).send("c60.0")` / `run_gobannos_v3.usb_cmd(ip, "x")`). The Arduino **never acknowledges**: HTTP 200 only means the phone got the command.

### 6.3 Telemetry (Arduino -> phone) and decoding
- **v4 firmware**: one 4-byte frame every `TELEM_PERIOD_MS` (250 = old 4 Hz, 10 = 100 Hz, 5 = 200 Hz; minimum 2 ms): `b1 = ADC(A1) & 0xFF`, `b2 = (ADC(A1) >> 4) | 0x80`, `b3 = ADC(A2) & 0xFF`, `b4 = ADC(A2) >> 4`. **Bit 7 of b2 is a start-of-frame flag**, so the PC finds the frames from the bytes alone for any period. Decode: `adc = ((b2 & 0x7F) << 4) | (b1 & 0x0F)`, check `(b1 >> 4) == (b2 & 0x0F)` (a failed check gives NaN). Use the **same number** in the session (asked at start) and `chipiron_decode_v4.TELEM_PERIOD_MS`.
- **Old firmware (v3, 250 ms, no flag)**: bytes at ~0 / 10 / 50 / 60 ms of each cycle, 4 bytes every 250 ms; b2 = adc >> 4 without flag. Recognised automatically by the decoder (framing `gap`: cycles found from time gaps > 100 ms; cycles with only 2 rows ~50 ms apart are `low8`: only the low 8 bits known, ADC counts NaN, aliased).
- Flagless fast data (framing `nibble`): A1/A2 order is ambiguous from the bytes alone; it is guessed (A1 = speed = the smoother channel) with a warning; flash `PWM_full_v4.ino` (flag) to avoid it. `phase=0|2` forces the other alignment.
- Time axis: `phone` timestamps, or `uniform` (frame number x period) when the phone stamps bytes in batches (chosen automatically at high rates; assumes no frame was lost, the number of lost frames is estimated and printed).
- Why two bytes: the ADC is 10 bit (0-1023, "12 bits" in an old note was wrong) but USB sends bytes; bits 4-7 are sent twice (consistency check). Worked example: adc 700 -> b1 = 188, b2 = 43; 43*16 + (188 & 15 = 12) = 700.
- Row parity is **not** a byte role: if the recording starts mid-cycle or a byte is lost, "even rows" become b2/b4. Use "bytes by role" in the plotting script, not even/odd.
- Conversion chain: voltage on pin -> ADC (0-1023) -> 2 bytes -> USB -> Usb csv -> decode -> `adc * Vref / 1023` = V -> `.edc` AO conversion -> rpm / A.
- Old (v3-era) description for reference: speed comes first, then current after ~50 ms, waiting ~150 ms+ for the next loop (timestamps `t0 speed`, `t0+50ms current`, `t1 = t0+~200 ms speed`); values following within < 50 ms are the bytes of the same output and need post-processing to rebuild the value.

---

## 7. Saved data and file formats

```
Phonefleet_data/<date>/<time>/
    <date>T..-Usb-1-...csv (Usb-2... if split), accelerometer / gyroscope / magnetic_field / gps csv
    input_<date>_<time>.csv      commands sent (laptop clock)
    edc_used.txt                 .edc used (line 1 = file name, line 2 = full path)
    sampling_period_ms.txt       telemetry period
    figures/                     saved plots
Phonefleet_data/<date>/metadataDDMMYYYY.txt, EXELOG*/, grain_size/, sampling_height/   (day files)
Phonefleet_data/LOG_days_media_Chipiron.xlsx    (day -> medium + Notes time ranges)
Sorties_UMA/LOG_UMA_Chipiron.xlsx               (station -> days processed by Chipiron)
Phonefleet_data/data_summary/chipiron_exp_summary.csv   (written by scan_chipiron_experiments.py)
```
- `<date>` = `YYYY-MM-DD`, `<time>` = `HH-MM-SS`, read from the Gobannos file names; the session saves them under the root folder (default `Phonefleet_data`). `mobile_v3.py pull` saves files as `P<phone>_D<original name>` in `Phonefleet/`.
- **Long recordings are split by the phone** (`...Usb-1-...`, `...Usb-2-...`; same for the other sensors); v5 plotting and the expedition loader concatenate all parts.
- Usb csv: no header; timestamp in microseconds (phone clock), byte as signed int8 (-43 = 213), empty last column; rows with timestamp <= 0 are placeholders (pre-allocated phone buffers) and removed. Other sensors: accelerometer `t_us, ax, ay, az` (m/s2, Android convention), gyroscope `gx, gy, gz` (rad/s), magnetic field `mx, my, mz` (microtesla), gps 4 columns (meaning unknown, empty so far). `t_s` = seconds since the first sample **of that file**: different files are not on a common time axis. Still to confirm: the timestamp clock and the GPS columns.
- **Input log** `input_<date>_<time>.csv`: columns `t_unix_s, datetime, command, ok, source, note`; every command sent by the session / scripts (not by browser / curl / `mobile_v3`), plus `[acquisition /start]` and `[acquisition /stop]`. The laptop clock is **not** aligned with the phone clock automatically (`tsync_v3.py` optional, not integrated): align on the start of the motion. Keep `programs/*.csv` with the data.
- Decoded Usb columns: `t_us, kind (full|low8), A1_lo, A2_lo, A1_adc, A2_adc, b1..b4, t_s, A1_V, A2_V`, plus `A1_speed_rpm, A2_current_A` after `apply_to_frame` with an `.edc`.
- `chipiron_io` data structure (not tested on real data): `<stem>_meta.txt` (hub: `key = value | unit | uncertainty`, `[links]` section naming the files), `<stem>_rheo.h5` (raw rheology, attrs `unit`/`uncertainty`), `<sample>.csv` (sample characterisation, long format).

---

## 8. ESCON `.edc` files and controller settings

Saved configurations to be loaded in ESCON Studio for a given motor. **Always give the scripts the `.edc` that is loaded in the controller** (nothing verifies it; its name is saved as `edc_used.txt`).

- `parameters.edc`: ESCON motor/controller configuration for Chipiron (DC motor). Saved on 13 december 26. (Stéphane? Bicwin 26?) *(date to check: it is later than today)*
- `parameters_SEP26.edc`: EC motor (PMMH lab), saved September 2026.
- `parameters_bipolar800_EC_ramp2000_offset-11_Oct26.edc`: EC motor (PMMH lab): bipolar 10 % > -800 rpm, 90 % > +800 rpm, offset -11 rpm, ramp 2000 rpm/s. Oct. 2026.
- `parameters_bipolar800_EC_ramp2000_offset-11_Oct26_ao_shift.edc`: same, with the analog outputs shifted to 0-4 V (**lab: AO1 = current, AO2 = speed**), speed P gain 1 (Speed 1 regulator). The -11 rpm offset adjust is stored but **not enabled** in this file (ignored by the scripts; enable it in ESCON Studio if wanted). Patched copies with another gain are named `..._gainN.edc`.
- `parameters_unipolar800_EC_ramp10_offset-11_Oct26.edc`: EC motor (PMMH lab): unipolar 10 % > 0 rpm, 90 % > +800 rpm, offset -11 rpm, ramp 10 rpm/s. Oct. 2026.
- `parameters_unipolar800_EC_ramp2000_offset-11_Oct26.edc`: EC motor (PMMH lab): unipolar 10 % > 0 rpm, 90 % > +800 rpm, offset -11 rpm, ramp 2000 rpm/s. Oct. 2026.
- `parameters_Chipiron_bipolar2640_...edc` and `parameters_Chipiron_unipolar_...edc`: Chipiron configurations used for the expedition data (**AO1 = actual speed, AO2 = actual current**). *(add exact file names, dates and ranges)*

### 8.1 Reading and patching
- Read: `py escon_edc_v3.py file.edc --no-save` prints motor type, operation mode, PWM mapping (10/50/90 %, duty giving 0 rpm, polarity), speed ramp, speed controller gains (the regulator in use is marked), AO conversions and warnings (e.g. offset adjust stored but not enabled).
- Patch a copy with `chipiron_patch_edc_v3.py` (`--in`, `--out`, `--force`): `--ao-shift / --no-ao-shift` (default on; speed -R..+R rpm -> 0-4 V with `--speed-range`, current with `--current-range`; **required with the bipolar mapping**, because a negative voltage reads as 0 on the Arduino and below about -0.5 V is outside the ATmega ratings), `--bipolar-setpoint R`, `--ramp RPM_PER_S`, `--motor-type EC|DC` (written only when given; a wrong value can make the controller drive the motor wrongly), `--offset RPM` (only acts if enabled in ESCON Studio). **Always open the result in ESCON Studio and check every value before writing it to the controller.** *(The `--speed-p-gain` option described in the Oct. 2026 README / `Running_Chipiron` is not in the `chipiron_patch_edc_v3.py` version currently in the project: see section 12.)*
- **Offset**: should be tuned so that the rest command (`c50.0` bipolar) gives ~0 rpm: tune by hand in ESCON Studio, store with `--offset`, check that offset adjust is enabled. The offset is recorded and reported but **not** added to the PWM mapping (rest duty stays 10 % unipolar / 50 % bipolar).
- If in ESCON Studio the "PWM value" is stuck at an out-of-range value and the motor is in disabled/error mode, check that the Arduino is plugged in and running (a reset Arduino boots with the motor disabled until `mo`, `mc` or `e`).
- Analog-output functions: AO1/AO2 function set in `PARAM_AnalogOutputsConfiguration` (lab file 0x42: AO1 = function 2 = Actual Current, AO2 = function 4 = Actual Speed; averaged versions are functions 3 and 5).

### 8.2 Speed controller P gain
- In the lab `.edc` the mode of controller is "Current Limiter", so ESCON uses the **Speed 1 regulator** (`PARAM_Speed1RegulatorP_Gain`, 0x4201:01); with "Inner Current Control Loop" it would be Speed 2 (0x4202:01). The P gain is a plain integer (no unit). It is **1** in `..._Oct26_ao_shift.edc` (ESCON default is 200), which is low and **skews the sine** (section 9). Integral time 11.2 ms in the lab file.
- Raise it step by step (e.g. 1, 5, 20, 50): too high gives noise / oscillation. Then open the result in ESCON Studio, check the value, write it to the controller.
- `/!\` To be checked on Chipiron: the gain and mode of controller of its own `.edc` (can differ from the lab one) and the largest gain without noise.

---

## 9. Sine skewing and speed-controller gain (lab results)

Lab motor and load, 21-22 October 2026 (details and full test log in `Running_Chipiron`, section 4).
- A commanded sine 0 +/- 600 rpm at 0.1 Hz came out with ~200 rpm amplitude (also on the live ESCON display: not a Python/plotting artefact). Cause: speed P gain = 1, so the closed speed loop is a first-order low-pass with tau ~ 4.6 s (corner 0.034 Hz). Neither the PWM mapping, the current limit (< 1 A) nor the ramp.
- Model: `tau * dy/dt + y = K u`; sine amplitude ratio `1/sqrt(1 + (2 pi f tau)^2)`, phase lag `atan(2 pi f tau)`. `tau * gain ~ 4.2-4.9 s` for gains 1 to 20 (tau ∝ 1/gain). At gain 20 (tau ~ 0.25 s): 0.1 Hz sine x0.99, 0.5 Hz sine x0.79 with ~38 deg lag, as predicted from the step.
- For >= 90 % amplitude: `tau <= 0.077/f` (gain ~ 6 at 0.1 Hz, ~18 at 0.3 Hz, ~30 at 0.5 Hz, ~60 at 1 Hz: rows above gain 20 are extrapolations, to be confirmed with sines).
- Peak acceleration `A * 2 pi f` must also stay below the ESCON ramp (e.g. 600 rpm at 1 Hz = 3770 rpm/s > 2000: reduce to <= 318 rpm or raise the ramp).
- Telemetry at 4 Hz cannot resolve tau below ~0.5 s (gains 50 and 100 were not characterised); with the 100 Hz firmware tau >~ 0.02 s is accessible. Use sines or the ESCON data recorder (~20 ms sampling) for fast loops.
- **Raw speed is quantised in steps of ~78 rpm** (probably one encoder edge per 0.187 ms control cycle; hypothesis, not verified; not the ADC, 1 count = 2.4 rpm). Use the **ESCON "averaged" analog outputs for the data** (AO1 = Actual Current Averaged, function 3; AO2 = Actual Speed Averaged, function 5, in the lab wiring): same mean, same dynamics (no measurable extra lag, <~ 10-20 ms at 0.5 Hz), noise divided by ~4. It cannot be done afterwards by averaging the 4 Hz samples. Keep raw outputs for diagnosing the loop itself.
- **Low-amplitude dead zone** (lab motor): with +/-1.5 % duty (+/-30 rpm) the speed never changed sign; reversal from about +/-4 % duty (~80 rpm). To be checked on Chipiron.
- Procedure for the real motor and load (summary): fix fmax and amplitude; set the measurement chain first (A1/A2, Vref, AO shift, averaged outputs, note the `.edc`); calibrate the load with a hold step at gain 1 then 5, fit tau; raise the gain by factors 2-3 (steps while tau >~ 0.5 s, sines at fmax and lower frequencies after); stop at the first gain meeting the target, add a margin x1.5-2; stop earlier on overshoot / ringing / buzzing / current near limit; adjust integral time only if overshoot; one parameter at a time; validate with the real protocol and **use the measured speed, not the commanded one, in the rheology analysis**; write gain, Ti, ramp, current limit, `.edc` name, raw or averaged outputs, load and sample in the experiment metadata. ESCON auto-tuning moves the motor: only with a free rotor.
- Plotting: `Chipiron_plotting_v5.py` has the step-fit cell (tau, rise time, corner frequency, predicted sine gain/lag, tau-vs-gain figure, needs scipy). The sine-fit cell described in the PDF was in `chipiron_plotting_v3.py` (section 12).
- Open points: can mean/raw files be exported at the same time through ESCON ? Is the 2000 rpm/s ramp active (steps at gain 50-100 are compatible with a ramp limit; test a 100-200 rpm step at high gain) ? Exact filter of the averaged outputs; origin of the 78 rpm quantisation; behaviour above gain 20 and with the real load.

### 9.1 Physical target for oscillatory rheology
The target is a sinusoidal **strain** `gamma(t) = gamma0 sin(2 pi f t)`, not a sinusoidal speed. If `gamma = K theta` (vane angle, to be calibrated on the real geometry), the angular speed is a **cosine** of amplitude `theta_dot_0 = 2 pi f gamma0 / K`, i.e. `rpm0 = 60 f gamma0 / K`. The scripts command `mean + amp sin(2 pi f t)` (firmware sine starts at the mean going up; no phase option), so a strain-controlled test is **not yet possible**.

---

## 10. Post-processing for CASCADES26

### 10.1 Experiment loading and cross-checks (`CASCADES26_expedition.py`)
- One experiment = one `<date>/<HH-MM-SS>` folder = one `Experiment` object: identity (date, time, station), day metadata (`metadataDDMMYYYY.txt`: medium, blades, temperatures, heights), test description read from the input log (type: sine, sweep, sweep_stepped, ramp, steps, hold, rest; direction, frequency, amplitudes, levels...), lazy raw data (usb raw / decoded / physical, sensors), matched Exelog temperature, images (grain size, sampling height; processed `vertical_and_grain_size_day<xxxxxx>_station<xxx>.csv`).
- Media: `air-empty`, `saltwater`, `frazil`, `labgrown_frazil`. The two xlsx day logs are **cross-checked** with the metadata txt (every difference raises a "CROSS-CHECK" warning and is listed in the summary and `exp.issues`). A day may hold several media/stations: the Notes column of the media log gives the time ranges (e.g. "frazil 101: from 10am to 11:30am. Empty; from 1:30pm to 2:50pm.") and each test gets the medium and station of the range containing its start time. Overrides: `MEDIUM_OVERRIDES`, `STATION_OVERRIDES`.
- Usage: `exps = load_expedition()`, `exps.table()`, `exps.filter(medium="frazil", test_type="sine", freq_hz=between(0.1, 1))`, `e.usb(physical=True)`, `e.temperature()`, `e.inputs`. Spyder-style cells: LOAD (console log `console_log_*.txt`, reference report, decodability check), PLOT sections, DEBUG match plots.
- **Reference days** (`REFERENCE_DAYS`): **13/11 = empty motor (`empty_ref`), 14/11 = labgrown slush (`ninja_ref`), 16/11 = salt water (`saltwater_ref`)**. Every test of those days is a reference test with a `test_key` built from its input log; every other test gets `ref_available = {ref_type: [matching reference tests]}` (same type and parameters within a tight tolerance of 1 duty step; usable-usb references first). Plot modes: `no_ref | empty_ref | ninja_ref | saltwater_ref`; `exps.pairs("empty_ref")` gives `[(test, ref or None)]`.
- `usb_usable(e)` / `decodability_report(E)`: list every test whose usb cannot be decoded and which references that affects.

### 10.2 Data audit (`chipiron_files_debug.py`)
`py chipiron_files_debug.py [--root ...] [--days ...] [--no-decode] [--only-errors]` reuses the loader and writes `audit_report_<stamp>.txt`, `audit_issues_<stamp>.csv`, `audit_tests_<stamp>.csv`, `audit_edc_<stamp>.csv` in `<root>/audit`. Severities ERROR / WARN / INFO; checks include missing/empty/short/undecodable usb, gaps, ADC saturation, missing edc, wrong AO wiring, commanded speed beyond the AO range, rest offset, plateau/peak ratio vs the command (wrong `.edc` scaling or Vref), medium/station unknown, missing temperature, unusable references, day-level cross-checks. `Exceptions.xlsx` (next to the script; one column per date, cells `HH:MM:SS -- instruction`) suppresses named checks or skips a test/day. The report ends with a **manual checklist** (phone files complete, what ESCON Studio really had loaded, `.edc` timeline per day, ADC_VREF, metadata typos, Exelog clock vs laptop clock, images, reference days).

### 10.3 Other scripts
- `scan_chipiron_experiments.py "<Phonefleet_data>"`: groups tests by their parameters then by date into `data_summary/chipiron_exp_summary.csv` (one row per test type x date: category, test, mean/amplitude rpm, frequency, date, nb_runs, hours, comments, input_files) and warns about time folders without input csv.
- `match_temperature_usb.py`: matches Exelog files (`yy-mm-dd/DA-HH-MM.csv`, 1 min resolution) to Usb acquisitions (pass 1 start match within `START_THRESHOLD_S` = 300 s, pass 2 overlap), lists unmatched ones, writes `usb_temperature_match.csv`, plots/saves selected temperatures. `TEMP_CLOCK_OFFSET_S` corrects a shifted logger clock.
- `plot_norm_111126.py`: **sine** (frazil run, + air run in `ref` mode: figs speed/current vs time, current vs speed loop, and fig 4 normalisation by the air: ratio, diff, diff_scaled, diff_speed, regression) and **ramp_test** (first ramp and first opposite ramp of a `chipiron_run` program, up/down in two colours; ref mode: raw time-matched frazil-minus-air current, no filtering/binning/averaging). Needs the `.edc` for rpm and A.
- `quick_data_visu.py`: one Beamer slide per acquisition (A1 vs t, A2 vs t, A2 vs A1) -> `slides_overview/slides.tex` + zip for Overleaf.

### 10.4 Stations processed with Chipiron
101 (little), 102 (little), 103, 202, 301, 400.

### 10.5 Frazil tank
30 x 30 x 35.5 cm.

---

## 11. Parallel frazil characterization

`solid_fraction_frazil_MausDeLaRosa_v1.py` (first version `solid_fraction_frazil_MausDeLaRosa.py`, tested with the synthetic input `frazil_synthetic_dataset_test.csv`): solid fraction and salinity of grease ice, Methods 1 to 4 of Maus & De la Rosa (2012), "Salinity and solid fraction of frazil and grease ice", J. Glaciol. 58(209), 594-612, plus Method 0 (PSS-78, conductivity + temperature -> salinity).
- Input: one raw CSV row per **bucket**; rows with the same date + time are the drained / undrained buckets of one sample (station checked, note decides drained/undrained); optional `X_unc` uncertainty columns.
- Each sample is computed in **"direct"** (salinities in psu given) and **"conductivity"** (PSS-78) mode (`conductivity_ref_T = 25` for a temperature-compensated meter, `None` for raw conductivity at the real temperature); one results csv per sample and mode (`grease_ice_results_<stamp>_<mode>.csv`, long format) + a summary table with the reason for every skipped method. Air correction per method (`include_air`, default phi_a = 0.01 +/- 0.005).
- Run on the first lab data (3/11/2026, station 101) and on the CASCADES26 samples. `plot_phi_cascades.py` plots phi_s (figure 1: one subplot per method, conductivity vs salinity-based; figure 2: all methods, top conductivity, bottom salinity; x = day / station / sample / datetime; colour = characterization iteration read from the raw CSV notes).
- `chipiron_io.py`: shared I/O for the general Chipiron data structure. **NOT TESTED ON REAL h5py / REAL DATA YET.**

---

## 12. Warnings, limitations and things to check

**Check on the real Chipiron**
- The imposed current does not depend on the exact battery voltage; the motor range is properly configured for the tests; we can decode and re-code the USB data properly.
- A1/A2 wiring (speed/current; AO2 on pin A3), Vref (multimeter GND-5V, powered as in the experiment), 2640 vs 2460 rpm range of the `.edc`.
- The range where the sine is not skewed (amplitude/frequency domain; ESCON speed ramp vs `A 2 pi f`; current limits); low-amplitude dead zone.
- If the ramp acceleration is fixed to a low 10 rpm/s value (slow `.edc`), spend enough time per measuring point. Max current 9 A, use a 4-5 A limit.

**Consistency / documentation issues**
- `ADC_VREF` differs between scripts (decoder / plotting / `plot_norm`: 5.06 V; `CASCADES26_expedition.py`: 4.85 V; lab 4.88 V). **Make them consistent.**
- Firmware reads the second telemetry channel on pin `A3` although everything calls it A2.
- The firmware v4 comments (A1 = speed, A2 = current) agree with the expedition `.edc`; the lab v3 results (A1 = current, A2 = speed) apply only to the lab `..._ao_shift.edc` wiring.
- `Running_Chipiron` was written for the v3 scripts (session asks phone ID + network number, 250 ms telemetry, `chipiron_plotting_v3.py` with step **and sine** fit cells, `--speed-p-gain` in the patch script). Current code: session v6 asks the full phone IP and the sampling period, plotting is v5 (step fit only), and the `chipiron_patch_edc_v3.py` in the project has no `--speed-p-gain` option. Re-add or update accordingly.
- `chipiron_oscillate_v4.py` / `chipiron_run_v4.py` docstrings still mention `PWM_full_v3.ino`; `chipiron_amp_sweep_v4.py` continuous sweep needs v4 firmware.
- `parameters.edc` "saved on 13 december 26" is later than today: fix the date.
- `chipiron_files_debug.py` docstring mentions `CASCADES26_expedition_v2.py`; the file name used is `CASCADES26_expedition.py`.
- Paths at the top of the analysis scripts (`C:/Users/Gulliver/...`, `D:/CASCADES26/...`) are hard-coded: edit them for each machine.
- Telemetry period: the number typed in the session must equal `TELEM_PERIOD_MS` of the flashed firmware (else lost-frame warnings and a wrong time scale).

**Script limitations**
- `/!\` **LoRa not tested** (`--link lora`, `LoraLink`, `USE_LORA 1`, bridge Arduino).
- No acknowledgement from the Arduino; no log of commands sent from browser/curl/`mobile_v3`; laptop and phone clocks not aligned automatically.
- Timing: commands over Wi-Fi have tens of ms jitter; ramps are staircases; stream sine is coarse; duty resolution 0.39 % (ramps can have several equal steps).
- ESCON ramp checks assume a start at rest and a simple limiter model; the `.edc` given to the scripts must be the one loaded in the controller (nothing verifies it); an offset stored in the `.edc` is ignored by the mapping unless enabled in the file.
- Only the session sends the rest command (0 rpm, then 30 s) before the start; standalone run/oscillate do not; the 30 s wait and the 30 s default rest before a program add up in the recording.
- No dry run apart from `--build-only` / `--check-only`; a run with the phone unreachable fails at the first command. Command buffer is 64 characters.
- Phone/software: `--song` needs the `play` command and `test.mp3`; file transfers time out after 60 s; HTTP routes depend on the Gobannos version; `tsync_v3.py` needs UDP ports 5000/5001 open on the laptop firewall; `usb <command>` in `mobile_v3.py` goes to all active phones.
- Lost bytes give `low8` cycles (old firmware) whose ADC value is only known modulo 256.
- `chipiron_io.py` not tested on real h5py / real data.

