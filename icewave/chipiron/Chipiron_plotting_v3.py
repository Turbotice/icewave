"""
chipiron_plotting_v3.py -- plot the Phonefleet logs of one or several acquisitions (date/time folders).


  - first cell : the root folder is hardcoded (ROOT_FOLDER below). The script lists the available dates, you
    pick one or several (e.g. "1 2"), then for each date it lists the available times and you pick one or
    several. It then lists the file types available (usb / accelerometer / gyroscope / magnetic_field / gps)
    and you pick the ones to plot: only those get questions below. If more than one acquisition is selected, you choose: all on the SAME figure, or ONE figure per
    acquisition.
  - then one cell per sensor, starting with Usb:   Usb -> accelerometer -> gyroscope -> magnetic_field -> gps
    each one prints the first 5 lines (of the first selected acquisition), asks which y / x to plot, the axis
    labels and linear / semilogx / semilogy / loglog. The same settings are applied to every selected
    acquisition. The log-scale lines are written explicitly in each cell (look for 'LOG SCALE').
  - last cell  : 'Save plots?' (none / all / numbers) -> <folder>/figures/

DEPENDENCIES
  chipiron_decode_v3.py   decoding of the Usb bytes (gap-based cycles, handles lost bytes) - REQUIRED
  escon_edc_v3.py         optional: converts the Usb volts to rpm / A with the .edc (asked in the Usb cell)

A1 / A2 = the two Arduino analog PINS. Wiring found on 2026-10-20 (ESCON Studio, motor held at 600 rpm, compared
with what the Arduino measured):   A1 <- ESCON AO1 = actual CURRENT,   A2 <- ESCON AO2 = actual SPEED.
(older versions of this script said A1 = speed, A2 = current: that was wrong). It is set once in WIRING below.
In the Usb cell you choose what to plot:  ADC counts (decoded, NOT converted) / volts / physical units (.edc),
or the 4 raw bytes by role (b1..b4). ADC counts are the safest way to check a measurement.
  numpy, pandas, matplotlib
INPUT   ROOT_FOLDER/YYYY-MM-DD/HH-MM-SS/ as written by chipiron_session_v3.py (or a mobile_v3.py
        download folder: files named P<phone>_D<original name>)

EXAMPLE (answers typed at the prompts)
  Dates:  [1] 2026-10-14  [2] 2026-10-15  -> Which ones? 2
  Times:  [1] 08-45-19  [2] 09-10-02  [3] 09-30-44  -> Which ones? 1 3
  Same figure (1) or one per acquisition (2)? 1
  Usb data: 1) decoded  -> Decoded: 3) A1 AND A2 -> .edc file for physical units (empty = none): parameters_Sep26.edc
            -> y = A2_speed_rpm (and A1_current_A), x = t_s, scale 1 (linear)  -> figures made, saved at the end ('all').
  Usb, ADC counts instead:  1) decoded -> units 1) ADC counts -> 3) A1 AND A2 -> y = A1_adc / A2_adc

What is known about the file contents (the csv files have NO header, so the column names below are an
interpretation, check them with whoever wrote Gobannos):
  - timestamp (column 0): microseconds. Only differences are used (t_s = seconds since the first sample).
  - accelerometer / gyroscope / magnetic_field: x, y, z in the Android convention (m/s^2, rad/s, microtesla).
  - Usb: ONE BYTE PER ROW (signed int8 as logged). The Arduino (PWM_full_v3.ino) sends 4 bytes every 250 ms:
    b1, b2 = A1 ADC (b1 = low byte, b2 = adc >> 4), b3, b4 = A2 ADC; chipiron_decode_v3 rebuilds A1_adc / A2_adc
    (+ A1_V / A2_V) and keeps the bytes as columns b1..b4.
    Row parity is NOT a byte role: if the recording starts in the middle of a cycle (or a byte is lost), "even rows"
    are b2/b4 instead of b1/b3. Use the "bytes by role" option, not the even/odd one, to look at raw bytes.
  - gps: 4 columns, meaning unknown (no data in the files seen so far).
"""
# %% Setup, helpers and folder selection
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))   # find chipiron_decode_v3 next to this file
from chipiron_decode_v3 import read_usb_rows, decode_usb

plt.ion()   # figures do not block the script

# --- ROOT FOLDER (hardcoded: edit here). Contains one sub-folder per date, each with one sub-folder per time. ---
ROOT_FOLDER = Path("C:/Users/matil/Desktop/PMMH/Chipiron/CASCADES26/tests_labo")

# key: (text to find in the file name, column names). GPS columns are generic names.
SENSORS = {
    "usb":            ("usb",            ["t_us", "usb_value"]),
    "accelerometer":  ("accelerometer",  ["t_us", "ax", "ay", "az"]),
    "gyroscope":      ("gyroscope",      ["t_us", "gx", "gy", "gz"]),
    "magnetic_field": ("magnetic_field", ["t_us", "mx", "my", "mz"]),
    "gps":            ("gps",            ["t_us", "gps_1", "gps_2", "gps_3"]),
}
ADC_VREF = 5.0              # volts at full scale (1023) of the Arduino ADC: measure the 5V pin of YOUR Arduino (USB ~4.9 V)
WIRING = {"A1": "AO1", "A2": "AO2"}   # Arduino pin -> ESCON analog output (A1 = AO1 current, A2 = AO2 speed). Verified 2026-10-20.
TIME_UNIT_US = 1e6          # microseconds per second (t_s = (t_us - t0) / TIME_UNIT_US)
FIG_SUBFOLDER = "figures"   # figures are saved in <data folder>/figures
FIG_FORMAT = "png"
FIG_DPI = 200
TS_RE = re.compile(r"(\d{4}-\d{2}-\d{2})T(\d{2})_(\d{2})_(\d{2})")
FIGS = []                   # registry of the figures made in this session
ACQS = []                   # selected acquisitions: dicts with label, folder, files, stamp
SAME_FIG = True             # several acquisitions: all on the same figure (True) or one figure each (False)


def ask(prompt, default=None):
  suffix = f" [{default}]" if default not in (None, "") else ""
  sys.stdout.write(f"{prompt}{suffix}: ")
  sys.stdout.flush()
  raw = sys.stdin.readline().strip()
  return raw if raw else (default if default is not None else "")

def find_files(folder):
    """Return {sensor_key: Path} for the csv files found in folder."""
    csvs = sorted(Path(folder).glob("*.csv"))
    print(f"{len(csvs)} csv file(s) in {folder}")
    found = {}
    for key, (pattern, _) in SENSORS.items():
        matches = [p for p in csvs if pattern in p.name.lower()]
        if len(matches) == 1:
            found[key] = matches[0]
            print(f"  [ok]      {key:15s} <- {matches[0].name}")
        elif len(matches) > 1:
            found[key] = matches[0]
            print(f"  [warning] {key:15s} {len(matches)} matches, using {matches[0].name}")
        else:
            print(f"  [missing] {key:15s} no file found")
    return found


def acquisition_stamp(files, folder):
    """'YYYY-MM-DD_HH-MM-SS' from the csv file names (else from the folder names)."""
    for p in files.values():
        m = TS_RE.search(p.name)
        if m:
            d, h, mi, s = m.groups()
            return f"{d}_{h}-{mi}-{s}"
    return "_".join(Path(folder).resolve().parts[-2:]).replace(" ", "")


def load(key, acq):
    """Load one sensor csv of one acquisition into a DataFrame (None if missing / unreadable)."""
    if key not in acq["files"]:
        print(f"\n[{key}] {acq['label']}: no file, skipping.")
        return None
    path = acq["files"][key]
    names = SENSORS[key][1]
    try:
        df = pd.read_csv(path, header=None, skipinitialspace=True)
    except Exception as exc:
        print(f"\n[{key}] could not read {path.name}: {exc}")
        return None
    if df.shape[1] > 1 and df.iloc[:, -1].isna().all():
        df = df.iloc[:, :-1]  # trailing comma -> empty last column
    if df.shape[1] == len(names):
        df.columns = names
    else:
        print(f"[{key}] {df.shape[1]} columns found, expected {len(names)}: using generic names")
        df.columns = [f"col{i}" for i in range(df.shape[1])]
    # The phone pre-allocates its buffers: unused rows are placeholders ("-1, NaN, ...").
    n_raw = len(df)
    valid = df.iloc[:, 0] > 0
    n_pad = int((~valid).sum())
    df = df[valid].reset_index(drop=True)
    print(f"\n[{key}] {acq['label']}: loaded {path.name}: {n_raw} rows read, {n_pad} placeholder row(s) "
          f"(timestamp <= 0) removed, {len(df)} kept")
    if len(df) == 0:
        print(f"[{key}] no valid data in this file, skipping.")
        return None
    t = df.iloc[:, 0]
    if t.nunique() > 1:
        df["t_s"] = (t - t.iloc[0]) / TIME_UNIT_US
        duration = df["t_s"].iloc[-1]
        print(f"[{key}] {df.shape[1] - 1} columns, duration {duration:.2f} s, "
              f"~{(len(df) - 1) / duration:.1f} Hz on average")
        if not t.is_monotonic_increasing:
            print(f"[{key}] warning: timestamps are not monotonic, x=t_s may look jumbled")
    else:
        print(f"[{key}] timestamp column is constant, no t_s added")
    n_nan = int(df.isna().sum().sum())
    if n_nan:
        print(f"[{key}] {n_nan} NaN value(s) left (they will not be drawn)")
    return df


def load_all(key):
    """Load one sensor for every selected acquisition -> list of (acq, DataFrame); empty if nothing usable."""
    out = []
    for acq in ACQS:
        df = load(key, acq)
        if df is not None:
            out.append((acq, df))
    if not out:
        print(f"\n[{key}] no usable data in the selected acquisition(s), skipping this section.")
    return out


def pick_column(df, prompt, default):
    """Ask for a column by index or name; re-ask until valid."""
    cols = list(df.columns)
    while True:
        raw = ask(prompt, default)
        if raw in cols:
            return raw
        if raw.lstrip("-").isdigit() and -len(cols) <= int(raw) < len(cols):
            return cols[int(raw)]
        print(f"  '{raw}' is not a valid column, use an index (0-{len(cols) - 1}) or a name.")


def ask_plot_settings(df, key, default_y=None):
    """Print head + columns, ask y, x, labels and scale. Returns a dict."""
    print(f"\n--- {key}: first 5 lines ---")
    print(df.head(5).to_string())
    print("\nColumns: " + ", ".join(f"{i}={c}" for i, c in enumerate(df.columns)))
    default_x = "t_s" if "t_s" in df.columns else df.columns[0]
    if default_y is None:
        default_y = df.columns[1] if df.shape[1] > 1 else df.columns[0]
    y = pick_column(df, "Column to plot (y), index or name", default_y)
    x = pick_column(df, "As a function of (x), index or name", default_x)
    pretty = lambda c: "time (s)" if c == "t_s" else c
    xlabel = ask("x axis label", pretty(x))
    ylabel = ask("y axis label", pretty(y))
    print("Scale:  1) linear   2) semilogx   3) semilogy   4) loglog")
    scale = ask("Choice", "1")
    logx = scale in ("2", "4")
    logy = scale in ("3", "4")
    for axis, col, is_log in (("x", x, logx), ("y", y, logy)):
        if is_log:
            n_bad = int((df[col] <= 0).sum())
            if n_bad:
                print(f"  warning: {n_bad} value(s) <= 0 in '{col}', hidden on a log {axis} axis")
    tag = ("log" if logx else "lin") + ("log" if logy else "lin")      # x scale then y scale
    print(f"Plotting y={y} vs x={x}, logx={logx}, logy={logy}")
    return dict(x=x, y=y, xlabel=xlabel, ylabel=ylabel, logx=logx, logy=logy, scale_tag=tag)


def has_cols(df, s, acq):
    """True if the chosen x / y columns exist in this acquisition's DataFrame (else say so and skip it)."""
    missing = [c for c in (s["x"], s["y"]) if c not in df.columns]
    if missing:
        print(f"  [{acq['label']}] column(s) {missing} not found, this acquisition is skipped on this figure")
        return False
    return True


def make_axes(group):
    """group = list of (acq, df). Returns [(fig, ax, sub_group)]: one figure for all (SAME_FIG) or one each."""
    if SAME_FIG or len(group) == 1:
        fig, ax = plt.subplots(figsize=(10, 5))
        return [(fig, ax, group)]
    out = []
    for item in group:
        fig, ax = plt.subplots(figsize=(10, 5))
        out.append((fig, ax, [item]))
    return out


def fig_title(group, key, variant=""):
    """File name in the title for one acquisition, a count for several."""
    if len(group) == 1:
        acq = group[0][0]
        base = acq["files"][key].name
        return f"{base} ({variant})" if variant else base
    return f"{key}{(' ' + variant) if variant else ''}: {len(group)} acquisitions"


def register_figure(fig, key, s, group, variant=""):
    """Remember the figure (for the final 'save plots?' prompt) and show it without blocking."""
    stamp = "+".join(a["stamp"] for a, _ in group)
    FIGS.append(dict(fig=fig, key=key, variant=variant, x=s["x"], y=s["y"], tag=s["scale_tag"], stamp=stamp))
    print(f"  -> figure #{len(FIGS)}: {describe(FIGS[-1])}")
    plt.show(block=False)
    fig.canvas.draw_idle()
    plt.pause(0.1)


def describe(f):
    return f"{f['stamp']} {f['key']}{(' ' + f['variant']) if f['variant'] else ''}: {f['y']} vs {f['x']} [{f['tag']}]"


def fig_filename(f):
    """e.g. 2026-10-15_08-45-19_usb_A1_A1_adc_vs_t_s_linlin.png (several acquisitions: stamps joined by '+')"""
    parts = [f["stamp"], f["key"], f["variant"], f["y"], "vs", f["x"], f["tag"]]
    name = "_".join(p for p in parts if p)
    return re.sub(r"[^A-Za-z0-9._+-]+", "-", name) + f".{FIG_FORMAT}"


def parse_selection(raw, n):
    """'all', 'none', '1,3', '1 3', '2-4' -> sorted list of numbers (1-based)."""
    raw = raw.strip().lower()
    if raw in ("", "none", "n", "no", "0"):
        return []
    if raw in ("all", "a"):
        return list(range(1, n + 1))
    picked = set()
    for part in re.split(r"[,\s]+", raw):
        m = re.fullmatch(r"(\d+)-(\d+)", part)
        if m:
            picked.update(range(int(m.group(1)), int(m.group(2)) + 1))
        elif part.isdigit():
            picked.add(int(part))
        elif part:
            print(f"  ignored '{part}'")
    bad = [i for i in picked if not 1 <= i <= n]
    if bad:
        print(f"  ignored out-of-range number(s): {sorted(bad)}")
    return sorted(i for i in picked if 1 <= i <= n)


def choose(title, names):
    """List names as [1] .. [n], ask which ones (e.g. '1 2 3', '2-4', 'all'); re-ask until valid. Returns names."""
    print(f"\n{title}")
    for i, n in enumerate(names, start=1):
        print(f"  [{i}] {n}")
    while True:
        picked = parse_selection(ask("Which ones? (numbers e.g. 1 2 3 / 2-4 / all)"), len(names))
        if picked:
            return [names[i - 1] for i in picked]
        print("  select at least one valid number.")


if not ROOT_FOLDER.is_dir():
    raise SystemExit(f"Root folder not found: {ROOT_FOLDER.resolve()} (edit ROOT_FOLDER at the top of the script)")
date_dirs = sorted(p for p in ROOT_FOLDER.iterdir() if p.is_dir())
if not date_dirs:
    raise SystemExit(f"No date sub-folder in {ROOT_FOLDER.resolve()}")
dates = choose(f"Dates available in {ROOT_FOLDER}:", [p.name for p in date_dirs])
for d in dates:
    time_dirs = sorted(p for p in (ROOT_FOLDER / d).iterdir() if p.is_dir())
    if not time_dirs:
        print(f"\nNo time sub-folder in {ROOT_FOLDER / d}, skipping this date.")
        continue
    times = choose(f"Times available on {d}:", [p.name for p in time_dirs])
    for t in times:
        folder = ROOT_FOLDER / d / t
        print(f"\n=== {d} {t} ===")
        files = find_files(folder)
        ACQS.append(dict(label=f"{d} {t}", folder=folder, files=files, stamp=acquisition_stamp(files, folder)))
if not ACQS:
    raise SystemExit("No acquisition selected.")
print(f"\n{len(ACQS)} acquisition(s) selected: " + ", ".join(a["label"] for a in ACQS))
available = [k for k in SENSORS if any(k in a["files"] for a in ACQS)]   # sensors with a file in at least one acquisition
if not available:
    raise SystemExit("No sensor file found in the selected acquisition(s).")
SELECTED = choose("File types available in the selected acquisition(s):", available)
print("File types to plot: " + ", ".join(SELECTED))
if len(ACQS) > 1:
    print("\nSeveral acquisitions selected:  1) all on the SAME figure   2) ONE figure per acquisition")
    SAME_FIG = ask("Choice", "1") != "2"

# %% Usb
df_usb_raw = [(a, read_usb_rows(a["files"]["usb"])) for a in ACQS if "usb" in SELECTED and "usb" in a["files"]]
df_usb_raw = [(a, d) for a, d in df_usb_raw if d is not None and len(d)]
if df_usb_raw:
    # The Usb file holds ONE BYTE PER ROW: 4 rows = A1 low/high byte + A2 low/high byte.
    print("\nUsb data:")
    print("  1) decoded ADC channels A1 / A2 (chipiron_decode_v3, as PWM_full_v3.ino sends them)")
    print("  2) raw rows as logged (one byte per row; even/odd = ROW parity, may not match the byte roles)")
    print("  3) the 4 bytes by ROLE: b1 / b2 (A1 low byte, A1 >> 4) and b3 / b4 (A2 low byte, A2 >> 4)")
    kind = ask("Choice", "1")
    jobs = []   # (default y, variant name, [(acq, DataFrame), ...])
    if kind == "2":
        print("Raw rows: 1) all   2) even rows only   3) odd rows only   4) even AND odd (two figures)")
        print("  (row parity depends on where the recording started: for a fixed role use option 3 of the Usb menu)")
        sel = ask("Choice", "1")
        even = [(a, d.iloc[0::2].reset_index(drop=True)) for a, d in df_usb_raw]
        odd = [(a, d.iloc[1::2].reset_index(drop=True)) for a, d in df_usb_raw]
        if sel == "2":
            jobs = [("usb_value", "raw-even", even)]
        elif sel == "3":
            jobs = [("usb_value", "raw-odd", odd)]
        elif sel == "4":
            jobs = [("usb_value", "raw-even", even), ("usb_value", "raw-odd", odd)]
        else:
            jobs = [("usb_value", "raw-all", df_usb_raw)]
        for (a, d), (_, e), (_, o) in zip(df_usb_raw, even, odd):
            print(f"{a['label']}: {len(d)} raw rows ({len(e)} even, {len(o)} odd)")
    elif kind == "3":
        # bytes by role, taken from the decoded cycles (not from the row parity)
        dec_b = [(a, decode_usb(d, vref=ADC_VREF)) for a, d in df_usb_raw]
        print("Bytes: 1) b1 (A1 low byte)  2) b2 (A1 >> 4)  3) b3 (A2 low byte)  4) b4 (A2 >> 4)  5) all four (4 figures)")
        sel = ask("Choice", "5")
        names = {"1": ["b1"], "2": ["b2"], "3": ["b3"], "4": ["b4"]}.get(sel, ["b1", "b2", "b3", "b4"])
        jobs = [(nm, nm, dec_b) for nm in names]
    else:
        # What to plot: ADC counts (decoded, not converted), volts at the pin, or physical units from the .edc.
        print("Units:  1) ADC counts (decoded, NOT converted)   2) volts at the Arduino pin   3) physical units (.edc)")
        units = ask("Choice", "1")
        edc_path = ""
        if units == "3":   # volts -> rpm / A with the ESCON file, using WIRING (A1 <- AO1 current, A2 <- AO2 speed)
            edc_path = ask(".edc file loaded in ESCON Studio when this was recorded").strip('"').strip("'")
        dec = []   # (acq, decoded DataFrame, name of A1 column, name of A2 column)
        for a, d in df_usb_raw:
            df_dec = decode_usb(d, vref=ADC_VREF)
            d1, d2 = ("A1_V", "A2_V") if units == "2" else ("A1_adc", "A2_adc")
            if edc_path:
                try:
                    from escon_edc_v3 import load_edc
                    load_edc(edc_path).apply_to_frame(df_dec, wiring=WIRING)
                    d1 = next((c for c in df_dec.columns if c.startswith("A1_") and c.endswith(("_rpm", "_A"))), d1)
                    d2 = next((c for c in df_dec.columns if c.startswith("A2_") and c.endswith(("_rpm", "_A"))), d2)
                except Exception as exc:
                    print(f"[usb] {a['label']}: could not use the .edc ({exc}): ADC counts only")
            if df_dec["A1_adc"].isna().all() and (df_dec["kind"] == "low8").any():
                print(f"[usb] {a['label']}: only low bytes were logged: plotting A1_lo / A2_lo (adc mod 256, aliased)")
                d1, d2 = "A1_lo", "A2_lo"
            dec.append((a, df_dec, d1, d2))
        print("Decoded: 1) A1 only   2) A2 only   3) A1 AND A2 (two figures)")
        sel = ask("Choice", "3")
        # the default y column is taken from the first acquisition
        g = [(a, df) for a, df, _, _ in dec]
        if sel == "1":
            jobs = [(dec[0][2], "A1", g)]
        elif sel == "2":
            jobs = [(dec[0][3], "A2", g)]
        else:
            jobs = [(dec[0][2], "A1", g), (dec[0][3], "A2", g)]

    for default_y, variant, group in jobs:
        s = ask_plot_settings(group[0][1], f"usb {variant}", default_y)

        for fig, ax, sub in make_axes(group):
            for acq, df_plot in sub:
                if not has_cols(df_plot, s, acq):
                    continue
                ax.plot(df_plot[s["x"]], df_plot[s["y"]], lw=0.8, label=acq["label"])
            ax.set_xlabel(s["xlabel"])
            ax.set_ylabel(s["ylabel"])
            # --- LOG SCALE (hardcoded: edit here if needed) ---
            if s["logx"]:
                ax.set_xscale("log")
            if s["logy"]:
                ax.set_yscale("log")
            ax.set_title(fig_title(sub, "usb", variant), fontsize=8)
            if len(sub) > 1:
                ax.legend(fontsize=8)
            ax.grid(True, which="both", alpha=0.3)
            fig.tight_layout()
            register_figure(fig, "usb", s, sub, variant)
elif "usb" in SELECTED:
    print("\n[usb] no usable file in the selected acquisition(s), skipping this section.")

# %% Accelerometer
dfs_acc = load_all("accelerometer") if "accelerometer" in SELECTED else []
if dfs_acc:
    s = ask_plot_settings(dfs_acc[0][1], "accelerometer")

    for fig, ax, sub in make_axes(dfs_acc):
        for acq, df_acc in sub:
            if not has_cols(df_acc, s, acq):
                continue
            ax.plot(df_acc[s["x"]], df_acc[s["y"]], lw=0.8, label=acq["label"])
        ax.set_xlabel(s["xlabel"])
        ax.set_ylabel(s["ylabel"])
        # --- LOG SCALE (hardcoded: edit here if needed) ---
        if s["logx"]:
            ax.set_xscale("log")
        if s["logy"]:
            ax.set_yscale("log")
        ax.set_title(fig_title(sub, "accelerometer"), fontsize=8)
        if len(sub) > 1:
            ax.legend(fontsize=8)
        ax.grid(True, which="both", alpha=0.3)
        fig.tight_layout()
        register_figure(fig, "accelerometer", s, sub)

# %% Gyroscope
dfs_gyro = load_all("gyroscope") if "gyroscope" in SELECTED else []
if dfs_gyro:
    s = ask_plot_settings(dfs_gyro[0][1], "gyroscope")

    for fig, ax, sub in make_axes(dfs_gyro):
        for acq, df_gyro in sub:
            if not has_cols(df_gyro, s, acq):
                continue
            ax.plot(df_gyro[s["x"]], df_gyro[s["y"]], lw=0.8, label=acq["label"])
        ax.set_xlabel(s["xlabel"])
        ax.set_ylabel(s["ylabel"])
        # --- LOG SCALE (hardcoded: edit here if needed) ---
        if s["logx"]:
            ax.set_xscale("log")
        if s["logy"]:
            ax.set_yscale("log")
        ax.set_title(fig_title(sub, "gyroscope"), fontsize=8)
        if len(sub) > 1:
            ax.legend(fontsize=8)
        ax.grid(True, which="both", alpha=0.3)
        fig.tight_layout()
        register_figure(fig, "gyroscope", s, sub)

# %% Magnetic field
dfs_mag = load_all("magnetic_field") if "magnetic_field" in SELECTED else []
if dfs_mag:
    s = ask_plot_settings(dfs_mag[0][1], "magnetic_field")

    for fig, ax, sub in make_axes(dfs_mag):
        for acq, df_mag in sub:
            if not has_cols(df_mag, s, acq):
                continue
            ax.plot(df_mag[s["x"]], df_mag[s["y"]], lw=0.8, label=acq["label"])
        ax.set_xlabel(s["xlabel"])
        ax.set_ylabel(s["ylabel"])
        # --- LOG SCALE (hardcoded: edit here if needed) ---
        if s["logx"]:
            ax.set_xscale("log")
        if s["logy"]:
            ax.set_yscale("log")
        ax.set_title(fig_title(sub, "magnetic_field"), fontsize=8)
        if len(sub) > 1:
            ax.legend(fontsize=8)
        ax.grid(True, which="both", alpha=0.3)
        fig.tight_layout()
        register_figure(fig, "magnetic_field", s, sub)

# %% GPS
dfs_gps = load_all("gps") if "gps" in SELECTED else []
# keep only the acquisitions with at least one non-NaN value besides the timestamp
dfs_gps_ok = [(a, d) for a, d in dfs_gps if not d.drop(columns=d.columns[0]).isna().all().all()]
if dfs_gps and not dfs_gps_ok:
    print("[gps] all values are NaN (no fix during this acquisition): nothing to plot.")
elif dfs_gps_ok:
    s = ask_plot_settings(dfs_gps_ok[0][1], "gps")

    for fig, ax, sub in make_axes(dfs_gps_ok):
        for acq, df_gps in sub:
            if not has_cols(df_gps, s, acq):
                continue
            ax.plot(df_gps[s["x"]], df_gps[s["y"]], ".-", lw=0.8, label=acq["label"])
        ax.set_xlabel(s["xlabel"])
        ax.set_ylabel(s["ylabel"])
        # --- LOG SCALE (hardcoded: edit here if needed) ---
        if s["logx"]:
            ax.set_xscale("log")
        if s["logy"]:
            ax.set_yscale("log")
        ax.set_title(fig_title(sub, "gps"), fontsize=8)
        if len(sub) > 1:
            ax.legend(fontsize=8)
        ax.grid(True, which="both", alpha=0.3)
        fig.tight_layout()
        register_figure(fig, "gps", s, sub)

# %% Save plots?
if not FIGS:
    print("No figure was made, nothing to save.")
else:
    print(f"\n{len(FIGS)} figure(s) made:")
    for i, f in enumerate(FIGS, start=1):
        print(f"  {i}) {describe(f)}")
    chosen = parse_selection(ask("Save plots? none / all / numbers (e.g. 1,3 or 2-4)", "none"), len(FIGS))
    if chosen:
        # one acquisition -> its own folder/figures ; several -> ROOT_FOLDER/figures
        default_out = (ACQS[0]["folder"] if len(ACQS) == 1 else ROOT_FOLDER) / FIG_SUBFOLDER
        out_dir = Path(ask("Save in folder", str(default_out))).expanduser()
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"Saving {len(chosen)} figure(s) in {out_dir.resolve()}")
        for i in chosen:
            f = FIGS[i - 1]
            path = out_dir / fig_filename(f)
            n = 2
            while path.exists():                      # never overwrite an existing figure
                path = path.with_name(f"{path.stem}_{n}{path.suffix}")
                n += 1
            try:
                f["fig"].savefig(path, dpi=FIG_DPI, bbox_inches="tight")
                print(f"  #{i} saved {path.name}")
            except Exception as exc:
                print(f"  #{i} could not be saved: {exc}")
    else:
        print("Nothing saved.")

# Keep the figure windows open when this file is run as a plain script (not in Spyder / VS Code cells).
try:
    get_ipython()   # noqa: F821  (defined only inside IPython / Spyder / VS Code interactive)
except NameError:
    if FIGS:
        print("Figures stay open: close the windows to end the script.")
        plt.ioff()
        plt.show()