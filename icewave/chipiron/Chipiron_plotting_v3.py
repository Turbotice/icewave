"""
chipiron_plotting_v3.py -- plot the Phonefleet logs of ONE acquisition (one date/time folder).


  - first cell : asks for the folder, finds and lists the files
  - then one cell per sensor, starting with Usb:   Usb -> accelerometer -> gyroscope -> magnetic_field -> gps
    each one prints the first 5 lines, asks which y / x to plot, the axis labels and linear / semilogx /
    semilogy / loglog. The log-scale lines are written explicitly in each cell (look for 'LOG SCALE').
  - last cell  : 'Save plots?' (none / all / numbers) -> <folder>/figures/

DEPENDENCIES
  chipiron_decode_v3.py   decoding of the Usb bytes (gap-based cycles, handles lost bytes) - REQUIRED
  escon_edc_v3.py         optional: converts the Usb volts to rpm / A with the .edc (asked in the Usb cell)
  numpy, pandas, matplotlib
INPUT   the folder written by chipiron_session_v3.py:  Phonefleet_data/YYYY-MM-DD/HH-MM-SS/ (or a mobile_v3.py
        download folder: files named P<phone>_D<original name>)

EXAMPLE (answers typed at the prompts)
  Folder containing the log files of one date/time: Phonefleet_data/2026-10-15/08-45-19
  Usb data: 1) decoded  -> Decoded: 3) A1 AND A2 -> .edc file for physical units (empty = none): parameters_Sep26.edc
            -> y = A1_speed_rpm, x = t_s, scale 1 (linear)  -> figures made, saved at the end ('all').

What is known about the file contents (the csv files have NO header, so the column names below are an
interpretation, check them with whoever wrote Gobannos):
  - timestamp (column 0): microseconds. Only differences are used (t_s = seconds since the first sample).
  - accelerometer / gyroscope / magnetic_field: x, y, z in the Android convention (m/s^2, rad/s, microtesla).
  - Usb: ONE BYTE PER ROW (signed int8 as logged). The Arduino (PWM_full_v3.ino) sends 4 bytes every 250 ms:
    b1, b2 = A1 ADC, b3, b4 = A2 ADC; chipiron_decode_v3 rebuilds A1_adc / A2_adc (+ A1_V / A2_V).
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

# key: (text to find in the file name, column names). GPS columns are generic names.
SENSORS = {
    "usb":            ("usb",            ["t_us", "usb_value"]),
    "accelerometer":  ("accelerometer",  ["t_us", "ax", "ay", "az"]),
    "gyroscope":      ("gyroscope",      ["t_us", "gx", "gy", "gz"]),
    "magnetic_field": ("magnetic_field", ["t_us", "mx", "my", "mz"]),
    "gps":            ("gps",            ["t_us", "gps_1", "gps_2", "gps_3"]),
}
TIME_UNIT_US = 1e6          # microseconds per second (t_s = (t_us - t0) / TIME_UNIT_US)
FIG_SUBFOLDER = "figures"   # figures are saved in <data folder>/figures
FIG_FORMAT = "png"
FIG_DPI = 200
TS_RE = re.compile(r"(\d{4}-\d{2}-\d{2})T(\d{2})_(\d{2})_(\d{2})")
FIGS = []                   # registry of the figures made in this session


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


def load(key):
    """Load one sensor csv into a DataFrame (None if missing / unreadable)."""
    if key not in FILES:
        print(f"\n[{key}] no file, skipping this section.")
        return None
    path = FILES[key]
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
    print(f"\n[{key}] loaded {path.name}: {n_raw} rows read, {n_pad} placeholder row(s) "
          f"(timestamp <= 0) removed, {len(df)} kept")
    if len(df) == 0:
        print(f"[{key}] no valid data in this file, skipping this section.")
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


def register_figure(fig, key, s, variant=""):
    """Remember the figure (for the final 'save plots?' prompt) and show it without blocking."""
    FIGS.append(dict(fig=fig, key=key, variant=variant, x=s["x"], y=s["y"], tag=s["scale_tag"]))
    print(f"  -> figure #{len(FIGS)}: {describe(FIGS[-1])}")
    plt.show(block=False)
    fig.canvas.draw_idle()
    plt.pause(0.1)


def describe(f):
    return f"{f['key']}{(' ' + f['variant']) if f['variant'] else ''}: {f['y']} vs {f['x']} [{f['tag']}]"


def fig_filename(f):
    """e.g. 2026-10-15_08-45-19_usb_A1_A1_adc_vs_t_s_linlin.png"""
    parts = [STAMP, f["key"], f["variant"], f["y"], "vs", f["x"], f["tag"]]
    name = "_".join(p for p in parts if p)
    return re.sub(r"[^A-Za-z0-9._-]+", "-", name) + f".{FIG_FORMAT}"


def parse_selection(raw, n):
    """'all', 'none', '1,3', '2-4' -> sorted list of figure numbers (1-based)."""
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


folder_in = ask("Folder containing the log files of one date/time").strip('"').strip("'")
FOLDER = Path(folder_in).expanduser()
if not FOLDER.is_dir():
    raise SystemExit(f"Not a folder: {FOLDER}")
FILES = find_files(FOLDER)
STAMP = acquisition_stamp(FILES, FOLDER)
print(f"Acquisition stamp used in figure file names: {STAMP}")

# %% Usb
df_usb_raw = read_usb_rows(FILES["usb"]) if "usb" in FILES else None
if df_usb_raw is not None and len(df_usb_raw):
    # The Usb file holds ONE BYTE PER ROW: 4 rows = A1 low/high byte + A2 low/high byte.
    print("\nUsb data:")
    print("  1) decoded ADC channels A1 / A2 (chipiron_decode_v3, as PWM_full_v3.ino sends them)")
    print("  2) raw rows as logged (one byte per row)")
    kind = ask("Choice", "1")
    jobs = []   # (DataFrame, default y, variant name)
    if kind == "2":
        print("Raw rows: 1) all   2) even rows only   3) odd rows only   4) even AND odd (two figures)")
        sel = ask("Choice", "1")
        even = df_usb_raw.iloc[0::2].reset_index(drop=True)
        odd = df_usb_raw.iloc[1::2].reset_index(drop=True)
        if sel == "2":
            jobs = [(even, "usb_value", "raw-even")]
        elif sel == "3":
            jobs = [(odd, "usb_value", "raw-odd")]
        elif sel == "4":
            jobs = [(even, "usb_value", "raw-even"), (odd, "usb_value", "raw-odd")]
        else:
            jobs = [(df_usb_raw, "usb_value", "raw-all")]
        print(f"{len(df_usb_raw)} raw rows ({len(even)} even, {len(odd)} odd)")
    else:
        df_dec = decode_usb(df_usb_raw)
        # optional conversion volts -> rpm / A with the ESCON file (A1 = speed = AO2, A2 = current = AO1: CHECK wiring)
        edc_path = ask(".edc file for physical units (rpm, A); empty = ADC counts only").strip('"').strip("'")
        d1, d2 = "A1_adc", "A2_adc"
        if edc_path:
            try:
                from escon_edc_v3 import load_edc
                load_edc(edc_path).apply_to_frame(df_dec)
                d1 = next((c for c in df_dec.columns if c.startswith("A1_") and c.endswith(("_rpm", "_A"))), d1)
                d2 = next((c for c in df_dec.columns if c.startswith("A2_") and c.endswith(("_rpm", "_A"))), d2)
            except Exception as exc:
                print(f"[usb] could not use the .edc ({exc}): ADC counts only")
        if df_dec["A1_adc"].isna().all() and (df_dec["kind"] == "low8").any():
            print("[usb] only low bytes were logged: plotting A1_lo / A2_lo (adc mod 256, aliased)")
            d1, d2 = "A1_lo", "A2_lo"
        print("Decoded: 1) A1 only   2) A2 only   3) A1 AND A2 (two figures)")
        sel = ask("Choice", "3")
        if sel == "1":
            jobs = [(df_dec, d1, "A1")]
        elif sel == "2":
            jobs = [(df_dec, d2, "A2")]
        else:
            jobs = [(df_dec, d1, "A1"), (df_dec, d2, "A2")]

    for df_plot, default_y, variant in jobs:
        s = ask_plot_settings(df_plot, f"usb {variant}", default_y)

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(df_plot[s["x"]], df_plot[s["y"]], lw=0.8)
        ax.set_xlabel(s["xlabel"])
        ax.set_ylabel(s["ylabel"])
        # --- LOG SCALE (hardcoded: edit here if needed) ---
        if s["logx"]:
            ax.set_xscale("log")
        if s["logy"]:
            ax.set_yscale("log")
        ax.set_title(f"{FILES['usb'].name} ({variant})", fontsize=8)
        ax.grid(True, which="both", alpha=0.3)
        fig.tight_layout()
        register_figure(fig, "usb", s, variant)

# %% Accelerometer
df_acc = load("accelerometer")
if df_acc is not None:
    s = ask_plot_settings(df_acc, "accelerometer")

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(df_acc[s["x"]], df_acc[s["y"]], lw=0.8)
    ax.set_xlabel(s["xlabel"])
    ax.set_ylabel(s["ylabel"])
    # --- LOG SCALE (hardcoded: edit here if needed) ---
    if s["logx"]:
        ax.set_xscale("log")
    if s["logy"]:
        ax.set_yscale("log")
    ax.set_title(FILES["accelerometer"].name, fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    register_figure(fig, "accelerometer", s)

# %% Gyroscope
df_gyro = load("gyroscope")
if df_gyro is not None:
    s = ask_plot_settings(df_gyro, "gyroscope")

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(df_gyro[s["x"]], df_gyro[s["y"]], lw=0.8)
    ax.set_xlabel(s["xlabel"])
    ax.set_ylabel(s["ylabel"])
    # --- LOG SCALE (hardcoded: edit here if needed) ---
    if s["logx"]:
        ax.set_xscale("log")
    if s["logy"]:
        ax.set_yscale("log")
    ax.set_title(FILES["gyroscope"].name, fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    register_figure(fig, "gyroscope", s)

# %% Magnetic field
df_mag = load("magnetic_field")
if df_mag is not None:
    s = ask_plot_settings(df_mag, "magnetic_field")

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(df_mag[s["x"]], df_mag[s["y"]], lw=0.8)
    ax.set_xlabel(s["xlabel"])
    ax.set_ylabel(s["ylabel"])
    # --- LOG SCALE (hardcoded: edit here if needed) ---
    if s["logx"]:
        ax.set_xscale("log")
    if s["logy"]:
        ax.set_yscale("log")
    ax.set_title(FILES["magnetic_field"].name, fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    register_figure(fig, "magnetic_field", s)

# %% GPS
df_gps = load("gps")
if df_gps is not None:
    if df_gps.drop(columns=df_gps.columns[0]).isna().all().all():
        print("[gps] all values are NaN (no fix during this acquisition): nothing to plot.")
    else:
        s = ask_plot_settings(df_gps, "gps")

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(df_gps[s["x"]], df_gps[s["y"]], ".-", lw=0.8)
        ax.set_xlabel(s["xlabel"])
        ax.set_ylabel(s["ylabel"])
        # --- LOG SCALE (hardcoded: edit here if needed) ---
        if s["logx"]:
            ax.set_xscale("log")
        if s["logy"]:
            ax.set_yscale("log")
        ax.set_title(FILES["gps"].name, fontsize=8)
        ax.grid(True, which="both", alpha=0.3)
        fig.tight_layout()
        register_figure(fig, "gps", s)

# %% Save plots?
if not FIGS:
    print("No figure was made, nothing to save.")
else:
    print(f"\n{len(FIGS)} figure(s) made:")
    for i, f in enumerate(FIGS, start=1):
        print(f"  {i}) {describe(f)}")
    chosen = parse_selection(ask("Save plots? none / all / numbers (e.g. 1,3 or 2-4)", "none"), len(FIGS))
    if chosen:
        out_dir = Path(ask("Save in folder", str(FOLDER / FIG_SUBFOLDER))).expanduser()
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