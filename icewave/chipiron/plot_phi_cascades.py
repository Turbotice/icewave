"""
Plot the solid fraction phi of grease ice / frazil / slush computed by solid_fraction_frazil_MausDeLaRosa_v1.py

WHAT IT DOES
    1. Looks into the results folder (frazil_carac) and loads every  grease_ice_results_*_<direct|conductivity>.csv
       (one file per sample and per salinity mode; the long-format files written by the solid-fraction script).
    2. Asks how ALL plots must be drawn (the choice is applied to every figure):
           x axis   : day | station | sample (chronological) | datetime
           error bars on/off, y-limits   (saving is set by save_plot / plot_output_folder at the top)
    3. Figure 1: phi vs chosen x, ONE SUBPLOT PER METHOD (1-4). In each subplot the conductivity-based
       (PSS-78, circles) and salinity-based (direct input in psu, squares) results are compared.
    4. Figure 2: phi vs chosen x, ALL METHODS on the same axes (one marker per method), in TWO SUBPLOTS:
       top = conductivity-based, bottom = salinity-based.
    5. Every point is labelled with its station number (label_stations), except with x = station where the station
       is written once on the x axis and ALL series of a station (methods / conductivity vs salinity) share exactly
       the same x, so they only differ vertically. The different characterization iterations of a station also share
       that x (colour tells them apart); set spread_iterations = True to offset them slightly instead.
       A missing station in a results file is filled in from the raw CSV note (matched by date + time).
    6. Colour = CHARACTERIZATION ITERATION, read from the 'note' column of the RAW CSV (raw_csv_path), e.g.
       "... station 301 -2e caractérisation ..." -> iteration 2. Samples whose note has no iteration number
       (characterized once) get the colour of iteration 1. One discrete colour bar for every figure.
       Raw rows are matched to the results files by (date, time, station); the station is read from the note
       ("station 301", "(station 400)"). If the raw CSV cannot be found, iterations are inferred instead from the
       results files (same station, different date/time -> 1, 2, 3 in chronological order).

NOTE ON DATES / TIMES
    * A date with year < 1900 (e.g. '0001-11-11', produced when the raw CSV had '11/11' without year) is read as
      the same day/month in `default_year`.
    * The 'Time' column of the results files comes from the 'time' column of the raw CSV. A sample without time is
      placed at 12:00 of its day. Old results files without time for a day+station that now has timed results
      are ignored (delete them from the folder).
    * Raw dates may be 'dd/mm/yyyy' or 'dd-mon' (e.g. '11-nov', year = default_year).
"""

import functools
import re
import unicodedata
from datetime import datetime
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D

# --- User configuration ---------------------------------------------------
results_folder = r"D:\CASCADES26\chipiron\frazil_carac"   # None = ask when the script runs
raw_csv_path = r"D:\CASCADES26\chipiron\CASCADES26_frazil_characterization_excel_FR.csv"   # raw CSV with the 'note' column
default_year = 2026                                       # used when a date has an impossible year (e.g. 0001)
cmap_name = "turbo"                                       # sampled at one colour per iteration; try 'viridis', 'coolwarm'
save_plot = True                                          # True = save the figures (PNG) in plot_output_folder
plot_output_folder = r"D:\CASCADES26\chipiron\frazil_carac\plots"   # hardcoded folder for the saved figures
interactive = True                                        # False = use the defaults below without asking
label_stations = True                                     # write the station number next to every point (not needed with x = station)
spread_iterations = False                                 # x = station only: True = repeated characterizations of a
                                                          # station side by side (small offset); False = exactly the same x
PHI = r"$\phi$"                                           # symbol used in titles / axis labels
defaults = dict(xaxis="day",           # day | station | sample | datetime
                errorbars=True,
                ylim=(-0.05, 1.05))    # None = automatic

METHODS = (1, 2, 3, 4)
METHOD_NAMES = {1: "Method 1 (S_g, S_b)", 2: "Method 2 (rho_g, S_b)",
                3: "Method 3 (S_g, rho_g)", 4: "Method 4 (drainage)"}
METHOD_MARKER = {1: "o", 2: "s", 3: "^", 4: "D"}
MODES = ("conductivity", "direct")
MODE_LABEL = {"conductivity": "conductivity-based (PSS-78)", "direct": "salinity-based (direct, psu)"}
MODE_MARKER = {"conductivity": "o", "direct": "s"}


# %% Loading -----------------------------------------------------------------
def parse_date(s):
    """'YYYY-MM-DD' -> datetime.date-like (year < 1900 replaced by default_year); None if unreadable."""
    m = re.match(r"\s*(\d{1,4})-(\d{1,2})-(\d{1,2})", str(s))
    if not m:
        return None
    y, mo, d = (int(g) for g in m.groups())
    if y < 1900:
        y = default_year
    try:
        return datetime(y, mo, d)
    except ValueError:
        return None


def parse_time(s):
    """'HH:MM' (also 'HHhMM', 'HH:MM:SS') -> hours (float); NaN if empty."""
    m = re.match(r"\s*(\d{1,2})\s*[:hH.]\s*(\d{2})(?:\s*[:.]\s*(\d{2}))?", str(s))
    return int(m.group(1)) + int(m.group(2)) / 60 + (int(m.group(3)) / 3600 if m.group(3) else 0) if m else np.nan

_MONTHS = {"jan": 1, "feb": 2, "fev": 2, "mar": 3, "apr": 4, "avr": 4, "may": 5, "mai": 5, "jun": 6, "jui": 6,
           "jul": 7, "aug": 8, "aou": 8, "sep": 9, "oct": 10, "nov": 11, "dec": 12}


def parse_raw_date(s):
    """Raw-CSV date: 'dd/mm/yyyy', 'dd-mm-yyyy', 'dd-mon' ('11-nov'), 'yyyy-mm-dd' -> datetime or None."""
    s = unicodedata.normalize("NFKD", str(s)).encode("ascii", "ignore").decode().strip().lower()
    m = re.match(r"(\d{1,2})[/.-](\d{1,2})[/.-](\d{2,4})", s)
    if m:
        d, mo, y = (int(g) for g in m.groups())
        y += 2000 if y < 100 else 0
    elif re.match(r"\d{4}-\d{1,2}-\d{1,2}", s):
        return parse_date(s)
    else:
        m = re.match(r"(\d{1,2})[-/ ]([a-z]{3})", s)
        if not m or m.group(2) not in _MONTHS:
            return None
        d, mo, y = int(m.group(1)), _MONTHS[m.group(2)], default_year
    try:
        return datetime(y, mo, d)
    except ValueError:
        return None


@functools.lru_cache(maxsize=None)
def read_note_iterations(path):
    """Raw CSV -> ({(day, minute_of_day, station): iteration}, {(day, station): iteration or None if ambiguous},
                   {(day, minute_of_day): station or None if ambiguous}).
    Iteration from the note ('-2e caractérisation' -> 2); a note without a number -> 1."""
    raw = pd.read_csv(path, sep=None, engine="python", dtype=str, encoding="utf-8-sig")
    raw.columns = [c.strip().lower() for c in raw.columns]
    if not {"date", "time", "note"} <= set(raw.columns):
        raise ValueError("raw CSV needs 'date', 'time' and 'note' columns")
    by_dts, by_ds, by_dm = {}, {}, {}
    for r in raw.itertuples():
        note = str(r.note) if isinstance(r.note, str) else ""
        ms = re.search(r"station\s*(\d+)", note, re.I)
        day = parse_raw_date(r.date)
        if not ms or day is None:
            continue
        mi = re.search(r"-\s*(\d+)\s*(?:e|er|ere|re|st|nd|rd|th)\b\s*caract", unicodedata.normalize("NFKD", note)
                       .encode("ascii", "ignore").decode(), re.I)
        it = int(mi.group(1)) if mi else 1
        h = parse_time(r.time)
        st = ms.group(1)
        if not np.isnan(h):
            by_dts[(day, int(round(h * 60)), st)] = it
            k = (day, int(round(h * 60)))
            by_dm[k] = st if by_dm.get(k, st) == st else None
        by_ds[(day, st)] = it if by_ds.get((day, st), it) == it else None
    return by_dts, by_ds, by_dm


def fill_missing_stations(df):
    """Results rows with an empty Station: take the station from the raw-CSV note (same date + time)."""
    empty = df["station"] == ""
    path = Path(raw_csv_path).expanduser() if raw_csv_path else None
    if not empty.any() or path is None or not path.is_file():
        return df
    try:
        by_dm = read_note_iterations(str(path))[2]
    except Exception:
        return df
    new = [by_dm.get((d, int(round(h * 60)))) if e and not np.isnan(h) else None
           for d, h, e in zip(df["day"], df["hour"], empty)]
    df["station"] = [n or s for n, s in zip(new, df["station"])]
    print(f"  Station filled in from the raw CSV for {int((empty & (df['station'] != '')).sum())} row(s)")
    return df


def assign_iterations(df):
    """Characterization iteration of every row: from the raw-CSV notes if available, else chronological per station."""
    samp = df.drop_duplicates("sid")[["sid", "day", "hour", "station", "ts"]].sort_values(["ts", "sid"])
    it = None
    path = Path(raw_csv_path).expanduser() if raw_csv_path else None
    if path is not None and path.is_file():
        try:
            by_dts, by_ds, _ = read_note_iterations(str(path))
            it, miss = {}, []
            for r in samp.itertuples():
                v = None if np.isnan(r.hour) else by_dts.get((r.day, int(round(r.hour * 60)), r.station))
                if v is None:
                    v = by_ds.get((r.day, r.station))
                if v is None:
                    miss.append(f"{r.day:%d/%m} st {r.station}")
                    v = 1
                it[r.sid] = v
            print(f"Iterations read from the notes of {path.name}")
            if miss:
                print("  ! no matching raw row (set to iteration 1): " + "; ".join(miss))
        except Exception as e:
            print(f"  ! could not read iterations from {path} ({e}); inferring from the results files")
            it = None
    elif raw_csv_path:
        print(f"  ! raw CSV not found ({raw_csv_path}); inferring iterations from the results files")
    if it is None:
        samp = samp.copy()
        samp["iter"] = samp.groupby("station").cumcount() + 1
        samp.loc[samp["station"] == "", "iter"] = 1
        it = dict(zip(samp["sid"], samp["iter"]))
    return df["sid"].map(it).astype(int)


def load_results(folder):
    """-> DataFrame, one row per (sample, mode, method): phi, unc, ts (datetime), hour, station, sid, iter."""
    folder = Path(folder)
    files = sorted(folder.glob("grease_ice_results_*.csv"))
    if not files:
        raise FileNotFoundError(f"no grease_ice_results_*.csv in {folder}")
    parts = []
    for f in files:
        d = pd.read_csv(f, dtype={"Date": str, "Time": str, "Station": str}, keep_default_na=True)
        d = d[(d["Label"] == "derived") & (d["Notation"] == "phi_s")]
        if d.empty:
            continue
        parts.append(d[["Date", "Time", "Station", "Salinity_Mode", "Method", "Value", "Uncertainty"]])
    if not parts:
        raise ValueError("the results files contain no derived phi_s")
    df = pd.concat(parts, ignore_index=True)
    df.columns = ["date", "time", "station", "mode", "method", "phi", "unc"]
    df["station"] = df["station"].fillna("").str.strip()
    df["method"] = df["method"].astype(int)
    df["day"] = df["date"].map(parse_date)
    bad = df["day"].isna()
    if bad.any():
        print(f"  ! {int(bad.sum())} row(s) with an unreadable date dropped")
        df = df[~bad]
    df["hour"] = df["time"].map(parse_time)
    df = fill_missing_stations(df)
    # timestamp (set below): real time if given, otherwise 12:00 of that day
    df["unc"] = df["unc"].fillna(0.0)
    # stale files: results written before the 'time' column existed (empty time) for a day+station that now also
    # has timed results would be counted twice -> keep only the timed ones
    timed = set(zip(df.loc[df["hour"].notna(), "day"], df.loc[df["hour"].notna(), "station"]))
    stale = df["hour"].isna() & pd.Series([(d, s) in timed for d, s in zip(df["day"], df["station"])], index=df.index)
    if stale.any():
        print(f"  ! {int(stale.sum())} row(s) from old results files WITHOUT time ignored (same day+station now "
              f"exists with a time): delete the old files from the folder")
        df = df[~stale]
    df["ts"] = df["day"] + pd.to_timedelta(df["hour"].fillna(12.0), unit="h")
    df["sid"] = df["day"].dt.strftime("%Y-%m-%d") + "|" + df["time"].fillna("") + "|" + df["station"]

    df["iter"] = assign_iterations(df)

    n_s = df["sid"].nunique()
    print(f"Loaded {len(files)} file(s): {n_s} sample(s), {df['station'].nunique()} station(s), "
          f"{df['day'].nunique()} day(s)")
    multi = df.drop_duplicates("sid").groupby("station")["iter"].max()
    multi = multi[(multi > 1) & (multi.index != "")]
    if len(multi):
        print("  Stations characterized several times: " + "; ".join(f"{s} (x{n})" for s, n in multi.items()))
    else:
        print("  No station characterized more than once (all points use the colour of iteration 1)")
    if df["hour"].isna().all():
        print("  ! no sampling time in the files: all samples placed at 12:00 of their day")
    elif df["hour"].isna().any():
        no_t = df[df["hour"].isna()].drop_duplicates("sid")
        print("  ! sample(s) WITHOUT a time (placed at 12:00): "
              + "; ".join(f"{r.day:%d/%m} st {r.station}" for r in no_t.itertuples()))
        print("    (if the first data row lost its time when the header was renamed, put it back in the raw CSV)")
    return df.sort_values(["ts", "station"]).reset_index(drop=True)


# %% x positions, colours ----------------------------------------------------
def build_axis(df, xaxis):
    """Adds base x position 'x0' and a group key 'grp' (samples sharing the same x). Returns (df, ticks, labels, xlabel)."""
    df = df.copy()
    samples = df.drop_duplicates("sid")[["sid", "ts", "day", "station"]].sort_values(["ts", "station"])
    if xaxis == "day":
        days = sorted(samples["day"].unique())
        pos = {d: i for i, d in enumerate(days)}
        df["x0"], df["grp"] = df["day"].map(pos), df["day"]
        ticks, labels, xl = list(range(len(days))), [pd.Timestamp(d).strftime("%d/%m") for d in days], "Day (dd/mm)"
    elif xaxis == "station":
        order = list(dict.fromkeys(samples["station"]))        # stations in order of first characterisation
        pos = {s: i for i, s in enumerate(order)}
        df["x0"], df["grp"] = df["station"].map(pos), df["station"]
        ticks, labels, xl = list(range(len(order))), [s or "?" for s in order], "Station"
    elif xaxis == "sample":
        order = list(samples["sid"])
        pos = {s: i for i, s in enumerate(order)}
        df["x0"], df["grp"] = df["sid"].map(pos), df["sid"]
        lab = {r.sid: f"{r.station or '?'}\n{r.ts:%d/%m}" for r in samples.itertuples()}
        ticks, labels, xl = list(range(len(order))), [lab[s] for s in order], "Sample (chronological: station, dd/mm)"
    else:  # datetime
        df["x0"], df["grp"] = mdates.date2num(df["ts"]), df["sid"]
        ticks, labels, xl = None, None, "Date"
    return df, ticks, labels, xl


def dodge(df, xaxis, series_col, series_order):
    """x = x0 + offset so that overlapping points (several samples in the same x slot, several series) stay readable."""
    df = df.copy()
    n_ser = len(series_order)
    sj = df[series_col].map({s: j for j, s in enumerate(series_order)})
    if xaxis == "station":                      # one station = one x: series (methods / modes) are NOT spread in x
        n_ser, sj = 1, sj * 0
    if xaxis == "datetime":
        df["x"] = df["x0"] + (sj - (n_ser - 1) / 2) * 0.04            # days
        return df
    # rank of each sample inside its x slot (chronological)
    samp = df.drop_duplicates("sid")[["sid", "grp", "ts"]].sort_values("ts")
    samp["rank"] = samp.groupby("grp").cumcount()
    samp["n"] = samp.groupby("grp")["sid"].transform("count")
    if xaxis == "station" and not spread_iterations:
        samp["rank"], samp["n"] = 0, 1
    df = df.merge(samp[["sid", "rank", "n"]], on="sid", how="left")
    n_slots = df["n"] * n_ser
    slot = df["rank"] * n_ser + sj.values
    df["x"] = df["x0"] + (slot - (n_slots - 1) / 2) * (0.7 / n_slots)
    return df


def make_iteration_colours(df):
    """Discrete colour map: one colour per characterization iteration (1 .. max). -> cmap, norm, n_iter."""
    n_it = int(df["iter"].max())
    base = plt.get_cmap(cmap_name)
    cols = base(np.linspace(0.1, 0.9, n_it)) if n_it > 1 else [base(0.1)]
    return ListedColormap(cols), BoundaryNorm(np.arange(0.5, n_it + 1.5), n_it), n_it


# %% Plot helpers ------------------------------------------------------------
def draw(ax, d, marker, cmap, norm, cvals, errorbars, annotate=True):
    d = d.dropna(subset=["phi"])
    if d.empty:
        return
    c = cvals.loc[d.index].to_numpy()
    if errorbars:
        for x, y, u, cv in zip(d["x"], d["phi"], d["unc"], c):
            if u > 0:
                ax.errorbar(x, y, yerr=u, fmt="none", ecolor=cmap(norm(cv)), elinewidth=1.0, capsize=2,
                            alpha=0.55, zorder=2)
    ax.scatter(d["x"], d["phi"], c=c, cmap=cmap, norm=norm, marker=marker, s=60,
               edgecolor="k", linewidth=0.7, zorder=3)
    if label_stations and annotate:
        for x, y, st in zip(d["x"], d["phi"], d["station"]):
            if st:
                ax.annotate(st, (x, y), xytext=(0, 7), textcoords="offset points", ha="center", va="bottom",
                            fontsize=7, zorder=4, clip_on=True)


def finish_axes(axes, ticks, labels, xl, ylim, xaxis):
    for ax in np.ravel(axes):
        ax.grid(alpha=0.3)
        ax.set_ylabel(f"Solid fraction {PHI}")
        if ylim is not None:
            ax.set_ylim(*ylim)
        if ticks is not None:
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels, rotation=45 if xaxis != "sample" else 0, ha="right" if xaxis != "sample" else "center")
            ax.set_xlim(-0.6, len(ticks) - 0.4)
        else:
            ax.xaxis_date()
            loc = mdates.AutoDateLocator()
            ax.xaxis.set_major_locator(loc)
            ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))
        ax.set_xlabel(xl)
        ax.tick_params(axis="x", labelbottom=True)       # tick labels (stations...) on every subplot, not only the bottom


def add_colorbar(fig, axes, cmap, norm, n_it):
    sm = ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=np.ravel(axes).tolist(), pad=0.02, aspect=35, ticks=list(range(1, n_it + 1)))
    cb.set_label("Characterization iteration")
    cb.ax.set_yticklabels([str(i) for i in range(1, n_it + 1)])
    return cb


# %% Figures -----------------------------------------------------------------
def figure_modes_per_method(df, ctx):
    """Fig 1: one subplot per method; conductivity vs salinity-based on the same axes."""
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True, sharey=True, constrained_layout=True)
    dd = dodge(df, ctx["xaxis"], "mode", MODES)
    cvals = ctx["cvals"]
    for ax, m in zip(axes.ravel(), METHODS):
        sub = dd[dd["method"] == m]
        for mode in MODES:
            draw(ax, sub[sub["mode"] == mode], MODE_MARKER[mode], ctx["cmap"], ctx["norm"], cvals, ctx["errorbars"],
                 ctx["xaxis"] != "station")
        ax.set_title(METHOD_NAMES[m])
        if sub.empty:
            ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center", va="center", color="0.5")
    finish_axes(axes, ctx["ticks"], ctx["labels"], ctx["xl"], ctx["ylim"], ctx["xaxis"])
    handles = [Line2D([], [], marker=MODE_MARKER[k], ls="", mfc="0.7", mec="k", ms=8, label=MODE_LABEL[k]) for k in MODES]
    axes[0, 0].legend(handles=handles, loc="best", fontsize=8)
    add_colorbar(fig, axes, ctx["cmap"], ctx["norm"], ctx["n_it"])
    fig.suptitle(f"{PHI}: conductivity- vs salinity-based, per method")
    return fig


def figure_methods_per_mode(df, ctx):
    """Fig 2: all methods on the same axes; top = conductivity-based, bottom = salinity-based."""
    fig, axes = plt.subplots(2, 1, figsize=(11, 9), sharex=True, sharey=True, constrained_layout=True)
    dd = dodge(df, ctx["xaxis"], "method", METHODS)
    cvals = ctx["cvals"]
    for ax, mode in zip(axes, MODES):
        sub = dd[dd["mode"] == mode]
        for m in METHODS:
            draw(ax, sub[sub["method"] == m], METHOD_MARKER[m], ctx["cmap"], ctx["norm"], cvals, ctx["errorbars"],
                 ctx["xaxis"] != "station")
        ax.set_title(f"All methods, {MODE_LABEL[mode]}")
        if sub.empty:
            ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center", va="center", color="0.5")
    finish_axes(axes, ctx["ticks"], ctx["labels"], ctx["xl"], ctx["ylim"], ctx["xaxis"])
    handles = [Line2D([], [], marker=METHOD_MARKER[m], ls="", mfc="0.7", mec="k", ms=8, label=METHOD_NAMES[m]) for m in METHODS]
    axes[0].legend(handles=handles, loc="best", fontsize=8, ncol=2)
    add_colorbar(fig, axes, ctx["cmap"], ctx["norm"], ctx["n_it"])
    fig.suptitle(f"{PHI}: all methods, per salinity type")
    return fig


# %% Interactive choices -----------------------------------------------------
def ask_choice(prompt, options, default):
    keys = list(options)
    print(f"\n{prompt}")
    for i, k in enumerate(keys, 1):
        print(f"  [{i}] {k}  - {options[k]}" + ("   (default)" if k == default else ""))
    a = input("Choice [Enter = default]: ").strip()
    if not a:
        return default
    if a.isdigit() and 1 <= int(a) <= len(keys):
        return keys[int(a) - 1]
    return a if a in options else default


def get_settings():
    s = dict(defaults)
    if not interactive:
        return s
    s["xaxis"] = ask_choice("X axis for ALL plots:", {
        "day": "one slot per day", "station": "one slot per station (repeated characterisations side by side)",
        "sample": "one slot per sample, chronological", "datetime": "continuous date/time axis"}, s["xaxis"])
    s["errorbars"] = input("\nShow error bars (1 sigma)? [Y/n]: ").strip().lower() not in ("n", "no")
    a = input("y-limits as 'min max', 'auto', or Enter for 0..1 : ").strip().lower()
    if a == "auto":
        s["ylim"] = None
    elif a:
        try:
            lo, hi = (float(v) for v in a.replace(",", " ").split()[:2])
            s["ylim"] = (lo, hi)
        except ValueError:
            print("  could not read y-limits, using 0..1")
    return s


# %% Main --------------------------------------------------------------------
def main(folder=None):
    folder = folder or results_folder
    if folder is None or interactive and not Path(str(folder)).is_dir():
        folder = input("Folder with the grease_ice_results_*.csv files: ").strip().strip('"').strip("'")
    folder = Path(folder).expanduser()
    df = load_results(folder)
    st = get_settings()

    df, ticks, labels, xl = build_axis(df, st["xaxis"])
    cmap, norm, n_it = make_iteration_colours(df)
    ctx = dict(xaxis=st["xaxis"], ticks=ticks, labels=labels, xl=xl, ylim=st["ylim"], errorbars=st["errorbars"],
               cmap=cmap, norm=norm, n_it=n_it, cvals=df["iter"].astype(float))

    figs = {"phi_modes_per_method": figure_modes_per_method(df, ctx),
            "phi_methods_per_mode": figure_methods_per_mode(df, ctx)}
    if save_plot:
        out = Path(plot_output_folder).expanduser()
        out.mkdir(parents=True, exist_ok=True)
        for name, fig in figs.items():
            fn = out / f"{name}_{st['xaxis']}.png"
            fig.savefig(fn, dpi=200)
            print(f"saved {fn}")
    plt.show()
    return figs


if __name__ == "__main__":
    main()