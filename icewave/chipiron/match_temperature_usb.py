"""
match_temperature_usb.py -- match the Exelog temperature files to the Phonefleet Usb acquisitions by time,
list what is NOT matched, plot / save selected temperature files, and write a match CSV for quick_data_visu.py.

FOLDERS READ
    USB    ROOT_FOLDER/YYYY-MM-DD/HH-MM-SS/*Usb*.csv     (start = folder name, duration = phone timestamps, us)
    TEMP   TEMP_FOLDER/yy-mm-dd/DA-HH-MM.csv             (start = folder + file name, 1 min resolution; if the csv
                                                          holds real timestamps, those are used instead)

MATCHING (every Usb acquisition gets at most ONE temperature file)
    pass 1  start match   |temp start - usb start| <= START_THRESHOLD_S, closest pairs first, one-to-one
    pass 2  overlap       Usb acquisitions still without temperature: the temperature file whose time span
                          covers most of the acquisition (a long temperature file can serve several acquisitions)
    Anything left over is listed: Usb without temperature / temperature files without Usb.

OUTPUT   MATCH_CSV (one row per Usb acquisition and per unmatched temperature file), columns:
    status (matched / usb_only / temp_only), usb_date, usb_time, usb_folder, usb_start, usb_end, usb_duration_s,
    temp_path, temp_file, temp_column, temp_start, temp_end, delta_start_s (temp - usb), overlap_s, match_type

Then you are asked which temperature files to plot (x = time since the Usb start when matched, shaded = Usb window)
and which figures to save (matched: <usb folder>/figures/, otherwise ORPHAN_FIG_FOLDER).

The temperature csv format is auto-detected (delimiter, header, decimal comma, datetime / date + time / time-of-day
column). What was detected is printed as [temp] lines: check them once. quick_data_visu.py imports
read_temperature / plot_temperature from this file.

Run:  py match_temperature_usb.py
"""
import csv
import re
import sys
import warnings
from collections import Counter
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ======================= USER CONFIGURATION =======================
ROOT_FOLDER = Path("C:/Users/Gulliver/Desktop/Chipiron/git_clone_turbotice_icewave/icewave/icewave/chipiron/Phonefleet_data")
TEMP_FOLDER = ROOT_FOLDER / "Exelog_T_backup13112026"
MATCH_CSV = Path("usb_temperature_match.csv")      # read by quick_data_visu.py
START_THRESHOLD_S = 300         # |temp start - usb start| accepted as "same experiment" [s]
TEMP_CLOCK_OFFSET_S = 0.0       # added to the temperature times if the logger clock is shifted vs the phone [s]
TEMP_COLUMNS = None             # list of temperature column names to use, e.g. ["Channel 2", "Channel 3"]; None = all "Channel k" with data
TEMP_FALLBACK_DT_S = 1.0        # only if a csv has no time column: sample spacing assumed [s]
FIG_SUBFOLDER = "figures"       # matched figures go to <usb folder>/figures
ORPHAN_FIG_FOLDER = Path("temperature_figures")
FIG_DPI = 200
# ==================================================================

DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")
TIME_RE = re.compile(r"^\d{2}-\d{2}-\d{2}$")
NAME_RE = re.compile(r"(\d{2})-(\d{2})(?:-(\d{2}))?$")
TOD_RE = r"\d{1,2}:\d{2}(?::\d{2}(?:[.,]\d+)?)?"
CSV_COLUMNS = ["status", "usb_date", "usb_time", "usb_folder", "usb_start", "usb_end", "usb_duration_s",
               "temp_path", "temp_file", "temp_column", "temp_start", "temp_end", "delta_start_s", "overlap_s",
               "match_type"]


def ask(prompt, default=""):
    raw = input(f"{prompt}" + (f" [{default}]" if default else "") + ": ").strip()
    return raw or default


# ----------------------------------------------------------------------------- Usb acquisitions
def scan_usb(root):
    """List of dicts: date, time, folder, start, end, duration_s (from the phone timestamps of the Usb files)."""
    out = []
    for dd in sorted(p for p in Path(root).iterdir() if p.is_dir() and DATE_RE.match(p.name)):
        for td in sorted(p for p in dd.iterdir() if p.is_dir() and TIME_RE.match(p.name)):
            files = sorted(f for f in td.glob("*.csv") if "usb" in f.name.lower())
            if not files:
                continue
            start = datetime.strptime(f"{dd.name} {td.name}", "%Y-%m-%d %H-%M-%S")
            lo, hi = None, None
            for f in files:
                try:
                    t = pd.to_numeric(pd.read_csv(f, header=None, usecols=[0]).iloc[:, 0], errors="coerce")
                except Exception as exc:
                    print(f"[usb] {f}: unreadable ({exc})")
                    continue
                t = t[t > 0]
                if len(t):
                    lo = t.min() if lo is None else min(lo, t.min())
                    hi = t.max() if hi is None else max(hi, t.max())
            dur = float((hi - lo) / 1e6) if lo is not None else 0.0
            if lo is None:
                print(f"[usb] {td}: no valid timestamp, duration set to 0")
            out.append(dict(date=dd.name, time=td.name, folder=td, start=pd.Timestamp(start),
                            end=pd.Timestamp(start) + pd.Timedelta(seconds=dur), duration_s=dur))
    return out


# ----------------------------------------------------------------------------- temperature files
def _num(x):
    try:
        float(x.replace(",", "."))
        return True
    except ValueError:
        return False


def _read_rows(path):
    """Cells of the csv: delimiter guessed, rows with an odd number of fields (metadata lines) dropped."""
    text = None
    for enc in ("utf-8-sig", "latin-1"):
        try:
            text = Path(path).read_text(encoding=enc)
            break
        except UnicodeDecodeError:
            continue
    lines = [l for l in (text or "").splitlines() if l.strip()]
    if not lines:
        return []
    head = lines[:50]
    delim = ","
    for d in (";", "\t", ","):
        if sum(d in l for l in head) >= 0.7 * len(head):
            delim = d
            break
    rows = [[c.strip() for c in r] for r in csv.reader(lines, delimiter=delim)]
    n = Counter(len(r) for r in rows).most_common(1)[0][0]
    return [r for r in rows if len(r) == n]


def pick_temp_columns(df, cols):
    """Temperature columns: TEMP_COLUMNS if given, else the Excelog 'Channel k' columns (an open input reads 'High' and
    is not numeric: it is dropped), else any column with 'temp' in its name, else the first non-counter column."""
    if TEMP_COLUMNS:
        sel = [c for c in TEMP_COLUMNS if c in cols]
        if sel:
            return sel
    sel = [c for c in cols if re.match(r"(?i)^(channel|ch\s*\d|temp)", c)]
    if sel:
        return sel
    for c in cols:                                   # skip a record counter 1, 2, 3...
        v = df[c].dropna()
        if not (len(v) > 2 and (v.diff().dropna() == 1).all()):
            return [c]
    return cols[:1]


def read_temperature(path, date_hint=None, name_start=None):
    """Temperature csv -> DataFrame with 'time' (datetime64) + one float column per numeric column.
    df.attrs['columns'] = list of temperature columns, df.attrs['time_source'] = 'file' or 'index'.
    Handles the Calex Excelog layout: title / 'Date:' lines, repeated header blocks, a units row, 'High' = open input.
    date_hint / name_start (datetime) are only used when the csv holds a time of day without date, or no time."""
    rows = _read_rows(path)
    if len(rows) < 2:
        raise ValueError("no table found")
    nonempty = [c for c in rows[0] if c]
    header = sum(any(ch.isalpha() for ch in c) for c in nonempty) > len(nonempty) / 2
    names, body = (rows[0], rows[1:]) if header else ([f"col{i}" for i in range(len(rows[0]))], rows)
    seen, uniq = Counter(), []
    for i, n in enumerate(names):
        n = n or f"col{i}"
        seen[n] += 1
        uniq.append(n if seen[n] == 1 else f"{n}_{seen[n]}")
    df = pd.DataFrame(body, columns=uniq)
    df = df.loc[:, (df != "").any()]

    numeric = {}
    for c in df.columns:
        v = pd.to_numeric(df[c].str.replace(",", ".", regex=False), errors="coerce")
        if v.notna().mean() > 0.5:
            numeric[c] = v
    time, tod, dcol = None, None, None
    for c in df.columns:
        if c in numeric:
            continue
        s = df[c]
        if s.str.fullmatch(TOD_RE).mean() > 0.9:
            tod = tod or c
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p = pd.to_datetime(s, errors="coerce", dayfirst=True)
        if p.notna().mean() > 0.7:
            if (p.dt.hour.eq(0) & p.dt.minute.eq(0) & p.dt.second.eq(0)).all():
                dcol = dcol or p
            elif time is None:
                time = p
    if time is None and tod is not None:
        t = pd.to_timedelta(df[tod].str.replace(",", ".", regex=False).map(lambda x: x if x.count(":") == 2 else x + ":00"))
        if dcol is not None:
            time = dcol.dt.normalize() + t
        else:
            base = pd.Timestamp((date_hint or name_start or datetime.now()).date())
            day = (t.diff() < pd.Timedelta(0)).cumsum()          # midnight roll-over
            time = base + t + pd.to_timedelta(day, unit="D")
    source = "file"
    if time is None:
        base = pd.Timestamp(name_start or date_hint or datetime.now())
        time = base + pd.to_timedelta(np.arange(len(df)) * TEMP_FALLBACK_DT_S, unit="s")
        source = "index"
    out = pd.DataFrame({"time": pd.to_datetime(time).reset_index(drop=True) + pd.Timedelta(seconds=TEMP_CLOCK_OFFSET_S)})
    for c, v in numeric.items():
        out[c] = v.reset_index(drop=True)
    out = out.dropna(subset=["time"]).sort_values("time").reset_index(drop=True)   # junk rows (units, repeated headers) have no time
    cols = [c for c in numeric if out[c].notna().any()]
    if not cols or out.empty:
        raise ValueError("no numeric column / no valid time found")
    out.attrs.update(columns=pick_temp_columns(out, cols), time_source=source)
    return out


def temp_name_start(path):
    """Start time from <yy-mm-dd or yyyy-mm-dd folder>/<DA-HH-MM[-SS]>.csv, or None."""
    parent = path.parent.name
    try:
        m = re.fullmatch(r"(\d{2})-(\d{2})-(\d{2})", parent)
        if m:
            d = datetime(2000 + int(m[1]), int(m[2]), int(m[3]))
        elif DATE_RE.match(parent):
            d = datetime.strptime(parent, "%Y-%m-%d")
        else:
            return None
        m = NAME_RE.search(path.stem)
        return d.replace(hour=int(m[1]), minute=int(m[2]), second=int(m[3] or 0)) if m else None
    except ValueError:
        return None


def scan_temp(folder):
    """List of dicts: path, name_start, start, end, df, column, error."""
    out = []
    for p in sorted(Path(folder).rglob("*.csv")):
        ns = temp_name_start(p)
        rec = dict(path=p, name_start=ns, start=None, end=None, df=None, column="", columns=[], error="")
        try:
            df = read_temperature(p, date_hint=ns, name_start=ns)
            rec.update(df=df, columns=df.attrs["columns"], column="|".join(df.attrs["columns"]),
                       start=df["time"].min(), end=df["time"].max())
            msg = (f"[temp] {p.parent.name}/{p.name}: {len(df)} rows, column(s) {df.attrs['columns']}, "
                   f"times from the {'file' if df.attrs['time_source'] == 'file' else 'FILE NAME + row index (no time column found)'}, "
                   f"{rec['start']:%Y-%m-%d %H:%M:%S} -> {rec['end']:%H:%M:%S}")
            if ns is not None and abs((rec["start"] - pd.Timestamp(ns)).total_seconds()) > 3600:
                msg += f"   !! differs from the name ({ns:%Y-%m-%d %H:%M}) by more than 1 h"
            print(msg)
        except Exception as exc:
            rec["error"] = str(exc)
            if ns is not None:
                rec["start"] = pd.Timestamp(ns)
                rec["end"] = pd.Timestamp(ns)
            print(f"[temp] {p.parent.name}/{p.name}: cannot read ({exc})" + ("; name start only" if ns else ""))
        if rec["start"] is not None:
            out.append(rec)
    return out


# ----------------------------------------------------------------------------- matching
def overlap_s(a0, a1, b0, b1):
    return max(0.0, (min(a1, b1) - max(a0, b0)).total_seconds())


def match(usbs, temps):
    """{usb index: (temp index, 'start' | 'overlap')}"""
    pairs = []
    for i, u in enumerate(usbs):
        for j, t in enumerate(temps):
            d = abs((t["start"] - u["start"]).total_seconds())
            if d <= START_THRESHOLD_S:
                pairs.append((d, i, j))
    pairs.sort()
    res, used = {}, set()
    for _, i, j in pairs:
        if i not in res and j not in used:
            res[i] = (j, "start")
            used.add(j)
    for i, u in enumerate(usbs):
        if i in res:
            continue
        best = None
        for j, t in enumerate(temps):
            ov = overlap_s(u["start"], u["end"], t["start"], t["end"])
            if ov > 0 and (best is None or ov > best[0]):
                best = (ov, j)
        if best:
            res[i] = (best[1], "overlap")
    return res


def build_rows(usbs, temps, res):
    rows, matched_t = [], set()
    for i, u in enumerate(usbs):
        r = dict.fromkeys(CSV_COLUMNS, "")
        r.update(usb_date=u["date"], usb_time=u["time"], usb_folder=str(u["folder"]), usb_start=u["start"].isoformat(),
                 usb_end=u["end"].isoformat(), usb_duration_s=f"{u['duration_s']:.1f}", status="usb_only")
        if i in res:
            j, kind = res[i]
            t = temps[j]
            matched_t.add(j)
            r.update(status="matched", temp_path=str(t["path"]), temp_file=t["path"].name, temp_column=t["column"],
                     temp_start=t["start"].isoformat(), temp_end=t["end"].isoformat(), match_type=kind,
                     delta_start_s=f"{(t['start'] - u['start']).total_seconds():.0f}",
                     overlap_s=f"{overlap_s(u['start'], u['end'], t['start'], t['end']):.0f}")
        rows.append(r)
    for j, t in enumerate(temps):
        if j not in matched_t:
            r = dict.fromkeys(CSV_COLUMNS, "")
            r.update(status="temp_only", temp_path=str(t["path"]), temp_file=t["path"].name, temp_column=t["column"],
                     temp_start=t["start"].isoformat(), temp_end=t["end"].isoformat())
            rows.append(r)
    return rows


# ----------------------------------------------------------------------------- plotting
def plot_temperature(df, cols, usb_start=None, usb_end=None, title="", figsize=(10, 3.2)):
    """Temperature vs time. With usb_start: x = seconds since the Usb start, Usb window shaded. Returns the figure."""
    fig, ax = plt.subplots(figsize=figsize)
    if usb_start is not None:
        usb_start = pd.Timestamp(usb_start)
        x = (df["time"] - usb_start).dt.total_seconds()
        ax.set_xlabel("time since Usb start (s)")
        if usb_end is not None:
            ax.axvspan(0, (pd.Timestamp(usb_end) - usb_start).total_seconds(), color="tab:orange", alpha=0.2,
                       label="Usb acquisition")
    else:
        x = df["time"]
        ax.set_xlabel("time")
        fig.autofmt_xdate()
    cols = [cols] if isinstance(cols, str) else list(cols)
    for c in cols:
        ax.plot(x, df[c], ".-", ms=3, lw=0.8, label=c)
    ax.set_ylabel("temperature (°C)" if all(c.lower().startswith("ch") for c in cols) else ", ".join(cols))
    if len(cols) > 1 or usb_start is not None:
        ax.legend(fontsize=7, ncol=4)
    if title:
        ax.set_title(title, fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig


def parse_selection(raw, n):
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
    return sorted(i for i in picked if 1 <= i <= n)


def save_fig(fig, folder, name):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / name
    k = 2
    while path.exists():                                # never overwrite
        path = folder / f"{Path(name).stem}_{k}{Path(name).suffix}"
        k += 1
    fig.savefig(path, dpi=FIG_DPI, bbox_inches="tight")
    return path


# ----------------------------------------------------------------------------- main
def main():
    for p in (ROOT_FOLDER, TEMP_FOLDER):
        if not Path(p).is_dir():
            raise SystemExit(f"Folder not found: {p} (edit the configuration at the top of the script)")
    usbs = scan_usb(ROOT_FOLDER)
    print(f"\n{len(usbs)} Usb acquisition(s) found under {ROOT_FOLDER}\n")
    temps = scan_temp(TEMP_FOLDER)
    print(f"\n{len(temps)} temperature file(s) found under {TEMP_FOLDER}")
    if not usbs or not temps:
        raise SystemExit("Nothing to match.")

    res = match(usbs, temps)
    rows = build_rows(usbs, temps, res)

    print(f"\n=== MATCHED ({len(res)}) === (threshold on the start: {START_THRESHOLD_S} s)")
    for i, (j, kind) in sorted(res.items()):
        u, t = usbs[i], temps[j]
        print(f"  usb {u['date']} {u['time']} ({u['duration_s']:.0f} s)  <->  {t['path'].parent.name}/{t['path'].name}"
              f"   start delta {(t['start'] - u['start']).total_seconds():+.0f} s   [{kind}]")
    no_t = [u for i, u in enumerate(usbs) if i not in res]
    print(f"\n=== Usb WITHOUT temperature ({len(no_t)}) ===")
    for u in no_t:
        print(f"  {u['date']} {u['time']}  ({u['duration_s']:.0f} s)")
    used_t = {j for j, _ in res.values()}
    no_u = [t for j, t in enumerate(temps) if j not in used_t]
    print(f"\n=== temperature files WITHOUT Usb ({len(no_u)}) ===")
    for t in no_u:
        print(f"  {t['path'].parent.name}/{t['path'].name}  start {t['start']:%Y-%m-%d %H:%M:%S}"
              + (f"   (unreadable: {t['error']})" if t["error"] else ""))

    pd.DataFrame(rows, columns=CSV_COLUMNS).to_csv(MATCH_CSV, index=False)
    print(f"\nMatch file written: {Path(MATCH_CSV).resolve()}")

    # ---- choose temperature files to plot
    readable = [(j, t) for j, t in enumerate(temps) if t["df"] is not None]
    if not readable:
        return
    usb_of = {j: i for i, (j, _) in res.items()}
    print("\nTemperature files:")
    for k, (j, t) in enumerate(readable, 1):
        u = usbs[usb_of[j]] if j in usb_of else None
        print(f"  [{k}] {t['path'].parent.name}/{t['path'].name}   "
              + (f"-> usb {u['date']} {u['time']}" if u else "(no usb)"))
    chosen = parse_selection(ask("Plot which ones? numbers e.g. 1 3 / 2-4 / all / none", "none"), len(readable))
    figs = []
    plt.ion()
    for k in chosen:
        j, t = readable[k - 1]
        u = usbs[usb_of[j]] if j in usb_of else None
        title = f"{t['path'].parent.name}/{t['path'].name}" + (f"   (usb {u['date']} {u['time']})" if u else "")
        fig = plot_temperature(t["df"], t["columns"], u["start"] if u else None, u["end"] if u else None, title)
        figs.append((k, fig, u, t))
        plt.show(block=False)
        plt.pause(0.1)
    if figs:
        s = parse_selection(ask("Save which figures? numbers (as above) / all / none", "none"), len(readable))
        for k, fig, u, t in figs:
            if k not in s:
                continue
            if u:
                path = save_fig(fig, u["folder"] / FIG_SUBFOLDER, f"{u['date']}_{u['time']}_temperature_vs_t_s.png")
            else:
                path = save_fig(fig, ORPHAN_FIG_FOLDER, f"temperature_{t['path'].parent.name}_{t['path'].stem}.png")
            print(f"  saved {path}")
        try:
            get_ipython()      # noqa: F821
        except NameError:
            print("Figures stay open: close the windows to end the script.")
            plt.ioff()
            plt.show()


if __name__ == "__main__":
    main()