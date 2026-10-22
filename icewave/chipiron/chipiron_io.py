"""
chipiron_io.py -- shared I/O for general data structure.   NOT TESTED ON REAL h5py / REAL DATA YET.

THREE KINDS OF FILES per experiment, tied together by the METADATA txt file:

    <stem>_meta.txt   human-readable "hub": experiment info, parameters (with units/uncertainties)
                      and, in its [links] section, the names of the 3 files:
                          meta_file   = <stem>_meta.txt
                          rheo_file   = <stem>_rheo.h5
                          sample_file = <sample1>.csv, <sample2>.csv   (or 'pending')
    <stem>_rheo.h5    raw rheology (HDF5): one dataset per quantity, attrs 'unit' and 'uncertainty'
                      (a scalar attr if constant, else a companion dataset '<name>_unc')
    <sample>.csv      sample characterisation (phi, salinity, beam properties, ...), "long" format
                      Parameter Name | Notation | Value | Uncertainty | Units | Label | Method
                      (the format written by the phi-estimation script); a "wide" CSV also works.

Metadata line format (anything after '#' is a comment):
        key = value | unit | uncertainty          e.g.   freq = 0.2 | Hz | 0.001
        key = lo .. hi | unit                      e.g.   omega = 10 .. 800 | rpm     (a varied parameter)
        key = text                                 e.g.   campaign = CASCADES26
Sections ([links], [experiment], [parameters], ...) are only for readability: keys must be unique.
"""
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

FORMAT_VERSION = "1"
LINK_KINDS = ("meta", "rheo", "sample")
PENDING = {"", "pending", "none", "nan"}


# ----------------------------------------------------------------------------- basics
@dataclass
class Q:
    """A quantity: value (float / str / numpy array), unc (None, scalar, or array like value), unit."""
    value: object
    unc: object = None
    unit: str = ""

    def unc_array(self):
        v = np.asarray(self.value, float)
        if self.unc is None:
            return np.full(v.shape, np.nan)
        return np.broadcast_to(np.asarray(self.unc, float), v.shape)


def slug(s, maxlen=60):
    """Safe for file names: letters, digits, '-', '_', '.' only."""
    s = re.sub(r"[^\w\-.]+", "-", str(s).strip()).strip("-_.")
    return s[:maxlen]


def _num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return str(s).strip()


# ----------------------------------------------------------------------------- parameter strings
_P = re.compile(r"^\s*(?P<name>[A-Za-z_]\w*)\s*(?:\[(?P<unit>[^\]]*)\])?\s*=\s*(?P<val>.+?)\s*$")


def parse_params(text):
    """Parse 'omega[rpm]=10..800, freq[Hz]=0.2+-0.001, T[degC]=-1.8' ->
    [dict(name, unit, value (float | str | (lo, hi)), unc)].  Separators: ',' or ';'."""
    out = []
    for chunk in re.split(r"[;,]", text or ""):
        if not chunk.strip():
            continue
        m = _P.match(chunk)
        if not m:
            print(f"  ! cannot parse parameter '{chunk.strip()}' (expected name[unit]=value) -> ignored")
            continue
        val, unc = m["val"], None
        if ".." in val:
            lo, hi = (x.strip() for x in val.split(".."))
            value = (_num(lo), _num(hi))
            if not all(isinstance(v, float) for v in value):
                value = val
        else:
            parts = re.split(r"\+-|\+/-|±", val)
            if len(parts) == 2:
                value, unc = _num(parts[0]), _num(parts[1])
                unc = unc if isinstance(unc, float) else None
            else:
                value = _num(val)
        out.append(dict(name=m["name"], unit=(m["unit"] or "").strip(), value=value, unc=unc))
    return out


def compact_params(params):
    """Short string for file names: omega10-800rpm_freq0.2Hz"""
    bits = []
    for p in params:
        v = p["value"]
        s = f"{v[0]:g}-{v[1]:g}" if isinstance(v, tuple) else (f"{v:g}" if isinstance(v, float) else str(v))
        bits.append(f"{p['name']}{s}{p['unit']}")
    return slug("_".join(bits), 80)


# ----------------------------------------------------------------------------- metadata txt
def meta_line(key, value, unit="", unc=None):
    if isinstance(value, tuple):
        vs = f"{value[0]:.10g} .. {value[1]:.10g}"
    elif isinstance(value, (int, float, np.number)) and not isinstance(value, bool):
        vs = f"{float(value):.10g}"
    else:
        vs = str(value)
    parts = [vs]
    if unit or unc is not None:
        parts.append(unit)
    if unc is not None:
        parts.append(f"{float(unc):.6g}")
    return f"{key} = " + " | ".join(parts)


def write_meta(path, sections):
    """sections: {section: [(key, value, unit, unc), ...]} (unit/unc optional, tuples may be shorter)."""
    lines = ["# chipiron metadata v1",
             "# key = value | unit | uncertainty      (range: lo .. hi)",
             "# [links] = names of the meta / rheo / sample files of this experiment "
             "(kept up to date by chipiron_io.register_file)"]
    for sec, items in sections.items():
        lines += ["", f"[{sec}]"] + [meta_line(*it) for it in items]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_meta(path):
    """-> {key: dict(value, unit, unc, section)}  (value: float | str | (lo, hi))."""
    items, sec = {}, ""
    for ln in Path(path).read_text(encoding="utf-8").splitlines():
        ln = ln.split("#", 1)[0].strip()
        if not ln:
            continue
        if ln.startswith("[") and ln.endswith("]"):
            sec = ln[1:-1].strip().lower()
            continue
        if "=" not in ln:
            continue
        key, rest = ln.split("=", 1)
        parts = [p.strip() for p in rest.split("|")]
        v = parts[0]
        if re.fullmatch(r"\S+\s*\.\.\s*\S+", v) and isinstance(_num(v.split("..")[0]), float):
            lo, hi = (x.strip() for x in v.split(".."))
            value = (_num(lo), _num(hi))
        else:
            value = _num(v)
        unit = parts[1] if len(parts) > 1 else ""
        unc = _num(parts[2]) if len(parts) > 2 and parts[2] else None
        if key.strip() in items:
            print(f"  ! duplicate key '{key.strip()}' in {Path(path).name}: last one kept")
        items[key.strip()] = dict(value=value, unit=unit, unc=unc if isinstance(unc, float) else None, section=sec)
    return items


def meta_to_row(items):
    """Flatten metadata for a table: ranges -> key_min, key_max (+ key = midpoint). Returns (row, units)."""
    row, units = {}, {}
    for k, d in items.items():
        v = d["value"]
        if isinstance(v, tuple):
            row[k + "_min"], row[k + "_max"] = v
            row[k] = 0.5 * (v[0] + v[1])
            units[k + "_min"] = units[k + "_max"] = d["unit"]
        else:
            row[k] = v
        units[k] = d["unit"]
        if d["unc"] is not None:
            row[k + "_unc"] = d["unc"]
    return row, units


def register_file(meta_path, kind, filename, verbose=True):
    """Make sure `filename` is listed in the [links] section of the metadata file as `<kind>_file`.
    kind in meta / rheo / sample. 'sample' accepts several files (comma list); the others are single:
    a different existing value is NOT overwritten (warning).  Returns True if the link is in place."""
    assert kind in LINK_KINDS, kind
    key, filename, p = f"{kind}_file", Path(filename).name, Path(meta_path)
    if not p.exists():
        print(f"[link] ERROR: metadata file not found: {p}")
        return False
    lines = p.read_text(encoding="utf-8").splitlines()

    def find_section():
        s = next((i for i, l in enumerate(lines) if l.strip().lower() == "[links]"), None)
        if s is None:
            return None, None
        e = next((i for i in range(s + 1, len(lines)) if lines[i].strip().startswith("[")), len(lines))
        return s, e

    s, e = find_section()
    if s is None:                                   # no [links] section yet: create it after the header comments
        top = next((i for i, l in enumerate(lines) if l.strip() and not l.strip().startswith("#")), len(lines))
        lines[top:top] = ["[links]"] + [f"{k}_file = pending" for k in LINK_KINDS] + [""]
        s, e = find_section()
    idx = next((i for i in range(s + 1, e) if lines[i].split("=")[0].strip() == key), None)
    cur = lines[idx].split("=", 1)[1].split("#")[0].strip() if idx is not None else ""
    entries = [x.strip() for x in cur.split(",") if x.strip().lower() not in PENDING]
    if filename in entries:
        if verbose:
            print(f"[link] {p.name}: {key} already contains {filename}")
        return True
    if entries and kind != "sample":
        print(f"[link] WARNING {p.name}: {key} = {', '.join(entries)} (not overwritten by {filename})")
        return False
    new = f"{key} = " + ", ".join(entries + [filename])
    if idx is None:
        lines.insert(s + 1, new)
    else:
        lines[idx] = new
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    if verbose:
        print(f"[link] {p.name}: added {key} = {filename}")
    return True


# ----------------------------------------------------------------------------- HDF5 rheology
def _s(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def write_rheo_h5(path, quantities, attrs=None):
    """quantities: {name: Q}. One dataset per quantity; attrs unit + uncertainty (scalar attr if constant,
    else companion dataset '<name>_unc'). MATLAB: h5read(f,'/current_A'); h5readatt(f,'/current_A','unit')."""
    import h5py
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        f.attrs["format_version"] = FORMAT_VERSION
        for k, v in (attrs or {}).items():
            f.attrs[k] = v
        for name, q in quantities.items():
            val = np.asarray(q.value)
            ds = f.create_dataset(name, data=val, compression="gzip" if val.ndim > 0 and val.size > 1 else None)
            ds.attrs["unit"] = q.unit
            if q.unc is not None:
                u = np.asarray(q.unc, float)
                if u.ndim == 0 or (u.size > 0 and np.all(u == u.flat[0])):
                    ds.attrs["uncertainty"] = float(u.flat[0])
                else:
                    f.create_dataset(name + "_unc", data=u, compression="gzip")
                    ds.attrs["uncertainty_dataset"] = name + "_unc"


def read_rheo_h5(path, names=None):
    """-> ({name: Q}, file_attrs). `names` restricts what is read (large files)."""
    import h5py
    out = {}
    with h5py.File(path, "r") as f:
        attrs = {k: (_s(v) if isinstance(v, (bytes, str)) else v) for k, v in f.attrs.items()}
        companions = {_s(d.attrs["uncertainty_dataset"]) for d in f.values()
                      if isinstance(d, h5py.Dataset) and "uncertainty_dataset" in d.attrs}
        for name, d in f.items():
            if not isinstance(d, h5py.Dataset) or name in companions or (names and name not in names):
                continue
            unc = None
            if "uncertainty_dataset" in d.attrs:
                unc = f[_s(d.attrs["uncertainty_dataset"])][()]
            elif "uncertainty" in d.attrs:
                unc = float(d.attrs["uncertainty"])
            val = d[()]
            out[name] = Q(val, unc, _s(d.attrs.get("unit", "")))
    return out, attrs


# ----------------------------------------------------------------------------- sample CSV
def read_sample_csv(path, prefer_method=None):
    """-> {name: Q(float, unc, unit)}.  Flexible on purpose:
    * long format (columns Value + Notation/Parameter Name/Name, optional Uncertainty, Units, Method)
    * wide format (one column per quantity, header 'name [unit]', optional '<name>_unc' column; row 0)
    Quantities that come with a 'Method' are stored as 'phi_s@M1', 'phi_s@M2', ...; the plain name 'phi_s'
    is an alias of `prefer_method` (or of the only method present)."""
    import pandas as pd
    df = pd.read_csv(path)
    cols = {c.lower().strip(): c for c in df.columns}
    out, methods = {}, {}
    if "value" in cols:
        ncol = cols.get("notation") or cols.get("name") or cols.get("parameter name") or df.columns[0]
        ucol, ucol2, mcol = cols.get("uncertainty"), cols.get("units") or cols.get("unit"), cols.get("method")
        for _, r in df.iterrows():
            nm = str(r[ncol]).strip()
            val = r[cols["value"]]
            val = float(val) if _is_num(val) else str(val)
            unc = float(r[ucol]) if ucol and _is_num(r[ucol]) and not pd.isna(r[ucol]) else None
            q = Q(val, unc, "" if not ucol2 or pd.isna(r[ucol2]) else str(r[ucol2]))
            if mcol and not pd.isna(r[mcol]):
                m = int(r[mcol])
                out[f"{nm}@M{m}"] = q
                methods.setdefault(nm, []).append(m)
            else:
                out[nm] = q
        for nm, ms in methods.items():
            if nm in out:
                continue
            if prefer_method is not None and prefer_method in ms:
                out[nm] = out[f"{nm}@M{prefer_method}"]
            elif len(ms) == 1:
                out[nm] = out[f"{nm}@M{ms[0]}"]
            else:
                print(f"  ! {Path(path).name}: '{nm}' exists for methods {ms}: use '{nm}@M<k>' or set prefer_method")
    else:
        for c in df.columns:
            m = re.match(r"^(.*?)\s*\[([^\]]*)\]\s*$", c)
            nm, unit = (m.group(1), m.group(2)) if m else (c.strip(), "")
            if nm.endswith("_unc"):
                continue
            v = df[c].to_numpy()
            v = v.astype(float) if np.issubdtype(v.dtype, np.number) else v
            uc = next((df[k].to_numpy(float) for k in (nm + "_unc", nm + "_unc [" + unit + "]") if k in df.columns), None)
            if getattr(v, "size", 0) == 1 or len(v) == 1:
                v, uc = v[0], (uc[0] if uc is not None else None)
            out[nm] = Q(v, uc, unit)
    return out


def _is_num(x):
    try:
        float(x)
        return True
    except (TypeError, ValueError):
        return False


# ----------------------------------------------------------------------------- file lookup
def index_files(root):
    """{basename: Path} of every file under root (warns about duplicated names)."""
    idx = {}
    for p in Path(root).rglob("*"):
        if p.is_file():
            if p.name in idx:
                print(f"  ! duplicate file name {p.name}: {idx[p.name]} and {p}")
            idx[p.name] = p
    return idx


def choose_meta_file(root, ask=input):
    """List metadata files (newest first) and let the user pick one (index or path). None = skip."""
    metas = sorted(Path(root).rglob("*_meta.txt"), key=lambda p: p.stat().st_mtime, reverse=True)
    if not metas:
        print(f"No *_meta.txt found under {root}")
    for i, m in enumerate(metas):
        print(f"  [{i}] {m.name}")
    a = ask("Link to which metadata file? index / path / Enter = skip: ").strip()
    if not a:
        return None
    if a.isdigit() and int(a) < len(metas):
        return metas[int(a)]
    return Path(a) if Path(a).exists() else (print(f"  ! not found: {a}") or None)
