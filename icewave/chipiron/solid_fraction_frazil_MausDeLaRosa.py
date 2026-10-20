"""
Solid fraction and salinity of grease ice: Methods 1 to 4 of
Maus & De la Rosa (2012), 'Salinity and solid fraction of frazil and grease ice', J. Glaciol. 58(209), 594-612,
plus Method 0 (PSS-78: conductivity + temperature -> salinity).

WHAT IT DOES
    Asks for a folder, lets you pick a raw CSV in it (one row = one sample), runs every method the available
    data allow on every row (missing/NaN values simply switch off the methods that need them), writes one results
    CSV per sample (+ one summary table) and optionally links the results files to an experiment's *_meta.txt.

EXPECTED INPUT CSV (the "template")
    * one row per sample; header names are case-sensitive and use the paper's notation (underscore for subscript,
      "'" for primes). Separator , or ;  and decimal point or comma are both accepted (French Excel export is fine).
    * for every quantity X there is an OPTIONAL uncertainty column  X_unc  (1 standard deviation, same unit as X).
      Empty/NaN X_unc = no uncertainty (it then contributes 0 to the propagated uncertainty).
    * empty cells, NaN, NA, n/a, '-' are all read as "not measured". Extra columns (e.g. notes) are ignored.

      column            unit       meaning
      ----------------  ---------  ----------------------------------------------------------------------------
      date              text       sampling date, e.g. 2026-09-22 or 22/09/2026 (day first)      [not in the paper]
      time              text       sampling time, e.g. 14:30 or 14:30:15                           [not in the paper]
      T_cond_g          degC       temperature of the MELTED grease when sigma_g was measured      [not in the paper]
      sigma_g           mS/cm      electrical conductivity of the melted grease                    [not in the paper]
      T_cond_b          degC       temperature of the brine when sigma_b was measured              [not in the paper]
      sigma_b           mS/cm      electrical conductivity of the brine                            [not in the paper]
      S_g               -          bulk grease salinity, if you already have it (replaces sigma_g/T_cond_g)
      S_b               -          brine salinity, if you already have it (replaces sigma_b/T_cond_b)
      S_w               -          sea-water salinity: used as S_b ~ S_w ONLY if neither S_b nor sigma_b/T_cond_b
                                   is available (the paper's sea-water approximation)
      V_g               m^3        grease sample volume before drainage
      M_g               kg         grease sample mass before drainage        (rho_g = M_g / V_g)
      M'_g              kg         mass of the drained grease ice (Method 4)
      S'_g              -          salinity of the drained grease ice (Method 4)
      phi_a             -          air volume fraction (used only where include_air is True)
    Salinities are dimensionless (practical salinity ~ g/kg), as in the paper.

    What each method needs (anything missing -> that method is skipped for that row, the others still run):
      Method 1: S_g, S_b                     (+ M_g, V_g, phi_a if its air correction is on)
      Method 2: S_b, M_g, V_g                (+ phi_a if air on)
      Method 3: S_g, M_g, V_g                (+ phi_a if air on)
      Method 4: S_b, V_g, M'_g, S'_g         (+ phi_a if air on)
      S_g comes from  S_g  or else PSS-78(sigma_g, T_cond_g);
      S_b comes from  S_b  or else PSS-78(sigma_b, T_cond_b) or else S_w.
    If include_air is on for a method but phi_a is empty for that row, phi_a_default is used (flagged in the output).

OUTPUT
    grease_ice_results_<YYYYMMDD_HHMMSS>.csv   one per sample, long format (read by chipiron_io.read_sample_csv):
        Date | Time | Parameter Name | Notation | Value | Uncertainty | Units | Label | Method | Source |
        Saved_Method_Flag | Air_Included | Air_Fraction_Used | Air_Source
        Label: 'measured' (your inputs), 'constant', 'computed' (S_g, S_b, rho_g derived from inputs),
               'derived' (Methods 1-4, with the Method number).
        Source: where S_g / S_b came from ('input', 'PSS-78', 'S_w approximation').
        Uncertainty 0 = no input uncertainty was given.
    grease_ice_summary_<input file name>.csv    one wide row per sample, with the reason for every skipped method.
"""

# %% Initialization
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from uncertainties import ufloat
from uncertainties.umath import sqrt

try:
    import chipiron_io as cio
except ImportError:
    cio = None

# --- User configuration ---------------------------------------------------
raw_data_folder = None      # e.g. r"C:\Users\me\PMMH\raw_grease" ; None = ask when the script runs
output_folder = None        # None = same folder as the input CSV
saved_method = 'all'        # 1, 2, 3, 4 or 'all': which methods go in the per-sample results CSV
link_to_metadata = True     # offer to register the results files in an experiment *_meta.txt
meta_root = "PhD_data"      # where chipiron_io looks for *_meta.txt

# Air correction per method (True = air-corrected equation, False = air-free equation), or one bool for all.
# NOTE Method 4: by our own derivation Eq 24 is already exact when V_g is the total sample volume, so the paper's
# correction is implemented as published but defaults to False. Set True to follow the paper literally.
include_air = {1: True, 2: True, 3: True, 4: False}
phi_a_default = (0.01, 0.005)   # (value, uncertainty) used when a row has no phi_a but air is switched on; None = skip instead

# Physical constants (Maus & De la Rosa 2012)
rho_i = 917.0     # Pure ice density (kg m^-3)
rho_0 = 1000.0    # Reference fresh-water density (kg m^-3)
beta = 0.81       # Haline contraction coefficient (kg m^-3 per unit salinity); 0.77 for NaCl solutions
# (the paper prints "g m^-3" for beta, which is a typo for kg m^-3)

# Input columns: (name, unit, description)
INPUTS = [
    ("T_cond_g", "°C", "Temperature of melted grease at conductivity measurement"),
    ("sigma_g", "mS cm^-1", "Electrical conductivity of melted grease"),
    ("T_cond_b", "°C", "Temperature of brine at conductivity measurement"),
    ("sigma_b", "mS cm^-1", "Electrical conductivity of brine"),
    ("S_g", "psu", "Grease bulk salinity (direct input)"),
    ("S_b", "psu", "Brine salinity (direct input)"),
    ("S_w", "psu", "Sea-water salinity (S_b ~ S_w approximation)"),
    ("V_g", "m^3", "Grease sample volume before drainage"),
    ("M_g", "kg", "Grease sample mass before drainage"),
    ("M'_g", "kg", "Drained grease ice mass"),
    ("S'_g", "psu", "Drained grease ice salinity"),
    ("phi_a", "dimensionless", "Air volume fraction"),
]
POSITIVE = ["sigma_g", "sigma_b", "S_g", "S_b", "S_w", "V_g", "M_g", "M'_g", "S'_g"]


# %% Helpers
def _nv(x):
    return getattr(x, 'nominal_value', x)


def _sd(x):
    return getattr(x, 'std_dev', 0.0)


def _mk(v, u):
    """ufloat if an uncertainty is given, else a plain float (avoids std_dev == 0 ufloats)."""
    return ufloat(v, u) if (u is not None and np.isfinite(u) and u > 0) else float(v)


def _val(row, name):
    """Value (+ optional '<name>_unc') from a row dict -> float | ufloat | None (None = missing/NaN)."""
    try:
        v = float(row.get(name, np.nan))
    except (TypeError, ValueError):
        return None
    if not np.isfinite(v):
        return None
    try:
        u = float(row.get(name + "_unc", np.nan))
    except (TypeError, ValueError):
        u = np.nan
    return _mk(v, u)


def _require(cond, msg):
    if not cond:
        raise ValueError(msg)


def _warn_range(name, x, lo, hi):
    if not (lo <= _nv(x) <= hi):
        warnings.warn(f"{name} = {_nv(x):.4g} outside expected range [{lo}, {hi}]")


def _air_flags(inc):
    if isinstance(inc, bool):
        return {m: inc for m in (1, 2, 3, 4)}
    if set(inc) != {1, 2, 3, 4}:
        raise ValueError("include_air must be a bool or a dict with keys 1, 2, 3, 4")
    return dict(inc)


def slug(s, maxlen=60):
    return re.sub(r"[^\w\-.]+", "-", str(s).strip()).strip("-_.")[:maxlen]


# %% Method 0: Salinity from Conductivity and Temperature (PSS-78)
# Valid for seawater and sea ice brines. NaCl solutions: PSS-78 is a seawater scale, so the
# result is only approximately the NaCl mass fraction in g/kg (the paper treats the two as equal).

def pss78_salinity(C, T):
    """
    Calculate Practical Salinity from conductivity (mS/cm) and temperature (°C)
    at standard atmospheric pressure (p=0) using UNESCO 1983 (PSS-78) polynomials.
    """
    # Standard conductivity of seawater at S=35, T=15°C, p=0 dbar
    C35150 = 42.914

    # Conductivity ratio R
    R = C / C35150

    # Temperature compensation polynomial r_T
    c0 = 0.6766097
    c1 = 2.00564e-2
    c2 = 1.104259e-4
    c3 = -6.9698e-7
    c4 = 1.0031e-9

    r_T = c0 + c1 * T + c2 * T**2 + c3 * T**3 + c4 * T**4

    # Temperature-corrected conductivity ratio R_T (for atmospheric pressure, R_p = 1)
    R_T = R / r_T
    sqrt_R_T = sqrt(R_T)

    # Salinity polynomial coefficients
    a0, a1, a2, a3, a4, a5 = 0.0080, -0.1692, 25.3851, 14.0941, -7.0261, 2.7081
    b0, b1, b2, b3, b4, b5 = 0.0005, -0.0056, -0.0066, -0.0375, 0.0636, -0.0144
    k = 0.0162

    # Main term
    S_main = a0 + a1 * sqrt_R_T + a2 * R_T + a3 * R_T * sqrt_R_T + a4 * R_T**2 + a5 * R_T**2 * sqrt_R_T

    # Temperature correction term
    dS = (T - 15.0) / (1.0 + k * (T - 15.0))
    S_corr = dS * (b0 + b1 * sqrt_R_T + b2 * R_T + b3 * R_T * sqrt_R_T + b4 * R_T**2 + b5 * R_T**2 * sqrt_R_T)

    return S_main + S_corr


# %% Methods 1 to 4 (one function each; they raise ValueError when the numbers are not physical)
def method1(S_g, S_b, rho_g, phi_a, air):
    """Eq 9 (air off) | Eq 12 with measured rho_g (air on)."""
    rho_b = rho_0 + beta * S_b                                            # Eq 10
    if air:
        phi_s = 1.0 - (S_g * rho_g) / (S_b * rho_b) - phi_a               # Eq 12
    else:
        phi_s = 1.0 - (1.0 + (S_b / S_g - 1.0) * (rho_b / rho_i)) ** (-1.0)   # Eq 9
    return dict(rho_b=rho_b, phi_s=phi_s, phi_b=1.0 - phi_s - phi_a)     # Eq 11


def method2(S_b, rho_g, phi_a, air):
    """Eq 14 (air off) | Eq 15 (air on)."""
    rho_b = rho_0 + beta * S_b
    _require(_nv(rho_b) != rho_i, "Method 2: rho_b equals rho_i")
    phi_s = (rho_b - rho_g - rho_b * phi_a) / (rho_b - rho_i)
    return dict(rho_b=rho_b, phi_s=phi_s, phi_b=1.0 - phi_s - phi_a)


def method3(S_g, rho_g, phi_a, air):
    """Eq 17-19 (air off) | Eq 20 -> Eq 17-18 with [rho_g - rho_i(1-phi_a)], then Eq 12 (air on)."""
    d_rho = rho_g - rho_i * (1.0 - phi_a)          # = (rho_g - rho_i) when phi_a = 0
    a = beta * d_rho
    b = rho_0 * d_rho - S_g * beta * rho_g
    c = S_g * rho_g * (rho_i - rho_0)
    _require(_nv(a) > 0, "Method 3: rho_g <= rho_i(1-phi_a), quadratic has no physical root")
    disc = b ** 2 - 4.0 * a * c
    _require(_nv(disc) >= 0, "Method 3: negative discriminant")
    S_b = (sqrt(disc) - b) / (2.0 * a)                                    # Eq 18
    rho_b = rho_0 + beta * S_b
    phi_s = 1.0 - (S_g * rho_g) / (S_b * rho_b) - phi_a                   # Eq 19 / Eq 12
    _warn_range("Method 3 brine salinity S_b", S_b, 20.0, 60.0)
    return dict(S_b=S_b, rho_b=rho_b, phi_s=phi_s, phi_b=1.0 - phi_s - phi_a)


def method4(S_b, V_g, M_prime_g, S_prime_g, phi_a, air):
    """Eq 24-25 (air off) | paper's air correction of Eq 24: V_g -> V_g(1-phi_a), minus phi_a (air on).
    Sg then uses the air-consistent form of Eq 25, Sg = Sb rho_b phi_b / (rho_i phi_s + rho_b phi_b)
    (reduces exactly to Eq 25 when phi_a = 0; this extension is not in the paper)."""
    rho_b = rho_0 + beta * S_b
    phi_s = (M_prime_g / (V_g * (1.0 - phi_a) * rho_i)) * (1.0 - S_prime_g / S_b) - phi_a   # Eq 24
    phi_b = 1.0 - phi_s - phi_a
    _require(_nv(phi_b) > 0 and _nv(phi_s) < 1.0, "Method 4: phi_s out of range (phi_b must be > 0)")
    S_g = S_b * rho_b * phi_b / (rho_i * phi_s + rho_b * phi_b)                              # Eq 25
    return dict(rho_b=rho_b, phi_s=phi_s, S_g=S_g, phi_b=phi_b)


# %% Analysis of ONE sample (a dict of raw values) with NaN handling
def analyse_sample(row):
    """-> (records, summary, log). records = list of dict rows for the long-format CSV."""
    air = _air_flags(include_air)
    log, records = [], []
    x = {n: _val(row, n) for n, _, _ in INPUTS}

    for n in POSITIVE:                                   # non-physical inputs are treated as missing
        if x[n] is not None and _nv(x[n]) <= 0:
            log.append(f"{n} = {_nv(x[n]):g} is not positive: ignored")
            x[n] = None
    if x["phi_a"] is not None and not (0.0 <= _nv(x["phi_a"]) < 1.0):
        log.append(f"phi_a = {_nv(x['phi_a']):g} outside [0, 1): ignored")
        x["phi_a"] = None

    def add(name, nt, val, unit, label, method=np.nan, source="", flag=np.nan, air_used=np.nan, air_src=""):
        records.append({'Parameter Name': name, 'Notation': nt, 'Value': _nv(val), 'Uncertainty': _sd(val),
                        'Units': unit, 'Label': label, 'Method': method if label == 'derived' else np.nan,
                        'Source': source, 'Saved_Method_Flag': (saved_method == 'all' or str(saved_method) == str(method)) if label == 'derived' else True,
                        'Air_Included': flag, 'Air_Fraction_Used': _nv(air_used) if flag is not np.nan else np.nan,
                        'Air_Source': air_src})

    # --- Method 0 / inputs -> S_g, S_b, rho_g
    S_g = S_b = None
    src_g = src_b = ""
    try:
        if x["S_g"] is not None:
            S_g, src_g = x["S_g"], "input"
        elif x["sigma_g"] is not None and x["T_cond_g"] is not None:
            _warn_range("T_cond_g", x["T_cond_g"], -2.0, 35.0)
            S_g, src_g = pss78_salinity(x["sigma_g"], x["T_cond_g"]), "PSS-78"
            _warn_range("S_g (PSS-78)", S_g, 2.0, 42.0)
        else:
            log.append("S_g unavailable (need S_g or sigma_g + T_cond_g)")
        if x["S_b"] is not None:
            S_b, src_b = x["S_b"], "input"
        elif x["sigma_b"] is not None and x["T_cond_b"] is not None:
            _warn_range("T_cond_b", x["T_cond_b"], -2.0, 35.0)
            S_b, src_b = pss78_salinity(x["sigma_b"], x["T_cond_b"]), "PSS-78"
            _warn_range("S_b (PSS-78)", S_b, 2.0, 42.0)
        elif x["S_w"] is not None:
            S_b, src_b = x["S_w"], "S_w approximation"
        else:
            log.append("S_b unavailable (need S_b, sigma_b + T_cond_b, or S_w)")
    except ValueError as e:
        log.append(f"Method 0: {e}")
    rho_g = x["M_g"] / x["V_g"] if (x["M_g"] is not None and x["V_g"] is not None) else None
    if rho_g is None:
        log.append("rho_g unavailable (need M_g and V_g)")
    else:
        _warn_range("rho_g", rho_g, 900.0, 1100.0)

    # --- measured / constant / computed rows
    for n, unit, desc in INPUTS:
        if x[n] is not None:
            add(desc, n, x[n], unit, 'measured')
    for nm, nt, v, u in [('Pure ice density', 'rho_i', rho_i, 'kg m^-3'), ('Reference fresh water density', 'rho_0', rho_0, 'kg m^-3'),
                         ('Haline contraction coefficient', 'beta', beta, 'kg m^-3 / unit salinity')]:
        add(nm, nt, v, u, 'constant')
    if S_g is not None and src_g != "input":
        add('Grease bulk salinity (PSS-78 from conductivity)', 'S_g', S_g, 'psu', 'computed', source=src_g)
    if S_b is not None and src_b != "input":
        add('Brine salinity', 'S_b', S_b, 'psu', 'computed', source=src_b)
    if rho_g is not None:
        add('Grease density (M_g / V_g)', 'rho_g', rho_g, 'kg m^-3', 'computed')

    # --- Methods 1-4
    def air_for(m):
        if not air[m]:
            return 0.0, "off"
        if x["phi_a"] is not None:
            return x["phi_a"], "input"
        if phi_a_default is not None:
            return _mk(*phi_a_default), "default"
        return None, "missing"

    summary = {}
    spec = {1: dict(need={"S_g": S_g, "S_b": S_b}, airneed={"rho_g": rho_g}),
            2: dict(need={"S_b": S_b, "rho_g": rho_g}, airneed={}),
            3: dict(need={"S_g": S_g, "rho_g": rho_g}, airneed={}),
            4: dict(need={"S_b": S_b, "V_g": x["V_g"], "M'_g": x["M'_g"], "S'_g": x["S'_g"]}, airneed={})}
    nice = {'rho_b': ('Brine density', 'kg m^-3'), 'phi_s': ('Solid volume fraction', 'dimensionless'),
            'phi_b': ('Brine volume fraction', 'dimensionless'), 'S_b': ('Derived brine salinity', 'psu'),
            'S_g': ('Derived in situ bulk salinity', 'psu')}
    order = {1: ['rho_b', 'phi_s', 'phi_b'], 2: ['rho_b', 'phi_s', 'phi_b'], 3: ['S_b', 'rho_b', 'phi_s', 'phi_b'],
             4: ['rho_b', 'phi_s', 'S_g', 'phi_b']}
    for m in (1, 2, 3, 4):
        need = dict(spec[m]["need"])
        if air[m]:
            need.update(spec[m]["airneed"])
        phi_a_m, air_src = air_for(m)
        if phi_a_m is None:
            need["phi_a"] = None
        missing = [k for k, v in need.items() if v is None]
        if missing:
            log.append(f"Method {m} skipped: missing {', '.join(missing)}")
            continue
        try:
            if m == 1:
                res = method1(S_g, S_b, rho_g, phi_a_m, air[m])
            elif m == 2:
                res = method2(S_b, rho_g, phi_a_m, air[m])
            elif m == 3:
                res = method3(S_g, rho_g, phi_a_m, air[m])
            else:
                res = method4(S_b, x["V_g"], x["M'_g"], x["S'_g"], phi_a_m, air[m])
            _warn_range(f"Method {m} phi_s", res["phi_s"], 0.0, 1.0)
        except (ValueError, ZeroDivisionError) as e:
            log.append(f"Method {m} skipped: {e}")
            continue
        for k in order[m]:
            add(*nice[k][:1], k, res[k], nice[k][1], 'derived', m, flag=air[m], air_used=phi_a_m, air_src=air_src)
            summary[f"{k}_M{m}"], summary[f"{k}_M{m}_unc"] = _nv(res[k]), _sd(res[k])
        summary[f"air_M{m}"] = air_src
    summary.update(S_g=_nv(S_g) if S_g is not None else np.nan, S_g_unc=_sd(S_g) if S_g is not None else np.nan, S_g_source=src_g,
                   S_b=_nv(S_b) if S_b is not None else np.nan, S_b_unc=_sd(S_b) if S_b is not None else np.nan, S_b_source=src_b,
                   rho_g=_nv(rho_g) if rho_g is not None else np.nan, rho_g_unc=_sd(rho_g) if rho_g is not None else np.nan,
                   notes="; ".join(log))
    return records, summary, log


# %% Reading the raw CSV
KNOWN = {n for n, _, _ in INPUTS} | {n + "_unc" for n, _, _ in INPUTS}


def read_raw_csv(path):
    """Tolerant reader: auto-detects ; or , separator, decimal comma, BOM; NaN-like cells -> NaN."""
    for enc in ("utf-8-sig", "latin-1"):
        try:
            df = pd.read_csv(path, sep=None, engine="python", dtype=str, keep_default_na=False, encoding=enc)
            break
        except UnicodeDecodeError:
            continue
    df.columns = [c.strip() for c in df.columns]
    for c in df.columns:
        if c in KNOWN:
            df[c] = pd.to_numeric(df[c].astype(str).str.strip().str.replace(",", ".", regex=False), errors="coerce")
    ignored = [c for c in df.columns if c not in KNOWN and c not in ("date", "time")]
    absent = [n for n, _, _ in INPUTS if n not in df.columns]
    print(f"Read {len(df)} sample(s) from {Path(path).name}")
    if ignored:
        print(f"  columns ignored: {ignored}")
    if absent:
        print(f"  input columns absent (treated as not measured): {absent}")
    for c in ("date", "time"):
        if c not in df.columns:
            print(f"  ! no '{c}' column: file names will use the row number")
    return df


def make_stamp(row, i):
    d, t = str(row.get("date", "")).strip(), str(row.get("time", "")).strip()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ts = pd.to_datetime(f"{d} {t}".strip(), dayfirst=True)
        return ts.strftime("%Y%m%d_%H%M%S") if t else ts.strftime("%Y%m%d")
    except (ValueError, TypeError):
        return slug(f"{d}_{t}") if (d or t) else f"row{i + 1}"


# %% Interactive choice of folder / file
def _clean_path(s):
    return Path(s.strip().strip('"').strip("'")).expanduser()


def choose_raw_csv(ask=None):
    ask = ask or input
    p = _clean_path(raw_data_folder) if raw_data_folder else _clean_path(ask("Folder (or CSV file) with the raw grease-ice data [Enter = current folder]: ") or ".")
    if p.is_file():
        return p
    if not p.is_dir():
        raise FileNotFoundError(f"not found: {p}")
    files = [f for f in sorted(p.glob("*.csv")) if not f.name.startswith(("grease_ice_results_", "grease_ice_summary_"))]
    if not files:
        raise FileNotFoundError(f"no raw CSV in {p}")
    for i, f in enumerate(files):
        print(f"  [{i}] {f.name}")
    a = ask(f"Which file? index [Enter = 0]: ").strip() or "0"
    return files[int(a)]


# %% Main
def main(ask=None):
    ask = ask or input
    csv_path = choose_raw_csv(ask)
    df = read_raw_csv(csv_path)
    out_dir = Path(output_folder) if output_folder else csv_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    stamps = [make_stamp(r, i) for i, r in enumerate(df.to_dict("records"))]
    seen = {}
    for i, s in enumerate(stamps):                       # identical date+time -> make file names unique
        seen[s] = seen.get(s, 0) + 1
        if seen[s] > 1:
            stamps[i] = f"{s}_r{i + 1}"
            print(f"  ! duplicate date/time {s}: row {i + 1} saved as {stamps[i]}")

    written, summaries = [], []
    for i, (row, stamp) in enumerate(zip(df.to_dict("records"), stamps)):
        records, summary, log = analyse_sample(row)
        d, t = str(row.get("date", "")), str(row.get("time", ""))
        print(f"\n--- sample {i + 1}: {d} {t} ---")
        for m in (1, 2, 3, 4):
            if f"phi_s_M{m}" in summary:
                u = summary[f"phi_s_M{m}_unc"]
                print(f"  Method {m} (air: {summary[f'air_M{m}']}): phi_s = {summary[f'phi_s_M{m}']:.3f} +/- {u:.3f}")
        for msg in log:
            print(f"  - {msg}")
        df_s = pd.DataFrame(records)
        df_s.insert(0, "Time", t)
        df_s.insert(0, "Date", d)
        if saved_method != 'all':
            df_s = df_s[(df_s['Label'] != 'derived') | (df_s['Method'] == int(saved_method))]
        fn = out_dir / f"grease_ice_results_{stamp}.csv"
        df_s.to_csv(fn, index=False)
        written.append(fn)
        summaries.append({"date": d, "time": t, "file": fn.name, **summary})

    summ = out_dir / f"grease_ice_summary_{csv_path.stem}.csv"
    pd.DataFrame(summaries).to_csv(summ, index=False)
    print(f"\n{len(written)} results file(s) + summary written to {out_dir}\n  summary: {summ.name}")

    if link_to_metadata and cio is not None and written:
        meta = cio.choose_meta_file(meta_root, ask=ask)
        if meta:
            for fn in written:
                cio.register_file(meta, "sample", fn.name)
    elif link_to_metadata and cio is None:
        print("chipiron_io not found: results not linked to any metadata file")
    return written, summ


if __name__ == "__main__":
    main()