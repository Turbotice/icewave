"""
Solid fraction and salinity of grease ice: Methods 1 to 4 of
Maus & De la Rosa (2012), 'Salinity and solid fraction of frazil and grease ice', J. Glaciol. 58(209), 594-612,
plus Method 0 (PSS-78: conductivity + temperature -> salinity).

WHAT IT DOES
    Asks for a folder, lets you pick a raw CSV in it, then
      1. checks which columns were recognised (and suggests the right name for near-misses); completely empty
         rows (e.g. ',,,,,') are dropped;
      2. GROUPS the rows into SAMPLES: rows with the same date + time (and, as a check, the same station) are the
         drained / undrained buckets of ONE sample and are merged into one sample;
      3. runs every method the available data allow, TWICE per sample:
            mode "direct"       salinities S_g, S_b, S'_g, S_w as given (psu) in the input file
            mode "conductivity" salinities from PSS-78 (sigma_g, sigma_b, sigma'_g, sigma_w; see conductivity_ref_T)
         (the two modes are never mixed: a method that lacks a salinity in one mode is skipped in that mode);
      4. writes one results CSV per sample AND per mode (+ one summary table), each with a 'Salinity_Mode' and a
         'Comment' column saying which salinities were used, and prints the same comment on screen;
      5. optionally links the results files to an experiment's *_meta.txt.

EXPECTED INPUT CSV (the "template")
    * one row per BUCKET; header names are case-sensitive and use the paper's notation (underscore for subscript,
      "'" for primes). Separator , ; or TAB (auto-detected from the header) ; decimal point, or decimal comma with ';'.
    * for every quantity X there is an OPTIONAL uncertainty column  X_unc  (1 standard deviation, same unit as X).
      Empty/NaN X_unc = no uncertainty (it then contributes 0 to the propagated uncertainty).
    * empty cells, NaN, NA, n/a, '-' are all read as "not measured".
    * 'note' (free text) and optional 'station' columns: they decide which bucket is drained / undrained and which
      rows belong together (see SAMPLES below).

      column            unit       meaning
      ----------------  ---------  ----------------------------------------------------------------------------
      date, time        text       sampling date (day first: 3/11/2026) and time; SAME date+time = same sample
      T_cond_g          degC       REAL temperature of the MELTED (undrained) grease when sigma_g was measured
      sigma_g           mS/cm      electrical conductivity of the melted (undrained) grease
      T_cond_b          degC       real temperature of the brine when sigma_b was measured
      sigma_b           mS/cm      electrical conductivity of the brine
      T'_cond_g         degC       real temperature of the melted DRAINED grease when sigma'_g was measured
      sigma'_g          mS/cm      electrical conductivity of the melted DRAINED grease
      T_cond_w          degC       real temperature of the sea water when sigma_w was measured
      sigma_w           mS/cm      electrical conductivity of the sea water
      S_g, S_b, S'_g    -          salinities given directly (used in mode "direct")
      S_w               -          sea-water salinity (direct input, mode "direct"); in mode "conductivity" it is
                                   computed from sigma_w. Used as S_b ~ S_w ONLY if S_b is not available in that mode
      V_g               m^3        grease volume BEFORE drainage (one value per bucket)
      M_g               kg         grease mass before drainage (one value per bucket)
      M'_g              kg         mass of the drained grease ice (drained bucket)
      phi_a             -          air volume fraction (used only where include_air is True)
      note, station     text       'undrained sample - station 101' / 'drained sample - station 101'

CONDUCTIVITY AND TEMPERATURE COMPENSATION
    If the conductivity meter is temperature-compensated (reports the conductivity the liquid WOULD have at a
    reference temperature, typically 25 degC), PSS-78 must be given that reference temperature, NOT the real
    temperature of the sample. Set  conductivity_ref_T = 25.0  (default here). The T_cond_* columns are then only
    recorded in the output (they are not needed for the salinity). Set  conductivity_ref_T = None  if the meter
    reports the raw conductivity at the real temperature: PSS-78 then uses the T_cond_* columns.
    Note: a meter's built-in compensation (often a linear ~2 %/degC) is not identical to PSS-78, so for samples far
    from 25 degC (cold brine) the compensated reading carries an extra, unquantified systematic error.

SAMPLES (drained / undrained buckets)
    A row is "drained" if its note contains 'drained' (but not 'undrained'), or, without a note, if it has M'_g,
    S'_g, sigma'_g or T'_cond_g. Otherwise it is "undrained". Rows with the same date + time form one sample.
    The station (column 'station', or 'station <id>' in the note) is checked: two rows with the same date+time but
    different stations are NOT merged; the same station on the same date at different times gives a warning.
    Because the two buckets hold different volumes/masses:
        Methods 1-3 use rho_g = M_g / V_g of the UNDRAINED bucket (rho_from = "undrained"), and S_g of that bucket;
        Method 4 uses V_g of the DRAINED bucket (volume before drainage), M'_g and S'_g;
        S_b (drained brine) is shared by all methods.
    The two bucket densities are compared and a warning is printed if they differ by more than rho_mismatch_tol.

    What each method needs (anything missing -> that method is skipped for that sample/mode):
      Method 1: S_g, S_b                     (+ M_g, V_g, phi_a if its air correction is on)
      Method 2: S_b, M_g, V_g                (+ phi_a if air on)
      Method 3: S_g, M_g, V_g                (+ phi_a if air on)
      Method 4: S_b, V_g (drained bucket), M'_g, S'_g   (+ phi_a if air on)
    If include_air is on for a method but phi_a is empty, phi_a_default is used (0.01 +/- 0.005: the paper's
    'typical air porosity', from Nakawo 1983: 0.005 < phi_a < 0.015) and flagged in the output.

OUTPUT
    grease_ice_results_<stamp>_<mode>.csv   one per sample and mode, long format (read by chipiron_io.read_sample_csv):
        Date | Time | Station | Salinity_Mode | Comment | Parameter Name | Notation | Value | Uncertainty | Units |
        Label | Method | Source | Saved_Method_Flag | Air_Included | Air_Fraction_Used | Air_Source
        Label: 'measured' (your inputs), 'constant', 'computed' (S_w, S_g, S_b, S'_g, rho_g), 'derived' (Methods 1-4).
        Uncertainty 0 = no input uncertainty was given.
    grease_ice_summary_<input file name>.csv    one wide row per sample AND mode, with the comment and the reason
                                                for every skipped method.
    Files go to output_folder (set it below), or next to the input CSV if output_folder is None.
"""

# %% Initialization
import difflib
import io
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
output_folder = r"C:\Users\Gulliver\Desktop\Chipiron\git_clone_turbotice_icewave\icewave\icewave\chipiron\frazil_characterization"       # e.g. r"C:\Users\me\PMMH\results_grease" ; None = same folder as the input CSV
saved_method = 'all'        # 1, 2, 3, 4 or 'all': which methods go in the per-sample results CSV
link_to_metadata = True     # offer to register the results files in an experiment *_meta.txt
meta_root = "PhD_data"      # where chipiron_io looks for *_meta.txt

salinity_modes = ("direct", "conductivity")   # every sample is computed once per mode (order = order of output)
rho_from = "undrained"      # bucket whose M_g / V_g gives rho_g for Methods 1-3: "undrained" or "drained"
rho_mismatch_tol = 0.05     # warn if the two buckets' densities differ by more than this (relative)

# Conductivity meter temperature compensation (see docstring):
#   25.0 -> the meter is compensated to 25 degC: PSS-78 is evaluated at 25 degC; T_cond_* are only recorded
#   None -> raw conductivity at the real temperature: PSS-78 uses the T_cond_* columns
conductivity_ref_T = 25.0

# Air correction per method (True = air-corrected equation, False = air-free equation), or one bool for all.
# NOTE Method 4: by our own derivation Eq 24 is already exact when V_g is the total sample volume, so the paper's
# correction is implemented as published but defaults to False. Set True to follow the paper literally.
include_air = {1: True, 2: True, 3: True, 4: False}
# Paper: "we may consider phi_a = 0.01 as a typical air porosity" (Nakawo 1983: 0.005-0.015) -> +/- 0.005
phi_a_default = (0.01, 0.005)   # (value, uncertainty) used when a sample has no phi_a but air is switched on; None = skip instead

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
    ("T'_cond_g", "°C", "Temperature of melted drained grease at conductivity measurement"),
    ("sigma'_g", "mS cm^-1", "Electrical conductivity of melted drained grease"),
    ("T_cond_w", "°C", "Temperature of sea water at conductivity measurement"),
    ("sigma_w", "mS cm^-1", "Electrical conductivity of sea water"),
    ("S_g", "psu", "Grease bulk salinity (direct input)"),
    ("S_b", "psu", "Brine salinity (direct input)"),
    ("S_w", "psu", "Sea-water salinity (direct input; S_b ~ S_w approximation)"),
    ("V_g", "m^3", "Grease sample volume before drainage"),
    ("M_g", "kg", "Grease sample mass before drainage"),
    ("M'_g", "kg", "Drained grease ice mass"),
    ("S'_g", "psu", "Drained grease ice salinity (direct input)"),
    ("phi_a", "dimensionless", "Air volume fraction"),
]
# values that only exist after the two buckets of a sample are merged (V_g / M_g of the drained bucket)
EXTRA = [
    ("V_g_drained", "m^3", "Grease sample volume before drainage (drained bucket)"),
    ("M_g_drained", "kg", "Grease sample mass before drainage (drained bucket)"),
]
POSITIVE = ["sigma_g", "sigma_b", "sigma_w", "sigma'_g", "S_g", "S_b", "S_w", "V_g", "M_g", "M'_g", "S'_g",
            "V_g_drained", "M_g_drained"]

STATION_RE = re.compile(r"station\s*[:#n°.]*\s*([A-Za-z0-9_-]+)", re.I)


def mode_text(mode):
    if mode == "direct":
        return "Salinities used: DIRECT (S_g, S_b, S'_g, S_w as given in psu in the input file)"
    if conductivity_ref_T is not None:
        return (f"Salinities used: CONDUCTIVITY (PSS-78 from sigma at {conductivity_ref_T:g} degC, "
                f"meter temperature-compensated: S_g, S_b, S'_g, S_w)")
    return "Salinities used: CONDUCTIVITY (PSS-78 from sigma and real T_cond: S_g, S_b, S'_g, S_w)"


# %% Helpers
def _nv(x):
    return getattr(x, 'nominal_value', x)


def _sd(x):
    return getattr(x, 'std_dev', 0.0)


def _mk(v, u):
    """ufloat if an uncertainty is given, else a plain float (avoids std_dev == 0 ufloats)."""
    return ufloat(v, u) if (u is not None and np.isfinite(u) and u > 0) else float(v)


def _f(x):
    """float or NaN."""
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def _s(x):
    """clean text ('' for NaN / None)."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return ""
    s = str(x).strip()
    return "" if s.lower() == "nan" else s


def _val(row, name):
    """Value (+ optional '<name>_unc') from a row dict -> float | ufloat | None (None = missing/NaN)."""
    v = _f(row.get(name, np.nan))
    if not np.isfinite(v):
        return None
    return _mk(v, _f(row.get(name + "_unc", np.nan)))


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


# %% Analysis of ONE sample (a merged dict of raw values) in ONE salinity mode
def resolve_salinities(x, mode, log):
    """S_w, S_g, S_b, S'_g for this mode only ('direct' = the psu columns, 'conductivity' = PSS-78). Never mixes modes.
    In conductivity mode PSS-78 is evaluated at conductivity_ref_T (temperature-compensated meter) or, if that is
    None, at the real T_cond_* of each measurement.
    Returns ({name: value or None}, {name: source text})."""
    out, src = {}, {}
    tag = f"[{mode}]"
    for key, sig, T in (("S_w", "sigma_w", "T_cond_w"), ("S_g", "sigma_g", "T_cond_g"),
                        ("S_b", "sigma_b", "T_cond_b"), ("S'_g", "sigma'_g", "T'_cond_g")):
        val, s = None, ""
        try:
            if mode == "direct":
                if x[key] is not None:
                    val, s = x[key], "input"
            elif x[sig] is not None:
                T_use = conductivity_ref_T if conductivity_ref_T is not None else x[T]
                if T_use is not None:
                    if conductivity_ref_T is None:
                        _warn_range(T, x[T], -2.0, 35.0)
                    val = pss78_salinity(x[sig], T_use)
                    s = "PSS-78" if conductivity_ref_T is None else f"PSS-78 @ {conductivity_ref_T:g} degC"
                    _warn_range(f"{key} (PSS-78)", val, 2.0, 42.0)
        except ValueError as e:
            log.append(f"{tag} {key}: {e}")
        if val is None and key == "S_b" and out["S_w"] is not None:
            val, s = out["S_w"], "S_w approximation"
        if val is None and key != "S_w":
            if mode == "direct":
                need = key
            else:
                need = sig if conductivity_ref_T is not None else f"{sig} + {T}"
            log.append(f"{tag} {key} unavailable (need {need}" + (" or S_w" if key == "S_b" else "") + ")")
        out[key], src[key] = val, s
    return out, src


def analyse_sample(row, mode):
    """-> (records, summary, log). records = list of dict rows for the long-format CSV."""
    air = _air_flags(include_air)
    log, records = [], []
    x = {n: _val(row, n) for n, _, _ in INPUTS + EXTRA}

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

    # --- salinities of this mode, densities
    sal, src = resolve_salinities(x, mode, log)
    S_g, S_b, S_pg = sal["S_g"], sal["S_b"], sal["S'_g"]

    def dens(M, V):
        return M / V if (M is not None and V is not None) else None
    rho_U, rho_D = dens(x["M_g"], x["V_g"]), dens(x["M_g_drained"], x["V_g_drained"])
    for lab, r_ in (("undrained bucket", rho_U), ("drained bucket", rho_D)):
        if r_ is not None and not (900.0 <= _nv(r_) <= 1100.0):
            log.append(f"rho_g ({lab}) = {_nv(r_):.0f} kg/m3 is outside 900-1100: check M_g / V_g of that bucket")
    if rho_U is not None and rho_D is not None and abs(_nv(rho_U) - _nv(rho_D)) > rho_mismatch_tol * _nv(rho_U):
        log.append(f"rho_g differs between buckets: undrained {_nv(rho_U):.0f}, drained {_nv(rho_D):.0f} kg/m3 "
                   f"(> {100 * rho_mismatch_tol:.0f} %); Methods 1-3 use the {rho_from} bucket")
    rho_g = (rho_U, rho_D) if rho_from == "undrained" else (rho_D, rho_U)
    rho_g = rho_g[0] if rho_g[0] is not None else rho_g[1]
    if rho_g is None:
        log.append("rho_g unavailable (need M_g and V_g)")

    # --- measured / constant / computed rows
    for n, unit, desc in INPUTS + EXTRA:
        if x[n] is not None:
            add(desc, n, x[n], unit, 'measured')
    for nm, nt, v, u in [('Pure ice density', 'rho_i', rho_i, 'kg m^-3'), ('Reference fresh water density', 'rho_0', rho_0, 'kg m^-3'),
                         ('Haline contraction coefficient', 'beta', beta, 'kg m^-3 / unit salinity')]:
        add(nm, nt, v, u, 'constant')
    for key, desc in (("S_w", "Sea-water salinity"), ("S_g", "Grease bulk salinity"), ("S_b", "Brine salinity"),
                      ("S'_g", "Drained grease ice salinity")):
        if sal[key] is not None and not (key == "S_w" and src[key] == "input"):   # direct S_w is already a 'measured' row
            add(f"{desc} [{mode} mode]", key, sal[key], 'psu', 'computed' if src[key] != "input" else 'measured',
                source=src[key])
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
            4: dict(need={"S_b": S_b, "V_g (drained bucket)": x["V_g_drained"], "M'_g": x["M'_g"], "S'_g": S_pg}, airneed={})}
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
            log.append(f"[{mode}] Method {m} skipped: missing {', '.join(missing)}")
            continue
        try:
            if m == 1:
                res = method1(S_g, S_b, rho_g, phi_a_m, air[m])
            elif m == 2:
                res = method2(S_b, rho_g, phi_a_m, air[m])
            elif m == 3:
                res = method3(S_g, rho_g, phi_a_m, air[m])
            else:
                res = method4(S_b, x["V_g_drained"], x["M'_g"], S_pg, phi_a_m, air[m])
            _warn_range(f"Method {m} phi_s", res["phi_s"], 0.0, 1.0)
        except (ValueError, ZeroDivisionError) as e:
            log.append(f"[{mode}] Method {m} skipped: {e}")
            continue
        for k in order[m]:
            add(*nice[k][:1], k, res[k], nice[k][1], 'derived', m, flag=air[m], air_used=phi_a_m, air_src=air_src)
            summary[f"{k}_M{m}"], summary[f"{k}_M{m}_unc"] = _nv(res[k]), _sd(res[k])
        summary[f"air_M{m}"] = air_src

    # --- comment: which salinities / air values were used
    comment = mode_text(mode) + " [" + ", ".join(f"{k}: {v}" for k, v in src.items() if v) + "]"
    dflt = [m for m in (1, 2, 3, 4) if summary.get(f"air_M{m}") == "default"]
    if dflt:
        comment += (f"; phi_a not measured -> paper default {phi_a_default[0]:g} +/- {phi_a_default[1]:g} "
                    f"used in Method(s) {', '.join(map(str, dflt))}")

    def _n(v):
        return _nv(v) if v is not None else np.nan

    def _u(v):
        return _sd(v) if v is not None else np.nan

    summary.update(salinity_mode=mode, comment=comment,
                   S_w=_n(sal["S_w"]), S_w_unc=_u(sal["S_w"]), S_w_source=src["S_w"],
                   S_g=_n(S_g), S_g_unc=_u(S_g), S_g_source=src["S_g"],
                   S_b=_n(S_b), S_b_unc=_u(S_b), S_b_source=src["S_b"],
                   S_prime_g=_n(S_pg), S_prime_g_unc=_u(S_pg), S_prime_g_source=src["S'_g"],
                   rho_g=_n(rho_g), rho_g_unc=_u(rho_g),
                   notes="; ".join(log))
    return records, summary, log


# %% Reading the raw CSV
KNOWN = {n for n, _, _ in INPUTS} | {n + "_unc" for n, _, _ in INPUTS}
TEXT_COLS = ("date", "time", "note", "station")


def read_raw_csv(path):
    """Tolerant reader: separator (, ; or TAB) taken from the header line, decimal comma, BOM; NaN-like cells -> NaN;
    completely empty rows (only separators) are dropped.
    Prints which columns were recognised, which have an uncertainty column, which were ignored (with a hint if
    the name looks like a typo of a known one) and which inputs are absent."""
    text = None
    for enc in ("utf-8-sig", "latin-1"):
        try:
            text = Path(path).read_text(encoding=enc)
            break
        except UnicodeDecodeError:
            continue
    lines = [ln for ln in text.splitlines() if ln.strip()]
    sep = max((",", ";", "\t"), key=lines[0].count)
    df = pd.read_csv(io.StringIO("\n".join(lines)), sep=sep, dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    blank = df.apply(lambda r: not any(str(v).strip() for v in r), axis=1)   # rows like ',,,,,,'
    if blank.any():
        print(f"  {int(blank.sum())} completely empty row(s) dropped")
    df = df[~blank].reset_index(drop=True)
    for c in df.columns:
        if c in KNOWN:
            df[c] = pd.to_numeric(df[c].astype(str).str.strip().str.replace(",", ".", regex=False), errors="coerce")
    names = [n for n, _, _ in INPUTS]
    recognised = [n for n in names if n in df.columns]
    with_unc = [n for n in recognised if n + "_unc" in df.columns]
    ignored = [c for c in df.columns if c not in KNOWN and c not in TEXT_COLS]
    absent = [n for n in names if n not in df.columns]
    print(f"Read {len(df)} row(s) from {Path(path).name} (separator {'TAB' if sep == chr(9) else repr(sep)})")
    print(f"  recognised inputs ({len(recognised)}): {recognised}")
    print(f"  with an _unc column: {with_unc}; _unc without its value column: "
          f"{[c for c in df.columns if c.endswith('_unc') and c[:-4] not in names]}")
    print(f"  text columns found: {[c for c in TEXT_COLS if c in df.columns]}")
    if ignored:
        print(f"  columns ignored: {ignored}")
        for c in ignored:
            near = difflib.get_close_matches(c, sorted(KNOWN), n=1, cutoff=0.7)
            if near:
                print(f"    ! '{c}' looks like '{near[0]}' (typo?)")
    if absent:
        print(f"  input columns absent (treated as not measured): {absent}")
    for c in ("date", "time"):
        if c not in df.columns:
            print(f"  ! no '{c}' column: rows cannot be grouped into samples by date/time")
    return df


# %% Grouping the rows into samples (drained / undrained buckets)
def _norm_dt(row):
    """(YYYY-MM-DD, HH:MM) from the raw date/time text (day first); falls back to the raw text."""
    d, t = _s(row.get("date")), _s(row.get("time"))
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ts = pd.to_datetime(f"{d} {t}".strip(), dayfirst=True)
        return ts.strftime("%Y-%m-%d"), (ts.strftime("%H:%M") if t else "")
    except (ValueError, TypeError):
        return d, t


def _row_station(row):
    s = _s(row.get("station"))
    if s:
        return s
    m = STATION_RE.search(_s(row.get("note")))
    return m.group(1) if m else ""


def _row_role(row):
    note = _s(row.get("note")).lower()
    if "undrained" in note or "not drained" in note:
        return "undrained"
    if "drained" in note:
        return "drained"
    has = any(np.isfinite(_f(row.get(n))) for n in ("M'_g", "S'_g", "sigma'_g", "T'_cond_g"))
    return "drained" if has else "undrained"


def _first(rows, name):
    for r in rows:
        v = _f(r.get(name))
        if np.isfinite(v):
            return v, _f(r.get(name + "_unc"))
    return None


def merge_group(rs):
    """Merge the bucket rows of ONE sample -> (merged dict, log)."""
    order = sorted(rs, key=lambda r: (r["_role"] != "undrained", r["_i"]))
    U = [r for r in order if r["_role"] == "undrained"]
    D = [r for r in order if r["_role"] == "drained"]
    log = []
    for role, lst in (("undrained", U), ("drained", D)):
        if len(lst) > 1:
            log.append(f"{len(lst)} {role} rows for the same sample: values merged (first non-empty kept)")
    if not U:
        log.append("no undrained row: S_g / rho_g come from the drained bucket if present")
    if not D:
        log.append("no drained row: Method 4 and S_b (drained brine) may be missing")
    m = {"date": order[0]["_date"], "time": order[0]["_time"],
         "station": next((r["_station"] for r in order if r["_station"]), ""),
         "note": " | ".join(_s(r.get("note")) for r in order if _s(r.get("note"))),
         "roles": ", ".join(f"row {r['_i']} ({r['_role']})" for r in order)}
    per_bucket = ("V_g", "M_g")
    for n, _, _ in INPUTS:
        if n in per_bucket:
            continue
        found = [f for f in (_first([r], n) for r in order) if f]
        if not found:
            continue
        m[n], m[n + "_unc"] = found[0]
        other = [f[0] for f in found[1:] if not np.isclose(f[0], found[0][0], rtol=1e-6)]
        if other:
            log.append(f"{n} given in several rows with different values: {found[0][0]:g} kept (also {other})")
    for n in per_bucket:                                  # V_g, M_g: one value per bucket
        fu, fd = _first(U, n), _first(D, n)
        if fu:
            m[n], m[n + "_unc"] = fu
        elif fd:
            m[n], m[n + "_unc"] = fd
            log.append(f"{n}: no undrained-bucket value, the drained bucket's is used for Methods 1-3")
        if n == "V_g":
            if fd:
                m["V_g_drained"], m["V_g_drained_unc"] = fd
            elif fu and D:
                m["V_g_drained"], m["V_g_drained_unc"] = fu
                log.append("V_g_drained: no drained-bucket value, the undrained bucket's is used for Method 4")
        elif fd:
            m["M_g_drained"], m["M_g_drained_unc"] = fd
    return m, log


def group_samples(df):
    """Rows -> samples. Same date+time = same sample (station checked). Returns [(merged dict, log), ...]."""
    rows = df.to_dict("records")
    for i, r in enumerate(rows, start=1):
        r["_i"] = i
        r["_date"], r["_time"] = _norm_dt(r)
        r["_station"] = _row_station(r)
        r["_role"] = _row_role(r)
    buckets = {}
    for r in rows:
        buckets.setdefault((r["_date"], r["_time"]), []).append(r)
    groups, extra = [], {}
    for key, rs in buckets.items():
        stations = sorted({r["_station"] for r in rs if r["_station"]})
        if len(stations) > 1:
            extra[key] = f"different stations {stations} at the same date/time: rows NOT merged across stations"
            for st in stations + [""]:
                sub = [r for r in rs if r["_station"] == st]
                if sub:
                    groups.append(sub)
        else:
            groups.append(rs)
    out = []
    for g in groups:
        merged, log = merge_group(g)
        key = (g[0]["_date"], g[0]["_time"])
        if key in extra:
            log.insert(0, extra[key])
        out.append((merged, log))
    # same station on the same date at different times: possible split sample?
    seen = {}
    for merged, log in out:
        if merged["station"]:
            seen.setdefault((merged["date"], merged["station"]), []).append(merged)
    for (d, st), lst in seen.items():
        if len(lst) > 1:
            msg = (f"station {st} on {d} appears at several times ({', '.join(m['time'] for m in lst)}): "
                   "different samples? (not merged: only identical date+time are merged)")
            for m in lst:
                next(lg for mm, lg in out if mm is m).append(msg)
    return out


def make_stamp(merged, i):
    d, t = merged.get("date", ""), merged.get("time", "")
    base = (d.replace("-", "") + (("_" + t.replace(":", "")) if t else "")) if d else f"row{i + 1}"
    return base + (f"_st{slug(merged['station'])}" if merged.get("station") else "")


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
def main(ask=None, csv_path=None):
    ask = ask or input
    csv_path = Path(csv_path) if csv_path else choose_raw_csv(ask)
    df = read_raw_csv(csv_path)
    out_dir = _clean_path(output_folder) if output_folder else csv_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    samples = group_samples(df)
    print(f"\n{len(df)} row(s) -> {len(samples)} sample(s):")
    for k, (m, glog) in enumerate(samples, start=1):
        print(f"  sample {k}: {m['date']} {m['time']}  station {m['station'] or '-'}  <- {m['roles']}")
        for msg in glog:
            print(f"      - {msg}")

    stamps = [make_stamp(m, i) for i, (m, _) in enumerate(samples)]
    seen = {}
    for i, s in enumerate(stamps):                       # identical stamps -> make file names unique
        seen[s] = seen.get(s, 0) + 1
        if seen[s] > 1:
            stamps[i] = f"{s}_s{i + 1}"
            print(f"  ! duplicate stamp {s}: sample {i + 1} saved as {stamps[i]}")

    written, summaries, compare = [], [], []
    for k, ((m, glog), stamp) in enumerate(zip(samples, stamps), start=1):
        print(f"\n=== sample {k}: {m['date']} {m['time']}  station {m['station'] or '-'} ===")
        printed = set()
        for mode in salinity_modes:
            records, summary, log = analyse_sample(m, mode)
            print(f"\n  [{mode.upper()}] {summary['comment']}")
            print(f"    S_w = {summary['S_w']:.2f} ({summary['S_w_source'] or 'n/a'}), "
                  f"S_g = {summary['S_g']:.2f} ({summary['S_g_source'] or 'n/a'}), S_b = {summary['S_b']:.2f} "
                  f"({summary['S_b_source'] or 'n/a'}), S'_g = {summary['S_prime_g']:.2f} ({summary['S_prime_g_source'] or 'n/a'}), "
                  f"rho_g = {summary['rho_g']:.0f} kg/m3")
            for meth in (1, 2, 3, 4):
                if f"phi_s_M{meth}" in summary:
                    u = summary[f"phi_s_M{meth}_unc"]
                    print(f"    Method {meth} (air: {summary[f'air_M{meth}']}): phi_s = {summary[f'phi_s_M{meth}']:.3f} +/- {u:.3f}"
                          f"   [{mode}]")
                    compare.append((k, m['date'], m['time'], meth, mode, summary[f'phi_s_M{meth}'], u))
            for msg in log:
                if msg not in printed:                    # sample-level warnings (same in both modes) shown once
                    print(f"    - {msg}")
                    printed.add(msg)
            df_s = pd.DataFrame(records)
            for col, val in (("Comment", summary["comment"]), ("Salinity_Mode", mode), ("Station", m["station"]),
                             ("Time", m["time"]), ("Date", m["date"])):
                df_s.insert(0, col, val)
            if saved_method != 'all':
                df_s = df_s[(df_s['Label'] != 'derived') | (df_s['Method'] == int(saved_method))]
            fn = out_dir / f"grease_ice_results_{stamp}_{mode}.csv"
            df_s.to_csv(fn, index=False)
            written.append(fn)
            summaries.append({"date": m["date"], "time": m["time"], "station": m["station"], "rows_used": m["roles"],
                              "file": fn.name, **summary})

    summ = out_dir / f"grease_ice_summary_{csv_path.stem}.csv"
    pd.DataFrame(summaries).to_csv(summ, index=False)

    if compare and len(salinity_modes) > 1:               # side-by-side overview of the two modes
        print("\n--- phi_s by method: " + " vs ".join(salinity_modes) + " ---")
        tab = pd.DataFrame(compare, columns=["sample", "date", "time", "method", "mode", "phi_s", "unc"])
        for (k, d, t), g in tab.groupby(["sample", "date", "time"]):
            print(f"  sample {k} ({d} {t})")
            for meth, gg in g.groupby("method"):
                vals = {r.mode: f"{r.phi_s:.3f} +/- {r.unc:.3f}" for r in gg.itertuples()}
                print(f"    Method {meth}: " + "   ".join(f"{md}: {vals.get(md, 'skipped')}" for md in salinity_modes))
    print(f"\n{len(written)} results file(s) ({len(samples)} sample(s) x {len(salinity_modes)} mode(s)) + summary written to {out_dir}"
          f"\n  summary: {summ.name}")

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