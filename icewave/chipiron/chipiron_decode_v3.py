"""
chipiron_decode_v3.py -- read / decode the "Usb" csv logged by Phonefleet from the Arduino.

DEPENDENCIES   numpy, pandas.  No other Chipiron script (escon_edc_v3 can then convert volts to rpm / A).
USED BY        chipiron_plotting_v3.py

EXAMPLES
    from chipiron_decode_v3 import read_raw_usb
    d = read_raw_usb("P101_D2026-10-15T08_45_19-Usb-1-....csv")     # decoded, one row per 250 ms cycle
    raw = read_raw_usb(path, decode=False)                          # rows exactly as logged
    d[["t_s", "A1_adc", "A1_V"]].head()
    # physical units with the ESCON file:
    from escon_edc_v3 import load_edc;  load_edc("parameters_Sep26.edc").apply_to_frame(d)   # A1_speed_rpm, A2_current_A
    # other Arduino reference / resolution:  read_raw_usb(path, vref=3.3, adc_bits=12)

What the file contains (csv, no header): timestamp in microseconds (phone clock, only differences matter), one
byte per row (signed int8 as logged, e.g. -43 = 213), and an empty last column. Unused rows of the phone buffer
are placeholders (timestamp -1) and are dropped.

What the Arduino sends (PWM_full_v3.ino), every 250 ms:
    b1 = ADC(A1) & 0xFF   (t = 0 ms)      b2 = ADC(A1) >> 4   (t ~ 10 ms)
    b3 = ADC(A2) & 0xFF   (t ~ 50 ms)     b4 = ADC(A2) >> 4   (t ~ 60 ms)
    decode: adc = (b2 << 4) | (b1 & 0x0F)         check: (b1 >> 4) == (b2 & 0x0F)
A1 = speed output, A2 = current output of the ESCON in the saved .edc (CHECK the wiring).
The cycle start is found from the time gaps (> 100 ms between cycles, <= 40 ms inside a cycle), NOT from a fixed
"every 4th row", so a dropped byte only loses its own cycle.

Cycles are classified as
    "full"  : 4 rows -> ADC counts rebuilt (A1_adc, A2_adc), NaN where the consistency check fails
    "low8"  : 2 rows ~50 ms apart -> only b1 and b3 were logged (b2/b4 lost by the phone):
              only the LOW 8 BITS of the ADC are known (A1_lo, A2_lo), A*_adc stay NaN
    anything else is dropped and counted.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

CYCLE_GAP_US = 100_000          # gap above this = new telemetry cycle
LOW8_DT_US = (30_000, 80_000)   # spacing of b1 -> b3 when only two bytes per cycle were logged
ADC_VREF = 5.0                  # volts at full scale (Arduino Mega default)
ADC_BITS = 10


def read_usb_rows(path, verbose=True):
    """Raw rows as logged. Columns: t_us, usb_value (signed as logged), byte (0..255), t_s."""
    df = pd.read_csv(path, header=None, skipinitialspace=True)
    df = df.iloc[:, :2]
    df.columns = ["t_us", "usb_value"]
    df = df.apply(pd.to_numeric, errors="coerce").dropna()
    n_raw = len(df)
    df = df[df["t_us"] > 0].astype("int64").reset_index(drop=True)
    df["byte"] = df["usb_value"] & 0xFF
    if len(df):
        df["t_s"] = (df["t_us"] - df["t_us"].iloc[0]) / 1e6
    if verbose:
        print(f"[usb] {path}: {n_raw} rows read, {n_raw - len(df)} placeholder row(s) removed, {len(df)} kept")
    return df


def decode_usb(raw, vref=ADC_VREF, adc_bits=ADC_BITS, gap_us=CYCLE_GAP_US, verbose=True):
    """Decode raw rows (from read_usb_rows) into one row per telemetry cycle."""
    t = raw["t_us"].to_numpy(np.int64)
    b = raw["byte"].to_numpy(np.int64)
    if len(t) < 2:
        raise ValueError("not enough Usb rows to decode")
    starts = np.r_[0, np.where(np.diff(t) > gap_us)[0] + 1]
    counts = np.diff(np.r_[starts, len(t)])

    cols = {k: [] for k in ("t_us", "kind", "A1_lo", "A2_lo", "A1_adc", "A2_adc")}

    # --- full cycles: 4 rows
    s = starts[counts == 4]
    b1, b2, b3, b4 = b[s], b[s + 1], b[s + 2], b[s + 3]
    hi_bits = adc_bits - 4                                   # b2 = adc >> 4 must fit in adc_bits-4 bits
    ok1 = (((b1 >> 4) & 0xF) == (b2 & 0xF)) & ((b2 >> hi_bits) == 0)
    ok2 = (((b3 >> 4) & 0xF) == (b4 & 0xF)) & ((b4 >> hi_bits) == 0)
    a1 = np.where(ok1, (b2 << 4) | (b1 & 0xF), np.nan)
    a2 = np.where(ok2, (b4 << 4) | (b3 & 0xF), np.nan)
    for k, v in zip(cols, (t[s], np.full(len(s), "full"), b1, b3, a1, a2)):
        cols[k].append(v)

    # --- low-byte-only cycles: 2 rows ~50 ms apart (b1, b3)
    s2 = starts[counts == 2]
    dt = t[s2 + 1] - t[s2]
    s2 = s2[(dt >= LOW8_DT_US[0]) & (dt <= LOW8_DT_US[1])]
    nan = np.full(len(s2), np.nan)
    for k, v in zip(cols, (t[s2], np.full(len(s2), "low8"), b[s2], b[s2 + 1], nan, nan)):
        cols[k].append(v)

    out = pd.DataFrame({k: np.concatenate(v) for k, v in cols.items()}).sort_values("t_us").reset_index(drop=True)
    out["t_us"] = out["t_us"].astype("int64")
    out["t_s"] = (out["t_us"] - t[0]) / 1e6
    out["A1_V"] = out["A1_adc"] * vref / (2 ** adc_bits - 1)
    out["A2_V"] = out["A2_adc"] * vref / (2 ** adc_bits - 1)

    if verbose:
        n_full, n_low8 = int((out.kind == "full").sum()), int((out.kind == "low8").sum())
        n_drop = len(starts) - n_full - n_low8
        print(f"[usb] {len(starts)} cycle(s): {n_full} full, {n_low8} low-byte-only, {n_drop} dropped "
              f"(row counts per cycle: {dict(zip(*np.unique(counts, return_counts=True)))})")
        if n_full:
            print(f"[usb] consistency check failed (A -> NaN): A1 {int((~ok1).sum())}, A2 {int((~ok2).sum())} of {n_full}")
            print(f"[usb] A1_adc {np.nanmin(out.A1_adc):.0f}..{np.nanmax(out.A1_adc):.0f}, "
                  f"A2_adc {np.nanmin(out.A2_adc):.0f}..{np.nanmax(out.A2_adc):.0f}; "
                  f"_V columns assume {vref} V at {2 ** adc_bits - 1}")
        if n_low8:
            print("[usb] WARNING: only b1/b3 (low bytes) were logged in some cycles: ADC counts cannot be rebuilt there "
                  "(A*_lo = adc mod 256, aliased), the upper bits b2/b4 are missing from the file.")
        if n_full == 0 and n_low8 == 0:
            print("[usb] WARNING: no usable cycle found (wrong file, or a different byte layout).")
    return out


def read_raw_usb(path, decode=True, vref=ADC_VREF, adc_bits=ADC_BITS, verbose=True):
    """Read a Usb csv. decode=True -> one row per cycle (see module doc); decode=False -> rows as logged."""
    raw = read_usb_rows(path, verbose=verbose)
    return decode_usb(raw, vref=vref, adc_bits=adc_bits, verbose=verbose) if decode else raw