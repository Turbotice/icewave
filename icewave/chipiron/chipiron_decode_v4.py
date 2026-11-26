"""
chipiron_decode_v4.py -- read / decode the "Usb" csv logged by Phonefleet from the Arduino, for ANY telemetry
sampling period (4 Hz ... a few hundred Hz) and for files recorded with the old 250 ms firmware.

THE ONE NUMBER TO TUNE
    TELEM_PERIOD_MS (below) = the TELEM_PERIOD_MS #define of the PWM_full_v4.ino that was flashed [ms].
        250 -> 4 Hz (old rate)      10 -> 100 Hz      5 -> 200 Hz
    It can also be passed per call (period_ms=...), or stored next to the data in <folder>/sampling_period_ms.txt
    (chipiron_session_v6.py writes it; Chipiron_plotting_v4.py reads it, else asks once).

DEPENDENCIES   numpy, pandas.  No other Chipiron script (escon_edc_v3 can then convert volts to rpm / A).
USED BY        Chipiron_plotting_v4.py, chipiron_session_v6.py (period file only)

EXAMPLES
    from chipiron_decode_v4 import read_raw_usb
    d = read_raw_usb("P101_D2026-10-15T08_45_19-Usb-1-....csv", period_ms=10)   # one row per telemetry frame
    d = read_raw_usb(path, period_ms=250)                                       # old firmware / old files
    raw = read_raw_usb(path, decode=False)                                      # rows exactly as logged
    d[["t_s", "A1_adc", "A1_V"]].head()
    from escon_edc_v3 import load_edc
    load_edc("parameters_Sep26.edc").apply_to_frame(d, wiring={"A1": "AO1", "A2": "AO2"})   # A1_speed_rpm, A2_current_A
    # the 4 bytes of each frame, by role (not by row parity): d[["b1", "b2", "b3", "b4"]]
    # other Arduino reference / resolution:  read_raw_usb(path, vref=3.3, adc_bits=12)

What the file contains (csv, no header): timestamp in microseconds (phone clock, only differences matter), one
byte per row (signed int8 as logged, e.g. -43 = 213), and an empty last column. Unused rows of the phone buffer
are placeholders (timestamp -1) and are dropped.

What the Arduino sends (PWM_full_v4.ino), one 4-byte frame every TELEM_PERIOD_MS:
    b1 = ADC(A1) & 0xFF      b2 = (ADC(A1) >> 4) | 0x80      b3 = ADC(A2) & 0xFF      b4 = ADC(A2) >> 4
    decode: adc = ((b2 & 0x7F) << 4) | (b1 & 0x0F)     check: (b1 >> 4) == (b2 & 0x0F)
    bit 7 of b2 = START-OF-FRAME FLAG: frames are found from the bytes alone, whatever the period.
A1 / A2 are the Arduino PINS (analogRead(A1), analogRead(A2)); what they measure depends on the wiring.
Wiring: A1 <- ESCON AO1, A2 <- ESCON AO2. What AO1 / AO2 MEASURE is set by the .edc (function of each analog output);
in every .edc used so far (parameters_Chipiron_bipolar2640_..., parameters_Chipiron_unipolar_...) it is
    AO1 = actual SPEED      AO2 = actual CURRENT      ->  A1 = speed, A2 = current
(checked on 2026-11-08 on the staircase run, the 250 ms sine run and the 10 ms sine run of 2026-11-06: the sine of
500 rpm commanded comes out as 498 rpm on A1 and as 6000+ rpm if A1 / A2 are swapped).

THREE FRAMINGS, chosen automatically (framing="auto"):
    "nibble": NO flag, fast period (< 160 ms), no 250 ms gaps: frames found from the two nibble checks (a1lo>>4 ==
             a1hi&0xF, same for A2). A1/A2 order is AMBIGUOUS from the bytes alone (see _scan_nibble, phase=): it is
             guessed (A1 = speed = the smoother channel) and a warning is printed. Flash PWM_full_v4.ino (flag) to avoid it.
    "flag" : firmware v4 (flag in b2). Works for every period.
    "gap"  : OLD firmware (no flag, bytes at 0 / 10 / 50 / 60 ms of each 250 ms cycle): cycles are found from the
             time gaps (> 100 ms between cycles). Only valid for that old ~250 ms layout.
    Auto = "flag" if the flag pattern is found in the file, else "gap" if the rows come in clusters of 4 separated by
    long gaps (old 250 ms firmware: this is recognised from the TIMESTAMPS, whatever period you typed), else "nibble"
    (period < 160 ms) or "gap" (period >= 160 ms).

TIME AXIS (time_mode="auto"):
    "phone"   t_s from the phone timestamps of the first byte of each frame.
    "uniform" t_s = frame_number * period. Chosen automatically when the phone stamps the bytes in batches (many
              frames with identical / almost identical timestamps, which happens at high rates because USB reads
              are chunked). It assumes no frame was lost; the number of lost frames is estimated and printed.

Output columns of decode_usb (one row per frame):
    A1_adc, A2_adc  decoded ADC counts (0..1023), NOT converted;  A1_V, A2_V  volts at the pin (counts * vref / 1023)
    b1..b4          the 4 raw bytes of the frame, by ROLE (b1/b3 = low bytes, b2/b4 = adc >> 4; b2 without the flag)
    kind            "full" (all 4 bytes), or "low8" (old firmware only: just b1 and b3 were logged, so only the low
                    8 bits of the ADC are known: A1_lo / A2_lo, A*_adc = NaN)
    df.attrs        period_ms, framing, time_mode actually used
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

# ======================= USER CONFIGURATION =======================
TELEM_PERIOD_MS = 10          # = TELEM_PERIOD_MS in PWM_full_v4.ino [ms]  (250 = old 4 Hz, 10 = 100 Hz)
ADC_VREF = 5.06 # 8/11 measured on the Amundsen, real Chipiron #4.88 PMMH lab 
# volts at full scale (measured on Arduino Mega)
ADC_BITS = 10
# ==================================================================
MIN_PERIOD_MS = 2.0             # 4 bytes at 115200 baud + 2 analogRead need ~1 ms: below this the firmware cannot keep up
PERIOD_FILE = "sampling_period_ms.txt"
FLAG_BIT = 0x80                 # bit 7 of b2 = start-of-frame
CYCLE_GAP_US = 100_000          # "gap" framing: gap above this = new telemetry cycle (old 250 ms firmware)
LOW8_DT_US = (30_000, 80_000)   # "gap" framing: spacing of b1 -> b3 when only two bytes per cycle were logged
FLAG_MIN_FRACTION = 0.2         # auto: "flag" if at least this fraction of (rows / 4) are valid flagged frames
LOST_TOL = 0.02                 # lost-frame warning above this fraction


# ----------------------------------------------------------------------------- sampling period
def check_period(period_ms=None):
    """period_ms or the module default -> float ms; ValueError if impossible."""
    p = TELEM_PERIOD_MS if period_ms is None else period_ms
    try:
        p = float(p)
    except (TypeError, ValueError):
        raise ValueError(f"sampling period must be a number of ms, got {period_ms!r}")
    if not np.isfinite(p) or p < MIN_PERIOD_MS:
        raise ValueError(f"sampling period {p:g} ms is not possible (minimum {MIN_PERIOD_MS:g} ms: 4 bytes at "
                         f"115200 baud plus two analogRead). Use the TELEM_PERIOD_MS of PWM_full_v4.ino.")
    return p


def read_period_file(folder):
    """Period [ms] stored in <folder>/sampling_period_ms.txt, or None if missing / unreadable."""
    try:
        v = float((Path(folder) / PERIOD_FILE).read_text(encoding="utf-8").split()[0].replace(",", "."))
    except (OSError, ValueError, IndexError):
        return None
    return v if np.isfinite(v) and v >= MIN_PERIOD_MS else None


def write_period_file(folder, period_ms):
    """Store the period next to the data (read back by read_period_file / the plotting script)."""
    p = check_period(period_ms)
    (Path(folder) / PERIOD_FILE).write_text(
        f"{p:g}\n# telemetry sampling period [ms] = TELEM_PERIOD_MS of PWM_full_v4.ino when this folder was recorded\n",
        encoding="utf-8")


# ----------------------------------------------------------------------------- raw rows
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


# ----------------------------------------------------------------------------- framing "flag"
def _scan_flagged(b, adc_bits):
    """Start index of every valid flagged frame, scanning sequentially (after a frame, the next one is 4 rows later;
    after a lost byte the scan resynchronises on the next valid pattern)."""
    n, hb = len(b), adc_bits - 4
    if n < 4:
        return np.array([], dtype=np.int64)
    i = np.arange(n - 3)
    good = ((((b[i] >> 4) & 0xF) == (b[i + 1] & 0xF))              # A1 nibble check
            & ((b[i + 1] & FLAG_BIT) != 0)                          # start-of-frame flag
            & (((b[i + 1] & 0x7F) >> hb) == 0)                      # b2 fits in adc_bits - 4 bits
            & ((((b[i + 2] >> 4) & 0xF)) == (b[i + 3] & 0xF))       # A2 nibble check
            & ((b[i + 3] >> hb) == 0))                              # b4 has no flag and fits
    starts, nxt = [], 0
    for p in np.flatnonzero(good):
        if p >= nxt:
            starts.append(p)
            nxt = p + 4
    return np.asarray(starts, dtype=np.int64)


def _scan_nibble(b, adc_bits, phase="auto", verbose=True):
    """FLAGLESS frames (firmware without the start-of-frame flag, e.g. a v3-type layout run at 10 ms).
    A frame is a1lo, a1hi, a2lo, a2hi with the two nibble checks and a1hi, a2hi < 2**(adc_bits-4). Those checks
    cannot tell (a1lo, a1hi, a2lo, a2hi) from (a2lo, a2hi, a1lo', a1hi'): both alignments (2 rows apart) are valid
    everywhere, so which pair is A1 is decided by `phase`:
        "auto" : the alignment where A1 (ESCON AO1 = SPEED, smooth) is the SMOOTHER channel and A2 (AO2 = CURRENT,
                 noisier) the other one. This is only a guess: it can fail on a constant-speed hold, where both
                 channels are flat. The only safe fix is the flag of PWM_full_v4.ino.
        0 / 2  : first valid alignment / the one 2 rows later (force it if auto guesses wrong)
    After a lost byte the scan resynchronises on the valid pattern 3..5 rows after the previous frame."""
    n, hb = len(b), adc_bits - 4
    if n < 8:
        return np.array([], dtype=np.int64)
    i = np.arange(n - 3)
    good = ((((b[i] >> 4) & 0xF) == (b[i + 1] & 0xF)) & ((b[i + 1] >> hb) == 0)
            & (((b[i + 2] >> 4) & 0xF) == (b[i + 3] & 0xF)) & ((b[i + 3] >> hb) == 0))

    def run(first):
        starts, p = [], first
        while p is not None:
            starts.append(p)
            if p + 4 < len(good) and good[p + 4]:
                p += 4
                continue
            nxt = np.flatnonzero(good[p + 3:])
            p = None if len(nxt) == 0 else p + 3 + int(nxt[0])
        return np.asarray(starts, dtype=np.int64)

    cand = np.flatnonzero(good)
    if len(cand) == 0:
        return np.array([], dtype=np.int64)
    p0 = int(cand[0])
    if phase in (0, "0"):
        return run(p0)
    if phase in (2, "2"):
        return run(p0 + 2)
    s0, s2 = run(p0), run(p0 + 2)

    def rough(st, lo):                                       # mean |step| of a channel [ADC counts]
        v = (((b[st + 1] << 4) | (b[st] & 0xF)) if lo == 1 else ((b[st + 3] << 4) | (b[st + 2] & 0xF))).astype(float)
        return float(np.mean(np.abs(np.diff(v)))) if len(v) > 2 else np.nan
    r0 = (rough(s0, 1), rough(s0, 2))                        # (A1, A2) roughness if alignment 0 is right
    r2 = (rough(s2, 1), rough(s2, 2))
    # alignment 0 is right if A1 (speed) is smoother than A2 (current) there, and more so than in alignment 2
    score0, score2 = r0[0] / max(r0[1], 1e-9), r2[0] / max(r2[1], 1e-9)
    pick0 = score0 <= score2
    weak = max(score0, score2) / max(min(score0, score2), 1e-9) < 2.0 or min(score0, score2) > 0.8
    if verbose:
        print(f"[usb] flagless frames: A1/A2 cannot be told apart from the bytes alone. Roughness A1 | A2 [counts/frame]: "
              f"alignment 0 -> {r0[0]:.1f} | {r0[1]:.1f};  alignment 2 -> {r2[0]:.1f} | {r2[1]:.1f}. "
              f"Using alignment {0 if pick0 else 2} (A1 = the smoother channel = speed, A2 = current). "
              "If A1/A2 look swapped, pass phase=0 or phase=2 (decode_usb / read_raw_usb) to force the other one.")
        if weak:
            print("[usb] WARNING: the two channels are about equally smooth: A1 / A2 are NOT reliably identified "
                  "(constant-speed hold?). Check the speed column against the command, or flash PWM_full_v4.ino (flag).")
    return s0 if pick0 else s2


def _decode_flag(t, b, s, adc_bits, period, time_mode, verbose, label="flag"):
    b1, b2, b3, b4 = b[s], b[s + 1] & 0x7F, b[s + 2], b[s + 3]
    a1 = ((b2 << 4) | (b1 & 0xF)).astype(float)
    a2 = ((b4 << 4) | (b3 & 0xF)).astype(float)
    t_f = t[s]
    # ---- time axis
    used, n = "phone", len(s)
    med = np.nan
    if n >= 3:
        d_ms = np.diff(t_f) / 1e3
        med = float(np.median(d_ms))
        batched = float(np.mean(d_ms < 0.2 * period)) > 0.2
        if time_mode == "uniform" or (time_mode == "auto" and batched):
            used = "uniform"
        if verbose and batched:
            print(f"[usb] the phone stamps the bytes in batches ({100 * np.mean(d_ms < 0.2 * period):.0f} % of the "
                  f"frame spacings are < {0.2 * period:.3g} ms)" + (": time axis rebuilt as frame_number x period."
                  if used == "uniform" else "."))
        if verbose and abs(med - period) > 0.3 * period and not batched:
            print(f"[usb] WARNING: median frame spacing {med:.3g} ms but the sampling period is {period:g} ms. "
                  "Either TELEM_PERIOD_MS here differs from the one in the flashed PWM_full_v4.ino, or many frames "
                  "were lost.")
    elif time_mode == "uniform":
        used = "uniform"
    t_s = np.arange(n) * period / 1e3 if used == "uniform" else (t_f - t[0]) / 1e6
    out = pd.DataFrame({"t_us": t_f, "kind": "full", "A1_lo": b1, "A2_lo": b3, "A1_adc": a1, "A2_adc": a2,
                        "b1": b1, "b2": b2, "b3": b3, "b4": b4})
    out["t_s"] = t_s
    if verbose:
        print(f"[usb] framing '{label}': {n} frame(s) from {len(b)} rows, {len(b) - 4 * n} byte(s) skipped "
              f"while resynchronising; sampling period {period:g} ms ({1000 / period:.3g} Hz); time axis '{used}'")
        dur = (t[-1] - t[0]) / 1e6
        expected = dur * 1000 / period + 1
        if n and n < (1 - LOST_TOL) * expected - 1:
            print(f"[usb] WARNING: ~{expected - n:.0f} frame(s) missing ({100 * (1 - n / expected):.1f} % of the "
                  f"{expected:.0f} expected in {dur:.1f} s at {period:g} ms). Wrong period, or the phone lost bytes."
                  + (" The uniform time axis assumes none is lost: it is compressed." if used == "uniform" else ""))
    return out, used


# ----------------------------------------------------------------------------- framing "gap" (old firmware)
def _looks_gapped(t, gap_us=CYCLE_GAP_US):
    """True if the rows come in clusters of 4 separated by gaps > gap_us (old 250 ms firmware). Independent of the
    period typed by the user, so an old file is never decoded as 'nibble' (whose A1/A2 order is only a guess)."""
    if len(t) < 40:
        return False
    starts = np.r_[0, np.where(np.diff(t) > gap_us)[0] + 1]
    counts = np.diff(np.r_[starts, len(t)])
    return len(starts) > 10 and float(np.mean(counts == 4)) > 0.8


def _decode_gap(t, b, adc_bits, gap_us, verbose):
    starts = np.r_[0, np.where(np.diff(t) > gap_us)[0] + 1]
    counts = np.diff(np.r_[starts, len(t)])
    cols = {k: [] for k in ("t_us", "kind", "A1_lo", "A2_lo", "A1_adc", "A2_adc", "b1", "b2", "b3", "b4")}

    s = starts[counts == 4]                                  # full cycles: 4 rows
    b1, b2, b3, b4 = b[s], b[s + 1], b[s + 2], b[s + 3]
    hi_bits = adc_bits - 4
    ok1 = (((b1 >> 4) & 0xF) == (b2 & 0xF)) & ((b2 >> hi_bits) == 0)
    ok2 = (((b3 >> 4) & 0xF) == (b4 & 0xF)) & ((b4 >> hi_bits) == 0)
    a1 = np.where(ok1, (b2 << 4) | (b1 & 0xF), np.nan)
    a2 = np.where(ok2, (b4 << 4) | (b3 & 0xF), np.nan)
    for k, v in zip(cols, (t[s], np.full(len(s), "full"), b1, b3, a1, a2, b1, b2, b3, b4)):
        cols[k].append(v)

    s2 = starts[counts == 2]                                 # low-byte-only cycles: 2 rows ~50 ms apart
    dt = t[s2 + 1] - t[s2]
    s2 = s2[(dt >= LOW8_DT_US[0]) & (dt <= LOW8_DT_US[1])]
    nan = np.full(len(s2), np.nan)
    for k, v in zip(cols, (t[s2], np.full(len(s2), "low8"), b[s2], b[s2 + 1], nan, nan, b[s2], nan, b[s2 + 1], nan)):
        cols[k].append(v)

    out = pd.DataFrame({k: np.concatenate(v) for k, v in cols.items()}).sort_values("t_us").reset_index(drop=True)
    out["t_us"] = out["t_us"].astype("int64")
    out["t_s"] = (out["t_us"] - t[0]) / 1e6
    if verbose:
        n_full, n_low8 = int((out.kind == "full").sum()), int((out.kind == "low8").sum())
        n_drop = len(starts) - n_full - n_low8
        print(f"[usb] framing 'gap' (old firmware): {len(starts)} cycle(s): {n_full} full, {n_low8} low-byte-only, "
              f"{n_drop} dropped (row counts per cycle: {dict(zip(*np.unique(counts, return_counts=True)))})")
        if n_full:
            print(f"[usb] consistency check failed (A -> NaN): A1 {int((~ok1).sum())}, A2 {int((~ok2).sum())} of {n_full}")
        if n_low8:
            print("[usb] WARNING: only b1/b3 (low bytes) were logged in some cycles: ADC counts cannot be rebuilt there "
                  "(A*_lo = adc mod 256, aliased), the upper bits b2/b4 are missing from the file.")
    return out


# ----------------------------------------------------------------------------- main entry point
def decode_usb(raw, vref=ADC_VREF, adc_bits=ADC_BITS, period_ms=None, framing="auto", time_mode="auto",
               gap_us=CYCLE_GAP_US, verbose=True, phase="auto"):
    """Decode raw rows (from read_usb_rows) into one row per telemetry frame (see the module docstring).
    period_ms : TELEM_PERIOD_MS of the firmware (None = module default TELEM_PERIOD_MS)
    framing   : "auto" | "flag" | "nibble" | "gap"        time_mode : "auto" | "phone" | "uniform"
    phase     : "nibble" framing only: "auto" | 0 | 2 (which byte pair is A1, see _scan_nibble)
    Returns a (possibly EMPTY) DataFrame: check len(df) before plotting."""
    if framing not in ("auto", "flag", "nibble", "gap"):
        raise ValueError("framing must be 'auto', 'flag', 'nibble' or 'gap'")
    if time_mode not in ("auto", "phone", "uniform"):
        raise ValueError("time_mode must be 'auto', 'phone' or 'uniform'")
    if not 4 <= adc_bits <= 11:
        raise ValueError("adc_bits must be between 4 and 11 (the start-of-frame flag needs bit 7 of b2 free)")
    period = check_period(period_ms)
    t = raw["t_us"].to_numpy(np.int64)
    b = raw["byte"].to_numpy(np.int64)
    if len(t) < 4:
        raise ValueError("not enough Usb rows to decode (need at least 4)")

    used, s = framing, np.array([], dtype=np.int64)
    if framing in ("auto", "flag"):
        s = _scan_flagged(b, adc_bits)
        if framing == "auto":
            if len(s) >= FLAG_MIN_FRACTION * (len(b) // 4):
                used = "flag"
            else:                                   # no flag: old firmware (slow, time gaps) or flagless fast data
                used = "gap" if _looks_gapped(t) else ("nibble" if period < 160 else "gap")
                if verbose:
                    print(f"[usb] no start-of-frame flag found: file recorded WITHOUT the v4 flag "
                          f"(firmware other than PWM_full_v4.ino as given) -> framing '{used}'.")
    if used == "nibble":
        s = _scan_nibble(b, adc_bits, phase, verbose)
    tmode = "phone"
    if used in ("flag", "nibble"):
        out, tmode = _decode_flag(t, b, s, adc_bits, period, time_mode, verbose, label=used)
    else:
        if verbose and framing == "auto":
            print("[usb] no start-of-frame flag found in this file: decoding it as OLD-firmware data (time gaps).")
        if verbose and period < 160 and not _looks_gapped(t):
            print(f"[usb] WARNING: gap framing only works for the old ~250 ms layout, but the sampling period is set "
                  f"to {period:g} ms. If this file comes from the old firmware, call with period_ms=250; if it comes "
                  "from PWM_full_v4.ino, the period or the flashed firmware is not what you think.")
        out = _decode_gap(t, b, adc_bits, gap_us, verbose)
    out["A1_V"] = out["A1_adc"] * vref / (2 ** adc_bits - 1)
    out["A2_V"] = out["A2_adc"] * vref / (2 ** adc_bits - 1)
    out.attrs.update(period_ms=period, framing=used, time_mode=tmode)

    if verbose:
        if len(out) == 0:
            print(f"[usb] WARNING: no usable frame found with framing '{used}' and period {period:g} ms "
                  "(wrong file, wrong period, or a different byte layout).")
        elif out["A1_adc"].notna().any():
            print(f"[usb] A1_adc {np.nanmin(out.A1_adc):.0f}..{np.nanmax(out.A1_adc):.0f}, "
                  f"A2_adc {np.nanmin(out.A2_adc):.0f}..{np.nanmax(out.A2_adc):.0f}; "
                  f"_V columns assume {vref} V at {2 ** adc_bits - 1}")
    return out


def read_raw_usb(path, decode=True, vref=ADC_VREF, adc_bits=ADC_BITS, period_ms=None, framing="auto",
                 time_mode="auto", verbose=True, phase="auto"):
    """Read a Usb csv. decode=True -> one row per frame (see module doc); decode=False -> rows as logged."""
    raw = read_usb_rows(path, verbose=verbose)
    if not decode:
        return raw
    return decode_usb(raw, vref=vref, adc_bits=adc_bits, period_ms=period_ms, framing=framing,
                      time_mode=time_mode, verbose=verbose, phase=phase)