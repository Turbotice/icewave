#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

"""
Unified geophone data unloading script.

Replaces the fragmented Move_data.py and prepare_smartsolo.py scripts with one interactive, transparent, two-phase workflow:

  Phase 1 — Metadata: Copy LOG files from geophone SD cards,
            3 geophones at a time, renaming them by geophone number.

  Phase 2 — Seismic data: Rename .miniseed files exported by Sololite
            (serial number → geophone number) and sort them into
            per-acquisition folders.

Dependencies:
    pip install prompt_toolkit pyyaml

Usage:
    python unload_geophones.py

Author: jonas + claude
Date:   2026-02-13
"""

import os
import sys
import re
import shutil
import textwrap
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("ERROR: PyYAML is required.  Install with:  pip install pyyaml")

try:
    from prompt_toolkit import prompt as pt_prompt
    from prompt_toolkit.completion import PathCompleter
except ImportError:
    sys.exit(
        "ERROR: prompt_toolkit is required.  Install with:  pip install prompt_toolkit"
    )

# ──────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────
CONFIG_PATH = Path.home() / "geophone_unload_config.yaml"
DEFAULT_DRIVES = ["H:\\", "G:\\", "D:\\"]
MINISEED_EXT = ".miniseed"


# Colours for terminal output — auto-disabled when not supported
def _supports_ansi() -> bool:
    """Return True if the terminal likely supports ANSI escape codes."""
    if not hasattr(sys.stdout, "isatty") or not sys.stdout.isatty():
        return False
    if sys.platform == "win32":
        # Windows 10 build 14393+ supports ANSI via Virtual Terminal Processing.
        # Try to enable it; if it fails, fall back to no-colour.
        try:
            import ctypes

            kernel32 = ctypes.windll.kernel32  # type: ignore[attr-defined]
            handle = kernel32.GetStdHandle(-11)  # STD_OUTPUT_HANDLE
            mode = ctypes.c_ulong()
            kernel32.GetConsoleMode(handle, ctypes.byref(mode))
            # ENABLE_VIRTUAL_TERMINAL_PROCESSING = 0x0004
            if not (mode.value & 0x0004):
                kernel32.SetConsoleMode(handle, mode.value | 0x0004)
            return True
        except Exception:
            return False
    return True  # Unix / macOS


if _supports_ansi():
    C_RESET = "\033[0m"
    C_BOLD = "\033[1m"
    C_GREEN = "\033[92m"
    C_YELLOW = "\033[93m"
    C_RED = "\033[91m"
    C_CYAN = "\033[96m"
else:
    C_RESET = C_BOLD = C_GREEN = C_YELLOW = C_RED = C_CYAN = ""


# ──────────────────────────────────────────────────────────────────────
# Utility helpers
# ──────────────────────────────────────────────────────────────────────
def info(msg: str) -> None:
    print(f"{C_GREEN}[INFO]{C_RESET}  {msg}")


def warn(msg: str) -> None:
    print(f"{C_YELLOW}[WARN]{C_RESET}  {msg}")


def error(msg: str) -> None:
    print(f"{C_RED}[ERROR]{C_RESET} {msg}")


def action(msg: str) -> None:
    """Log a file-operation line (copy / rename / move)."""
    print(f"  {C_CYAN}→{C_RESET} {msg}")


def header(title: str) -> None:
    width = 60
    print()
    print(f"{C_BOLD}{'═' * width}{C_RESET}")
    print(f"{C_BOLD}  {title}{C_RESET}")
    print(f"{C_BOLD}{'═' * width}{C_RESET}")
    print()


def confirm(msg: str, default_yes: bool = False) -> bool:
    """Ask a yes/no question. Returns True for yes."""
    suffix = "[Y/n]" if default_yes else "[y/N]"
    answer = input(f"{msg} {suffix} ").strip().lower()
    if answer == "":
        return default_yes
    return answer in ("y", "yes")


# ──────────────────────────────────────────────────────────────────────
# Config management  (YAML, persisted to ~/geophone_unload_config.yaml)
# ──────────────────────────────────────────────────────────────────────
def load_config() -> dict:
    if CONFIG_PATH.exists():
        try:
            with open(CONFIG_PATH, "r") as f:
                data = yaml.safe_load(f) or {}
            return data
        except Exception as e:
            warn(f"Could not read config file {CONFIG_PATH}: {e}")
    return {}


def save_config(cfg: dict) -> None:
    try:
        with open(CONFIG_PATH, "w") as f:
            yaml.dump(cfg, f, default_flow_style=False)
    except Exception as e:
        warn(f"Could not save config file {CONFIG_PATH}: {e}")


# ──────────────────────────────────────────────────────────────────────
# Interactive path input  (with tab-completion via prompt_toolkit)
# ──────────────────────────────────────────────────────────────────────
_path_completer_dirs = PathCompleter(expanduser=True, only_directories=True)
_path_completer_all = PathCompleter(expanduser=True, only_directories=False)


def ask_path(
    label: str,
    default: str = "",
    only_directories: bool = True,
    must_exist: bool = False,
) -> str:
    """Prompt the user for a filesystem path with tab-completion.

    Parameters
    ----------
    label : str
        The prompt text shown to the user.
    default : str
        Pre-filled value (from config). User can press Enter to accept.
    only_directories : bool
        If True, tab-completer only shows directories.
    must_exist : bool
        If True, keep asking until the path actually exists on disk.
    """
    completer = _path_completer_dirs if only_directories else _path_completer_all

    while True:
        raw = pt_prompt(
            f"{label}: ",
            default=default,
            completer=completer,
        ).strip()

        if not raw:
            warn("Empty path — please enter something.")
            continue

        expanded = os.path.expanduser(raw)
        if must_exist and not os.path.exists(expanded):
            warn(f"Path does not exist: {expanded}")
            if not confirm("Try again?", default_yes=True):
                return expanded  # user insists
            continue

        return expanded


def normalise_drives(drives: list[str]) -> list[str]:
    """Normalise drive paths: ensure trailing backslash on Windows drive letters."""
    normalised = []
    for d in drives:
        d = d.strip()
        if not d:
            continue
        if len(d) >= 2 and d[1] == ":":
            if not d.endswith("\\"):
                d = d + "\\"
        normalised.append(d)
    return normalised


# ──────────────────────────────────────────────────────────────────────
# Geophones table  (maps 9-digit serial number → 2-digit geophone id)
# ──────────────────────────────────────────────────────────────────────
def load_geophones_table(path: str) -> dict[str, str]:
    """Return {serial_number: geophone_id} dict, e.g. {'453024765': '01'}.

    The file is TSV with a header line 'Num\\tgeo_SN'.
    """
    table: dict[str, str] = {}
    with open(path, "r") as f:
        for lineno, line in enumerate(f, 1):
            stripped = line.strip()
            if not stripped:
                continue
            # Support both CSV (comma) and TSV (tab / whitespace)
            if "," in stripped:
                parts = [p.strip() for p in stripped.split(",")]
            else:
                parts = stripped.split()
            if len(parts) < 2:
                warn(f"geophones_table line {lineno}: malformed → {stripped!r}")
                continue
            num, serial = parts[0], parts[1]
            # Skip header
            if num.lower() == "num" or serial.lower().startswith("geo"):
                continue
            table[serial] = num
    if not table:
        error(f"Geophones table at {path} appears empty or unparseable.")
        sys.exit(1)
    info(f"Loaded geophones table: {len(table)} entries from {path}")
    for serial, num in sorted(table.items(), key=lambda kv: kv[1]):
        print(f"    geophone {num}  ←  serial {serial}")
    return table


def serial_to_geophone(serial: str, table: dict[str, str]) -> str | None:
    """Look up a serial number (9 digits or last-4/5 digits) in the table."""
    # Exact 9-digit match
    if serial in table:
        return table[serial]
    # Try matching by last 4 or 5 digits (some filenames use shorter serials)
    for full_serial, geo_num in table.items():
        if full_serial.endswith(serial):
            return geo_num
    return None


# ──────────────────────────────────────────────────────────────────────
# Phase 1:  Metadata extraction  (replaces Move_data.py)
# ──────────────────────────────────────────────────────────────────────
def parse_digisolo(filepath: str) -> dict:
    """Parse a DigiSolo .txt / .LOG file and extract metadata.

    Returns a dict with keys:
        serial_number : str   (9 digits)
        acquisitions  : list of dicts, each with date/time fields
    """
    result = {"serial_number": None, "acquisitions": []}

    try:
        with open(filepath, "r", errors="replace") as f:
            lines = f.readlines()
    except Exception as e:
        error(f"Cannot read {filepath}: {e}")
        return result

    current_acq: dict | None = None

    for i, raw_line in enumerate(lines):
        line = raw_line.strip()

        # Serial number  (e.g.  "Serial Number = 453024765")
        m = re.search(r"Serial\s+Number\s*=\s*(\d{7,9})", line)
        if m:
            result["serial_number"] = m.group(1)
            continue

        # Start of an acquisition record
        if "Start Acquisition FileName" in line:
            current_acq = {
                "filename": line,
                "date": None,
                "time": None,
                "latitudes": [],
                "longitudes": [],
            }
            result["acquisitions"].append(current_acq)
            # The date is typically on the *previous* line
            if i > 0:
                prev = lines[i - 1].strip()
                # Pattern like: Date = "2025/02/10"  or   2025/02/10, HH:MM:SS
                dm = re.search(r"(\d{4})[/\-](\d{2})[/\-](\d{2})", prev)
                if dm:
                    current_acq["date"] = f"{dm.group(1)}/{dm.group(2)}/{dm.group(3)}"
                # Time often follows a comma or is on the same line
                tm = re.search(r"(\d{2}):(\d{2}):(\d{2})", prev)
                if tm:
                    current_acq["time"] = f"{tm.group(1)}:{tm.group(2)}:{tm.group(3)}"
            continue

        # Notify blocks carry date+time too  (legacy format used by Move_data.py)
        if "[Notify" in line and current_acq is None:
            # Next line typically has the timestamp
            if i + 1 < len(lines):
                nxt = lines[i + 1].strip()
                dm = re.search(r"(\d{4})[/\-](\d{2})[/\-](\d{2})", nxt)
                tm = re.search(r"(\d{2}):(\d{2}):(\d{2})", nxt)
                acq = {
                    "filename": line,
                    "date": None,
                    "time": None,
                    "latitudes": [],
                    "longitudes": [],
                }
                if dm:
                    acq["date"] = f"{dm.group(1)}/{dm.group(2)}/{dm.group(3)}"
                if tm:
                    acq["time"] = f"{tm.group(1)}:{tm.group(2)}:{tm.group(3)}"
                result["acquisitions"].append(acq)
            continue

        # GPS coordinates
        if current_acq is not None:
            lat_m = re.search(r"Latitude\s*=\s*([-\d.]+)", line)
            if lat_m:
                try:
                    current_acq["latitudes"].append(float(lat_m.group(1)))
                except ValueError:
                    pass
            lon_m = re.search(r"Longitude\s*=\s*([-\d.]+)", line)
            if lon_m:
                try:
                    current_acq["longitudes"].append(float(lon_m.group(1)))
                except ValueError:
                    pass

    return result


def run_phase1(cfg: dict, geo_table: dict[str, str]) -> str:
    """Phase 1: copy metadata from geophone drives to the output folder.

    Returns the output_path for use in Phase 2.
    """
    header("PHASE 1 — Metadata extraction")
    print(
        textwrap.dedent(
            """\
        Plug in 3 geophones at a time. This phase will:
          • scan the mounted drives for .LOG files
          • extract serial number → rename by geophone number → copy to output
    """
        )
    )

    # --- Ask for output path ---
    output_path = ask_path(
        "Output folder (metadata + data will go here)",
        default=cfg.get("all_data_output_path", ""),
        only_directories=True,
        must_exist=False,
    )
    cfg["all_data_output_path"] = output_path
    save_config(cfg)

    # Create output dir
    if os.path.isdir(output_path):
        contents = os.listdir(output_path)
        if contents:
            warn(f"Output folder already exists and contains {len(contents)} items:")
            for item in contents[:10]:
                print(f"    {item}")
            if len(contents) > 10:
                print(f"    … and {len(contents) - 10} more")
            if not confirm("Continue anyway?"):
                info("Aborted by user.")
                sys.exit(0)
    else:
        os.makedirs(output_path, exist_ok=True)
        info(f"Created output folder: {output_path}")

    metadata_dir = os.path.join(output_path, "metadata")
    os.makedirs(metadata_dir, exist_ok=True)

    # --- Batch loop: 3 geophones at a time ---
    drives = normalise_drives(cfg.get("geophone_drives", DEFAULT_DRIVES))
    info(f"Using geophone drives from config: {', '.join(drives)}")
    batch_number = 0
    total_files_copied = 0

    while True:
        batch_number += 1
        header(f"Metadata batch #{batch_number}")

        batch_copied = 0
        for drive in drives:
            info(f"Scanning drive: {drive}")
            if not os.path.isdir(drive):
                warn(f"  Drive {drive} is not accessible — skipping.")
                continue

            found_any = False
            for fname in os.listdir(drive):
                if not fname.upper().endswith(".LOG"):
                    continue

                found_any = True
                src = os.path.join(drive, fname)
                meta = parse_digisolo(src)
                serial = meta["serial_number"]
                geo_num = serial_to_geophone(serial, geo_table) if serial else None

                if geo_num:
                    dest_name = f"LOG_{geo_num}{os.path.splitext(fname)[1]}"
                    dest = os.path.join(metadata_dir, dest_name)
                    if os.path.exists(dest):
                        warn(f"LOG_{geo_num} already exists")
                        if not confirm("  Overwrite?"):
                            continue
                    shutil.copy2(src, dest)
                    n_acq = len(meta.get("acquisitions", []))
                    info(
                        f"Extracted log from geophone {geo_num} (SN: {serial}, {n_acq} acq)"
                    )
                    batch_copied += 1
                else:
                    shutil.copy2(src, os.path.join(metadata_dir, fname))
                    warn(f"Unknown serial in {fname} — copied as-is")
                    batch_copied += 1

            if not found_any:
                warn(f"No .LOG file found on {drive}")

        # Warn if fewer than 3 geophones found in this batch
        geophones_found = sum(
            1
            for d in drives
            if os.path.isdir(d)
            and any(f.upper().endswith(".LOG") for f in os.listdir(d))
        )
        if geophones_found < 3:
            warn(f"Only {geophones_found} geophone(s) with .LOG data (expected 3)")

        total_files_copied += batch_copied
        info(
            f"Batch #{batch_number}: {batch_copied} log(s) copied ({total_files_copied} total)"
        )

        print()
        choice = (
            pt_prompt(
                "Plug in next set of geophones and press Enter, "
                "or type 'done' to proceed to data phase: ",
            )
            .strip()
            .lower()
        )

        if choice in ("done", "d", "quit", "q", "exit"):
            break

    info(f"Phase 1 complete — {total_files_copied} log file(s) saved to {metadata_dir}")
    return output_path


# ──────────────────────────────────────────────────────────────────────
# Phase 2:  Seismic data organisation  (replaces prepare_smartsolo.py)
# ──────────────────────────────────────────────────────────────────────
def rename_miniseed_files(directory: str, geo_table: dict[str, str]) -> int:
    """Step 2a: rename .miniseed files from serial-number prefix to
    2-digit geophone number.

    Returns the number of files renamed.
    """
    info("Step 2a — Renaming .miniseed files (serial → geophone number)…")
    count = 0
    skipped = 0

    files = sorted(os.listdir(directory))
    for fname in files:
        if not fname.lower().endswith(MINISEED_EXT):
            continue

        # The first 9 characters are the serial number
        if len(fname) < 10:
            warn(f"  Filename too short to contain serial: {fname}")
            skipped += 1
            continue

        serial_candidate = fname[:9]
        if not serial_candidate.isdigit():
            # Already renamed (starts with 2-digit geophone number)?
            if fname[:2].isdigit() and fname[2] == ".":
                continue  # already processed
            warn(f"  Cannot extract serial from: {fname} — skipping")
            skipped += 1
            continue

        geo_num = serial_to_geophone(serial_candidate, geo_table)
        if geo_num is None:
            warn(
                f"  Serial {serial_candidate} not in geophones table "
                f"— skipping {fname}"
            )
            skipped += 1
            continue

        new_fname = geo_num + fname[9:]
        old_path = os.path.join(directory, fname)
        new_path = os.path.join(directory, new_fname)

        if os.path.exists(new_path):
            warn(f"  Target already exists: {new_fname} — skipping")
            skipped += 1
            continue

        os.rename(old_path, new_path)
        action(f"RENAME  {fname}  →  {new_fname}")
        count += 1

    info(f"  Renamed {count} file(s), skipped {skipped}.")
    return count


def fix_malformed_filenames(directory: str) -> int:
    """Step 2b: fix common malformations in .miniseed filenames.

    1. Replace double dots  '..' → '.'
    2. Zero-pad single-digit acquisition numbers
       e.g. '01.1.2024...' → '01.001.2024...'

    Returns the number of files fixed.
    """
    info("Step 2b — Fixing malformed filenames…")
    count = 0

    files = sorted(os.listdir(directory))
    for fname in files:
        if not fname.lower().endswith(MINISEED_EXT):
            continue

        new_fname = fname

        # Fix double dots (at most 2 replacements, like original code)
        if ".." in new_fname:
            new_fname = new_fname.replace("..", ".", 2)

        # Zero-pad acquisition number:
        # Pattern: "NN.D." where D is a single digit → "NN.00D."
        # e.g. "01.1.2024" → "01.001.2024"
        m = re.match(r"^(\d{2})\.(\d)\.(.+)$", new_fname)
        if m:
            geo = m.group(1)
            acq_digit = m.group(2)
            rest = m.group(3)
            new_fname = f"{geo}.00{acq_digit}.{rest}"

        # Also handle 2-digit acq numbers: "01.12." → "01.012."
        m2 = re.match(r"^(\d{2})\.(\d{2})\.(.+)$", new_fname)
        if m2:
            geo = m2.group(1)
            acq_digits = m2.group(2)
            rest = m2.group(3)
            # Only pad if it looks like a short acq number (not a year like 20xx)
            if not acq_digits.startswith("20"):
                new_fname = f"{geo}.0{acq_digits}.{rest}"

        if new_fname != fname:
            old_path = os.path.join(directory, fname)
            new_path = os.path.join(directory, new_fname)
            if os.path.exists(new_path):
                warn(f"  Target already exists: {new_fname} — skipping fix for {fname}")
                continue
            os.rename(old_path, new_path)
            action(f"FIX     {fname}  →  {new_fname}")
            count += 1

    info(f"  Fixed {count} filename(s).")
    return count


def sort_into_acquisition_folders(source_dir: str, dest_base: str) -> int:
    """Step 2c: move .miniseed files into per-acquisition sub-folders.

    The acquisition number is the dot-separated field immediately before
    the 4-digit year (20xx) in the filename:
        01.000.1.2026.02.13...E.miniseed  →  acq '001'
        01.000.2.2026.02.13...E.miniseed  →  acq '002'

    Files are moved from source_dir into dest_base/<acq_number>/.

    Returns the number of files moved.
    """
    info("Step 2c — Sorting files into acquisition folders…")
    count = 0

    files = sorted(os.listdir(source_dir))
    for fname in files:
        if not fname.lower().endswith(MINISEED_EXT):
            continue

        parts = fname.split(".")
        if len(parts) < 3:
            warn(f"  Cannot parse: {fname} — skipping")
            continue

        # Find the acquisition number: field just before the 4-digit year
        acq_number = None
        for idx in range(2, len(parts)):
            if len(parts[idx]) == 4 and parts[idx].startswith("20"):
                acq_number = parts[idx - 1].zfill(3)
                break

        if acq_number is None:
            # Fallback: use parts[1]
            acq_number = parts[1]

        acq_folder = os.path.join(dest_base, acq_number)
        os.makedirs(acq_folder, exist_ok=True)

        src_path = os.path.join(source_dir, fname)
        dst_path = os.path.join(acq_folder, fname)

        if os.path.exists(dst_path):
            warn(f"  Already exists: {dst_path} — skipping")
            continue

        shutil.move(src_path, dst_path)
        count += 1

    # Report what acquisition folders were created
    acq_folders = sorted(
        d
        for d in os.listdir(dest_base)
        if os.path.isdir(os.path.join(dest_base, d)) and d != "metadata"
    )
    acq_with_data = [
        (
            d,
            sum(
                1
                for f in os.listdir(os.path.join(dest_base, d))
                if f.lower().endswith(MINISEED_EXT)
            ),
        )
        for d in acq_folders
    ]
    acq_with_data = [(d, n) for d, n in acq_with_data if n > 0]
    if acq_with_data:
        info(f"Sorted {count} file(s) into {len(acq_with_data)} acquisition folder(s):")
        for folder, n in acq_with_data:
            print(f"    {folder}/  — {n} file(s)")
    else:
        info(f"{count} file(s) moved.")
    return count


def run_phase2(cfg: dict, geo_table: dict[str, str], output_path: str) -> None:
    """Phase 2: rename and sort the Sololite-exported .miniseed files."""
    header("PHASE 2 — Seismic data organisation")
    print(
        textwrap.dedent(
            """\
        This phase will:
          1. Rename .miniseed files (serial number → geophone number)
          2. Fix malformed filenames (double dots, acquisition number padding)
          3. Sort files into per-acquisition sub-folders
    """
        )
    )

    # Ask where Sololite exported the data
    sololite_path = ask_path(
        "Path to Sololite-exported .miniseed files",
        default=cfg.get("sololite_data_export_path", ""),
        only_directories=True,
        must_exist=True,
    )
    cfg["sololite_data_export_path"] = sololite_path
    save_config(cfg)

    # Count miniseed files
    all_files = os.listdir(sololite_path)
    miniseed_files = [f for f in all_files if f.lower().endswith(MINISEED_EXT)]
    info(f"Found {len(miniseed_files)} .miniseed file(s) in {sololite_path}")

    if not miniseed_files:
        warn("No .miniseed files found — nothing to do.")
        return

    # We work on a copy in the output folder to avoid modifying the Sololite
    # export directory.  The user can also choose to work in-place.
    print()
    work_in_place = confirm(
        f"Copy files to output folder ({output_path}) before processing?\n"
        f"  If No, files will be renamed/moved IN-PLACE in {sololite_path}.",
        default_yes=True,
    )

    if work_in_place:
        work_dir = output_path
        info(f"Copying {len(miniseed_files)} file(s) to {work_dir} …")
        for fname in miniseed_files:
            src = os.path.join(sololite_path, fname)
            dst = os.path.join(work_dir, fname)
            if os.path.exists(dst):
                continue  # already there
            shutil.copy2(src, dst)
        info("Copy complete.")
    else:
        work_dir = sololite_path

    print()

    # Step 2a: rename
    rename_miniseed_files(work_dir, geo_table)
    print()

    # Step 2b: fix filenames
    fix_malformed_filenames(work_dir)
    print()

    # Step 2c: sort into folders
    sort_into_acquisition_folders(work_dir, work_dir)


# ──────────────────────────────────────────────────────────────────────
# Summary
# ──────────────────────────────────────────────────────────────────────
def print_summary(output_path: str) -> None:
    header("SUMMARY")

    metadata_dir = os.path.join(output_path, "metadata")
    if os.path.isdir(metadata_dir):
        meta_files = os.listdir(metadata_dir)
        info(f"Metadata files: {len(meta_files)} in {metadata_dir}")
    else:
        info("Metadata folder: not found (Phase 1 may not have run)")

    # Count acquisition folders and miniseed files
    acq_folders = []
    total_miniseed = 0
    for item in sorted(os.listdir(output_path)):
        item_path = os.path.join(output_path, item)
        if os.path.isdir(item_path) and item not in ("metadata",):
            ms_count = sum(
                1 for f in os.listdir(item_path) if f.lower().endswith(MINISEED_EXT)
            )
            if ms_count > 0:
                acq_folders.append((item, ms_count))
                total_miniseed += ms_count

    # Also count any unsorted miniseed files still in the root
    root_ms = sum(
        1 for f in os.listdir(output_path) if f.lower().endswith(MINISEED_EXT)
    )
    if root_ms:
        warn(f"{root_ms} .miniseed file(s) remain unsorted in the output root.")

    if acq_folders:
        info(
            f"Acquisitions: {len(acq_folders)},  total .miniseed files: {total_miniseed}"
        )
        for folder, count in acq_folders:
            print(f"    {folder}/  — {count} file(s)")

    print()
    info("Done. Config saved to: " + str(CONFIG_PATH))


# ──────────────────────────────────────────────────────────────────────
# Main entry point
# ──────────────────────────────────────────────────────────────────────
def main() -> None:
    header("GEOPHONE DATA UNLOADING TOOL")
    print(
        textwrap.dedent(
            """\
        This script consolidates the geophone data unloading workflow:
          Phase 1 — Copy metadata (LOG files) from geophone drives
          Phase 2 — Rename & sort .miniseed files exported by Sololite

        Paths are saved to a config file for next time. Tab-completion is
        available when entering paths.
    """
        )
    )

    cfg = load_config()
    if cfg:
        info(f"Loaded previous config from {CONFIG_PATH}")

    # ── Geophones table (must be in the same directory as this script) ──
    script_dir = os.path.dirname(os.path.abspath(__file__))
    geo_table_candidates = [
        os.path.join(script_dir, "geophones_table.csv"),
        os.path.join(script_dir, "geophones_table"),
    ]
    geo_table_path = None
    for candidate in geo_table_candidates:
        if os.path.isfile(candidate):
            geo_table_path = candidate
            break
    if geo_table_path is None:
        error(
            f"geophones_table not found in {script_dir}\n"
            f"  Expected one of: {', '.join(os.path.basename(c) for c in geo_table_candidates)}"
        )
        sys.exit(1)

    geo_table = load_geophones_table(geo_table_path)

    # ── Phase 1: Metadata ──
    output_path = run_phase1(cfg, geo_table)

    # ── Phase 2: Seismic data ──
    run_phase2(cfg, geo_table, output_path)

    if output_path and os.path.isdir(output_path):
        print_summary(output_path)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n")
        info("Interrupted by user.")
        sys.exit(0)
    except EOFError:
        print("\n")
        info("Input closed.")
        sys.exit(0)
