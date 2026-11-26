
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
scan_chipiron_experiments.py

CHIPIRON EXPERIMENT DATABASE
============================

Recursively scans Phonefleet_data for all CSV files containing
"input" in their filename.

Groups experiments by their actual test parameters, then by date.

ADDITIONAL DATA-INTEGRITY CHECK
-------------------------------
Compares the number of time-stamped experiment subfolders with
the number of time-stamped folders containing at least one
"input" CSV file.

If some time-stamped folders have no input CSV, a WARNING is
printed listing their full relative paths.

This check is independent of the experiment analysis.

OUTPUT
------
Phonefleet_data/
    data_summary/
        chipiron_exp_summary.csv

CSV ORGANIZATION
----------------
One row per TEST TYPE x DATE.

Columns:
    category
    test
    mean_rpm
    amplitude_rpm
    frequency_hz
    date
    nb_runs
    hours
    comments
    input_files

MO:
    Sines are grouped by mean RPM, amplitude RPM and frequency.
    Sweeps are grouped by their sweep parameters.

MC:
    Tests are grouped by their descriptive protocol notes.
    If no useful protocol note is available, the commanded PWM
    sequence is used as a fallback identifier.

OTHER:
    Tests are grouped by their available descriptive notes.

Each input file counts as one run.

The summary CSV is overwritten at every execution.

DEPENDENCIES
------------
Python standard library only.

USAGE
-----
python scan_chipiron_experiments.py "D:/icewave/chipiron/Phonefleet_data"
"""

import csv
import re
import sys

from pathlib import Path
from collections import defaultdict


# ================================================================
# CONFIGURATION
# ================================================================

SUMMARY_FOLDER = "data_summary"
SUMMARY_FILENAME = "chipiron_exp_summary.csv"

PARAMETER_PRECISION = 6

# A time-stamped folder must contain HH-MM-SS or HH:MM:SS.
# The same format is accepted in input filenames.
TIME_FOLDER_RE = re.compile(
    r"(?<!\d)(\d{2})[-:](\d{2})[-:](\d{2})(?!\d)"
)


# ================================================================
# REGULAR EXPRESSIONS
# ================================================================

NUMBER = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"

DATE_RE = re.compile(
    r"\b(20\d{2}-\d{2}-\d{2})\b"
)

TIME_RE = re.compile(
    r"\b(\d{2})[-:](\d{2})[-:](\d{2})\b"
)

MEAN_RE = re.compile(
    rf"\bmean\s*[:=]?\s*({NUMBER})\s*(?:rpm)?",
    re.IGNORECASE
)

AMP_RE = re.compile(
    rf"\bamp(?:litude)?\s*[:=]?\s*({NUMBER})\s*(?:rpm)?",
    re.IGNORECASE
)

FREQ_RE = re.compile(
    rf"\b({NUMBER})\s*Hz"
    rf"|\bfreq(?:uency)?\s*[:=]?\s*({NUMBER})",
    re.IGNORECASE
)

CYCLES_RE = re.compile(
    rf"\b({NUMBER})\s*cycles?"
    rf"|\bcycles?\s*[:=]?\s*({NUMBER})",
    re.IGNORECASE
)

O_COMMAND_RE = re.compile(
    rf"^o\s*({NUMBER})\s*,\s*"
    rf"({NUMBER})\s*,\s*"
    rf"({NUMBER})$",
    re.IGNORECASE
)

C_COMMAND_RE = re.compile(
    rf"^c\s*({NUMBER})$",
    re.IGNORECASE
)


# ================================================================
# GENERIC UTILITIES
# ================================================================

def safe_float(value):
    """Convert a value to float, or return None."""

    if value is None:
        return None

    try:
        return float(value)
    except (ValueError, TypeError):
        return None


def extract_first(regex, text):
    """Extract the first numerical capture group."""

    match = regex.search(text)

    if not match:
        return None

    for group in match.groups():

        if group is not None:
            return safe_float(group)

    return None


def format_number(value):
    """Format a number without unnecessary trailing zeros."""

    if value is None:
        return "unknown"

    return f"{value:g}"


def normalized_number(value):
    """Round a parameter for consistent grouping."""

    if value is None:
        return None

    return round(
        float(value),
        PARAMETER_PRECISION
    )


# ================================================================
# DATE AND TIME EXTRACTION
# ================================================================

def extract_date(path):
    """
    Extract date from filename or parent directories.

    Priority:
        1. Filename
        2. Closest parent directory
        3. Other parent directories
    """

    candidates = [path.name] + [
        p.name for p in path.parents
    ]

    for candidate in candidates:

        match = DATE_RE.search(candidate)

        if match:
            return match.group(1)

    return "Unknown"


def extract_time(path, rows=None):
    """
    Extract experiment time.

    Priority:
        1. Session folder / input filename.
        2. First datetime in the input CSV.

    Returns HH:MM:SS.
    """

    candidates = [path.stem] + [
        p.name for p in path.parents
    ]

    for candidate in candidates:

        match = TIME_RE.search(candidate)

        if match:

            return ":".join(match.groups())

    # Fallback to CSV timestamp.
    if rows:

        for row in rows:

            timestamp = str(
                row.get("datetime", "") or ""
            ).strip()

            match = TIME_RE.search(timestamp)

            if match:

                return ":".join(match.groups())

    return "Unknown"


# ================================================================
# NEW: TIME-STAMPED FOLDER / INPUT FILE CONSISTENCY CHECK
# ================================================================

def normalize_time_folder(path):
    """
    Return the HH:MM:SS time extracted from a folder name.

    Returns None if the folder does not contain a valid time.
    """

    match = TIME_FOLDER_RE.search(path.name)

    if not match:
        return None

    hour, minute, second = map(
        int,
        match.groups()
    )

    if hour > 23 or minute > 59 or second > 59:
        return None

    return f"{hour:02d}:{minute:02d}:{second:02d}"


def check_time_folders(root, input_files):
    """
    Compare time-stamped experiment folders with input CSV files.

    A time-stamped folder is considered represented if at least
    one input CSV is located inside it or one of its descendants.

    Each time-stamped folder is counted only once.

    The check uses folder paths rather than just timestamps,
    because two experiments can potentially have identical times
    in different parent directories.

    Prints warnings for:
        - Time-stamped folders without input CSV files.
        - Input CSV files not associated with any time-stamped
          folder.
        - A mismatch between the number of time-stamped folders
          and the number of input CSV files.
    """

    print("\n" + "=" * 80)
    print("TIME-FOLDER / INPUT-FILE CONSISTENCY CHECK")
    print("=" * 80)

    # ------------------------------------------------------------
    # 1. Find all time-stamped directories.
    # ------------------------------------------------------------

    time_folders = {}

    for path in root.rglob("*"):

        if not path.is_dir():
            continue

        relative_parts = path.relative_to(root).parts

        if SUMMARY_FOLDER in relative_parts:
            continue

        time_value = normalize_time_folder(path)

        if time_value is None:
            continue

        time_folders[path] = time_value

    # ------------------------------------------------------------
    # 2. Associate input files with their closest time folder.
    # ------------------------------------------------------------

    represented_folders = set()

    unassociated_files = []

    for input_file in input_files:

        matching_folders = [
            folder
            for folder in time_folders
            if folder in input_file.parents
        ]

        if matching_folders:

            # Use the closest (deepest) time-stamped folder.
            closest_folder = max(
                matching_folders,
                key=lambda p: len(p.parts)
            )

            represented_folders.add(
                closest_folder
            )

        else:

            unassociated_files.append(
                input_file
            )

    # ------------------------------------------------------------
    # 3. Identify missing folders.
    # ------------------------------------------------------------

    missing_folders = [
        folder
        for folder in time_folders
        if folder not in represented_folders
    ]

    nb_folders = len(time_folders)
    nb_input_files = len(input_files)
    nb_represented = len(represented_folders)

    print(
        f"\nTime-stamped folders found: {nb_folders}"
    )

    print(
        f"Input CSV files found:      {nb_input_files}"
    )

    print(
        f"Folders containing input:   {nb_represented}"
    )

    print(
        f"Folders without input:      {len(missing_folders)}"
    )

    print(
        f"Input files without folder: {len(unassociated_files)}"
    )

    # ------------------------------------------------------------
    # 4. Report discrepancies.
    # ------------------------------------------------------------

    discrepancy = (
        nb_folders != nb_input_files
        or bool(missing_folders)
        or bool(unassociated_files)
    )

    if discrepancy:

        print("\n" + "!" * 80)
        print("[WARNING] TIME-FOLDER / INPUT-FILE MISMATCH")
        print("!" * 80)

        print(
            f"\nThere are {nb_folders} time-stamped folders "
            f"but {nb_input_files} input CSV files."
        )

        # --------------------------------------------------------
        # List folders requiring recovery.
        # --------------------------------------------------------

        if missing_folders:

            print(
                "\n[TIME FOLDERS WITH NO INPUT CSV FILE]"
            )

            print(
                "These folders should be checked / recovered:"
            )

            for folder in sorted(missing_folders):

                relative_path = folder.relative_to(root)

                date = extract_date(folder)

                time_value = time_folders[folder]

                print(
                    f"  [MISSING INPUT] "
                    f"{date} | {time_value} | {relative_path}"
                )

        else:

            print(
                "\n[OK] Every time-stamped folder contains "
                "at least one input CSV."
            )

        # --------------------------------------------------------
        # List unassociated input files.
        # --------------------------------------------------------

        if unassociated_files:

            print(
                "\n[INPUT FILES NOT ASSOCIATED WITH "
                "A TIME-STAMPED FOLDER]"
            )

            for path in sorted(unassociated_files):

                print(
                    f"  [UNASSOCIATED] {path.relative_to(root)}"
                )

        print("\n" + "!" * 80)

    else:

        print(
            "\n[OK] Number of time-stamped folders matches "
            "number of input CSV files."
        )

        print(
            "[OK] Every time-stamped folder contains "
            "exactly one associated input CSV."
        )

    print("=" * 80 + "\n")

    return {
        "nb_time_folders": nb_folders,
        "nb_input_files": nb_input_files,
        "nb_represented_folders": nb_represented,
        "missing_folders": missing_folders,
        "unassociated_files": unassociated_files
    }


# ================================================================
# READ INPUT CSV
# ================================================================

def read_input_csv(path):
    """Read input CSV with encoding and delimiter detection."""

    encodings = [
        "utf-8-sig",
        "utf-8",
        "cp1252",
        "latin-1"
    ]

    for encoding in encodings:

        try:

            with open(
                path,
                "r",
                encoding=encoding,
                newline=""
            ) as f:

                sample = f.read(4096)

                f.seek(0)

                try:

                    dialect = csv.Sniffer().sniff(
                        sample,
                        delimiters=",;\t"
                    )

                except csv.Error:

                    dialect = csv.excel

                reader = csv.DictReader(
                    f,
                    dialect=dialect
                )

                if not reader.fieldnames:
                    return []

                reader.fieldnames = [
                    h.strip().lower()
                    for h in reader.fieldnames
                ]

                return list(reader)

        except UnicodeDecodeError:
            continue

        except Exception as exc:

            print(
                f"[WARNING] Cannot read {path}: {exc}"
            )

            return []

    print(
        f"[WARNING] Unknown encoding: {path}"
    )

    return []


# ================================================================
# COMMAND EXTRACTION
# ================================================================

def get_commands(rows):
    """Return all non-empty commands in chronological order."""

    return [

        str(row.get("command", "") or "").strip()

        for row in rows

        if str(row.get("command", "") or "").strip()

    ]


def get_mode(commands):
    """Identify MC or MO mode."""

    for command in commands:

        cmd = command.lower()

        if cmd == "mc":
            return "mc"

        if cmd == "mo":
            return "mo"

    return "unknown"


def extract_o_commands(commands):
    """
    Extract oscillation commands.

    Format:
        o<mean_duty>,<amplitude_duty>,<frequency>

    Example:
        o50.00,7.58,0.5000

    NOTE:
        The first two values are duty parameters, not RPM.
    """

    results = []

    for command in commands:

        match = O_COMMAND_RE.match(command)

        if match:

            results.append({

                "mean_duty": float(match.group(1)),

                "amp_duty": float(match.group(2)),

                "frequency": float(match.group(3))

            })

    return results


def extract_c_commands(commands):
    """Extract constant PWM commands."""

    values = []

    for command in commands:

        match = C_COMMAND_RE.match(command)

        if match:

            values.append(
                float(match.group(1))
            )

    return values


def compress_sequence(values):
    """
    Remove consecutive duplicates from a command sequence.

    Example:
        [10, 10, 20, 20, 30, 20, 20]
        becomes:
        [10, 20, 30, 20]
    """

    result = []

    for value in values:

        if not result or value != result[-1]:

            result.append(value)

    return result


# ================================================================
# NOTE EXTRACTION
# ================================================================

def extract_notes(rows):
    """Collect unique, non-empty notes in original order."""

    notes = []

    for row in rows:

        note = str(
            row.get("note", "") or ""
        ).strip()

        if note and note.lower() not in [
            "nan", "none"
        ]:

            if note not in notes:

                notes.append(note)

    return notes


def extract_physical_parameters(notes):
    """
    Extract physical parameters from human-readable notes.

    Example:
        sine mean 0 rpm, amp 500 rpm, 0.5 Hz, 10 cycles
    """

    text = " | ".join(notes)

    return {

        "mean_rpm": extract_first(
            MEAN_RE, text
        ),

        "amp_rpm": extract_first(
            AMP_RE, text
        ),

        "frequency_hz": extract_first(
            FREQ_RE, text
        ),

        "cycles": extract_first(
            CYCLES_RE, text
        )

    }


# ================================================================
# NOTE CLEANING / PROTOCOL IDENTIFICATION
# ================================================================

def is_boilerplate_note(note):
    """
    Identify notes that do not describe the actual experiment.

    These notes are retained in the original input files but
    excluded from test identification.
    """

    text = note.lower().strip()

    boilerplate = [

        "arduino mode",

        "rpm before start",

        "before start",

        "rest (end)",

        "rest after experiment",

        "rest before experiment",

        "acquisition",

        "starting acquisition",

        "stopping acquisition",

        "start acquisition",

        "stop acquisition"

    ]

    return any(
        phrase in text
        for phrase in boilerplate
    )


def get_protocol_notes(notes):
    """
    Return meaningful notes for identifying a test protocol.

    Excludes generic setup/rest messages.
    """

    return [

        note

        for note in notes

        if not is_boilerplate_note(note)

    ]


def normalize_description(text):
    """
    Normalize a description for grouping.

    Case-insensitive, removes redundant spaces and trailing
    punctuation.
    """

    text = text.lower().strip()

    text = re.sub(
        r"\s+",
        " ",
        text
    )

    text = re.sub(
        r"[.,;]+$",
        "",
        text
    )

    return text


def get_mc_description(protocol_notes):
    """
    Identify an MC test from its descriptive notes.

    Prefer notes explicitly describing a plateau/ramp protocol.

    Otherwise use the first meaningful note.
    """

    if not protocol_notes:
        return None

    preferred_keywords = [

        "plateau",

        "plateaux",

        "ramp",

        "ramps",

        "step",

        "steps",

        "rpm",

        "speed",

        "velocity",

        "test"

    ]

    for note in protocol_notes:

        lower = note.lower()

        if any(
            keyword in lower
            for keyword in preferred_keywords
        ):

            return note.strip()

    return protocol_notes[0].strip()


# ================================================================
# TEST CLASSIFICATION
# ================================================================

def is_sweep_test(notes, o_commands):
    """Identify an oscillation sweep."""

    text = " ".join(notes).lower()

    keywords = [

        "sweep",

        "frequency sweep",

        "amplitude sweep",

        "frequency scan",

        "parameter sweep",

        "multi-frequency"

    ]

    return (

        any(keyword in text for keyword in keywords)

        or len(o_commands) > 1

    )


def classify_experiment(
    mode,
    notes,
    o_commands,
    c_commands
):
    """Classify an experiment."""

    if mode == "mc":

        if c_commands:

            return "MC plateaux-ramps"

        return "MC other"

    if mode == "mo":

        if is_sweep_test(notes, o_commands):

            return "MO oscillation sweep"

        if o_commands:

            return "MO oscillations"

        return "MO other"

    # Fallback when mode command is missing.

    if o_commands:

        if is_sweep_test(notes, o_commands):

            return "MO oscillation sweep"

        return "MO oscillations"

    if c_commands:

        return "MC plateaux-ramps"

    return "Other"


# ================================================================
# IDENTIFY THE ACTUAL TEST
# ================================================================

def identify_test(
    category,
    physical,
    o_commands,
    c_commands,
    protocol_notes
):
    """
    Return:
        test label
        grouping signature

    The grouping signature determines which experiments
    belong to the same test.
    """

    mean_rpm = physical["mean_rpm"]

    amp_rpm = physical["amp_rpm"]

    freq = physical["frequency_hz"]

    # ------------------------------------------------------------
    # MO: SINGLE SINE
    # ------------------------------------------------------------

    if category == "MO oscillations":

        # Frequency can also be extracted from the o command.
        if freq is None and o_commands:

            freq = o_commands[0]["frequency"]

        mean_key = normalized_number(mean_rpm)

        amp_key = normalized_number(amp_rpm)

        freq_key = normalized_number(freq)

        test_label = (
            f"Sine: "
            f"mean={format_number(mean_rpm)} rpm, "
            f"A={format_number(amp_rpm)} rpm, "
            f"f={format_number(freq)} Hz"
        )

        signature = (
            "sine",
            mean_key,
            amp_key,
            freq_key
        )

        return test_label, signature

    # ------------------------------------------------------------
    # MO: OSCILLATION SWEEP
    # ------------------------------------------------------------

    if category == "MO oscillation sweep":

        if freq is None and len(o_commands) == 1:

            freq = o_commands[0]["frequency"]

        frequencies = sorted(set(

            normalized_number(item["frequency"])

            for item in o_commands

        ))

        if len(frequencies) > 1:

            frequency_label = ", ".join(

                format_number(v)

                for v in frequencies

            ) + " Hz"

        elif freq is not None:

            frequency_label = (
                format_number(freq) + " Hz"
            )

        else:

            frequency_label = "unknown"

        mean_key = normalized_number(mean_rpm)

        amp_key = normalized_number(amp_rpm)

        # Use descriptive notes to distinguish sweep protocols.
        description = get_mc_description(
            protocol_notes
        )

        test_label = (
            f"Oscillation sweep: "
            f"mean={format_number(mean_rpm)} rpm, "
            f"A={format_number(amp_rpm)} rpm, "
            f"frequencies={frequency_label}"
        )

        signature = (
            "sweep",
            mean_key,
            amp_key,
            tuple(frequencies),
            normalize_description(description or "")
        )

        return test_label, signature

    # ------------------------------------------------------------
    # MC: PLATEAUX AND RAMPS
    # ------------------------------------------------------------

    if category == "MC plateaux-ramps":

        description = get_mc_description(
            protocol_notes
        )

        duty_sequence = compress_sequence(
            c_commands
        )

        if description:

            # Descriptive note is the primary identifier.
            # PWM sequence is retained as a fallback / verification
            # detail, but does not split tests with the same note.
            test_label = description

            signature = (
                "mc_description",
                normalize_description(description)
            )

        elif duty_sequence:

            # Fallback if no descriptive note exists.
            duty_text = " -> ".join(
                format_number(v)
                for v in duty_sequence
            )

            test_label = (
                f"MC PWM sequence: {duty_text} %"
            )

            signature = (
                "mc_pwm",
                tuple(
                    normalized_number(v)
                    for v in duty_sequence
                )
            )

        else:

            test_label = "MC test (protocol unknown)"

            signature = (
                "mc_unknown",
            )

        return test_label, signature

    # ------------------------------------------------------------
    # OTHER TESTS
    # ------------------------------------------------------------

    description = (
        " | ".join(protocol_notes)
        if protocol_notes
        else "Protocol unknown"
    )

    test_label = description

    signature = (
        "other",
        normalize_description(description)
    )

    return test_label, signature


# ================================================================
# ANALYSE ONE INPUT FILE
# ================================================================

def analyse_file(path, root):

    rows = read_input_csv(path)

    if not rows:
        return None

    commands = get_commands(rows)

    notes = extract_notes(rows)

    protocol_notes = get_protocol_notes(notes)

    mode = get_mode(commands)

    o_commands = extract_o_commands(commands)

    c_commands = extract_c_commands(commands)

    physical = extract_physical_parameters(notes)

    category = classify_experiment(
        mode,
        notes,
        o_commands,
        c_commands
    )

    test_label, signature = identify_test(
        category,
        physical,
        o_commands,
        c_commands,
        protocol_notes
    )

    # Extract frequency from command if not present in notes.
    frequency = physical["frequency_hz"]

    if frequency is None and o_commands:

        if len(o_commands) == 1:

            frequency = o_commands[0]["frequency"]

    # Time from session folder / filename, with CSV fallback.
    experiment_time = extract_time(
        path,
        rows
    )

    # Keep all notes as comments.
    comments = list(notes)

    # Explicitly add cycles if parsed but not already apparent
    # in the original notes.
    cycles = physical["cycles"]

    if cycles is not None:

        cycle_text = (
            f"{format_number(cycles)} cycles"
        )

        if not any(
            "cycle" in note.lower()
            for note in comments
        ):

            comments.append(cycle_text)

    return {

        "date": extract_date(path),

        "time": experiment_time,

        "category": category,

        "test": test_label,

        "signature": signature,

        "mean_rpm": physical["mean_rpm"],

        "amplitude_rpm": physical["amp_rpm"],

        "frequency_hz": frequency,

        "cycles": cycles,

        "comments": comments,

        "relative_path": str(
            path.relative_to(root)
        )

    }


# ================================================================
# SCAN ALL INPUT FILES
# ================================================================

def scan_root(root):
    """
    Recursively find all CSV files containing 'input'.

    Excludes data_summary to avoid scanning generated outputs.
    """

    files = sorted(

        path for path in root.rglob("*")

        if path.is_file()

        and "input" in path.name.lower()

        and path.suffix.lower() == ".csv"

        and SUMMARY_FOLDER not in path.relative_to(root).parts

    )

    print("\n" + "=" * 80)

    print("CHIPIRON EXPERIMENT SCANNER")

    print("=" * 80)

    print(f"\nRoot folder: {root}")

    print(f"Input files found: {len(files)}")

    # NEW: Check folder/input consistency before analysing files.
    check_time_folders(
        root,
        files
    )

    experiments = []

    for i, path in enumerate(files, start=1):

        result = analyse_file(
            path,
            root
        )

        if result is not None:

            experiments.append(result)

        else:

            print(
                f"[WARNING] Empty/unreadable: {path}"
            )

        if i % 50 == 0:

            print(
                f"Scanned {i}/{len(files)} files..."
            )

    experiments.sort(
        key=lambda x: (
            x["category"],
            str(x["signature"]),
            x["date"],
            x["time"]
        )
    )

    print(
        f"\nSuccessfully analysed "
        f"{len(experiments)} input files."
    )

    return experiments


# ================================================================
# GROUP BY TEST TYPE AND DATE
# ================================================================

def build_test_summary(experiments):
    """
    Group experiments by:

        1. Test signature
        2. Date

    Each input file counts as one run.

    Hours are collected from the session folder / filename.

    Comments are combined without duplicates.
    """

    groups = defaultdict(list)

    for exp in experiments:

        key = (

            exp["signature"],

            exp["date"]

        )

        groups[key].append(exp)

    summary_rows = []

    for (signature, date), items in groups.items():

        first = items[0]

        # All session times for this test on this date.
        hours = [

            item["time"]

            for item in items

            if item["time"] != "Unknown"

        ]

        # Keep unique times in chronological order.
        hours = sorted(set(hours))

        # Collect all comments.
        comments = []

        for item in items:

            for comment in item["comments"]:

                if comment not in comments:

                    comments.append(comment)

        # Keep all input file paths for traceability.
        input_files = [

            item["relative_path"]

            for item in items

        ]

        summary_rows.append({

            "category": first["category"],

            "test": first["test"],

            "mean_rpm": (
                first["mean_rpm"]
                if first["mean_rpm"] is not None
                else ""
            ),

            "amplitude_rpm": (
                first["amplitude_rpm"]
                if first["amplitude_rpm"] is not None
                else ""
            ),

            "frequency_hz": (
                first["frequency_hz"]
                if first["frequency_hz"] is not None
                else ""
            ),

            "date": date,

            "nb_runs": len(items),

            "hours": "; ".join(hours),

            "comments": " || ".join(comments),

            "input_files": " || ".join(input_files)

        })

    # Sort by category, test and date.
    summary_rows.sort(
        key=lambda row: (
            row["category"],
            row["test"],
            row["date"]
        )
    )

    return summary_rows


# ================================================================
# SAVE CSV
# ================================================================

def save_summary(root, summary_rows):
    """
    Create data_summary if it does not exist.

    Save / overwrite chipiron_exp_summary.csv.
    """

    summary_dir = root / SUMMARY_FOLDER

    summary_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    output_path = (
        summary_dir / SUMMARY_FILENAME
    )

    fieldnames = [

        "category",

        "test",

        "mean_rpm",

        "amplitude_rpm",

        "frequency_hz",

        "date",

        "nb_runs",

        "hours",

        "comments",

        "input_files"

    ]

    with open(
        output_path,
        "w",
        newline="",
        encoding="utf-8-sig"
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames
        )

        writer.writeheader()

        writer.writerows(summary_rows)

    print("\n" + "=" * 80)

    print("[SAVED] CHIPIRON EXPERIMENT SUMMARY")

    print(output_path)

    print("=" * 80)

    return output_path


# ================================================================
# DISPLAY SUMMARY IN TERMINAL
# ================================================================

def print_summary(summary_rows, experiments):

    print("\n")

    print("=" * 90)

    print("CHIPIRON TEST CATALOGUE")

    print("=" * 90)

    if not summary_rows:

        print("No experiments found.")

        return

    current_test = None

    for row in summary_rows:

        test_key = (
            row["category"],
            row["test"]
        )

        if test_key != current_test:

            current_test = test_key

            print()

            print(
                f"[{row['category']}] {row['test']}"
            )

            print("-" * 85)

        print(
            f"    {row['date']} | "
            f"{row['nb_runs']} run(s) | "
            f"hours: {row['hours'] or 'Unknown'}"
        )

        if row["comments"]:

            print(
                f"        Comments: {row['comments']}"
            )

    # ------------------------------------------------------------
    # GLOBAL STATISTICS
    # ------------------------------------------------------------

    print("\n")

    print("=" * 90)

    print("GLOBAL STATISTICS")

    print("=" * 90)

    print(
        f"Total input files / runs: {len(experiments)}"
    )

    print(
        f"Number of experiment days: "
        f"{len(set(e['date'] for e in experiments))}"
    )

    print(
        f"Number of distinct test/date combinations: "
        f"{len(summary_rows)}"
    )

    test_signatures = set(
        exp["signature"]
        for exp in experiments
    )

    print(
        f"Number of distinct test protocols: "
        f"{len(test_signatures)}"
    )

    categories = defaultdict(int)

    for exp in experiments:

        categories[exp["category"]] += 1

    print("\nRuns by category:")

    for category, count in sorted(categories.items()):

        print(
            f"    {category:<30}: {count}"
        )

    print("=" * 90)


# ================================================================
# MAIN
# ================================================================

def main():

    if len(sys.argv) < 2:

        print(
            "\nUsage:\n"
            'python scan_chipiron_experiments.py "PATH_TO_Phonefleet_data"'
        )

        sys.exit(1)

    root = Path(
        sys.argv[1]
    ).expanduser()

    if not root.exists():

        print(
            f"[ERROR] Folder does not exist: {root}"
        )

        sys.exit(1)

    if not root.is_dir():

        print(
            f"[ERROR] Not a directory: {root}"
        )

        sys.exit(1)

    # 1. Scan all input files and check missing time folders.
    experiments = scan_root(root)

    # 2. Group by test type and date.
    summary_rows = build_test_summary(
        experiments
    )

    # 3. Save the single summary CSV.
    output_path = save_summary(
        root,
        summary_rows
    )

    # 4. Display the results.
    print_summary(
        summary_rows,
        experiments
    )

    print("\nSCAN COMPLETED")

    print(f"Summary: {output_path}")


if __name__ == "__main__":

    main()
