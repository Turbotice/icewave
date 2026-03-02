# Geophone Data Unloading Tool

1. Open SoloLite software AS ADMINISTRATOR
2. Click `Open Project` and open the corresponding project
3. Connect first batch of geophones
4. Launch `python unload_smartsolo_geophones.py` in a terminal
5. Follow the steps to copy metadata (set path, press enter to copy metadata and repeat for next set of geophones)
6. Type `done` to pass to the seismic data phase 
7. In SoloLite, click the blue circle "Export seismic data"
Set: 
- A job name
- Output Data Component `All Components in separate file`
- Output Data Format `MiniSEED`
- Sample Interval: `1ms`
- Start and End time corresponding to your acquisition date
8. Click `Prepare`, make sure everything is as expected (No. of acquisitions, geophones, etc.) then click `Run`
9. When finished, find the output folder (probably something like `C:\SOLODATA\project_name\job_name`) and copy the path to the terminal 
10. The script will move the files to the path that you first set and rename them appropriately. Then move them onto the server for safekeeping.



**Script:** `unload_smartsolo_geophones.py`

A unified, interactive script that replaces the old fragmented workflow (`Move_data.py`, `prepare_smartsolo.py`, and their many variants) with a single two-phase pipeline for unloading SmartSolo geophone data after a day of field work.

---

## Quick Start

```bash
pip install prompt_toolkit pyyaml
python unload_smartsolo_geophones.py
```

The script is fully interactive — it will ask for every path it needs, with **tab-completion** and **pre-filled defaults** from the last run.

The **geophones table** file (`geophones_table.csv` or `geophones_table`) must be in the **same directory** as the script.

---

## What It Does

### Phase 1 — Metadata Extraction

Copies **DigiSolo / LOG files** from the geophone SD cards to a central folder, renaming them by geophone number.

| Step | What happens |
|------|-------------|
| 1 | You specify an **output folder** (created if it doesn't exist). |
| 2 | The script reads the **geophone drives** from `~/geophone_unload_config.yaml` (edit the YAML to change them). |
| 3 | The script scans each drive for DigiSolo/LOG files, parses the serial number, and copies them to `<output>/metadata/DigiSolo_<geophone_number>.txt`. |
| 4 | A **warning** is shown if fewer than 3 geophones with data are found across the configured drives. |
| 5 | You're prompted to plug in the next batch of 3 geophones and press Enter, or type `done` when all geophones are processed. |

Every file operation is printed to the console:
```
→ COPY  H:\DigiSolo.txt
       → C:\Canada_2026\metadata\DigiSolo_01.txt
       (serial 453024765 → geophone 01)
```

### Phase 2 — Seismic Data Organisation

Renames and sorts the `.miniseed` files that were exported by the Sololite software.

| Sub-step | What happens |
|----------|-------------|
| **2a — Rename** | Replaces the 9-digit serial number prefix with the 2-digit geophone number (e.g. `453024765.001.…` → `01.001.…`). |
| **2b — Fix filenames** | Cleans up common issues: double dots (`..` → `.`) and single-digit acquisition numbers (`01.1.…` → `01.001.…`). |
| **2c — Sort into folders** | Moves each file into a sub-folder named by acquisition number (e.g. `0001/`, `0002/`). |

You can choose to copy files to the output folder first (safe — Sololite export untouched) or work in-place.

---

## Workflow Diagram

```
  ┌─────────────────────┐
  │  Plug in 3 geophones │
  │  (mounted as drives) │
  └────────┬────────────┘
           │
           ▼
  ┌─────────────────────────┐
  │  Phase 1: Metadata      │◄─── repeats for each batch of 3
  │  • scan drives          │
  │  • parse DigiSolo files │
  │  • copy + rename        │
  └────────┬────────────────┘
           │  type 'done'
           ▼
  ┌──────────────────────────┐
  │  Export via Sololite     │  (separate vendor software)
  └────────┬─────────────────┘
           │
           ▼
  ┌──────────────────────────────────────┐
  │  Phase 2: Seismic data               │
  │  2a. Rename (serial → num)           │
  │  2b. Fix malformed filenames         │
  │  2c. Sort into acquisition folders   │
  └────────┬─────────────────────────────┘
           │
           ▼
  ┌──────────────────┐
  │  Summary report  │
  └──────────────────┘
```

---

## Config File

All paths you enter are saved to `~/geophone_unload_config.yaml`. Next time you run the script, they're pre-filled as defaults — just press Enter to accept.

Example config:
```yaml
all_data_output_path: C:\Canada_2026\0213\Geophones
geophone_drives:
- H:\
- G:\
- D:\
sololite_data_export_path: C:\SOLODATA\Bicwin2024_12dB\BicWin20260213
```

The **`geophone_drives`** list is used directly — the script does not prompt for drive letters. Edit this file to change them.

---

## Geophones Table

A tab-separated file placed **in the same directory as the script**, named `geophones_table.csv` (or `geophones_table`). The script auto-detects it — no path prompt needed. If the file is missing, the script exits with an error.

Format — header line, then 2-digit geophone IDs mapped to 9-digit serial numbers:

```
Num	geo_SN
01	453024765
02	453024177
03	453025132
...
16	453025291
```

The script skips the header automatically.

---

## Output Folder Structure

After both phases, the output folder looks like:

```
<output_path>/
├── metadata/
│   ├── DigiSolo_01.txt
│   ├── DigiSolo_02.txt
│   └── ...
├── 0001/
│   ├── 01.0001.2024.02.11.14.30.00.000.E.miniseed
│   ├── 01.0001.2024.02.11.14.30.00.000.N.miniseed
│   ├── 01.0001.2024.02.11.14.30.00.000.Z.miniseed
│   ├── 02.0001.…
│   └── ...
├── 0002/
│   └── ...
└── ...
```

- `metadata/` — DigiSolo/LOG files (one per geophone)
- `0001/`, `0002/`, … — one folder per acquisition, containing all the `.miniseed` files for that recording

---

## Menu Options

When launched, you can choose:

| Option | Description |
|--------|-------------|
| **1** | Run both phases (full workflow) |
| **2** | Phase 1 only (metadata extraction — useful when still plugging in geophones) |
| **3** | Phase 2 only (rename + sort — useful when metadata was already done) |

---

## Dependencies

| Package | Purpose | Install |
|---------|---------|---------|
| `prompt_toolkit` | Cross-platform tab-completion for path input (works on Windows) | `pip install prompt_toolkit` |
| `pyyaml` | Persist config between runs | `pip install pyyaml` |

Both are pure-Python with no system dependencies.

---

## What This Replaces

This script consolidates everything from:

| Old script | What it did | Status |
|------------|-------------|--------|
| `Move_data.py` | Copy DigiSolo metadata from SD cards | 3 nearly-identical variants with hardcoded paths and bugs |
| `Move_data (2).py`, `Move_data (3).py` | Same thing, different computers | Redundant |
| `prepare_smartsolo.py` | Rename miniseed files + sort into folders | 4 variants, one with unresolved git merge conflict |
| `prepare_smartsolo_seb.py`, `prepare_smartsolo_old.py`, `prepare_smartsolo (2).py` | Same thing, different authors | Redundant |

**Key fixes over the old code:**
- No hardcoded paths — everything is prompted with tab-completion
- Paths are persisted in a YAML config for next time
- Serial number parsing uses regex instead of brittle character-position slicing
- Fresh metadata dict per geophone (old code reused a global dict, leaking stale values)
- Every file operation is logged transparently
- Warns before overwriting existing files
- Handles missing/inaccessible drives gracefully
