"""
quick_data_visu.py -- scan Phonefleet_data/<date>/<time>/figures for the three overview plots of each acquisition
and build a Beamer deck (one slide per acquisition, three panels side by side).

Panels picked (linear scales only, nothing else in the figures folders is touched):
    <date>_<time>_usb_A1_A1_V_vs_t_s_linlin.png     A1 (V) vs time
    <date>_<time>_usb_A2_A2_V_vs_t_s_linlin.png     A2 (V) vs time
    <date>_<time>_usb_A2_A2_V_vs_A1_V_linlin.png    A2 (V) vs A1 (V)
(a trailing _2, _3... added by Chipiron_plotting_v5 when a file already existed is accepted; the newest file wins)

OUTPUT (in OUT_FOLDER, NOT inside ROOT_FOLDER, which the plotting script reads as a list of dates)
    slides_overview/slides.tex     compile with pdflatex (Overleaf: upload the zip, set main file to slides.tex)
    slides_overview/figs/          copies of the figures, flat names <date>_<time>_<panel>.png
    slides_overview.zip            the two above, ready to upload to Overleaf

Run:  py make_slides.py
"""
import re
import shutil
from pathlib import Path

# ======================= USER CONFIGURATION =======================
ROOT_FOLDER = Path("C:/Users/Gulliver/Desktop/Chipiron/git_clone_turbotice_icewave/icewave/icewave/chipiron/Phonefleet_data")
OUT_FOLDER = Path("slides_overview")
EXCLUDE_DATES = {"2026-11-08", "2026-10-22"}     # different experiments, not in the comparison
FIG_SUBFOLDER = "figures"
PANEL_WIDTH = 0.325                               # fraction of the text width for each panel
# ==================================================================

# (short key used in the copied file name, core of the original file name, caption)
PANELS = [
    ("A1_V_vs_t", "A1_A1_V_vs_t_s", "A1 (V) vs time"),
    ("A2_V_vs_t", "A2_A2_V_vs_t_s", "A2 (V) vs time"),
    ("A2_vs_A1", "A2_A2_V_vs_A1_V", "A2 (V) vs A1 (V)"),
]
STAMP = r"(?P<stamp>\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})"


def tex_escape(s):
    for a, b in (("\\", r"\textbackslash{}"), ("_", r"\_"), ("%", r"\%"), ("&", r"\&"), ("#", r"\#")):
        s = s.replace(a, b)
    return s


def find_panel(fig_dir, core):
    """Newest PNG of fig_dir matching the exact figure name for this panel, or None."""
    pat = re.compile(rf"^{STAMP}_usb_{re.escape(core)}_linlin(?:_\d+)?\.png$", re.I)
    hits = [p for p in fig_dir.glob("*.png") if pat.match(p.name)]
    if len(hits) > 1:
        hits.sort(key=lambda p: p.stat().st_mtime)
        print(f"  note: {len(hits)} files for '{core}' in {fig_dir}: using the newest ({hits[-1].name})")
    return hits[-1] if hits else None


def first_line(path):
    try:
        return path.read_text(encoding="utf-8").splitlines()[0].strip()
    except (OSError, IndexError):
        return None


def main():
    if not ROOT_FOLDER.is_dir():
        raise SystemExit(f"Root folder not found: {ROOT_FOLDER} (edit ROOT_FOLDER at the top of the script)")
    figs_out = OUT_FOLDER / "figs"
    if OUT_FOLDER.exists():
        shutil.rmtree(OUT_FOLDER)                    # regenerated at every run
    figs_out.mkdir(parents=True)

    frames, n_skipped = [], 0
    for date_dir in sorted(p for p in ROOT_FOLDER.iterdir() if p.is_dir()):
        if date_dir.name in EXCLUDE_DATES:
            print(f"{date_dir.name}: excluded")
            continue
        for time_dir in sorted(p for p in date_dir.iterdir() if p.is_dir()):
            fig_dir = time_dir / FIG_SUBFOLDER
            if not fig_dir.is_dir():
                n_skipped += 1
                continue
            found = [find_panel(fig_dir, core) for _, core, _ in PANELS]
            if not any(found):
                n_skipped += 1
                continue
            label = f"{date_dir.name} {time_dir.name}"
            tag = f"{date_dir.name}_{time_dir.name}"
            boxes = []
            for (key, _, caption), src in zip(PANELS, found):
                if src is None:
                    print(f"{label}: MISSING {caption}")
                    boxes.append((None, caption))
                else:
                    dest = figs_out / f"{tag}_{key}.png"
                    shutil.copy2(src, dest)
                    boxes.append((f"figs/{dest.name}", caption))
            notes = []
            edc = first_line(time_dir / "edc_used.txt")
            if edc:
                notes.append(f"ESCON: {edc}")
            per = first_line(time_dir / "sampling_period_ms.txt")
            if per:
                notes.append(f"sampling period: {per} ms")
            frames.append((label, boxes, " | ".join(notes)))
            print(f"{label}: {sum(f is not None for f in found)}/3 panels")

    if not frames:
        raise SystemExit("No acquisition with the expected figures found.")

    L = [r"\documentclass[aspectratio=169,10pt]{beamer}",
         r"\usepackage{graphicx}",
         r"\setbeamertemplate{navigation symbols}{}",
         r"\setbeamersize{text margin left=0.3cm,text margin right=0.3cm}",
         r"\setbeamerfont{frametitle}{size=\normalsize}",
         r"\begin{document}", ""]
    for label, boxes, note in frames:
        L.append(rf"\begin{{frame}}[t]{{{tex_escape(label)}}}")
        L.append(r"\centering\vspace{0.3cm}")
        for i, (path, caption) in enumerate(boxes):
            L.append(rf"\begin{{minipage}}[t]{{{PANEL_WIDTH}\textwidth}}\centering")
            if path:
                L.append(rf"\includegraphics[width=\linewidth]{{{path}}}")
            else:
                L.append(rf"\fbox{{\parbox[c][3cm][c]{{0.9\linewidth}}{{\centering\scriptsize missing}}}}")
            L.append(rf"\\ {{\scriptsize {tex_escape(caption)}}}")
            L.append(r"\end{minipage}" + (r"\hfill" if i < len(boxes) - 1 else ""))
        if note:
            L.append(rf"\vfill{{\tiny {tex_escape(note)}}}")
        L += [r"\end{frame}", ""]
    L.append(r"\end{document}")
    (OUT_FOLDER / "slides.tex").write_text("\n".join(L), encoding="utf-8")

    zip_path = shutil.make_archive(str(OUT_FOLDER), "zip", root_dir=OUT_FOLDER)
    print(f"\n{len(frames)} slide(s), {n_skipped} acquisition folder(s) without these figures.")
    print(f"Written: {(OUT_FOLDER / 'slides.tex').resolve()}")
    print(f"Overleaf zip: {Path(zip_path).resolve()}")


if __name__ == "__main__":
    main()