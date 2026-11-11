"""
balloon_runner.py — draws the favourable-allele balloon heatmap of the Favourable Allele Analysis screen.

The figure is made by balloon_heatmap.R, which is the R script supplied with the Heatmap_Final folder (data preparation,
quality control and balloon heatmap unchanged). This module only writes the genes and lines the user chose in the
two-header-row layout of that script's Excel sheet, runs Rscript and reads back what it wrote.
"""

import csv
import subprocess
import tempfile
import traceback
from pathlib import Path

from tools.common.r_locator import find_rscript, get_r_env, rscript_not_found_message
from tools.common.subprocess_utils import no_window_kwargs

SCRIPT = Path(__file__).resolve().parent / "balloon_heatmap.R"
META = ["Gene_ID", "Gene name", "Gene_label", "Broad trait category", "Sub category", "Marker_in_8k",
        "CHROM", "POS", "REF", "ALT", "Unfavourable allele", "Favourable allele"]
DPIS = (150, 300, 600)

_COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}


def _revcomp(s: str) -> str:
    return "".join(_COMP.get(c, "N") for c in reversed(s.upper()))


def _alleles(ref: str, alt: str, fav: str, unfav: str):
    """(favourable, unfavourable) written on the strand of REF and ALT, as in the sheet the script was made for.
    An allele that matches neither REF nor ALT (not even on the opposite strand) is passed on as written, so the
    script's own quality control drops that marker as contradicting."""
    for a, b in ((ref, alt), (alt, ref)):
        if a and fav and (fav == a or (fav != "-" and fav == _revcomp(a))):
            return a, b
    return fav, unfav


def _unique(names):
    seen, out = {}, []
    for n in names:
        n = (n or "").strip() or "unnamed"
        if n in seen:
            seen[n] += 1
            n = f"{n} ({seen[n]})"
        else:
            seen[n] = 1
        out.append(n)
    return out


def _read_csv(path: Path):
    if not path.is_file():
        return []
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def draw(payload: dict) -> dict:
    try:
        snps = payload.get("snps") or []
        lines = payload.get("varieties") or []
        calls = payload.get("calls") or []
        if not snps or not lines:
            return {"status": "error", "message": "Choose at least one gene and one line."}
        if len(calls) != len(lines) or any(len(c) != len(snps) for c in calls):
            return {"status": "error", "message": "The genotype table does not match the chosen genes and lines."}
        dpi = int(payload.get("dpi") or 300)
        dpi = dpi if dpi in DPIS else 300

        rscript = find_rscript()
        if not rscript:
            return {"status": "error", "message": rscript_not_found_message()}

        n_genes = len({s.get("name") or s.get("gene") for s in snps})
        # the size of the original figure (340 x 200 mm for 15 genes and 20 lines), grown or shrunk with the data
        width = min(900, max(170, 130 + 14 * n_genes))
        height = min(1200, max(110, 60 + 7 * len(lines)))

        work = Path(tempfile.mkdtemp(prefix="terratools_balloon_"))
        names = _unique([v.get("name") for v in lines])
        groups = ["Landrace" if v.get("group") == "Landrace" else "Variety" for v in lines]
        sheet = work / "input.csv"
        with open(sheet, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(META + groups)
            w.writerow([""] * len(META) + names)
            for i, s in enumerate(snps):
                ref, alt = (s.get("ref") or "").strip(), (s.get("alt") or "").strip()
                fav, unfav = _alleles(ref, alt, (s.get("fav") or "").strip(), (s.get("unfav") or "").strip())
                gname = s.get("name") or s.get("gene")
                w.writerow([s.get("gene"), gname, gname, s.get("trait") or "Unspecified", s.get("sub") or "",
                            s.get("marker"), s.get("chrom"), s.get("pos"), ref, alt, unfav, fav]
                           + [-1 if c[i] == "x" else c[i] for c in calls])

        out = work / "out"
        proc = subprocess.run([rscript, str(SCRIPT), str(sheet), str(out), str(width), str(height), str(dpi), str(payload.get("group_label") or "Lines")],
                              capture_output=True, text=True, encoding="utf-8", errors="replace",
                              env=get_r_env(rscript), cwd=str(work), **no_window_kwargs())
        log = ((proc.stdout or "") + "\n" + (proc.stderr or "")).strip()
        png = out / "gene_balloon.png"
        if proc.returncode != 0 or not png.is_file():
            tail = "\n".join(log.splitlines()[-12:])
            if "there is no package called" in log:
                tail = "An R package the heatmap needs is missing.\n" + tail
            else:
                tail = ("This usually means no marker is left after the quality control (a marker is dropped when it is missing in 80 % "
                        "or more of the chosen lines, or when its favourable allele matches neither REF nor ALT). "
                        "Choose more lines or other genes.\n" + tail)
            return {"status": "error", "message": f"The heatmap could not be drawn.\n{tail}"}

        qc = _read_csv(out / "qc_marker_filter_report.csv")
        dropped = [r for r in qc if r.get("Dropped") == "TRUE"]
        counts = _read_csv(out / "locus_counts.csv")
        return {
            "status": "success",
            "image": str(png),
            "outdir": str(out),
            "n_lines": len(lines),
            "n_genes": n_genes,
            "n_markers_in": len(qc),
            "n_markers_used": len(qc) - len(dropped),
            "dropped": [{"marker": r["Marker_in_8k"], "gene": r["Gene"], "missing": r["Missing_fraction"],
                         "reason": ("favourable allele matches neither REF nor ALT" if r["Contradicting"] == "TRUE"
                                    else "missing in 80 % or more of the lines")}
                        for r in dropped],
            "load": [{"line": r["Variety"], "favourable": int(r["Favourable"]), "total": int(r["Total"])} for r in counts],
            "files": [{"label": "Gene frequency matrix (CSV)", "path": str(out / "gene_freq_matrix.csv"), "name": "gene_freq_matrix.csv"},
                      {"label": "Marker QC report (CSV)", "path": str(out / "qc_marker_filter_report.csv"), "name": "qc_marker_filter_report.csv"},
                      {"label": "Favourable counts per line (CSV)", "path": str(out / "locus_counts.csv"), "name": "locus_counts.csv"}],
            "width_mm": width, "height_mm": height, "dpi": dpi,
        }
    except Exception as e:  # noqa: BLE001
        return {"status": "error", "message": f"{e}\n{traceback.format_exc()[-600:]}"}
