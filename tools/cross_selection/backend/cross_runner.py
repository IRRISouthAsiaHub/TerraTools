from __future__ import annotations

from pathlib import Path

from tools.cross_selection.backend.ocs_optimizer import OCSError, run_ocs


BACKEND_DIR = Path(__file__).resolve().parent


def run_cross_selection(
    config: dict,
    backend_dir: Path | None = None,
    results_dir: Path | None = None,
) -> dict:

    if results_dir is None:

        results_dir = BACKEND_DIR / "results"

    results_dir = Path(results_dir)

    try:
        ocs = run_ocs(config, results_dir)
    except OCSError as e:
        raise RuntimeError(str(e)) from e

    return {
        "status": "success",
        "ocs": ocs,
    }


# ================================================================
# FAST PREVIEW
# Reads and aligns the files so the screen can report the number of
# genotypes and markers before the (heavier) optimisation is run.
# ================================================================

def preview_markers(geno_file, effects_file, map_file, map_kind=None, cm_per_mb=None) -> dict:
    """Counts shown before a run: dosage matrix + marker effects + genetic map (or VCF)."""
    from tools.cross_selection.backend import ocs_usefulness as U

    if not (geno_file and effects_file and map_file):
        return {"status": "error", "message": "Select the allele dosage matrix, the marker effects and the genetic map (or the VCF) to preview."}
    try:
        res = U.preview(geno_file, effects_file, map_file, map_kind, cm_per_mb)
    except U.UsefulnessError as e:
        return {"status": "error", "message": str(e)}
    res["status"] = "success"
    return res
