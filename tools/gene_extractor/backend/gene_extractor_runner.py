import json
import shutil
import subprocess
from pathlib import Path

from tools.common.r_locator import find_rscript, rscript_not_found_message
from tools.common.subprocess_utils import no_window_kwargs


# ============================================================
# PATHS
# ============================================================

BACKEND_DIR = Path(__file__).resolve().parent
R_SCRIPT = BACKEND_DIR / "gene_extractor.R"


def _error(message):
    return {"status": "error", "message": message}


# ============================================================
# PUBLIC ENTRYPOINT
#
# Mirrors the shape of gene_runner.run_gene(config, results_dir) from
# the gene mapping tool, so the Qt bridge/worker pattern is identical.
# ============================================================

def run_extraction(config, results_dir):
    """
    config: {
        "extracted_markers_path": str,
        "gene_markers_path": str,
        "rice_gene_path": str
    }

    results_dir: Path or str, directory the R script will write
                 Markers_with_physical_position.xlsx and
                 Genes_in_gene_regions.xlsx into.

    Returns a dict:
        {
          "status": "success",
          "markers": [ {Markers, Chromosome, Position}, ... ],
          "genes":   [ {Chromosome, Start, Stop, Strand, Gene_id,
                         Description, gene_start, gene_stop}, ... ],
          "markers_file": "<path>",
          "genes_file": "<path>"
        }
    or:
        { "status": "error", "message": "..." }
    """

    extracted_markers_path = config.get("extracted_markers_path")
    gene_markers_path = config.get("gene_markers_path")
    rice_gene_path = config.get("rice_gene_path")

    # Same validation the original Shiny observeEvent(input$extract, ...) did
    if not extracted_markers_path or not gene_markers_path or not rice_gene_path:
        return _error("Please upload all three required files.")

    for label, path in (
        ("Extracted markers", extracted_markers_path),
        ("gene markers", gene_markers_path),
        ("Rice gene list", rice_gene_path),
    ):
        if not Path(path).exists():
            return _error(f"{label} file not found: {path}")

    rscript = find_rscript()
    if not rscript:
        return _error(rscript_not_found_message())

    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        rscript,
        str(R_SCRIPT),
        str(extracted_markers_path),
        str(gene_markers_path),
        str(rice_gene_path),
        str(results_dir),
    ]

    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
            **no_window_kwargs(),
        )
    except subprocess.TimeoutExpired:
        return _error("gene extraction timed out.")
    except Exception as e:
        return _error(f"Failed to launch R backend: {e}")

    stdout = (proc.stdout or "").strip()
    stderr = (proc.stderr or "").strip()

    if not stdout:
        return _error(stderr or "R script produced no output.")

    try:
        result = json.loads(stdout)
    except json.JSONDecodeError:
        # If R printed warnings/messages before the JSON, try to recover
        # the JSON object from the last line of stdout.
        last_line = stdout.splitlines()[-1] if stdout.splitlines() else ""
        try:
            result = json.loads(last_line)
        except json.JSONDecodeError:
            return _error(f"Could not parse R output: {stdout[:500]}")

    if result.get("status") == "error":
        return result

    # normalize missing fields defensively
    result.setdefault("markers", [])
    result.setdefault("genes", [])
    result.setdefault("markers_file", str(results_dir / "Markers_with_physical_position.xlsx"))
    result.setdefault("genes_file", str(results_dir / "Genes_in_gene_regions.xlsx"))

    return result
