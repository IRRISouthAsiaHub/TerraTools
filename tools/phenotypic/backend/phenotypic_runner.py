"""
phenotypic_runner.py — Python side of the Phenotypic Analysis tool's R
backend. Calls phenotypic_cli.R via subprocess and reads back summary.json,
then converts each output CSV into a JSON array so the frontend doesn't
need its own CSV parser.

Same deadlock-safe / timeout-safe subprocess pattern as mendelian_runner.py:
capture_output=True, an explicit timeout, and file-based result handoff
instead of parsing stdout.
"""

import csv
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

from tools.phenotypic.config import PhenotypicConfig
from tools.common.r_locator import find_rscript, rscript_not_found_message

R_TIMEOUT_SECONDS = 900  # this pipeline fits many mixed models — allow up to 15 min


def _csv_to_records(path: str) -> list:
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        return [row for row in reader]


class PhenotypicEngine:
    def __init__(self, config: PhenotypicConfig):
        self.config = config
        self._backend_dir = Path(__file__).resolve().parent
        self._cli_path = self._backend_dir / "phenotypic_cli.R"

    def run(self) -> dict:
        rscript = find_rscript()
        if rscript is None:
            return {"status": "error", "message": rscript_not_found_message()}

        outdir = Path(tempfile.mkdtemp(prefix="terratools_phenotypic_"))

        cmd = [
            rscript,
            str(self._cli_path),
            "--input", self.config.input_file,
            "--trait", self.config.trait,
            "--interaction", "true" if self.config.include_interaction else "false",
            "--mode", self.config.analysis_mode,
            "--environments", "||".join(self.config.selected_environments),
            "--n_clusters", str(self.config.n_clusters),
            "--clustering_method", self.config.clustering_method,
            "--color_scheme", self.config.color_scheme,
            "--outdir", str(outdir),
        ]

        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=R_TIMEOUT_SECONDS,
            )
        except subprocess.TimeoutExpired:
            return {
                "status": "error",
                "message": f"R process did not finish within {R_TIMEOUT_SECONDS}s (timed out).",
            }
        except Exception as e:
            return {"status": "error", "message": f"Failed to launch Rscript: {e}"}

        summary_path = outdir / "summary.json"
        if not summary_path.exists():
            return {
                "status": "error",
                "message": (
                    "R process exited without producing a result "
                    f"(exit code {proc.returncode}). stderr:\n{proc.stderr[-2000:]}"
                ),
            }

        with open(summary_path, "r") as f:
            summary = json.load(f)

        if summary.get("status") != "success":
            return summary

        # Convert every table CSV into a JSON array for the frontend.
        tables_json = {}
        for name, path in (summary.get("tables") or {}).items():
            try:
                tables_json[name] = _csv_to_records(path)
            except Exception as e:
                tables_json[name] = []
                summary.setdefault("errors", {})[f"{name}_load"] = str(e)
        summary["tables_data"] = tables_json

        return summary
