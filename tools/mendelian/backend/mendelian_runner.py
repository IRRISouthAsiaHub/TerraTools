"""
mendelian_runner.py — Python side of the Mendelian tool's R backend.

Calls mendelian_cli.R via subprocess and reads back summary.json.
Design notes (learned from debugging the GWAS tool's hang earlier):
  - capture_output=True on subprocess.run() so R's stdout/stderr can never
    fill a pipe buffer and deadlock.
  - An explicit timeout, so a genuinely stuck R process fails loudly
    instead of hanging the app forever.
  - Result handoff via a JSON file written by R (summary.json), not by
    parsing subprocess stdout — R may print warnings/messages that would
    corrupt naive stdout-JSON parsing.
"""

import json
import subprocess
import tempfile
from pathlib import Path

from tools.mendelian.config import MendelianConfig
from tools.common.r_locator import find_rscript, rscript_not_found_message
from tools.common.subprocess_utils import no_window_kwargs

R_TIMEOUT_SECONDS = 600  # 10 min ceiling; adjust if your datasets are larger


class MendelianEngine:
    def __init__(self, config: MendelianConfig):
        self.config = config
        self._backend_dir = Path(__file__).resolve().parent
        self._cli_path = self._backend_dir / "mendelian_cli.R"

    def run(self) -> dict:
        rscript = find_rscript()
        if rscript is None:
            return {
                "status": "error",
                "message": rscript_not_found_message(),
            }

        outdir = Path(tempfile.mkdtemp(prefix="terratools_mendelian_"))

        cmd = [
            rscript,
            str(self._cli_path),
            "--input", self.config.input_file,
            "--pop", self.config.pop,
            "--pthr", str(self.config.pthr),
            "--missing", self.config.missing_codes,
            "--drop_unplaced", "true" if self.config.drop_unplaced else "false",
            "--outdir", str(outdir),
        ]

        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=R_TIMEOUT_SECONDS,
                stdin=subprocess.DEVNULL,
                **no_window_kwargs(),
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
            # R crashed before it could even write summary.json (e.g. a
            # missing package). Surface stderr so the real cause is visible.
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
            return summary  # already {"status": "error", "message": ...}

        # Load the per-marker table for the frontend.
        results_json_path = Path(summary["results_json"])
        with open(results_json_path, "r") as f:
            table_rows = json.load(f)

        summary["table"] = table_rows
        return summary
