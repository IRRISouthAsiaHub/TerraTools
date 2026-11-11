"""
diversity_runner.py — thin wrapper around diversity_engine.run_diversity(),
mirroring GWASEngine's shape (a .run() method returning a result dict) so
diversity_window.py's threading follows the exact same pattern as the
already-working GWAS tool.
"""

import tempfile
import traceback
from pathlib import Path

from tools.diversity.backend.diversity_engine import run_diversity
from tools.diversity.config import DiversityConfig


class DiversityEngine:
    def __init__(self, config: DiversityConfig):
        self.config = config

    def run(self) -> dict:
        try:
            outdir = Path(tempfile.mkdtemp(prefix="terratools_diversity_"))
            result = run_diversity(
                input_file=self.config.input_file,
                n_clusters=self.config.n_clusters,
                top_blocks=self.config.top_blocks,
                radar_normalized=self.config.radar_normalized,
                outdir=outdir,
            )
            return result
        except Exception as e:
            return {
                "status": "error",
                "message": (
                    f"{e}\n\n"
                    "Expected CSV format:\n"
                    "- Column 1: Variety/Cultivar name\n"
                    "- Column 2: Sample ID\n"
                    "- Columns 3+: Haplotype blocks (numeric)\n\n"
                    f"Debug info:\n{traceback.format_exc()}"
                ),
            }
