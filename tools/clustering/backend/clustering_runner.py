"""
clustering_runner.py — thin wrapper around
clustering_engine.cluster_split_data(), mirroring the shape of the other
pure-Python engines (GWAS, Diversity): a .run() method returning a result
dict, so clustering_window.py's threading follows the same pattern.
"""

import tempfile
import traceback
from pathlib import Path

from tools.clustering.backend.clustering_engine import cluster_split_data
from tools.clustering.config import ClusteringConfig


class ClusteringEngine:
    def __init__(self, config: ClusteringConfig):
        self.config = config

    def run(self) -> dict:
        try:
            outdir = Path(tempfile.mkdtemp(prefix="terratools_clustering_"))
            return cluster_split_data(
                geno_file=self.config.input_file,
                n_clusters=self.config.n_clusters,
                test_size=self.config.test_size,
                outdir=outdir,
                pheno_file=self.config.pheno_file,
                keep_test_pheno=self.config.keep_test_pheno,
                n_pcs=self.config.n_pcs,
            )
        except Exception as e:
            return {
                "status": "error",
                "message": (
                    f"{e}\n\n"
                    "Expected input: a genotype matrix (CSV or TSV) with sample IDs "
                    "in the first column and SNP/marker columns after that. "
                    "Missing genotypes should be coded as -1. The optional phenotype "
                    "file has sample IDs in the first column and the trait in the second.\n\n"
                    f"Debug info:\n{traceback.format_exc()}"
                ),
            }
