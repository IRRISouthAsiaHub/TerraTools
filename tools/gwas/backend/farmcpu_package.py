"""
FarmCPU through the rMVP R package.

TerraTools writes the QC-passed genotypes (0/1/2 dosage, mean-imputed), the map and the phenotype, runs rMVP's FarmCPU in the
bundled R with the package's own settings (static bin search), and reads back its effect, standard error and p-value for every
marker. No statistic is computed here.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import pandas as pd

from tools.gwas.backend import r_packages
from tools.gwas.backend.r_packages import RPackageError, run_r

RUNNER = Path(__file__).with_name("farmcpu_run.R")


class FarmcpuError(RPackageError):
    pass


def package_version():
    return r_packages.package_version("rMVP")


def run(G, ids, snp, chrom, pos, y, max_loops=10, log=print):
    """Returns {'p', 'beta', 'se', 'version'} (arrays in marker order). G: samples x markers dosage without missing values."""
    with tempfile.TemporaryDirectory(prefix="terratools_farmcpu_") as tmp:
        n, m = r_packages.write_scan_inputs(tmp, G, ids, snp, chrom, pos, y)
        log(f"[GWAS] rMVP package: FarmCPU on {m:,} SNPs and {n:,} samples")
        r = run_r([RUNNER, tmp, n, m, int(max_loops)])
        err = os.path.join(tmp, "error.txt")
        if r.returncode != 0 or os.path.isfile(err):
            msg = Path(err).read_text(errors="replace") if os.path.isfile(err) else (r.stderr or r.stdout)[-1500:]
            raise FarmcpuError(f"The rMVP R package stopped: {msg.strip()}")
        out = pd.read_csv(os.path.join(tmp, "farmcpu.csv"))
        version = Path(tmp, "version.txt").read_text().strip()
    return {"p": out["p"].to_numpy(float), "beta": out["effect"].to_numpy(float), "se": out["se"].to_numpy(float), "version": version}
