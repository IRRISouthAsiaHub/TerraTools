"""
BLINK through the GAPIT R package (the R implementation of BLINK).

TerraTools writes the QC-passed genotypes (0/1/2 dosage, mean-imputed), the map and the phenotype, runs GAPIT with model "Blink" in
the bundled R with the package's own settings, and reads back each marker's effect and p-value. No statistic is computed here.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from tools.gwas.backend import r_packages
from tools.gwas.backend.r_packages import RPackageError, run_r

RUNNER = Path(__file__).with_name("blink_run.R")


class BlinkError(RPackageError):
    pass


def package_version():
    return r_packages.package_version("GAPIT")


def run(G, ids, snp, chrom, pos, y, log=print):
    """Returns {'p', 'beta', 'se', 'version'} in marker order (se is missing: the package does not report it)."""
    snp = np.asarray(snp).astype(str)
    with tempfile.TemporaryDirectory(prefix="terratools_blink_") as tmp:
        n, m = r_packages.write_scan_inputs(tmp, G, ids, snp, chrom, pos, y)
        log(f"[GWAS] GAPIT package: BLINK on {m:,} SNPs and {n:,} samples")
        r = run_r([RUNNER, tmp, n, m])
        err = os.path.join(tmp, "error.txt")
        if r.returncode != 0 or os.path.isfile(err):
            msg = Path(err).read_text(errors="replace") if os.path.isfile(err) else (r.stderr or r.stdout)[-1500:]
            raise BlinkError(f"The GAPIT R package stopped: {msg.strip()}")
        out = pd.read_csv(os.path.join(tmp, "blink.csv"))
        version = Path(tmp, "version.txt").read_text().strip()
    got = out.drop_duplicates("SNP").set_index("SNP")
    p = got["p"].reindex(snp).to_numpy(float)
    beta = got["effect"].reindex(snp).to_numpy(float)
    return {"p": p, "beta": beta, "se": np.full(len(snp), np.nan), "version": version}
