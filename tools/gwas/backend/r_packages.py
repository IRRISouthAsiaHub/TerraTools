"""Calling the bundled R for the GWAS models that come from R packages (mrMLM: mrMLM and FASTmrMLM; rMVP: FarmCPU)."""

from __future__ import annotations

import subprocess
import sys

from tools.common.r_locator import find_rscript, get_r_env, rscript_not_found_message

TIMEOUT = 6 * 3600


class RPackageError(ValueError):
    pass


def run_r(args, timeout=TIMEOUT):
    rscript = find_rscript()
    if rscript is None:
        raise RPackageError(rscript_not_found_message())
    flags = {"creationflags": 0x08000000} if sys.platform.startswith("win") else {}           # no console window
    return subprocess.run([rscript, *map(str, args)], capture_output=True, text=True, timeout=timeout, env=get_r_env(rscript), **flags)


def package_version(pkg):
    """Installed version of an R package in the bundled R, or None when R or the package is missing."""
    try:
        r = run_r(["-e", f'if (requireNamespace("{pkg}", quietly=TRUE)) cat(as.character(packageVersion("{pkg}")))'], timeout=180)
    except Exception:
        return None
    return (r.stdout or "").strip().splitlines()[-1] if (r.stdout or "").strip() else None


def write_scan_inputs(folder, G, ids, snp, chrom, pos, y):
    """Files read by farmcpu_run.R and blink_run.R: geno.bin (float32, sample-major), map.csv and phenotype.csv."""
    import os

    import numpy as np
    import pandas as pd

    G = np.asarray(G, dtype="<f4")
    snp = np.asarray(snp).astype(str)
    if len(set(snp)) != G.shape[1]:
        raise RPackageError("Marker names must be unique for the R packages.")
    G.tofile(os.path.join(folder, "geno.bin"))                       # sample-major = column-major of the markers x samples matrix
    pd.DataFrame({"SNP": snp, "CHROM": np.asarray(chrom, dtype=int), "POS": np.asarray(pos, dtype=np.int64)}).to_csv(
        os.path.join(folder, "map.csv"), index=False)
    pd.DataFrame({"Taxa": [str(i) for i in ids], "trait1": np.asarray(y, float)}).to_csv(os.path.join(folder, "phenotype.csv"), index=False)
    return G.shape
