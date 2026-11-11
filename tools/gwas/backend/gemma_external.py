"""
Optional comparison with an external GEMMA executable.

Nothing here runs unless a GEMMA binary is found (configuration `gemma_binary`, the GEMMA environment variable, or `gemma` on the
PATH). The same QC-passed genotypes, phenotype, covariates and relationship matrix are written in GEMMA's formats, GEMMA is run with
`-lmm 1` (Wald) or `-lmm 2` (likelihood ratio), and its p-values are compared with the in-program GEMMA-style scan on the same
markers. Agreement is reported as numbers (rank correlation, largest difference in -log10 p, overlap of the top markers). The
in-program model is described as "GEMMA-style"; it is only said to agree with GEMMA when this comparison was run and the numbers
below the stated tolerances.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile

import numpy as np
from scipy import stats

INSTALL = ("GEMMA was not found, so no comparison was run. Install it from https://github.com/genetics-statistics/GEMMA/releases "
           "(download the build for your system, unzip it, and either put it on the PATH as `gemma` or set the 'GEMMA executable' "
           "option to its full path). Then run the analysis again with the GEMMA comparison switched on.")

TOLERANCE = {"spearman_min": 0.99, "max_abs_delta_neglog10p": 0.5, "top_overlap_min": 0.8}


def find_gemma(path=""):
    for cand in (path, os.environ.get("GEMMA", "")):
        if cand and os.path.isfile(cand):
            return cand
    return shutil.which("gemma") or shutil.which("gemma.exe")


def write_inputs(folder, G, y, snp_ids, X0=None, K=None):
    """BIMBAM mean-genotype file, phenotype, covariates (with intercept column) and optional relationship matrix."""
    n, m = G.shape
    geno = os.path.join(folder, "geno.txt")
    with open(geno, "w") as fh:
        for j in range(m):
            fh.write(f"{snp_ids[j]},A,T," + ",".join(f"{v:.6g}" for v in G[:, j]) + "\n")
    pheno = os.path.join(folder, "pheno.txt")
    np.savetxt(pheno, y, fmt="%.10g")
    cov = os.path.join(folder, "covariates.txt")
    np.savetxt(cov, np.ones((n, 1)) if X0 is None else np.asarray(X0, float), fmt="%.10g")
    kin = None
    if K is not None:
        kin = os.path.join(folder, "kinship.txt")
        np.savetxt(kin, K, fmt="%.10g")
    return geno, pheno, cov, kin


def compare(p_ours, p_gemma, tol=TOLERANCE):
    """Agreement numbers between two p-value vectors over the same markers."""
    p1, p2 = np.asarray(p_ours, float), np.asarray(p_gemma, float)
    ok = np.isfinite(p1) & np.isfinite(p2) & (p1 > 0) & (p2 > 0)
    if ok.sum() < 3:
        return {"status": "not comparable", "n": int(ok.sum())}
    a, b = -np.log10(p1[ok]), -np.log10(p2[ok])
    k = max(int(round(0.01 * ok.sum())), 5)
    top = len(set(np.argsort(-a)[:k]) & set(np.argsort(-b)[:k])) / k
    out = {"n": int(ok.sum()), "spearman": float(stats.spearmanr(a, b)[0]), "pearson_neglog10p": float(np.corrcoef(a, b)[0, 1]),
           "max_abs_delta_neglog10p": float(np.max(np.abs(a - b))), "median_abs_delta_neglog10p": float(np.median(np.abs(a - b))),
           "top_overlap": float(top), "top_k": k}
    out["agrees"] = bool(out["spearman"] >= tol["spearman_min"] and out["max_abs_delta_neglog10p"] <= tol["max_abs_delta_neglog10p"]
                         and top >= tol["top_overlap_min"])
    out["tolerances"] = dict(tol)
    return out


def run_gemma(G, y, snp_ids, p_ours, X0=None, K=None, test="wald", binary="", workdir=None, timeout=3600):
    """Run GEMMA and compare. Returns {'status': 'not run', 'message': install instructions} when GEMMA is unavailable."""
    exe = find_gemma(binary)
    if not exe:
        return {"status": "not run", "message": INSTALL}
    own = workdir is None
    folder = workdir or tempfile.mkdtemp(prefix="terratools_gemma_")
    try:
        geno, pheno, cov, kin = write_inputs(folder, G, y, snp_ids, X0, K)
        env = dict(cwd=folder, capture_output=True, text=True, timeout=timeout)
        if kin is None:
            r = subprocess.run([exe, "-g", geno, "-p", pheno, "-gk", "1", "-o", "k"], **env)
            if r.returncode:
                return {"status": "failed", "message": (r.stderr or r.stdout)[-800:]}
            kin = os.path.join(folder, "output", "k.cXX.txt")
        r = subprocess.run([exe, "-g", geno, "-p", pheno, "-c", cov, "-k", kin, "-lmm", "1" if test == "wald" else "2", "-o", "assoc"], **env)
        if r.returncode:
            return {"status": "failed", "message": (r.stderr or r.stdout)[-800:]}
        import pandas as pd
        tab = pd.read_csv(os.path.join(folder, "output", "assoc.assoc.txt"), sep=r"\s+")
        col = "p_wald" if test == "wald" else "p_lrt"
        got = dict(zip(tab["rs"].astype(str), tab[col]))
        p_ref = np.array([got.get(str(s), np.nan) for s in snp_ids], float)
        rep = compare(p_ours, p_ref)
        rep.update(status="run", test=test, executable=exe)
        return rep
    except subprocess.TimeoutExpired:
        return {"status": "failed", "message": "GEMMA did not finish within the time limit."}
    finally:
        if own:
            shutil.rmtree(folder, ignore_errors=True)
