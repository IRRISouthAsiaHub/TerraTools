"""
mrMLM and FASTmrMLM through the mrMLM R package.

TerraTools writes the QC-passed genotypes, the phenotype and (optionally) principal components in the package's own file formats,
runs the package in the bundled R, and reads back its two result tables:
  "<trait>_intermediate result.csv"  first-step effect and -log10 p of every marker, for each method
  "<trait>_Final result.csv"         the QTNs each method selects, with effect, LOD, -log10 p, r2 (%), MAF
No statistic is computed here: the package decides which markers are QTNs and what their LOD is. Genotypes go in as -1/0/1
(dosage - 1), with the counted allele as "genotype for code 1".
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from tools.gwas.backend import r_packages
from tools.gwas.backend.r_packages import RPackageError, TIMEOUT, run_r

METHODS = ("mrMLM", "FASTmrMLM")
RUNNER = Path(__file__).with_name("mrmlm_run.R")


class MrmlmError(RPackageError):
    pass


def _run_r(args, timeout=TIMEOUT):
    return run_r(args, timeout)


def package_version():
    """Installed mrMLM version, or None when R or the package is missing."""
    return r_packages.package_version("mrMLM")


def _write_inputs(folder, G, ids, snp, chrom, pos, y, alleles, Q):
    ids = [str(i) for i in ids]
    gen = pd.DataFrame(np.asarray(G, float).T - 1.0, columns=ids)
    gen.insert(0, "genotype for code 1", alleles)
    gen.insert(0, "pos", np.asarray(pos, dtype=np.int64))
    gen.insert(0, "chrom", np.asarray(chrom, dtype=int))
    gen.insert(0, "rs#", [str(s) for s in snp])
    gen_path = os.path.join(folder, "genotype.csv")
    gen.to_csv(gen_path, index=False)
    phe_path = os.path.join(folder, "phenotype.csv")
    with open(phe_path, "w", newline="") as fh:
        fh.write("<Phenotype>,trait1\n")
        for i, v in zip(ids, y):
            fh.write(f"{i},{float(v):.10g}\n")
    ps_path = "none"
    if Q is not None and np.asarray(Q).shape[1] > 0:
        Q = np.asarray(Q, float)
        ps_path = os.path.join(folder, "popstr.csv")
        with open(ps_path, "w", newline="") as fh:
            fh.write("<PopStr>" + "," * Q.shape[1] + "\n")
            fh.write("<ID>," + ",".join(f"PC{k + 1}" for k in range(Q.shape[1])) + "\n")
            for i, row in zip(ids, Q):
                fh.write(f"{i}," + ",".join(f"{v:.8g}" for v in row) + "\n")
    return gen_path, phe_path, ps_path


def run(G, ids, snp, chrom, pos, y, alleles=None, Q=None, methods=METHODS, radius=20.0, critical_lod=3.0, log=print):
    """Run the package. Returns {method: {'neglog10p', 'beta', 'qtn': DataFrame}} plus 'version'.

    G: samples x markers dosage 0/1/2 without missing values (mean-imputed beforehand). Q: samples x k structure covariates or None."""
    methods = [m for m in methods if m in METHODS]
    snp = np.asarray(snp).astype(str)
    if len(set(snp)) != len(snp):
        raise MrmlmError("Marker names must be unique for the mrMLM package.")
    if alleles is None:
        alleles = ["A"] * len(snp)
    else:
        counted = [str(a).split("/")[-1] for a in alleles]
        alleles = [a if a and a != "NA" else "A" for a in counted]
    with tempfile.TemporaryDirectory(prefix="terratools_mrmlm_") as tmp:
        gen, phe, ps = _write_inputs(tmp, G, ids, snp, chrom, pos, y, alleles, Q)
        out = os.path.join(tmp, "out")
        log(f"[GWAS] mrMLM package: {', '.join(methods)} on {len(snp):,} SNPs and {len(ids):,} samples")
        r = _run_r([RUNNER, gen, phe, ps, out, ",".join(methods), f"{float(radius):g}", f"{float(critical_lod):g}"])
        err = os.path.join(out, "error.txt")
        if r.returncode != 0 or os.path.isfile(err):
            msg = Path(err).read_text(errors="replace") if os.path.isfile(err) else (r.stderr or r.stdout)[-1500:]
            raise MrmlmError(f"The mrMLM R package stopped: {msg.strip()}")
        inter = _read(out, "intermediate result")
        final = _read(out, "Final result", required=False)
        version = Path(out, "version.txt").read_text().strip() if os.path.isfile(os.path.join(out, "version.txt")) else ""
    if inter is None:
        raise MrmlmError("The mrMLM package finished without writing its intermediate result file.")
    index = {s: i for i, s in enumerate(snp)}
    rows = np.array([index.get(str(s), -1) for s in inter["RS#"]])
    res = {"version": version}
    for m in methods:
        pcol = next((c for c in inter.columns if "log10(P)" in c and f"({m})" in c), None)
        ecol = next((c for c in inter.columns if c.startswith("SNP effect") and f"({m})" in c), None)
        if pcol is None:
            raise MrmlmError(f"The package output has no first-step column for {m}.")
        nlp = np.full(len(snp), np.nan)
        beta = np.full(len(snp), np.nan)
        ok = rows >= 0
        nlp[rows[ok]] = pd.to_numeric(inter[pcol], errors="coerce").to_numpy()[ok]
        if ecol:
            beta[rows[ok]] = pd.to_numeric(inter[ecol], errors="coerce").to_numpy()[ok]
        q = final[final["Method"] == m].copy() if final is not None and len(final) else pd.DataFrame()
        if len(q):
            q["index"] = [index.get(str(s), -1) for s in q["RS#"]]
            q = q[q["index"] >= 0]
        res[m] = {"neglog10p": nlp, "beta": beta, "qtn": q}
    return res


def _read(folder, kind, required=True):
    hits = [f for f in os.listdir(folder) if kind in f and f.endswith(".csv")] if os.path.isdir(folder) else []
    if not hits:
        return None
    return pd.read_csv(os.path.join(folder, sorted(hits)[0]))
