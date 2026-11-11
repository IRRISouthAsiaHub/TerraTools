"""
One GWAS analysis from genotype data to result tables (no plots, no files).

Order of the steps, and why:
  1. Genotypes are converted to a single internal scale (dosage 0..2, NaN = missing) according to the encoding chosen by the user.
  2. Quality control (duplicate positions, SNP/sample call rate, monomorphic, MAF) on the raw matrix, with a ledger per stage.
     The analysis stops with an explanation if fewer than two SNPs remain.
  3. Relationship matrix and principal components are computed from the OBSERVED calls (a missing call contributes 0 to the
     centred genotype; allele frequencies use observed calls only). Imputation therefore does not change K or the PCs: it is
     applied after them, and only to the copy of the matrix that the marker scans need to be complete.
  4. Each requested model is run on the same marker set; each reports its own effect, SE, statistic and p-value.
  5. Significance: Bonferroni, BH and (for single-locus models) a simulation threshold from the model's own null fit are
     reported separately; one chosen rule sets the category. Multi-locus models keep their own QTN rule.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from tools.gwas.backend import blink_package, farmcpu_package, mrmlm_package, genotype_io, inference as inf, mixed, population as pop, qc as qcmod
from tools.gwas.backend.inference import ModelResult

NEEDS_MAP = ("FarmCPU", "BLINK", "mrMLM", "FASTmrMLM")
MIXED = {"MLM": "mlm", "EMMAX": "emmax", "GEMMA": "gemma"}
SIMULATED = {"OLS": "ols", "MLM": "mlm", "EMMAX": "emmax", "GEMMA": "gemma"}
MTA_COLUMNS = ["Model", "Category", "Rank", "SNP", "Chr", "Position", "Alleles", "MAF", "Effect", "SE", "P_value", "P_Bonferroni",
               "P_BH", "neg_log10P", "LOD", "PVE_%", "PVE_basis", "Criterion"]
SCAN_COLUMNS = ["SNP", "Chromosome", "Position", "Alleles", "MAF", "Effect", "SE", "P_value", "P_Bonferroni", "P_BH", "LOD",
                "PVE_percent", "Category"]


@dataclass
class Analysis:
    cfg: object
    data: genotype_io.GenotypeData
    qc: qcmod.QCResult
    ids: np.ndarray
    y: np.ndarray
    snp: np.ndarray
    chrom: np.ndarray
    pos: np.ndarray
    alleles: np.ndarray
    maf: np.ndarray
    G: np.ndarray
    G_observed: np.ndarray
    K: np.ndarray
    kinship_info: dict
    pc: dict
    structure: dict
    results: dict = field(default_factory=dict)
    errors: dict = field(default_factory=dict)
    pc_compare: list = field(default_factory=list)
    pc_compare_p: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)
    store: dict = field(default_factory=dict)


# ============================================================================================ helpers

def read_genotypes(cfg, path, generated):
    """GenotypeData from the GWAS-ready file. Files built by the program are dosage by construction."""
    enc = "dosage" if generated else cfg.genotype_encoding
    return genotype_io.read_gwas_file(path=path, encoding=enc, extra_missing=cfg.missing_code_list())


def _design(n, pc, k):
    k = int(max(k, 0))
    if k == 0:
        return np.ones((n, 1)), ["Intercept"], []
    k = min(k, pc["scores"].shape[1], max(n - 6, 0))
    names = [f"PC{i + 1}" for i in range(k)]
    return mixed.design_matrix(n, pc["scores"][:, :k], names)


def _describe_covariates(used, dropped):
    base = "no covariates" if len(used) == 1 else "covariates: " + ", ".join(used[1:])
    return base + (f" (dropped as redundant: {', '.join(d[0] for d in dropped)})" if dropped else "")


def pcs_from_eigen(w, U, n_snps, k):
    """Leading eigenvectors of the relationship matrix as principal components (eigenvalues in descending order)."""
    order = np.argsort(w)[::-1]
    ws, Vs = w[order], U[:, order]
    pos_w = np.clip(ws, 0, None)
    return {"scores": Vs[:, :max(int(k), 0)], "eigenvalues": ws, "variance_explained": pos_w / max(pos_w.sum(), 1e-300),
            "n_snps": int(n_snps), "ld_pruned": False, "marker_idx": None}


def _eig_info(K):
    w, U, info = pop.eig_psd(K)
    return (w, U, info)


def _pvalue_curve(p, cap=4000):
    p = np.sort(np.asarray(p, float)[inf.tested(p)])
    if len(p) > cap:                                    # thin the plain bulk, keep the extreme tail whole
        idx = np.unique(np.concatenate([np.linspace(0, len(p) - 1, cap - 200).astype(int), np.arange(min(200, len(p)))]))
        return p[idx], len(p), idx
    return p, len(p), np.arange(len(p))


# ============================================================================================ models

def _run_single(name, cfg, A, X0, used, dropped, eig, nf_global, log):
    y, G, K = A.y, A.G, A.K
    cov_text = _describe_covariates(used, dropped)
    if name == "OLS":
        r = mixed.ols_scan(y, X0, G)
        label = ("OLS - uncorrected for structure (" + cov_text + ")") if len(used) == 1 else \
            f"OLS - corrected with {len(used) - 1} principal component(s) ({cov_text})"
        res = ModelResult("OLS", r["p"], r["beta"], r["se"], "single", label, r["info"]["stat_name"], info=r["info"])
        res.pve = inf.pve_ols_marginal(y, X0, G)
        res.pve_definition = "OLS marginal: (RSS0 - RSS1) / SS_total, covariates in both models"
        return res
    method = MIXED[name]
    test = cfg.gemma_test if method == "gemma" else "wald"
    if cfg.loco:
        r = mixed.mixed_scan_loco(y, X0, G, A.chrom, pop.loco_kinship(A.G_observed, A.chrom, A.kinship_info.get("_cols")), method, test)
        kin = "leave-one-chromosome-out kinship"
    else:
        r = mixed.mixed_scan(y, X0, G, K, method, test, eig=eig)
        kin = f"kinship ({A.kinship_info['kinship_method']}, {A.kinship_info['n_snps']:,} SNPs)"
    how = {"mlm": "variance ratio from the null model, fixed; sigma^2 re-estimated per marker; t test",
           "emmax": "both variance components from the null model, fixed for all markers; normal (chi-square) test",
           "gemma": f"variance ratio re-estimated per marker; {'likelihood-ratio' if test == 'lrt' else 'Wald'} test"}[method]
    res = ModelResult(name, r["p"], r["beta"], r["se"], "single", f"{name} - {kin}; {how}; {cov_text}", r["info"]["stat_name"],
                      info=r["info"])
    res.pve = inf.pve_mixed_marginal(r["beta"], G, nf_global, float(np.mean(np.diag(K))))
    res.pve_definition = "mixed-model marginal: b^2 v / (b^2 v + s2_g mean(diag K) + s2_e), null-model variance components"
    if test == "lrt":
        res.notes.append("The effect and its SE come from the maximum-likelihood fit of the alternative model; the p-value is the "
                         "likelihood-ratio test.")
    return res


def _run_iterative(name, cfg, A, log):
    """FarmCPU (rMVP package) and BLINK (GAPIT package): the package computes the scan; nothing is re-implemented here."""
    y, G = A.y, A.G
    n = len(y)
    if name == "FarmCPU":
        out = farmcpu_package.run(G, A.ids, A.snp, A.chrom, A.pos, y, max_loops=int(cfg.farmcpu_iterations), log=log)
        label = f"FarmCPU - run by the rMVP R package {out['version']} (iterative model with pseudo-QTNs as covariates, no kinship matrix)"
        package = "rMVP"
    else:
        out = blink_package.run(G, A.ids, A.snp, A.chrom, A.pos, y, log=log)
        label = f"BLINK - run by the GAPIT R package {out['version']} (iterative model with LD pruning and BIC selection, no kinship matrix)"
        package = "GAPIT"
    p = np.asarray(out["p"], dtype=float)
    res = ModelResult(name, p, np.asarray(out["beta"], float), np.asarray(out["se"], float), "single", label, "t",
                      info={"package": package, "package_version": out["version"]})
    one = np.ones((n, 1))
    res.pve = inf.pve_ols_marginal(y, one, G)
    res.pve_definition = "OLS marginal, computed by TerraTools (the package does not return its selected markers)"
    res.notes.append("Computed by the " + package + " R package. Significance applies the chosen rule to the package's p-values; no "
                     "permutation threshold exists for this model." + (" The package gives no standard error, so none is shown." if name == "BLINK" else ""))
    return res


def _run_multilocus(name, cfg, A, log):
    """mrMLM and FASTmrMLM: both come from one run of the mrMLM R package (nothing is re-implemented here)."""
    if "mrmlm" not in A.store:
        wanted = [m for m in ("mrMLM", "FASTmrMLM") if m in cfg.models]
        k_pc = max(min(int(cfg.mrmlm_n_pca), A.pc["scores"].shape[1], len(A.y) - 6), 0)
        Q = A.pc["scores"][:, :k_pc] if k_pc > 0 else None
        out = mrmlm_package.run(A.G, A.ids, A.snp, A.chrom, A.pos, A.y, alleles=A.alleles, Q=Q, methods=wanted,
                                radius=float(cfg.mrmlm_radius_kb), critical_lod=float(cfg.mrmlm_lod), log=log)
        A.store["mrmlm"] = (out, k_pc, Q)
    out, k_pc, Q = A.store["mrmlm"]
    y, G, n = A.y, A.G, len(A.y)
    r = out[name]
    with np.errstate(over="ignore"):
        p = np.where(np.isfinite(r["neglog10p"]), np.power(10.0, -r["neglog10p"]), np.nan)
    q = r["qtn"]
    qtn = q["index"].to_numpy(dtype=int) if len(q) else np.zeros(0, dtype=int)
    lod = q["LOD score"].to_numpy(dtype=float) if len(q) else np.zeros(0)
    beta = np.asarray(r["beta"], float)
    qbeta = q["QTN effect"].to_numpy(dtype=float) if len(q) else np.zeros(0)
    beta[qtn] = qbeta
    res = ModelResult(name, p, beta, np.full(len(p), np.nan), "multilocus",
                      f"{name} - multi-locus model run by the mrMLM R package {out['version']}; structure covariates: {k_pc} PC(s)", "Wald",
                      qtn_idx=qtn, qtn_lod=lod, qtn_beta=qbeta, qtn_criterion=f"LOD >= {cfg.mrmlm_lod:g}, selected by the mrMLM package",
                      info={"package": "mrMLM", "package_version": out["version"]})
    X0 = np.ones((n, 1)) if Q is None else np.column_stack([np.ones((n, 1)), Q])
    pv = np.full(G.shape[1], np.nan)
    if len(qtn):
        pv[qtn], _joint = inf.pve_conditional(y, X0, G, qtn)
    res.pve = pv
    res.pve_definition = "conditional on the other QTNs and the structure covariates (joint least squares, computed by TerraTools); unavailable for non-QTN markers"
    res.notes.append("Computed by the mrMLM R package. The p-value shown for every marker, QTNs included, is the package's first-step "
                     "value; a QTN is reported because of its LOD, which is not converted into a p-value. The package gives no standard "
                     "error, so none is shown.")
    return res


# ============================================================================================ tables

def _fmt(v, nd=None):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(v):
        return None
    return round(v, nd) if nd is not None else v


def _model_arrays(A, res):
    sig = res.info["sig"]
    p = np.asarray(res.p, float)
    lod = inf.chi2_lod(p)
    if res.kind == "multilocus" and len(res.qtn_idx):
        lod = lod.copy()
        lod[res.qtn_idx] = res.qtn_lod
    return sig, p, lod


def scan_columns(A, res):
    """Per-marker columns for the full-scan sheet, in map order."""
    sig, p, lod = _model_arrays(A, res)
    order = np.lexsort((A.pos, A.chrom))
    pve = None if res.pve is None else np.asarray(res.pve, float) * 100.0
    cols = {"SNP": A.snp.astype(str), "Chromosome": A.chrom.astype(int), "Position": A.pos.astype(np.int64),
            "Alleles": A.alleles, "MAF": A.maf, "Effect": res.beta, "SE": res.se, "P_value": p,
            "P_Bonferroni": sig["bonferroni_p"], "P_BH": sig["bh_p"], "LOD": lod,
            "PVE_percent": pve if pve is not None else np.full(len(p), np.nan), "Category": sig["category"]}
    return {k: np.asarray(v)[order] for k, v in cols.items()}


def mta_rows(A):
    cfg = A.cfg
    rows, dropped_pve = [], 0
    min_pve = float(cfg.minimum_pve) if cfg.show_pve else 0.0
    for name, res in A.results.items():
        sig, p, lod = _model_arrays(A, res)
        pve = None if res.pve is None else np.asarray(res.pve, float) * 100.0
        thr = sig["thresholds"]

        def row(i, category, rank, criterion):
            lp = inf.neglog10(p[i:i + 1])[0]
            return {"Model": name, "Category": category, "Rank": rank, "SNP": str(A.snp[i]), "Chr": int(A.chrom[i]),
                    "Position": int(A.pos[i]), "Alleles": str(A.alleles[i]), "MAF": _fmt(A.maf[i], 4), "Effect": _fmt(res.beta[i], 6),
                    "SE": _fmt(res.se[i], 6), "P_value": _fmt(p[i]), "P_Bonferroni": _fmt(sig["bonferroni_p"][i]),
                    "P_BH": _fmt(sig["bh_p"][i]), "neg_log10P": _fmt(lp, 4), "LOD": _fmt(lod[i], 3),
                    "PVE_%": None if pve is None else _fmt(pve[i], 3), "PVE_basis": res.pve_definition if pve is not None else "",
                    "Criterion": criterion}

        gw = np.flatnonzero(sig["genome_wide"])
        gw = gw[np.argsort(p[gw], kind="stable")]
        if min_pve > 0 and pve is not None:
            keep = np.array([(pve[i] >= min_pve) if np.isfinite(pve[i]) else False for i in gw], dtype=bool)
            dropped_pve += int((~keep).sum())
            gw = gw[keep]
        for rank, i in enumerate(gw, 1):
            if res.kind == "multilocus":
                k = int(np.flatnonzero(res.qtn_idx == i)[0])
                crit = f"{res.qtn_criterion} (LOD {res.qtn_lod[k]:.2f}); p-value shown is the first-step scan"
            else:
                rule = sig["rule_used"]
                crit = {"bonferroni": f"Bonferroni: p <= {thr['bonferroni_p']:.2e}", "bh": f"BH adjusted p <= {thr['alpha']:g}",
                        "permutation": f"simulated threshold: p <= {thr.get('permutation_p', float('nan')):.2e}"}[rule]
            rows.append(row(i, res.info["sig"]["category"][i], rank, crit))
    return rows, dropped_pve


def summary_rows(A):
    out = []
    for name, res in A.results.items():
        d = res.info["diag"]
        sig = res.info["sig"]
        thr = sig["thresholds"]
        none = d["n_genome_wide"] == 0
        out.append({"Model": name, "Description": res.label, "SNPs_Tested": d["snps_tested"], "Invalid_Statistics": d["invalid_statistics"],
                    "Genomic_Lambda": _fmt(d["lambda"], 3), "Min_P": _fmt(d["min_p"]), "Bonferroni_threshold_p": _fmt(thr["bonferroni_p"]),
                    "Bonferroni_hits": d["n_bonferroni"], "BH_hits": d["n_bh"],
                    "Permutation_threshold_p": _fmt(thr.get("permutation_p")), "Permutation_hits": d["n_permutation"],
                    "Category_rule": d["rule"], "Genome_wide": d["n_genome_wide"], "Suggestive": d["n_suggestive"],
                    "Result": inf.NO_HITS if none else f"{d['n_genome_wide']} genome-wide significant"})
    return out


# ============================================================================================ the analysis

def analyse(cfg, data: genotype_io.GenotypeData, labels=None, alleles=None, log=print):
    notes = list(data.notes)
    log(f"[GWAS] {data.encoding['message']}")
    settings = qcmod.QCSettings(min_maf=float(cfg.maf), max_snp_missing=float(cfg.max_snp_missing),
                                max_sample_missing=float(cfg.max_sample_missing),
                                drop_duplicate_positions=bool(cfg.drop_duplicate_positions), imputation=cfg.imputation)
    qc = qcmod.run_qc(data.G, data.snp_ids, data.chrom, data.pos, data.ids, settings)          # QCError stops the run here
    for st in qc.ledger:
        log(f"[GWAS] QC {st['stage']}: {st['samples_left']:,} samples, {st['snps_left']:,} SNPs ({st['rule']})")
    notes += qc.notes
    si, mi = qc.sample_idx, qc.snp_idx
    n, m = len(si), len(mi)
    y = data.y[si]
    if np.std(y) == 0:
        raise qcmod.QCError("The phenotype is constant in the samples that passed quality control.")
    chrom, pos = data.chrom[mi], data.pos[mi]
    if alleles is not None and len(alleles) == len(data.snp_ids):
        allele_txt = np.asarray(alleles, dtype=object)[mi]
    elif data.alt is not None:
        allele_txt = np.array([f"{a}/{b}" for a, b in zip(data.ref[mi], data.alt[mi])], dtype=object)
    else:
        allele_txt = np.array(["NA"] * m, dtype=object)

    # relationship matrix and principal components from the observed calls
    K, kinfo, kcols = pop.kinship(qc.G_observed, chrom, pos, cfg.kinship_method, ld_r2=float(cfg.ld_prune_r2))
    kinfo["_cols"] = kcols
    w, U, einfo = eig = _eig_info(K)
    kinfo["eigen"] = einfo
    if not einfo["psd"]:
        notes.append(f"The relationship matrix is not positive semi-definite (smallest eigenvalue {einfo['min_eigenvalue']:.3g}). "
                     "Negative eigenvalues were set to 0 and nothing else was changed.")
    need = max([int(cfg.ols_pcs), int(cfg.mlm_pcs), int(cfg.emmax_pcs), int(cfg.gemma_pcs), int(cfg.mrmlm_n_pca), 3])
    pc = pcs_from_eigen(w, U, kinfo["n_snps"], min(need, n - 6))          # from the decomposition of K already done: no extra pass over the genotypes
    structure = pop.structure_summary(pc, K, labels=None if labels is None else np.asarray(labels)[si])
    log(f"[GWAS] kinship: {kinfo['method']} on {kinfo['n_snps']:,} SNPs")

    A = Analysis(cfg=cfg, data=data, qc=qc, ids=data.ids[si], y=y, snp=data.snp_ids[mi], chrom=chrom, pos=pos, alleles=allele_txt,
                 maf=qc.maf, G=qc.G, G_observed=qc.G_observed, K=K, kinship_info=kinfo, pc=pc, structure=structure, notes=notes)

    map_ok = data.marker_map.valid
    for name in cfg.models:
        try:
            if name in NEEDS_MAP and not map_ok:
                raise ValueError("this model needs physical marker positions (marker names such as Chr3_123456); "
                                 + "; ".join(data.marker_map.notes))
            if name in ("OLS", *MIXED):
                k = {"OLS": cfg.ols_pcs, "MLM": cfg.mlm_pcs, "EMMAX": cfg.emmax_pcs, "GEMMA": cfg.gemma_pcs}[name]
                X0, used, dropped = _design(n, pc, k)
                nf = mixed.null_fit(U.T @ y, U.T @ X0, w)
                res = _run_single(name, cfg, A, X0, used, dropped, eig, nf, log)
                A.results[name] = res
                res._design = (X0, used, dropped)
            elif name in ("FarmCPU", "BLINK"):
                A.results[name] = _run_iterative(name, cfg, A, log)
            elif name in ("mrMLM", "FASTmrMLM"):
                A.results[name] = _run_multilocus(name, cfg, A, log)
            else:
                raise ValueError("unknown model")
            log(f"[GWAS] {name} done")
        except Exception as ex:                          # a failing model is reported, the others still run
            A.errors[name] = f"{type(ex).__name__}: {ex}" if not isinstance(ex, (ValueError, mixed.ModelError)) else str(ex)
            log(f"[GWAS] {name} failed: {A.errors[name]}")
    if not A.results:
        raise qcmod.QCError("No model completed. " + " | ".join(f"{k}: {v}" for k, v in A.errors.items()))

    # significance by simulation under each single-locus model's own null fit
    R = int(cfg.permutations or 0)
    for name, res in A.results.items():
        if name not in SIMULATED or R <= 0:
            continue
        if cfg.loco and name in MIXED:
            res.notes.append("The simulated threshold is not available together with leave-one-chromosome-out kinship (the null model "
                             "differs per chromosome); Bonferroni and BH are reported.")
            continue
        try:
            X0 = res._design[0]
            res.permutation = mixed.null_threshold(y, X0, A.G, K, SIMULATED[name], n_perm=R, alphas=(float(cfg.alpha),),
                                                   null=cfg.permutation_null, seed=int(cfg.permutation_seed), eig=eig)
            t = res.permutation["thresholds"][f"{float(cfg.alpha):g}"]
            log(f"[GWAS] {name}: simulated threshold p <= {t['p']:.2e} ({res.permutation['n_successful']}/{R} replicates)")
        except Exception as ex:
            res.notes.append(f"Simulated threshold failed: {ex}")

    for name, res in A.results.items():
        inf.evaluate(res, alpha=float(cfg.alpha), rule=cfg.significance_rule)
        res.info["diag"] = inf.diagnostics(res)

    # optional comparison of the GEMMA-style scan with an external GEMMA executable
    if getattr(cfg, "compare_gemma", False) and "GEMMA" in A.results:
        from tools.gwas.backend import gemma_external
        g = A.results["GEMMA"]
        rep = gemma_external.run_gemma(A.G, y, list(A.snp), g.p, X0=g._design[0], K=K, test=cfg.gemma_test, binary=cfg.gemma_binary)
        g.info["external_gemma"] = rep
        if rep["status"] == "run":
            g.notes.append(f"Compared with GEMMA on the same inputs: Spearman {rep['spearman']:.4f}, largest difference in -log10 p "
                           f"{rep['max_abs_delta_neglog10p']:.2f}, top-{rep['top_k']} overlap {rep['top_overlap']:.0%}; "
                           + ("within the stated tolerances." if rep["agrees"] else "outside the stated tolerances, so no equivalence is claimed."))
        else:
            g.notes.append("GEMMA comparison: " + rep.get("message", rep["status"]))

    return A
