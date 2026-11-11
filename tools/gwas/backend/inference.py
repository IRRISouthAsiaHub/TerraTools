"""
Multiple-testing, significance categories, diagnostics and variance explained (PVE) for GWAS results.

Significance
  Three criteria are reported separately and never merged:
    Bonferroni            p <= alpha / m                      (m = markers actually tested)
    BH                    Benjamini-Hochberg adjusted p <= alpha
    Permutation           p <= the model's own simulated genome-wide threshold (only where one was computed)
  The marker category follows ONE chosen rule (`rule`); the other two are shown as flags. "Suggestive" means p <= 1/m (one false
  positive expected per scan) and not genome-wide. A top-ranked marker is never called significant: markers that pass nothing can be
  listed as "Exploratory" and carry that label. When nothing passes, the message is "No genome-wide significant associations detected".
  A multi-locus model (FarmCPU, BLINK, mrMLM, FASTmrMLM) keeps its own detection rule: its QTNs are reported with the model's LOD
  and are not given a p-value derived from that LOD.

Genomic inflation (lambda)
  Median of the 1-df chi-square statistics implied by the p-values divided by the median of the chi-square(1) distribution. It is
  reported as a number. It is not turned into a verdict and no p-value is rescaled by it.

Variance explained (PVE), as a fraction of the phenotypic variance, always on one scale:
  ols_marginal      (RSS0 - RSS1) / SS_total, RSS0 from the covariate-only model, RSS1 with the marker added, SS_total the total sum of
                    squares of the phenotype. Equal to the squared partial correlation times the share of variance the covariates
                    leave. Defined for a plain linear model only.
  mixed_marginal    b^2 v / (b^2 v + s2_g (tr K / n) + s2_e), with b the marker effect from the mixed model, v the sample variance of
                    the marker, and s2_g, s2_e the null-model variance components (REML). It is a marginal share relative to the
                    variance model of the null fit; it ignores the covariate variance and is not a joint estimate.
  conditional_joint (RSS without marker - RSS with all selected markers) / SS_total from the joint linear model of the QTNs and
                    covariates: what the marker adds given the other QTNs. Used for multi-locus models.
  Where an estimate cannot be defined (invalid statistic, singular design, no model) the value is missing and reported as
  "unavailable", not as 0 or an approximation. No equivalence with another program's PVE is claimed.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import stats

NO_HITS = "No genome-wide significant associations detected"
RULES = ("bonferroni", "bh", "permutation")
CHI2_MEDIAN = float(stats.chi2.ppf(0.5, 1))
LOD_PER_CHI2 = 1.0 / (2.0 * np.log(10.0))


@dataclass
class ModelResult:
    name: str
    p: np.ndarray
    beta: np.ndarray
    se: np.ndarray
    kind: str = "single"                  # "single" or "multilocus"
    label: str = ""                       # plain description shown next to the name
    stat_name: str = ""
    qtn_idx: np.ndarray | None = None
    qtn_lod: np.ndarray | None = None
    qtn_beta: np.ndarray | None = None
    qtn_criterion: str = ""
    info: dict = field(default_factory=dict)
    pve: np.ndarray | None = None
    pve_definition: str = ""
    permutation: dict | None = None       # output of mixed.null_threshold
    notes: list = field(default_factory=list)


# ============================================================================================ basic statistics

def tested(p):
    p = np.asarray(p, dtype=float)
    return np.isfinite(p) & (p > 0) & (p <= 1)


def genomic_lambda(p):
    """Median chi-square(1) of the observed p-values over the median of chi-square(1)."""
    p = np.asarray(p, dtype=float)
    p = p[tested(p)]
    if len(p) < 10:
        return float("nan")
    return float(np.median(stats.chi2.isf(p, 1)) / CHI2_MEDIAN)


def adjust(p):
    """(Bonferroni, BH) adjusted p-values over the tested markers; untested markers stay NaN."""
    p = np.asarray(p, dtype=float)
    ok = tested(p)
    m = int(ok.sum())
    bonf = np.full(len(p), np.nan)
    bh = np.full(len(p), np.nan)
    if m == 0:
        return bonf, bh
    pv = p[ok]
    bonf[ok] = np.minimum(pv * m, 1.0)
    order = np.argsort(pv, kind="stable")
    ranked = pv[order] * m / (np.arange(m) + 1.0)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(ranked, 1.0)
    bh[ok] = out
    return bonf, bh


def chi2_lod(p):
    """LOD equivalent of a 1-df test: chi-square / (2 ln 10)."""
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    ok = tested(p)
    out[ok] = stats.chi2.isf(p[ok], 1) * LOD_PER_CHI2
    return out


def neglog10(p):
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p) & (p > 0)
    out[ok] = -np.log10(p[ok])
    return out


# ============================================================================================ significance

def evaluate(res: ModelResult, alpha=0.05, rule="bonferroni", suggestive_p=None):
    """Attach per-marker significance to `res` (in res.info['sig']) and return it.

    Returns dict: bonferroni_p, bh_p (arrays), flags (bool arrays), category (array of str), thresholds (dict), rule_used."""
    p = np.asarray(res.p, dtype=float)
    ok = tested(p)
    m = int(ok.sum())
    bonf_adj, bh_adj = adjust(p)
    thr = {"alpha": alpha, "m_tested": m, "bonferroni_p": (alpha / m) if m else float("nan"),
           "suggestive_p": float(suggestive_p) if suggestive_p else (1.0 / m if m else float("nan"))}
    flag_bonf = ok & (p <= thr["bonferroni_p"])
    flag_bh = ok & (bh_adj <= alpha)
    perm = None
    flag_perm = np.zeros(len(p), dtype=bool)
    if res.permutation:
        key = f"{alpha:g}"
        t = res.permutation["thresholds"].get(key) or next(iter(res.permutation["thresholds"].values()))
        perm = t["p"]
        flag_perm = ok & (p <= perm)
        thr["permutation_p"] = perm
        thr["permutation_interval_p"] = [t["p_stricter"], t["p_looser"]]
    rule_used = rule
    notes = []
    if res.kind == "multilocus":
        genome = np.zeros(len(p), dtype=bool)
        if res.qtn_idx is not None and len(res.qtn_idx):
            genome[np.asarray(res.qtn_idx, dtype=int)] = True
        rule_used = "model QTN criterion"
    else:
        if rule == "permutation" and perm is None:
            notes.append("No permutation threshold exists for this model; the Bonferroni criterion was used for its categories.")
            rule_used = "bonferroni"
        genome = {"bonferroni": flag_bonf, "bh": flag_bh, "permutation": flag_perm}[rule_used]
    sugg = ok & (p <= thr["suggestive_p"]) & ~genome
    cat = np.full(len(p), "Not significant", dtype=object)
    cat[sugg] = "Suggestive"
    cat[genome] = "QTN" if res.kind == "multilocus" else "Genome-wide significant"
    cat[~ok] = "Not tested"
    sig = {"bonferroni_p": bonf_adj, "bh_p": bh_adj, "flag_bonferroni": flag_bonf, "flag_bh": flag_bh, "flag_permutation": flag_perm,
           "genome_wide": genome, "suggestive": sugg, "category": cat, "thresholds": thr, "rule_used": rule_used, "notes": notes}
    res.info["sig"] = sig
    return sig


def exploratory_idx(res: ModelResult, top_n):
    """Top-ranked markers among those that are not genome-wide significant, to be shown separately and labelled exploratory."""
    sig = res.info["sig"]
    p = np.asarray(res.p, dtype=float)
    cand = np.flatnonzero(tested(p) & ~sig["genome_wide"])
    return cand[np.argsort(p[cand], kind="stable")][:max(int(top_n), 0)]


# ============================================================================================ diagnostics

def diagnostics(res: ModelResult, bins=20):
    """Numbers only: counts, lambda, smallest p, the p-value histogram, convergence and invalid statistics."""
    p = np.asarray(res.p, dtype=float)
    ok = tested(p)
    sig = res.info.get("sig") or evaluate(res)
    hist, _ = np.histogram(p[ok], bins=bins, range=(0, 1))
    d = {"model": res.name, "snps_tested": int(ok.sum()), "invalid_statistics": int((~ok).sum()),
         "lambda": genomic_lambda(p), "min_p": float(p[ok].min()) if ok.any() else float("nan"),
         "n_genome_wide": int(sig["genome_wide"].sum()), "n_suggestive": int(sig["suggestive"].sum()),
         "n_bonferroni": int(sig["flag_bonferroni"].sum()), "n_bh": int(sig["flag_bh"].sum()),
         "n_permutation": (int(sig["flag_permutation"].sum()) if res.permutation else None),
         "p_histogram": [int(x) for x in hist], "rule": sig["rule_used"]}
    if res.kind == "multilocus":
        d["lambda_note"] = "first-step random-effect scan; its p-values are not on the scale lambda assumes"
    null = res.info.get("null")
    if null:
        d["convergence"] = ("converged" if null.get("converged") else
                            f"variance ratio at the {null.get('at_bound')} edge of the search range" if null.get("at_bound")
                            else "optimiser did not report convergence")
        d["reml_check"] = "passed" if null.get("check_ok") else "failed"
        d["h2"] = null.get("h2")
    if "n_delta_at_grid_edge" in res.info:
        d["markers_with_variance_ratio_at_edge"] = res.info["n_delta_at_grid_edge"]
    return d


# ============================================================================================ PVE

def _qr(X):
    Q, _ = np.linalg.qr(X)
    return Q


def pve_ols_marginal(y, X0, G, chunk=2000):
    """(RSS0 - RSS1)/SS_total for every marker, in closed form (see the module description)."""
    y = np.asarray(y, float)
    sst = float(((y - y.mean()) ** 2).sum())
    out = np.full(G.shape[1], np.nan)
    if sst <= 0:
        return out
    Q = _qr(X0)
    yr = y - Q @ (Q.T @ y)
    for s in range(0, G.shape[1], chunk):
        Gb = G[:, s:s + chunk]
        Gr = Gb - Q @ (Q.T @ Gb)
        gg = (Gr ** 2).sum(0)
        gy = Gr.T @ yr
        ok = gg > 1e-10 * (Gb ** 2).sum(0)
        with np.errstate(all="ignore"):
            out[s:s + chunk] = np.where(ok, gy ** 2 / gg / sst, np.nan)
    return out


def pve_mixed_marginal(beta, G, null, mean_diag_k):
    """b^2 v / (b^2 v + s2_g * mean(diag K) + s2_e), v = marker sample variance, variance components from the null REML fit."""
    v = np.nanvar(G, axis=0, ddof=1)
    bv = np.asarray(beta, float) ** 2 * v
    base = null["sigma2_g"] * mean_diag_k + null["sigma2_e"]
    with np.errstate(all="ignore"):
        out = bv / (bv + base)
    out[~np.isfinite(out)] = np.nan
    return out


def pve_conditional(y, X0, G, idx):
    """Per selected marker: (RSS of the joint model without it - RSS of the joint model) / SS_total; and the joint share of all."""
    y = np.asarray(y, float)
    idx = np.asarray(idx, dtype=int)
    sst = float(((y - y.mean()) ** 2).sum())
    out = np.full(len(idx), np.nan)
    if len(idx) == 0 or sst <= 0:
        return out, float("nan")

    def rss(cols):
        X = np.column_stack([X0] + ([G[:, cols]] if len(cols) else []))
        b = np.linalg.lstsq(X, y, rcond=None)[0]
        r = y - X @ b
        return float(r @ r)

    full = rss(list(idx))
    for k in range(len(idx)):
        rest = [j for i, j in enumerate(idx) if i != k]
        out[k] = max(rss(rest) - full, 0.0) / sst
    joint = max(rss([]) - full, 0.0) / sst
    return out, joint


def to_percent(x):
    """Fractions to percent; NaN stays NaN (reported as unavailable)."""
    return None if x is None else np.asarray(x, float) * 100.0
