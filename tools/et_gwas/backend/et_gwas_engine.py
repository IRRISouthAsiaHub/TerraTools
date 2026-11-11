"""
et_gwas_engine.py - Extreme Trait GWAS (Et-GWAS).

Method: Gnanapragasam N, et al. (2024) "Extreme trait GWAS (Et-GWAS): Unraveling rare variants in the
3,000 rice genome", Life Science Alliance, doi:10.26508/lsa.202302352. The statistics follow the authors'
code (github.com/IRRI-South-Asia-Hub/Et-GWAS, after Yang et al. 2015 XP-GWAS):

1. Sampling   - the lines are sorted by phenotype; the lowest and the highest `bulk_pct` % form the low and
                high bulks, and an equal number of lines drawn at random from the remaining (middle) lines
                forms the random bulk.
2. Pooling    - for every SNP the reference and alternative alleles of the lines of each bulk are counted
                (two alleles per line; missing calls are skipped), which gives the pooled allele frequency.
3. Screening  - each SNP is tested with a binomial GLM (logit link) of the (ref, alt) counts of the three bulks
                on their rank (low = 1, random = 2, high = 3), i.e. a test for a monotone trend of the allele
                frequency. One pseudo-count is added to each count. The statistic is the deviance
                explained by the trend (1 d.f. chi-square). The statistics are divided by an inflation factor
                (lambda, regression of the sorted statistics on the expected chi-square quantiles; Devlin &
                Roeder 1999) before the p-value is taken; q-values are Benjamini-Hochberg FDR.
4. PVE        - PVE = 2 b^2 f (1-f) / (2 b^2 f (1-f) + 2 se(b)^2 N f (1-f)), b = slope, f = MAF, N = lines in the bulks.

Everything here is vectorised over SNPs (one IRLS run for the whole genome), so a 50,000-SNP panel takes seconds.
"""

from __future__ import annotations

import math
import os
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import expit, logit

BULKS = ("low", "random", "high")
RANK_X = np.array([1.0, 2.0, 3.0])            # the bulks as a linear sequence: low, random, high


class EtGwasError(ValueError):
    """A problem the user can fix; shown as is in the app."""


# ──────────────────────────── 1. sampling ────────────────────────────

def make_bulks(y, bulk_pct, seed=1):
    """Return (low_idx, random_idx, high_idx, thresholds). Ties at a bulk border stay in the bulk."""
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 20:
        raise EtGwasError(f"Only {n} lines have a phenotype and a genotype; Et-GWAS needs at least 20.")
    val = int(math.ceil(bulk_pct * n / 100.0))
    if val < 3:
        raise EtGwasError(f"A bulk of {bulk_pct:g}% of {n} lines has only {val} line(s); use a larger percentage.")
    if 2 * val > n * 0.8:
        raise EtGwasError(f"The two extreme bulks would take {2 * val} of the {n} lines; use a smaller percentage (<= 40).")
    order = np.argsort(y, kind="stable")
    low_thr = float(y[order[:val]].max())
    high_thr = float(y[order[::-1][:val]].min())
    if low_thr >= high_thr:
        raise EtGwasError("The phenotype has too many identical values to separate a low and a high bulk.")
    group = np.where(y <= low_thr, 1, np.where(y >= high_thr, 3, 2))
    low, high, mid = (np.where(group == g)[0] for g in (1, 3, 2))
    k = min(len(low), len(high))
    if len(mid) < k:
        raise EtGwasError("Too few lines are left outside the extreme bulks to draw a random bulk.")
    rng = np.random.default_rng(seed)
    rand = np.sort(rng.choice(mid, size=k, replace=False))
    return low, rand, high, {"low_threshold": low_thr, "high_threshold": high_thr, "group": group}


# ──────────────────────────── 2. pooling ────────────────────────────

def pooled_counts(G, idx):
    """Reference and alternative allele counts of the lines `idx` (G = lines x SNPs, dosage 0/1/2 of the
    alternative allele, NaN = missing)."""
    sub = G[idx]
    called = ~np.isnan(sub)
    alt = np.nansum(sub, axis=0)
    ref = 2.0 * called.sum(axis=0) - alt
    return ref, alt


# ──────────────────────────── 3. the GLM trend test ────────────────────────────

def _deviance(y, n, mu):
    return 2.0 * (y * np.log(y / (n * mu)) + (n - y) * np.log((n - y) / (n * (1.0 - mu)))).sum(axis=0)


def trend_test(ref, alt, x=RANK_X, pseudo=1.0, max_iter=25):
    """Binomial GLM  logit(P(ref)) ~ 1 + x  on the three bulks, for every SNP at once.

    ref, alt: (3, m) counts in the order of `x`. Returns a dict of (m,) arrays: stat (deviance explained
    by the slope, 1 d.f.), slope and se (log-odds of the REFERENCE allele per bulk step), null_dev, res_dev.
    The same numbers as R's  glm(cbind(ref, alt) ~ x, family = binomial)  and anova()."""
    y = np.asarray(ref, dtype=float) + pseudo
    n = y + np.asarray(alt, dtype=float) + pseudo
    xx = x.reshape(-1, 1)
    p0 = y.sum(axis=0) / n.sum(axis=0)
    null_dev = _deviance(y, n, np.broadcast_to(p0, y.shape))

    eta = logit((y + 0.5) / (n + 1.0))                        # glm's starting values for the binomial family
    dev_old = np.full(y.shape[1], np.inf)
    for _ in range(max_iter):
        mu = expit(eta)
        w = n * mu * (1.0 - mu)
        z = eta + (y - n * mu) / w
        s0, s1, s2 = w.sum(0), (w * xx).sum(0), (w * xx * xx).sum(0)
        t0, t1 = (w * z).sum(0), (w * xx * z).sum(0)
        det = s0 * s2 - s1 ** 2
        a = (s2 * t0 - s1 * t1) / det
        b = (s0 * t1 - s1 * t0) / det
        eta = a + b * xx
        dev = _deviance(y, n, expit(eta))
        if np.all(np.abs(dev - dev_old) < 1e-8 * (0.1 + np.abs(dev))):
            break
        dev_old = dev
    mu = expit(eta)
    w = n * mu * (1.0 - mu)
    s0, s1, s2 = w.sum(0), (w * xx).sum(0), (w * xx * xx).sum(0)
    se = np.sqrt(s0 / (s0 * s2 - s1 ** 2))
    return {"stat": np.maximum(null_dev - dev, 0.0), "slope": b, "se": se, "null_dev": null_dev, "res_dev": dev}


def estimate_lambda(stat):
    """Genomic-control inflation factor as in the authors' estlambda(method = "regression"):
    slope through the origin of the sorted statistics on the expected chi-square(1) quantiles."""
    s = np.asarray(stat, dtype=float)
    s = s[np.isfinite(s)]
    s = np.sort(s[np.abs(s) >= 1e-8])
    k = len(s)
    if k < 2:
        return 1.0
    a = 3.0 / 8.0 if k <= 10 else 0.5
    ppoi = np.sort(stats.chi2.isf((np.arange(1, k + 1) - a) / (k + 1 - 2 * a), 1))
    return float((s * ppoi).sum() / (ppoi ** 2).sum())


def bh_qvalues(p):
    """Benjamini-Hochberg adjusted p-values (R's p.adjust(method = "fdr"))."""
    p = np.asarray(p, dtype=float)
    m = len(p)
    order = np.argsort(p)
    q = p[order] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(q, 1.0)
    return out


def pve_formula(slope, se, maf, n_lines):
    f = maf * (1.0 - maf)
    num = 2.0 * slope ** 2 * f
    return 100.0 * num / (num + 2.0 * se ** 2 * n_lines * f)


# ──────────────────────────── regions ────────────────────────────

def merge_regions(table, flag, merge_bp):
    """Group the flagged SNPs of each chromosome into regions: SNPs closer than `merge_bp` belong together."""
    rows = []
    sig = table[flag].sort_values(["Chr", "Position"])
    for chrom, d in sig.groupby("Chr", sort=True):
        pos = d["Position"].to_numpy()
        start = 0
        for i in range(1, len(d) + 1):
            if i == len(d) or pos[i] - pos[i - 1] > merge_bp:
                part = d.iloc[start:i]
                peak = part.loc[part["P_value"].idxmin()]
                rows.append({"Chr": int(chrom), "Start": int(part.Position.min()), "End": int(part.Position.max()),
                             "Size_bp": int(part.Position.max() - part.Position.min() + 1), "N_SNPs": len(part),
                             "Peak_SNP": peak["SNP"], "Peak_P": float(peak["P_value"]),
                             "Peak_Position": int(peak["Position"]),
                             "Enriched_in": peak["Counted_allele_enriched_in"],
                             "Max_PVE_%": float(part["PVE_%"].max())})
                start = i
    out = pd.DataFrame(rows)
    if len(out):
        out = out.sort_values("Peak_P").reset_index(drop=True)
        out.insert(0, "Region", [f"ET{i + 1}" for i in range(len(out))])
    return out


# ──────────────────────────── the analysis ────────────────────────────

def analyse(G, y, ids, chrom, pos, snp_names, p):
    """G: lines x SNPs dosage (NaN = missing). p: parameter dict. Returns dict of tables and arrays."""
    n, m = G.shape
    low, rand, high, info = make_bulks(y, p["bulk_pct"], p["seed"])
    bulk_idx = {"low": low, "random": rand, "high": high}

    ref = np.zeros((3, m))
    alt = np.zeros((3, m))
    for k, name in enumerate(BULKS):
        ref[k], alt[k] = pooled_counts(G, bulk_idx[name])

    called = (~np.isnan(G)).sum(axis=0)
    freq_all = np.where(called > 0, np.nansum(G, axis=0) / np.maximum(2.0 * called, 1), np.nan)
    maf = np.minimum(freq_all, 1 - freq_all)

    min_alleles = p["min_alleles"] if p.get("min_alleles") else len(low)         # paper: depth ~ bulk size
    depth = ref + alt
    keep = (maf >= p["maf"]) & (depth.min(axis=0) >= min_alleles)
    keep &= np.isfinite(maf)
    if keep.sum() < 20:
        raise EtGwasError(f"Only {int(keep.sum())} SNPs pass the MAF ({p['maf']}) and depth ({min_alleles} alleles per bulk) filters.")

    res = trend_test(ref[:, keep], alt[:, keep])
    stat = res["stat"]
    lam = estimate_lambda(stat)
    lam_used = lam if p["use_gc"] else 1.0
    pval = stats.chi2.sf(stat / lam_used, 1)
    qval = bh_qvalues(pval)

    # orientation: the counted allele is the "alt" one; the GLM models the reference allele, so flip the sign
    slope_alt = -res["slope"]
    af = alt[:, keep] / np.maximum(depth[:, keep], 1)
    n_lines = len(low) + len(rand) + len(high)
    tbl = pd.DataFrame({
        "SNP": np.asarray(snp_names)[keep], "Chr": np.asarray(chrom)[keep].astype(int), "Position": np.asarray(pos)[keep].astype(np.int64),
        "MAF": maf[keep],
        "AF_low": af[0], "AF_random": af[1], "AF_high": af[2], "Delta_AF_high_minus_low": af[2] - af[0],
        "Alleles_low": depth[0, keep], "Alleles_random": depth[1, keep], "Alleles_high": depth[2, keep],
        "Chi2": stat, "P_value": pval, "Q_value": qval,
        "Slope_logodds": slope_alt, "SE": res["se"],
        "R2_deviance": np.where(res["null_dev"] > 0, stat / res["null_dev"], 0.0),
        "PVE_%": pve_formula(slope_alt, res["se"], maf[keep], n_lines),
    })
    tbl["Counted_allele_enriched_in"] = np.where(tbl["Delta_AF_high_minus_low"] >= 0, "high bulk", "low bulk")
    tbl = tbl.sort_values(["Chr", "Position"]).reset_index(drop=True)

    m_t = len(tbl)
    tbl["Suggestive"] = tbl["P_value"] < p["suggestive"]
    tbl["FDR_significant"] = tbl["Q_value"] < p["fdr"]
    tbl["Bonferroni_significant"] = tbl["P_value"] < 0.05 / m_t
    flag = tbl["FDR_significant"] if p["call_by"] == "fdr" else tbl["Suggestive"]
    tbl["MTA"] = flag
    regions = merge_regions(tbl, tbl["MTA"], int(p["merge_kb"] * 1000))
    return {"table": tbl, "regions": regions, "lambda": lam, "lambda_used": lam_used, "bulk_idx": bulk_idx, "info": info,
            "n_tested": m_t, "n_input": m, "min_alleles": int(min_alleles), "stat_all": stat, "keep": keep}


# ──────────────────────────── figures ────────────────────────────

CHR_COLORS = ["#2F6F9F", "#D98324"]


def _save(fig, path):
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    import matplotlib.pyplot as plt
    plt.close(fig)
    return str(path)


def _chrom_layout(tbl):
    chroms = sorted(tbl.Chr.unique())
    length = {c: float(tbl.loc[tbl.Chr == c, "Position"].max()) for c in chroms}
    gap = 0.012 * sum(length.values())
    off, cum = {}, 0.0
    for c in chroms:
        off[c] = cum
        cum += length[c] + gap
    x = tbl["Chr"].map(off).to_numpy() + tbl["Position"].to_numpy()
    ticks = [off[c] + length[c] / 2 for c in chroms]
    return x, chroms, ticks, cum - gap


def plot_bulks(y, res, trait):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 4.2))
    bins = np.histogram_bin_edges(y, bins=min(40, max(10, len(y) // 8)))
    ax.hist(y, bins=bins, color="#d9d9d9", edgecolor="white", label="whole panel")
    for name, col in (("low", "#c0392b"), ("random", "#2e9d5b"), ("high", "#2f6f9f")):
        ax.hist(y[res["bulk_idx"][name]], bins=bins, color=col, alpha=0.75, edgecolor="white",
                label=f"{name} bulk (n={len(res['bulk_idx'][name])})")
    ax.axvline(res["info"]["low_threshold"], color="#c0392b", ls="--", lw=1)
    ax.axvline(res["info"]["high_threshold"], color="#2f6f9f", ls="--", lw=1)
    ax.set_xlabel(trait); ax.set_ylabel("Number of lines"); ax.set_title("Phenotype distribution and the three bulks")
    ax.legend(frameon=False, fontsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return fig


def plot_manhattan(tbl, p):
    import matplotlib.pyplot as plt
    x, chroms, ticks, total = _chrom_layout(tbl)
    lp = -np.log10(np.clip(tbl["P_value"].to_numpy(), 1e-300, 1))
    fig, ax = plt.subplots(figsize=(13, 4.6))
    colors = np.array([CHR_COLORS[chroms.index(c) % 2] for c in tbl.Chr])
    ax.scatter(x, lp, s=5, c=colors, linewidths=0, alpha=0.75)
    mta = tbl["MTA"].to_numpy()
    ax.scatter(x[mta], lp[mta], s=22, facecolors="none", edgecolors="#c0392b", linewidths=1.0, label="MTA", zorder=3)
    sig_q = tbl.loc[tbl.Q_value < p["fdr"], "P_value"]
    if len(sig_q):
        ax.axhline(-np.log10(sig_q.max()), color="#c0392b", ls="--", lw=1, label=f"FDR {p['fdr']:g}")
    ax.axhline(-np.log10(p["suggestive"]), color="#e69f00", ls="-.", lw=1, label=f"suggestive (p<{p['suggestive']:g})")
    ax.axhline(-np.log10(0.05 / len(tbl)), color="#6a3d9a", ls=":", lw=1, label="Bonferroni (0.05/m)")
    ax.set_xticks(ticks); ax.set_xticklabels([str(c) for c in chroms], fontsize=8)
    ax.set_xlim(0, total); ax.set_xlabel("Chromosome"); ax.set_ylabel("-log10(P)")
    ax.set_title(f"Et-GWAS Manhattan plot  ({len(tbl):,} SNPs, inflation factor applied)")
    ax.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=0.95, fontsize=8, ncol=4,
              loc="upper center", bbox_to_anchor=(0.5, -0.16))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return fig


def plot_qq(stat, lam, use_gc):
    import matplotlib.pyplot as plt
    s = np.sort(stat[np.isfinite(stat)])
    k = len(s)
    exp = np.sort(stats.chi2.isf((np.arange(1, k + 1) - 0.5) / k, 1))
    step = max(1, k // 20000)
    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    ax.scatter(exp[::step], s[::step], s=6, color="#2f6f9f", linewidths=0, label="observed")
    if use_gc:
        ax.scatter(exp[::step], (s / lam)[::step], s=6, color="#2e9d5b", linewidths=0, label=f"after correction (λ={lam:.2f})")
    top = max(exp.max(), s.max())
    ax.plot([0, top], [0, top], color="black", lw=0.8)
    ax.set_xlabel("Expected χ² (1 d.f.)"); ax.set_ylabel("Observed χ²")
    ax.set_title(f"QQ plot  (inflation factor λ = {lam:.2f})")
    ax.legend(frameon=False, fontsize=8)
    for sd in ("top", "right"):
        ax.spines[sd].set_visible(False)
    return fig


def plot_allele_freq(tbl, n_top=12):
    import matplotlib.pyplot as plt
    top = tbl.sort_values("P_value").head(n_top)
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    cmap = plt.get_cmap("tab20")
    for i, (_, r) in enumerate(top.iterrows()):
        ax.plot([1, 2, 3], [r.AF_low, r.AF_random, r.AF_high], marker="o", lw=1.4, color=cmap(i % 20),
                label=f"{r.SNP}  (p={r.P_value:.1e})")
    ax.set_xticks([1, 2, 3]); ax.set_xticklabels(["low bulk", "random bulk", "high bulk"])
    ax.set_ylabel("Frequency of the counted allele"); ax.set_ylim(-0.02, 1.02)
    ax.set_title(f"Allele frequency across the bulks - the {len(top)} most significant SNPs")
    ax.legend(frameon=False, fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5))
    for sd in ("top", "right"):
        ax.spines[sd].set_visible(False)
    return fig


def plot_pca(G, res, rng_seed=1, max_snps=5000):
    import matplotlib.pyplot as plt
    n, m = G.shape
    rng = np.random.default_rng(rng_seed)
    cols = np.sort(rng.choice(m, size=min(m, max_snps), replace=False))
    X = G[:, cols].astype(float)
    mu = np.nanmean(X, axis=0)
    X = np.where(np.isnan(X), mu, X) - mu
    sd = X.std(axis=0)
    X = X[:, sd > 0] / sd[sd > 0]
    U, S, _ = np.linalg.svd(X, full_matrices=False)
    pc = U[:, :2] * S[:2]
    var = S ** 2 / (S ** 2).sum()
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    ax.scatter(pc[:, 0], pc[:, 1], s=14, color="#d0d0d0", label="other lines", linewidths=0)
    for name, col in (("low", "#c0392b"), ("random", "#2e9d5b"), ("high", "#2f6f9f")):
        idx = res["bulk_idx"][name]
        ax.scatter(pc[idx, 0], pc[idx, 1], s=26, color=col, alpha=0.85, label=f"{name} bulk", linewidths=0)
    ax.set_xlabel(f"PC1 ({var[0] * 100:.1f}%)"); ax.set_ylabel(f"PC2 ({var[1] * 100:.1f}%)")
    ax.set_title("Population structure of the lines in the bulks")
    ax.legend(frameon=False, fontsize=8)
    for sd_ in ("top", "right"):
        ax.spines[sd_].set_visible(False)
    return fig


# ──────────────────────────── the run ────────────────────────────

def _num(v, nd=None):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(v):
        return None
    return round(v, nd) if nd is not None else v


def _records(df, cols, limit=1000):
    out = []
    for r in df[cols].head(limit).to_dict("records"):
        out.append({k: (_num(v) if isinstance(v, (float, np.floating)) else (int(v) if isinstance(v, (np.integer,)) else
                                                                          (bool(v) if isinstance(v, (np.bool_,)) else v)))
                    for k, v in r.items()})
    return out


MTA_COLS = ["SNP", "Chr", "Position", "MAF", "AF_low", "AF_random", "AF_high", "Delta_AF_high_minus_low", "Chi2",
            "P_value", "Q_value", "Slope_logodds", "SE", "PVE_%", "Counted_allele_enriched_in", "Bonferroni_significant"]
REGION_COLS = ["Region", "Chr", "Start", "End", "Size_bp", "N_SNPs", "Peak_SNP", "Peak_Position", "Peak_P", "Enriched_in", "Max_PVE_%"]


def run_et_gwas(cfg, log=print):
    """cfg: an EtGwasConfig. Returns the result dict shown by the app (status 'success' or 'error')."""
    from tools.gwas.backend.gwas_inputs import InputError, prepare_input

    out_dir = Path(tempfile.mkdtemp(prefix="terratools_etgwas_"))
    try:
        log("[ET-GWAS] reading the genotypes and the phenotype...")
        path, in_info = prepare_input(cfg, str(out_dir), log=log)
    except InputError as e:
        return {"status": "error", "message": str(e)}
    except Exception as e:
        import traceback
        traceback.print_exc()
        return {"status": "error", "message": f"Could not prepare the input data: {e}"}

    try:
        from tools.gwas.backend.gwas_engine import _SimpleFileObj, parse_data
        df, snp_cols, chroms, positions = parse_data(file_obj=_SimpleFileObj(path), text_input=None,
                                                       encoding="dosage" if in_info else getattr(cfg, "genotype_encoding", "auto"))
        y = df["Phenotype"].to_numpy(dtype=float)
        ids = df["ID"].astype(str).to_numpy()
        G = df[snp_cols].to_numpy(dtype=float)
        del df
        p = {"bulk_pct": float(cfg.bulk_pct), "seed": int(cfg.seed), "maf": float(cfg.maf), "min_alleles": cfg.min_alleles,
             "use_gc": bool(cfg.use_gc), "fdr": float(cfg.fdr), "suggestive": float(cfg.suggestive),
             "call_by": cfg.call_by, "merge_kb": float(cfg.merge_kb)}
        log(f"[ET-GWAS] {len(y)} lines x {G.shape[1]:,} SNPs; building the bulks and testing every SNP...")
        res = analyse(G, y, ids, np.asarray(chroms), np.asarray(positions), np.asarray(snp_cols), p)
        tbl, regions = res["table"], res["regions"]

        # ---- bulk summary ----
        bulk_rows = []
        for name in BULKS:
            idx = res["bulk_idx"][name]
            bulk_rows.append({"Bulk": name, "Lines": len(idx), "Mean": _num(y[idx].mean(), 4), "Median": _num(np.median(y[idx]), 4),
                              "Min": _num(y[idx].min(), 4), "Max": _num(y[idx].max(), 4)})
        bulk_rows.append({"Bulk": "whole panel", "Lines": len(y), "Mean": _num(y.mean(), 4), "Median": _num(np.median(y), 4),
                          "Min": _num(y.min(), 4), "Max": _num(y.max(), 4)})
        bulk_membership = pd.DataFrame({"ID": ids, "Phenotype": y, "Bulk": ""})
        for name in BULKS:
            bulk_membership.loc[res["bulk_idx"][name], "Bulk"] = name

        # ---- figures ----
        import matplotlib
        matplotlib.use("Agg")
        plots = {}
        trait = str(getattr(cfg, "trait", "") or "Phenotype")
        figs = {
            "01_Phenotype_bulks": lambda: plot_bulks(y, res, trait),
            "02_Manhattan": lambda: plot_manhattan(tbl, p),
            "03_QQ": lambda: plot_qq(res["stat_all"], res["lambda"], p["use_gc"]),
            "04_Allele_frequency_top_SNPs": lambda: plot_allele_freq(tbl),
            "05_PCA_bulks": lambda: plot_pca(G, res, p["seed"]),
        }
        for name, make in figs.items():
            try:
                plots[name] = _save(make(), out_dir / f"{name}.png")
            except Exception as e:                              # one failing figure must not lose the analysis
                print(f"[ET-GWAS] figure {name} failed: {e}", flush=True)

        # ---- files ----
        mta = tbl[tbl["MTA"]].sort_values("P_value")
        all_csv = out_dir / "et_gwas_all_snps.csv"
        mta_csv = out_dir / "et_gwas_mtas.csv"
        bulks_csv = out_dir / "et_gwas_bulks.csv"
        tbl.drop(columns=["MTA"]).to_csv(all_csv, index=False, float_format="%.6g")
        mta[MTA_COLS].to_csv(mta_csv, index=False, float_format="%.6g")
        bulk_membership.to_csv(bulks_csv, index=False)

        crit = (f"FDR (q) < {p['fdr']:g}" if p["call_by"] == "fdr" else f"p < {p['suggestive']:g} (suggestive)")
        notes = pd.DataFrame([
            ("Method", "Et-GWAS: Gnanapragasam et al. 2024, Life Sci Alliance, doi:10.26508/lsa.202302352 (statistics follow the authors' code)."),
            ("Lines", f"{len(y)} lines with a phenotype and genotypes; bulk size {p['bulk_pct']:g}% = {len(res['bulk_idx']['low'])} lines per bulk (ties kept)."),
            ("Bulks", "low = lowest, high = highest phenotypes; random = lines drawn at random from the lines between the extremes (seed " + str(p["seed"]) + ")."),
            ("SNPs", f"{res['n_input']:,} read, {res['n_tested']:,} tested (MAF >= {p['maf']}, at least {res['min_alleles']} called alleles in each bulk)."),
            ("Test", "Binomial GLM (logit) of the pooled ref/alt allele counts on the bulk order low=1, random=2, high=3 (1 pseudo-count added); the statistic is the deviance explained by the trend (1 d.f.)."),
            ("Inflation factor", f"lambda = {res['lambda']:.3f} " + ("(statistics divided by it before the p-value)" if p["use_gc"] else "(estimated, NOT applied)")),
            ("MTA rule", crit + "; regions merge MTAs closer than " + f"{p['merge_kb']:g} kb."),
            ("Counted allele", "The allele whose dosage is counted (the alternative allele); 'Counted_allele_enriched_in' says in which extreme bulk it is more frequent."),
            ("PVE_%", "2 b^2 f(1-f) / (2 b^2 f(1-f) + 2 se(b)^2 N f(1-f)); b = slope (log-odds per bulk step), f = MAF, N = lines in the three bulks."),
            ("R2_deviance", "Share of the null deviance explained by the linear trend (TerraTools' definition)."),
        ], columns=["Item", "Details"])
        xlsx = out_dir / "et_gwas_results.xlsx"
        summary_row = pd.DataFrame([{"Lines": len(y), "Bulk_size_%": p["bulk_pct"], "Lines_per_bulk": len(res["bulk_idx"]["low"]),
                                     "SNPs_tested": res["n_tested"], "Lambda": round(res["lambda"], 4),
                                     "MTAs": int(len(mta)), "Regions": int(len(regions)),
                                     "FDR_significant": int(tbl.FDR_significant.sum()), "Suggestive": int(tbl.Suggestive.sum()),
                                     "Bonferroni_significant": int(tbl.Bonferroni_significant.sum())}])
        with pd.ExcelWriter(xlsx, engine="openpyxl") as xw:
            summary_row.to_excel(xw, sheet_name="Summary", index=False)
            mta[MTA_COLS].to_excel(xw, sheet_name="MTAs", index=False)
            (regions if len(regions) else pd.DataFrame(columns=REGION_COLS)).to_excel(xw, sheet_name="Regions", index=False)
            pd.DataFrame(bulk_rows).to_excel(xw, sheet_name="Bulks", index=False)
            bulk_membership.to_excel(xw, sheet_name="Bulk_members", index=False)
            tbl.drop(columns=["MTA"]).replace([np.inf, -np.inf], np.nan).to_excel(xw, sheet_name="All_SNPs", index=False)
            notes.to_excel(xw, sheet_name="Notes", index=False)

        zip_path = out_dir / "et_gwas_plots.zip"
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for name, f in plots.items():
                zf.write(f, f"{name}.png")

        summary = "\n".join([
            f"Lines: {len(y)}   bulks of {p['bulk_pct']:g}% = {len(res['bulk_idx']['low'])} lines each (low / random / high)",
            f"SNPs tested: {res['n_tested']:,} of {res['n_input']:,}   inflation factor lambda = {res['lambda']:.3f}"
            + ("" if p["use_gc"] else " (not applied)"),
            f"MTAs ({crit}): {len(mta)} in {len(regions)} region(s)",
            f"FDR < {p['fdr']:g}: {int(tbl.FDR_significant.sum())}   suggestive (p < {p['suggestive']:g}): {int(tbl.Suggestive.sum())}   "
            f"Bonferroni: {int(tbl.Bonferroni_significant.sum())}",
        ])
        if in_info:
            summary = f"Input: {in_info['source']} - trait '{in_info['trait']}', {in_info['matched']} of {in_info['phenotype_ids']} phenotype IDs matched\n" + summary
        return {
            "status": "success", "output_dir": str(out_dir), "summary": summary,
            "lambda": _num(res["lambda"], 4), "n_lines": int(len(y)), "n_snps_tested": int(res["n_tested"]),
            "n_mtas": int(len(mta)), "n_regions": int(len(regions)), "min_alleles": res["min_alleles"],
            "bulk_table": bulk_rows,
            "mta_table": _records(mta, MTA_COLS),
            "region_table": _records(regions, REGION_COLS) if len(regions) else [],
            "plots": plots, "excel": str(xlsx), "mta_csv": str(mta_csv), "all_csv": str(all_csv),
            "bulks_csv": str(bulks_csv), "zip": str(zip_path),
        }
    except EtGwasError as e:
        return {"status": "error", "message": str(e)}
    except Exception as e:
        import traceback
        traceback.print_exc()
        return {"status": "error", "message": f"Et-GWAS failed: {e}"}
