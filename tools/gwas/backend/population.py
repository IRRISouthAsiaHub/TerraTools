"""
Relationship matrices, LD pruning and principal components for GWAS.

Kinship (genomic relationship matrix)
    Z  = G - 2p                               G = 0/1/2 dosage, p = frequency of the counted allele (observed calls only)
    K  = Z Z' / (2 * sum p (1 - p))           so that the average diagonal is about 1 for outbred material
  - a missing call contributes 0 to Z (it is the mean), and p never comes from imputed values;
  - K is symmetric by construction; it is positive semi-definite up to rounding. `eig_psd` reports the smallest eigenvalue
    and only sets eigenvalues that are negative by rounding noise to 0 (nothing is added to the diagonal, nothing is rescaled);
  - the number of markers that entered K is returned with it.
  Variants: all QC-passed markers ("vanraden"), LD-pruned markers ("ld_pruned"), and leave-one-chromosome-out matrices
  (the markers of the chromosome being tested are left out of K, which avoids testing a marker against a relationship matrix
  that contains it).

Principal components
    z = (g - 2p) / sqrt(2p(1-p)) per marker, missing = 0, the sample-by-sample matrix is divided by the number of markers
    observed in both samples. LD-pruned markers are used by default because LD blocks would otherwise drive the leading
    components. Eigenvalues and the share of variance explained are returned for a scree plot.
"""

from __future__ import annotations

import numpy as np

KINSHIP_METHODS = ("vanraden", "ld_pruned")


# ============================================================================================ kinship

def _freq(G):
    with np.errstate(invalid="ignore"):
        p = np.nanmean(G, axis=0) / 2.0
    return np.nan_to_num(p, nan=0.0)


def centred(G):
    """Z = G - 2p with missing calls set to 0."""
    G = np.asarray(G, dtype=float)
    p = _freq(G)
    Z = G - 2.0 * p
    return np.where(np.isnan(Z), 0.0, Z), p


def vanraden(G):
    """(K, info). K = Z Z' / (2 sum p(1-p)) over the markers in G."""
    Z, p = centred(G)
    d = 2.0 * float(np.sum(p * (1.0 - p)))
    if d <= 0:
        raise ValueError("No marker varies: the relationship matrix cannot be computed.")
    K = Z @ Z.T / d
    K = (K + K.T) / 2.0
    info = {"method": "VanRaden (centred on twice the allele frequency, scaled by 2 sum p(1-p))", "n_snps": int(G.shape[1]),
            "mean_diagonal": float(np.mean(np.diag(K))),
            "mean_off_diagonal": float((K.sum() - np.trace(K)) / max(K.shape[0] * (K.shape[0] - 1), 1))}
    return K, info


def eig_psd(K, tol=1e-8):
    """Eigendecomposition of a kinship matrix: (eigenvalues ascending, eigenvectors, info).

    A symmetric matrix of the form Z Z' has no truly negative eigenvalues; those found are rounding noise and are set to 0.
    If a clearly negative eigenvalue (below -tol times the largest) is found the matrix is not positive semi-definite and
    info['psd'] is False, so the caller can warn; the eigenvalue is still not altered beyond the 0 floor."""
    Ks = (np.asarray(K) + np.asarray(K).T) / 2.0
    w, V = np.linalg.eigh(Ks)
    scale = max(float(w.max()), 1e-300)
    psd = bool(w.min() >= -tol * scale)
    info = {"min_eigenvalue": float(w.min()), "max_eigenvalue": float(w.max()), "psd": psd,
            "n_negative_set_to_zero": int((w < 0).sum()), "rank": int((w > tol * scale).sum())}
    return np.clip(w, 0.0, None), V, info


def kinship(G, chrom=None, pos=None, method="vanraden", ld_r2=0.2, ld_window=50, ld_step=5):
    """Relationship matrix and a description. method 'vanraden' uses every marker, 'ld_pruned' the LD-pruned ones."""
    if method not in KINSHIP_METHODS:
        raise ValueError(f"Unknown kinship method '{method}'.")
    cols = None
    if method == "ld_pruned":
        if chrom is None or pos is None:
            raise ValueError("LD pruning needs the chromosome and position of every marker.")
        cols = ld_prune(G, chrom, pos, r2=ld_r2, window=ld_window, step=ld_step)
        if len(cols) < 2:
            raise ValueError("Fewer than two markers remain after LD pruning; use a higher r² threshold or the full kinship.")
    K, info = vanraden(G if cols is None else G[:, cols])
    info["kinship_method"] = method
    info["n_snps_available"] = int(G.shape[1])
    if cols is not None:
        info["ld_prune"] = {"r2": ld_r2, "window_snps": ld_window, "step": ld_step}
    return K, info, cols


def loco_kinship(G, chrom, cols=None):
    """Iterator of (chromosome, K without that chromosome, info) for every chromosome that has markers.

    Uses S = Z Z' over all markers and subtracts the part of each chromosome, so the cost is one n x n product per
    chromosome. `cols` restricts the markers that make up K (e.g. the LD-pruned ones)."""
    G = np.asarray(G, dtype=float)
    cols = np.arange(G.shape[1]) if cols is None else np.asarray(cols)
    Z, p = centred(G[:, cols])
    ch = np.asarray(chrom)[cols]
    w = p * (1.0 - p)
    S = Z @ Z.T
    d = 2.0 * float(w.sum())
    for c in np.unique(ch):
        sel = ch == c
        Sc = Z[:, sel] @ Z[:, sel].T
        dc = 2.0 * float(w[sel].sum())
        if d - dc <= 0:
            continue                                   # a single-chromosome genome: nothing is left to build K from
        K = (S - Sc) / (d - dc)
        yield int(c), (K + K.T) / 2.0, {"n_snps": int((~sel).sum()), "left_out_chromosome": int(c)}


# ============================================================================================ LD pruning

def ld_prune(G, chrom, pos, r2=0.2, window=50, step=5):
    """Indices of markers that remain after greedy LD pruning (map order, window of `window` markers moved by `step`):
    a marker is dropped when its squared correlation with an earlier kept marker in the window exceeds `r2`."""
    G = np.asarray(G, dtype=float)
    chrom, pos = np.asarray(chrom), np.asarray(pos)
    X = np.where(np.isnan(G), np.nanmean(G, axis=0, keepdims=True), G)
    X = np.nan_to_num(X)
    keep_all = []
    for c in np.unique(chrom):
        idx = np.flatnonzero(chrom == c)
        idx = idx[np.argsort(pos[idx], kind="stable")]
        sd = X[:, idx].std(axis=0)
        idx = idx[sd > 1e-9]
        if len(idx) == 0:
            continue
        Z = X[:, idx] - X[:, idx].mean(axis=0)
        Z /= np.sqrt((Z ** 2).sum(axis=0))
        alive = np.ones(len(idx), dtype=bool)
        for s in range(0, len(idx), max(step, 1)):
            win = np.flatnonzero(alive[s:s + window]) + s
            if len(win) < 2:
                continue
            R2 = (Z[:, win].T @ Z[:, win]) ** 2
            kept = np.ones(len(win), dtype=bool)
            for i in range(1, len(win)):
                if (kept[:i] & (R2[:i, i] > r2)).any():
                    kept[i] = False
            alive[win[~kept]] = False
        keep_all.append(idx[alive])
    return np.sort(np.concatenate(keep_all)) if keep_all else np.zeros(0, dtype=int)


# ============================================================================================ principal components

def pca(G_observed, n_pc=10, chrom=None, pos=None, ld_prune_markers=True, ld_r2=0.2, block=4000):
    """Principal components of the standardised genotypes.

    G_observed keeps NaN for missing calls. Returns a dict with 'scores' (samples x n_pc, eigenvector scaled to unit norm),
    'eigenvalues' (all, descending), 'variance_explained' (fractions of the sum of the positive eigenvalues), 'n_snps' and the
    marker indices used."""
    G = np.asarray(G_observed, dtype=float)
    n, m = G.shape
    cols = np.arange(m)
    pruned = False
    if ld_prune_markers and chrom is not None and pos is not None and m >= 100:
        kept = ld_prune(np.where(np.isnan(G), np.nanmean(G, axis=0, keepdims=True), G), chrom, pos, r2=ld_r2)
        if len(kept) >= max(n_pc + 1, 50):
            cols, pruned = kept, True
    Gs = G[:, cols]
    p = _freq(Gs)
    ok = np.where((p > 0) & (p < 1))[0]
    if len(ok) < 2:
        raise ValueError("Fewer than two variable markers are available for the principal components.")
    S = np.zeros((n, n))
    C = np.zeros((n, n))
    for i in range(0, len(ok), block):
        j = ok[i:i + block]
        Z = (Gs[:, j] - 2 * p[j]) / np.sqrt(2 * p[j] * (1 - p[j]))
        obs = (~np.isnan(Z)).astype(float)
        Z = np.where(np.isnan(Z), 0.0, Z)
        S += Z @ Z.T
        C += obs @ obs.T
    M = S / np.maximum(C, 1.0)
    w, V = np.linalg.eigh((M + M.T) / 2.0)
    order = np.argsort(w)[::-1]
    w, V = w[order], V[:, order]
    pos_w = np.clip(w, 0, None)
    k = int(min(max(n_pc, 0), n))
    return {"scores": V[:, :k], "eigenvalues": w, "variance_explained": pos_w / max(pos_w.sum(), 1e-300),
            "n_snps": int(len(ok)), "ld_pruned": pruned, "marker_idx": cols[ok]}


def structure_summary(pc, K=None, labels=None):
    """Plain description of population structure and relatedness for the report.

    pc: output of pca(). labels: optional subpopulation label per sample (e.g. indica / japonica / aus); without labels the
    samples are grouped by k-means on the leading components (k = 3), which is stated as such."""
    ve = pc["variance_explained"]
    out = {"pc1_pc2_variance": float(ve[:2].sum()), "pc_variance": [float(x) for x in ve[:10]]}
    S = pc["scores"]
    if S.shape[1] >= 2:
        grp, source = None, ""
        if labels is not None:
            lab = np.asarray(labels).astype(str)
            if len(set(lab)) > 1:
                grp, source = lab, "labels supplied with the data"
        if grp is None and S.shape[0] >= 12:
            import os
            os.environ.setdefault("LOKY_MAX_CPU_COUNT", str(os.cpu_count() or 1))
            from sklearn.cluster import KMeans
            km = KMeans(n_clusters=3, n_init=10, random_state=1).fit(S[:, :min(S.shape[1], 3)])
            grp, source = np.array([f"group{c + 1}" for c in km.labels_]), "k-means (k = 3) on the leading components, not a classification"
        if grp is not None:
            r2 = []
            for k in range(min(S.shape[1], 3)):
                tot = ((S[:, k] - S[:, k].mean()) ** 2).sum()
                within = sum(((S[grp == g, k] - S[grp == g, k].mean()) ** 2).sum() for g in np.unique(grp))
                r2.append(float(1 - within / tot) if tot > 0 else 0.0)
            out["groups"] = {"source": source, "sizes": {str(g): int((grp == g).sum()) for g in np.unique(grp)},
                             "pc_variance_between_groups": r2, "assignment": np.asarray(grp)}
    if K is not None:
        off = K[np.triu_indices_from(K, 1)]
        out["relatedness"] = {"pairs_above_0.25": int((off > 0.25).sum()), "pairs_above_0.5": int((off > 0.5).sum()),
                              "max_pair": float(off.max()), "mean_pair": float(off.mean())}
    return out
