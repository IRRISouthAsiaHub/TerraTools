"""
clustering_engine.py — Genotypic Data Clustering / Genomic Selection tool.

PCA -> K-Means clustering -> stratified train/test split (the logic of the Gradio app's
cluster_split_data()), restructured to write files to a given output directory.

Optional phenotype file. When it is given, the
training file also gets a `phenotypes` column (sample_id, phenotypes, markers..., Cluster),
so it is a phenotyped training set. Lines without a phenotype cannot be
trained on, so they are put in the test set.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless rendering
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.model_selection import StratifiedShuffleSplit

# Sample IDs are written on the plot only when there are few enough to read.
MAX_LABELLED = 150
# K-selection aid: K = 2 ... MAX_K_TRIED; the silhouette is computed on at most this many lines
MAX_K_TRIED = 10
SILHOUETTE_SAMPLE = 2000


def read_table(path) -> pd.DataFrame:
    """CSV, TSV/TXT or Excel, decided by the extension and the first line."""
    p = str(path)
    if p.lower().endswith((".xlsx", ".xlsm", ".xls")):
        return pd.read_excel(p)
    with open(p, "r", encoding="utf-8-sig", errors="replace") as fh:
        first = fh.readline()
    return pd.read_csv(p, sep="\t" if "\t" in first else ",", encoding="utf-8-sig")


def load_phenotypes(path):
    """Two columns: sample IDs first, then the trait (a column called 'phenotypes' wins, else the second
    column). Returns (DataFrame[sample_id, phenotypes], column name used, counts)."""
    ph = read_table(path)
    if ph.shape[1] < 2:
        raise ValueError("The phenotype file needs at least two columns: sample IDs, then the trait values.")
    col = next((c for c in ph.columns[1:] if str(c).strip().lower() in ("phenotypes", "phenotype")), ph.columns[1])
    out = pd.DataFrame({"sample_id": ph.iloc[:, 0].astype(str).str.strip(),
                        "phenotypes": pd.to_numeric(ph[col], errors="coerce")})
    n_rows = len(out)
    out = out.dropna(subset=["phenotypes"])
    n_blank = n_rows - len(out)
    n_dup = int(out["sample_id"].duplicated().sum())
    out = out.drop_duplicates("sample_id", keep="first")
    return out, str(col), {"rows": n_rows, "blank": n_blank, "duplicates": n_dup}


def k_selection(Z: np.ndarray, outdir: Path, chosen_k: int):
    """Silhouette width (higher = better separated clusters) and within-cluster sum of squares (elbow) for K = 2 ... 10,
    on the space the clustering is done in. Returns (table, suggested K by silhouette, png path) or (None, None, None)."""
    n = Z.shape[0]
    ks = list(range(2, min(MAX_K_TRIED, n - 1) + 1))
    if len(ks) < 2:
        return None, None, None
    rng = np.random.default_rng(1)
    idx = rng.choice(n, SILHOUETTE_SAMPLE, replace=False) if n > SILHOUETTE_SAMPLE else np.arange(n)
    rows = []
    for k in ks:
        km = KMeans(n_clusters=k, random_state=60, n_init=10).fit(Z)
        lab = km.labels_
        sil = float(silhouette_score(Z[idx], lab[idx])) if len(set(lab[idx])) > 1 else float("nan")
        rows.append({"K": k, "Silhouette": round(sil, 4), "Within_SS": round(float(km.inertia_), 2)})
    tab = pd.DataFrame(rows)
    ok = tab.dropna(subset=["Silhouette"])
    best = int(ok.loc[ok["Silhouette"].idxmax(), "K"]) if len(ok) else None
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.8))
    a1.plot(tab["K"], tab["Silhouette"], "o-", color="#2E6B35")
    a1.set_xlabel("Number of clusters K"); a1.set_ylabel("Mean silhouette width"); a1.set_title("Separation (higher is better)")
    a2.plot(tab["K"], tab["Within_SS"], "o-", color="#6c63ff")
    a2.set_xlabel("Number of clusters K"); a2.set_ylabel("Within-cluster sum of squares"); a2.set_title("Elbow")
    for a in (a1, a2):
        if chosen_k in ks:
            a.axvline(chosen_k, color="#CC0000", ls="--", lw=1)
        if best:
            a.axvline(best, color="#d9a52b", ls=":", lw=1.5)
        a.grid(alpha=0.25)
    fig.suptitle(f"Choice of K: used {chosen_k} (red), highest silhouette {best} (gold)", fontsize=10)
    plt.tight_layout()
    path = outdir / "k_selection.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return tab, best, path


def cluster_split_data(geno_file: str, n_clusters: int, test_size: float, outdir: Path,
                       pheno_file: str = "", keep_test_pheno: bool = False, n_pcs: int = 0) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    notes = []

    # ── 1. Load data ────────────────────────────────────────────────────
    try:
        geno = read_table(geno_file)
    except Exception as e:
        raise ValueError(f"Could not read file: {e}")

    id_col = geno.columns[0]
    sample_ids = geno[id_col].astype(str)
    keys = sample_ids.str.strip()
    geno_data = geno.iloc[:, 1:].copy()

    # ── 1b. Phenotypes (optional) ───────────────────────────────────────
    pmap, pheno_col, has = None, None, None
    pheno_counts = {}
    if pheno_file:
        ph, pheno_col, c = load_phenotypes(pheno_file)
        pmap = ph.set_index("sample_id")["phenotypes"]
        has = keys.isin(pmap.index).to_numpy()
        pheno_counts = {
            "n_pheno_matched": int(has.sum()),
            "n_geno_without_pheno": int((~has).sum()),
            "n_pheno_without_geno": int((~pmap.index.isin(set(keys))).sum()),
        }
        if c["duplicates"]:
            notes.append(f"{c['duplicates']} duplicated sample ID(s) in the phenotype file: the first value of each was used.")
        if c["blank"]:
            notes.append(f"{c['blank']} line(s) in the phenotype file have no numeric value and were ignored.")
        if has.sum() < 4:
            raise ValueError(
                f"Only {int(has.sum())} line(s) of the genotype file have a phenotype. Check that the sample IDs of "
                "the two files are written the same way.")
        if (~has).any():
            notes.append(f"{int((~has).sum())} genotyped line(s) have no phenotype: they are all placed in the test set.")

    # ── 2. Preprocess ────────────────────────────────────────────────────
    geno_data.replace(-1, np.nan, inplace=True)
    geno_data = geno_data.apply(
        lambda col: col.fillna(col.mode()[0] if not col.mode().empty else 0)
    )
    geno_data = geno_data.apply(pd.to_numeric, errors="coerce").fillna(0)

    n_samples, n_features = geno_data.shape

    # ── 3. PCA ───────────────────────────────────────────────────────────
    max_components = min(n_samples, n_features)
    pca_components = max(2, min(max_components, max_components))

    pca = PCA(n_components=pca_components, random_state=42)
    geno_reduced = pca.fit_transform(geno_data)

    pc1 = geno_reduced[:, 0]
    pc2 = geno_reduced[:, 1] if geno_reduced.shape[1] > 1 else np.zeros(n_samples)

    explained = pca.explained_variance_ratio_
    pc1_var = round(explained[0] * 100, 1)
    pc2_var = round(explained[1] * 100, 1) if len(explained) > 1 else 0.0

    # ── 4. Clustering ────────────────────────────────────────────────────
    n_clusters = int(n_clusters)
    n_clusters = max(2, min(n_clusters, n_samples - 1))

    # K-means on the first n_pcs principal components (0 = all of them, which is clustering in the full genotype space)
    n_use = geno_reduced.shape[1] if not n_pcs else max(2, min(int(n_pcs), geno_reduced.shape[1]))
    space = geno_reduced[:, :n_use]
    kmeans = KMeans(n_clusters=n_clusters, random_state=60, n_init=20)
    clusters = kmeans.fit_predict(space)
    try:
        k_tab, k_best, k_png = k_selection(space, outdir, n_clusters)
    except Exception as e:                      # an aid only: never stops the run
        k_tab, k_best, k_png = None, None, None
        notes.append(f"The K-selection plot could not be made: {e}")
    if n_use < geno_reduced.shape[1]:
        notes.append(f"K-means used the first {n_use} of {geno_reduced.shape[1]} principal components.")

    geno["Cluster"] = clusters

    # ── 5. Plot ───────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 6))
    scatter = ax.scatter(pc1, pc2, c=clusters, cmap="viridis", alpha=0.7, s=80 if n_samples <= MAX_LABELLED else 22,
                         edgecolors="k", linewidths=0.4)

    if n_samples <= MAX_LABELLED:
        for i, sid in enumerate(sample_ids):
            ax.annotate(sid, (pc1[i], pc2[i]), fontsize=7, ha="left", va="bottom", alpha=0.75)

    ax.set_xlabel(f"PC1 ({pc1_var}% variance)")
    ax.set_ylabel(f"PC2 ({pc2_var}% variance)")
    ax.set_title(f"Genotypic Clusters (K={n_clusters})")
    plt.colorbar(scatter, ax=ax, label="Cluster")
    plt.tight_layout()

    plot_path = outdir / "cluster_plot.png"
    fig.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ── 6. Train / Test split ────────────────────────────────────────────
    # With phenotypes only the lines that have one can be split; the others go to the test set.
    all_idx = np.arange(n_samples)
    elig = all_idx[has] if has is not None else all_idx
    extra_test = all_idx[~has] if has is not None else np.array([], dtype=int)

    cl = clusters[elig]
    cluster_counts = pd.Series(cl).value_counts()
    train_idx = test_idx = None
    split_method = ""
    if (cluster_counts >= 2).all() and len(elig) >= 4:
        try:
            sss = StratifiedShuffleSplit(n_splits=1, test_size=test_size, random_state=42)
            a, b = next(sss.split(elig, cl))
            train_idx, test_idx, split_method = elig[a], elig[b], "stratified"
        except ValueError:
            pass            # too few lines per cluster for this test fraction
    if train_idx is None:
        rng = np.random.default_rng(42)
        perm = rng.permutation(elig)
        cut = max(1, int(len(elig) * (1 - test_size)))
        train_idx, test_idx = perm[:cut], perm[cut:]
        split_method = "random (fallback — some clusters had <2 members)"
    if len(extra_test):
        test_idx = np.concatenate([test_idx, extra_test])

    geno_train_df = geno.iloc[train_idx].copy()
    geno_test_df = geno.iloc[test_idx].copy()

    test_has_pheno = False
    if pmap is not None:
        for d in (geno_train_df, geno_test_df):
            if "phenotypes" in d.columns:
                d.drop(columns=["phenotypes"], inplace=True)
        geno_train_df.insert(1, "phenotypes", keys.iloc[train_idx].map(pmap).to_numpy())
        if keep_test_pheno and len(extra_test) == 0:
            geno_test_df.insert(1, "phenotypes", keys.iloc[test_idx].map(pmap).to_numpy())
            test_has_pheno = True
        elif keep_test_pheno:
            notes.append("The test file has no phenotypes column, because some test lines have no phenotype "
                         "(a test set is either fully phenotyped or not at all).")

    train_path = outdir / "Training_data_Cluster.csv"
    test_path = outdir / "Testing_data_Cluster.csv"
    geno_train_df.to_csv(train_path, index=False)
    geno_test_df.to_csv(test_path, index=False)

    # ── 7. Membership table ──────────────────────────────────────────────
    which = np.full(n_samples, "", dtype=object)
    which[train_idx] = "Train"
    which[test_idx] = "Test"
    table = pd.DataFrame({id_col: sample_ids.to_numpy(), "Cluster": clusters, "Set": which})
    if pmap is not None:
        table["Phenotype"] = keys.map(pmap).to_numpy()
    cluster_table = table.sort_values("Cluster", kind="stable").reset_index(drop=True)
    cluster_table_path = outdir / "cluster_membership.csv"
    cluster_table.to_csv(cluster_table_path, index=False)

    # JSON-safe records (NaN -> None) with the columns in the order written above
    cols = list(cluster_table.columns)
    records = [{c: (None if (isinstance(v, float) and np.isnan(v)) else (v.item() if hasattr(v, "item") else v))
                for c, v in zip(cols, row)} for row in cluster_table.itertuples(index=False, name=None)]

    return {
        "status": "success",
        "n_samples": int(n_samples),
        "n_features": int(n_features),
        "n_clusters_used": int(n_clusters),
        "n_train": int(len(train_idx)),
        "n_test": int(len(test_idx)),
        "split_method": split_method,
        "pc1_variance": pc1_var,
        "pc2_variance": pc2_var,
        "pheno_used": pmap is not None,
        "pheno_column": pheno_col,
        "train_has_pheno": pmap is not None,
        "test_has_pheno": test_has_pheno,
        **pheno_counts,
        "notes": notes,
        "plot_png": str(plot_path),
        "n_pcs_used": int(n_use),
        "n_pcs_available": int(geno_reduced.shape[1]),
        "k_selection_png": str(k_png) if k_png else "",
        "k_selection": [] if k_tab is None else k_tab.to_dict("records"),
        "k_suggested": k_best,
        "train_csv": str(train_path),
        "test_csv": str(test_path),
        "cluster_table_csv": str(cluster_table_path),
        "cluster_table": records,
    }
