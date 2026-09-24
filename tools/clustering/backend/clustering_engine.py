"""
clustering_engine.py — Genotypic Data Clustering / Genomic Selection tool.

Same logic as the Gradio app's cluster_split_data(): PCA -> K-Means
clustering -> stratified train/test split. No math has been changed —
only restructured to write files to a given output directory instead of
returning Gradio component values.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless rendering
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.model_selection import StratifiedShuffleSplit


def cluster_split_data(geno_file: str, n_clusters: int, test_size: float, outdir: Path) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)

    # ── 1. Load data ────────────────────────────────────────────────────
    try:
        geno = pd.read_csv(geno_file, sep="\t")
        if geno.shape[1] == 1:  # probably comma-separated after all
            geno = pd.read_csv(geno_file, sep=",")
    except Exception as e:
        raise ValueError(f"Could not read file: {e}")

    id_col = geno.columns[0]
    sample_ids = geno[id_col].astype(str)
    geno_data = geno.iloc[:, 1:].copy()

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

    kmeans = KMeans(n_clusters=n_clusters, random_state=60, n_init=20)
    clusters = kmeans.fit_predict(geno_reduced)

    geno["Cluster"] = clusters

    # ── 5. Plot ───────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 6))
    scatter = ax.scatter(pc1, pc2, c=clusters, cmap="viridis", alpha=0.7, s=80, edgecolors="k", linewidths=0.4)

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
    cluster_counts = pd.Series(clusters).value_counts()
    can_stratify = (cluster_counts >= 2).all() and (n_samples >= 4)

    if can_stratify:
        sss = StratifiedShuffleSplit(n_splits=1, test_size=test_size, random_state=42)
        train_idx, test_idx = next(sss.split(geno, geno["Cluster"]))
        split_method = "stratified"
    else:
        rng = np.random.default_rng(42)
        idx = rng.permutation(n_samples)
        split = max(1, int(n_samples * (1 - test_size)))
        train_idx, test_idx = idx[:split], idx[split:]
        split_method = "random (fallback — some clusters had <2 members)"

    geno_train_df = geno.iloc[train_idx].copy()
    geno_test_df = geno.iloc[test_idx].copy()

    train_path = outdir / "Training_data_Cluster.csv"
    test_path = outdir / "Testing_data_Cluster.csv"
    geno_train_df.to_csv(train_path, index=False)
    geno_test_df.to_csv(test_path, index=False)

    cluster_table = geno[[id_col, "Cluster"]].sort_values("Cluster").reset_index(drop=True)
    cluster_table_path = outdir / "cluster_membership.csv"
    cluster_table.to_csv(cluster_table_path, index=False)

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
        "plot_png": str(plot_path),
        "train_csv": str(train_path),
        "test_csv": str(test_path),
        "cluster_table_csv": str(cluster_table_path),
        "cluster_table": cluster_table.to_dict(orient="records"),
    }
