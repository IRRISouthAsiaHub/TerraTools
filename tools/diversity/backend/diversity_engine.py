"""
diversity_engine.py — Haplotypic diversity analysis engine.

The diversity metrics, PCA, clustering and the custom neighbour-joining
implementation are the same computations as before. What changed is the output:
every figure is now a PNG on a white background (diversity_plots.py), so the app
shows it at its natural proportions, opens it in a pop-up and saves it as an
image. The tree is a circular iTOL-style tree with bootstrap support.
"""

import re
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from tools.diversity.backend import diversity_plots as P

warnings.filterwarnings('ignore')

# =============================================================================
# DIVERSITY METRICS COMPUTATION
# =============================================================================
#
# Haplotype block values are CODES for haplotypes, not abundances. Each
# sample's metrics describe how its non-missing blocks are distributed over
# those codes:
#   - a block coded -1 (or left blank) is MISSING and is left out
#   - 0 is a real haplotype code and counts like any other
#   - p_i = share of the sample's non-missing blocks that carry code i
#
# This follows the corrected manual script (app_ENHANCED_NonNormalized.py).

METRIC_COLUMNS = [
    'Richness', 'Load', 'Shannon_Diversity', 'Simpson_Diversity', 'Gini_Simpson',
    'Effective_Alleles', 'Heterozygosity', 'Gene_Diversity', 'PIC_Index',
    'Nucleotide_Diversity', 'Allelic_Entropy', 'Richness_Normalized',
]


def compute_diversity_metrics(df, block_cols):
    n = len(df)
    # -1 = missing block call
    values = df[block_cols].replace(-1, np.nan).to_numpy(dtype=float)
    out = {name: np.zeros(n) for name in METRIC_COLUMNS}

    out['Richness'] = (~np.isnan(values)).sum(axis=1).astype(float)
    out['Load'] = np.nansum(values, axis=1)          # missing counts as 0

    for i in range(n):
        valid = values[i][~np.isnan(values[i])]
        if valid.size == 0:
            out['Effective_Alleles'][i] = 1.0
            continue

        _, counts = np.unique(valid, return_counts=True)
        probs = counts / counts.sum()
        sum_sq = float(np.sum(probs ** 2))

        out['Shannon_Diversity'][i] = -np.sum(probs * np.log(probs + 1e-10))
        out['Simpson_Diversity'][i] = 1 - sum_sq
        out['Gene_Diversity'][i] = 1 - sum_sq
        out['Gini_Simpson'][i] = 1 - sum_sq
        out['Effective_Alleles'][i] = 1 / sum_sq

        alt_freq = valid.sum() / (2 * valid.size)
        he = 2 * alt_freq * (1 - alt_freq)
        out['Heterozygosity'][i] = he if not np.isnan(he) else 0

        max_freq = probs.max()
        min_freq = 1 - max_freq
        pic = 2 * max_freq * min_freq * (1 - 2 * max_freq * min_freq)
        out['PIC_Index'][i] = max(0.0, min(pic, 1.0))

        # Extra columns (not in the manual script), on the same cleaned data
        out['Allelic_Entropy'][i] = -np.sum(probs * np.log2(probs + 1e-10))
        out['Nucleotide_Diversity'][i] = (
            float(np.mean(np.abs(np.diff(valid)))) if valid.size > 1 else 0.0
        )

    max_richness = len(block_cols)
    out['Richness_Normalized'] = out['Richness'] / max_richness if max_richness > 0 else np.zeros(n)

    return pd.DataFrame(out, index=df.index)



# =============================================================================
# NEIGHBOUR JOINING (same algorithm, the pair search is vectorised)
# =============================================================================

def neighbour_joining_simple(dist_matrix):
    """Neighbour Joining — returns list of (parent, child, branch_length) edges."""
    n = dist_matrix.shape[0]
    clusters = list(range(n))
    D = dist_matrix.copy().astype(float)
    edges = []
    next_id = n

    while len(clusters) > 2:
        m = len(clusters)
        row_sums = D.sum(axis=1)
        Q = (m - 2) * D - row_sums[:, None] - row_sums[None, :]
        np.fill_diagonal(Q, np.inf)
        flat_idx = np.argmin(Q)
        i, j = divmod(flat_idx, m)

        li = 0.5 * (D[i, j] + (row_sums[i] - row_sums[j]) / (m - 2))
        lj = D[i, j] - li

        new_node = next_id
        next_id += 1
        edges.append((new_node, clusters[i], max(li, 0)))
        edges.append((new_node, clusters[j], max(lj, 0)))

        old = [k for k in range(m) if k != i and k != j]
        new_dists = (D[i, old] + D[j, old] - D[i, j]) / 2

        size = len(old) + 1
        new_D = np.zeros((size, size))
        new_D[:-1, :-1] = D[np.ix_(old, old)]
        new_D[:-1, -1] = new_dists
        new_D[-1, :-1] = new_dists

        D = new_D
        clusters = [clusters[k] for k in old] + [new_node]

    if len(clusters) == 2:
        edges.append((clusters[0], clusters[1], D[0, 1] / 2))
        edges.append((clusters[1], clusters[0], D[0, 1] / 2))

    return edges


def nj_structure(dist_matrix):
    """Neighbour-joining tree as a rooted structure: children, post-order, leaf order, distance from the root."""
    edges = neighbour_joining_simple(dist_matrix)
    adj = defaultdict(list)
    all_nodes = set()
    for parent, child, dist in edges:
        adj[parent].append((child, dist))
        adj[child].append((parent, dist))
        all_nodes.add(parent)
        all_nodes.add(child)
    root = max(all_nodes)

    children = defaultdict(list)
    visited = {root}
    stack = [root]
    while stack:                                   # iterative, so a long chain cannot hit the recursion limit
        node = stack.pop()
        for nb, d in adj[node]:
            if nb not in visited:
                visited.add(nb)
                children[node].append((nb, d))
                stack.append(nb)

    leaf_order, order = [], []
    stack = [(root, False)]
    while stack:
        node, done = stack.pop()
        if done:
            order.append(node)
            continue
        stack.append((node, True))
        ch = children[node]
        if not ch:
            leaf_order.append(node)
        for c, _ in reversed(ch):
            stack.append((c, False))

    radius = {root: 0.0}
    stack = [root]
    while stack:
        node = stack.pop()
        for child, dist in children[node]:
            radius[child] = radius[node] + max(dist, 0.0)
            stack.append(child)
    return dict(children=children, order=order, leaf_order=leaf_order, radius=radius, root=root)


def _splits(tree, n):
    """{node: bipartition} for every branch; the side without sample 0 names the split."""
    leaves = {}
    for node in tree["order"]:
        ch = tree["children"][node]
        leaves[node] = frozenset([node]) if not ch else frozenset().union(*(leaves[c] for c, _ in ch))
    full = frozenset(range(n))
    out = {}
    for node, s in leaves.items():
        if node == tree["root"] or len(s) < 2 or len(s) > n - 2:
            continue
        out[node] = s if 0 not in s else full - s
    return out


def bootstrap_support(tree, X, n_boot, seed=1):
    """Share of bootstrap trees (blocks resampled with replacement) that contain each branch of `tree`."""
    n = X.shape[0]
    base = _splits(tree, n)
    if n_boot <= 0 or not base:
        return {}
    rng = np.random.default_rng(seed)
    wanted = {}
    for node, s in base.items():
        wanted.setdefault(s, []).append(node)
    hits = {s: 0 for s in wanted}
    for _ in range(n_boot):
        cols = rng.integers(0, X.shape[1], X.shape[1])
        rep = nj_structure(squareform(pdist(X[:, cols], metric='euclidean')))
        for s in set(_splits(rep, n).values()):
            if s in hits:
                hits[s] += 1
    return {node: hits[s] / n_boot for s, nodes in wanted.items() for node in nodes}


# =============================================================================
# CORE ANALYSIS FUNCTION
# =============================================================================

def run_diversity(input_file, n_clusters, top_blocks, outdir, radar_normalized=True):
    outdir.mkdir(parents=True, exist_ok=True)

    df_raw = pd.read_csv(input_file)

    # A second "Sample" column is read by pandas as "Sample.1"; it is a copy of
    # the sample IDs, not a haplotype block.
    if 'Sample' in df_raw.columns:
        dup = [c for c in df_raw.columns[2:] if re.fullmatch(r"Sample\.\d+", str(c))]
        df_raw = df_raw.drop(columns=dup)

    variety_col = df_raw.columns[0]
    sample_col = df_raw.columns[1] if len(df_raw.columns) > 1 else None
    block_cols = list(df_raw.columns[2:]) if len(df_raw.columns) > 2 else []

    if len(block_cols) == 0:
        raise ValueError(
            "CSV must have at least 3 columns: Variety, Sample_ID, and Haplotype Blocks"
        )

    df = df_raw.copy()
    # Blank / non-numeric block calls stay missing (NaN); -1 is also missing.
    # The metrics below need to see them as missing, not as zeros.
    df[block_cols] = df[block_cols].apply(pd.to_numeric, errors="coerce")
    n_missing_blocks = int((df[block_cols].isna() | (df[block_cols] == -1)).to_numpy().sum())
    pct_missing = 100 * n_missing_blocks / max(df[block_cols].size, 1)

    # ── 1. DIVERSITY METRICS ────────────────────────────────────────────────
    metrics = compute_diversity_metrics(df, block_cols)
    df = pd.concat([df, metrics], axis=1)


    # ── 2. BLOCK VARIANCE BAR ───────────────────────────────────────────────
    block_var = df[block_cols].var()   # missing (NaN) calls are skipped
    # PCA / clustering / tree need a complete matrix, as in the manual script.
    df[block_cols] = df[block_cols].fillna(0)
    block_var = block_var[block_var > 1e-4].sort_values(ascending=False)
    top_blocks = min(top_blocks, len(block_var))
    if top_blocks == 0:
        top_blocks = len(block_var)

    plot_paths = {}
    top_var = block_var.head(top_blocks)
    plot_paths["block_variance"] = None
    if len(top_var) > 0:
        plot_paths["block_variance"] = P.plot_block_variance(
            [str(b) for b in top_var.index], top_var.to_numpy(float),
            f"Top {top_blocks} Informative Haplotypic Blocks", outdir / "block_variance.png")

    # ── 3. DIVERSITY HEATMAP ────────────────────────────────────────────────
    metric_cols = [
        'Richness', 'Shannon_Diversity', 'Simpson_Diversity',
        'Gene_Diversity', 'PIC_Index', 'Heterozygosity'
    ]
    heatmap_data = df[[variety_col, sample_col] + metric_cols].copy()
    for col in metric_cols:
        max_val = heatmap_data[col].max()
        heatmap_data[col] = heatmap_data[col] / max_val if max_val > 0 else 0
    plot_paths["heatmap"] = P.plot_heatmap(
        heatmap_data[metric_cols].T.values, [m.replace("_", " ") for m in metric_cols],
        [str(v) for v in heatmap_data[variety_col]], [str(s) for s in heatmap_data[sample_col]],
        "Diversity Metrics Heatmap (Normalised)", outdir / "heatmap.png")

    # ── 4. VIOLIN PLOT ──────────────────────────────────────────────────────
    plot_paths["violin"] = None
    if df[variety_col].nunique() > 1:
        names = sorted(df[variety_col].unique(), key=str)
        plot_paths["violin"] = P.plot_violin(
            [df[df[variety_col] == v]['Shannon_Diversity'].to_numpy(float) for v in names], [str(v) for v in names],
            "Shannon Diversity Distribution by Variety", "Shannon Diversity Index", outdir / "violin.png")

    # ── 5. PCA ───────────────────────────────────────────────────────────────
    X = StandardScaler().fit_transform(df[block_cols].fillna(0))
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X)
    df["PC1"] = pcs[:, 0]
    df["PC2"] = pcs[:, 1]
    xl = f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)"
    yl = f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)"

    varieties = list(df[variety_col].unique())
    color_map = {v: P.color_for(i) for i, v in enumerate(varieties)}
    plot_paths["pca"] = P.plot_scatter(
        [(v, df.loc[df[variety_col] == v, "PC1"], df.loc[df[variety_col] == v, "PC2"], color_map[v])
         for v in sorted(varieties, key=str)],
        "PCA Diversity — coloured by variety", xl, yl, "Variety", outdir / "pca.png")

    # ── 6. HIERARCHICAL CLUSTERING ──────────────────────────────────────────
    Z = linkage(pdist(X), method="ward")
    df["Cluster"] = fcluster(Z, n_clusters, criterion="maxclust").astype(str)
    cluster_order = sorted(df["Cluster"].unique(), key=lambda x: int(x))
    plot_paths["clusters_pca"] = P.plot_scatter(
        [(f"Cluster {cl}", df.loc[df["Cluster"] == cl, "PC1"], df.loc[df["Cluster"] == cl, "PC2"], P.color_for(i))
         for i, cl in enumerate(cluster_order)],
        f"Hierarchical clusters ({n_clusters} clusters)", xl, yl, "Cluster", outdir / "clusters_pca.png")

    # ── 7. RADAR CHART ──────────────────────────────────────────────────────
    # radar_normalized=True: each metric scaled 0-1
    # radar_normalized=False: the actual mean value of each metric
    plot_paths["radar"] = None
    radar_metrics = [
        'Richness_Normalized', 'Shannon_Diversity', 'Simpson_Diversity',
        'Gene_Diversity', 'PIC_Index', 'Heterozygosity'
    ]
    if df[variety_col].nunique() > 1:
        varieties_sorted = sorted(df[variety_col].unique(), key=str)
        actual = {v: df[df[variety_col] == v][radar_metrics].mean().tolist() for v in varieties_sorted}
        if radar_normalized:
            maxima = {m: df[m].max() for m in radar_metrics}
            shown = {v: [(a / maxima[m] if maxima[m] > 0 else a) for m, a in zip(radar_metrics, vals)]
                     for v, vals in actual.items()}
            radial_max, radial_title = 1, "Normalised (0–1)"
        else:
            shown = actual
            radial_max = max(max(vals) for vals in actual.values()) * 1.1
            radial_title = "Actual mean values"
        plot_paths["radar"] = P.plot_radar(
            radar_metrics, [(v, shown[v], P.color_for(i)) for i, v in enumerate(varieties_sorted)],
            radial_max, radial_title, "Multi-metric diversity profile by variety", outdir / "radar.png")

    # ── 8. CIRCULAR NJ TREE ──────────────────────────────────────────────
    n_samples = len(df)
    dist_matrix = squareform(pdist(X, metric='euclidean'))
    sample_labels_list = [str(s) for s in df[sample_col].tolist()]
    variety_labels_list = [str(v) for v in df[variety_col].tolist()]
    plot_paths["nj_tree"] = None
    n_boot = 0
    if n_samples >= 3:
        try:
            tree = nj_structure(dist_matrix)
            n_boot = 100 if n_samples <= 150 else 50 if n_samples <= 400 else 0
            support = bootstrap_support(tree, X, n_boot) if n_boot else {}
            plot_paths["nj_tree"] = P.plot_nj_tree(
                tree, sample_labels_list, variety_labels_list, outdir / "nj_tree.png", support=support, n_boot=n_boot)
        except Exception as tree_error:                    # the other results are still returned
            plot_paths["nj_tree"] = None
            tree_note = str(tree_error)[:160]
            print(f"[DIVERSITY] tree failed: {tree_note}")

    # ── 9. SUMMARY (built directly as HTML — no markdown parser needed,
    #    so it renders with zero network dependency in QWebEngineView) ───────
    def esc(v):
        return str(v).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

    summary = f"""
<h2>Comprehensive Diversity Analysis Summary</h2>

<h3>Dataset Overview</h3>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td><strong>Total Varieties/Cultivars</strong></td><td>{df[variety_col].nunique()}</td></tr>
<tr><td><strong>Total Samples</strong></td><td>{len(df)}</td></tr>
<tr><td><strong>Total Haplotype Blocks</strong></td><td>{len(block_cols)}</td></tr>
<tr><td><strong>Informative Blocks</strong></td><td>{len(block_var)}</td></tr>
<tr><td><strong>Missing block calls (&minus;1 / blank)</strong></td><td>{n_missing_blocks:,} ({pct_missing:.2f}%)</td></tr>
</table>

<h3>Diversity Metrics (Mean &plusmn; Std)</h3>
<table>
<tr><th>Metric</th><th>Mean &plusmn; Std</th><th>Min &mdash; Max</th><th>Interpretation</th></tr>
<tr><td><strong>Richness</strong></td><td>{df['Richness'].mean():.2f} &plusmn; {df['Richness'].std():.2f}</td><td>{df['Richness'].min():.0f} &mdash; {df['Richness'].max():.0f}</td><td>Higher = more allelic variants</td></tr>
<tr><td><strong>Shannon Diversity</strong></td><td>{df['Shannon_Diversity'].mean():.3f} &plusmn; {df['Shannon_Diversity'].std():.3f}</td><td>{df['Shannon_Diversity'].min():.3f} &mdash; {df['Shannon_Diversity'].max():.3f}</td><td>0=no diversity, higher=greater evenness</td></tr>
<tr><td><strong>Simpson Diversity</strong></td><td>{df['Simpson_Diversity'].mean():.3f} &plusmn; {df['Simpson_Diversity'].std():.3f}</td><td>{df['Simpson_Diversity'].min():.3f} &mdash; {df['Simpson_Diversity'].max():.3f}</td><td>Probability 2 random alleles differ</td></tr>
<tr><td><strong>Gene Diversity</strong></td><td>{df['Gene_Diversity'].mean():.3f} &plusmn; {df['Gene_Diversity'].std():.3f}</td><td>{df['Gene_Diversity'].min():.3f} &mdash; {df['Gene_Diversity'].max():.3f}</td><td>Expected heterozygosity</td></tr>
<tr><td><strong>PIC Index</strong></td><td>{df['PIC_Index'].mean():.3f} &plusmn; {df['PIC_Index'].std():.3f}</td><td>{df['PIC_Index'].min():.3f} &mdash; {df['PIC_Index'].max():.3f}</td><td>&gt;0.5 = highly informative</td></tr>
<tr><td><strong>Effective Alleles</strong></td><td>{df['Effective_Alleles'].mean():.2f} &plusmn; {df['Effective_Alleles'].std():.2f}</td><td>{df['Effective_Alleles'].min():.1f} &mdash; {df['Effective_Alleles'].max():.1f}</td><td>Equally-frequent alleles equiv.</td></tr>
<tr><td><strong>Gini-Simpson</strong></td><td>{df['Gini_Simpson'].mean():.3f} &plusmn; {df['Gini_Simpson'].std():.3f}</td><td>{df['Gini_Simpson'].min():.3f} &mdash; {df['Gini_Simpson'].max():.3f}</td><td>0=low, 1=high diversity</td></tr>
<tr><td><strong>Heterozygosity</strong></td><td>{df['Heterozygosity'].mean():.3f} &plusmn; {df['Heterozygosity'].std():.3f}</td><td>{df['Heterozygosity'].min():.3f} &mdash; {df['Heterozygosity'].max():.3f}</td><td>Within-sample genetic variation</td></tr>
<tr><td><strong>Nucleotide Diversity</strong></td><td>{df['Nucleotide_Diversity'].mean():.3f} &plusmn; {df['Nucleotide_Diversity'].std():.3f}</td><td>{df['Nucleotide_Diversity'].min():.3f} &mdash; {df['Nucleotide_Diversity'].max():.3f}</td><td>Pairwise sequence divergence</td></tr>
<tr><td><strong>Allelic Entropy</strong></td><td>{df['Allelic_Entropy'].mean():.3f} &plusmn; {df['Allelic_Entropy'].std():.3f}</td><td>{df['Allelic_Entropy'].min():.3f} &mdash; {df['Allelic_Entropy'].max():.3f}</td><td>Information entropy of allele freq.</td></tr>
<tr><td><strong>Load</strong></td><td>{df['Load'].mean():.2f} &plusmn; {df['Load'].std():.2f}</td><td>{df['Load'].min():.0f} &mdash; {df['Load'].max():.0f}</td><td>Cumulative haplotype burden</td></tr>
</table>

<h3>Variety-Level Diversity</h3>
"""
    for variety in sorted(df[variety_col].unique()):
        sub = df[df[variety_col] == variety]
        summary += (
            f"<p><strong>{esc(variety)}</strong> ({len(sub)} samples)<br>"
            f"Avg Richness: {sub['Richness'].mean():.2f} &nbsp;|&nbsp; "
            f"Avg Shannon: {sub['Shannon_Diversity'].mean():.3f} &nbsp;|&nbsp; "
            f"Avg PIC: {sub['PIC_Index'].mean():.3f}</p>\n"
        )

    summary += f"""
<h3>Clustering Analysis ({n_clusters} clusters)</h3>
<table>
<tr><th>Cluster</th><th>Samples</th><th>Avg Richness</th><th>Avg Shannon</th><th>Primary Varieties</th></tr>
"""
    for cl in cluster_order:
        sub = df[df["Cluster"] == cl]
        primary_vars = sub[variety_col].value_counts().head(2).index.tolist()
        primary_vars_str = ", ".join(esc(v) for v in primary_vars) if primary_vars else "Mixed"
        summary += (
            f"<tr><td>Cluster {esc(cl)}</td><td>{len(sub)}</td>"
            f"<td>{sub['Richness'].mean():.2f}</td>"
            f"<td>{sub['Shannon_Diversity'].mean():.3f}</td>"
            f"<td>{primary_vars_str}</td></tr>\n"
        )
    summary += "</table>\n"

    summary += """
<h3>Interpretation Guide</h3>
<p><strong>Shannon Diversity</strong>: Measures both richness and evenness.</p>
<ul>
<li>Low (&lt;1.0): Limited allelic variation</li>
<li>Moderate (1.0&ndash;2.0): Good diversity</li>
<li>High (&gt;2.0): Excellent genetic diversity</li>
</ul>
<p><strong>PIC Index</strong>: Polymorphic Information Content.</p>
<ul>
<li>&lt;0.25: Low information</li>
<li>0.25&ndash;0.5: Intermediate</li>
<li>&gt;0.5: High information (highly polymorphic)</li>
</ul>
<p><strong>Circular NJ Tree</strong>: genetic relatedness of the samples.</p>
<ul>
<li>The root is at the centre; distance from the centre is genetic distance along the branches (scale bar at the bottom left)</li>
<li>Pale sectors and dots mark the variety/cultivar of every sample</li>
<li>Red triangles: branches found in at least 70% of bootstrap trees (blocks resampled with replacement)</li>
</ul>
"""

    # ── 10. RESULTS TABLE ───────────────────────────────────────────────────
    results_df = df[[
        variety_col, sample_col, 'Richness', 'Shannon_Diversity',
        'Simpson_Diversity', 'Gene_Diversity', 'PIC_Index',
        'Heterozygosity', 'Nucleotide_Diversity', 'Allelic_Entropy',
        'Effective_Alleles', 'Load', 'Cluster'
    ]].copy().round(4)

    results_csv_path = outdir / "results.csv"
    results_df.to_csv(results_csv_path, index=False)

    results_json_path = outdir / "results.json"
    results_df.to_json(results_json_path, orient="records")

    return {
        "status": "success",
        "n_varieties": int(df[variety_col].nunique()),
        "n_samples": int(len(df)),
        "n_blocks": int(len(block_cols)),
        "n_informative_blocks": int(len(block_var)),
        "n_clusters": int(n_clusters),
        "radar_normalized": bool(radar_normalized),
        "pct_missing": round(pct_missing, 3),
        "summary_html": summary,
        "results_csv": str(results_csv_path),
        "results_json": str(results_json_path),
        "plots": plot_paths,
    }
