"""
diversity_engine.py — Haplotypic diversity analysis engine.

This is the same computation logic as the Hugging Face app.py (Gradio)
version: compute_diversity_metrics(), the custom neighbour-joining
implementation, build_nj_tree_plotly(), and the analyze_haplotypes()
pipeline. No math has been changed.

The only thing that's different is the "shape" of the output: instead of
returning Plotly Figure objects + a DataFrame for Gradio to render, this
writes each figure as a self-contained interactive HTML file (Plotly.js
embedded, so it works fully offline in QWebEngineView) plus a results CSV
and JSON, and returns a dict of paths + the summary markdown string.
"""

import warnings
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings('ignore')

# =============================================================================
# COLOR PALETTES (unchanged from app.py)
# =============================================================================
DARK_COLORS = [
    "#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD",
    "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
]

VARIETY_CLADE_COLORS = [
    "#C0392B", "#E67E22", "#F1C40F", "#27AE60", "#1ABC9C",
    "#2980B9", "#8E44AD", "#E91E63", "#00BCD4", "#FF5722",
    "#9C27B0", "#4CAF50",
]


def hex_to_rgba(hex_color, alpha=0.25):
    hex_color = hex_color.lstrip('#')
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    return f"rgba({r}, {g}, {b}, {alpha})"


# =============================================================================
# DIVERSITY METRICS COMPUTATION (10+ INDICES) — unchanged
# =============================================================================

def compute_diversity_metrics(df, block_cols):
    metrics = pd.DataFrame(index=df.index)
    data = df[block_cols].fillna(0)

    metrics['Richness'] = (data != 0).sum(axis=1)
    metrics['Load'] = data.sum(axis=1)

    data_norm = data.div(data.sum(axis=1), axis=0).replace([np.inf, -np.inf], 0).fillna(0)
    shannon = -(data_norm * np.log(data_norm + 1e-10)).sum(axis=1)
    metrics['Shannon_Diversity'] = shannon

    simpson = 1 - (data_norm ** 2).sum(axis=1)
    metrics['Simpson_Diversity'] = simpson
    metrics['Gini_Simpson'] = (data_norm * (1 - data_norm)).sum(axis=1)
    metrics['Effective_Alleles'] = 1 / (data_norm ** 2).sum(axis=1).replace(0, np.nan)
    metrics['Effective_Alleles'] = metrics['Effective_Alleles'].fillna(0)

    nucleotide_div = []
    for idx, row in data.iterrows():
        if row.sum() > 0:
            diff = np.abs(np.diff(row.values))
            nuc_div = np.mean(diff) if len(diff) > 0 else 0
        else:
            nuc_div = 0
        nucleotide_div.append(nuc_div)
    metrics['Nucleotide_Diversity'] = nucleotide_div

    het = []
    for idx, row in data.iterrows():
        unique_vals = len(row[row > 0].unique())
        max_het = unique_vals / len(block_cols) if len(block_cols) > 0 else 0
        het.append(max_het)
    metrics['Heterozygosity'] = het

    entropy = []
    for idx, row in data.iterrows():
        non_zero = row[row > 0]
        if len(non_zero) > 0:
            probs = non_zero.value_counts(normalize=True).values
            ent = -np.sum(probs * np.log2(probs + 1e-10))
        else:
            ent = 0
        entropy.append(ent)
    metrics['Allelic_Entropy'] = entropy

    pic = []
    for idx, row in data.iterrows():
        if row.sum() > 0:
            counts = row[row > 0].value_counts()
            freq = counts / counts.sum()
            pic_val = 1 - (freq ** 2).sum()
        else:
            pic_val = 0
        pic.append(pic_val)
    metrics['PIC_Index'] = pic

    metrics['Gene_Diversity'] = 1 - (data_norm ** 2).sum(axis=1)
    max_richness = len(block_cols)
    metrics['Richness_Normalized'] = metrics['Richness'] / max_richness if max_richness > 0 else 0

    return metrics


# =============================================================================
# PLOT STYLING — unchanged
# =============================================================================

def apply_dark_style(fig, title="", xaxis_title="", yaxis_title="",
                      height=480, margin=None, extra=None):
    if margin is None:
        margin = dict(l=60, r=40, t=65, b=60)

    layout_kw = dict(
        font=dict(family="Consolas, 'Courier New', monospace", color="#333333"),
        paper_bgcolor="#FFFFFF",
        plot_bgcolor="#F8F9FA",
        title=dict(
            text=title,
            font=dict(size=15, color="#2C3E50", family="Consolas, monospace"),
            x=0.03,
        ),
        legend=dict(
            bgcolor="rgba(255,255,255,0.9)",
            bordercolor="#CCCCCC",
            borderwidth=1,
            font=dict(color="#333333"),
        ),
        xaxis=dict(
            title=xaxis_title,
            gridcolor="#E0E0E0",
            linecolor="#CCCCCC",
            tickfont=dict(size=10, color="#555555"),
            title_font=dict(color="#333333"),
        ),
        yaxis=dict(
            title=yaxis_title,
            gridcolor="#E0E0E0",
            linecolor="#CCCCCC",
            tickfont=dict(size=10, color="#555555"),
            title_font=dict(color="#333333"),
        ),
        coloraxis=dict(colorbar=dict(tickfont=dict(size=10, color="#555555"),
                                      title_font=dict(color="#333333"))),
        height=height,
        margin=margin,
    )

    if extra:
        layout_kw.update(extra)

    fig.update_layout(**layout_kw)
    return fig


# =============================================================================
# NEIGHBOUR JOINING ALGORITHM — unchanged
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
        Q = np.zeros((m, m))
        for i in range(m):
            for j in range(m):
                if i != j:
                    Q[i, j] = (m - 2) * D[i, j] - row_sums[i] - row_sums[j]

        np.fill_diagonal(Q, np.inf)
        flat_idx = np.argmin(Q)
        i, j = divmod(flat_idx, m)

        li = 0.5 * (D[i, j] + (row_sums[i] - row_sums[j]) / (m - 2))
        lj = D[i, j] - li

        new_node = next_id
        next_id += 1
        edges.append((new_node, clusters[i], max(li, 0)))
        edges.append((new_node, clusters[j], max(lj, 0)))

        old_indices = [k for k in range(m) if k != i and k != j]
        new_dists = [(D[i, k] + D[j, k] - D[i, j]) / 2 for k in old_indices]

        new_clusters = [clusters[k] for k in old_indices] + [new_node]
        new_size = len(new_clusters)
        new_D = np.zeros((new_size, new_size))
        for a, ka in enumerate(old_indices):
            for b, kb in enumerate(old_indices):
                new_D[a, b] = D[ka, kb]
        for a in range(len(old_indices)):
            new_D[a, new_size - 1] = new_dists[a]
            new_D[new_size - 1, a] = new_dists[a]

        D = new_D
        clusters = new_clusters

    if len(clusters) == 2:
        edges.append((clusters[0], clusters[1], D[0, 1] / 2))
        edges.append((clusters[1], clusters[0], D[0, 1] / 2))

    return edges


# =============================================================================
# RECTANGULAR CLADOGRAM NJ TREE — unchanged
# =============================================================================

def build_nj_tree_plotly(dist_matrix, sample_labels, variety_labels):
    """
    Rectangular cladogram NJ tree:
      - Root on LEFT, leaves extend to RIGHT
      - Each variety's entire subtree rendered in a distinct solid color
      - Grey backbone connects variety clades at deep divergence
      - Colored clade label bars on the far right margin
      - Branch lines are thin horizontal + vertical L-shapes (no diagonals)
      - Leaf labels shown for <=120 samples
    """
    n = len(sample_labels)
    if n < 3:
        fig = go.Figure()
        apply_dark_style(fig, title="NJ Tree (need >= 3 samples)")
        return fig

    edges = neighbour_joining_simple(dist_matrix)

    adj = defaultdict(list)
    for parent, child, dist in edges:
        adj[parent].append((child, dist))
        adj[child].append((parent, dist))

    all_nodes = set()
    for p, c, _ in edges:
        all_nodes.add(p)
        all_nodes.add(c)
    root = max(all_nodes)

    children_map = defaultdict(list)
    parent_map = {}
    visited = set()

    def dfs_root(node, par):
        visited.add(node)
        for nb, d in adj[node]:
            if nb not in visited:
                children_map[node].append((nb, d))
                parent_map[nb] = node
                dfs_root(nb, node)

    dfs_root(root, -1)

    leaf_order = []

    def collect_leaves(node):
        ch = children_map[node]
        if not ch:
            leaf_order.append(node)
        else:
            for c, _ in ch:
                collect_leaves(c)

    collect_leaves(root)

    leaf_y = {leaf: i for i, leaf in enumerate(leaf_order)}

    node_x = {root: 0.0}
    stack = [root]
    while stack:
        node = stack.pop()
        for child, dist in children_map[node]:
            node_x[child] = node_x[node] + dist
            stack.append(child)

    node_y = {}

    def assign_y(node):
        ch = children_map[node]
        if not ch:
            node_y[node] = leaf_y[node]
        else:
            for c, _ in ch:
                assign_y(c)
            ys = [node_y[c] for c, _ in ch]
            node_y[node] = (min(ys) + max(ys)) / 2.0

    assign_y(root)

    max_x = max(node_x.values()) if node_x else 1.0
    if max_x == 0:
        max_x = 1.0
    node_x = {k: v / max_x for k, v in node_x.items()}

    unique_varieties = []
    seen = set()
    for v in variety_labels:
        if v not in seen:
            unique_varieties.append(v)
            seen.add(v)
    var_color = {v: VARIETY_CLADE_COLORS[i % len(VARIETY_CLADE_COLORS)]
                 for i, v in enumerate(unique_varieties)}

    def get_leaf_varieties(node):
        ch = children_map[node]
        if not ch:
            return {variety_labels[node]}
        result = set()
        for c, _ in ch:
            result |= get_leaf_varieties(c)
        return result

    node_variety = {}

    def assign_node_variety(node):
        ch = children_map[node]
        if not ch:
            node_variety[node] = variety_labels[node]
            return
        for c, _ in ch:
            assign_node_variety(c)
        child_vars = [node_variety[c] for c, _ in ch]
        most_common = Counter(child_vars).most_common(1)[0][0]
        node_variety[node] = most_common

    assign_node_variety(root)

    def is_pure_variety(node):
        vs = get_leaf_varieties(node)
        return list(vs)[0] if len(vs) == 1 else None

    purity_cache = {}

    def get_purity(node):
        if node not in purity_cache:
            purity_cache[node] = is_pure_variety(node)
        return purity_cache[node]

    BACKBONE_COLOR = "#666666"

    def edge_color(child_node):
        p = get_purity(child_node)
        if p is not None:
            return var_color[p]
        return BACKBONE_COLOR

    color_traces = defaultdict(lambda: {"x": [], "y": []})

    def add_edge(px, py, cx, cy, color):
        color_traces[color]["x"] += [px, cx, None]
        color_traces[color]["y"] += [cy, cy, None]
        if abs(py - cy) > 0.001:
            color_traces[color]["x"] += [px, px, None]
            color_traces[color]["y"] += [py, cy, None]

    def draw_edges(node):
        ch = children_map[node]
        px = node_x[node]
        py = node_y[node]
        for child, _ in ch:
            cx = node_x[child]
            cy = node_y[child]
            color = edge_color(child)
            add_edge(px, py, cx, cy, color)
            draw_edges(child)

    draw_edges(root)

    fig = go.Figure()

    for color, coords in color_traces.items():
        is_backbone = (color == BACKBONE_COLOR)
        fig.add_trace(go.Scatter(
            x=coords["x"],
            y=coords["y"],
            mode="lines",
            line=dict(color=color, width=0.9 if not is_backbone else 1.1),
            hoverinfo="none",
            showlegend=False,
        ))

    for variety in unique_varieties:
        idx_list = [i for i in range(n) if variety_labels[i] == variety]
        if not idx_list:
            continue
        lx = [node_x[i] for i in idx_list if i in node_x]
        ly = [node_y[i] for i in idx_list if i in node_y]
        hover = [f"<b>{sample_labels[i]}</b><br>Variety: {variety}" for i in idx_list if i in node_x]
        fig.add_trace(go.Scatter(
            x=lx, y=ly,
            mode="markers",
            name=variety,
            marker=dict(
                size=4,
                color=var_color[variety],
                line=dict(color="rgba(0,0,0,0.3)", width=0.3),
            ),
            text=hover,
            hovertemplate="%{text}<extra></extra>",
        ))

    show_labels = n <= 120
    if show_labels:
        for i in range(n):
            if i not in node_x:
                continue
            fig.add_annotation(
                x=node_x[i] + 0.012,
                y=node_y[i],
                text=sample_labels[i],
                showarrow=False,
                xanchor="left",
                font=dict(
                    size=7,
                    color=var_color.get(variety_labels[i], "#AAAAAA"),
                    family="Consolas, monospace",
                ),
                yshift=0,
            )

    bar_x_right = 1.08 if show_labels else 1.04
    label_x = bar_x_right + 0.025

    i = 0
    while i < len(leaf_order):
        leaf = leaf_order[i]
        var = variety_labels[leaf]
        j = i
        while j < len(leaf_order) and variety_labels[leaf_order[j]] == var:
            j += 1
        y_min = i - 0.4
        y_max = j - 1 + 0.4
        y_mid = (i + j - 1) / 2.0
        color = var_color[var]
        fig.add_shape(
            type="rect",
            x0=bar_x_right, x1=bar_x_right + 0.018,
            y0=y_min, y1=y_max,
            fillcolor=color,
            line=dict(width=0),
            xref="x", yref="y",
        )
        fig.add_annotation(
            x=label_x + 0.018,
            y=y_mid,
            text=f"<b>{var}</b>",
            showarrow=False,
            xanchor="left",
            font=dict(size=9, color=color, family="Consolas, monospace"),
        )
        i = j

    fig.add_trace(go.Scatter(
        x=[node_x[root]],
        y=[node_y[root]],
        mode="markers",
        marker=dict(size=6, color=BACKBONE_COLOR, symbol="diamond"),
        hoverinfo="text",
        text=["Root"],
        showlegend=False,
    ))

    plot_height = max(600, n * 10)
    x_range_right = label_x + 0.18
    if show_labels:
        x_range_right += 0.06

    apply_dark_style(
        fig,
        title="Neighbour Joining Phylogenetic Tree (Rectangular Cladogram)",
        height=plot_height,
        margin=dict(l=40, r=20, t=65, b=40),
        extra=dict(
            xaxis=dict(
                visible=True,
                title="Genetic Distance (normalized)",
                range=[-0.05, x_range_right],
                gridcolor="#2A2A4A",
                tickfont=dict(size=9, color="#A0A0C0"),
                title_font=dict(color="#C0C0E0"),
                zeroline=False,
            ),
            yaxis=dict(
                visible=False,
                range=[-1, n],
            ),
            legend=dict(
                title=dict(text="Variety / Cultivar", font=dict(color="#C0C0E0")),
                yanchor="top", y=0.99, xanchor="left", x=0.01,
                itemsizing="constant",
            ),
        ),
    )

    return fig


# =============================================================================
# CORE ANALYSIS FUNCTION
# =============================================================================
# Same as app.py's analyze_haplotypes(), but writes outputs to disk instead
# of returning Gradio component values.

def run_diversity(input_file, n_clusters, top_blocks, outdir):
    outdir.mkdir(parents=True, exist_ok=True)

    df_raw = pd.read_csv(input_file)

    variety_col = df_raw.columns[0]
    sample_col = df_raw.columns[1] if len(df_raw.columns) > 1 else None
    block_cols = list(df_raw.columns[2:]) if len(df_raw.columns) > 2 else []

    if len(block_cols) == 0:
        raise ValueError(
            "CSV must have at least 3 columns: Variety, Sample_ID, and Haplotype Blocks"
        )

    df = df_raw.copy()
    df[block_cols] = df[block_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

    # ── 1. DIVERSITY METRICS ────────────────────────────────────────────────
    metrics = compute_diversity_metrics(df, block_cols)
    df = pd.concat([df, metrics], axis=1)

    # ── 2. BLOCK VARIANCE BAR ───────────────────────────────────────────────
    block_var = df[block_cols].var()
    block_var = block_var[block_var > 1e-4].sort_values(ascending=False)
    top_blocks = min(top_blocks, len(block_var))
    if top_blocks == 0:
        top_blocks = len(block_var)

    top_var = block_var.head(top_blocks).reset_index()
    fig_block = None
    if len(top_var) > 0:
        top_var.columns = ["Block", "Variance"]
        top_var["Label"] = top_var["Block"].astype(str)
        fig_block = go.Figure(go.Bar(
            x=top_var["Label"],
            y=top_var["Variance"],
            marker=dict(
                color=top_var["Variance"],
                colorscale=[[0.0, "#1A3A5C"], [0.5, "#1F77B4"], [1.0, "#00BFFF"]],
                line=dict(color="#000000", width=0.6),
                showscale=True,
                colorbar=dict(
                    title=dict(text="Variance", font=dict(color="#C0C0E0")),
                    thickness=12,
                    tickfont=dict(color="#C0C0E0"),
                ),
            ),
            hovertemplate="<b>%{x}</b><br>Variance: %{y:.4f}<extra></extra>",
        ))
        apply_dark_style(
            fig_block,
            title=f"Top {top_blocks} Informative Haplotypic Blocks",
            xaxis_title="Haplotype Block",
            yaxis_title="Variance (log scale)",
            margin=dict(l=60, r=40, t=65, b=110),
            extra=dict(yaxis_type="log", xaxis_tickangle=-45),
        )

    # ── 3. DIVERSITY HEATMAP ────────────────────────────────────────────────
    metric_cols = [
        'Richness', 'Shannon_Diversity', 'Simpson_Diversity',
        'Gene_Diversity', 'PIC_Index', 'Heterozygosity'
    ]
    heatmap_data = df[[variety_col, sample_col] + metric_cols].copy()
    for col in metric_cols:
        max_val = heatmap_data[col].max()
        heatmap_data[col] = heatmap_data[col] / max_val if max_val > 0 else 0

    heatmap_matrix = heatmap_data[metric_cols].T.values
    hover_text = []
    for metric in metric_cols:
        row_hover = []
        for idx, sample in enumerate(heatmap_data[sample_col]):
            variety = heatmap_data.iloc[idx][variety_col]
            value = df[metric].iloc[idx]
            row_hover.append(
                f"<b>{metric}</b><br>Variety: {variety}<br>Sample: {sample}<br>Value: {value:.4f}")
        hover_text.append(row_hover)

    fig_heatmap = go.Figure(data=go.Heatmap(
        z=heatmap_matrix,
        x=[f"{v}<br>{s}" for v, s in zip(heatmap_data[variety_col], heatmap_data[sample_col])],
        y=metric_cols,
        colorscale=[[0.0, "#0D1B2A"], [0.4, "#1F4E79"], [0.7, "#2874A6"], [1.0, "#1ABC9C"]],
        hovertemplate="%{customdata}<extra></extra>",
        customdata=hover_text,
    ))
    apply_dark_style(
        fig_heatmap,
        title="Diversity Metrics Heatmap (Normalized)",
        xaxis_title="Variety — Sample",
        yaxis_title="Diversity Metric",
        margin=dict(l=120, r=40, t=65, b=120),
        extra=dict(xaxis_tickangle=-45),
    )

    # ── 4. VIOLIN PLOT ──────────────────────────────────────────────────────
    fig_violin = None
    if df[variety_col].nunique() > 1:
        fig_violin = go.Figure()
        for i, variety in enumerate(sorted(df[variety_col].unique())):
            sub = df[df[variety_col] == variety]
            fig_violin.add_trace(go.Violin(
                y=sub['Shannon_Diversity'],
                name=variety,
                box_visible=True,
                meanline_visible=True,
                opacity=0.85,
                points=False,
                line_color=DARK_COLORS[i % len(DARK_COLORS)],
                fillcolor=DARK_COLORS[i % len(DARK_COLORS)],
            ))
        apply_dark_style(
            fig_violin,
            title="Shannon Diversity Distribution by Variety",
            xaxis_title="Variety",
            yaxis_title="Shannon Diversity Index",
            height=500,
            margin=dict(l=60, r=40, t=65, b=80),
        )

    # ── 5. PCA ───────────────────────────────────────────────────────────────
    X = StandardScaler().fit_transform(df[block_cols].fillna(0))
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X)
    df["PC1"] = pcs[:, 0]
    df["PC2"] = pcs[:, 1]

    varieties = df[variety_col].unique()
    color_map = {v: DARK_COLORS[i % len(DARK_COLORS)] for i, v in enumerate(varieties)}

    fig_pca = go.Figure()
    for variety in sorted(varieties):
        sub = df[df[variety_col] == variety]
        fig_pca.add_trace(go.Scatter(
            x=sub["PC1"], y=sub["PC2"],
            mode="markers",
            name=variety,
            marker=dict(
                size=10,
                color=color_map[variety],
                line=dict(color="#000000", width=0.6),
                opacity=0.9,
            ),
            text=[f"<b>{v}</b><br>Sample: {s}<br>Shannon: {sh:.3f}"
                  for v, s, sh in zip(sub[variety_col], sub[sample_col], sub['Shannon_Diversity'])],
            hovertemplate="%{text}<br>PC1: %{x:.3f}<br>PC2: %{y:.3f}<extra></extra>",
        ))
    apply_dark_style(
        fig_pca,
        title="PCA Diversity — Colored by Variety",
        xaxis_title=f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)",
        yaxis_title=f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)",
        extra=dict(legend=dict(title="Variety", yanchor="top", y=0.99, xanchor="left", x=0.01)),
    )

    # ── 6. HIERARCHICAL CLUSTERING ──────────────────────────────────────────
    Z = linkage(pdist(X), method="ward")
    df["Cluster"] = fcluster(Z, n_clusters, criterion="maxclust").astype(str)
    cluster_order = sorted(df["Cluster"].unique(), key=lambda x: int(x))

    fig_cpca = go.Figure()
    for i, cl in enumerate(cluster_order):
        sub = df[df["Cluster"] == cl]
        fig_cpca.add_trace(go.Scatter(
            x=sub["PC1"], y=sub["PC2"],
            mode="markers",
            name=f"Cluster {cl}",
            marker=dict(
                size=10,
                color=DARK_COLORS[i % len(DARK_COLORS)],
                line=dict(color="#000000", width=0.6),
                opacity=0.9,
                symbol="circle",
            ),
            text=[f"<b>Cluster {cl}</b><br>{v} — {s}<br>Richness: {r:.0f}"
                  for v, s, r in zip(sub[variety_col], sub[sample_col], sub['Richness'])],
            hovertemplate="%{text}<br>PC1: %{x:.3f}<br>PC2: %{y:.3f}<extra></extra>",
        ))
    apply_dark_style(
        fig_cpca,
        title=f"Hierarchical Clusters ({n_clusters} Clusters) with Varieties",
        xaxis_title=f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)",
        yaxis_title=f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)",
        extra=dict(legend=dict(title="Cluster", yanchor="top", y=0.99, xanchor="right", x=0.99)),
    )

    # ── 7. RADAR CHART ──────────────────────────────────────────────────────
    fig_radar = None
    if df[variety_col].nunique() <= 8:
        fig_radar = go.Figure()
        radar_metrics = [
            'Richness_Normalized', 'Shannon_Diversity', 'Simpson_Diversity',
            'Gene_Diversity', 'PIC_Index', 'Heterozygosity'
        ]
        norm_df = df.copy()
        for metric in radar_metrics:
            max_val = norm_df[metric].max()
            if max_val > 0:
                norm_df[metric] = norm_df[metric] / max_val

        for i, variety in enumerate(sorted(df[variety_col].unique())):
            sub = norm_df[norm_df[variety_col] == variety]
            avg_vals = sub[radar_metrics].mean().tolist()
            fill_color = hex_to_rgba(DARK_COLORS[i % len(DARK_COLORS)], 0.25)
            fig_radar.add_trace(go.Scatterpolar(
                r=avg_vals,
                theta=radar_metrics,
                fill='toself',
                name=variety,
                line=dict(color=DARK_COLORS[i % len(DARK_COLORS)], width=2),
                fillcolor=fill_color,
                opacity=0.85,
            ))

        fig_radar.update_layout(
            polar=dict(
                radialaxis=dict(
                    visible=True, range=[0, 1],
                    gridcolor="#2A2A4A",
                    tickfont=dict(color="#A0A0C0", size=9),
                ),
                angularaxis=dict(
                    gridcolor="#2A2A4A",
                    tickfont=dict(color="#C0C0E0", size=10),
                ),
                bgcolor="#16213E",
            ),
            paper_bgcolor="#1A1A2E",
            font=dict(family="Consolas, monospace", color="#E0E0E0"),
            title=dict(
                text="Multi-Metric Diversity Profile by Variety",
                font=dict(size=15, color="#A0A8D0"),
            ),
            legend=dict(
                bgcolor="rgba(20,20,40,0.85)",
                bordercolor="#444466",
                borderwidth=1,
                font=dict(color="#C0C0E0"),
            ),
            height=580,
            margin=dict(l=80, r=80, t=80, b=80),
        )

    # ── 8. RECTANGULAR NJ TREE ──────────────────────────────────────────────
    dist_matrix = squareform(pdist(X, metric='euclidean'))
    sample_labels_list = df[sample_col].tolist()
    variety_labels_list = df[variety_col].tolist()

    try:
        fig_nj = build_nj_tree_plotly(dist_matrix, sample_labels_list, variety_labels_list)
    except Exception as tree_error:
        fig_nj = go.Figure()
        fig_nj.add_annotation(
            text=f"Tree error: {str(tree_error)[:120]}",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(color="#FF6666")
        )
        apply_dark_style(fig_nj, title="NJ Tree (Error)")

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
<p><strong>Rectangular NJ Tree</strong>: Cladogram showing genetic relatedness.</p>
<ul>
<li>Root on LEFT, leaves extend RIGHT</li>
<li>Branch colors match variety/cultivar grouping</li>
<li>Grey backbone = mixed-variety nodes (deep divergence)</li>
<li>Colored bars on right margin = variety clade labels</li>
<li>Branch length proportional to genetic distance</li>
</ul>
"""

    # ── 10. RESULTS TABLE ───────────────────────────────────────────────────
    results_df = df[[
        variety_col, sample_col, 'Richness', 'Shannon_Diversity',
        'Simpson_Diversity', 'Gene_Diversity', 'PIC_Index',
        'Heterozygosity', 'Nucleotide_Diversity', 'Allelic_Entropy',
        'Effective_Alleles', 'Load', 'Cluster'
    ]].copy().round(4)

    # ── WRITE OUTPUTS ────────────────────────────────────────────────────────
    figures = {
        "block_variance": fig_block,
        "pca": fig_pca,
        "clusters_pca": fig_cpca,
        "heatmap": fig_heatmap,
        "violin": fig_violin,
        "radar": fig_radar,
        "nj_tree": fig_nj,
    }

    plot_paths = {}
    for name, fig in figures.items():
        if fig is None:
            plot_paths[name] = None
            continue
        path = outdir / f"{name}.html"
        fig.write_html(str(path), full_html=True, include_plotlyjs=True)
        plot_paths[name] = str(path)

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
        "summary_html": summary,
        "results_csv": str(results_csv_path),
        "results_json": str(results_json_path),
        "plots": plot_paths,
    }
