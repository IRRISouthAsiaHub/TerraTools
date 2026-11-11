"""
diversity_plots.py — figures of the Diversity Analysis tool.

Every figure is drawn with matplotlib on a white background and saved as a PNG,
so the app can show it at its natural proportions (never squeezed), open it in
a pop-up and save it as an image.

The circular neighbour-joining tree follows the look of an iTOL export: black
branches, a ring of pale coloured sectors (one colour per variety) behind the
sample names, light concentric distance circles and red triangles on the nodes
with bootstrap support.
"""

from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
import numpy as np                                    # noqa: E402
from matplotlib import cm                             # noqa: E402
from matplotlib.collections import LineCollection     # noqa: E402
from matplotlib.colors import LogNorm, Normalize, to_rgb   # noqa: E402
from matplotlib.lines import Line2D                   # noqa: E402
from matplotlib.patches import Circle, Wedge          # noqa: E402

DPI = 200

# Distinct, print-friendly colours (varieties, clusters)
PALETTE = [
    "#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD",
    "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
    "#393B79", "#AD494A", "#637939", "#CE6DBD", "#8C6D31",
    "#3182BD", "#E6550D", "#31A354", "#756BB1", "#636363",
]

INK = "#2C3E50"
GRID = "#E3E6EA"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.edgecolor": "#B8BEC6",
    "axes.labelcolor": "#333333",
    "axes.titlecolor": INK,
    "xtick.color": "#555555",
    "ytick.color": "#555555",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
})


def color_for(i):
    return PALETTE[i % len(PALETTE)]


def _save(fig, path):
    fig.savefig(str(path), dpi=DPI, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)
    return str(path)


def _style_axes(ax, title, xlabel="", ylabel=""):
    ax.set_title(title, loc="left", fontsize=14, fontweight="bold", color=INK, pad=14)
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def _legend_outside(ax, title, ncol=1):
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return
    ax.legend(handles, labels, title=title, loc="upper left", bbox_to_anchor=(1.02, 1.0),
              frameon=False, ncol=ncol, fontsize=10, title_fontsize=11, borderaxespad=0)


# =============================================================================
# Block variance bar chart
# =============================================================================

def plot_block_variance(labels, variances, title, path):
    n = len(labels)
    fig, ax = plt.subplots(figsize=(max(9.5, 0.36 * n + 3.0), 6.4))
    norm = LogNorm(vmin=max(min(variances), 1e-6), vmax=max(variances)) if max(variances) > min(variances) else None
    cmap = plt.get_cmap("YlGnBu")
    colors = [cmap(0.35 + 0.6 * (norm(v) if norm else 0.5)) for v in variances]
    ax.bar(range(n), variances, color=colors, edgecolor="#4A5568", linewidth=0.5)
    ax.set_yscale("log")
    ax.set_xticks(range(n))
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=9)
    ax.set_xlim(-0.8, n - 0.2)
    _style_axes(ax, title, "Haplotype block", "Variance (log scale)")
    ax.grid(axis="x", visible=False)
    if norm:
        sm = cm.ScalarMappable(norm=Normalize(0.35, 0.95), cmap=cmap)
        sm.set_array([])
        cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=matplotlib.colors.LinearSegmentedColormap.from_list(
            "yg", [cmap(0.35), cmap(0.95)])), ax=ax, pad=0.015, fraction=0.03)
        cb.set_label("Variance", fontsize=10)
        cb.outline.set_visible(False)
    return _save(fig, path)


# =============================================================================
# Metrics heatmap
# =============================================================================

def plot_heatmap(matrix, metric_names, varieties, samples, title, path):
    n = matrix.shape[1]
    width = min(17.0, max(9.5, 0.12 * n + 4.0))
    fig, ax = plt.subplots(figsize=(width, 5.6 + (1.2 if n <= 80 else 0.8)))
    im = ax.imshow(matrix, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1, interpolation="nearest")
    ax.set_yticks(range(len(metric_names)))
    ax.set_yticklabels(metric_names, fontsize=10)
    if n <= 80:
        ax.set_xticks(range(n))
        ax.set_xticklabels([f"{v} · {s}" for v, s in zip(varieties, samples)], rotation=90, fontsize=max(5, 9 - n // 20))
    else:                                              # one label per run of the same variety
        starts, labs = [], []
        prev = None
        for i, v in enumerate(varieties):
            if v != prev:
                starts.append(i)
                labs.append(v)
                prev = v
        ends = starts[1:] + [n]
        ax.set_xticks([(a + b - 1) / 2 for a, b in zip(starts, ends)])
        ax.set_xticklabels(labs, rotation=60, ha="right", fontsize=9)
        for s in starts[1:]:
            ax.axvline(s - 0.5, color="white", linewidth=1.0)
    ax.set_title(title, loc="left", fontsize=14, fontweight="bold", color=INK, pad=14)
    ax.set_xlabel("Sample" if n > 80 else "Variety · sample", fontsize=11)
    for side in ax.spines.values():
        side.set_visible(False)
    cb = fig.colorbar(im, ax=ax, pad=0.012, fraction=0.025)
    cb.set_label("Normalised value (0–1)", fontsize=10)
    cb.outline.set_visible(False)
    return _save(fig, path)


# =============================================================================
# Violin
# =============================================================================

def plot_violin(groups, names, title, ylabel, path):
    k = len(names)
    fig, ax = plt.subplots(figsize=(max(8.0, 1.3 * k + 3.0), 6.0))
    parts = ax.violinplot(groups, positions=range(k), widths=0.8, showextrema=False)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(color_for(i))
        body.set_edgecolor(color_for(i))
        body.set_alpha(0.55)
    bp = ax.boxplot(groups, positions=range(k), widths=0.14, patch_artist=True, showfliers=False,
                    medianprops=dict(color="white", linewidth=1.6), whiskerprops=dict(color="#444444"),
                    capprops=dict(color="#444444"))
    for patch in bp["boxes"]:
        patch.set(facecolor="#444444", edgecolor="#444444")
    for i, g in enumerate(groups):
        ax.hlines(np.mean(g), i - 0.2, i + 0.2, color=color_for(i), linewidth=2.2, zorder=4)
    ax.set_xticks(range(k))
    ax.set_xticklabels(names, rotation=0 if k <= 6 else 45, ha="center" if k <= 6 else "right", fontsize=10)
    _style_axes(ax, title, "Variety", ylabel)
    ax.grid(axis="x", visible=False)
    return _save(fig, path)


# =============================================================================
# Scatter (PCA)
# =============================================================================

def plot_scatter(series, title, xlabel, ylabel, legend_title, path):
    """series: list of (name, x, y, color)"""
    fig, ax = plt.subplots(figsize=(9.2, 6.4))
    for name, x, y, color in series:
        ax.scatter(x, y, s=62, color=color, edgecolor="white", linewidth=0.8, alpha=0.95, label=str(name), zorder=3)
    ax.axhline(0, color="#C9CED4", linewidth=0.9, zorder=1)
    ax.axvline(0, color="#C9CED4", linewidth=0.9, zorder=1)
    _style_axes(ax, title, xlabel, ylabel)
    ax.margins(0.07)
    _legend_outside(ax, legend_title, ncol=1 if len(series) <= 16 else 2)
    return _save(fig, path)


# =============================================================================
# Radar
# =============================================================================

def plot_radar(metrics, series, radial_max, radial_title, title, path):
    """series: list of (name, values, color)"""
    k = len(metrics)
    theta = np.linspace(0, 2 * np.pi, k, endpoint=False)
    fig = plt.figure(figsize=(9.2, 7.0))
    ax = fig.add_subplot(111, polar=True)
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    for name, vals, color in series:
        v = list(vals) + [vals[0]]
        t = list(theta) + [theta[0]]
        ax.plot(t, v, color=color, linewidth=2.0, label=str(name))
        ax.fill(t, v, color=color, alpha=0.16)
    ax.set_ylim(0, radial_max)
    ax.set_xticks(theta)
    ax.set_xticklabels([m.replace("_", " ") for m in metrics], fontsize=10, color="#333333")
    ax.tick_params(axis="x", pad=10)
    ax.tick_params(axis="y", labelsize=8, labelcolor="#777777")
    ax.grid(color="#D5D9DE", linewidth=0.9)
    ax.spines["polar"].set_color("#C9CED4")
    ax.set_facecolor("#FAFBFC")
    ax.set_title(title, fontsize=14, fontweight="bold", color=INK, pad=30)
    ax.text(0.5, -0.12, f"Radial axis: {radial_title}", transform=ax.transAxes, ha="center", fontsize=10, color="#666666")
    ax.legend(loc="upper left", bbox_to_anchor=(1.08, 1.02), frameon=False, fontsize=10, title="Variety", title_fontsize=11)
    return _save(fig, path)


# =============================================================================
# Circular neighbour-joining tree (iTOL style)
# =============================================================================

def plot_nj_tree(tree, sample_labels, variety_labels, path, support=None, n_boot=0, support_min=0.7):
    """
    tree: dict from diversity_engine.nj_structure (children, order, leaf_order, radius)
    support: {node: share of bootstrap replicates that contain the branch above the node}
    """
    n = len(sample_labels)
    children, order, leaf_order = tree["children"], tree["order"], tree["leaf_order"]
    radius = dict(tree["radius"])
    r_max = max(radius.values()) or 1.0
    radius = {k: v / r_max for k, v in radius.items()}

    gap = 2 * np.pi * 0.012
    step = (2 * np.pi - gap) / n
    angle = {leaf: i * step for i, leaf in enumerate(leaf_order)}
    for node in order:
        ch = children[node]
        if ch:
            angs = [angle[c] for c, _ in ch]
            angle[node] = (min(angs) + max(angs)) / 2.0

    varieties = list(dict.fromkeys(variety_labels))
    vcol = {v: color_for(i) for i, v in enumerate(varieties)}

    size = 11.0 if n <= 60 else 13.0 if n <= 150 else 16.0 if n <= 300 else 19.0
    font = 9.5 if n <= 40 else 8 if n <= 90 else 6.8 if n <= 160 else 5.6 if n <= 260 else 4.6
    show_labels = n <= 450
    ring_in = 1.03
    label_r = 1.06
    longest = max((len(str(s)) for s in sample_labels), default=1) if show_labels else 0
    label_in = longest * font * 0.58 / 72.0                         # inches the longest name needs
    lim = (label_r + 0.05) / max(1.0 - label_in / (0.46 * size), 0.35)

    fig_w = size * 1.22                                   # a margin on the right holds the legend
    fig = plt.figure(figsize=(fig_w, size))
    ax = fig.add_axes([0.01, 0.02, 0.96 * size / fig_w, 0.92])
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axis("off")

    # pale sector behind the names of every run of the same variety
    r_out = lim * 0.985
    i = 0
    while i < n:
        v = variety_labels[leaf_order[i]]
        j = i
        while j < n and variety_labels[leaf_order[j]] == v:
            j += 1
        a0 = np.degrees(angle[leaf_order[i]] - step * 0.5)
        a1 = np.degrees(angle[leaf_order[j - 1]] + step * 0.5)
        ax.add_patch(Wedge((0, 0), r_out, a0, a1, width=r_out - ring_in,
                           facecolor=to_rgb(vcol[v]), alpha=0.20, edgecolor="none", zorder=0))
        i = j

    # concentric distance circles
    dmax = r_max
    target = dmax * 0.15
    mag = 10 ** np.floor(np.log10(target)) if target > 0 else 1.0
    bar = min((m for m in (1, 2, 5, 10) if m * mag >= target), default=10) * mag
    k = 1
    while k * bar / dmax <= 1.0 + 1e-9:
        ax.add_patch(Circle((0, 0), k * bar / dmax, fill=False, edgecolor="#9A9AEE", linewidth=0.7, zorder=1))
        k += 1
    ax.add_patch(Circle((0, 0), ring_in - 0.012, fill=False, edgecolor="#6C6CD8", linewidth=1.0, zorder=1))

    # branches
    segs = []
    for node in order:
        ch = children[node]
        if not ch:
            continue
        angs = [angle[c] for c, _ in ch]
        if len(ch) > 1:
            lo, hi = min(angs), max(angs)
            t = np.linspace(lo, hi, max(2, int((hi - lo) / np.radians(1.5)) + 2))
            segs.append(np.column_stack([radius[node] * np.cos(t), radius[node] * np.sin(t)]))
        for c, _ in ch:
            a = angle[c]
            segs.append(np.array([[radius[node] * np.cos(a), radius[node] * np.sin(a)],
                                  [radius[c] * np.cos(a), radius[c] * np.sin(a)]]))
    ax.add_collection(LineCollection(segs, colors="#111111", linewidths=0.9 if n <= 150 else 0.7, zorder=3,
                                     capstyle="round", joinstyle="round"))

    # coloured dots on the tips (a dotted guide from the tip to the ring keeps short branches readable)
    tips = np.array([[radius[l] * np.cos(angle[l]), radius[l] * np.sin(angle[l])] for l in leaf_order])
    ax.scatter(tips[:, 0], tips[:, 1], s=8 if n <= 150 else 3, c=[vcol[variety_labels[l]] for l in leaf_order],
               zorder=5, linewidths=0)

    # bootstrap support
    n_sup = 0
    if support:
        xs, ys, ss = [], [], []
        for node, s in support.items():
            if s >= support_min and node in angle:
                xs.append(radius[node] * np.cos(angle[node]))
                ys.append(radius[node] * np.sin(angle[node]))
                ss.append(s)
        n_sup = len(xs)
        if xs:
            ax.scatter(xs, ys, marker="^", s=[(10 + 34 * (s - support_min) / (1 - support_min + 1e-9)) * (1.0 if n <= 150 else 0.6)
                                                for s in ss],
                       c="#E02424", edgecolors="none", zorder=6)

    # sample names
    if show_labels:
        for leaf in leaf_order:
            deg = np.degrees(angle[leaf])
            flip = 90 < deg < 270
            ax.text(label_r * np.cos(angle[leaf]), label_r * np.sin(angle[leaf]), str(sample_labels[leaf]),
                    rotation=(deg + 180) if flip else deg, rotation_mode="anchor",
                    ha="right" if flip else "left", va="center", fontsize=font, color="#1A1A1A", zorder=7)

    # scale bar + legend
    ax.plot([-lim * 0.97, -lim * 0.97 + bar / dmax], [-lim * 0.97, -lim * 0.97], color="#222222", linewidth=2.0, zorder=8)
    ax.text(-lim * 0.97, -lim * 0.97 - 0.015 * lim, f"{bar:g} = genetic distance (circles every {bar:g})",
            ha="left", va="top", fontsize=9, color="#444444", zorder=8)

    handles = [Line2D([0], [0], marker="s", linestyle="", markersize=10, markerfacecolor=to_rgb(vcol[v]) + (0.45,),
                      markeredgecolor=vcol[v], label=str(v)) for v in varieties]
    if n_sup:
        handles.append(Line2D([0], [0], marker="^", linestyle="", markersize=9, markerfacecolor="#E02424",
                              markeredgecolor="none", label=f"Bootstrap ≥ {int(support_min * 100)}%  ({n_boot} replicates)"))
    leg = fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.8, 0.9), frameon=True, framealpha=0.92,
                    edgecolor="#CCCCCC", fontsize=9 if len(handles) <= 14 else 7.5, title="Variety / Cultivar",
                    title_fontsize=10, ncol=1 if len(handles) <= 24 else 2)
    leg.set_zorder(20)
    fig.text(0.5 * 0.96 / 1.22 + 0.01, 0.965, "Neighbour-joining tree", ha="center", va="center", fontsize=17, fontweight="bold", color=INK)
    fig.savefig(str(path), dpi=DPI if size <= 13 else 170, facecolor="white")
    plt.close(fig)
    return str(path)
