"""
Figures for the GWAS tool: Manhattan, QQ, p-value histograms, circos, principal components, scree and the PC-count comparison.

Plots only draw what the analysis has already decided: a marker is highlighted when it passes its model's criterion, and the
top-ranked markers of a scan that pass nothing are drawn as open circles labelled exploratory.
"""

import io
import os
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from pycirclize import Circos
from scipy import stats

# ─────────────────────── DEFAULT COLOR PALETTES ──────────────────────────────────
BG_DARK     = "#0D1117"
PANEL_DARK  = "#161B22"
PANEL_MID   = "#21262D"
ACCENT_TEAL = "#2DD4BF"
ACCENT_GOLD = "#F59E0B"
ACCENT_CORAL= "#F87171"
ACCENT_BLUE = "#60A5FA"
ACCENT_LIME = "#84CC16"
ACCENT_PURP = "#A78BFA"
ACCENT_PINK = "#F472B6"
ACCENT_CYAN = "#22D3EE"
ACCENT_ORG  = "#FB923C"
TEXT_MAIN   = "#E6EDF3"
TEXT_DIM    = "#8B949E"
GRID_LINE   = "#30363D"
SIG_COLOR   = "#F59E0B"

# Light-theme equivalents -- used for ALL plots (PCA must NOT default to dark)
BG_LIGHT      = "#FFFFFF"
PANEL_LIGHT   = "#FAFBFF"
TEXT_DARK     = "#2D3142"
GRID_LIGHT    = "#E8ECF4"
ACCENT_TEAL_D = "#0F766E"
ACCENT_GOLD_D = "#B45309"
ACCENT_PINK_D = "#BE185D"

MODEL_COLORS = {
    "OLS":        "#60A5FA",
    "MLM":        "#2DD4BF",
    "EMMAX":      "#F472B6",
    "FarmCPU":    "#84CC16",
    "GEMMA":      "#E879F9",
    "BLINK":      "#FB923C",
    "mrMLM":      "#A78BFA",
    "FASTmrMLM":  "#F59E0B",
}

CHR_PALETTE = [
    "#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
    "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC",
    "#499894","#86BCB6","#D4A6C8","#FFBE7D","#8CD17D",
]

ALL_MODELS = ["OLS", "MLM", "EMMAX", "FarmCPU", "GEMMA", "BLINK", "mrMLM", "FASTmrMLM"]

def style_fig_white(fig, font_size=9, font_color="#2D3142", grid_color="#E8ECF4"):
    fig.patch.set_facecolor("#FFFFFF")
    for ax in fig.get_axes():
        ax.set_facecolor("#FAFBFF")
        ax.tick_params(colors=font_color, labelsize=font_size)
        ax.xaxis.label.set_color(font_color)
        ax.yaxis.label.set_color(font_color)
        ax.title.set_color(font_color)
        for spine in ax.spines.values():
            spine.set_edgecolor("#CCCCCC")
        ax.grid(color=grid_color, linewidth=0.5, linestyle="--", alpha=0.7)
    return fig

def apply_custom_style(fig, cfg):
    bg = cfg.get("bg_color", "#FFFFFF")
    panel = cfg.get("panel_color", "#FAFBFF")
    fc = cfg.get("font_color", "#2D3142")
    fs = cfg.get("font_size", 9)
    gc = cfg.get("grid_color", "#E8ECF4")
    fig.patch.set_facecolor(bg)
    for ax in fig.get_axes():
        ax.set_facecolor(panel)
        ax.tick_params(colors=fc, labelsize=fs)
        ax.xaxis.label.set_color(fc)
        ax.yaxis.label.set_color(fc)
        ax.title.set_color(fc)
        for spine in ax.spines.values():
            spine.set_edgecolor("#CCCCCC")
        ax.grid(color=gc, linewidth=0.5, linestyle="--", alpha=0.7)
        ax.xaxis.label.set_fontsize(fs)
        ax.yaxis.label.set_fontsize(fs)
        ax.title.set_fontsize(fs + 2)
    return fig


RENDER_DPI = 300  # publication-quality
EXPORT_DPI = 400  # all downloadable PNG/figure files must exceed 300 dpi

def fig_to_png_bytes(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=RENDER_DPI, bbox_inches="tight", facecolor=fig.get_facecolor())
    buf.seek(0)
    return buf.read()

def fig_to_pil(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=RENDER_DPI, bbox_inches="tight", facecolor=fig.get_facecolor())
    buf.seek(0)
    from PIL import Image
    return Image.open(buf).copy()

def build_cfg(w, h, fs, fc, gc="#E8ECF4", bg="#FFFFFF", panel="#FAFBFF",
              sig_col="#E15759", sug_col="#4E79A7", marker_size=12, alpha=0.75,
              font_color=None, dot_col=None, line_col=None, fill_col=None,
              axis_col=None):
    """
    Build a plot configuration dictionary.
    'fc' is the required positional font-color argument.
    'font_color' is an optional alias (takes precedence if provided).
    """
    effective_fc = font_color if font_color is not None else fc
    return {
        "width":        w,
        "height":       h,
        "font_size":    fs,
        "font_color":   effective_fc,
        "grid_color":   gc,
        "bg_color":     bg,
        "panel_color":  panel,
        "sig_color":    sig_col,
        "sug_color":    sug_col,
        "marker_size":  marker_size,
        "alpha":        alpha,
        "dot_color":    dot_col  or "#4E79A7",
        "line_color":   line_col or "#E15759",
        "fill_color":   fill_col or "#A78BFA",
        "axis_color":   axis_col or "#2D3142",
    }

def _nature_qq_ticks(max_val):
    """Pick clean, evenly-spaced tick values (Nature-style QQ plots use round
    steps like 0,3,6,9,12 rather than whatever matplotlib defaults to)."""
    if not np.isfinite(max_val) or max_val <= 0:
        return [0, 1]
    for step in (1, 2, 3, 5, 10, 15, 20, 25, 50, 100):
        n_ticks = max_val / step
        if n_ticks <= 5:
            top = int(np.ceil(max_val / step)) * step
            return list(np.arange(0, top + step, step))
    top = int(np.ceil(max_val / 100)) * 100
    return list(np.arange(0, top + 100, 100))

def neglog10p_to_lod(neglog10p):
    """Convert a -log10(p) value (1 d.f. test) to its equivalent LOD score.

    LOD = chi2 / (2*ln(10)), where chi2 is the 1-d.f. chi-square value
    whose upper-tail p-value equals 10**(-neglog10p). This is the
    standard linkage/QTL LOD scale, so a p-based threshold can be read
    directly as a LOD threshold as well.
    """
    neglog10p = np.atleast_1d(np.asarray(neglog10p, dtype=float))
    p = np.clip(10.0 ** (-neglog10p), 1e-300, 1.0)
    chi2_val = stats.chi2.isf(p, df=1)
    lod = chi2_val / (2.0 * np.log(10.0))
    return lod if lod.shape != (1,) else lod[0]

def lod_to_neglog10p(lod):
    """Exact inverse of neglog10p_to_lod: given a LOD score, return the
    -log10(p) value it corresponds to, so a LOD threshold can be drawn on
    a -log10(p)-scaled Manhattan/Circos axis."""
    lod = np.atleast_1d(np.asarray(lod, dtype=float))
    chi2_val = lod * 2.0 * np.log(10.0)
    p = stats.chi2.sf(chi2_val, df=1)
    neglog10p = -np.log10(np.clip(p, 1e-300, 1.0))
    return neglog10p if neglog10p.shape != (1,) else neglog10p[0]



# ─────────────────────── LAYOUT HELPERS ──────────────────────────────────────────
def chromosome_layout(chroms, positions):
    """x coordinate of every marker on a genome-wide axis, plus tick positions and labels."""
    chroms = np.asarray(chroms)
    positions = np.asarray(positions, dtype=float)
    chrom_list = sorted(set(chroms.tolist()))
    gap = positions.max() * 0.018 + 1
    x = np.zeros(len(chroms))
    offset = 0.0
    ticks, labels = [], []
    for ch in chrom_list:
        mask = chroms == ch
        pos = positions[mask]
        x[mask] = pos - pos.min() + offset
        ticks.append(float(x[mask].mean()))
        labels.append(f"Chr{ch}")
        offset += float(pos.max() - pos.min()) + gap
    return x, ticks, labels, offset, chrom_list


def _neglog10(p):
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p) & (p > 0)
    out[ok] = -np.log10(p[ok])
    return out


def _draw_manhattan(ax, x, chroms, chrom_list, p, color, thr, genome, explor, cfg, show_bonferroni=True, show_suggestive=True,
                    alternate=True):
    lp = _neglog10(p)
    for i, ch in enumerate(chrom_list):
        mask = (chroms == ch) & np.isfinite(lp)
        c = color if not alternate else (color if i % 2 == 0 else cfg.get("alt_color", "#7C8AA5"))
        ax.scatter(x[mask], lp[mask], c=c, s=cfg["marker_size"], alpha=cfg["alpha"], linewidths=0)
    if show_suggestive and np.isfinite(thr.get("suggestive_p", np.nan)):
        v = -np.log10(thr["suggestive_p"])
        ax.axhline(v, color=cfg["sug_color"], ls=":", lw=1.4, alpha=0.9, label=f"Suggestive (1/m = {thr['suggestive_p']:.1e})")
    if show_bonferroni and np.isfinite(thr.get("bonferroni_p", np.nan)):
        v = -np.log10(thr["bonferroni_p"])
        ax.axhline(v, color=cfg["sig_color"], ls="--", lw=1.8, alpha=0.95, label=f"Bonferroni (0.05/m = {thr['bonferroni_p']:.1e})")
    if thr.get("permutation_p"):
        ax.axhline(-np.log10(thr["permutation_p"]), color=ACCENT_LIME, ls="-", lw=1.6, alpha=0.95,
                   label=f"Simulated null threshold (p = {thr['permutation_p']:.1e})")
    if explor is not None and explor.any():
        ax.scatter(x[explor], lp[explor], facecolors="none", edgecolors="#555555", s=max(cfg["marker_size"] * 3, 30), linewidths=0.9,
                   zorder=8, label=f"Exploratory, not significant ({int(explor.sum())})")
    if genome is not None and genome.any():
        ax.scatter(x[genome], lp[genome], c=ACCENT_GOLD, s=max(cfg["marker_size"] * 4, 55), marker="D", zorder=10,
                   edgecolors="#333333", linewidth=0.7, label=f"Genome-wide / QTN ({int(genome.sum())})")
    elif genome is not None:
        ax.text(0.01, 0.97, "No genome-wide significant associations detected", transform=ax.transAxes, va="top", ha="left",
                fontsize=cfg["font_size"], color="#8A1F1F", fontweight="bold")


PASTEL_CHR = ["#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
              "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"]


def _axis_ticks(vmax, n=4):
    """Round tick values 0..top for an axis that has to reach vmax; top is the first round value >= vmax."""
    vmax = max(float(vmax), 1e-9)
    raw = vmax / n
    mag = 10 ** np.floor(np.log10(raw))
    step = next(m * mag for m in (1, 2, 5, 10) if raw <= m * mag * 1.0001)
    top = np.ceil(vmax / step - 1e-9) * step
    return [float(v) for v in np.arange(0, top + step / 2, step)], float(top)


def _fmt_tick(v):
    return f"{v:g}"


def _spread_labels(xs, lo, hi, sep):
    """Horizontal positions for labels that must be at least `sep` apart and stay inside [lo, hi] (xs sorted)."""
    out = list(xs)
    for k in range(1, len(out)):
        out[k] = max(out[k], out[k - 1] + sep)
    if out and out[-1] > hi:
        out[-1] = hi
        for k in range(len(out) - 2, -1, -1):
            out[k] = min(out[k], out[k + 1] - sep)
    return [max(v, lo) for v in out]


def plot_manhattan(p, chroms, positions, model_name="", thresholds=None, genome=None, exploratory=None, cfg=None,
                   show_bonferroni=True, show_suggestive=True, title="Manhattan plot", snp_names=None, max_labels=12):
    """-log10(p) of every marker, one colour per chromosome, a marker-density strip with the chromosome names underneath, the
    significance lines, and the name of each genome-wide / QTN marker above the plot. The y axis ends at the highest point or
    threshold of THIS model, so every model is shown at its own scale. Markers that pass nothing are never highlighted; the
    top-ranked ones that pass nothing are open circles."""
    cfg = cfg or build_cfg(14, 5.6, 10, "#222222")
    thr = thresholds or {}
    chroms = np.asarray(chroms)
    x, ticks, labels, total, chrom_list = chromosome_layout(chroms, positions)
    lp = _neglog10(p)
    fs = cfg["font_size"]
    lines = []
    if show_suggestive and np.isfinite(thr.get("suggestive_p", np.nan)):
        lines.append((-np.log10(thr["suggestive_p"]), "--", "#555555", 1.0, "suggestive"))
    if show_bonferroni and np.isfinite(thr.get("bonferroni_p", np.nan)):
        lines.append((-np.log10(thr["bonferroni_p"]), "-", "#222222", 1.0, "Bonferroni"))
    if thr.get("permutation_p"):
        lines.append((-np.log10(thr["permutation_p"]), "-", "#3A8F3A", 1.0, "simulated null"))
    top_data = float(np.nanmax(lp)) if np.isfinite(lp).any() else 1.0
    yticks, ytop = _axis_ticks(max([top_data] + [v for v, *_ in lines]) * 1.02, 4)

    fig, (ax, st) = plt.subplots(2, 1, figsize=(cfg["width"], cfg["height"]), sharex=True,
                                 gridspec_kw={"height_ratios": [1, 0.085], "hspace": 0.035})
    fig.patch.set_facecolor(cfg["bg_color"])
    n = int(np.isfinite(lp).sum())
    s = float(np.clip(60000.0 / max(n, 1), 5, 26))
    for i, ch in enumerate(chrom_list):
        m = (chroms == ch) & np.isfinite(lp)
        ax.scatter(x[m], lp[m], c=PASTEL_CHR[i % len(PASTEL_CHR)], s=s, alpha=0.95, linewidths=0, rasterized=True)
    for v, ls, col, lw, name in lines:
        ax.axhline(v, color=col, ls=ls, lw=lw)
        ax.text(total * 1.006, v, name, ha="left", va="center", fontsize=fs - 3, color=col, clip_on=False)
    if exploratory is not None and exploratory.any():
        ax.scatter(x[exploratory], lp[exploratory], facecolors="none", edgecolors="#555555", s=max(s * 2.5, 30), linewidths=0.9, zorder=8)
    hit = np.zeros(len(lp), dtype=bool) if genome is None else np.asarray(genome, dtype=bool) & np.isfinite(lp)
    if hit.any():
        ax.scatter(x[hit], lp[hit], c="#E4572E", s=max(s * 1.8, 26), edgecolors="#333333", linewidths=0.6, zorder=10)
    ax.set_xlim(-total * 0.004, total * 1.004)
    ax.set_ylim(0, ytop)
    ax.set_yticks(yticks)
    ax.set_yticklabels([_fmt_tick(v) for v in yticks], fontsize=fs, fontweight="bold", color=cfg["font_color"])
    ax.set_ylabel(f"{model_name}\n−log₁₀(p)" if model_name else "−log₁₀(p)", fontsize=fs + 2, color=cfg["font_color"])
    ax.set_facecolor("none")
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.tick_params(axis="x", length=0, labelbottom=False)
    ax.grid(False)

    if hit.any() and snp_names is not None:
        # one name per cluster of hits: the strongest marker of each neighbourhood
        chosen = []
        for i in np.flatnonzero(hit)[np.argsort(-lp[hit])]:
            if all(abs(x[i] - x[j]) > total * 0.03 for j in chosen):
                chosen.append(i)
            if len(chosen) >= max_labels:
                break
        idx = np.array(sorted(chosen, key=lambda i: x[i]), dtype=int)
        lx = _spread_labels([x[i] for i in idx], 0, total, total * 0.02 * fs / 10)
        for i, lab_x in zip(idx, lx):
            ax.plot([x[i], x[i], lab_x], [lp[i], ytop * 1.035, ytop * 1.075], color="#444444", lw=0.6, clip_on=False, zorder=1)
            ax.text(lab_x, ytop * 1.08, str(snp_names[i]), rotation=90, ha="center", va="bottom", fontsize=fs, color="#111111")
    elif not hit.any() and genome is not None:
        ax.text(0.01, 0.97, "No genome-wide significant associations detected", transform=ax.transAxes, va="top", ha="left",
                fontsize=fs, color="#8A1F1F", fontweight="bold")

    st.set_facecolor("none")
    for side in ("top", "right", "left", "bottom"):
        st.spines[side].set_visible(False)
    st.set_yticks([])
    st.set_ylim(0, 1)
    for ch in chrom_list:
        m = chroms == ch
        xm = x[m]
        st.add_patch(plt.Rectangle((xm.min(), 0), max(xm.max() - xm.min(), total * 0.002), 1, fill=False, ec="#111111", lw=0.9))
        st.vlines(xm, 0.04, 0.96, color="#111111", lw=0.3)
    st.set_xticks(ticks)
    st.set_xticklabels([l.replace("Chr", "chr") for l in labels], fontsize=fs - 1 if len(chrom_list) <= 16 else fs - 3,
                       rotation=0 if len(chrom_list) <= 16 else 90, color=cfg["font_color"])
    st.tick_params(axis="x", length=0, pad=4)
    return fig


def plot_model_comparison_manhattan(items, chroms, positions, cfg=None):
    """One panel per model. items: list of (name, p, thresholds, genome_mask, exploratory_mask)."""
    cfg = cfg or build_cfg(16, 3.5, 8, "#2D3142")
    chroms = np.asarray(chroms)
    x, ticks, labels, total, chrom_list = chromosome_layout(chroms, positions)
    n = len(items)
    fig, axes = plt.subplots(n, 1, figsize=(cfg["width"], cfg["height"] * n), sharex=True)
    fig.patch.set_facecolor(cfg["bg_color"])
    axes = [axes] if n == 1 else list(axes)
    for k, (name, p, thr, genome, explor) in enumerate(items):
        ax = axes[k]
        ax.set_facecolor(cfg["panel_color"])
        col = MODEL_COLORS.get(name, CHR_PALETTE[k % len(CHR_PALETTE)])
        _draw_manhattan(ax, x, chroms, chrom_list, p, col, thr, genome, explor, cfg)
        ax.set_ylabel("−log₁₀(p)", fontsize=cfg["font_size"], color=cfg["font_color"])
        ax.set_title(f"{name}", color=col, fontsize=cfg["font_size"] + 1, pad=3, fontweight="bold")
        ax.set_ylim(bottom=0)
        ax.tick_params(colors=cfg["font_color"], labelsize=cfg["font_size"] - 1)
        h, _ = ax.get_legend_handles_labels()
        if h:
            ax.legend(fontsize=cfg["font_size"] - 2, framealpha=0.85, facecolor=cfg["bg_color"], edgecolor="#CCCCCC",
                      labelcolor=cfg["font_color"], ncol=3, loc="upper right")
        for sp in ax.spines.values():
            sp.set_edgecolor("#CCCCCC")
        ax.grid(color=cfg["grid_color"], linewidth=0.3, linestyle="--", alpha=0.5)
    axes[-1].set_xticks(ticks)
    axes[-1].set_xticklabels(labels, rotation=40, fontsize=cfg["font_size"] - 1, color=cfg["font_color"])
    axes[-1].set_xlim(-total * 0.01, total * 1.01)
    fig.suptitle("Manhattan plots — all models", color=cfg["font_color"], fontsize=cfg["font_size"] + 3, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    return fig


# ─────────────────────── QQ PLOT ─────────────────────────────────────────────────
def plot_qq(curves, lambdas=None, title="QQ plot", cfg=None):
    """Observed against expected -log10(p). curves: {name: p array}; lambdas: {name: lambda}. The genomic inflation factor is
    printed beside each name as a number; nothing is rescaled. Each axis stops at its own largest value and the figure takes
    the proportions of the data (equal units on both axes, so the diagonal stays at 45 degrees), so there is no empty margin
    from a fixed square."""
    from matplotlib.ticker import MaxNLocator
    cfg = cfg or build_cfg(6, 6, 10, "#000000")
    prepared = []
    for name, p in curves.items():
        p = np.asarray(p, float)
        p = np.sort(p[np.isfinite(p) & (p > 0) & (p <= 1)])
        if len(p) == 0:
            continue
        n = len(p)
        prepared.append((name, -np.log10((np.arange(1, n + 1) - 0.5) / n), -np.log10(p)))
    xmax = max([float(e.max()) for _, e, _ in prepared] + [1.0])
    ymax = max([float(o.max()) for _, _, o in prepared] + [1.0])
    xlim, ylim = xmax * 1.04, ymax * 1.04
    base = 5.6
    w, h = (base, base * ylim / xlim) if xlim >= ylim else (base * xlim / ylim, base)
    equal = min(w, h) >= 2.8
    w, h = max(w, 2.8), max(h, 2.8)
    fig, ax = plt.subplots(figsize=(w + 1.0, h + 0.9))
    fig.patch.set_facecolor("#FFFFFF")
    ax.set_facecolor("#FFFFFF")
    d = min(xlim, ylim)
    ax.plot([0, d], [0, d], color="#000000", lw=1.2, zorder=3)
    handles = []
    for k, (name, exp, obs) in enumerate(prepared):
        col = MODEL_COLORS.get(name, CHR_PALETTE[k % len(CHR_PALETTE)])
        ax.scatter(exp, obs, s=max(cfg["marker_size"] * 0.45, 6), color=col, alpha=0.9, linewidths=0, zorder=5)
        lam = None if lambdas is None else lambdas.get(name)
        handles.append(Line2D([0], [0], marker="o", color="none", markerfacecolor=col, markersize=7,
                              label=name + (f"  λ = {lam:.3f}" if lam is not None and np.isfinite(lam) else "")))
    ax.set_xlim(0, xlim)
    ax.set_ylim(0, ylim)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
    if equal:
        ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Expected −log$_{10}$($P$ value)", fontsize=cfg["font_size"] + 1)
    ax.set_ylabel("Observed −log$_{10}$($P$ value)", fontsize=cfg["font_size"] + 1)
    ax.set_title(title, fontsize=cfg["font_size"] + 1, color="#333333")
    if handles:
        ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=cfg["font_size"] - 1, handletextpad=0.2,
                  borderaxespad=0.2)
    ax.tick_params(colors="#000000", labelsize=cfg["font_size"], length=4, width=1.1)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(False)
    fig.tight_layout(pad=0.6)
    return fig


def plot_pvalue_histograms(items, cfg=None):
    """items: list of (name, p). One histogram per model with the uniform expectation drawn as a reference line."""
    cfg = cfg or build_cfg(4, 3, 9, "#2D3142")
    n = len(items)
    ncols = min(4, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(cfg["width"] * ncols, cfg["height"] * nrows), squeeze=False)
    fig.patch.set_facecolor("#FFFFFF")
    for k, (name, p) in enumerate(items):
        ax = axes[k // ncols][k % ncols]
        p = np.asarray(p, float)
        p = p[np.isfinite(p) & (p >= 0) & (p <= 1)]
        col = MODEL_COLORS.get(name, "#4E79A7")
        ax.hist(p, bins=20, range=(0, 1), color=col, alpha=0.85, edgecolor="white")
        ax.axhline(len(p) / 20.0, color="#333333", lw=1.0, ls="--")
        ax.set_title(f"{name}  (n = {len(p):,})", fontsize=cfg["font_size"] + 1, color=col, fontweight="bold")
        ax.set_xlabel("p-value", fontsize=cfg["font_size"])
        ax.set_ylabel("Markers", fontsize=cfg["font_size"])
        ax.tick_params(labelsize=cfg["font_size"] - 1)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    for k in range(n, nrows * ncols):
        axes[k // ncols][k % ncols].set_visible(False)
    fig.suptitle("p-value distributions (dashed line = uniform expectation)", fontsize=cfg["font_size"] + 2, y=1.0)
    fig.tight_layout()
    return fig


# ─────────────────────── POPULATION STRUCTURE ────────────────────────────────────
def plot_scree(pc, n_used=0, n_show=20, cfg=None):
    cfg = cfg or build_cfg(7, 4.2, 9, TEXT_DARK)
    ve = np.asarray(pc["variance_explained"], float)[:n_show] * 100
    k = np.arange(1, len(ve) + 1)
    fig, ax = plt.subplots(figsize=(cfg["width"], cfg["height"]))
    fig.patch.set_facecolor("#FFFFFF")
    ax.bar(k, ve, color=[ACCENT_TEAL_D if i <= n_used else "#B7C0CC" for i in k], edgecolor="white")
    ax2 = ax.twinx()
    ax2.plot(k, np.cumsum(ve), color=ACCENT_GOLD_D, marker="o", ms=3.5, lw=1.4)
    ax2.set_ylabel("Cumulative (%)", color=ACCENT_GOLD_D, fontsize=cfg["font_size"])
    ax2.set_ylim(0, 100)
    ax.set_xlabel("Principal component", fontsize=cfg["font_size"])
    ax.set_ylabel("Variance explained (%)", fontsize=cfg["font_size"])
    sub = f"{pc['n_snps']:,} markers" + (", LD-pruned" if pc.get("ld_pruned") else "")
    ax.set_title(f"Scree plot ({sub})" + (f" — dark bars: used as covariates by some model" if n_used else ""),
                 fontsize=cfg["font_size"] + 1, color="#333333")
    ax.set_xticks(k[::max(len(k) // 10, 1)])
    for side in ("top",):
        ax.spines[side].set_visible(False)
        ax2.spines[side].set_visible(False)
    fig.tight_layout()
    return fig


def plot_pca_scatter(pc, y, groups=None, group_source="", cfg=None):
    """PC1 against PC2 and PC2 against PC3, coloured by subpopulation label (or k-means group) when there is one, otherwise by
    the phenotype."""
    cfg = cfg or build_cfg(13, 5.4, 9, TEXT_DARK, marker_size=45, alpha=0.85)
    S = pc["scores"]
    ve = np.asarray(pc["variance_explained"], float) * 100
    pairs = [(0, 1)] + ([(1, 2)] if S.shape[1] >= 3 else [])
    fig, axes = plt.subplots(1, len(pairs), figsize=(cfg["width"] * len(pairs) / 2, cfg["height"]), squeeze=False)
    fig.patch.set_facecolor("#FFFFFF")
    for ax, (a, b) in zip(axes[0], pairs):
        ax.set_facecolor(PANEL_LIGHT)
        if groups is not None:
            for k, g in enumerate(sorted(set(groups.tolist()))):
                m = groups == g
                ax.scatter(S[m, a], S[m, b], s=cfg["marker_size"], alpha=cfg["alpha"], color=CHR_PALETTE[k % len(CHR_PALETTE)],
                           edgecolors="white", linewidths=0.4, label=f"{g} ({int(m.sum())})")
        else:
            sc = ax.scatter(S[:, a], S[:, b], c=y, cmap="viridis", s=cfg["marker_size"], alpha=cfg["alpha"], edgecolors="white",
                            linewidths=0.4)
            fig.colorbar(sc, ax=ax).set_label("Phenotype", fontsize=cfg["font_size"] - 1)
        ax.set_xlabel(f"PC{a + 1} ({ve[a]:.1f}%)", fontsize=cfg["font_size"])
        ax.set_ylabel(f"PC{b + 1} ({ve[b]:.1f}%)", fontsize=cfg["font_size"])
        ax.grid(color=GRID_LIGHT, linewidth=0.5, linestyle="--", alpha=0.7)
        for sp in ax.spines.values():
            sp.set_edgecolor(GRID_LIGHT)
    if groups is not None:
        axes[0][0].legend(fontsize=cfg["font_size"] - 1, frameon=False, title=group_source or None, title_fontsize=cfg["font_size"] - 2)
    fig.suptitle("Population structure", fontsize=cfg["font_size"] + 3, fontweight="bold", color=TEXT_DARK)
    fig.tight_layout()
    return fig


def plot_pc_comparison(rows, curves, cfg=None):
    """Left: genomic inflation factor against the number of PCs used as covariates. Right: QQ curves for each choice."""
    cfg = cfg or build_cfg(13, 5.2, 9, TEXT_DARK)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(cfg["width"], cfg["height"]))
    fig.patch.set_facecolor("#FFFFFF")
    for model, col in (("OLS", MODEL_COLORS["OLS"]), ("MLM", MODEL_COLORS["MLM"])):
        pts = sorted((r["PCs"], r["Lambda"]) for r in rows if r["Model"] == model and r["Lambda"] is not None)
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], marker="o", color=col, label=model + (" (K + PCs)" if model == "MLM" else ""))
    ax.axhline(1.0, color="#333333", lw=1.0, ls="--")
    ax.set_xlabel("Principal components used as covariates", fontsize=cfg["font_size"])
    ax.set_ylabel("Genomic inflation factor λ", fontsize=cfg["font_size"])
    ax.legend(frameon=False, fontsize=cfg["font_size"] - 1)
    ax.set_title("λ by number of PCs", fontsize=cfg["font_size"] + 1)
    lim = 1.0
    for k, (name, (p, n_total, idx)) in enumerate(curves.items()):
        if len(p) == 0:
            continue
        obs = -np.log10(p)
        exp = -np.log10((idx + 0.5) / n_total)
        lim = max(lim, float(obs.max()), float(exp.max()))
        bx.scatter(exp, obs, s=7, alpha=0.85, linewidths=0, color=CHR_PALETTE[k % len(CHR_PALETTE)], label=name)
    bx.plot([0, lim], [0, lim], color="#000000", lw=1.0)
    bx.set_xlabel("Expected −log$_{10}$($P$ value)", fontsize=cfg["font_size"])
    bx.set_ylabel("Observed −log$_{10}$($P$ value)", fontsize=cfg["font_size"])
    bx.legend(frameon=False, fontsize=cfg["font_size"] - 2, markerscale=2)
    bx.set_title("QQ by covariate choice", fontsize=cfg["font_size"] + 1)
    for a in (ax, bx):
        for side in ("top", "right"):
            a.spines[side].set_visible(False)
    fig.tight_layout()
    return fig


# ─────────────────────── CIRCOS — MTAs (pyCirclize) ──────────────────────────────
def _nice_step(span, target=6):
    raw = span / target
    mag = 10 ** np.floor(np.log10(max(raw, 1e-9)))
    for m in (1, 2, 2.5, 5, 10):
        if raw <= m * mag:
            return m * mag
    return 10 * mag


# One colour per chromosome, cycled (the order of the reference figures).
CIRCOS_CHR_COLORS = ["#4A9BE0", "#F6C026", "#4B3F9F", "#4F5D2A", "#E0267A", "#EA7A1C", "#D77FB5", "#83D1AB", "#7A1230", "#1E6B52"]


def plot_circos_mta(all_results, chroms, positions, snp_names, ring_names, mtas,
                    thresholds=None, max_labels=12):
    """Circos figure of the associations: a black chromosome band, then one ring per model (outer ring first) showing -log10(p)
    of every marker in the colour of its chromosome. Each ring has its own y axis that ends at the highest value or threshold of
    that model, so a weak model is as readable as a strong one. Genome-wide significant markers (QTNs for multi-locus models)
    are red stars on their model's ring, the Bonferroni level is a dashed line, and the strongest markers are named outside.
    `mtas` are the genome-wide / QTN rows only; thresholds = {model: Bonferroni p}."""
    chroms = np.asarray(chroms).astype(int)
    positions = np.asarray(positions, dtype=float)
    chr_list = sorted(set(chroms.tolist()))
    pos = positions.copy()
    for c in chr_list:
        m = chroms == c
        pos[m] -= pos[m].min()
    length = {c: max(float(pos[chroms == c].max()), 1.0) for c in chr_list}
    models = [m for m in ring_names if m in all_results]
    circos = Circos({str(c): length[c] for c in chr_list}, space=2.0)
    thresholds = thresholds or {}

    for sector in circos.sectors:
        c = int(sector.name)
        t = sector.add_track((91.5, 100))
        t.axis(fc="#111111", ec="#111111", lw=0.8)
        t.text(f"Chr{c}", size=9 if len(chr_list) <= 14 else 7, color="#FFFFFF", weight="bold")

    n = max(len(models), 1)
    top, hole, gap = 88.5, 26.0, 4.0
    each = min((top - hole - gap * (n - 1)) / n, 24.0)
    mta_by_model = {m: [r for r in mtas if r["Model"] == m] for m in models}
    snp_pos = {str(s): i for i, s in enumerate(snp_names)}
    dot = float(np.clip(60000.0 / max(len(chroms), 1), 3, 14))
    chr_col = {c: CIRCOS_CHR_COLORS[i % len(CIRCOS_CHR_COLORS)] for i, c in enumerate(chr_list)}
    for k, name in enumerate(models):
        pvals = np.asarray(all_results[name], dtype=float)
        lp = np.nan_to_num(-np.log10(np.clip(pvals, 1e-300, 1)), nan=0.0)
        thr_negp = -np.log10(thresholds[name]) if thresholds.get(name) else None
        yticks, ymax = _axis_ticks(max(float(lp.max()), thr_negp or 0.0, 1.5) * 1.02, 4)
        r_top = top - k * (each + gap)
        r_bot = r_top - each
        hit = {r["SNP"] for r in mta_by_model[name]}
        for sector in circos.sectors:
            c = int(sector.name)
            idx = np.where(chroms == c)[0]
            tr = sector.add_track((r_bot, r_top))
            tr.axis(fc="#FFFFFF", ec="#CFCFCF", lw=0.4)
            for v in yticks[1:]:
                tr.line([0, length[c]], [v, v], vmin=0, vmax=ymax, color="#E2E2E2", lw=0.5)
            if sector is circos.sectors[0]:
                tr.yticks(yticks, labels=[_fmt_tick(v) for v in yticks], vmin=0, vmax=ymax, side="left", tick_length=1.2,
                          label_size=7, label_margin=0.8, text_kws=dict(color="#222222", fontweight="bold"),
                          line_kws=dict(color="#222222", lw=0.7))
            if idx.size:
                tr.scatter(pos[idx], lp[idx], vmin=0, vmax=ymax, s=dot, color=chr_col[c], ec="none", alpha=0.9)
            if thr_negp is not None and thr_negp <= ymax:
                tr.line([0, length[c]], [thr_negp, thr_negp], vmin=0, vmax=ymax, color="#C0392B", lw=0.8, ls=(0, (2, 2)))
            hit_idx = [i for i in idx if str(snp_names[i]) in hit]
            if hit_idx:
                tr.scatter(pos[hit_idx], lp[hit_idx], vmin=0, vmax=ymax, marker="*", s=40,
                           color="#E4572E", ec="#5A1A0C", lw=0.4, zorder=5)

    best = {}
    for r in mtas:
        if snp_pos.get(r["SNP"]) is None:
            continue
        cur = best.get(r["SNP"])
        if cur is None or (r["neg_log10P"] or 0) > (cur["neg_log10P"] or 0):
            best[r["SNP"]] = r
    placed = []
    for r in sorted(best.values(), key=lambda r: -(r["neg_log10P"] or 0)):
        i = snp_pos[r["SNP"]]
        c, p = int(chroms[i]), float(pos[i])
        if any(pc == c and abs(pp - p) < 0.08 * length[c] for pc, pp in placed):
            continue
        placed.append((c, p))
        circos.get_sector(str(c)).tracks[0].annotate(p, r["SNP"], label_size=7, shorten=22)
        if len(placed) >= max_labels:
            break

    if models:
        circos.text("Rings, outer to inner\n" + "\n".join(f"{k + 1}  {m}" for k, m in enumerate(models)), r=0, size=10)
    fig = circos.plotfig(figsize=(11, 11))
    fig.patch.set_facecolor("white")
    legend = [Line2D([0], [0], marker="*", color="none", markerfacecolor="#E4572E", markeredgecolor="#5A1A0C", markersize=11,
                     label="Genome-wide significant / QTN")]
    if thresholds:
        legend.append(Line2D([0], [0], color="#C0392B", lw=1.2, ls=(0, (2, 2)), label="Bonferroni level"))
    fig.legend(handles=legend, loc="upper left", ncol=1, frameon=False, fontsize=10, bbox_to_anchor=(0.0, 1.0))
    return fig


# ─────────────────────── R FIGURES (CMplot + base graphics) ──────────────────────────
R_PLOT_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gwas_plots.R")


def r_figures(chroms, positions, snp_names, results, thresholds, hits, lambdas, out_dir, show_bonferroni=True,
              show_suggestive=True, dpi=300, max_labels=12):
    """Draws the Manhattan (one per model), QQ and circular figures with R (gwas_plots.R) and returns {file name: path}.
    results: {model: p array}; thresholds: {model: {bonferroni_p, suggestive_p, permutation_p}}; hits: {model: [marker names]}.
    Raises when R, or a figure, is not available so the caller can fall back to the matplotlib figures."""
    import json
    import pandas as pd
    from tools.gwas.backend.r_packages import run_r

    models = list(results)
    df = pd.DataFrame({"SNP": np.asarray(snp_names).astype(str), "Chromosome": np.asarray(chroms).astype(int),
                       "Position": np.asarray(positions, dtype=np.int64)})
    for m in models:
        df[m] = np.asarray(results[m], dtype=float)
    csv = os.path.join(out_dir, "gwas_plot_input.csv")
    df.to_csv(csv, index=False)

    def _clean(t):
        return {k: float(v) for k, v in (t or {}).items() if v is not None and np.isfinite(v)}

    cfg = {"input_csv": csv.replace("\\", "/"), "out_dir": out_dir.replace("\\", "/"), "models": models,
           "thresholds": {m: _clean(thresholds.get(m)) for m in models},
           "hits": {m: [str(h) for h in hits.get(m, [])] for m in models},
           "lambdas": {m: float(v) for m, v in (lambdas or {}).items() if v is not None and np.isfinite(v)},
           "show_bonferroni": bool(show_bonferroni), "show_suggestive": bool(show_suggestive), "dpi": int(dpi),
           "max_labels": int(max_labels), "model_colors": [MODEL_COLORS.get(m, CHR_PALETTE[i % len(CHR_PALETTE)]) for i, m in enumerate(models)]}
    cfg_path = os.path.join(out_dir, "gwas_plot_config.json")
    with open(cfg_path, "w", encoding="utf-8") as f:
        json.dump(cfg, f)
    try:
        r = run_r([R_PLOT_SCRIPT, cfg_path], timeout=1800)
        if r.returncode != 0:
            raise RuntimeError((r.stderr or r.stdout or "R failed").strip().splitlines()[-1])
        line = [l for l in (r.stdout or "").splitlines() if l.startswith("{")][-1]
        files = json.loads(line)["files"]
    finally:
        for tmp in (csv, cfg_path):
            try:
                os.remove(tmp)
            except OSError:
                pass
    files = {os.path.basename(f): f for f in ([files] if isinstance(files, str) else files) if os.path.isfile(f)}
    from PIL import Image, ImageChops
    for f in files.values():                      # trim the empty margins the R devices leave round the figure
        try:
            im = Image.open(f).convert("RGB")
            box = ImageChops.difference(im, Image.new("RGB", im.size, (255, 255, 255))).getbbox()
            if box:
                pad = max(12, int(im.size[0] * 0.012))
                im.crop((max(box[0] - pad, 0), max(box[1] - pad, 0), min(box[2] + pad, im.size[0]), min(box[3] + pad, im.size[1]))).save(f, dpi=(dpi, dpi))
        except Exception:
            pass
    need = 1 + 1 + 1 if len(models) == 1 else len(models) + 2
    if len(files) < need:
        raise RuntimeError(f"R drew {len(files)} of {need} figures")
    return files


def save_png_zip(paths, out_dir):
    """gwas_plots.zip from figures that are already PNG files."""
    zip_path = os.path.join(out_dir, "gwas_plots.zip")
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, path in paths.items():
            zf.write(path, name if name.lower().endswith(".png") else name + ".png")
    return zip_path


def save_all_plots_zip(plot_dict, out_dir, dpi=EXPORT_DPI):
    zip_path = os.path.join(out_dir, "gwas_plots.zip")
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zf:
        for name, fig in plot_dict.items():
            if fig is not None:
                buf = io.BytesIO()
                fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight",
                            facecolor=fig.get_facecolor())
                buf.seek(0)
                zf.writestr(f"{name}.png", buf.read())
    return zip_path
