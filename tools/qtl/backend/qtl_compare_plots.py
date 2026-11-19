"""
qtl_compare_plots.py — Figures and JSON helpers for QTL "Compare Models"
and "Common QTL".

Drawn with matplotlib's object API (Figure + Agg canvas, no pyplot) so it
is safe to call from the app's worker threads.
"""

import json
import math
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

# One colour per model, the same in the ideogram and the bar chart of the comparison
MODEL_COLORS = {
    "SIM": "#FF69B4",       # pink
    "CIM": "#FF1A1A",       # bright red
    "MQM": "#FF8C00",       # orange
    "RF": "#8C8C8C",        # grey
    "SVR": "#7EC8F0",       # light blue
    "XGBOOST": "#7B2CBF",   # purple
}
BACKEND_DIR = Path(__file__).resolve().parent
CHR_FILL, CHR_EDGE, INK, MUTED = "#e3f1e9", "#2d6a4f", "#1a3a2a", "#5a7566"


# ───────────────────────── JSON safety ─────────────────────────

def to_native(obj):
    """numpy/pandas scalars -> plain Python; NaN/inf -> None. pandas hands
    back numpy int64 values, which json.dumps can't encode — that error is
    what stopped Compare Models / Common QTL results reaching the page."""
    if isinstance(obj, dict):
        return {str(k): to_native(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_native(v) for v in obj]
    if isinstance(obj, np.generic):
        obj = obj.item()
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    return obj


# ───────────────────────── data helpers ─────────────────────────

def _chr_key(c):
    try:
        return (0, float(c), "")
    except (TypeError, ValueError):
        return (1, 0.0, str(c))


def load_scans(model_results):
    """{model: full-scan DataFrame} from each model's qtl_full_scan_*.csv."""
    scans = {}
    for model, res in model_results.items():
        path = res.get("full_csv")
        if not path or not Path(path).exists():
            continue
        df = pd.read_csv(path)
        if {"Marker", "Chromosome", "Position", "LOD"} <= set(df.columns):
            df["LOD"] = pd.to_numeric(df["LOD"], errors="coerce")
            df["Position"] = pd.to_numeric(df["Position"], errors="coerce")
            df["Chromosome"] = df["Chromosome"].astype(str)
            scans[model] = df
    return scans


def _chromosome_lengths(scans):
    allpos = pd.concat([d[["Chromosome", "Position"]] for d in scans.values()], ignore_index=True)
    lengths = allpos.groupby("Chromosome")["Position"].max().fillna(0)
    order = sorted(lengths.index, key=_chr_key)
    return order, lengths


def _save(fig, path):
    FigureCanvasAgg(fig)
    fig.savefig(path, dpi=150, facecolor="white", bbox_inches="tight")
    return str(Path(path).resolve()).replace("\\", "/")


def _style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color("#b9d3c4")
    ax.spines["bottom"].set_color("#b9d3c4")
    ax.tick_params(colors=MUTED, labelsize=9)


def _draw_chromosomes(ax, order, lengths, width=0.56):
    for i, ch in enumerate(order):
        ax.bar(i, lengths[ch], width=width, bottom=0, color=CHR_FILL, edgecolor=CHR_EDGE, linewidth=1.2, zorder=1)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"Chr {c}" for c in order])
    ax.set_xlim(-0.7, len(order) - 0.3)
    ymax = max(float(lengths.max()), 1.0)
    ax.set_ylim(ymax * 1.04, -ymax * 0.04)   # 0 cM at the top, like a karyotype
    ax.set_ylabel("Position", color=MUTED)



# ───────────────────────── Compare Models ─────────────────────────

def _regions_by_model(model_results):
    from tools.qtl.backend.qtl_regions import load_regions
    reg = load_regions(model_results)
    reg["Chromosome"] = reg["Chromosome"].astype(str)
    return reg


def _r_ideogram(model_results, out_dir, plot_cfg, trait):
    """The chromosome map of the single-model run (same R function), with a lane and a colour per model.
    Returns the PNG path, or None when R is not available / fails (the caller then draws a plain version)."""
    from tools.common.r_locator import find_rscript, get_r_env
    from tools.common.subprocess_utils import no_window_kwargs

    rscript = find_rscript()
    if rscript is None:
        return None
    models = [{"name": m, "qtl_csv": r.get("qtl_csv"), "full_csv": r.get("full_csv"), "color": MODEL_COLORS.get(str(m).upper(), "#CC0000")}
              for m, r in model_results.items() if r.get("full_csv")]
    if not models:
        return None
    out = Path(out_dir) / "compare_ideogram.png"
    with tempfile.TemporaryDirectory() as tmp:
        cfg = Path(tmp) / "ideogram.json"
        cfg.write_text(json.dumps({"models": models, "plot": plot_cfg or {}, "trait": trait or ""}), encoding="utf-8")
        p = subprocess.run([str(rscript), str(BACKEND_DIR / "qtl_compare_ideogram.R"), str(cfg), str(out)],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding="utf-8", errors="replace",
                           cwd=str(BACKEND_DIR), env=get_r_env(rscript), timeout=600, **no_window_kwargs())
    if p.returncode != 0 or not out.exists():
        print("[QTL] compare ideogram (R) failed:", (p.stderr or p.stdout)[-600:], flush=True)
        return None
    return str(out.resolve()).replace("\\", "/")


def plot_comparison(model_results, out_dir, plot_cfg=None, trait=""):
    """Per-model QTL regions (bar = left..right flanking marker, dot = peak)
    and the number of QTL each model found."""
    scans = load_scans(model_results)
    if not scans:
        return None, None
    out_dir = Path(out_dir)
    ideogram = _r_ideogram(model_results, out_dir, plot_cfg, trait)
    if ideogram:
        return ideogram, _counts_bar(model_results, scans, out_dir)
    order, lengths = _chromosome_lengths(scans)
    models = list(scans)
    n = len(models)
    regions = _regions_by_model(model_results)

    fig = Figure(figsize=(max(7, 1.3 * len(order) + 2), 6.5))
    ax = fig.add_subplot(111)
    _draw_chromosomes(ax, order, lengths)
    lane = 0.5 / n
    for k, model in enumerate(models):
        sub = regions[regions["Model"] == model]
        x0 = -0.25 + lane * (k + 0.5)
        col = MODEL_COLORS.get(str(model).upper(), INK)
        for i, ch in enumerate(order):
            for _, r in sub[sub["Chromosome"] == ch].iterrows():
                ax.vlines(i + x0, r["Left_Position"], r["Right_Position"], color=col, linewidth=5, alpha=0.45, zorder=3)
                ax.plot([i + x0], [r["Peak_Position"]], marker="o", markersize=5, color=col, zorder=4)
    ax.legend(handles=[Line2D([0], [0], color=MODEL_COLORS.get(m, INK), lw=3, label=m) for m in models],
              loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False, title="Model", fontsize=9)
    ax.set_title("QTL regions by model (bar = flanking markers, dot = peak)", color=INK, fontsize=12, fontweight="bold")
    _style(ax)
    ideogram = _save(fig, out_dir / "compare_ideogram.png")
    return ideogram, _counts_bar(model_results, scans, out_dir)


def _counts_bar(model_results, scans, out_dir):
    """Bar chart: number of QTL regions found by each model."""
    models = list(scans)
    n = len(models)
    regions = _regions_by_model(model_results)
    counts = [int((regions["Model"] == m).sum()) for m in models]
    fig = Figure(figsize=(max(5.5, 1.1 * n + 2), 4.2))
    ax = fig.add_subplot(111)
    bars = ax.bar(models, counts, color=[MODEL_COLORS.get(m, INK) for m in models], width=0.6, zorder=2)
    for b, c in zip(bars, counts):
        ax.annotate(str(c), (b.get_x() + b.get_width() / 2, b.get_height()), ha="center", va="bottom",
                    xytext=(0, 3), textcoords="offset points", color=INK, fontsize=10, fontweight="bold")
    ax.set_ylabel("QTL regions", color=MUTED)
    ax.set_title("QTL regions detected per model", color=INK, fontsize=12, fontweight="bold")
    ax.grid(axis="y", color="#e6efe9", zorder=0)
    ax.set_ylim(0, max(counts + [1]) * 1.18)
    _style(ax)
    return _save(fig, out_dir / "compare_counts.png")


# ───────────────────────── Common QTL ─────────────────────────

def plot_common(common_rows, model_results, min_models, out_dir):
    scans = load_scans(model_results)
    if not scans:
        return None, None
    out_dir = Path(out_dir)
    order, lengths = _chromosome_lengths(scans)
    common = pd.DataFrame(common_rows)
    max_models = len(scans)
    cmap = {k: c for k, c in zip(range(2, 7), ["#95d5b2", "#52b788", "#2d6a4f", "#1b4332", "#081c15"])}

    fig = Figure(figsize=(max(7, 1.3 * len(order) + 2), 6.5))
    ax = fig.add_subplot(111)
    _draw_chromosomes(ax, order, lengths)
    if not common.empty:
        common["Chromosome"] = common["Chromosome"].astype(str)
        for _, r in common.iterrows():
            if r["Chromosome"] not in order:
                continue
            i = order.index(r["Chromosome"])
            col = cmap.get(int(r["N_Models_Detected"]), INK)
            ax.vlines(i, r["Left_Position"], r["Right_Position"], color=col, linewidth=9, alpha=0.55, zorder=3)
            ax.hlines(r["Peak_Position"], i - 0.26, i + 0.26, color=col, linewidth=2.2, zorder=4)
    levels = range(max(2, min_models), max_models + 1)
    ax.legend(handles=[Line2D([0], [0], color=cmap.get(k, INK), lw=6, label=f"{k} models") for k in levels],
              loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False, title="Detected by", fontsize=9)
    ax.set_title(f"Common QTL regions (found by ≥ {min_models} models)", color=INK, fontsize=12, fontweight="bold")
    _style(ax)
    ideogram = _save(fig, out_dir / "common_ideogram.png")

    fig = Figure(figsize=(max(5.5, 0.9 * len(order) + 2), 4.2))
    ax = fig.add_subplot(111)
    per_chr = common.groupby("Chromosome").size() if not common.empty else pd.Series(dtype=int)
    counts = [int(per_chr.get(c, 0)) for c in order]
    bars = ax.bar([f"Chr {c}" for c in order], counts, color="#52b788", edgecolor=CHR_EDGE, width=0.6, zorder=2)
    for b, c in zip(bars, counts):
        if c:
            ax.annotate(str(c), (b.get_x() + b.get_width() / 2, b.get_height()), ha="center", va="bottom",
                        xytext=(0, 3), textcoords="offset points", color=INK, fontsize=10, fontweight="bold")
    ax.set_ylabel("Common QTL regions", color=MUTED)
    ax.set_title("Common QTL regions per chromosome", color=INK, fontsize=12, fontweight="bold")
    ax.grid(axis="y", color="#e6efe9", zorder=0)
    ax.set_ylim(0, max(counts + [1]) * 1.18)
    _style(ax)
    bar = _save(fig, out_dir / "common_counts.png")
    return ideogram, bar