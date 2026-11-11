"""Figures (matplotlib, saved as PNG and PDF) and the PDF report."""

from __future__ import annotations

import datetime as _dt
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402
import numpy as np                                   # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.patches import Patch                  # noqa: E402

from tools.qtl_seq.backend import stats              # noqa: E402

LEVEL_COLOR = {90: "#f4a261", 95: "#2a6f97", 99: "#c1121f"}
CHROM_COLORS = ("#7a9e7e", "#3d5a40")
DPI = 300
MAX_POINTS = 250_000


def safe_name(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(s))


def save_figure(fig, stem):
    """Writes <stem>.png (300 dpi) and <stem>.pdf. Returns the png path."""
    stem = Path(stem)
    fig.savefig(stem.with_suffix(".png"), dpi=DPI, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    return str(stem.with_suffix(".png"))


def layout(chrom_lengths, chroms):
    """Offsets that put the chromosomes one after the other on a single axis. -> ({chrom: offset}, total)"""
    off, x = {}, 0
    gap = 0.01 * sum(chrom_lengths[c] for c in chroms)
    for c in chroms:
        off[c] = x
        x += chrom_lengths[c] + gap
    return off, x - gap


def _decimate(df, n=MAX_POINTS):
    return df if len(df) <= n else df.sample(n, random_state=1).sort_values(["Chromosome", "Position"])


def _axis_chroms(ax, off, chrom_lengths, chroms, total):
    ticks, labels = [], []
    for i, c in enumerate(chroms):
        ticks.append(off[c] + chrom_lengths[c] / 2)
        labels.append(re.sub(r"^(chromosome|chrom|chr)[_ ]?", "", str(c), flags=re.I) or str(c))
        if i:
            ax.axvline(off[c] - 0.005 * total, color="#cccccc", lw=0.6, zorder=0)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_xlim(0, total)
    ax.set_xlabel("Chromosome")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def fig_bulk_index(snps, chrom_lengths, which, bulk_name):
    """SNP-index of one bulk (which = 'High' or 'Low') along the genome."""
    chroms = [c for c in chrom_lengths if c in set(snps["Chromosome"])]
    off, total = layout(chrom_lengths, chroms)
    fig, ax = plt.subplots(figsize=(13, 3.6))
    d = _decimate(snps)
    for i, c in enumerate(chroms):
        s = d[d["Chromosome"] == c]
        ax.scatter(s["Position"] + off[c], s[f"SNP_index_{which}"], s=2.2, color=CHROM_COLORS[i % 2], alpha=0.55,
                   linewidths=0, rasterized=True)
    ax.set_ylim(-0.02, 1.02)
    ax.set_ylabel("SNP-index")
    ax.set_title(f"SNP-index of the {bulk_name}", fontsize=11)
    _axis_chroms(ax, off, chrom_lengths, chroms, total)
    fig.tight_layout()
    return fig


def _draw_delta(ax, snps, win, qtls, off, chroms, levels=stats.LEVELS, points=True):
    if points:
        d = _decimate(snps)
        for i, c in enumerate(chroms):
            s = d[d["Chromosome"] == c]
            ax.scatter(s["Position"] + off[c], s["Delta_SNP_index"], s=2.2, color=CHROM_COLORS[i % 2], alpha=0.35,
                       linewidths=0, rasterized=True)
    ax.axhline(0, color="#555555", lw=0.6)
    for c in chroms:
        w = win[win["Chromosome"] == c]
        if w.empty:
            continue
        x = w["Window_Mid"] + off[c]
        for lv in levels:
            for side in ("Upper", "Lower"):
                ax.plot(x, w[f"CI{lv}_{side}"], color=LEVEL_COLOR[lv], lw=1.0, alpha=0.95, zorder=3)
        ax.plot(x, w["Mean_Delta_SNP_index"], color="#1b1b1b", lw=1.5, zorder=4)
        # the stretches of the curve that are outside the limits
        for lv in levels:
            sig = w[f"Significant_{lv}"].to_numpy()
            if sig.any():
                y = np.where(sig, w["Mean_Delta_SNP_index"].to_numpy(), np.nan)
                ax.plot(x, y, color=LEVEL_COLOR[lv], lw=3.0, zorder=5, solid_capstyle="butt")
    if qtls is not None and len(qtls):
        for q in qtls.itertuples(index=False):
            if q.Chromosome in off:
                ax.axvspan(q.Start + off[q.Chromosome], q.End + off[q.Chromosome], color="#52b788", alpha=0.18, zorder=1)
                ax.plot([q.Peak_Position + off[q.Chromosome]], [q.Peak_Delta_SNP_index], marker="v" if q.Peak_Delta_SNP_index > 0 else "^",
                        color="#1b1b1b", ms=6, zorder=6)
    ax.set_ylim(-1.05, 1.05)
    ax.set_ylabel("Δ SNP-index")


def _legend(ax, levels=stats.LEVELS):
    h = [plt.Line2D([0], [0], color="#1b1b1b", lw=1.5, label="Sliding-window mean")]
    h += [plt.Line2D([0], [0], color=LEVEL_COLOR[lv], lw=1.2, label=f"{lv}% confidence limit") for lv in levels]
    h += [Patch(color="#52b788", alpha=0.3, label="QTL interval")]
    ax.legend(handles=h, loc="lower left", ncol=5, fontsize=8, frameon=False, bbox_to_anchor=(0, 1.0))


def fig_genomewide(snps, win, qtls, chrom_lengths, high_name="High", low_name="Low"):
    chroms = [c for c in chrom_lengths if c in set(snps["Chromosome"])]
    off, total = layout(chrom_lengths, chroms)
    fig, ax = plt.subplots(figsize=(13, 5.2))
    _draw_delta(ax, snps, win, qtls, off, chroms)
    _axis_chroms(ax, off, chrom_lengths, chroms, total)
    ax.set_title(f"Δ SNP-index ({high_name} − {low_name})", fontsize=12, pad=26)
    _legend(ax)
    fig.tight_layout()
    return fig


def fig_chromosome(chrom, snps, win, qtls, chrom_lengths, high_name="High", low_name="Low"):
    """High-bulk index, low-bulk index and Δ SNP-index of one chromosome, with limits, QTL intervals and peaks."""
    s = snps[snps["Chromosome"] == chrom]
    w = win[win["Chromosome"] == chrom]
    q = qtls[qtls["Chromosome"] == chrom] if qtls is not None and len(qtls) else qtls
    L = chrom_lengths[chrom]
    fig, axes = plt.subplots(3, 1, figsize=(11, 8.2), sharex=True, gridspec_kw={"height_ratios": [1, 1, 1.6]})
    for ax, key, name in ((axes[0], "SNP_index_High", high_name), (axes[1], "SNP_index_Low", low_name)):
        ax.scatter(s["Position"] / 1e6, s[key], s=3, color=CHROM_COLORS[0], alpha=0.55, linewidths=0, rasterized=True)
        if len(w):
            ax.plot(w["Window_Mid"] / 1e6, w["Mean_SNP_index_High" if key.endswith("High") else "Mean_SNP_index_Low"],
                    color="#1b1b1b", lw=1.4)
        ax.set_ylim(-0.02, 1.02)
        ax.set_ylabel(f"SNP-index\n{name}", fontsize=9)
    ax = axes[2]
    ax.scatter(s["Position"] / 1e6, s["Delta_SNP_index"], s=3, color=CHROM_COLORS[0], alpha=0.4, linewidths=0, rasterized=True)
    ax.axhline(0, color="#555555", lw=0.6)
    if len(w):
        x = w["Window_Mid"] / 1e6
        for lv in stats.LEVELS:
            for side in ("Upper", "Lower"):
                ax.plot(x, w[f"CI{lv}_{side}"], color=LEVEL_COLOR[lv], lw=1.0)
            sig = w[f"Significant_{lv}"].to_numpy()
            if sig.any():
                ax.plot(x, np.where(sig, w["Mean_Delta_SNP_index"], np.nan), color=LEVEL_COLOR[lv], lw=3)
        ax.plot(x, w["Mean_Delta_SNP_index"], color="#1b1b1b", lw=1.5, zorder=4)
    if q is not None and len(q):
        for r in q.itertuples(index=False):
            for a in axes:
                a.axvspan(r.Start / 1e6, r.End / 1e6, color="#52b788", alpha=0.18, zorder=0)
            ax.plot([r.Peak_Position / 1e6], [r.Peak_Delta_SNP_index], marker="v" if r.Peak_Delta_SNP_index > 0 else "^",
                    color="#1b1b1b", ms=7, zorder=6)
            ax.annotate(f"{r.QTL_ID}  peak {r.Peak_Position / 1e6:.2f} Mb", (r.Peak_Position / 1e6, r.Peak_Delta_SNP_index),
                        textcoords="offset points", xytext=(6, 6 if r.Peak_Delta_SNP_index > 0 else -14), fontsize=8)
    ax.set_ylim(-1.05, 1.05)
    ax.set_ylabel("Δ SNP-index", fontsize=9)
    ax.set_xlim(0, L / 1e6)
    ax.set_xlabel(f"Position on {chrom} (Mb)")
    for a in axes:
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
    axes[0].set_title(f"{chrom}", fontsize=11)
    fig.tight_layout()
    return fig


def fig_simulation(limits, ci_method, depth_example=None, seed=1):
    """Left: simulated null distribution of Δ SNP-index for one depth with the 90/95/99 % limits. Right: how the
    limits change with read depth."""
    meta = limits.meta
    fr = limits.frame
    if depth_example is None:
        # the most common depth pair of the table: the one closest to the median depth
        row = stats.typical_pair(fr)
    else:
        row = fr.iloc[0]
    dh, dl = int(row["Depth_High"]), int(row["Depth_Low"])
    rng = np.random.default_rng(seed)
    sim = stats.simulate_null(dh, dl, meta["population"], meta["n_high"], meta["n_low"], int(meta["replicates"]), rng,
                              meta.get("bc_recurrent", "non_index"))
    fig, (a, b) = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"width_ratios": [1.15, 1]})
    span = max(abs(float(row["CI99_Lower"])), abs(float(row["CI99_Upper"]))) * 1.35
    a.hist(sim, bins=32, range=(-span, span), color="#b7d9c4", edgecolor="white")
    a.set_xlim(-span, span)
    for lv in stats.LEVELS:
        for side in ("Lower", "Upper"):
            v = row[f"CI{lv}_{side}"]
            a.axvline(v, color=LEVEL_COLOR[lv], lw=1.6, label=f"{lv}%" if side == "Lower" else None)
    a.set_xlabel("Simulated Δ SNP-index under the null (no QTL)")
    a.set_ylabel("Replicates")
    a.set_title(f"Null distribution: {meta['population']}, {meta['n_high']} + {meta['n_low']} plants, "
                f"depth {dh} / {dl}, {int(meta['replicates']):,} replicates", fontsize=10)
    a.legend(title="Confidence limit", frameon=False, fontsize=8)
    eq = fr[fr["Depth_High"] == fr["Depth_Low"]].sort_values("Depth_High")
    use = eq if len(eq) >= 2 else fr.assign(_m=np.minimum(fr["Depth_High"], fr["Depth_Low"])).sort_values("_m")
    xs = use["Depth_High"] if use is eq else use["_m"]
    for lv in stats.LEVELS:
        b.plot(xs, use[f"CI{lv}_Upper"], color=LEVEL_COLOR[lv], lw=1.6, label=f"{lv}%")
        b.plot(xs, use[f"CI{lv}_Lower"], color=LEVEL_COLOR[lv], lw=1.6)
    b.axhline(0, color="#555555", lw=0.6)
    b.set_xlabel("Read depth" + (" (both bulks)" if use is eq else " (smaller of the two bulks)"))
    b.set_ylabel("Confidence limit of Δ SNP-index")
    b.set_title("Limits by read depth", fontsize=10)
    b.set_ylim(-1.02, 1.02)
    for ax in (a, b):
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    note = ("Window limit = mean of the per-SNP limits in the window" if ci_method == "snp_average"
            else "Window limit = the whole window simulated (number of SNPs and depths of the window)")
    fig.text(0.5, -0.02, note, ha="center", fontsize=8, color="#555555")
    fig.tight_layout()
    return fig


# ============================================================
# REPORT
# ============================================================

def _text_page(lines, title, size=9):
    fig = plt.figure(figsize=(8.27, 11.69))
    fig.text(0.07, 0.96, title, fontsize=15, weight="bold", va="top")
    y = 0.925
    for kind, text in lines:
        if kind == "h":
            y -= 0.012
            fig.text(0.07, y, text, fontsize=11, weight="bold", va="top")
            y -= 0.027
        else:
            for chunk in _wrap(text, 105):
                fig.text(0.07, y, chunk, fontsize=size, va="top", family="DejaVu Sans")
                y -= 0.0175
        if y < 0.05:
            break
    return fig


def _wrap(text, width):
    out, cur = [], ""
    for word in str(text).split(" "):
        if len(cur) + len(word) + 1 > width:
            out.append(cur)
            cur = word
        else:
            cur = (cur + " " + word).strip()
    return out + [cur]


def write_report(path, ctx):
    """QTLseq_analysis_report.pdf. ctx: params, ledger, warnings, snps, windows, qtls, genes, chrom_lengths, limits."""
    p, led = ctx["params"], ctx["ledger"]
    qtls, snps, win = ctx["qtls"], ctx["snps"], ctx["windows"]
    counts = stats.level_counts(qtls)
    with PdfPages(path) as pdf:
        lines = [("h", "Analysis"),
                 ("p", f"Date: {_dt.datetime.now().strftime('%Y-%m-%d %H:%M')}   Input: {p['mode'].upper()}   "
                       f"Marker source: {ctx['marker_mode_text']}"),
                 ("p", f"Population: {p['population']}   {p['high_name']}: {p['n_high']} plants   {p['low_name']}: {p['n_low']} plants"),
                 ("p", f"SNP-index counts the allele of: {ctx['index_text']}"),
                 ("h", "SNP filtering"),
                 ("p", f"SNPs identified: {led.get('total_snps', 0):,}"),
                 ("p", f"Passing parental polymorphism filtering: {led.get('passing_parental', 0):,}  "
                       f"(removed for parental criteria: {led.get('removed_parental', 0):,})"),
                 ("p", f"Removed for quality: {led.get('removed_quality', 0):,}   Removed for low depth (< {p['min_depth']}): "
                       f"{led.get('removed_low_depth', 0):,}   for high depth: {led.get('removed_high_depth', 0):,}   "
                       f"for SNP-index < {p['min_snp_index']} in both bulks: {led.get('removed_snp_index', 0):,}"),
                 ("p", f"SNPs retained for the analysis: {led.get('retained', 0):,}"),
                 ("h", "Sliding window and confidence limits"),
                 ("p", f"Window {p['window_size']:,} bp, step {p['window_step']:,} bp, at least {p['min_snps_per_window']} SNPs per window. "
                       f"{ctx['n_windows']:,} windows."),
                 ("p", f"Simulation: {p['replicates']:,} replicates, seed {p['seed']}, limit method '{p['ci_method']}', "
                       f"depth {'observed per SNP' if p['depth_mode'] == 'auto' else 'fixed at ' + str(p['fixed_depth'])}."),
                 ("p", "Limits are two-sided 90 / 95 / 99 % quantiles of the simulated null distribution. There is no correction "
                       "for the many overlapping windows tested."),
                 ("h", "QTL"),
                 ("p", f"QTL intervals: {len(qtls)}   reaching 99%: {counts['99']}   95%: {counts['95']}   90%: {counts['90']}"),
                 ]
        if ctx["warnings"]:
            lines.append(("h", "Warnings"))
            lines += [("p", "- " + w) for w in ctx["warnings"][:14]]
        pdf.savefig(_text_page(lines, "QTL-seq analysis report"))
        for fig in (fig_genomewide(snps, win, qtls, ctx["chrom_lengths"], p["high_name"], p["low_name"]),
                    fig_bulk_index(snps, ctx["chrom_lengths"], "High", p["high_name"]),
                    fig_bulk_index(snps, ctx["chrom_lengths"], "Low", p["low_name"])):
            pdf.savefig(fig)
            plt.close(fig)
        lines = [("h", "Significant QTL intervals")]
        if len(qtls):
            for q in qtls.itertuples(index=False):
                lines.append(("p", f"{q.QTL_ID}  {q.Chromosome}:{q.Start:,}-{q.End:,}  peak {q.Peak_Position:,}  "
                                   f"Δ {q.Peak_Delta_SNP_index:+.3f}  {q.Significance_Level}  SNPs {q.SNP_Count}"))
        else:
            lines.append(("p", "No window crossed the selected confidence limit."))
        genes = ctx.get("genes")
        if genes is not None and len(genes):
            lines.append(("h", "Genes overlapping the intervals (position only; not evidence of causality)"))
            for g in genes.head(40).itertuples(index=False):
                lines.append(("p", f"{g.QTL_ID}  {g.Gene_ID}  {g.Gene_Start:,}-{g.Gene_End:,}  {g.Strand}  {str(g.Annotation)[:70]}"))
        pdf.savefig(_text_page(lines, "Results"))
        fig = fig_simulation(ctx["limits"], p["ci_method"]) if ctx.get("limits") is not None else None
        if fig is not None:
            pdf.savefig(fig)
            plt.close(fig)
        for chrom in list(dict.fromkeys(qtls["Chromosome"]))[:12] if len(qtls) else []:
            fig = fig_chromosome(chrom, snps, win, qtls, ctx["chrom_lengths"], p["high_name"], p["low_name"])
            pdf.savefig(fig)
            plt.close(fig)
