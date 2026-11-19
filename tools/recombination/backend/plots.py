"""Figures of the Recombination Analysis (matplotlib, PNG)."""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402
import numpy as np                                   # noqa: E402
from matplotlib.colors import ListedColormap         # noqa: E402
from matplotlib.patches import Patch                 # noqa: E402

from tools.recombination.backend.recombination import RICE_CENTROMERES  # noqa: E402

C_A, C_B, C_H, C_U = "#FFD400", "#FF0052", "#4393C3", "#E8E8E8"     # recurrent parent, donor, het, missing
CHR_COLORS = ["#E63946", "#F4A261", "#E9C46A", "#2A9D8F", "#457B9D", "#9B5DE5", "#F15BB5", "#00BBF9", "#00F5A0", "#FF6B6B",
              "#A8DADC", "#6A4C93", "#8D99AE", "#BC6C25", "#606C38", "#8338EC"]
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.spines.top": False, "axes.spines.right": False, "font.size": 11})


def _save(fig, out_dir, name, dpi=150):
    p = os.path.join(out_dir, name)
    fig.savefig(p, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


def _chr_color(c):
    return CHR_COLORS[(int(c) - 1) % len(CHR_COLORS)]


def fig_crossovers(a, out):
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1.1, 1]})
    per = a["per_ril"]["Crossovers"].to_numpy()
    top = max(int(per.max()), 1)
    ax[0].hist(per, bins=np.arange(-0.5, top + 1.5, 1), color="#4DAF4A", edgecolor="white")
    m = per.mean()
    ax[0].axvline(m, color="red", ls="--", lw=1.2)
    ax[0].text(m, ax[0].get_ylim()[1] * 0.96, f" mean {m:.1f}", color="red", va="top", fontweight="bold")
    ax[0].set(xlabel="Crossovers per line", ylabel="Number of lines", title="Crossovers per line")
    pc = a["per_chr"]
    ax[1].bar([f"Chr{c}" for c in pc["Chromosome"]], pc["Mean_COs_per_RIL"], color=[_chr_color(c) for c in pc["Chromosome"]])
    ax[1].set(ylabel="Mean crossovers per line", title="Mean crossovers per line and chromosome")
    ax[1].tick_params(axis="x", rotation=90)
    return _save(fig, out, "crossovers_per_line.png")


def fig_rate(a, out):
    chroms, L, bins = a["chroms"], a["chr_len_mb"], a["bins"]
    gap = 2.0
    off, offs = 0.0, {}
    for c in chroms:
        offs[c] = off
        off += L[c] + gap
    fig, ax = plt.subplots(figsize=(15, 4.6))
    ymax = max(float(bins["COs_per_RIL"].max()) * 1.12, 0.05)
    for k, c in enumerate(chroms):
        o = offs[c]
        ax.axvspan(o, o + L[c], color="#FFF8DC" if k % 2 == 0 else "#FFFFFF", zorder=0)
        if a["rice_centromeres"] and c in RICE_CENTROMERES and RICE_CENTROMERES[c][1] <= L[c] + 1:
            s, e = RICE_CENTROMERES[c]
            ax.axvspan(o + s, o + e, color="#E8A33D", alpha=0.35, zorder=1)
        b = bins[bins["Chromosome"] == c]
        x = o + (b["Bin_start_Mb"].to_numpy() + b["Bin_end_Mb"].to_numpy()) / 2
        y = b["COs_per_RIL"].to_numpy()
        ax.fill_between(x, 0, y, color="#F1948A", alpha=0.65, zorder=2)
        ax.plot(x, y, color="#D90000", lw=1.1, zorder=3)
        ax.text(o + L[c] / 2, ymax * 1.01, f"Chr{c}", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.set_xlim(0, off - gap)
    ax.set_ylim(0, ymax)
    ax.set_xticks([])
    ax.set_ylabel(f"Crossovers per line per {a['bin_mb']:g} Mb")
    ax.set_title("Genome-wide crossover rate", pad=22)
    if a["rice_centromeres"]:
        ax.legend(handles=[Patch(color="#E8A33D", alpha=0.35, label="Centromere (rice IRGSP-1.0)")], loc="upper right", frameon=False)
    return _save(fig, out, "crossover_rate_genome.png")


def fig_donor(a, out):
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    d = a["per_ril"]["Donor_pct"].to_numpy()
    ax[0].hist(d, bins=np.arange(0, max(55, d.max() + 5), 2.5), color="#FF6B6B", edgecolor="white")
    ax[0].axvline(d.mean(), color="black", ls="--", lw=1.2)
    ax[0].text(d.mean(), ax[0].get_ylim()[1] * 0.96, f" mean {d.mean():.1f} %", va="top", fontweight="bold")
    ax[0].set(xlabel="Donor genome per line (%)", ylabel="Number of lines", title="Donor content per line")
    dc = a["donor_chr"]
    bp = ax[1].boxplot([dc[f"Chr{c}"].to_numpy() for c in a["chroms"]], patch_artist=True, showfliers=False,
                       tick_labels=[f"Chr{c}" for c in a["chroms"]], medianprops={"color": "black"})
    for patch, c in zip(bp["boxes"], a["chroms"]):
        patch.set_facecolor(_chr_color(c))
    ax[1].set(ylabel="Donor genome (%)", title="Donor content by chromosome")
    ax[1].tick_params(axis="x", rotation=90)
    return _save(fig, out, "donor_content.png")


def _order(a):
    return np.argsort(a["per_ril"]["Donor_pct"].to_numpy(), kind="stable")


def fig_mosaic(a, out):
    chroms = a["chroms"]
    pm = a["poly_markers"]
    cols = np.concatenate([np.flatnonzero(pm["Chrom"].to_numpy() == c) for c in chroms])
    order = _order(a)
    M = a["smooth"][cols][:, order].T
    nr, nm = M.shape
    cmap = ListedColormap([C_U, C_A, C_B, C_H])
    h = min(max(4.5, nr * 0.045 + 1.5), 14)
    fig, ax = plt.subplots(figsize=(15, h))
    ax.imshow(M, aspect="auto", cmap=cmap, vmin=-0.5, vmax=3.5, interpolation="nearest")
    start = 0
    ticks, labels = [], []
    for c in chroms:
        n = int((pm["Chrom"] == c).sum())
        ax.axvline(start - 0.5, color="black", lw=0.8)
        ticks.append(start + n / 2)
        labels.append(f"Chr{c}")
        start += n
    ax.axvline(nm - 0.5, color="black", lw=0.8)
    ax.set_xticks(ticks, labels)
    ax.set_yticks([])
    ax.set_ylabel(f"{nr} lines (sorted by donor %)")
    ax.set_title("Haplotype mosaic after smoothing (markers equally spaced)")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.legend(handles=[Patch(color=C_A, label="Recurrent parent (A)"), Patch(color=C_B, label="Donor (B)"),
                       Patch(color=C_H, label="Heterozygous"), Patch(color=C_U, label="Missing")],
              loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=4, frameon=False)
    return _save(fig, out, "haplotype_mosaic.png", dpi=130)


def fig_circos(a, out):
    chroms, L = a["chroms"], a["chr_len_mb"]
    pm = a["poly_markers"]
    order = _order(a)
    sm = a["smooth"][:, order]
    nr = sm.shape[1]
    gap = np.deg2rad(2.0)
    total = sum(L[c] for c in chroms)
    scale = (2 * np.pi - gap * len(chroms)) / total
    r0, r1 = 0.30, 1.0
    redges = np.linspace(r0, r1, nr + 1)
    cmap = ListedColormap([C_U, C_A, C_B, C_H])
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    th = 0.0
    for c in chroms:
        ix = np.flatnonzero(pm["Chrom"].to_numpy() == c)
        if len(ix) >= 2:
            p = pm["POS"].to_numpy()[ix] / 1e6
            mid = (p[1:] + p[:-1]) / 2
            e = np.concatenate([[0.0], mid, [L[c]]])
            ang = th + e * scale
            ax.pcolormesh(ang, redges, sm[ix].T, cmap=cmap, vmin=-0.5, vmax=3.5, shading="flat", rasterized=True)
        mid_ang = th + L[c] * scale / 2
        ax.text(mid_ang, r1 + 0.06, f"Chr{c}", ha="center", va="center", fontsize=11, fontweight="bold", color=_chr_color(c), clip_on=False)
        ax.plot(np.linspace(th, th + L[c] * scale, 60), np.full(60, r1 + 0.015), color=_chr_color(c), lw=4, solid_capstyle="butt")
        th += L[c] * scale + gap
    ax.set_ylim(0, r1 + 0.12)
    ax.axis("off")
    ax.set_title("Circos haplotype plot (innermost ring = lowest donor %)", pad=18)
    fig.legend(handles=[Patch(color=C_A, label="Recurrent parent (A)"), Patch(color=C_B, label="Donor (B)"),
                        Patch(color=C_H, label="Heterozygous"), Patch(color=C_U, label="Missing")], loc="lower center", ncol=4, frameon=False)
    return _save(fig, out, "circos_haplotype.png", dpi=170)


def fig_intra(a, out):
    co, pc = a["co"], a["per_chr"]
    ni = int((co["CO_type"] == "intragenic").sum())
    ne = int((co["CO_type"] == "intergenic").sum())
    n = max(ni + ne, 1)
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1, 1.7]})
    b = ax[0].bar(["Intragenic", "Intergenic"], [ni, ne], color=["#E91E8C", "#27AE60"])
    for r, v in zip(b, (ni, ne)):
        ax[0].text(r.get_x() + r.get_width() / 2, v, f"{v:,}\n({100 * v / n:.1f}%)", ha="center", va="bottom", fontweight="bold")
    ax[0].set(ylabel="Crossovers", title="A. Overall crossover location")
    ax[0].set_ylim(0, max(ni, ne) * 1.2)
    x = np.arange(len(pc))
    ax[1].bar(x, pc["Intragenic"], color="#E91E8C", label="Intragenic")
    ax[1].bar(x, pc["Intergenic"], bottom=pc["Intragenic"], color="#27AE60", label="Intergenic")
    ax[1].set_xticks(x, [f"Chr{c}" for c in pc["Chromosome"]], rotation=90)
    ax[1].set(ylabel="Crossovers", title="B. Crossover location per chromosome")
    ax[1].legend(frameon=False)
    return _save(fig, out, "intragenic_vs_intergenic.png")


def fig_genes(a, out):
    g = a["genes"].head(10)
    labels = [(n.split("|")[0] if isinstance(n, str) and n else l) for n, l in zip(g["Gene_Name"], g["LOC_ID"])]
    fig, ax = plt.subplots(figsize=(11, 5.2))
    cols = ["#C0392B", "#E67E22", "#F1C40F", "#27AE60", "#1ABC9C", "#2980B9", "#8E44AD", "#E91E63", "#00BCD4", "#FF5722"]
    bars = ax.bar(range(len(g)), g["Intragenic_COs"], color=cols[:len(g)])
    for r, v in zip(bars, g["Intragenic_COs"]):
        ax.text(r.get_x() + r.get_width() / 2, v, str(int(v)), ha="center", va="bottom", fontweight="bold")
    ax.set_xticks(range(len(g)), labels, rotation=45, ha="right")
    ax.set(ylabel="Intragenic crossovers (all lines)", title=f"Top {len(g)} genes with the most intragenic crossovers")
    return _save(fig, out, "top_genes.png")


def make_all(a: dict, out_dir: str) -> list:
    """-> [{'key', 'label', 'path'}] of the figures that apply to this run."""
    items = [("crossovers", "Crossovers per line and per chromosome", fig_crossovers),
             ("rate", "Genome-wide crossover rate", fig_rate),
             ("donor", "Donor content per line and per chromosome", fig_donor),
             ("mosaic", "Haplotype mosaic of all lines", fig_mosaic),
             ("circos", "Circos haplotype plot", fig_circos)]
    if a["has_genes"] and len(a["co"]):
        items.append(("intra", "Intragenic versus intergenic crossovers", fig_intra))
        if len(a["genes"]):
            items.append(("genes", "Genes with the most intragenic crossovers", fig_genes))
    out = []
    for key, label, fn in items:
        try:
            out.append({"key": key, "label": label, "path": fn(a, out_dir)})
        except Exception:                                     # noqa: BLE001  one failed figure must not lose the analysis
            import traceback
            traceback.print_exc()
    return out
