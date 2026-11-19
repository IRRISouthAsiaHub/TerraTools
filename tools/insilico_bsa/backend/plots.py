"""Figures of the In silico BSA tool (matplotlib, PNG)."""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402
import numpy as np                                   # noqa: E402

C_LOW, C_HIGH, C_THR = "#2B7BBA", "#D7263D", "#222222"
C_GENO = {0: "#FFD400", 1: "#4393C3", 2: "#FF0052"}
CHR_PAIR = ("#2A9D8F", "#264653")
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.spines.top": False, "axes.spines.right": False, "font.size": 11})


def _save(fig, out_dir, name, dpi=150):
    p = os.path.join(out_dir, name)
    fig.savefig(p, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


def _genome_x(mk):
    """Cumulative x position (Mb) of every marker, chromosome offsets and the centre of each chromosome."""
    x = np.zeros(len(mk))
    centres, offsets, off = {}, {}, 0.0
    for c in sorted(mk["Chromosome"].unique()):
        sel = (mk["Chromosome"] == c).to_numpy()
        pos = mk["Position"].to_numpy()[sel] / 1e6
        x[sel] = pos - pos.min() + off
        offsets[c] = off
        centres[c] = off + (pos.max() - pos.min()) / 2
        off += pos.max() - pos.min() + max(1.0, 0.02 * (mk["Position"].max() / 1e6))
    return x, centres, offsets, off


def _win_x(mk, tab):
    """x (Mb, cumulative over the chromosomes as in _genome_x) of the centre of every window of `tab`."""
    _, centres, offsets, _ = _genome_x(mk)
    mn = mk.groupby("Chromosome")["Position"].min()
    cen = (tab["Window_start"].to_numpy() + tab["Window_end"].to_numpy()) / 2.0
    return np.array([offsets[c] + (m - mn[c]) / 1e6 for c, m in zip(tab["Chromosome"], cen)]), centres


def _win_panel(ax, mk, tab, col, thr, ylabel, size_kb, title=None, regions=None, sig_col="Sig_Gw"):
    x, centres = _win_x(mk, tab)
    v = tab[col].to_numpy()
    sig = tab[sig_col].to_numpy().astype(bool)
    width = size_kb / 1000.0 * 0.92
    for c in sorted(tab["Chromosome"].unique()):
        sel = (tab["Chromosome"] == c).to_numpy() & ~sig & np.isfinite(v)
        ax.bar(x[sel], v[sel], width=width, color=CHR_PAIR[int(c) % 2], linewidth=0)
    sel = sig & np.isfinite(v)
    ax.bar(x[sel], v[sel], width=width, color=C_HIGH, linewidth=0, zorder=3)
    if np.isfinite(thr):
        ax.axhline(thr, color=C_THR, ls="--", lw=1.1)
        ax.text(x.max(), thr, f" threshold {thr:.3g}", va="bottom", ha="right", fontsize=9, color=C_THR)
    if regions is not None and len(regions):
        for r in regions.itertuples():
            i = np.flatnonzero((tab["Chromosome"].to_numpy() == r.Chromosome) & (tab["Window_start"].to_numpy() == r.Peak_window_start))
            if len(i) and np.isfinite(v[i[0]]):
                ax.annotate(r.Region.split("-")[-1], (x[i[0]], v[i[0]]), textcoords="offset points", xytext=(0, 3), ha="center", fontsize=8, color=C_HIGH, fontweight="bold")
    ax.set_xticks(list(centres.values()))
    ax.set_xticklabels([f"{c}" for c in centres])
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=0)
    if title:
        ax.set_title(title, loc="left", fontweight="bold", fontsize=11)


def fig_gw(a, out):
    """The same genome in the three window sizes, one above the other: where the sizes agree and where they differ."""
    mk, wins = a["markers"], a["windows"]
    fig, axs = plt.subplots(len(wins), 1, figsize=(13, 3.6 * len(wins) + 0.6), sharex=True)
    for ax, (lab, w) in zip(np.atleast_1d(axs), wins.items()):
        _win_panel(ax, mk, w["table"], "Gw", w["thr_use"]["Gw"], "Window Gw", w["kb"], regions=w["regions"],
                   title=f"{lab} windows — {len(w['regions'])} region{'s' if len(w['regions']) != 1 else ''}")
    np.atleast_1d(axs)[-1].set_xlabel("Chromosome")
    fig.suptitle(f"Window Gw along the genome — {a['trait']}, {a['pool_size']} lines per pool (threshold: permutations, alpha {a['alpha']:g})", x=0.01, ha="left", fontsize=12)
    fig.tight_layout()
    return _save(fig, out, "gw_windows_compared.png")


def fig_methods(a, out):
    mk, w = a["markers"], a["windows"][a["primary"]]
    fig, axs = plt.subplots(3, 1, figsize=(13, 9), sharex=True)
    items = (("Gw", "Gw", "Sig_Gw", "Gw — alleles and genotypes"), ("G1", "G1", "Sig_G1", "G1 — allele counts only (as G′ / ED / ΔSNP-index)"),
             ("abs_dAF", "dAF", "Sig_dAF", "|Δ allele frequency| (the Δ SNP-index)"))
    for ax, (col, key, sig, title) in zip(axs, items):
        _win_panel(ax, mk, w["table"], col, w["thr_use"][key], f"Window value ({a['primary']})", w["kb"], title=title,
                   regions=w["regions"] if key == "Gw" else None, sig_col=sig)
    axs[-1].set_xlabel("Chromosome")
    fig.tight_layout()
    return _save(fig, out, "methods_compared.png")


def fig_loci(a, out):
    """Every locus (rows) and the interval each window size gives for it (bars): resolution against sensitivity."""
    loci, wins = a["loci"], a["windows"]
    if not len(loci):
        raise ValueError("no loci")
    labs = list(wins)
    top = loci.head(24)
    cols = {lab: c for lab, c in zip(labs, ("#2A9D8F", "#E9C46A", "#E76F51", "#6A4C93", "#264653", "#8AB17D"))}
    fig, ax = plt.subplots(figsize=(11, 0.55 * len(top) * len(labs) / 2 + 2.4))
    reg = a["regions"].set_index("Region")
    y = 0
    ticks, names = [], []
    for r in top.itertuples():
        base = y
        for j, lab in enumerate(labs):
            rid = getattr(r, f"Region_{lab}")
            if rid:
                g = reg.loc[rid]
                ax.barh(y, (g.End_bp - g.Start_bp) / 1e6, left=g.Start_bp / 1e6, height=0.7, color=cols[lab], label=lab if r.Index == 0 else None)
                ax.plot(g.Peak_bp / 1e6, y, "|", color="black", ms=9)
            y += 1
        ticks.append((base + y - 1) / 2)
        names.append(f"{r.Locus}  chr{r.Chromosome}")
        y += 0.6
    ax.set_yticks(ticks)
    ax.set_yticklabels(names, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Position along the chromosome (Mb); | marks the peak marker")
    h, l = ax.get_legend_handles_labels()
    seen = dict(zip(l, h))
    ax.legend(seen.values(), seen.keys(), frameon=False, title="Window", ncol=len(seen), loc="lower center", bbox_to_anchor=(0.5, 1.0))
    ax.set_title("Loci found by the three window sizes", loc="left", fontweight="bold", pad=34)
    return _save(fig, out, "loci_compared.png")


def fig_pools(a, out):
    fig, ax = plt.subplots(figsize=(8.5, 4.4))
    y = a["y"]
    lo, hi = a["pool_cuts"]
    bins = np.histogram_bin_edges(y, bins=min(40, max(10, len(y) // 8)))
    ax.hist(y, bins=bins, color="#BBBBBB", edgecolor="white")
    ax.hist(y[a["lo_i"]], bins=bins, color=C_LOW, edgecolor="white", label=f"Low pool ({a['pool_size']} lines)")
    ax.hist(y[a["hi_i"]], bins=bins, color=C_HIGH, edgecolor="white", label=f"High pool ({a['pool_size']} lines)")
    ax.set(xlabel=a["trait"], ylabel="Number of lines", title="Trait distribution and the two in silico pools")
    ax.legend(frameon=False)
    return _save(fig, out, "pools.png")


def fig_genotypes(a, out):
    reg = a["windows"][a["primary"]]["regions"].head(6)
    if not len(reg):
        raise ValueError("no regions")
    mk = a["markers"].set_index("Marker")
    n = len(reg)
    fig, axs = plt.subplots(1, n, figsize=(2.6 * n + 1, 4.2), sharey=True)
    axs = np.atleast_1d(axs)
    for ax, r in zip(axs, reg.itertuples()):
        row = mk.loc[r.Peak_marker]
        for j, pool in enumerate(("Low", "High")):
            c = np.array([row[f"{pool}_0"], row[f"{pool}_1"], row[f"{pool}_2"]], dtype=float)
            f = c / max(c.sum(), 1)
            bottom = 0.0
            for g in (0, 1, 2):
                if f[g] > 0:
                    ax.bar(j, f[g], bottom=bottom, color=C_GENO[g], edgecolor="white", width=0.7)
                    if f[g] > 0.07:
                        ax.text(j, bottom + f[g] / 2, f"{f[g] * 100:.0f}%", ha="center", va="center", fontsize=9)
                bottom += f[g]
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Low", "High"])
        ax.set_title(f"{r.Region}  chr{r.Chromosome}\n{r.Peak_bp / 1e6:.2f} Mb", fontsize=10)
    axs[0].set_ylabel("Genotype frequency in the pool")
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color=C_GENO[g], label=l) for g, l in ((0, "0  parent 1"), (1, "1  heterozygous"), (2, "2  parent 2"))],
               loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.06))
    fig.suptitle("Genotype classes in the low and high pool at the peak marker of each region", y=1.02)
    return _save(fig, out, "genotypes_at_peaks.png")


def fig_null(a, out):
    if a["n_perm"] < 20:
        raise ValueError("no permutations")
    wins = a["windows"]
    fig, axs = plt.subplots(1, len(wins), figsize=(4.8 * len(wins), 4.2), sharey=True)
    for ax, (lab, w) in zip(np.atleast_1d(axs), wins.items()):
        mx = w["maxima"]["Gw"]
        mx = mx[np.isfinite(mx)]
        ax.hist(mx, bins=40, color="#9AA5B1", edgecolor="white")
        for al, ls in ((0.10, ":"), (0.05, "--"), (0.01, "-")):
            v = w["thr"]["Gw"].get(al)
            if v is not None and np.isfinite(v):
                ax.axvline(v, color=C_THR, ls=ls, lw=1.2, label=f"alpha {al:g}: {v:.3g}")
        obs = np.nanmax(w["table"]["Gw"].to_numpy()) if np.isfinite(w["table"]["Gw"].to_numpy()).any() else np.nan
        if np.isfinite(obs):
            ax.axvline(obs, color=C_HIGH, lw=2, label=f"observed maximum {obs:.3g}")
        ax.set(xlabel="Largest window Gw of the genome", title=f"{lab} windows")
        ax.legend(frameon=False, fontsize=8)
    np.atleast_1d(axs)[0].set_ylabel("Permutations")
    fig.suptitle(f"Null distribution of the largest window Gw ({a['n_perm']:,} random low / high pools)", x=0.01, ha="left")
    fig.tight_layout()
    return _save(fig, out, "permutation_null.png")


def make_all(a: dict, out_dir: str) -> list:
    """-> [{'key', 'label', 'path'}] of the figures that apply to this run."""
    items = [("gw", "Window Gw along the genome, 250 kb / 500 kb / 1 Mb compared", fig_gw), ("loci", "Loci found by each window size", fig_loci),
             ("methods", "Gw compared with the allele-frequency statistics", fig_methods),
             ("genotypes", "Genotype classes in the two pools at each peak", fig_genotypes), ("pools", "Trait distribution and the two pools", fig_pools),
             ("null", "Permutation null distributions", fig_null)]
    out = []
    for key, label, fn in items:
        try:
            out.append({"key": key, "label": label, "path": fn(a, out_dir)})
        except ValueError:                                        # nothing to draw (no regions, no permutations): skipped on purpose
            continue
        except Exception:                                         # noqa: BLE001  one failed figure must not lose the analysis
            import traceback
            traceback.print_exc()
    return out
