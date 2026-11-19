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


# ───────────────────────────── SNP-index figure (QTL-seq style) ─────────────────────────────

C_CHR_PAIR = ("#7CCBA2", "#F5AE7E")           # pastel mint / peach, alternating chromosomes
C_CHR_BAND = ("#F1FAF5", "#FEF6EE")           # the faint stripe behind each chromosome
C_REGION, C_CI99, C_CI95, C_LINE = "#E58C9B", "#E8738A", "#6F9BE0", "#33384D"
C_RULER = "#EDE7F6"


def _roll_sd(pos, val, centres, size):
    """Standard deviation of the marker values inside the sliding window around every centre (NaN with fewer than 3 markers)."""
    ok = np.isfinite(val)
    p, v = np.asarray(pos, dtype=float)[ok], np.asarray(val, dtype=float)[ok]
    o = np.argsort(p)
    p, v = p[o], v[o]
    c1 = np.concatenate([[0.0], np.cumsum(v)])
    c2 = np.concatenate([[0.0], np.cumsum(v * v)])
    lo = np.searchsorted(p, np.asarray(centres, dtype=float) - size / 2.0, side="left")
    hi = np.searchsorted(p, np.asarray(centres, dtype=float) + size / 2.0, side="right")
    n = (hi - lo).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        m = (c1[hi] - c1[lo]) / n
        var = (c2[hi] - c2[lo]) / n - m * m
    return np.where(n >= 5, np.sqrt(np.clip(var, 0, None)), np.nan)


def _track_x(mk, tr):
    """x (Mb, the cumulative axis of _genome_x) of the centre of every sliding window of the track."""
    _, centres, offsets, _ = _genome_x(mk)
    mn = mk.groupby("Chromosome")["Position"].min()
    return np.array([offsets[c] + (p - mn[c]) / 1e6 for c, p in zip(tr["Chromosome"], tr["Position"])]), centres, offsets, mn


def _nice_step(span_mb):
    for st in (0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100):
        if span_mb / st <= 12:
            return st
    return 100


def _region_x(offsets, mn, c, bp):
    return offsets[c] + (bp - mn[c]) / 1e6


def fig_snp_index(a, out):
    """A: SNP-index of the high pool, B: of the low pool, C: Delta SNP-index with the 95 / 99 % permutation intervals, D: the strongest region."""
    from matplotlib.gridspec import GridSpec
    from matplotlib.patches import Patch, Rectangle, ConnectionPatch
    from matplotlib.lines import Line2D
    mk, tr, regs = a["markers"], a["snp_track"], a["snp_regions"]
    if not tr["Delta_SNP_index"].notna().any():
        raise ValueError("no track")
    parent = a["snp_parent"]
    gw_regs = a["windows"][a["primary"]]["regions"]
    mkx, centres, offsets, _ = _genome_x(mk)
    tx, _c, _o, mn = _track_x(mk, tr)
    chroms = sorted(mk["Chromosome"].unique())
    colour = {c: C_CHR_PAIR[i % 2] for i, c in enumerate(chroms)}
    target = regs.iloc[0] if len(regs) else (gw_regs.iloc[0] if len(gw_regs) else None)
    with_d = target is not None
    fig = plt.figure(figsize=(13, 17.5 if with_d else 10.5))
    gs = GridSpec(4 if with_d else 3, 1, height_ratios=[3, 3, 3.6, 5.0] if with_d else [3, 3, 3.6], hspace=0.38, figure=fig)
    axA = fig.add_subplot(gs[0])
    axB = fig.add_subplot(gs[1], sharex=axA)
    axC = fig.add_subplot(gs[2], sharex=axA)
    fig.subplots_adjust(left=0.085, right=0.975, top=0.95, bottom=0.11 if with_d else 0.12)
    scored = mk["Scored"].to_numpy().astype(bool)
    prim_kb = a["snp_window_kb"]
    k = a["pool_size"]
    lo_cut, hi_cut = a["pool_cuts"]
    panels = ((axA, "A", "SNP_index_high", f"High pool SNP-index\n({k} lines, {a['trait']} ≥ {hi_cut:.4g})"),
              (axB, "B", "SNP_index_low", f"Low pool SNP-index\n({k} lines, {a['trait']} ≤ {lo_cut:.4g})"),
              (axC, "C", "Delta_SNP_index", "ΔSNP-index\n(high pool − low pool)"))
    x_lo, x_hi = 0.0, float(mkx.max())
    for ax, letter, col, ylab in panels:
        v = mk[col].to_numpy()
        for i, c in enumerate(chroms):                                   # a faint stripe behind every chromosome, as in the QTL-seq figures
            selc = (mk["Chromosome"] == c).to_numpy()
            ax.axvspan(mkx[selc].min() - 0.5, mkx[selc].max() + 0.5, color=C_CHR_BAND[i % 2], lw=0, zorder=0)
        for c in chroms:
            selc = (mk["Chromosome"] == c).to_numpy()
            sel = selc & scored & np.isfinite(v)
            t = tr[tr["Chromosome"] == c]
            mu = t[col].to_numpy()
            sd = _roll_sd(mk["Position"].to_numpy()[sel], v[sel], t["Position"].to_numpy(), prim_kb * 1000.0)
            ax.fill_between(tx[t.index], mu - sd, mu + sd, color=colour[c], alpha=0.30, lw=0, zorder=1)          # spread of the markers in the window
            ax.scatter(mkx[sel], v[sel], s=30, color=colour[c], alpha=0.6, edgecolors="white", linewidths=0.3, rasterized=True, zorder=2)
            ax.plot(tx[t.index], mu, color=C_LINE, lw=1.7, zorder=5, solid_capstyle="round")
        ax.set_ylabel(ylab, fontsize=10.5, fontweight="bold", color=C_LINE)
        ax.grid(axis="y", color="#E6E6EE", lw=0.7, zorder=0)
        ax.set_axisbelow(True)
        ax.text(-0.085, 1.0, letter, transform=ax.transAxes, fontsize=20, fontweight="bold", ha="left", va="top")
        ax.set_xticks(list(centres.values()))
        ax.set_xticklabels([f"chr{c}" for c in centres], fontsize=9)
    axA.set_ylim(-0.03, 1.03)
    axB.set_ylim(-0.03, 1.03)
    axC.set_ylim(-1.05, 1.05)
    allt = np.concatenate([[0], mkx])
    for ax in (axA, axB, axC):
        ax.set_xlim(allt.min() - 1, mkx.max() + 1)
    axC.axhline(0, color="#888888", lw=0.8)
    ci_ok = bool(tr["CI99_upper"].notna().any())
    if ci_ok:
        for c in chroms:
            t = tr[tr["Chromosome"] == c]
            for col, colr in (("CI95", C_CI95), ("CI99", C_CI99)):
                for side in ("upper", "lower"):
                    axC.plot(tx[t.index], t[f"{col}_{side}"], color=colr, lw=1.6, zorder=4)
    # regions: red = Delta SNP-index outside the 99 % interval; blue hatch = Gw region that the Delta SNP-index does not show
    span_min = 0.004 * (mkx.max() - allt.min())
    for r in regs.itertuples():
        x0, x1 = _region_x(offsets, mn, r.Chromosome, r.Left_bp), _region_x(offsets, mn, r.Chromosome, r.Right_bp)
        x1 = max(x1, x0 + span_min)
        for ax in (axA, axB, axC):
            ax.axvspan(x0, x1, color=C_REGION, alpha=0.38, lw=0, zorder=1)
    for r in gw_regs.itertuples():
        if len(regs) and ((regs["Chromosome"] == r.Chromosome) & (regs["Left_bp"] <= r.End_bp) & (regs["Right_bp"] >= r.Start_bp)).any():
            continue
        x0, x1 = _region_x(offsets, mn, r.Chromosome, r.Start_bp), _region_x(offsets, mn, r.Chromosome, r.End_bp)
        x1 = max(x1, x0 + span_min)
        for ax in (axA, axB, axC):
            ax.axvspan(x0, x1, facecolor="none", edgecolor=C_CI95, hatch="////", alpha=0.5, lw=0, zorder=1)
    # the position of the strongest regions, written above panel C
    for i, r in enumerate(regs.head(4).itertuples()):
        xm = (_region_x(offsets, mn, r.Chromosome, r.Left_bp) + _region_x(offsets, mn, r.Chromosome, r.Right_bp)) / 2
        axC.annotate(f"chr{r.Chromosome}: {r.Left_bp / 1e3:,.0f}–{r.Right_bp / 1e3:,.0f} Kb", xy=(xm, 1.0), xycoords=("data", "axes fraction"),
                     xytext=(xm, 1.07 + 0.09 * (i % 2)), textcoords=("data", "axes fraction"), ha="center", va="bottom", fontsize=8.5, fontweight="bold",
                     arrowprops=dict(arrowstyle="-|>", color="black", lw=0.9, shrinkA=0, shrinkB=0), annotation_clip=False)
    win_txt = f"{prim_kb / 1000:g} Mb" if prim_kb >= 1000 else f"{prim_kb:g} kb"
    handles = [Line2D([], [], color=C_LINE, lw=1.7, label=f"{win_txt} sliding window (shaded: ± 1 SD of the markers in it)")]
    if ci_ok:
        handles += [Line2D([], [], color=C_CI99, lw=1.6, label=f"99 % interval ({a['snp_null_perms']:,} permutations)"),
                    Line2D([], [], color=C_CI95, lw=1.6, label="95 % interval")]
    handles += [Patch(facecolor=C_REGION, alpha=0.4, label="peak region (outside the 99 % interval)")]
    if len(gw_regs):
        handles += [Patch(facecolor="none", edgecolor=C_CI95, hatch="////", label="Gw region (genotype counts) without a ΔSNP-index peak")]
    axC.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False, fontsize=8.5)

    if with_d:
        c = int(target["Chromosome"])
        left = int(target["Left_bp"]) if "Left_bp" in target else int(target["Start_bp"])
        right = int(target["Right_bp"]) if "Right_bp" in target else int(target["End_bp"])
        peak_m = target["Peak_marker"]
        chr_mk = mk[(mk["Chromosome"] == c) & mk["Scored"]]
        chr_len = max(chr_mk["Position"].max() / 1e6, right / 1e6)
        pk_row = chr_mk[chr_mk["Marker"] == peak_m]
        peak_bp = float(pk_row["Position"].iloc[0]) if len(pk_row) else (left + right) / 2.0
        # the enlarged view: at most two windows around the peak marker, inside the region
        half = prim_kb * 1000.0
        z_lo, z_hi = max(0.0, peak_bp - half), min(chr_len * 1e6, peak_bp + half)
        if z_hi - z_lo < 1000:
            z_lo, z_hi = peak_bp - 5000, peak_bp + 5000
        bb = gs[3].get_position(fig)
        sub = fig.add_gridspec(3, 1, left=0.13, right=0.94, top=bb.y1 - 0.03, bottom=bb.y0 + 0.012, height_ratios=[0.75, 0.9, 2.9], hspace=0.4)
        axR, axZ, axT = fig.add_subplot(sub[0]), fig.add_subplot(sub[1]), fig.add_subplot(sub[2])
        fig.text(0.0, bb.y1 - 0.004, "D", fontweight="bold", fontsize=20, ha="left", va="top")
        fig.text(0.04, bb.y1 - 0.009, f"Strongest region: chr{c}, {left / 1e3:,.0f}–{right / 1e3:,.0f} Kb", fontweight="bold", fontsize=11, ha="left", va="top")
        # ruler of the whole chromosome
        axR.set_xlim(-0.5, chr_len + 1.0)
        axR.set_ylim(-1.3, 1.4)
        axR.add_patch(Rectangle((0, 0), chr_len, 0.9, facecolor=C_RULER, edgecolor="#8A84A8", lw=1.2, zorder=2))
        axR.add_patch(Rectangle((left / 1e6, 0), max((right - left) / 1e6, chr_len * 0.004), 0.9, facecolor=C_REGION, edgecolor="none", alpha=0.55, zorder=3))
        axR.add_patch(Rectangle((z_lo / 1e6, 0), max((z_hi - z_lo) / 1e6, chr_len * 0.003), 0.9, facecolor=C_REGION, edgecolor="none", zorder=4))
        st = _nice_step(chr_len)
        for t in np.arange(0, chr_len + 1e-9, st):
            axR.plot([t, t], [0, -0.22], color="black", lw=0.9)
            axR.text(t, -0.34, f"{t:g}", ha="center", va="top", fontsize=8)
        axR.text(-0.5, 1.05, f"chr{c} (Mb)", ha="left", va="bottom", fontweight="bold", fontsize=10)
        axR.text(chr_len + 1.0, 1.05, f"{chr_len:.2f} Mb (to the last marker)", ha="right", va="bottom", fontsize=8.5)
        axR.axis("off")
        # the enlarged view
        span = z_hi - z_lo
        z0, z1 = (z_lo - 0.03 * span) / 1e3, (z_hi + 0.03 * span) / 1e3
        axZ.set_xlim(z0, z1)
        axZ.set_ylim(-1.0, 1.8)
        axZ.add_patch(Rectangle((z_lo / 1e3, 0), span / 1e3, 0.9, facecolor=C_RULER, edgecolor="#8A84A8", lw=1.0, zorder=1))
        r0, r1 = max(left, z_lo), min(right, z_hi)
        if r1 > r0:
            axZ.add_patch(Rectangle((r0 / 1e3, 0), (r1 - r0) / 1e3, 0.9, facecolor=C_REGION, alpha=0.85, edgecolor="none", zorder=2))
        inz = chr_mk[(chr_mk["Position"] >= z_lo) & (chr_mk["Position"] <= z_hi)]
        for p in inz["Position"]:
            axZ.plot([p / 1e3, p / 1e3], [0, 0.9], color="white" if left <= p <= right else "#BDB8D4", lw=0.9, zorder=3)
        axZ.plot([peak_bp / 1e3] * 2, [-0.15, 0.9], color=C_LINE, lw=2.4, zorder=4)
        axZ.text(z_lo / 1e3, 1.3, f"{z_lo / 1e3:,.0f} Kb", ha="left", va="bottom", fontsize=9, fontweight="bold")
        axZ.text(z_hi / 1e3, 1.3, f"{z_hi / 1e3:,.0f} Kb", ha="right", va="bottom", fontsize=9, fontweight="bold")
        axZ.text(peak_bp / 1e3, -0.25, f"peak {peak_bp / 1e3:,.0f} Kb", ha="center", va="top", fontsize=8.5, color="black")
        axZ.axis("off")
        for xe in (z_lo, z_hi):
            fig.add_artist(ConnectionPatch(xyA=(xe / 1e6, 0), coordsA=axR.transData, xyB=(xe / 1e3, 0.9), coordsB=axZ.transData, color="black", ls=":", lw=0.9))
        # marker names: the markers nearest to the peak, spread evenly below with a line to their position
        pool = inz.assign(_d=(inz["Position"] - peak_bp).abs()).sort_values("_d")
        chosen, min_gap = [], max(span / 40.0, 1.0)
        for r in pool.itertuples():
            if all(abs(r.Position - q.Position) >= min_gap for q in chosen):
                chosen.append(r)
            if len(chosen) == 7:
                break
        chosen.sort(key=lambda q: q.Position)
        axT.set_xlim(z0, z1)
        axT.set_ylim(0, 1)
        if chosen:
            xs = np.linspace(z0 + 0.05 * (z1 - z0), z1 - 0.05 * (z1 - z0), len(chosen)) if len(chosen) > 1 else [(z0 + z1) / 2]
            for r, xl in zip(chosen, xs):
                bold = r.Marker == peak_m
                axT.plot([r.Position / 1e3, xl], [1.0, 0.82], color="black" if bold else "#777777", lw=0.9 if bold else 0.7)
                axT.text(xl, 0.8, str(r.Marker), rotation=90, ha="center", va="top", fontsize=8.5, fontweight="bold" if bold else "normal")
        axT.axis("off")
        if chosen:                                                       # long marker names: shrink the font so that they stay above the info box
            h_in = axT.get_position().height * fig.get_figheight()
            longest = max(len(str(r.Marker)) for r in chosen)
            fs = float(np.clip(0.85 * h_in * 72.0 / (0.7 * longest), 6.5, 8.5))
            for tx_ in axT.texts:
                tx_.set_fontsize(fs)
        if len(regs) and "Delta_SNP_index_peak" in target:
            info = (f"peak region = ΔSNP-index outside the 99 % interval of the permutation test\n"
                    f"peak marker {peak_m}: SNP-index high {target['SNP_index_high_peak']:.2f}, low {target['SNP_index_low_peak']:.2f}, "
                    f"Δ {target['Delta_SNP_index_peak']:+.2f} ({target['Direction']})")
        else:
            info = ("Gw region without a ΔSNP-index peak (typical of an over-dominant locus: the allele frequency does not differ, the genotype counts do)\n"
                    f"peak marker {peak_m}")
        fig.text(0.13, bb.y0 - 0.032, info, fontsize=9, ha="left", va="top", bbox=dict(boxstyle="round,pad=0.6", facecolor="#FFF5F6", edgecolor=C_REGION, lw=1.3))
    fig.text(0.0, 0.004, f"SNP-index = frequency of the allele of parent {parent} (genotype {parent}) in a pool.   ΔSNP-index = SNP-index of the high pool − SNP-index of the low pool: "
             f"positive means that allele goes with the higher trait value.", fontsize=8.5, ha="left", va="bottom", color="#333333")
    return _save(fig, out, "snp_index_tracks.png", dpi=200)


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
    items = [("snpindex", "SNP-index of the high and the low pool and the ΔSNP-index along the genome", fig_snp_index), ("gw", "Window Gw along the genome, 250 kb / 500 kb / 1 Mb compared", fig_gw), ("loci", "Loci found by each window size", fig_loci),
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
