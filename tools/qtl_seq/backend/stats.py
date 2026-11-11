"""
Statistics of the SNP-index analysis (numpy / pandas only, no I/O).

SNP table
    One row per informative SNP with the columns
    Chromosome, Position, Ref, Alt, Ref_High, Alt_High, Ref_Low, Alt_Low
    (read counts of the two bulks). "Alt" is the allele the SNP-index counts (the index allele).

    SNP-index_H = Alt_High / (Ref_High + Alt_High)       SNP-index_L = Alt_Low / (Ref_Low + Alt_Low)
    Delta SNP-index = SNP-index_H - SNP-index_L           (-1 ... +1)

    Positive Delta: the index allele is enriched in the HIGH bulk. Negative: enriched in the LOW bulk.

Null distribution (no QTL)
    For one bulk of n plants the share of index alleles in the bulk is
        F2         k ~ Binomial(2n, 1/2),  f = k / 2n          (every plant carries two alleles)
        RIL        k ~ Binomial(n, 1/2),   f = k / n           (every plant is homozygous)
        backcross  k ~ Binomial(n, 1/2),   f = k / 2n  (recurrent parent without the index allele)
                                           f = (n + k) / 2n    (recurrent parent with the index allele)
    The sequencing then samples d reads from that bulk: alt reads ~ Binomial(d, f). Doing this for both bulks
    and subtracting gives one simulated Delta SNP-index; `replicates` of them (10,000 by default) give the null
    distribution for that pair of read depths. The 90/95/99 % limits are its two-sided quantiles.

    Two ways to turn that into limits for a sliding window (selectable, both reported in the parameter file):
    "snp_average"  the limits of every SNP in the window (each from its own pair of read depths) are averaged,
                   the same way the Delta values are. Low coverage widens the limits; the number of SNPs does not
                   narrow them, which is the conservative choice.
    "window"       the window is simulated as a whole: one bulk allele frequency per bulk (the SNPs of a window are
                   tightly linked) and the reads of all m SNPs of the window pooled. More SNPs narrow the limits.

    Assumptions that cannot be recovered from a description of the method and are therefore explicit here:
    random sampling of the plants of each bulk under the null (no selection), no residual heterozygosity in RILs,
    independence between the two bulks, no sequencing/alignment error in the simulation, the per-SNP lookup uses
    the (binned) read depths of both bulks, and no multiple-testing correction across windows.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from tools.qtl_seq.config import LEVELS, MAX_SIM_DEPTH

SNP_COLUMNS = ["Chromosome", "Position", "Ref", "Alt", "Ref_High", "Alt_High", "Ref_Low", "Alt_Low"]


# ============================================================
# SNP-INDEX
# ============================================================

def snp_index(alt, ref):
    """alt / (alt + ref); NaN where nothing covers the SNP."""
    alt = np.asarray(alt, dtype=float)
    depth = alt + np.asarray(ref, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(depth > 0, alt / depth, np.nan)


def delta_snp_index(index_high, index_low):
    """SNP-index(high bulk) - SNP-index(low bulk)."""
    return np.asarray(index_high, dtype=float) - np.asarray(index_low, dtype=float)


def add_indices(df):
    """Adds Depth_High/Low, SNP_index_High/Low and Delta_SNP_index to a SNP table (returns a copy)."""
    out = df.copy()
    out["Depth_High"] = out["Ref_High"] + out["Alt_High"]
    out["Depth_Low"] = out["Ref_Low"] + out["Alt_Low"]
    out["SNP_index_High"] = snp_index(out["Alt_High"], out["Ref_High"])
    out["SNP_index_Low"] = snp_index(out["Alt_Low"], out["Ref_Low"])
    out["Delta_SNP_index"] = delta_snp_index(out["SNP_index_High"], out["SNP_index_Low"])
    return out


def apply_filters(df, min_depth=7, max_depth=0, min_snp_index=0.3, use_snp_index=True):
    """Depth and SNP-index filters, applied in this order; every removal is counted once.

    1. read depth below `min_depth` in either bulk
    2. read depth above `max_depth` in either bulk (0 = no upper limit)
    3. SNP-index below `min_snp_index` in BOTH bulks. A SNP is kept when at least one bulk reaches the threshold:
       the SNP-index of the other bulk is kept as it is, even when it is below the threshold.
    Returns (filtered table with indices, ledger)."""
    d = add_indices(df)
    ledger = {"input": int(len(d))}
    keep = (d["Depth_High"] >= min_depth) & (d["Depth_Low"] >= min_depth)
    ledger["removed_low_depth"] = int((~keep).sum())
    d = d[keep]
    if max_depth and max_depth > 0:
        keep = (d["Depth_High"] <= max_depth) & (d["Depth_Low"] <= max_depth)
        ledger["removed_high_depth"] = int((~keep).sum())
        d = d[keep]
    else:
        ledger["removed_high_depth"] = 0
    if use_snp_index:
        keep = np.maximum(d["SNP_index_High"], d["SNP_index_Low"]) >= min_snp_index
        ledger["removed_snp_index"] = int((~keep).sum())
        d = d[keep]
    else:
        ledger["removed_snp_index"] = 0
    d = d.reset_index(drop=True)
    ledger["retained"] = int(len(d))
    return d, ledger


# ============================================================
# NULL DISTRIBUTION
# ============================================================

def bulk_frequency(population, n_plants, reps, rng, bc_recurrent="non_index"):
    """Index-allele share of a bulk of `n_plants` random plants, `reps` times."""
    n = int(n_plants)
    if n < 1:
        raise ValueError("A bulk needs at least one plant.")
    pop = str(population).upper()
    if pop == "F2":
        return rng.binomial(2 * n, 0.5, reps) / (2.0 * n)
    if pop == "RIL":
        return rng.binomial(n, 0.5, reps) / float(n)
    if pop in ("BC", "BACKCROSS"):
        k = rng.binomial(n, 0.5, reps)
        return ((n + k) if bc_recurrent == "index" else k) / (2.0 * n)
    raise ValueError(f"Unknown population type: {population}")


def simulate_null(depth_high, depth_low, population, n_high, n_low, reps, rng, bc_recurrent="non_index"):
    """`reps` simulated Delta SNP-index values for one SNP covered by depth_high / depth_low reads."""
    fh = bulk_frequency(population, n_high, reps, rng, bc_recurrent)
    fl = bulk_frequency(population, n_low, reps, rng, bc_recurrent)
    ah = rng.binomial(int(depth_high), fh)
    al = rng.binomial(int(depth_low), fl)
    return ah / float(depth_high) - al / float(depth_low)


def level_quantiles(levels=LEVELS):
    """[(level, lower quantile, upper quantile)] — two-sided limits."""
    return [(lv, (1 - lv / 100.0) / 2.0, 1 - (1 - lv / 100.0) / 2.0) for lv in levels]


def bin_depth(depth):
    """Depths up to 50 are simulated exactly, up to 100 in steps of 5, above in steps of 10, capped at MAX_SIM_DEPTH.
    (The null distribution hardly changes between neighbouring high depths; binning keeps the number of
    simulated depth pairs small.)"""
    d = np.minimum(np.asarray(depth, dtype=float), MAX_SIM_DEPTH)
    d = np.where(d <= 50, np.rint(d), np.where(d <= 100, np.rint(d / 5) * 5, np.rint(d / 10) * 10))
    return np.maximum(d, 1).astype(np.int64)


class LimitTable:
    """Simulated limits for a set of (depth_high, depth_low) pairs."""

    def __init__(self, frame, meta):
        self.frame = frame
        self.meta = meta
        self._keys = (frame["Depth_High"].to_numpy(np.int64) * 100_000 + frame["Depth_Low"].to_numpy(np.int64))
        order = np.argsort(self._keys)
        self._order = order
        self._sorted = self._keys[order]

    def lookup(self, depth_high, depth_low):
        """Row positions of the table for SNPs with these (unbinned) depths."""
        key = bin_depth(depth_high) * 100_000 + bin_depth(depth_low)
        i = np.searchsorted(self._sorted, key)
        return self._order[np.clip(i, 0, len(self._sorted) - 1)]

    def limits_for(self, depth_high, depth_low):
        """{column: array} of the CIxx_Lower / CIxx_Upper limits of each SNP."""
        rows = self.lookup(depth_high, depth_low)
        return {c: self.frame[c].to_numpy()[rows] for c in self.frame.columns if c.startswith("CI")}


def simulate_limits(depth_high, depth_low, population, n_high, n_low, reps=10_000, seed=12345,
                    bc_recurrent="non_index", levels=LEVELS, progress=None, cancelled=None):
    """Run the null simulation for every distinct (binned) pair of read depths in the data.

    Returns a LimitTable; its frame has Depth_High, Depth_Low, Replicates and CI<level>_Lower / CI<level>_Upper."""
    dh, dl = bin_depth(depth_high), bin_depth(depth_low)
    pairs = np.unique(np.stack([dh, dl], axis=1), axis=0)          # sorted -> reproducible for one seed
    rng = np.random.default_rng(int(seed))
    qs = level_quantiles(levels)
    cols = {"Depth_High": pairs[:, 0], "Depth_Low": pairs[:, 1], "Replicates": np.full(len(pairs), int(reps))}
    for lv, _lo, _up in qs:
        cols[f"CI{lv}_Lower"] = np.empty(len(pairs))
        cols[f"CI{lv}_Upper"] = np.empty(len(pairs))
    probs = [q for _lv, lo, up in qs for q in (lo, up)]
    for i, (a, b) in enumerate(pairs):
        if cancelled is not None and cancelled():
            raise InterruptedError("cancelled")
        sim = simulate_null(a, b, population, n_high, n_low, int(reps), rng, bc_recurrent)
        vals = np.quantile(sim, probs)
        for j, (lv, _lo, _up) in enumerate(qs):
            cols[f"CI{lv}_Lower"][i] = vals[2 * j]
            cols[f"CI{lv}_Upper"][i] = vals[2 * j + 1]
        if progress is not None and (i % 200 == 0):
            progress(i / max(len(pairs), 1))
    meta = {"population": population, "n_high": int(n_high), "n_low": int(n_low), "replicates": int(reps),
            "seed": int(seed), "bc_recurrent": bc_recurrent, "pairs": int(len(pairs))}
    return LimitTable(pd.DataFrame(cols), meta)


def fixed_depth_limits(depth, *args, **kw):
    """One depth for every SNP (the user typed the sequencing depth instead of using the observed depths)."""
    d = np.array([int(depth)])
    return simulate_limits(d, d, *args, **kw)


def window_limits(counts, depth_high, depth_low, population, n_high, n_low, reps=10_000, seed=12345,
                  bc_recurrent="non_index", levels=LEVELS, cancelled=None):
    """Limits of windows simulated as a whole ("window" method).

    counts = SNPs in each window, depth_high / depth_low = mean read depth of the window's SNPs. Per window: one
    index-allele share per bulk, then the reads of all SNPs pooled (m SNPs of depth d = Binomial(m*d, f) / (m*d)).
    Windows with similar (m, d, d) are simulated once. Returns {column: array}."""
    m = np.maximum(np.asarray(counts, dtype=float), 1)
    # SNP counts are grouped in ~3 % steps and depths rounded, which keeps the number of simulations small
    mb = np.maximum(np.rint(np.exp(np.rint(np.log(m) / 0.03) * 0.03)), 1).astype(np.int64)
    dh = np.maximum(np.rint(np.asarray(depth_high, dtype=float)), 1).astype(np.int64)
    dl = np.maximum(np.rint(np.asarray(depth_low, dtype=float)), 1).astype(np.int64)
    keys = np.stack([mb, dh, dl], axis=1)
    uniq, inverse = np.unique(keys, axis=0, return_inverse=True)
    rng = np.random.default_rng(int(seed))
    qs = level_quantiles(levels)
    probs = [q for _lv, lo, up in qs for q in (lo, up)]
    table = np.empty((len(uniq), len(probs)))
    for i, (mm, a, b) in enumerate(uniq):
        if cancelled is not None and cancelled():
            raise InterruptedError("cancelled")
        fh = bulk_frequency(population, n_high, int(reps), rng, bc_recurrent)
        fl = bulk_frequency(population, n_low, int(reps), rng, bc_recurrent)
        nh, nl = int(mm * a), int(mm * b)
        sim = rng.binomial(nh, fh) / float(nh) - rng.binomial(nl, fl) / float(nl)
        table[i] = np.quantile(sim, probs)
    inverse = np.asarray(inverse).reshape(-1)
    out = {}
    for j, (lv, _lo, _up) in enumerate(qs):
        out[f"CI{lv}_Lower"] = table[inverse, 2 * j]
        out[f"CI{lv}_Upper"] = table[inverse, 2 * j + 1]
    return out


# ============================================================
# SLIDING WINDOWS
# ============================================================

def window_starts(chrom_length, size, step):
    """1-based start positions of the windows of one chromosome: every `step` bp while the window fits."""
    chrom_length, size, step = int(chrom_length), int(size), int(step)
    if chrom_length <= size:
        return np.array([1], dtype=np.int64)
    return np.arange(1, chrom_length - size + 2, step, dtype=np.int64)


def sliding_windows(positions, values, chrom_length, size=2_000_000, step=10_000, min_snps=10):
    """Mean of every column of `values` over the SNPs inside each window.

    positions: sorted 1-based SNP positions of one chromosome. values: {name: array aligned to positions}.
    A window [s, s + size - 1] holds the SNPs inside it. Windows with fewer than `min_snps` SNPs are dropped.
    Returns a dict: Window_Start, Window_End, SNP_Count and one array per value name (all aligned)."""
    pos = np.asarray(positions, dtype=np.int64)
    starts = window_starts(chrom_length, size, step)
    ends = np.minimum(starts + size - 1, int(chrom_length)) if chrom_length > size else np.array([int(chrom_length)])
    lo = np.searchsorted(pos, starts, side="left")
    hi = np.searchsorted(pos, ends, side="right")
    n = hi - lo
    ok = n >= max(int(min_snps), 1)
    out = {"Window_Start": starts[ok], "Window_End": ends[ok], "SNP_Count": n[ok]}
    for name, v in values.items():
        v = np.asarray(v, dtype=float)
        cs = np.concatenate([[0.0], np.cumsum(v)])
        with np.errstate(divide="ignore", invalid="ignore"):
            mean = (cs[hi] - cs[lo]) / np.maximum(n, 1)
        out[name] = mean[ok]
    return out


def level_label(level):
    return f"{int(level)}%"


def p_threshold(level):
    return round(1 - level / 100.0, 4)


def build_windows(snps, chrom_lengths, p, limits=None, progress=None, cancelled=None):
    """Sliding windows of the whole genome with confidence limits and significance flags.

    snps: filtered SNP table with indices (apply_filters). chrom_lengths: {chrom: length}.
    limits: LimitTable for the "snp_average" method (the limits of each SNP are averaged in each window) —
    for the "window" method pass None and the window limits are simulated from the window's SNP count and depths.
    Returns a DataFrame with the all_sliding_windows columns (and the depth columns used by the window method)."""
    levels = LEVELS
    method = p.get("ci_method", "snp_average")
    frames = []
    snp_lim = None
    if method == "snp_average":
        # the per-SNP limits were looked up in the simulation table beforehand (columns CI<level>_Lower / _Upper)
        cols = [f"CI{lv}_{s}" for lv in levels for s in ("Lower", "Upper")]
        if any(c not in snps.columns for c in cols):
            if limits is None:
                raise ValueError("snp_average needs per-SNP limits or a LimitTable")
            snps = snps.assign(**limits.limits_for(snps["Depth_High"].to_numpy(), snps["Depth_Low"].to_numpy()))
        snp_lim = {c: snps[c].to_numpy() for c in cols}
    present = set(snps["Chromosome"])
    chroms = [c for c in chrom_lengths if c in present]
    for ci, chrom in enumerate(chroms):
        if cancelled is not None and cancelled():
            raise InterruptedError("cancelled")
        idx = np.flatnonzero((snps["Chromosome"] == chrom).to_numpy())
        sub = snps.iloc[idx]
        order = np.argsort(sub["Position"].to_numpy(), kind="stable")
        sub, idx = sub.iloc[order], idx[order]
        vals = {"Mean_SNP_index_High": sub["SNP_index_High"].to_numpy(),
                "Mean_SNP_index_Low": sub["SNP_index_Low"].to_numpy(),
                "Mean_Delta_SNP_index": sub["Delta_SNP_index"].to_numpy(),
                "Mean_Depth_High": sub["Depth_High"].to_numpy(float),
                "Mean_Depth_Low": sub["Depth_Low"].to_numpy(float)}
        if snp_lim is not None:
            for c, arr in snp_lim.items():
                vals[c] = arr[idx]
        w = sliding_windows(sub["Position"].to_numpy(), vals, chrom_lengths[chrom], p["window_size"],
                            p["window_step"], p["min_snps_per_window"])
        if len(w["Window_Start"]) == 0:
            continue
        f = pd.DataFrame(w)
        f.insert(0, "Chromosome", chrom)
        frames.append(f)
        if progress is not None:
            progress((ci + 1) / max(len(chroms), 1))
    if not frames:
        return pd.DataFrame(columns=window_columns())
    win = pd.concat(frames, ignore_index=True)
    if method == "window":
        lim = window_limits(win["SNP_Count"], win["Mean_Depth_High"], win["Mean_Depth_Low"], p["population"],
                            p["n_high"], p["n_low"], p["replicates"], p["seed"], p.get("bc_recurrent", "non_index"),
                            cancelled=cancelled)
        for c, arr in lim.items():
            win[c] = arr
    return flag_windows(win)


def window_columns():
    cols = ["Chromosome", "Window_Start", "Window_End", "SNP_Count", "Mean_SNP_index_High", "Mean_SNP_index_Low",
            "Mean_Delta_SNP_index"]
    for lv in LEVELS:
        cols += [f"CI{lv}_Upper", f"CI{lv}_Lower"]
    return cols + [f"Significant_{lv}" for lv in LEVELS]


def flag_windows(win):
    """Significant_<level> = the window mean lies outside the lower/upper limit of that level."""
    win = win.copy()
    m = win["Mean_Delta_SNP_index"]
    for lv in LEVELS:
        win[f"Significant_{lv}"] = (m > win[f"CI{lv}_Upper"]) | (m < win[f"CI{lv}_Lower"])
    win["Window_Mid"] = ((win["Window_Start"] + win["Window_End"]) // 2).astype(np.int64)
    return win


# ============================================================
# QTL INTERVALS
# ============================================================

QTL_COLUMNS = ["Chromosome", "Start", "End", "Peak_Position", "Peak_Delta_SNP_index", "Mean_Delta_SNP_index",
               "Minimum_Delta_SNP_index", "Maximum_Delta_SNP_index", "SNP_Count", "Significance_Level",
               "Confidence_Threshold"]


def significant_windows(win, min_level=90):
    """The windows that cross the limits of `min_level` or a stricter one, with direction and highest level."""
    w = win.copy()
    w["Sign"] = np.where(w["Mean_Delta_SNP_index"] > 0, 1, -1)
    w["Level"] = 0
    for lv in sorted(LEVELS):
        w.loc[w[f"Significant_{lv}"], "Level"] = lv
    return w[w["Level"] >= int(min_level)].copy()


def call_qtls(win, snp_positions, window_step=10_000, min_level=90, merge_gap=100_000, basis="window_extent"):
    """Merge neighbouring significant windows into QTL intervals.

    Rule: significant windows (at `min_level` or stricter) of the same chromosome and the same sign belong to one
    QTL when no more than `merge_gap` bp of non-significant (or missing) windows lie between them. A window that
    is not significant, or a window of the opposite sign, ends the interval.
    basis (what Start/End report; the outer extent of all merged windows is always kept in Window_Extent_*):
      "top_level"       the outer ends of the first/last window that is significant at the HIGHEST level the QTL reaches
                        (a QTL that reaches 99 % is bounded by the stretch above the 99 % limit)
      "window_extent"   the outer ends of the first/last window significant at the lowest level used to call QTLs
      "window_midpoints" the midpoints of those windows (the stretch where the plotted curve is outside the limits)
    The peak is the significant window with the largest |mean Delta|; it is a position inside the interval, not its
    boundary.
    snp_positions: {chrom: sorted array of SNP positions} used for SNP_Count (SNPs inside Start..End).
    Returns (qtl table, significant-window table with a QTL_ID column)."""
    cols = QTL_COLUMNS + ["QTL_ID", "Direction", "Window_Count", "Window_Extent_Start", "Window_Extent_End",
                          "Peak_Window_Start", "Peak_Window_End", "P_Threshold", "Interval_Basis"]
    sig = significant_windows(win, min_level).sort_values(["Chromosome", "Window_Start"], kind="stable")
    if sig.empty:
        return pd.DataFrame(columns=cols), sig.assign(QTL_ID=pd.Series(dtype=str))
    qid = np.zeros(len(sig), dtype=np.int64)
    chrom = sig["Chromosome"].to_numpy()
    start = sig["Window_Start"].to_numpy()
    sign = sig["Sign"].to_numpy()
    cur = 0
    for i in range(len(sig)):
        if i == 0 or chrom[i] != chrom[i - 1] or sign[i] != sign[i - 1] or start[i] - start[i - 1] > window_step + merge_gap:
            cur += 1
        qid[i] = cur
    sig = sig.assign(_q=qid)
    rows, ids = [], []
    # a stable chromosome order for the numbering follows the order of the genome in `win`
    for n, (q, g) in enumerate(sig.groupby("_q", sort=True), start=1):
        peak = g.loc[g["Mean_Delta_SNP_index"].abs().idxmax()]
        top = int(g["Level"].max())
        positive = bool(peak["Mean_Delta_SNP_index"] > 0)
        ext_s, ext_e = int(g["Window_Start"].min()), int(g["Window_End"].max())
        if basis == "window_midpoints":
            s, e = int(g["Window_Mid"].min()), int(g["Window_Mid"].max())
        elif basis == "top_level":
            tl = g[g["Level"] == top]
            s, e = int(tl["Window_Start"].min()), int(tl["Window_End"].max())
        else:
            s, e = ext_s, ext_e
        span = win[(win["Chromosome"] == peak["Chromosome"]) & (win["Window_Start"] >= g["Window_Start"].min())
                   & (win["Window_Start"] <= g["Window_Start"].max())]["Mean_Delta_SNP_index"]
        pos = np.asarray(snp_positions.get(peak["Chromosome"], []))
        n_snps = int(np.searchsorted(pos, e, side="right") - np.searchsorted(pos, s, side="left")) if len(pos) else 0
        thr_col = f"CI{top}_Upper" if positive else f"CI{top}_Lower"
        rows.append({
            "Chromosome": peak["Chromosome"], "Start": s, "End": e,
            "Peak_Position": int(peak["Window_Mid"]), "Peak_Delta_SNP_index": float(peak["Mean_Delta_SNP_index"]),
            "Mean_Delta_SNP_index": float(span.mean()), "Minimum_Delta_SNP_index": float(span.min()),
            "Maximum_Delta_SNP_index": float(span.max()), "SNP_Count": n_snps,
            "Significance_Level": level_label(top), "Confidence_Threshold": float(peak[thr_col]),
            "QTL_ID": f"QTL{n}",
            "Direction": ("Index allele enriched in the high bulk (positive Delta)" if positive
                          else "Index allele enriched in the low bulk (negative Delta)"),
            "Window_Count": int(len(g)), "Window_Extent_Start": ext_s, "Window_Extent_End": ext_e,
            "Peak_Window_Start": int(peak["Window_Start"]), "Peak_Window_End": int(peak["Window_End"]),
            "P_Threshold": p_threshold(top), "Interval_Basis": basis,
        })
        ids.extend([f"QTL{n}"] * len(g))
    sig = sig.drop(columns="_q").assign(QTL_ID=ids)
    return pd.DataFrame(rows)[cols], sig


def level_counts(qtls):
    """Number of QTLs that reach each confidence level (a 99 % QTL also counts for 95 % and 90 %)."""
    lv = qtls["Significance_Level"].str.rstrip("%").astype(int) if len(qtls) else pd.Series([], dtype=int)
    return {str(L): int((lv >= L).sum()) for L in LEVELS}


def typical_pair(frame):
    """The (depth_high, depth_low) row of a limit table closest to the median depth (used to draw one null distribution)."""
    dm = float(np.median(frame["Depth_High"].to_numpy()))
    d = (frame["Depth_High"] - dm).abs() + (frame["Depth_Low"] - dm).abs()
    return frame.loc[d.idxmin()]


def null_histogram(limits, bins=32, seed=1):
    """Histogram of the simulated null distribution at the typical read depth, with that depth's limits."""
    meta, row = limits.meta, typical_pair(limits.frame)
    dh, dl = int(row["Depth_High"]), int(row["Depth_Low"])
    sim = simulate_null(dh, dl, meta["population"], meta["n_high"], meta["n_low"], int(meta["replicates"]),
                        np.random.default_rng(seed), meta.get("bc_recurrent", "non_index"))
    span = max(abs(float(row["CI99_Lower"])), abs(float(row["CI99_Upper"]))) * 1.35
    counts, edges = np.histogram(sim, bins=bins, range=(-span, span))
    return {"depth_high": dh, "depth_low": dl, "edges": [round(float(e), 4) for e in edges], "counts": counts.tolist(),
            "limits": {str(lv): [round(float(row[f"CI{lv}_Lower"]), 4), round(float(row[f"CI{lv}_Upper"]), 4)] for lv in LEVELS}}
