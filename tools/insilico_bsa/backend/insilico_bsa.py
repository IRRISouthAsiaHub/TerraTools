"""
In silico BSA — bulked segregant analysis made from the genotypes and phenotypes of every line of a mapping population.

No sequencing of pools is needed: the lines are sorted by the trait, the lowest and the highest fraction form the two
"pools" (in silico), and the genotype counts of every marker in the two pools are compared. This is how Li et al. (2022)
analysed their real rice IMF2 data with dQTG-seq (Plant Communications 3:100319): with the genotype of every individual known,
the allele and genotype counts are counts of individuals, not read counts.

Statistics per marker (pool x genotype contingency tables; no factor 2 in G, as in the paper):

    G1 = sum_ij  nQ_ij ln( nQ_ij / E(nQ_ij) )       allele counts  (2 x 2)   - the allele-frequency G statistic
    G2 = sum_ij  nG_ij ln( nG_ij / E(nG_ij) )       genotype counts (2 x 3)
    Gw = G1/(G1+G2) * G1 + G2/(G1+G2) * G2          weighted mean (dQTG-seq2): sees over-dominant loci that G1 cannot
    dAF = AF(high pool) - AF(low pool)              allele frequency difference (the Delta SNP-index of QTL-seq)

Window analysis (no per-marker testing, no kernel smoothing): the genome is cut into non-overlapping physical windows of 250 kb,
500 kb and 1 Mb. The statistic of a window is the mean of the marker statistic over the scored markers it holds (a window with fewer
than min_window_markers scored markers is not tested; automatic = 3, fewer on a sparse map). The threshold of each window size comes from permutation: in every permutation
a low and a high pool of the same size are drawn at random from all lines, the statistics are recomputed and the largest window value of
the genome is recorded, T_b = max_w Gw_wb. The threshold is the (1 - alpha) quantile of the T_b of that window size, and the adjusted
p-value of window w is  P_adj,w = (1 + #{b : T_b >= Gw_w}) / (B + 1). The same permutations serve all window sizes. Significant windows
of a chromosome that touch (or lie within the merge distance) form a region; regions of the three window sizes are then compared locus
by locus (found at which size, how long, how significant, where the peak lies).

Populations (chosen explicitly, see POPULATIONS): F2 / F3 / IMF2 and BC1 carry two alleles per line and use the heterozygous
calls (three genotype classes, two for a backcross); RIL and DH carry one allele per line, heterozygous calls are treated as
missing and Gw = G1. The observed heterozygosity is checked against what the chosen population should show.

Markers that cannot be scored (too few called lines in a pool) get no statistic (empty, never zero) and do not enter a window mean.

SNP-index view (as in QTL-seq): the SNP-index of a pool is the frequency of the allele of the chosen parent (parent 2 by default) in that pool,
so every marker has a SNP-index in the high pool and in the low pool, and Delta SNP-index = SNP-index(high pool) - SNP-index(low pool): positive
when that parent's allele is more frequent in the high pool. The three tracks are averaged in a sliding window (the primary window size, moved in
steps of a fifth of it). The confidence intervals of the Delta SNP-index come from the same permutations as the Gw threshold: at every position
the 2.5 / 97.5 % (95 % interval) and the 0.5 / 99.5 % (99 % interval) quantiles of the smoothed Delta of random pools. Positions outside the 99 %
interval form the Delta SNP-index peak regions.
"""

from __future__ import annotations

import os
import re
import warnings
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

MAX_ROWS = 400
MAX_PERMUTATIONS = 20000
WINDOW_SIZES_KB = (250, 500, 1000)
STATS = ("Gw", "G1", "dAF")
SLIDE_STEPS = 5                     # sliding window = primary window, moved in steps of window / SLIDE_STEPS
MAX_NULL_PERMUTATIONS = 2000        # permutations kept for the pointwise Delta SNP-index intervals
SLIDE_MIN_MARKERS = 5               # the SNP-index window is widened until at least SLIDE_FILL of the windows hold this many markers ...
SLIDE_FILL = 0.75
SLIDE_MAX_KB = 10000.0              # ... but never beyond 10 Mb
SLIDE_STEPS_KB = (250.0, 500.0, 750.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 8000.0, 10000.0)
STAT_LABEL = {"Gw": "Gw (alleles + genotypes)", "G1": "G1 (alleles only)", "dAF": "|Δ allele frequency|"}


from tools.insilico_bsa.backend.json_util import json_safe        # noqa: E402,F401


class BsaError(ValueError):
    """A problem the user can fix; shown as it is in the app."""


# Segregation models. diploid: every line carries two alleles and heterozygous calls are used; otherwise one allele per line and
# heterozygous calls are missing. het / tol: the heterozygosity expected from the cross (and how far the data may differ before
# the choice is questioned); classes: the genotype classes the population can show.
POPULATIONS = {
    "f2": {"label": "F2 / F3 / IMF2 (three genotype classes, 1:2:1)", "short": "F2 / IMF2", "diploid": True, "het": 0.50, "tol": 0.20},
    "bc1": {"label": "Backcross BC1 (two genotype classes, 1:1)", "short": "backcross", "diploid": True, "het": 0.50, "tol": 0.20},
    "ril": {"label": "RIL (inbred lines, one allele per line)", "short": "RIL", "diploid": False, "het": 0.0, "tol": 0.10},
    "dh": {"label": "DH (doubled haploid, one allele per line)", "short": "DH", "diploid": False, "het": 0.0, "tol": 0.05},
}
LEGACY_POPULATION = {"diploid": "f2", "inbred": "ril"}
AUTO_HET_CUTOFF = 0.15


def class_fractions(G: np.ndarray) -> list:
    n = max(int((G >= 0).sum()), 1)
    return [float((G == c).sum()) / n for c in (0, 1, 2)]


def suggest_population(G: np.ndarray) -> str:
    """A hint for the user, never applied silently: heterozygosity separates inbred lines from F2 / backcross; one empty homozygous class
    separates a backcross from an F2."""
    f0, f1, f2 = class_fractions(G)
    if f1 <= AUTO_HET_CUTOFF:
        return "ril"
    return "bc1" if min(f0, f2) < 0.2 * max(f0, f2) else "f2"


def resolve_population(selected, G: np.ndarray):
    """-> (key, notes). An explicit choice is used as it is; 'auto' (kept for old settings) falls back on a 15 % heterozygosity cut-off and says so."""
    sel = LEGACY_POPULATION.get(str(selected or "auto").lower(), str(selected or "auto").lower())
    notes = []
    if sel == "auto":
        sel = "f2" if class_fractions(G)[1] > AUTO_HET_CUTOFF else "ril"
        notes.append(f"The population type was not chosen; it was guessed from the heterozygosity (above {AUTO_HET_CUTOFF * 100:.0f} % = F2, otherwise RIL) "
                     f"and set to '{POPULATIONS[sel]['short']}'. Choose the population type yourself.")
    if sel not in POPULATIONS:
        raise BsaError(f"Unknown population type '{selected}'. Choose one of: " + ", ".join(POPULATIONS) + ".")
    return sel, notes


def population_check(G: np.ndarray, key: str) -> list:
    """Compares the genotype classes of the data with what the chosen population should show."""
    pop = POPULATIONS[key]
    f0, f1, f2 = class_fractions(G)
    notes = []
    if abs(f1 - pop["het"]) > pop["tol"]:
        notes.append(f"{f1 * 100:.1f} % of the calls are heterozygous, but a {pop['short']} population should show about {pop['het'] * 100:.0f} %. "
                     "Check that the population type is the right one.")
    if key == "bc1" and min(f0, f2) > 0.10:
        notes.append(f"Both homozygous classes are common ({f0 * 100:.0f} % and {f2 * 100:.0f} %); a backcross carries only one of them. This looks like an F2: "
                     "choose F2 / F3 / IMF2.")
    if key == "f2" and min(f0, f2) < 0.05 and f1 > 0.3:
        notes.append("One homozygous class is almost absent (an F2 shows 1:2:1); this looks like a backcross: choose BC1.")
    return notes


# ───────────────────────────── reading ─────────────────────────────

def _read_pheno_table(path, sheet=None) -> pd.DataFrame:
    low = str(path).lower()
    try:
        if low.endswith((".xlsx", ".xlsm", ".xls")):
            return pd.read_excel(path, sheet_name=sheet or 0)
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            head = fh.readline()
        sep = max((",", "\t", ";"), key=head.count)
        return pd.read_csv(path, sep=sep)
    except Exception as e:                                       # noqa: BLE001
        raise BsaError(f"Could not read the phenotype file: {e}") from None


def read_phenotypes(path, sheet=None) -> dict:
    """-> {'table': DataFrame (first column = sample), 'traits': [numeric columns]}."""
    if not path or not os.path.isfile(str(path)):
        raise BsaError("The phenotype file was not found.")
    df = _read_pheno_table(path, sheet)
    if df.shape[1] < 2 or df.shape[0] < 5:
        raise BsaError("The phenotype file needs a sample column, at least one trait column and at least 5 rows.")
    df.columns = [str(c).strip() for c in df.columns]
    id_col = df.columns[0]
    df[id_col] = df[id_col].astype(str).str.strip()
    traits = []
    for c in df.columns[1:]:
        v = pd.to_numeric(df[c], errors="coerce")
        if v.notna().sum() >= 5 and v.notna().mean() >= 0.5 and v.nunique() > 2:
            df[c] = v
            traits.append(c)
    if not traits:
        raise BsaError("No numeric trait column was found (the first column is taken as the sample name).")
    return {"table": df, "traits": traits, "id_col": id_col}


def read_genotypes(path, sheet=None) -> dict:
    """-> dict(markers, samples, G int8 markers x samples (-1 missing), pheno {sample: value} or None, traits, layout, notes).

    Two layouts: a table with one row per marker (CHROM/POS or Marker/Chromosome, then one column per line; what the Genofile
    Converter writes) or the GWAS-ready table of one row per line (ID, Pheno, markers), which holds the trait itself."""
    from tools.recombination.backend import recombination as R
    path = str(path)
    if not os.path.isfile(path):
        raise BsaError("The genotype file was not found.")
    try:
        d = R.read_genotypes(path, sheet)
        notes = list(d["notes"])
        if str(d["pos_source"]).startswith("column"):          # the column is the truth; say so when the names disagree with it
            named = d["markers"]["ID"].map(lambda v: R._NAME_POS.match(str(v)))
            from_name = named.map(lambda mm: float(mm.group(2)) if mm else np.nan)
            differ = int((from_name.notna() & (from_name != d["markers"]["POS"])).sum())
            if differ:
                notes.append(f"For {differ:,} markers the position written in the marker name differs from the {d['pos_source']}; the column was used.")
        return {**d, "notes": notes, "pheno": None, "traits": [], "layout": "marker rows"}
    except R.RecombError as first:
        try:
            from tools.gwas.backend import genotype_io as gio
            g = gio.read_gwas_file(path)
        except Exception:                                        # noqa: BLE001
            raise BsaError(str(first)) from None
    G = np.where(np.isnan(g.G), -1, np.round(g.G)).T
    G = np.where(np.isin(G, (0, 1, 2)), G, -1).astype(np.int8)
    markers = pd.DataFrame({"ID": [str(s) for s in g.snp_ids], "Chrom": np.asarray(g.chrom, dtype=int), "POS": np.asarray(g.pos, dtype=np.int64)})
    notes = list(g.marker_map.notes) + list(g.notes)
    return {"markers": markers, "samples": [str(s) for s in g.ids], "G": G, "layout": "one row per line (with the trait)",
            "pos_source": g.marker_map.pos_source, "notes": notes, "pheno": {"Pheno": dict(zip([str(s) for s in g.ids], g.y))},
            "traits": ["Pheno"], "id_col": None, "n_rows": int(G.shape[0])}


def inspect_inputs(geno_path, pheno_path=None, sheet=None, pheno_sheet=None) -> dict:
    d = read_genotypes(geno_path, sheet)
    m, G, s = d["markers"], d["G"], d["samples"]
    het = float((G == 1).sum() / max(int((G >= 0).sum()), 1))
    out = {"status": "success", "file": os.path.basename(str(geno_path)), "layout": d["layout"], "n_markers": int(len(m)),
           "n_samples": len(s), "chromosomes": sorted(m["Chrom"].unique().tolist()), "pos_source": d["pos_source"],
           "het_pct": round(het * 100, 2), "missing_pct": round(float((G == -1).mean() * 100), 2),
           "suggest_population": suggest_population(G), "class_pct": [round(f * 100, 1) for f in class_fractions(G)], "notes": list(d["notes"]),
           "traits": list(d["traits"]), "pheno_from": "genotype file" if d["pheno"] else "", "n_matched": 0}
    prev = pd.DataFrame({"Marker": m["ID"][:8], "Chromosome": m["Chrom"][:8], "Position": m["POS"][:8]})
    for j, name in enumerate(s[:6]):
        prev[name] = G[:8, j]
    out["preview_columns"], out["preview_rows"] = list(prev.columns), prev.to_dict("records")
    if pheno_path:
        p = read_phenotypes(pheno_path, pheno_sheet)
        out["traits"], out["pheno_from"] = p["traits"], "phenotype file"
        out["n_matched"] = int(p["table"][p["id_col"]].isin(set(s)).sum())
        if out["n_matched"] < 10:
            raise BsaError(f"Only {out['n_matched']} sample names of the phenotype file match the genotype file "
                           f"(e.g. genotype: '{s[0]}', phenotype: '{p['table'][p['id_col']].iloc[0]}'). Check the names.")
    elif d["pheno"]:
        out["n_matched"] = len(s)
    return out


# ───────────────────────────── statistics ─────────────────────────────

def _g_table(O: np.ndarray) -> np.ndarray:
    """G = sum O ln(O/E) of a pool x class table, vectorised.  O: (2, k, ...) counts of pool x class; -> (...)."""
    rows = O.sum(axis=1, keepdims=True)                      # (2, 1, ...)
    cols = O.sum(axis=0, keepdims=True)                      # (1, k, ...)
    tot = O.sum(axis=(0, 1), keepdims=True)                  # (1, 1, ...)
    with np.errstate(divide="ignore", invalid="ignore"):
        E = rows * cols / tot
        t = np.where((O > 0) & (E > 0), O * np.log(np.where(E > 0, O / np.where(E > 0, E, 1), 1)), 0.0)
    return t.sum(axis=(0, 1))


def pool_stats(Z0, Z1, Z2, hi, lo, diploid: bool, min_calls: int):
    """Statistics of every marker for pools given as 0/1 sample indicator matrices hi, lo (samples x P).

    Z0, Z1, Z2: (markers x samples) indicator (float) of each genotype class among the called genotypes.
    Returns dict of (markers x P) arrays: Gw, G1, G2, dAF, plus counts needed for reporting (only for P == 1)."""
    n = {}
    for name, Z in (("0", Z0), ("1", Z1), ("2", Z2)):
        n["h" + name] = Z @ hi
        n["l" + name] = Z @ lo
    hn, ln_ = n["h0"] + n["h1"] + n["h2"], n["l0"] + n["l1"] + n["l2"]
    ok = (hn >= min_calls) & (ln_ >= min_calls)
    if diploid:
        a1 = [2 * n["h0"] + n["h1"], 2 * n["h2"] + n["h1"]]
        b1 = [2 * n["l0"] + n["l1"], 2 * n["l2"] + n["l1"]]
    else:
        a1, b1 = [n["h0"], n["h2"]], [n["l0"], n["l2"]]
    Q = np.stack([np.stack(a1), np.stack(b1)])               # (pool, allele, markers, P)
    Gn = np.stack([np.stack([n["h0"], n["h1"], n["h2"]]), np.stack([n["l0"], n["l1"], n["l2"]])])
    g1 = _g_table(Q)
    g2 = _g_table(Gn)
    s = g1 + g2
    with np.errstate(divide="ignore", invalid="ignore"):
        gw = np.where(s > 0, (g1 * g1 + g2 * g2) / s, 0.0)
        tot_h, tot_l = Q[0].sum(axis=0), Q[1].sum(axis=0)
        sdaf = np.where((tot_h > 0) & (tot_l > 0), Q[0][1] / np.where(tot_h > 0, tot_h, 1) - Q[1][1] / np.where(tot_l > 0, tot_l, 1), 0.0)
    daf = np.abs(sdaf)
    gw, g1, g2, daf, sdaf = (np.where(ok, x, 0.0) for x in (gw, g1, g2, daf, sdaf))
    return {"Gw": gw, "G1": g1, "G2": g2, "dAF": daf, "sdAF": sdaf, "n": n, "ok": ok}


# ───────────────────────────── analysis ─────────────────────────────

def _as_bool(v, default=False):
    if v is None:
        return default
    return str(v).strip().lower() in ("1", "true", "yes", "on")


def build_pools(y: np.ndarray, frac: float, n_each: int | None = None):
    """Indices of the low and high pool: the lowest and highest n lines by trait value (stable order for ties)."""
    n = len(y)
    k = int(n_each) if n_each else int(round(frac * n))
    k = max(2, min(k, n // 2))
    order = np.argsort(y, kind="stable")
    return order[:k], order[-k:], k


def analyse(config: dict, log=print) -> dict:
    geno = read_genotypes(config["input_file"], config.get("sheet") or None)
    samples, G, markers = geno["samples"], geno["G"], geno["markers"]
    notes = list(geno["notes"])
    trait = config.get("trait") or ""
    # ----- phenotype -----
    if config.get("pheno_file"):
        p = read_phenotypes(config["pheno_file"], config.get("pheno_sheet") or None)
        if trait not in p["traits"]:
            if not trait and p["traits"]:
                trait = p["traits"][0]
            else:
                raise BsaError(f"The trait '{trait}' is not a numeric column of the phenotype file.")
        t = p["table"].drop_duplicates(p["id_col"]).set_index(p["id_col"])[trait]
        y_all = np.array([t.get(s, np.nan) for s in samples], dtype=float)
    elif geno["pheno"]:
        trait = trait or geno["traits"][0]
        y_all = np.array([geno["pheno"][trait].get(s, np.nan) for s in samples], dtype=float)
    else:
        raise BsaError("Choose a phenotype file: this genotype file does not contain the trait.")
    have = np.isfinite(y_all)
    if int(have.sum()) < 20:
        raise BsaError(f"Only {int(have.sum())} lines have both a genotype and a trait value; at least 20 are needed.")
    if int(have.sum()) < len(samples):
        notes.append(f"{len(samples) - int(have.sum())} of {len(samples)} lines have no trait value and were left out.")
    keep_idx = np.flatnonzero(have)
    names = [samples[i] for i in keep_idx]
    y = y_all[keep_idx]
    G = G[:, keep_idx]

    # ----- population type -----
    het = float((G == 1).sum() / max(int((G >= 0).sum()), 1))
    popkey, pnotes = resolve_population(config.get("population"), G)
    pop = POPULATIONS[popkey]
    notes += pnotes + population_check(G, popkey)
    diploid = pop["diploid"]
    if not diploid:
        G = np.where(G == 1, -1, G).astype(np.int8)
        if het > 0.0:
            notes.append(f"{het * 100:.1f} % of the calls are heterozygous; they are treated as missing because the population is {pop['short']} "
                         "(one allele per line). Choose F2 / IMF2 or BC1 to use them.")
    parent = int(config.get("snp_index_parent", 2) or 2)
    if parent not in (1, 2):
        raise BsaError("The SNP-index can show the allele of parent 1 or of parent 2.")
    sgn = 1.0 if parent == 2 else -1.0          # allele frequencies are kept for parent 2 (code 2); parent 1 is their complement
    # ----- markers -----
    min_calls = int(config.get("min_calls", 5))
    pf = float(config.get("pool_fraction", 0.2))
    if not 0.05 <= pf <= 0.5:
        raise BsaError("The pool fraction must be between 5 % and 50 % of the lines.")
    lo_i, hi_i, k = build_pools(y, pf, config.get("pool_size"))
    if k < min_calls:
        raise BsaError(f"Each pool would hold only {k} lines; at least {min_calls} are needed (raise the pool fraction or add lines).")
    called = G >= 0
    n_called = called.sum(axis=1)
    maf_ok = np.zeros(len(markers), dtype=bool)
    cnt = np.stack([(G == 0).sum(axis=1), (G == 1).sum(axis=1) if diploid else np.zeros(len(markers)), (G == 2).sum(axis=1)])
    allele2 = (2 * cnt[2] + cnt[1]) if diploid else cnt[2]
    allele_n = (2 * n_called) if diploid else n_called
    with np.errstate(divide="ignore", invalid="ignore"):
        f2 = np.where(allele_n > 0, allele2 / allele_n, 0.0)
    maf = np.minimum(f2, 1 - f2)
    min_maf = float(config.get("min_maf", 0.05))
    miss_max = float(config.get("max_missing", 0.5))
    maf_ok = (maf >= min_maf) & (1 - n_called / G.shape[1] <= miss_max)
    if int(maf_ok.sum()) < 50:
        raise BsaError(f"Only {int(maf_ok.sum())} markers pass the filters (minor allele frequency >= {min_maf:g}, at most {miss_max * 100:g} % missing).")
    markers = markers.loc[maf_ok].reset_index(drop=True)
    G = G[maf_ok]
    m = len(markers)
    chrom, pos = markers["Chrom"].to_numpy(), markers["POS"].to_numpy()
    Z = [(G == c).astype(np.float64) for c in (0, 1, 2)]

    # ----- observed -----
    n_samp = G.shape[1]
    hi = np.zeros((n_samp, 1))
    lo = np.zeros((n_samp, 1))
    hi[hi_i, 0] = 1
    lo[lo_i, 0] = 1
    obs = pool_stats(Z[0], Z[1], Z[2], hi, lo, diploid, min_calls)
    scored = obs["ok"].reshape(-1)
    raw = {s: np.where(scored, obs[s].reshape(-1), np.nan) for s in STATS}      # unscored: empty, never a zero
    n_unscored = int((~scored).sum())
    if n_unscored:
        notes.append(f"{n_unscored:,} of {m:,} markers had fewer than {min_calls} called lines in a pool. They get no statistic "
                     "(they are not counted as zero) and do not enter any window.")

    # ----- windows -----
    sizes = _window_sizes(config)
    primary_kb = float(config.get("primary_window_kb", 500))
    if primary_kb not in sizes:
        primary_kb = sizes[len(sizes) // 2]
    min_cfg = int(config.get("min_window_markers", 0) or 0)           # 0 / empty = automatic
    if min_cfg < 0:
        raise BsaError("The minimum number of scored markers per window must not be negative.")
    wins = {}
    for kb in sizes:
        A, tab = build_windows(chrom, pos, kb * 1000.0)
        # automatic: 3 scored markers per window, fewer when the map is sparse (median markers per occupied window), never below 1 -
        # a fixed 3 would leave almost every window untested on a 300-marker map
        auto = int(min(3, max(1, np.floor(np.median(np.asarray(A.sum(axis=1)).reshape(-1))))))
        wins[window_label(kb)] = {"kb": kb, "A": A, "table": tab, "min_markers": min_cfg or auto, "min_auto": not min_cfg}
    min_wm = {lab: w["min_markers"] for lab, w in wins.items()}
    primary = window_label(primary_kb)

    # ----- sliding SNP-index grid (the primary window size, moved in steps of a fifth of it) -----
    # a sparse or clustered map (a few markers per Mb) would give a jagged, broken line: widen the window until it holds enough markers
    snp_kb = primary_kb
    for cand in [primary_kb] + [k for k in SLIDE_STEPS_KB if k > primary_kb]:
        snp_kb = cand
        cnt = np.asarray(build_sliding(chrom, pos, cand * 1000.0, cand * 1000.0 / SLIDE_STEPS)[0].sum(axis=1)).reshape(-1)
        if cnt.size == 0 or (cnt >= SLIDE_MIN_MARKERS).mean() >= SLIDE_FILL or cand >= SLIDE_MAX_KB:
            break
    snp_auto = snp_kb != primary_kb
    slide_S, slide_tab = build_sliding(chrom, pos, snp_kb * 1000.0, snp_kb * 1000.0 / SLIDE_STEPS)
    slide_min = max(wins[primary]["min_markers"], SLIDE_MIN_MARKERS if snp_auto else 1)
    if snp_auto:
        notes.append(f"The map is sparse ({len(chrom) / max(1.0, (pos.max() - pos.min()) / 1e6):.1f} markers per Mb on average), so the SNP-index tracks use a sliding window of "
                     f"{window_label(snp_kb)} instead of {window_label(primary_kb)} to hold enough markers per window (at least {slide_min}). The Gw windows are not affected.")

    # ----- permutations: random low / high pools of the same size, the largest window statistic of the genome per window size -----
    n_perm = int(config.get("permutations", 1000))
    if n_perm < 0 or n_perm > MAX_PERMUTATIONS:
        raise BsaError(f"The number of permutations must be between 0 and {MAX_PERMUTATIONS:,}.")
    n_null = min(n_perm, MAX_NULL_PERMUTATIONS)
    null_delta = np.full((len(slide_tab), n_null), np.nan, dtype=np.float32) if n_null >= 20 else None
    seed = int(config.get("seed", 2022))
    rng = np.random.default_rng(seed)
    for w in wins.values():
        w["maxima"] = {s: np.full(n_perm, np.nan) for s in STATS}
    batch = max(1, min(n_perm, int(2.5e7 // max(m * 3, 1)), 200))
    done = 0
    while done < n_perm:
        b = min(batch, n_perm - done)
        H = np.zeros((n_samp, b))
        L = np.zeros((n_samp, b))
        for j in range(b):
            perm = rng.permutation(n_samp)
            L[perm[:k], j] = 1
            H[perm[-k:], j] = 1
        st = pool_stats(Z[0], Z[1], Z[2], H, L, diploid, min_calls)
        if null_delta is not None and done < n_null:
            nb = min(b, n_null - done)
            null_delta[:, done:done + nb] = (sgn * sliding_mean(slide_S, st["sdAF"][:, :nb], st["ok"][:, :nb], slide_min)).astype(np.float32)
        for w in wins.values():
            means, valid = window_means(w["A"], st, w["min_markers"])
            for s in STATS:
                mx = np.where(valid, means[s], -np.inf).max(axis=0)
                w["maxima"][s][done:done + b] = np.where(np.isfinite(mx), mx, np.nan)
        done += b
        log(f"Permutations {done:,} of {n_perm:,}")
    alpha_levels = (0.10, 0.05, 0.01)
    alpha = float(config.get("alpha", 0.05))
    levels = tuple(sorted(set(alpha_levels + (alpha,)), reverse=True))

    # ----- window statistics, thresholds, adjusted p-values (observed pools) -----
    n_o = {k_: v.reshape(-1) for k_, v in obs["n"].items()}
    tot_h = n_o["h0"] + n_o["h1"] + n_o["h2"]
    tot_l = n_o["l0"] + n_o["l1"] + n_o["l2"]
    with np.errstate(divide="ignore", invalid="ignore"):
        def frac_(a, b):
            return np.where(b > 0, a / np.where(b > 0, b, 1), np.nan)
        af_h = frac_((2 * n_o["h2"] + n_o["h1"]) if diploid else n_o["h2"], (2 * tot_h) if diploid else tot_h)
        af_l = frac_((2 * n_o["l2"] + n_o["l1"]) if diploid else n_o["l2"], (2 * tot_l) if diploid else tot_l)
        het_h, het_l = frac_(n_o["h1"], tot_h), frac_(n_o["l1"], tot_l)
    mk = pd.DataFrame({
        "Marker": markers["ID"], "Chromosome": chrom, "Position": pos,
        "Low_0": n_o["l0"].astype(int), "Low_1": n_o["l1"].astype(int), "Low_2": n_o["l2"].astype(int),
        "High_0": n_o["h0"].astype(int), "High_1": n_o["h1"].astype(int), "High_2": n_o["h2"].astype(int),
        "AF2_low": af_l, "AF2_high": af_h, "dAF": af_h - af_l, "Het_low": het_l, "Het_high": het_h,
        "SNP_index_high": af_h if parent == 2 else 1 - af_h, "SNP_index_low": af_l if parent == 2 else 1 - af_l, "Delta_SNP_index": sgn * (af_h - af_l),
        "Scored": scored, "G1": raw["G1"], "G2": np.where(scored, obs["G2"].reshape(-1), np.nan), "Gw": raw["Gw"], "abs_dAF": raw["dAF"],
    })
    merge_kb = float(config.get("merge_gap_kb", 0))
    if merge_kb < 0:
        raise BsaError("The distance for merging neighbouring windows must not be negative.")
    all_regions = []
    for lab, w in wins.items():
        A = w["A"]
        means, valid = window_means(A, obs, w["min_markers"])
        tab = w["table"]
        tab["Markers"] = np.asarray(A.sum(axis=1)).reshape(-1).astype(int)
        tab["Scored_markers"] = np.asarray(A @ scored.astype(float)).reshape(-1).astype(int)
        w["thr"] = {s: {a_: (float(np.nanquantile(w["maxima"][s], 1 - a_)) if n_perm >= 20 and np.isfinite(w["maxima"][s]).any() else float("nan"))
                        for a_ in levels} for s in STATS}
        w["thr_use"] = {s: w["thr"][s][alpha] for s in STATS}
        for s in STATS:
            v = means[s].reshape(-1)
            col = {"Gw": "Gw", "G1": "G1", "dAF": "abs_dAF"}[s]
            tab[col] = v
            if n_perm >= 20:
                srt = np.sort(w["maxima"][s][np.isfinite(w["maxima"][s])])
                exceed = len(srt) - np.searchsorted(srt, np.nan_to_num(v), side="left")           # #{T_b >= value}
                tab["P_adj_" + s] = np.where(np.isfinite(v), (1 + exceed) / (n_perm + 1), np.nan)   # (1 + #{T_b >= S}) / (B + 1)
            else:
                tab["P_adj_" + s] = np.nan
            # significant = adjusted p-value <= alpha, the same rule the reported p-values follow (a window sitting exactly on the interpolated
            # quantile would otherwise be called significant with P_adj just above alpha)
            tab["Sig_" + s] = (np.isfinite(v) & (np.nan_to_num(tab["P_adj_" + s].to_numpy(dtype=float), nan=1.0) <= alpha + 1e-12)) if n_perm >= 20 else False
        n_skip = int((~np.isfinite(tab["Gw"].to_numpy())).sum())
        w["n_untested"] = n_skip
        if n_skip:
            notes.append(f"{lab}: {n_skip:,} of {len(tab):,} windows hold fewer than {w['min_markers']} scored markers and were not tested"
                         + (" (marker density is too low for this window size; use a larger window or fewer required markers)." if n_skip > 0.5 * len(tab) else "."))
        # per-marker view of the window it lies in
        wid = _window_index(chrom, pos, w["kb"] * 1000.0, tab)
        mk[f"Gw_win_{lab}"] = tab["Gw"].to_numpy()[wid]
        mk[f"P_adj_win_{lab}"] = tab["P_adj_Gw"].to_numpy()[wid]
        mk[f"Sig_win_{lab}"] = tab["Sig_Gw"].to_numpy()[wid].astype(bool)
        w["regions"] = window_regions(mk, tab, lab, w["kb"], merge_kb * 1000.0)
        w["near"] = near_regions(mk, tab, lab, w["kb"], merge_kb * 1000.0, w["thr_use"]["Gw"]) if not len(w["regions"]) else pd.DataFrame()
        all_regions.append(w["regions"])
    mk["P_adj_Gw"] = mk[f"P_adj_win_{primary}"]
    regions = pd.concat(all_regions, ignore_index=True) if all_regions else pd.DataFrame()
    # ----- SNP-index tracks: high pool, low pool and the difference, smoothed; permutation intervals of the difference -----
    track = slide_tab.copy()
    af_hi = np.where(scored, mk["SNP_index_high"].to_numpy(), np.nan)
    af_lo = np.where(scored, mk["SNP_index_low"].to_numpy(), np.nan)
    track["Markers"] = np.asarray(slide_S @ scored.astype(float)).reshape(-1).astype(int)
    track["SNP_index_high"] = sliding_mean(slide_S, af_hi, scored, slide_min).reshape(-1)
    track["SNP_index_low"] = sliding_mean(slide_S, af_lo, scored, slide_min).reshape(-1)
    track["Delta_SNP_index"] = track["SNP_index_high"] - track["SNP_index_low"]
    for lab_, q in (("CI95", 2.5), ("CI99", 0.5)):
        if null_delta is not None:
            with np.errstate(all="ignore"), warnings.catch_warnings():
                warnings.simplefilter("ignore")
                lo_q, hi_q = np.nanpercentile(null_delta, [q, 100 - q], axis=1)
            track[f"{lab_}_lower"], track[f"{lab_}_upper"] = lo_q.astype(float), hi_q.astype(float)
        else:
            track[f"{lab_}_lower"] = track[f"{lab_}_upper"] = np.nan
    snp_regs = snp_index_regions(track, mk, snp_kb * 1000.0, merge_kb * 1000.0, parent, wins[primary]["regions"])
    near_all = [w["near"] for w in wins.values() if len(w["near"])]
    near = pd.concat(near_all, ignore_index=True) if near_all else pd.DataFrame()
    comparison, loci = compare_windows(wins, sizes)
    pools = pd.DataFrame({"Sample": names, "Trait": y, "Pool": "-"})
    pools.loc[lo_i, "Pool"] = "Low"
    pools.loc[hi_i, "Pool"] = "High"
    pools = pools.sort_values("Trait").reset_index(drop=True)
    notes = list(dict.fromkeys(notes))
    pw = wins[primary]
    spacing = [np.diff(np.sort(pos[chrom == c])) for c in np.unique(chrom) if (chrom == c).sum() > 1]
    return {
        "markers": mk, "regions": regions, "near": near, "pools": pools, "windows": wins, "primary": primary, "comparison": comparison, "loci": loci,
        "thresholds": pw["thr"], "thr_use": pw["thr_use"], "alpha": alpha, "maxima": pw["maxima"], "sizes": sizes, "min_window_markers": min_wm, "min_auto": not min_cfg,
        "n_perm": n_perm, "diploid": diploid, "population": popkey, "population_label": pop["label"], "merge_gap_kb": merge_kb,
        "median_spacing_kb": round(float(np.median(np.concatenate(spacing))) / 1000.0, 2) if spacing else 0.0,
        "n_unscored": n_unscored, "pool_size": k, "n_lines": len(names), "trait": trait, "y": y, "lo_i": lo_i, "hi_i": hi_i,
        "notes": notes, "layout": geno["layout"], "pos_source": geno["pos_source"],
        "n_markers_in": int(len(geno["markers"])), "n_markers_used": m, "het_pct": het * 100, "min_calls": min_calls, "pool_fraction": pf,
        "min_maf": min_maf, "max_missing": miss_max, "seed": seed, "chroms": sorted(set(chrom.tolist())),
        "pool_cuts": (float(y[lo_i].max()), float(y[hi_i].min())),
        "snp_track": track, "snp_regions": snp_regs, "snp_parent": parent, "snp_window_kb": snp_kb, "snp_window_auto": snp_auto, "snp_null_perms": n_null if null_delta is not None else 0,
    }


# ───────────────────────────── sliding SNP-index tracks ─────────────────────────────

def build_sliding(chrom: np.ndarray, pos: np.ndarray, size_bp: float, step_bp: float):
    """Sliding windows of `size_bp` centred every `step_bp` along every chromosome.
    -> (sparse grid x markers 0/1 membership, DataFrame[Chromosome, Position] of the window centres)."""
    rows, cols, cc, xx = [], [], [], []
    g = 0
    for c in np.unique(chrom):
        idx = np.flatnonzero(chrom == c)
        order = idx[np.argsort(pos[idx], kind="stable")]
        ps = pos[order].astype(float)
        centres = np.arange(ps.min(), ps.max() + step_bp, step_bp)
        lo = np.searchsorted(ps, centres - size_bp / 2.0, side="left")
        hi = np.searchsorted(ps, centres + size_bp / 2.0, side="right")
        for i in range(len(centres)):
            if hi[i] > lo[i]:
                members = order[lo[i]:hi[i]]
                rows.append(np.full(len(members), g + i))
                cols.append(members)
        cc.append(np.full(len(centres), int(c)))
        xx.append(centres)
        g += len(centres)
    r = np.concatenate(rows) if rows else np.zeros(0, dtype=int)
    k = np.concatenate(cols) if cols else np.zeros(0, dtype=int)
    S = sparse.csr_matrix((np.ones(len(r)), (r, k)), shape=(g, len(chrom)))
    return S, pd.DataFrame({"Chromosome": np.concatenate(cc).astype(int), "Position": np.concatenate(xx).astype(np.int64)})


def sliding_mean(S: sparse.csr_matrix, values: np.ndarray, ok: np.ndarray, min_markers: int) -> np.ndarray:
    """Mean of `values` over the scored markers of every sliding window (NaN when fewer than `min_markers`); values: markers x P."""
    okf = ok.astype(float)
    if okf.ndim == 1:
        okf = okf[:, None]
    cnt = np.asarray(S @ okf)
    v = np.nan_to_num(np.asarray(values, dtype=float).reshape(len(okf), -1))
    return np.where(cnt >= min_markers, np.asarray(S @ (v * okf)) / np.where(cnt > 0, cnt, 1.0), np.nan)


def snp_index_regions(track: pd.DataFrame, mk: pd.DataFrame, size_bp: float, merge_gap_bp: float, parent: int, gw_regions: pd.DataFrame) -> pd.DataFrame:
    """Delta SNP-index peak regions: runs of sliding windows whose smoothed Delta lies outside the 99 % interval of the permutations.
    A region ends where the run is broken by more than `merge_gap_bp`; its peak is the window with the largest |Delta|."""
    cols = ["Region", "Chromosome", "Left_marker", "Peak_marker", "Right_marker", "Start_bp", "End_bp", "Left_bp", "Right_bp", "Length_kb", "Peak_bp",
            "Delta_SNP_index_peak", "SNP_index_high_peak", "SNP_index_low_peak", "CI99_at_peak", "Windows_outside_CI99", "Markers", "Direction", "Also_Gw_region"]
    out = []
    if not {"CI99_lower", "CI99_upper"} <= set(track.columns) or not track["CI99_upper"].notna().any():
        return pd.DataFrame(columns=cols)
    d = track["Delta_SNP_index"].to_numpy()
    sig = np.isfinite(d) & ((d > track["CI99_upper"].to_numpy()) | (d < track["CI99_lower"].to_numpy()))
    mk_sc = mk[mk["Scored"]]
    for c, g in track.groupby("Chromosome", sort=True):
        ii = np.flatnonzero(sig[g.index.to_numpy()])
        if not len(ii):
            continue
        x = g["Position"].to_numpy()[ii]
        cut = np.flatnonzero(np.diff(x) > max(merge_gap_bp, size_bp)) + 1      # runs closer than one window (or the merge distance) are one region
        for grp in np.split(np.arange(len(ii)), cut):
            sel = g.iloc[ii[grp]]
            pk = sel.iloc[int(np.argmax(np.abs(sel["Delta_SNP_index"].to_numpy())))]
            lo_b, hi_b = float(sel["Position"].min()) - size_bp / 2.0, float(sel["Position"].max()) + size_bp / 2.0
            inreg = mk_sc[(mk_sc["Chromosome"] == c) & (mk_sc["Position"] >= lo_b) & (mk_sc["Position"] <= hi_b)]
            if not len(inreg):
                continue
            pkm = inreg.iloc[int(np.argmin(np.abs(inreg["Position"].to_numpy() - float(pk["Position"]))))]
            left, right = inreg.loc[inreg["Position"].idxmin()], inreg.loc[inreg["Position"].idxmax()]
            gws = gw_regions[(gw_regions["Chromosome"] == c) & (gw_regions["Start_bp"] <= right.Position) & (gw_regions["End_bp"] >= left.Position)] if len(gw_regions) else gw_regions
            out.append({"Chromosome": int(c), "Left_marker": left.Marker, "Peak_marker": pkm.Marker, "Right_marker": right.Marker,
                        "Start_bp": int(max(lo_b, 0)), "End_bp": int(hi_b), "Left_bp": int(left.Position), "Right_bp": int(right.Position),
                        "Length_kb": round((right.Position - left.Position) / 1000.0, 1), "Peak_bp": int(pkm.Position),
                        "Delta_SNP_index_peak": round(float(pk["Delta_SNP_index"]), 4), "SNP_index_high_peak": round(float(pk["SNP_index_high"]), 4),
                        "SNP_index_low_peak": round(float(pk["SNP_index_low"]), 4),
                        "CI99_at_peak": round(float(pk["CI99_upper"] if pk["Delta_SNP_index"] > 0 else pk["CI99_lower"]), 4),
                        "Windows_outside_CI99": int(len(grp)), "Markers": int(len(inreg)),
                        "Direction": (f"high pool carries more of the parent {parent} allele" if pk["Delta_SNP_index"] > 0
                                      else f"low pool carries more of the parent {parent} allele"),
                        "Also_Gw_region": bool(len(gws))})
    df = pd.DataFrame(out, columns=[x for x in cols if x != "Region"])
    if len(df):
        df = df.reindex(df["Delta_SNP_index_peak"].abs().sort_values(ascending=False).index).reset_index(drop=True)
    df.insert(0, "Region", [f"D{i + 1}" for i in range(len(df))])
    return df


# ───────────────────────────── windows ─────────────────────────────

def window_label(kb: float) -> str:
    """250 -> '250kb', 1000 -> '1Mb'."""
    return f"{kb / 1000:g}Mb" if kb >= 1000 and float(kb) % 1000 == 0 else f"{kb:g}kb"


def _window_sizes(config: dict) -> list:
    raw = config.get("window_sizes_kb") or WINDOW_SIZES_KB
    try:
        sizes = sorted({float(x) for x in raw})
    except (TypeError, ValueError):
        raise BsaError("The window sizes must be numbers (kb).") from None
    if not sizes or sizes[0] <= 0 or len(sizes) > 6:
        raise BsaError("Give between 1 and 6 window sizes, each larger than 0 kb.")
    return sizes


def build_windows(chrom: np.ndarray, pos: np.ndarray, size_bp: float):
    """Non-overlapping physical windows of `size_bp` on every chromosome (window i covers [i*size, (i+1)*size)).
    -> (sparse windows x markers 0/1 membership, DataFrame[Chromosome, Window_start, Window_end] of the windows that hold a marker)."""
    b = np.floor(pos.astype(float) / size_bp).astype(np.int64)
    df = pd.DataFrame({"c": chrom, "b": b})
    gid = df.groupby(["c", "b"], sort=True).ngroup().to_numpy()
    nw = int(gid.max()) + 1
    A = sparse.csr_matrix((np.ones(len(gid)), (gid, np.arange(len(gid)))), shape=(nw, len(gid)))
    u = df.drop_duplicates().sort_values(["c", "b"])
    tab = pd.DataFrame({"Chromosome": u["c"].to_numpy().astype(int), "Window_start": (u["b"].to_numpy() * size_bp).astype(np.int64),
                        "Window_end": ((u["b"].to_numpy() + 1) * size_bp).astype(np.int64)}).reset_index(drop=True)
    return A, tab


def _window_index(chrom, pos, size_bp, tab) -> np.ndarray:
    """Row of `tab` (the windows of build_windows) that each marker lies in."""
    b = np.floor(pos.astype(float) / size_bp).astype(np.int64) * int(size_bp)
    key = pd.MultiIndex.from_arrays([tab["Chromosome"].to_numpy(), tab["Window_start"].to_numpy()])
    return key.get_indexer(pd.MultiIndex.from_arrays([np.asarray(chrom).astype(int), b]))


def window_means(A: sparse.csr_matrix, st: dict, min_markers: int):
    """Window statistic = mean of the marker statistic over the SCORED markers of the window (unscored markers are not zeros).
    A window with fewer than `min_markers` scored markers is not tested (NaN). st: pool_stats() result (markers x P).
    -> ({stat: windows x P}, valid mask windows x P)."""
    okf = st["ok"].astype(float)
    cnt = np.asarray(A @ okf)
    valid = cnt >= min_markers
    den = np.where(cnt > 0, cnt, 1.0)
    return {s: np.where(valid, np.asarray(A @ (st[s] * okf)) / den, np.nan) for s in STATS}, valid


def window_regions(mk: pd.DataFrame, tab: pd.DataFrame, label: str, size_kb: float, merge_gap_bp: float = 0.0) -> pd.DataFrame:
    """QTL regions of one window size: significant windows (window Gw above the permutation threshold of that window size) of a chromosome
    are one region when the physical gap between them is at most `merge_gap_bp` (0 = touching windows only). The peak window has the highest
    window Gw; the peak marker is its marker with the highest Gw."""
    cols = ["Region", "Chromosome", "Left_marker", "Peak_marker", "Right_marker", "P_adj_Gw_at_peak", "LOD", "Delta_SNP_index", "Window", "Window_kb", "Start_bp", "End_bp",
            "Left_bp", "Right_bp", "Length_kb", "Windows_above_threshold", "Peak_window_start", "Peak_bp", "Peak_value", "Markers", "AF2_low_at_peak", "AF2_high_at_peak", "Delta_SNP_index_window", "Het_low_at_peak",
            "Het_high_at_peak", "Pattern", "Also_significant_by_G1", "Also_significant_by_dAF"]
    out = []
    for c, g in tab.groupby("Chromosome", sort=True):
        sig = g[g["Sig_Gw"].astype(bool)].sort_values("Window_start")
        if not len(sig):
            continue
        ws, we = sig["Window_start"].to_numpy(), sig["Window_end"].to_numpy()
        cut = np.flatnonzero(ws[1:] - we[:-1] > merge_gap_bp) + 1
        for grp in np.split(np.arange(len(sig)), cut):
            s_bp, e_bp = int(ws[grp[0]]), int(we[grp[-1]])
            sub = sig.iloc[grp]
            pkw = sub.iloc[int(np.argmax(sub["Gw"].to_numpy()))]
            inmk = mk[(mk["Chromosome"] == c) & (mk["Position"] >= pkw.Window_start) & (mk["Position"] < pkw.Window_end) & mk["Scored"]]
            pk = inmk.iloc[int(np.argmax(inmk["Gw"].to_numpy()))]
            span = g[(g["Window_start"] >= s_bp) & (g["Window_end"] <= e_bp)]
            inreg = mk[(mk["Chromosome"] == c) & (mk["Position"] >= s_bp) & (mk["Position"] < e_bp) & mk["Scored"]]
            left, right = inreg.loc[inreg["Position"].idxmin()], inreg.loc[inreg["Position"].idxmax()]
            out.append({"Window": label, "Window_kb": size_kb, "Chromosome": int(c), "Start_bp": s_bp, "End_bp": e_bp,
                        "Length_kb": round((e_bp - s_bp) / 1000.0, 1), "Windows_above_threshold": int(len(grp)),
                        "Peak_window_start": int(pkw.Window_start), "Peak_marker": pk.Marker, "Peak_bp": int(pk.Position),
                        "Left_marker": left.Marker, "Left_bp": int(left.Position), "Right_marker": right.Marker, "Right_bp": int(right.Position),
                        "Peak_value": float(pkw.Gw), "LOD": round(float(pkw.Gw) / np.log(10.0), 3), "Markers": int(span["Markers"].sum()), "AF2_low_at_peak": float(pk.AF2_low),
                        "AF2_high_at_peak": float(pk.AF2_high), "Delta_SNP_index": round(float(getattr(pk, "Delta_SNP_index", pk.dAF)), 4), "Delta_SNP_index_window": round(float(pkw.abs_dAF), 4), "Het_low_at_peak": float(pk.Het_low),
                        "Het_high_at_peak": float(pk.Het_high), "Pattern": _pattern(pk), "P_adj_Gw_at_peak": float(pkw.P_adj_Gw),
                        "Also_significant_by_G1": bool(sub["Sig_G1"].any()), "Also_significant_by_dAF": bool(sub["Sig_dAF"].any())})
    df = pd.DataFrame(out, columns=[x for x in cols if x != "Region"])
    if len(df):
        df = df.sort_values("Peak_value", ascending=False).reset_index(drop=True)
    df.insert(0, "Region", [f"{label}-R{i + 1}" for i in range(len(df))])
    return df


def near_regions(mk: pd.DataFrame, tab: pd.DataFrame, label: str, size_kb: float, merge_gap_bp: float, thr: float, n_windows: int = 8, n_keep: int = 5) -> pd.DataFrame:
    """For a window size with NO significant region: the candidates closest to significance, as regions of the same form. The `n_windows`
    windows with the highest Gw are grouped like significant windows (touching ones form one candidate) and the `n_keep` strongest
    candidates are returned, with how far they are from the threshold. They are NOT significant; they are what to follow up first."""
    ok = tab[tab["Gw"].notna()]
    if not len(ok):
        return pd.DataFrame()
    t2 = tab.copy()
    t2["Sig_Gw"] = False
    t2.loc[ok.sort_values("Gw", ascending=False).head(n_windows).index, "Sig_Gw"] = True
    reg = window_regions(mk, t2, label, size_kb, merge_gap_bp).head(n_keep).copy()
    if not len(reg):
        return reg
    reg["Region"] = [f"{label}-N{i + 1}" for i in range(len(reg))]
    reg.insert(reg.columns.get_loc("LOD") + 1, "Fraction_of_threshold", [round(float(v) / thr, 3) if np.isfinite(thr) and thr > 0 else np.nan for v in reg["Peak_value"]])
    return reg


def compare_windows(wins: dict, sizes: list):
    """-> (per window size summary, loci). A locus is a group of regions of the different window sizes that overlap on a chromosome; for
    every window size the table says whether it found the locus, with which region, how long, how significant, and where its peak lies."""
    labs = [window_label(k) for k in sizes]
    rows = []
    for lab in labs:
        w = wins[lab]
        reg, tab = w["regions"], w["table"]
        tested = int(np.isfinite(tab["Gw"].to_numpy()).sum())
        rows.append({"Window": lab, "Window_kb": w["kb"], "Windows_with_markers": int(len(tab)), "Windows_tested": tested,
                     "Median_markers_per_window": float(tab["Scored_markers"].median()) if len(tab) else np.nan,
                     "Gw_threshold": w["thr_use"]["Gw"], "Significant_windows": int(tab["Sig_Gw"].sum()), "Regions": int(len(reg)),
                     "Total_region_length_Mb": round(float(reg["Length_kb"].sum()) / 1000.0, 2) if len(reg) else 0.0,
                     "Median_region_length_kb": float(reg["Length_kb"].median()) if len(reg) else np.nan,
                     "Smallest_P_adj_Gw": float(np.nanmin(tab["P_adj_Gw"].to_numpy())) if np.isfinite(tab["P_adj_Gw"].to_numpy()).any() else np.nan})
    comp = pd.DataFrame(rows)
    allr = pd.concat([w["regions"] for w in wins.values()], ignore_index=True)
    cols = ["Locus", "Chromosome", "Start_bp", "End_bp", "Span_kb", "Peak_marker", "Peak_bp", "Found_by", "Windows_that_found_it", "Peak_spread_kb"]
    for lab in labs:
        cols += [f"Region_{lab}", f"Length_kb_{lab}", f"P_adj_{lab}", f"Peak_bp_{lab}"]
    if not len(allr):
        return comp, pd.DataFrame(columns=cols)
    loci = []
    for c, g in allr.groupby("Chromosome", sort=True):
        g = g.sort_values("Start_bp")
        cur, cur_end = [], -1
        groups = []
        for r in g.itertuples():
            if cur and r.Start_bp > cur_end:
                groups.append(cur)
                cur, cur_end = [], -1
            cur.append(r)
            cur_end = max(cur_end, r.End_bp)
        groups.append(cur)
        for grp in groups:
            best = min(grp, key=lambda r: (r.P_adj_Gw_at_peak if np.isfinite(r.P_adj_Gw_at_peak) else 9, -r.Peak_value))
            row = {"Chromosome": int(c), "Start_bp": int(min(r.Start_bp for r in grp)), "End_bp": int(max(r.End_bp for r in grp)),
                   "Peak_marker": best.Peak_marker, "Peak_bp": int(best.Peak_bp)}
            row["Span_kb"] = round((row["End_bp"] - row["Start_bp"]) / 1000.0, 1)
            found = []
            for lab in labs:
                rs = [r for r in grp if r.Window == lab]
                if rs:
                    r = max(rs, key=lambda x: x.Peak_value)
                    found.append(lab)
                    row.update({f"Region_{lab}": r.Region, f"Length_kb_{lab}": r.Length_kb, f"P_adj_{lab}": r.P_adj_Gw_at_peak, f"Peak_bp_{lab}": int(r.Peak_bp)})
                else:
                    row.update({f"Region_{lab}": "", f"Length_kb_{lab}": np.nan, f"P_adj_{lab}": np.nan, f"Peak_bp_{lab}": np.nan})
            pks = [row[f"Peak_bp_{lab}"] for lab in found]
            row["Found_by"] = ", ".join(found)
            row["Windows_that_found_it"] = f"{len(found)} of {len(labs)}"
            row["Peak_spread_kb"] = round((max(pks) - min(pks)) / 1000.0, 1)
            loci.append(row)
    df = pd.DataFrame(loci)
    df["_n"] = df["Found_by"].str.count(",") + 1
    df = df.sort_values(["_n", "Chromosome", "Start_bp"], ascending=[False, True, True]).drop(columns="_n").reset_index(drop=True)
    df.insert(0, "Locus", [f"L{i + 1}" for i in range(len(df))])
    return comp, df[cols]


def _pattern(pk) -> str:
    """Reads the genotype classes at the peak: an allele shift (additive / dominant) or a genotype shift without it (over-dominant)."""
    d = abs(float(pk.dAF)) if np.isfinite(pk.dAF) else 0.0
    dh = (pk.Het_high - pk.Het_low) if np.isfinite(pk.Het_high) and np.isfinite(pk.Het_low) else 0.0
    if d >= 0.15:
        return "allele frequency shift (additive or dominant)"
    if abs(dh) >= 0.10:
        return "heterozygotes enriched in the high pool (over-dominant)" if dh > 0 else "heterozygotes enriched in the low pool (over-dominant)"
    return "small shift"


# ───────────────────────────── outputs ─────────────────────────────

def _records(df: pd.DataFrame, n=MAX_ROWS):
    out = df.head(n).copy()
    for c in out.columns:
        if out[c].dtype.kind == "f":
            out[c] = out[c].round(4)
    return json_safe(out.astype(object).where(pd.notna(out), None).to_dict("records"))


def _thr(v):
    return round(v, 4) if np.isfinite(v) else "(no permutations)"


def summary_table(a: dict) -> pd.DataFrame:
    wins, prim, loci = a["windows"], a["primary"], a["loci"]
    pr = wins[prim]["regions"]
    rows = [("Lines analysed", a["n_lines"]), ("Trait", a["trait"]), ("Population", a["population_label"]),
            ("Lines per pool", a["pool_size"]), ("Pool fraction", f"{a['pool_fraction'] * 100:g} %"),
            ("Low pool: trait up to", round(a["pool_cuts"][0], 4)), ("High pool: trait from", round(a["pool_cuts"][1], 4)),
            ("Markers analysed", a["n_markers_used"]), ("Window sizes compared", ", ".join(window_label(k) for k in a["sizes"])),
            ("Primary window (figures, headline numbers)", prim)]
    for lab, w in wins.items():
        rows.append((f"Significant regions, {lab} windows (alpha {a['alpha']:g})", int(len(w["regions"]))))
    n_all = int((loci["Windows_that_found_it"].str.split(" of ").str[0].astype(int) == len(a["sizes"])).sum()) if len(loci) else 0
    rows += [("Loci (regions of the different window sizes that overlap)", int(len(loci))), ("Loci found by every window size", n_all),
             ("Gw threshold", _thr(a["thr_use"]["Gw"])), ("G1 threshold", _thr(a["thr_use"]["G1"])), ("|dAF| threshold", _thr(a["thr_use"]["dAF"])),
             ("Permutations", a["n_perm"])]
    rows.append((f"Delta SNP-index peak regions (outside the 99 % interval, {window_label(a['snp_window_kb'])} sliding window)", int(len(a["snp_regions"]))))
    if len(pr):
        rows.append((f"Regions that G1 alone also finds ({prim})", int(pr["Also_significant_by_G1"].sum())))
        rows.append((f"Regions found only through genotype counts ({prim})", int((~pr["Also_significant_by_G1"]).sum())))
    return pd.DataFrame(rows, columns=["Metric", "Value"])


def settings_table(a: dict, config: dict) -> pd.DataFrame:
    rows = [("Genotype file", os.path.basename(str(config["input_file"]))), ("Layout", a["layout"]), ("Positions read from", a["pos_source"]),
            ("Phenotype file", os.path.basename(str(config.get("pheno_file") or "")) or "(trait in the genotype file)"), ("Trait", a["trait"]),
            ("Pool fraction (each side)", f"{a['pool_fraction'] * 100:g} %"), ("Minimum called lines per pool", a["min_calls"]),
            ("Minor allele frequency filter", a["min_maf"]), ("Maximum missing per marker", f"{a['max_missing'] * 100:g} %"),
            ("Window sizes", ", ".join(window_label(k) for k in a["sizes"]) + " (non-overlapping; window statistic = mean over scored markers)"),
            ("Minimum scored markers per window", ", ".join(f"{lab}: {v}" for lab, v in a["min_window_markers"].items()) + (" (automatic)" if a["min_auto"] else "")), ("Primary window", a["primary"]),
            ("Permutations", a["n_perm"]), ("Significance level", a["alpha"]),
            ("SNP-index", f"frequency of the allele of parent {a['snp_parent']} (code {a['snp_parent']}) in a pool; every marker has one value for the high pool and one for the low pool"),
            ("Delta SNP-index", f"SNP-index of the high pool minus SNP-index of the low pool (positive: the allele of parent {a['snp_parent']} is more frequent in the high pool); "
                                "Delta_SNP_index of a region is the value at its peak marker, Delta_SNP_index_window the mean |difference| of the peak window"),
            ("SNP-index tracks", f"sliding window of {window_label(a['snp_window_kb'])}{' (widened automatically because the map is sparse)' if a.get('snp_window_auto') else ''} moved in steps of {window_label(a['snp_window_kb'] / SLIDE_STEPS)}; the 95 % and 99 % intervals are the "
                                 f"2.5 / 97.5 % and 0.5 / 99.5 % quantiles of the smoothed Delta of {a['snp_null_perms']:,} random low / high pools at every position (pointwise)"
                                 if a["snp_null_perms"] else "no permutations: no intervals"),
            ("LOD of a region", "Gw of the peak window / ln(10): G here is the log likelihood ratio (no factor 2), so G / ln(10) is a LOD score"),
            ("Population model", a["population_label"]), ("Merge significant windows closer than", f"{a['merge_gap_kb']:g} kb"),
            ("Median distance between neighbouring markers", f"{a['median_spacing_kb']:g} kb"),
            ("Markers without a statistic (too few calls in a pool)", a["n_unscored"]), ("Random seed", a["seed"])]
    return pd.DataFrame(rows, columns=["Setting", "Value"])


def thresholds_table(a: dict) -> pd.DataFrame:
    rows = []
    for lab, w in a["windows"].items():
        for s in STATS:
            for al, v in sorted(w["thr"][s].items(), reverse=True):
                rows.append({"Statistic": f"{STAT_LABEL[s]}, {lab} windows", "Alpha": al, "Threshold": v})
    return pd.DataFrame(rows)


def window_statistics(a: dict) -> pd.DataFrame:
    """Every window of every size with its statistics, adjusted p-values and significance flags."""
    parts = []
    for lab, w in a["windows"].items():
        t = w["table"].copy()
        t.insert(0, "Window", lab)
        parts.append(t)
    return pd.concat(parts, ignore_index=True)


def run(config: dict, log=print) -> dict:
    """Entry point used by the bridge. Writes the tables, the workbook and the plots into config['output_dir']."""
    try:
        from tools.insilico_bsa.backend import plots
        a = analyse(config, log=log)
        out = Path(config["output_dir"])
        out.mkdir(parents=True, exist_ok=True)
        mk, reg, pools = a["markers"], a["regions"], a["pools"]
        summ, sett, thr = summary_table(a), settings_table(a, config), thresholds_table(a)
        comp, loci, wst, near = a["comparison"], a["loci"], window_statistics(a), a["near"]
        mk.to_csv(out / "marker_statistics.csv", index=False)
        reg.to_csv(out / "significant_regions.csv", index=False)
        comp.to_csv(out / "window_comparison.csv", index=False)
        loci.to_csv(out / "loci_comparison.csv", index=False)
        wst.to_csv(out / "window_statistics.csv", index=False)
        if len(near):
            near.to_csv(out / "closest_to_significant.csv", index=False)
        a["snp_track"].to_csv(out / "snp_index_track.csv", index=False)
        a["snp_regions"].to_csv(out / "snp_index_regions.csv", index=False)
        pools.to_csv(out / "pool_members.csv", index=False)
        thr.to_csv(out / "thresholds.csv", index=False)
        xlsx = out / "InSilico_BSA.xlsx"
        with pd.ExcelWriter(xlsx, engine="openpyxl") as xw:
            summ.to_excel(xw, sheet_name="Summary", index=False)
            sett.to_excel(xw, sheet_name="Settings", index=False)
            comp.to_excel(xw, sheet_name="Window comparison", index=False)
            loci.to_excel(xw, sheet_name="Loci", index=False)
            reg.to_excel(xw, sheet_name="Regions", index=False)
            a["snp_regions"].to_excel(xw, sheet_name="SNP-index regions", index=False)
            a["snp_track"].to_excel(xw, sheet_name="SNP-index track", index=False)
            thr.to_excel(xw, sheet_name="Thresholds", index=False)
            wst.to_excel(xw, sheet_name="Windows", index=False)
            if len(near):
                near.to_excel(xw, sheet_name="Closest to significant", index=False)
            pools.to_excel(xw, sheet_name="Pools", index=False)
            mk.to_excel(xw, sheet_name="Markers", index=False)
        pl = plots.make_all(a, str(out))
        zpath = out / "InSilico_BSA_all_files.zip"
        with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
            for f in sorted(out.iterdir()):
                if f.is_file() and f != zpath:
                    z.write(f, f.name)
        s = dict(zip(summ["Metric"], summ["Value"]))
        top = mk.sort_values("Gw", ascending=False)[["Marker", "Chromosome", "Position", "Gw", "G1", "G2", "dAF", "Het_low", "Het_high", "P_adj_Gw"]]
        return {
            "status": "success",
            "message": (f"{len(a['loci'])} {'locus' if len(a['loci']) == 1 else 'loci'} for {a['trait']} ("
                        + ", ".join(f"{lab}: {len(w['regions'])}" for lab, w in a["windows"].items()) + f" regions; {a['pool_size']} lines per pool)"),
            "summary": s, "summary_rows": _records(summ), "settings": _records(sett), "regions": _records(reg), "regions_total": int(len(reg)), "snp_regions": _records(a["snp_regions"], 100), "snp_parent": a["snp_parent"],
            "comparison": _records(comp), "loci": _records(loci, 200), "near_regions": _records(near, 60),
            "windows_top": _records(wst.sort_values("Gw", ascending=False), 60), "thresholds": _records(thr), "top_markers": _records(top, 40), "pools": _records(pools, 400), "pools_total": int(len(pools)),
            "notes": a["notes"], "plots": pl,
            "files": {"workbook": str(xlsx), "zip": str(zpath), "markers": str(out / "marker_statistics.csv"), "regions": str(out / "significant_regions.csv"), "comparison": str(out / "window_comparison.csv"),
                      "loci": str(out / "loci_comparison.csv"), "near": str(out / "closest_to_significant.csv") if len(near) else "", "windows": str(out / "window_statistics.csv"),
                      "pools": str(out / "pool_members.csv"), "snp_track": str(out / "snp_index_track.csv"), "snp_regions": str(out / "snp_index_regions.csv")},
            "counts": {"lines": a["n_lines"], "markers": a["n_markers_used"], "pool": a["pool_size"], "regions": int(len(reg)), "permutations": a["n_perm"]},
        }
    except BsaError as e:
        return {"status": "error", "message": str(e)}
    except Exception as e:                                       # noqa: BLE001
        import traceback
        traceback.print_exc()
        return {"status": "error", "message": f"The analysis failed: {e}"}
