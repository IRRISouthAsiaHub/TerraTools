"""
Recombination Analysis — crossovers, donor content and recombination hot / cold regions of a mapping population.

Input: a genotype table with one row per marker (chromosome, position in bp, marker name) and one column per sample,
coded 0 / 1 / 2 / -1 (-1 = missing). Two kinds of file are read:

  * a table that holds the recurrent parent as one of its samples (a VCF-style table with  CHROM, POS, ID, …, <sample>_GT);
    the parent's homozygous allele is "A", the other homozygous allele is "B" (the donor);
  * a table that is already coded against its two parents (what the Genofile Converter writes: Marker, Chromosome, then
    the individuals; 0 = Parent 1 = "A", 2 = Parent 2 = "B").

Method (the analysis the lab ran as R scripts, now one reproducible tool):
  1. markers that cannot show recombination are dropped: the parent must be homozygous and at least one line must carry the
     other homozygous allele (with an optional donor sample: the donor must be homozygous for the other allele);
  2. every line gets A / B / H (het) / U (missing) calls per marker;
  3. isolated calls are smoothed away: a run of consecutive A or B calls shorter than a minimum block (a % of the markers of
     the chromosome) is a genotyping error rather than a double crossover;
  4. a crossover is a change A <-> B between two neighbouring informative markers; its interval (BP_start, BP_end) is
     the region in which the breakpoint lies;
  5. donor content, crossover rate by bin, and (when the marker names carry the gene, chr1_134521_LOC_Os01g01280) the
     split into intragenic and intergenic crossovers.
"""

from __future__ import annotations

import io
import os
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
GENE_NAMES = HERE / "data" / "gene_names.csv"

# Centromere (start, end) in Mb of the 12 rice chromosomes (IRGSP-1.0) used for the shading of the rate plot.
RICE_CENTROMERES = {1: (16.0, 20.0), 2: (3.5, 7.0), 3: (14.0, 17.5), 4: (5.5, 9.0), 5: (9.5, 12.5), 6: (5.0, 8.0),
                    7: (6.5, 9.5), 8: (4.0, 7.5), 9: (7.0, 10.0), 10: (8.0, 11.5), 11: (13.0, 16.5), 12: (9.5, 13.0)}

HAP_U, HAP_A, HAP_B, HAP_H = 0, 1, 2, 3
HAP_LETTER = {HAP_A: "A", HAP_B: "B", HAP_H: "H", HAP_U: "U"}

MAX_ROWS = 400                      # rows of a table that are sent to the page (the full table is in the files)


class RecombError(ValueError):
    """A problem the user can fix; shown as it is in the app."""


# ───────────────────────────── reading ─────────────────────────────

_CHR_NAMES = {"chrom", "chr", "chromosome", "#chrom", "chromosomes", "linkagegroup", "lg"}
_POS_NAMES = {"pos", "position", "bp", "physicalposition", "basepair", "coordinate", "start"}
_ID_NAMES = {"id", "marker", "markerid", "markername", "snp", "snpid", "rs", "name", "locus"}
_SKIP_NAMES = {"ref", "alt", "qual", "filter", "info", "format", "strand", "alleles", "cm", "gene", "geneid", "end",
               "qccode", "assembly", "center", "protlsid", "assayLSID", "panelLSID".lower()}
_NAME_POS_PLAIN = re.compile(r"^(?:chr|chromosome)?[_\-]?0*(\d+)[_\-:.](\d+)(?:[_\-:.](.+))?$", re.I)
# SNP names that start with the assembly: MSU7_1_4797026_C-T, IRGSP1_C01_22891448, Os_Nipponbare_3_14933451
_NAME_POS_ASSEMBLY = re.compile(r"^[A-Za-z][A-Za-z]*\d*(?:_[A-Za-z]+)*?_C?0*(\d+)_(\d{3,})(?:_(.+))?$")


class _NamePosPattern:
    """Marker name -> match with group 1 = chromosome label, group 2 = position (chr1_134521, chr01_1044946, MSU7_1_4797026_C-T, IRGSP1_C01_22891448)."""

    @staticmethod
    def match(text):
        return _NAME_POS_PLAIN.match(text) or _NAME_POS_ASSEMBLY.match(text)


_NAME_POS = _NamePosPattern()
_CHR_LABEL = re.compile(r"^(?:chr|chromosome)?[_\-\s]?0*(\d+)$", re.I)


def _norm(s) -> str:
    return re.sub(r"[^0-9a-z]", "", str(s).lower())


def _read_table(path: str, sheet=None) -> pd.DataFrame:
    low = path.lower()
    try:
        if low.endswith((".xlsx", ".xlsm", ".xls")):
            return pd.read_excel(path, sheet_name=sheet or 0, dtype=str)
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            head = fh.readline()
        sep = max((",", "\t", ";"), key=head.count)
        return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, low_memory=False)
    except Exception as e:
        raise RecombError(f"Could not read the file: {e}") from None


def chromosome_number(label):
    """1 for 'Chr01', 'chr1', '1', 'chromosome 1'; None for scaffolds, organelles and anything else."""
    m = _CHR_LABEL.match(str(label).strip())
    return int(m.group(1)) if m else None


def _geno_columns(df: pd.DataFrame, skip: set) -> list:
    """Columns whose values are genotype codes 0 / 1 / 2 / -1 (or empty / NA)."""
    cols = []
    for c in df.columns:
        if c in skip:
            continue
        s = df[c].astype(str).str.strip()
        s = s[~s.str.upper().isin(["", "NA", "N/A", "NAN", ".", "./."])]
        if s.empty:
            continue
        v = pd.to_numeric(s, errors="coerce")
        if v.notna().mean() >= 0.98 and v.dropna().isin([-1, 0, 1, 2]).mean() >= 0.98:
            cols.append(c)
    return cols


def read_genotypes(path, sheet=None) -> dict:
    """-> dict(markers=DataFrame[ID, Chr, Chrom, POS], samples=[...], G=int8 matrix markers x samples, layout, notes)."""
    path = str(path)
    if not os.path.isfile(path):
        raise RecombError("The genotype file was not found.")
    df = _read_table(path, sheet)
    if df.shape[1] < 4 or df.shape[0] < 3:
        raise RecombError("The genotype file needs marker rows and at least three more columns (chromosome, markers, samples).")
    df.columns = [str(c).strip() for c in df.columns]
    # QTL layout (markers across): its positions are cM, not base pairs
    first = [_norm(v) for v in df.iloc[:2, 0]]
    if len(first) == 2 and first[0] in ("chromosome", "chr") and first[1] in ("position", "pos", "cm"):
        raise RecombError("This is the QTL-mapping layout (markers across, positions in cM). Recombination Analysis needs physical "
                          "positions in base pairs: use a table with one row per marker and a CHROM / POS (or Marker / Chromosome) "
                          "column, for example the output of the Genofile Converter.")
    cols = {_norm(c): c for c in df.columns}
    chr_col = next((cols[k] for k in cols if k in _CHR_NAMES), None)
    pos_col = next((cols[k] for k in cols if k in _POS_NAMES), None)
    id_col = next((cols[k] for k in cols if k in _ID_NAMES), None)
    skip = {c for c in (chr_col, pos_col, id_col) if c} | {cols[k] for k in cols if k in _SKIP_NAMES}
    gcols = _geno_columns(df, skip)
    if len(gcols) < 3:
        raise RecombError("Fewer than 3 columns hold genotype codes 0 / 1 / 2 / -1. Code the genotypes first (Genofile Converter).")
    notes = []
    ids = df[id_col].astype(str).str.strip() if id_col else pd.Series([f"m{i + 1}" for i in range(len(df))])
    parsed = ids.map(lambda s: _NAME_POS.match(s))
    if chr_col:
        chrom_label = df[chr_col].astype(str).str.strip()
    elif parsed.notna().mean() >= 0.5:
        chrom_label = parsed.map(lambda m: m.group(1) if m else "")
    else:
        raise RecombError("The file has no chromosome column and the marker names do not carry the chromosome (chr1_134521).")
    if pos_col:
        pos = pd.to_numeric(df[pos_col].astype(str).str.replace(",", ""), errors="coerce")
        pos_source = f"column '{pos_col}'"
    elif parsed.notna().mean() >= 0.5:                      # a few special markers (a QTL name, a gene marker) carry no position: they are left out
        pos = parsed.map(lambda m: float(m.group(2)) if m else np.nan)
        pos_source = "the marker names"
    else:
        raise RecombError("The file has no position column and the marker names do not end in a position (chr1_134521). "
                          "Recombination Analysis needs the physical position of every marker in base pairs.")
    chrom = chrom_label.map(chromosome_number)
    unplaced = int(chrom.isna().sum())
    bad_pos = int(pos.isna().sum())
    keep = chrom.notna() & pos.notna() & (pos > 0)
    if unplaced:
        notes.append(f"{unplaced:,} markers on scaffolds or non-numeric chromosomes were left out.")
    if bad_pos and bad_pos != unplaced:
        notes.append(f"{int(bad_pos):,} markers without a valid position were left out.")
    if keep.sum() < 20:
        raise RecombError("Fewer than 20 markers have a numeric chromosome and a position.")
    sub = df.loc[keep]
    G = (sub[gcols].apply(pd.to_numeric, errors="coerce").fillna(-1).to_numpy())
    G = np.where(np.isin(G, (0, 1, 2)), G, -1).astype(np.int8)
    markers = pd.DataFrame({"ID": ids[keep].to_numpy(), "Chrom": chrom[keep].astype(int).to_numpy(), "POS": pos[keep].astype(np.int64).to_numpy()})
    order = np.lexsort((markers["POS"].to_numpy(), markers["Chrom"].to_numpy()))
    markers, G = markers.iloc[order].reset_index(drop=True), G[order]
    samples = [re.sub(r"_GT$", "", c) for c in gcols]
    dup = sorted({s for s in samples if samples.count(s) > 1})
    if dup:
        raise RecombError("Sample names must be unique; repeated: " + ", ".join(dup[:5]))
    layout = "marker rows with a parent column" if chr_col and pos_col else "marker rows"
    return {"markers": markers, "samples": samples, "G": G, "layout": layout, "pos_source": pos_source, "notes": notes,
            "id_col": id_col, "n_rows": len(df)}


def guess_parent(samples: list, G: np.ndarray | None = None, has_pos_column: bool = True) -> str:
    """The recurrent parent: a sample named like one, else the first sample (the layout of the lab's 4k files)."""
    for s in samples:
        if re.search(r"parent|recurrent|(^|[^a-z])p1($|[^a-z0-9])|^p1$", s, re.I):
            return s
    return samples[0] if (has_pos_column and samples) else ""


def inspect_file(path, sheet=None) -> dict:
    d = read_genotypes(path, sheet)
    m, G, s = d["markers"], d["G"], d["samples"]
    has_pos_col = d["layout"] == "marker rows with a parent column"
    gene_ids = m["ID"].str.extract(r"^chr[0-9]+[_\-:.][0-9]+[_\-:.](.+)$")[0]
    n_gene = int(gene_ids.notna().sum())
    het = (G == 1).mean(axis=0)
    miss = (G == -1).mean(axis=0)
    prev_cols = ["Marker", "Chromosome", "Position"] + s[:6]
    prev = pd.DataFrame({"Marker": m["ID"][:8], "Chromosome": m["Chrom"][:8], "Position": m["POS"][:8]})
    for j, name in enumerate(s[:6]):
        prev[name] = G[:8, j]
    chrs = sorted(m["Chrom"].unique().tolist())
    return {
        "status": "success", "file": os.path.basename(str(path)), "layout": d["layout"], "samples": s, "n_samples": len(s),
        "n_markers": int(len(m)), "chromosomes": chrs, "pos_source": d["pos_source"],
        "guess_parent": guess_parent(s, G, has_pos_col), "coded": not has_pos_col or False,
        "gene_in_names": n_gene >= 0.5 * len(m), "n_gene_names": n_gene,
        "het_pct": round(float(het.mean() * 100), 2), "missing_pct": round(float(miss.mean() * 100), 2),
        "notes": d["notes"],
        "sample_stats": [{"sample": n, "het_pct": round(float(h * 100), 1), "missing_pct": round(float(x * 100), 1)}
                         for n, h, x in zip(s[:400], het[:400], miss[:400])],
        "preview_columns": prev_cols, "preview_rows": prev.to_dict("records"),
    }


# ───────────────────────────── smoothing ─────────────────────────────

def _runs(vals):
    """[(value, length), …] of a 1-D sequence."""
    out = []
    for v in vals:
        if out and out[-1][0] == v:
            out[-1][1] += 1
        else:
            out.append([v, 1])
    return out


def _expand(runs):
    return np.concatenate([np.full(n, v, dtype=np.int8) for v, n in runs]) if runs else np.empty(0, dtype=np.int8)


def smooth_legacy(ab: np.ndarray, min_block: int) -> np.ndarray:
    """The smoothing of the original R script: the first run shorter than the minimum takes the value of the run before it
    (the first run takes that of the run after it); repeated at most 50 times."""
    if len(ab) < 2:
        return ab
    runs = _runs(ab.tolist())
    for _ in range(50):
        if all(n >= min_block for _, n in runs):
            break
        i = next(i for i, (_, n) in enumerate(runs) if n < min_block)
        if i == 0:
            if len(runs) == 1:
                break
            rv = runs[1][0]
        else:
            rv = runs[i - 1][0]
        runs[i][0] = rv
        merged = []
        for r in runs:
            if merged and merged[-1][0] == r[0]:
                merged[-1][1] += r[1]
            else:
                merged.append(list(r))
        runs = merged
    return _expand(runs)


def smooth_robust(ab: np.ndarray, min_block: int) -> np.ndarray:
    """Shortest run first, into the longer of its two neighbours (the only neighbour at a chromosome end), until every
    run has the minimum length. No iteration cap, and a real double crossover is only removed when it is the shortest."""
    if len(ab) < 2:
        return ab
    runs = _runs(ab.tolist())
    while len(runs) > 1:
        lens = [n for _, n in runs]
        n_min = min(lens)
        if n_min >= min_block:
            break
        i = lens.index(n_min)
        if i == 0:
            j = 1
        elif i == len(runs) - 1:
            j = i - 1
        else:
            j = i - 1 if runs[i - 1][1] >= runs[i + 1][1] else i + 1
        runs[i][0] = runs[j][0]
        merged = []
        for r in runs:
            if merged and merged[-1][0] == r[0]:
                merged[-1][1] += r[1]
            else:
                merged.append(list(r))
        runs = merged
    return _expand(runs)


# ───────────────────────────── the analysis ─────────────────────────────

def _as_bool(v, default=False):
    if v is None:
        return default
    if isinstance(v, str):
        return v.strip().lower() in ("1", "true", "yes", "on")
    return bool(v)


def load_gene_names() -> dict:
    try:
        t = pd.read_csv(GENE_NAMES, dtype=str)
        return dict(zip(t["LOC_ID"], t["Gene_Name"]))
    except Exception:
        return {}


def analyse(config: dict, log=print) -> dict:
    """Everything except the files: returns a dict of DataFrames and arrays used by the writers and the page."""
    d = read_genotypes(config["input_file"], config.get("sheet") or None)
    markers, G, samples = d["markers"], d["G"], d["samples"]
    notes = list(d["notes"])
    parent = (config.get("parent") or "").strip()
    donor = (config.get("donor") or "").strip()
    exclude = {str(s) for s in (config.get("exclude_samples") or [])}
    for nm, label in ((parent, "parent"), (donor, "donor")):
        if nm and nm not in samples:
            raise RecombError(f"The {label} '{nm}' is not a sample of the file.")
    if parent and parent == donor:
        raise RecombError("The recurrent parent and the donor must be different samples.")
    smoothing = str(config.get("smoothing") or "legacy")
    if smoothing not in ("legacy", "robust", "none"):
        raise RecombError(f"Unknown smoothing method '{smoothing}'.")
    pct = float(config.get("min_block_pct", 5.0))
    floor = int(config.get("min_block_min", 3))
    min_donor_rils = max(1, int(config.get("min_donor_rils", 1)))
    bin_mb = float(config.get("bin_mb", 0.5))
    if not (0.05 <= bin_mb <= 20):
        raise RecombError("The bin size must be between 0.05 and 20 Mb.")
    if not (0 <= pct <= 50) or floor < 1:
        raise RecombError("The minimum block must be 0-50 % of the markers of a chromosome, and at least 1 marker.")

    drop = {parent, donor} | exclude
    ril_idx = [i for i, s in enumerate(samples) if s not in drop]
    ril_names = [samples[i] for i in ril_idx]
    if len(ril_names) < 5:
        raise RecombError(f"Only {len(ril_names)} lines are left to analyse; at least 5 are needed.")
    R = G[:, ril_idx]                                     # markers x lines

    # 1. polymorphic markers
    total = len(markers)
    if parent:
        pg = G[:, samples.index(parent)]
        a_allele = pg.astype(np.int8)
        ok = (pg == 0) | (pg == 2)
        reason_parent = int((~ok).sum())
        if donor:
            dg = G[:, samples.index(donor)]
            ok &= (dg == 0) | (dg == 2)
            ok &= dg != pg
        b_allele = np.where(a_allele == 0, 2, 0).astype(np.int8)
        n_donor_like = (R == b_allele[:, None]).sum(axis=1)
    else:
        a_allele = np.zeros(total, dtype=np.int8)
        b_allele = np.full(total, 2, dtype=np.int8)
        ok = np.ones(total, dtype=bool)
        reason_parent = 0
        n_donor_like = (R == 2).sum(axis=1)
    poly = ok & (n_donor_like >= min_donor_rils)
    n_poly = int(poly.sum())
    if n_poly < 20:
        raise RecombError(f"Only {n_poly} polymorphic markers remain (parent homozygous and at least {min_donor_rils} line(s) with the "
                          "other allele). Check the parent, or that the genotypes are coded 0 / 1 / 2.")
    pm = markers[poly].reset_index(drop=True)
    Rp, ap, bp = R[poly], a_allele[poly], b_allele[poly]
    log(f"[Recombination] {n_poly:,} polymorphic markers of {total:,}; {len(ril_names)} lines")

    # 2. haplotype calls
    hap = np.full(Rp.shape, HAP_U, dtype=np.int8)
    hap[Rp == ap[:, None]] = HAP_A
    hap[Rp == bp[:, None]] = HAP_B
    hap[Rp == 1] = HAP_H
    raw_counts = {k: int((hap == v).sum()) for k, v in (("A", HAP_A), ("B", HAP_B), ("H", HAP_H), ("U", HAP_U))}

    # 3. smoothing per line and chromosome
    chroms = sorted(pm["Chrom"].unique().tolist())
    cidx = {c: np.flatnonzero(pm["Chrom"].to_numpy() == c) for c in chroms}
    min_block = {c: max(floor, int(round(len(ix) * pct / 100.0))) for c, ix in cidx.items()}
    sm = hap.copy()
    changed = 0
    fn = {"legacy": smooth_legacy, "robust": smooth_robust}.get(smoothing)
    if fn is not None:
        for c, ix in cidx.items():
            for j in range(hap.shape[1]):
                h = hap[ix, j]
                ab = np.flatnonzero((h == HAP_A) | (h == HAP_B))
                if len(ab) < 2:
                    continue
                new = fn(h[ab], min_block[c])
                changed += int((new != h[ab]).sum())
                sm[ix[ab], j] = new

    # 4. crossovers
    pos = pm["POS"].to_numpy()
    pid = pm["ID"].to_numpy()
    co_rows = []
    for j, ril in enumerate(ril_names):
        for c, ix in cidx.items():
            h = sm[ix, j]
            keep = np.flatnonzero((h == HAP_A) | (h == HAP_B))
            if len(keep) < 2:
                continue
            hk = h[keep]
            for k in np.flatnonzero(hk[1:] != hk[:-1]):
                l, r = ix[keep[k]], ix[keep[k + 1]]
                co_rows.append((ril, c, int(pos[l]), int(pos[r]), (pos[l] + pos[r]) / 2.0, int(pos[r] - pos[l]), int(r - l - 1),
                                HAP_LETTER[int(hk[k])], HAP_LETTER[int(hk[k + 1])], pid[l], pid[r]))
    co = pd.DataFrame(co_rows, columns=["RIL", "Chrom", "BP_start", "BP_end", "BP_mid", "BP_size", "N_SNPs", "From", "To", "Marker_left", "Marker_right"])

    # 5. intragenic / intergenic (needs the gene in the marker names)
    loc = pm["ID"].str.extract(r"^chr[0-9]+[_\-:.][0-9]+[_\-:.](.+)$")[0]
    has_genes = bool(loc.notna().mean() >= 0.5)
    names = load_gene_names()
    if has_genes and len(co):
        by_id = dict(zip(pm["ID"], loc))
        co["LOC_left"] = co["Marker_left"].map(by_id)
        co["LOC_right"] = co["Marker_right"].map(by_id)
        same = co["LOC_left"].notna() & (co["LOC_left"] == co["LOC_right"])
        co["CO_type"] = np.where(same, "intragenic", "intergenic")
        co["LOC_ID"] = np.where(same, co["LOC_left"], "")
        co["Gene_Name"] = [names.get(g, "") if s else "" for g, s in zip(co["LOC_left"], same)]
    elif len(co):
        co["CO_type"] = "unknown"
        notes.append("The marker names do not carry a gene (chr1_134521_LOC_Os01g01280), so crossovers were not split into intragenic "
                     "and intergenic.")

    # 6. per line / per chromosome
    n_rils = len(ril_names)
    n_ab_chr = {c: np.zeros(n_rils) for c in chroms}
    n_b_chr = {c: np.zeros(n_rils) for c in chroms}
    for c, ix in cidx.items():
        h = sm[ix]
        n_ab_chr[c] = ((h == HAP_A) | (h == HAP_B)).sum(axis=0).astype(float)
        n_b_chr[c] = (h == HAP_B).sum(axis=0).astype(float)
    donor_chr = pd.DataFrame({f"Chr{c}": np.where(n_ab_chr[c] > 0, np.round(100 * n_b_chr[c] / np.maximum(n_ab_chr[c], 1), 2), 0.0) for c in chroms},
                             index=ril_names)
    tot_ab = sum(n_ab_chr.values())
    tot_b = sum(n_b_chr.values())
    co_ril = co.groupby("RIL").size().reindex(ril_names).fillna(0).astype(int) if len(co) else pd.Series(0, index=ril_names)
    per_ril = pd.DataFrame({
        "RIL": ril_names,
        "Crossovers": co_ril.to_numpy(),
        "Donor_pct": donor_chr.mean(axis=1).round(2).to_numpy(),
        "Donor_pct_all_markers": np.round(100 * tot_b / np.maximum(tot_ab, 1), 2),
        "Het_pct": np.round(100 * (hap == HAP_H).mean(axis=0), 2),
        "Missing_pct": np.round(100 * (hap == HAP_U).mean(axis=0), 2),
    })
    for c in chroms:
        per_ril[f"Chr{c}_donor_pct"] = donor_chr[f"Chr{c}"].to_numpy()
        per_ril[f"Chr{c}_COs"] = (co[co["Chrom"] == c].groupby("RIL").size().reindex(ril_names).fillna(0).astype(int).to_numpy()
                                  if len(co) else 0)

    span_mb = {c: float((markers.loc[markers["Chrom"] == c, "POS"].max() - markers.loc[markers["Chrom"] == c, "POS"].min()) / 1e6) for c in chroms}
    chr_len_mb = {c: float(markers.loc[markers["Chrom"] == c, "POS"].max() / 1e6) for c in chroms}
    rows = []
    for c in chroms:
        sub = co[co["Chrom"] == c] if len(co) else co
        row = {"Chromosome": c, "Polymorphic_markers": int(len(cidx[c])), "Length_Mb": round(chr_len_mb[c], 2),
               "Min_block_markers": int(min_block[c]), "Crossovers": int(len(sub)), "Mean_COs_per_RIL": round(len(sub) / n_rils, 3),
               "COs_per_RIL_per_100Mb": round(100 * len(sub) / n_rils / max(chr_len_mb[c], 1e-9), 2),
               "Mean_donor_pct": round(float(donor_chr[f"Chr{c}"].mean()), 2)}
        if has_genes and len(sub):
            row["Intragenic"] = int((sub["CO_type"] == "intragenic").sum())
            row["Intergenic"] = int((sub["CO_type"] == "intergenic").sum())
        rows.append(row)
    per_chr = pd.DataFrame(rows)

    # 7. crossover rate by bin
    bin_rows = []
    for c in chroms:
        L = chr_len_mb[c]
        edges = np.arange(0, L + bin_mb, bin_mb)
        mids = co.loc[co["Chrom"] == c, "BP_mid"].to_numpy() / 1e6 if len(co) else np.empty(0)
        cnt, _ = np.histogram(mids, bins=edges)
        for k, n in enumerate(cnt):
            bin_rows.append((c, round(edges[k], 4), round(edges[k + 1], 4), int(n), n / n_rils))
    bins = pd.DataFrame(bin_rows, columns=["Chromosome", "Bin_start_Mb", "Bin_end_Mb", "Crossovers", "COs_per_RIL"])

    # 8. genes
    genes = pd.DataFrame()
    if has_genes and len(co):
        ig = co[co["CO_type"] == "intragenic"]
        if len(ig):
            genes = (ig.groupby(["LOC_ID", "Gene_Name", "Chrom"]).agg(Intragenic_COs=("RIL", "size"), RILs=("RIL", "nunique")).reset_index()
                     .sort_values(["Intragenic_COs", "LOC_ID"], ascending=[False, True]).reset_index(drop=True))
            genes.insert(0, "Rank", np.arange(1, len(genes) + 1))

    summary = [
        ("Lines analysed", n_rils), ("Markers in file", total), ("Polymorphic markers used", n_poly),
        ("Markers removed (monomorphic / parent not homozygous)", total - n_poly),
        ("Total crossovers", int(len(co))), ("Mean crossovers per line", round(len(co) / n_rils, 2)),
        ("Min crossovers per line", int(co_ril.min())), ("Max crossovers per line", int(co_ril.max())),
        ("Lines with no crossover", int((co_ril == 0).sum())),
        ("Mean donor % (mean of the chromosomes)", round(float(per_ril["Donor_pct"].mean()), 2)),
        ("Mean donor % (all markers pooled)", round(float(per_ril["Donor_pct_all_markers"].mean()), 2)),
        ("Min donor %", round(float(per_ril["Donor_pct"].min()), 2)), ("Max donor %", round(float(per_ril["Donor_pct"].max()), 2)),
        ("Heterozygous calls %", round(100 * raw_counts["H"] / hap.size, 2)), ("Missing calls %", round(100 * raw_counts["U"] / hap.size, 2)),
        ("Calls changed by smoothing", changed),
    ]
    if has_genes and len(co):
        ni = int((co["CO_type"] == "intragenic").sum())
        summary += [("Intragenic crossovers", ni), ("Intergenic crossovers", int(len(co) - ni)),
                    ("Intragenic %", round(100 * ni / len(co), 1))]
    if smoothing == "none":
        notes.append("No smoothing was applied: every isolated call counts as two crossovers, so genotyping errors inflate the counts.")
    sparse = per_ril[per_ril["Missing_pct"] > 50]
    if len(sparse):
        notes.append(f"{len(sparse)} line(s) have more than 50 % missing calls ({', '.join(sparse['RIL'].astype(str)[:5])}"
                     f"{' …' if len(sparse) > 5 else ''}); their crossovers cannot be located. Leave them out under 'Samples to leave out'.")
    notes.append("A crossover is placed in the interval between its two flanking markers and plotted at the middle of that interval, so "
                 "the rate along the chromosome is only as fine as the marker spacing.")
    if (per_ril["Het_pct"] > 25).any():
        notes.append(f"{int((per_ril['Het_pct'] > 25).sum())} line(s) are more than 25 % heterozygous: they are early-generation lines, "
                     "and heterozygous calls are not used to place a crossover.")
    return {"markers": markers, "poly_markers": pm, "ril_names": ril_names, "hap": hap, "smooth": sm, "co": co, "per_ril": per_ril,
            "per_chr": per_chr, "bins": bins, "genes": genes, "summary": pd.DataFrame(summary, columns=["Metric", "Value"]),
            "chroms": chroms, "chr_len_mb": chr_len_mb, "has_genes": has_genes, "notes": notes, "min_block": min_block,
            "n_total": total, "n_poly": n_poly, "reason_parent": reason_parent, "parent": parent, "donor": donor,
            "smoothing": smoothing, "pct": pct, "floor": floor, "bin_mb": bin_mb, "raw_counts": raw_counts, "donor_chr": donor_chr,
            "layout": d["layout"], "pos_source": d["pos_source"], "rice_centromeres": _as_bool(config.get("centromeres"), True)}


# ───────────────────────────── outputs ─────────────────────────────

def _records(df: pd.DataFrame, n=MAX_ROWS):
    out = df.head(n).copy()
    for c in out.columns:
        if out[c].dtype.kind == "f":
            out[c] = out[c].round(4)
    return out.where(pd.notna(out), None).to_dict("records")


def _settings_table(a: dict, config: dict) -> pd.DataFrame:
    rows = [("Input file", os.path.basename(str(config["input_file"]))), ("Layout", a["layout"]), ("Positions read from", a["pos_source"]),
            ("Recurrent parent", a["parent"] or "(none: file is coded, 0 = Parent 1 = A, 2 = Parent 2 = B)"),
            ("Donor sample", a["donor"] or "(none: B = the other homozygous allele)"),
            ("Smoothing", {"legacy": "legacy (as in the original R script)", "robust": "shortest run into its longer neighbour",
                           "none": "none"}[a["smoothing"]]),
            ("Minimum block", f"{a['pct']:g} % of the markers of a chromosome, at least {a['floor']}"),
            ("Bin size (Mb)", a["bin_mb"])]
    return pd.DataFrame(rows, columns=["Setting", "Value"])


def run(config: dict) -> dict:
    """Entry point used by the bridge. Writes the tables, the workbook and the plots into config['output_dir']."""
    try:
        from tools.recombination.backend import plots
        a = analyse(config)
        out = Path(config["output_dir"])
        out.mkdir(parents=True, exist_ok=True)
        a["co"].to_csv(out / "crossover_breakpoints.csv", index=False)
        a["per_ril"].to_csv(out / "per_RIL_summary.csv", index=False)
        a["per_chr"].to_csv(out / "per_chromosome_summary.csv", index=False)
        a["summary"].to_csv(out / "recombination_summary.csv", index=False)
        a["bins"].to_csv(out / "crossover_rate_by_bin.csv", index=False)
        a["poly_markers"].rename(columns={"ID": "Marker", "Chrom": "Chromosome"}).to_csv(out / "polymorphic_markers.csv", index=False)
        if len(a["genes"]):
            a["genes"].to_csv(out / "genes_with_intragenic_crossovers.csv", index=False)
        xlsx = out / "Recombination_Analysis.xlsx"
        with pd.ExcelWriter(xlsx, engine="openpyxl") as xw:
            a["summary"].to_excel(xw, sheet_name="Summary", index=False)
            _settings_table(a, config).to_excel(xw, sheet_name="Settings", index=False)
            a["per_chr"].to_excel(xw, sheet_name="Per_chromosome", index=False)
            a["per_ril"].to_excel(xw, sheet_name="Per_RIL", index=False)
            a["co"].to_excel(xw, sheet_name="Crossovers", index=False)
            a["bins"].to_excel(xw, sheet_name="CO_rate_by_bin", index=False)
            if len(a["genes"]):
                a["genes"].to_excel(xw, sheet_name="Genes", index=False)
        pl = plots.make_all(a, str(out))
        zpath = out / "Recombination_Analysis_all_files.zip"
        with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
            for f in sorted(out.iterdir()):
                if f.is_file() and f != zpath:
                    z.write(f, f.name)
        co = a["co"]
        s = dict(zip(a["summary"]["Metric"], a["summary"]["Value"]))
        return {
            "status": "success", "message": f"{len(co):,} crossovers in {len(a['ril_names'])} lines",
            "summary": s, "summary_rows": _records(a["summary"]), "per_chr": _records(a["per_chr"]),
            "per_ril": _records(a["per_ril"][["RIL", "Crossovers", "Donor_pct", "Donor_pct_all_markers", "Het_pct", "Missing_pct"]].sort_values("Crossovers", ascending=False)),
            "per_ril_total": int(len(a["per_ril"])),
            "crossovers": _records(co.drop(columns=[c for c in ("Marker_left", "Marker_right") if c in co.columns])),
            "crossovers_total": int(len(co)), "genes": _records(a["genes"], 30), "has_genes": a["has_genes"],
            "settings": _records(_settings_table(a, config)), "notes": a["notes"], "plots": pl,
            "files": {"crossovers": str(out / "crossover_breakpoints.csv"), "per_ril": str(out / "per_RIL_summary.csv"),
                      "per_chr": str(out / "per_chromosome_summary.csv"), "bins": str(out / "crossover_rate_by_bin.csv"),
                      "workbook": str(xlsx), "zip": str(zpath), "polymorphic": str(out / "polymorphic_markers.csv")},
            "counts": {"markers": a["n_total"], "polymorphic": a["n_poly"], "lines": len(a["ril_names"]), "chromosomes": len(a["chroms"])},
        }
    except RecombError as e:
        return {"status": "error", "message": str(e)}
    except Exception as e:                                   # noqa: BLE001
        import traceback
        traceback.print_exc()
        return {"status": "error", "message": f"The analysis failed: {e}"}
