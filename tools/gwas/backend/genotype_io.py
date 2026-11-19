"""
Genotype input for GWAS: explicit encodings, strict VCF genotype parsing and marker-map validation.

Internal representation
    G        float matrix (samples x markers): diploid alternate-allele dosage 0 .. 2, NaN = missing.
    ids      sample IDs as written in the file
    snp_ids  marker IDs as written in the file
    chrom    integer chromosome per marker, pos  integer position (bp) per marker

Encodings
    "dosage"   0 / 1 / 2 (copies of the alternate allele; fractional values in [0, 2] are accepted as imputed dosages)
    "centred"  -1 / 0 / 1 (the same genotypes shifted by -1; fractional values in [-1, 1] are accepted)
    "auto"     decide only when the values leave no doubt:
                 a value above 1 (and none below 0)  -> dosage
                 a value below 0 (and none above 1)  -> centred
                 only 0 and 1 (or values in [0, 1])  -> ambiguous: the user must choose. 0 / 1 alone could be a
                 0/1/2 matrix of inbred lines scored 0 and 1, or a centred matrix without heterozygotes.
Missing values are NaN. The sentinels -9, -99 and -999 and any code the user lists are turned into NaN first; no other
value is ever treated as missing.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from io import StringIO

import numpy as np
import pandas as pd


class InputError(ValueError):
    """A problem the user can fix (shown as is in the app)."""


ENCODINGS = {"dosage": "0/1/2 alternate-allele dosage", "centred": "-1/0/1 centred genotypes"}
SENTINELS = (-9.0, -99.0, -999.0)
NA_TOKENS = ["NA", "N/A", "NaN", "nan", "NULL", "null", "", ".", "./.", ".|."]


# ============================================================================================ encodings

def _examples(values, limit=5):
    u = np.unique(values)
    return ", ".join(f"{v:g}" for v in u[:limit]) + (" …" if len(u) > limit else "")


def clean_missing(G, extra_codes=()):
    """Replace the missing-value sentinels (and the codes the user lists) by NaN. Returns (G, {code: count})."""
    G = np.array(G, dtype=float, copy=True)
    found = {}
    for code in tuple(SENTINELS) + tuple(float(c) for c in extra_codes):
        hit = G == code
        n = int(hit.sum())
        if n:
            found[code] = n
            G[hit] = np.nan
    return G, found


def detect_encoding(G):
    """(encoding or None, reason). None means the values do not decide between the encodings."""
    v = G[~np.isnan(G)]
    if v.size == 0:
        return None, "the matrix has no genotype values"
    lo, hi = float(v.min()), float(v.max())
    if lo < 0 and hi > 1:
        return None, (f"it holds both negative values and values above 1 (range {lo:g} to {hi:g}); that is neither "
                      "0/1/2 nor -1/0/1. If one of those numbers marks a missing call, list it as a missing code")
    if hi > 1:
        return "dosage", f"it holds values above 1 (range {lo:g} to {hi:g}), which only the 0/1/2 coding has"
    if lo < 0:
        return "centred", f"it holds values below 0 (range {lo:g} to {hi:g}), which only the -1/0/1 coding has"
    return None, (f"it holds only the values {_examples(v)}; that is a valid 0/1/2 matrix without heterozygotes or "
                  "without the allele count 2 and also a valid -1/0/1 matrix without the value -1")


def to_dosage(G, encoding="auto", extra_missing=()):
    """Validate `G` against the encoding and convert it to dosage 0..2 (NaN = missing).

    Returns (dosage matrix, report). Raises InputError when the encoding is ambiguous (auto) or the values do not fit it."""
    G, sentinels = clean_missing(G, extra_missing)
    requested = str(encoding or "auto").lower()
    if requested == "auto":
        # TerraTools files (Linkage Mapping, Mendelian, Genofile Converter) mark a missing call with -1 next to 0/1/2.
        # A -1/0/1 centred matrix cannot hold the value 2, so -1 together with values above 1 can only be a missing code.
        v0 = G[~np.isnan(G)]
        neg = np.unique(v0[v0 < 0])
        if neg.size == 1 and neg[0] == -1.0 and (v0 > 1).any():
            hit = G == -1.0
            sentinels = dict(sentinels)
            sentinels[-1.0] = int(hit.sum())
            G = G.copy()
            G[hit] = np.nan
    if requested not in ("auto", "dosage", "centred"):
        raise InputError(f"Unknown genotype encoding '{encoding}'. Use 0/1/2 dosage or -1/0/1 centred.")
    n_missing = int(np.isnan(G).sum())
    how = "chosen by the user"
    enc = requested
    if requested == "auto":
        enc, why = detect_encoding(G)
        if enc is None:
            raise InputError(
                "The genotype encoding cannot be told from the values: " + why + ".\n"
                "Choose the encoding explicitly: '0/1/2 dosage' (copies of the alternate allele) or '-1/0/1 centred'. "
                "TerraTools never guesses it from a matrix of 0 and 1.")
        how = "detected because " + why
    v = G[~np.isnan(G)]
    lo_ok, hi_ok = (0.0, 2.0) if enc == "dosage" else (-1.0, 1.0)
    bad = v[(v < lo_ok - 1e-9) | (v > hi_ok + 1e-9)]
    if bad.size:
        raise InputError(
            f"{bad.size:,} genotype values are not valid for the {ENCODINGS[enc]} encoding (allowed: {lo_ok:g} to {hi_ok:g}; "
            f"found e.g. {_examples(bad)}). Check the encoding, or list the code of missing calls under the missing codes.")
    out = G + 1.0 if enc == "centred" else G
    frac = int((np.abs(v - np.round(v)) > 1e-9).sum())
    report = {"requested": requested, "encoding": enc, "how": how, "converted": enc == "centred",
              "missing_calls": n_missing, "sentinels_to_missing": {f"{k:g}": c for k, c in sentinels.items()},
              "fractional_values": frac,
              "message": (f"Genotype encoding: {ENCODINGS[enc]} ({how})."
                          + (" Converted to 0/1/2 dosage by adding 1." if enc == "centred" else " Used as 0/1/2 dosage.")
                          + (f" {frac:,} non-integer values (imputed dosages) are kept as they are." if frac else "")
                          + (" Missing-value codes turned into NaN: " + ", ".join(f"{k} ({c:,})" for k, c in report_items(sentinels)) + "."
                             if sentinels else ""))}
    return out, report


def report_items(d):
    return [(f"{k:g}", c) for k, c in d.items()]


# ============================================================================================ marker map

# Chr3_123456, chr03-123456, 3:123456, and names that go on after the position: chr1_134521_LOC_Os01g01280, Chr3_123456_A_G
_NAME = re.compile(r"^(?:chr|chromosome)?[_\-]?0*(\d+)[_\-:.](\d+)(?:[_\-:.].*)?$", re.I)


@dataclass
class MarkerMap:
    chrom: np.ndarray
    pos: np.ndarray
    chrom_source: str          # "chromosome row", "marker names" or "none"
    pos_source: str            # "position row", "marker names" or "none" (positions are then the order of the markers)
    valid: bool                # real chromosomes and physical positions are known
    notes: list = field(default_factory=list)


def map_from_names(snp_ids, chrom_row=None):
    """Chromosome and position of every marker.

    Sources, in this order: an optional chromosome row of the file (then positions are unknown), marker names such as
    Chr3_123456, chr3-123456, 3:123456. Nothing is invented: if neither source works the markers are placed on one
    chromosome in file order and the map is marked invalid, so that models that need physical positions refuse to run."""
    m = len(snp_ids)
    notes = []
    if chrom_row is not None:
        chrom = np.asarray(chrom_row, dtype=int)
        # The marker names usually carry the physical position (chr1_134521, Chr3_123456_A_G): use it when EVERY name gives one
        # and the chromosome in the name is the chromosome of the row.
        parsed = [_NAME.match(str(n).strip()) for n in snp_ids]
        if all(parsed):
            nchr = np.array([int(t.group(1)) for t in parsed])
            npos = np.array([int(t.group(2)) for t in parsed], dtype=np.int64)
            if np.array_equal(nchr, chrom):
                notes.append("Chromosomes come from the chromosome row of the file and physical positions from the marker names.")
                firsts = [int(npos[chrom == c].min()) for c in dict.fromkeys(chrom.tolist())]
                lasts = [int(npos[chrom == c].max()) for c in dict.fromkeys(chrom.tolist())]
                if len(firsts) > 1 and all(firsts[i] > lasts[i - 1] for i in range(1, len(firsts))):
                    notes.append("The positions in the marker names run on from one chromosome to the next (whole-genome coordinates). "
                                 "Distances within a chromosome are still correct, which is what the models use.")
                return MarkerMap(chrom, npos, "chromosome row", "marker names", True, notes)
            notes.append(f"The chromosome in the marker names differs from the chromosome row for {int((nchr != chrom).sum()):,} markers, "
                         "so the positions in the names were not used.")
        counters, pos = {}, np.zeros(m, dtype=np.int64)
        for i, c in enumerate(chrom):
            counters[c] = counters.get(c, 0) + 1
            pos[i] = counters[c]
        notes.append("The file gives the chromosome of every marker but no positions: positions are the order of the markers "
                     "on each chromosome, not physical distances. Add a Position row, or name the markers Chr<number>_<position>.")
        return MarkerMap(chrom, pos, "chromosome row", "none", False, notes)
    chrom, pos = np.zeros(m, dtype=int), np.zeros(m, dtype=np.int64)
    ok = np.zeros(m, dtype=bool)
    for i, name in enumerate(snp_ids):
        mt = _NAME.match(str(name).strip())
        if mt:
            chrom[i], pos[i], ok[i] = int(mt.group(1)), int(mt.group(2)), True
    if ok.all():
        return MarkerMap(chrom, pos, "marker names", "marker names", True, notes)
    if ok.any():
        raise InputError(f"{int((~ok).sum()):,} of {m:,} marker names do not give a chromosome and position "
                         f"(e.g. '{snp_ids[int(np.flatnonzero(~ok)[0])]}'); the others do (e.g. '{snp_ids[int(np.flatnonzero(ok)[0])]}'). "
                         "Name every marker like Chr3_123456, or give none of them a position.")
    notes.append("The marker names carry no chromosome and position. All markers are placed on one chromosome in file order, so "
                 "the Manhattan plot shows marker order and the models that need physical positions (FarmCPU, mrMLM, FASTmrMLM) "
                 "are not available. Add a Position row under the header (and a Chromosome row), or name the markers Chr<number>_<position>, to enable them.")
    return MarkerMap(np.ones(m, dtype=int), np.arange(1, m + 1, dtype=np.int64), "none", "none", False, notes)


_POS_LABELS = ("pos", "position", "positions", "bp", "physicalposition", "physpos", "basepair", "coordinate")
_CHR_LABELS = ("chr", "chrom", "chromosome", "chromosomes", "chrs")
_NAME_CHR = re.compile(r"(?:^|[_\-:.])(?:chr|chromosome)?0*(\d+)[_\-:.]\d+", re.I)


def map_from_rows(snp_ids, chrom_row=None, pos_row=None):
    """Marker map from the rows of the file itself.

    chromosome row + position row  -> real physical map (the layout Linkage Mapping writes for QTL mapping: marker names,
                                      then a Position row and a Chromosome row, then one row per individual)
    chromosome row only            -> positions unknown (order of the markers)
    position row only              -> chromosomes taken from marker names such as Chr3_123456, otherwise an error"""
    m = len(snp_ids)
    if pos_row is None:
        return map_from_names(snp_ids, chrom_row)
    pos = np.asarray(pos_row, dtype=np.int64)
    if chrom_row is not None:
        chrom = np.asarray(chrom_row, dtype=int)
        return MarkerMap(chrom, pos, "chromosome row", "position row", True,
                         ["Marker positions (bp) were read from the position row of the file and chromosomes from its chromosome row."])
    chrom = np.zeros(m, dtype=int)
    for i, name in enumerate(snp_ids):
        mt = _NAME_CHR.search(str(name).strip())
        if not mt:
            raise InputError("The file has a position row but no chromosome row, and the marker name '%s' does not start with a "
                             "chromosome (Chr3_123456). Add a chromosome row below the position row." % name)
        chrom[i] = int(mt.group(1))
    return MarkerMap(chrom, pos, "marker names", "position row", True,
                     ["Marker positions (bp) were read from the position row of the file; chromosomes from the marker names."])


def check_map(chrom, pos):
    """Problems of a marker map: non-positive or non-integer positions."""
    problems = []
    if np.any(np.asarray(pos) <= 0):
        problems.append(f"{int((np.asarray(pos) <= 0).sum()):,} markers have a position of 0 or less")
    if np.any(np.asarray(chrom) <= 0):
        problems.append(f"{int((np.asarray(chrom) <= 0).sum()):,} markers have a chromosome of 0 or less")
    return problems


# ============================================================================================ GWAS-ready file

@dataclass
class GenotypeData:
    ids: np.ndarray
    y: np.ndarray
    snp_ids: np.ndarray
    chrom: np.ndarray
    pos: np.ndarray
    G: np.ndarray                       # dosage, NaN = missing (never altered afterwards: QC works on copies)
    encoding: dict
    marker_map: MarkerMap
    ref: np.ndarray | None = None
    alt: np.ndarray | None = None
    notes: list = field(default_factory=list)


def _read_table(path):
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        head = fh.readline()
    sep = max((",", "\t", ";"), key=head.count)
    if head.count(sep) < 2:
        return pd.read_csv(path, sep=None, engine="python", na_values=NA_TOKENS, keep_default_na=True), head, None
    return pd.read_csv(path, sep=sep, na_values=NA_TOKENS, keep_default_na=True), head, sep


def _norm_label(v):
    return re.sub(r"[^0-9a-z]", "", str(v).lower()) if not pd.isna(v) else ""


def read_gwas_file(path=None, text=None, encoding="auto", extra_missing=()):
    """Read a GWAS-ready file: column 1 = sample ID, column 2 = phenotype, then one column per marker.
    Optional rows under the header without a phenotype hold the marker positions (bp) and/or the chromosome of every marker."""
    if path is not None:
        raw, head, sep = _read_table(path)
    elif text and text.strip():
        head = text.splitlines()[0]
        sep = max((",", "\t", ";"), key=head.count)
        raw = pd.read_csv(StringIO(text), sep=sep, na_values=NA_TOKENS)
    else:
        raise InputError("No data provided.")
    if raw.shape[1] < 4:
        raise InputError("The genotype file needs an ID column, a phenotype column and at least two marker columns.")
    raw_cols = [c.strip().strip('"') for c in head.rstrip("\r\n").split(sep or ",")]
    dup_snp = sorted({c for c in raw_cols[2:] if raw_cols[2:].count(c) > 1})
    if dup_snp:
        raise InputError(f"{len(dup_snp)} marker IDs appear more than once in the header, e.g. {', '.join(dup_snp[:5])}. "
                         "Marker IDs must be unique.")

    cols = list(raw.columns)
    id_col, phen_col, snp_cols = cols[0], cols[1], cols[2:]
    notes = []

    chrom_row, pos_row = None, None
    map_rows = []                      # leading rows without a phenotype that hold one number per marker
    for k in range(min(2, len(raw) - 1)):
        row = raw.iloc[k]
        if not pd.isna(pd.to_numeric(row[phen_col], errors="coerce")):
            break
        vals = pd.to_numeric(row[snp_cols], errors="coerce")
        if vals.isna().any() or not (np.abs(vals - vals.round()) < 1e-9).all():
            break
        label = _norm_label(row[id_col])
        known = label in _POS_LABELS or label in _CHR_LABELS
        blank = label in ("", "nan")
        if not (known or (blank and (float(vals.min()) >= 1 or map_rows))):
            break
        map_rows.append((k, label, vals.astype(float).astype(np.int64).to_numpy()))
    for k, label, vals in map_rows:
        if label in _POS_LABELS:
            pos_row = vals
        elif label in _CHR_LABELS:
            chrom_row = vals
    unlabelled = [r for r in map_rows if r[1] not in _POS_LABELS and r[1] not in _CHR_LABELS]
    if len(unlabelled) == 2:           # no labels: the row with the larger numbers holds base-pair positions
        a, b = sorted(unlabelled, key=lambda r: int(r[2].max()))
        chrom_row, pos_row = a[2], b[2]
    elif len(unlabelled) == 1:
        only = unlabelled[0][2]
        if chrom_row is None and pos_row is None:
            # a lone unlabelled row: chromosome numbers are small, base-pair positions are not
            if int(only.max()) > 1000:
                pos_row = only
            else:
                chrom_row = only
        elif chrom_row is None:
            chrom_row = only
        else:
            pos_row = only
    if map_rows:
        raw = raw.drop(index=[r[0] for r in map_rows]).reset_index(drop=True)
        what = (["positions"] if pos_row is not None else []) + (["chromosomes"] if chrom_row is not None else [])
        notes.append("The first %d row(s) under the header were read as the marker %s." % (len(map_rows), " and ".join(what)))

    pheno = pd.to_numeric(raw[phen_col], errors="coerce")
    ok = pheno.notna()
    n_bad = int((~ok).sum())
    if n_bad:
        notes.append(f"{n_bad:,} samples without a numeric phenotype were left out.")
    raw, pheno = raw[ok].reset_index(drop=True), pheno[ok].reset_index(drop=True)
    ids = raw[id_col].astype(str).str.strip().to_numpy()
    dup_ids = sorted({i for i in ids if (ids == i).sum() > 1})
    if dup_ids:
        raise InputError(f"{len(dup_ids)} sample IDs appear more than once, e.g. {', '.join(dup_ids[:5])}. Sample IDs must be "
                         "unique (average repeated measurements before the analysis).")

    geno = raw[snp_cols]
    if not all(pd.api.types.is_numeric_dtype(t) for t in geno.dtypes):
        geno = geno.apply(pd.to_numeric, errors="coerce")
    G, enc = to_dosage(geno.to_numpy(dtype=float), encoding, extra_missing)
    snp_ids = np.array([str(c) for c in snp_cols])
    mmap = map_from_rows(snp_ids, chrom_row, pos_row)
    problems = check_map(mmap.chrom, mmap.pos)
    if problems:
        raise InputError("Invalid marker map: " + "; ".join(problems) + ".")
    return GenotypeData(ids=ids, y=pheno.to_numpy(dtype=float), snp_ids=snp_ids, chrom=mmap.chrom, pos=mmap.pos, G=G,
                        encoding=enc, marker_map=mmap, notes=notes + mmap.notes)


# ============================================================================================ VCF genotypes

_GT_FAST = {}
for _a in "01.":
    for _b in "01.":
        for _s in "/|":
            _t = f"{_a}{_s}{_b}"
            _GT_FAST[_t] = np.nan if "." in (_a, _b) else float((_a == "1") + (_b == "1"))


class VcfGenotypeError(InputError):
    pass


def parse_gt(token, where=""):
    """Dosage (0, 1 or 2 copies of ALT) of one diploid biallelic GT; NaN if any allele is missing.

    Phased (|) and unphased (/) calls are accepted. Other ploidies and allele numbers above 1 raise an error: they are
    never read as something else."""
    gt = token.split(":", 1)[0]
    v = _GT_FAST.get(gt)
    if v is not None:
        return v
    alleles = re.split(r"[/|]", gt)
    if len(alleles) != 2:
        raise VcfGenotypeError(f"Unsupported ploidy in the genotype '{gt}'{where}: only diploid calls (0/1, 1|1, ./.) are "
                               "supported. For haploid, polyploid or pooled-sequencing data use another tool.")
    if any(a not in ("0", "1", ".") for a in alleles):
        raise VcfGenotypeError(f"The genotype '{gt}'{where} names an allele other than 0 or 1; only biallelic SNPs are supported.")
    return np.nan


POOL_NAME = re.compile(r"bulk|pool", re.I)


def reject_pooled(sample_names):
    pooled = [s for s in sample_names if POOL_NAME.search(str(s))]
    if pooled:
        raise InputError(
            "Samples named like pooled sequencing bulks were found (" + ", ".join(pooled[:4]) + "). The genotype of a pooled bulk "
            "is not the diploid genotype of one individual, so it cannot go into a GWAS. Use QTL-Seq for bulked samples, or "
            "rename the samples if they are individual plants.")
