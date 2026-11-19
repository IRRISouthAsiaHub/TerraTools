"""
gff_genes.py — genes of a GFF3 annotation that lie in the regions of an input table.

The input is the QTL table of QTL Mapping (qtl_significant_*.csv, positions in cM) together with a marker position file
(a Marker, Chromosome, bp table, or a VCF), which is REQUIRED: the base-pair position of the left / peak / right marker of every QTL
is read from it, and the confidence interval is placed on the genome between them. A table of intervals that is already in base
pairs (Chromosome, Start, End) is not an input of this tool.
Columns are found by their NAME, in any order. Excel workbooks are read sheet by sheet.

Every gene line of the annotation inside a region is reported with ALL of its GFF details (name, biotype,
description/function, and any other attribute). With a reference FASTA each gene's genomic sequence is added,
reverse-complemented for genes on the minus strand, and a Gene_sequences.fasta is written (Excel cells hold at most
32,767 characters). The unique gene IDs are also written as a plain list, Gene_ID_list.txt.
"""

from __future__ import annotations

import gzip
import math
import re
from pathlib import Path

import pandas as pd

from tools.common.seqtools import Fasta, _attrs, _strip_prefix, norm_chrom, revcomp

CORE = ["Chromosome", "Start", "Stop", "Strand", "Gene_id", "Description"]
DETAIL = ["Name", "Feature", "Source"]
# attributes already shown in a dedicated column
_USED = {"ID", "Name", "description", "gene_id"}   # gene_id repeats the ID
EXCEL_CELL = 32000
MAX_ATTR_COLS = 40
MAX_REGIONS = 20000


class GeneExtractionError(Exception):
    pass


def _chr_key(c):
    """Chr01, chr1, 1 and the GWAS style S1 are the same chromosome."""
    m = re.fullmatch(r"[Ss](\d+)", str(c).strip())
    return norm_chrom(m.group(1)) if m else norm_chrom(c)


def _opener(path):
    return gzip.open if str(path).lower().endswith(".gz") else open


# ------------------------------------------------------------
# GFF3
# ------------------------------------------------------------

def read_gff_genes(path):
    """One dict per gene feature (type ends in 'gene'), keeping every attribute of its GFF line."""
    genes = []
    with _opener(path)(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 9 or not f[2].lower().endswith("gene"):
                continue
            try:
                start, stop = int(f[3]), int(f[4])
            except ValueError:
                continue
            a = _attrs(f[8])
            gid = _strip_prefix(a.get("ID", "")) or a.get("Name", "")
            if not gid:
                continue
            genes.append({
                "Chromosome": f[0], "Start": start, "Stop": stop,
                "Strand": f[6] if f[6] in ("+", "-") else ".",
                "Gene_id": gid,
                # "description" (Ensembl, RAP-DB) or "Note" (MSU, others) holds the function
                "Description": a.get("description") or a.get("Note", ""),
                "Name": a.get("Name", ""), "Feature": f[2], "Source": f[1],
                "attrs": {k: v for k, v in a.items() if k not in _USED},
            })
    if not genes:
        raise GeneExtractionError(
            "No gene features were found in the annotation file. It must be a GFF3 with 'gene' lines.")
    return genes


# ------------------------------------------------------------
# REGIONS  (any table of intervals or positions)
# ------------------------------------------------------------

_CHROM = {"chr", "chrom", "chromosome", "chromosomename", "chromo", "linkagegroup", "lg", "scaffold", "contig",
          "seqid", "seqname", "chromosomeid"}
_START = {"start", "qtlstart", "begin", "left", "lower", "from", "startbp", "startpos", "regionstart", "intervalstart",
          "leftpos", "leftposition", "cistart", "lb", "chromstart", "bpstart"}
_END = {"end", "qtlend", "stop", "right", "upper", "to", "endbp", "endpos", "regionend", "intervalend",
        "rightpos", "rightposition", "ciend", "ub", "chromend", "bpend"}
_POS = {"pos", "position", "bp", "bppos", "physicalposition", "snppos", "snpposition", "location", "peak", "peakpos",
        "peakposition", "positionbp", "coordinate", "basepair", "bpposition"}
_ID = {"qtlid", "id", "snp", "snpid", "marker", "markername", "markerid", "name", "variant", "variantid", "rs", "rsid",
       "regionid", "region", "peakid", "mta", "qtl", "locus"}
_PVAL = {"p", "pval", "pvalue", "pvalues", "pvalueadj"}
_KNOWN = _CHROM | _START | _END | _POS | _ID | _PVAL
# Chr1_1124393, chr3:456789, S5_1200, Chr1_12345_A_G
_MARKER = re.compile(r"^(?P<c>[^_:\-@\s]+?)[_:@\-](?P<p>\d+)(?:[_:@\-][A-Za-z,<>*]+(?:[_:@\-][A-Za-z,<>*]+)?)?$")


def _read_tables(path):
    """Every sheet of an Excel file (sheet name, DataFrame), or the one table of a CSV/TSV/TXT file."""
    p = str(path).lower()
    if p.endswith((".xlsx", ".xls")):
        return list(pd.read_excel(path, header=None, dtype=str, sheet_name=None).items())
    sep = "\t" if p.endswith((".tsv", ".txt", ".tab")) else ("," if p.endswith(".csv") else None)
    return [("", pd.read_csv(path, header=None, dtype=str, sep=sep, engine="python", encoding="utf-8-sig"))]


def _is_number(x):
    try:
        return not math.isnan(float(x))
    except (TypeError, ValueError):
        return False


def _norm(h):
    return re.sub(r"[^a-z0-9]", "", str(h).lower())


def _unit(h):
    """Scale to bp from a column name such as Left_Mb or Peak_kb."""
    h = str(h).lower()
    if re.search(r"(^|[^a-z])mbp?($|[^a-z])", h):
        return 1e6
    if re.search(r"(^|[^a-z])kbp?($|[^a-z])", h):
        return 1e3
    return 1.0


def chromosome_of(qtl_id):
    """'QTL_3.1' -> '3', 'QTL_Chr03.2' -> 'Chr03' (the part before the last '.')."""
    s = re.sub(r"^qtl[_\-. ]*", "", str(qtl_id).strip(), flags=re.I)
    return s.rsplit(".", 1)[0] if "." in s else s


def _numeric_share(body, i):
    v = body.iloc[:, i].dropna()
    return 0.0 if len(v) == 0 else sum(_is_number(x) for x in v) / len(v)


def _find(heads, raw, body, exact, subs=(), numeric=True):
    """Column whose normalised name is in `exact`, else contains one of `subs`; numeric ones only when asked.
    Names that carry a unit (bp, Mb, kb) win over unit-less ones (cM positions)."""
    cand = []
    for i, h in enumerate(heads):
        hit = 0 if h in exact else (1 if any(sub in h for sub in subs) else None)
        if hit is None:
            continue
        if numeric and _numeric_share(body, i) < 0.8:
            continue
        if re.search(r"(^|[^a-z])cm($|[^a-z])", str(raw[i]).lower()):
            continue                                  # a genetic-map position (CI_Start_cM): not base pairs
        has_unit = _unit(raw[i]) != 1.0 or h.endswith("bp")
        cand.append((0 if has_unit else 1, hit, i))
    return min(cand)[2] if cand else None


def _cell(row, i):
    v = row.iloc[i]
    return "" if v is None or str(v).strip().lower() in ("nan", "none") else str(v).strip()


def _marker_bp(name, chrom):
    """Physical position (bp) carried by a marker name (chr11_22044151, IRGSP1_C02_34852782, S2_1234567), or None."""
    from tools.qtl.backend.qtl_regions import physical_mb
    mb = physical_mb(name, chromosome_of(chrom) if not _is_number(chrom) else chrom)
    if mb is not None:
        return mb * 1e6
    m = _MARKER.match(str(name))
    return float(m.group("p")) if m and _chr_key(m.group("c")) == _chr_key(chrom) else None


# ------------------------------------------------------------
# MARKER POSITIONS  (a table of Marker / Chromosome / bp, or a VCF)
# ------------------------------------------------------------

_MARK_ID = {"marker", "markername", "markerid", "snp", "snpid", "id", "name", "variant", "variantid", "rs", "rsid",
            "probeset", "probesetid", "locus"}
_BP_NAMES = ("physicalbp", "bp", "bppos", "positionbp", "bpposition", "physicalposition", "basepair", "bpposition")


def _canon(name):
    """Marker name as a lookup key: lower case, without a leading 'chr' and without zeros in front of the chromosome number,
    so that chr01_18102213, Chr1_18102213 and 01_18102213 are the same marker."""
    k = str(name).strip().lower()
    k = re.sub(r"^chr(?:omosome)?[_\-]?", "", k)
    m = re.match(r"^0*(\d+)(?=[_:@])", k)
    return (m.group(1) + k[m.end():]) if m else k


def _positions_from_vcf(path, wanted):
    """{marker: (chromosome, bp)} from the CHROM / POS / ID columns of a VCF. A variant is known by its ID and also as
    CHROM_POS and CHROM:POS (the names a position-carrying marker name has). With `wanted`, only those names are kept."""
    out = {}
    with _opener(path)(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 4)
            if len(f) < 3:
                continue
            try:
                bp = float(int(f[1]))
            except ValueError:
                continue
            names = [n for n in f[2].split(";") if n and n != "."] + [f"{f[0]}_{f[1]}", f"{f[0]}:{f[1]}"]
            for n in names:
                k = _canon(n)
                if wanted is None or k in wanted:
                    out.setdefault(k, (f[0], bp))
    return out


def _positions_from_table(path, wanted):
    tables = _read_tables(path)
    tables.sort(key=lambda t: 0 if "physical" in t[0].lower() else 1)    # the 'Physical Map' sheet of Linkage Mapping first
    first_problem = None
    for _sheet, df in tables:
        if df.shape[0] < 2:
            continue
        raw = [str(c).strip() for c in df.iloc[0]]
        heads = [_norm(c) for c in raw]
        body = df.iloc[1:]
        i_m = next((i for i, h in enumerate(heads) if h in _MARK_ID), None)
        i_c = next((i for i, h in enumerate(heads) if h == "physicalchromosome"), None)
        if i_c is None:
            i_c = _find(heads, raw, body, _CHROM, numeric=False)
        i_b = next((i for n in _BP_NAMES for i, h in enumerate(heads) if h == n and _numeric_share(body, i) >= 0.8), None)
        if i_b is None:                                    # a unit-less 'Position' counts only when it cannot be cM
            i_b = _find(heads, raw, body, {"pos", "position", "location", "coordinate"}, ("position",))
            if i_b is not None:
                vals = [float(x) for x in body.iloc[:, i_b] if _is_number(x)]
                if not vals or max(vals) * _unit(raw[i_b]) < 1000:
                    i_b = None
        if None in (i_m, i_c, i_b) or len({i_m, i_c, i_b}) < 3:
            if first_problem is None:
                missing = [n for n, i in (("Marker", i_m), ("Chromosome", i_c), ("bp (Position_bp, bp, Physical_bp)", i_b)) if i is None]
                first_problem = (missing, raw)
            continue
        i_src = next((i for i, h in enumerate(heads) if h == "physicalsource"), None)
        k = _unit(raw[i_b])
        out = {}
        for _, row in body.iterrows():
            name, chrom = _cell(row, i_m), _cell(row, i_c)
            if not name or not chrom or not _is_number(row.iloc[i_b]):
                continue
            if i_src is not None and _cell(row, i_src).lower().startswith("cm x"):
                continue                                   # extrapolated from a bp-per-cM factor: not a real position
            key = _canon(name)
            if wanted is None or key in wanted:
                out.setdefault(key, (chrom, float(row.iloc[i_b]) * k))
        return out
    found = ""
    if first_problem:
        found = (f" No column was recognised for: {', '.join(first_problem[0])}. The header row of the file reads: "
                 f"{', '.join(first_problem[1][:8])}{' …' if len(first_problem[1]) > 8 else ''}. Check the spelling of these headers.")
    raise GeneExtractionError(
        "The marker position file needs a header row with a Marker column, a Chromosome column and a bp column "
        "(Position_bp, bp, Physical_bp …), or it must be a VCF." + found)


def read_marker_positions(path, wanted=None):
    """{lower-case marker name: (chromosome, bp)} from a table (CSV / TSV / Excel: Marker, Chromosome, bp — also the
    'Physical Map' sheet of a Linkage Mapping workbook) or a VCF (.vcf / .vcf.gz). `wanted`: lower-case names to keep."""
    if str(path).lower().endswith((".vcf", ".vcf.gz")):
        return _positions_from_vcf(path, wanted)
    return _positions_from_table(path, wanted)


from tools.qtl.backend.qtl_regions import physical_position  # noqa: E402


def _anchor(name, lg, positions, state):
    """(chromosome, bp) of a marker: from the marker position file, else from the position its name carries."""
    hit = positions.get(_canon(name)) if positions else None
    if hit:
        state["used"] = True
        return hit
    bp = _marker_bp(name, lg)
    if bp is not None:
        return lg, bp
    own = physical_position(name)                    # the name names another chromosome than the QTL table does
    if own:
        state["name_chr"] = True
        return str(own[0]), float(own[1])
    return None


def _qtl_bp_interval(row, cols, positions, mode, state):
    """Base-pair interval of one QTL of the QTL table of QTL Mapping (positions in cM, flanking markers by name).
    -> (chromosome, start, stop, basis), or (None, None, None, reason).

    mode 'ci'    : the confidence interval (CI_Start_cM, CI_End_cM) is placed between the markers of the QTL (left, peak,
                   right) whose cM and bp are both known, by linear interpolation.
          'flank': the span of the left and right flanking markers (wider)."""
    import numpy as np

    lg = _cell(row, cols["chr"])
    pts = {}
    for key, pkey in (("lm", "lp"), ("pm", "pp"), ("rm", "rp")):
        if cols.get(key) is None:
            continue
        a = _anchor(_cell(row, cols[key]), lg, positions, state)
        if a is None:
            continue
        cm = float(row.iloc[cols[pkey]]) if cols.get(pkey) is not None and _is_number(row.iloc[cols[pkey]]) else None
        pts[key] = (a[0], a[1], cm)
    if "lm" not in pts or "rm" not in pts:
        gone = [_cell(row, cols[k]) for k in ("lm", "rm") if k not in pts]
        return None, None, None, "no base-pair position for the marker(s) " + ", ".join(gone)
    if _chr_key(pts["lm"][0]) != _chr_key(pts["rm"][0]):
        return None, None, None, "the left and right markers are on different chromosomes"
    chrom = pts["lm"][0]
    pts = {k: v for k, v in pts.items() if _chr_key(v[0]) == _chr_key(chrom)}
    bps = [pts[k][1] for k in ("lm", "pm", "rm") if k in pts]
    odd = len(bps) == 3 and not (bps[0] <= bps[1] <= bps[2] or bps[0] >= bps[1] >= bps[2])
    flag = "; the map order of these markers differs from the genome, so the interval may be wide" if odd else ""

    ci = [row.iloc[cols[k]] for k in ("cs", "ce") if cols.get(k) is not None]
    anchors = sorted({(v[2], v[1]) for v in pts.values() if v[2] is not None})
    if mode == "ci" and len(ci) == 2 and all(_is_number(x) for x in ci) and len({a[0] for a in anchors}) >= 2:
        lo, hi = sorted(float(x) for x in ci)
        cms, bpv = np.array([a[0] for a in anchors]), np.array([a[1] for a in anchors])
        vals = [float(np.interp(lo, cms, bpv)), float(np.interp(hi, cms, bpv))] + [float(v) for v in bpv[(cms >= lo) & (cms <= hi)]]
        return chrom, min(vals), max(vals), "confidence interval placed between the QTL's markers" + flag
    return chrom, min(bps[0], bps[-1]), max(bps[0], bps[-1]), "flanking markers" + flag


def _parse_table(df, notes, positions=None, mode="ci", state=None):
    """-> (regions, "interval"); None when no region columns are recognised."""
    state = state if state is not None else {"used": False}
    if df.shape[0] < 2:
        return None
    raw = [str(c).strip() for c in df.iloc[0]]
    heads = [_norm(c) for c in raw]
    numeric_row = df.shape[1] >= 3 and _is_number(df.iloc[0, 1]) and _is_number(df.iloc[0, 2])
    has_header = not numeric_row and any(
        h in _KNOWN or any(k in h for k in ("start", "pos", "chrom", "marker", "snp")) for h in heads)

    if not has_header:
        if df.shape[1] >= 3 and _numeric_share(df, 1) >= 0.8 and _numeric_share(df, 2) >= 0.8:
            state["interval_table"] = True                     # ID, start, end without a header: an interval file, not a QTL table
        return None

    body = df.iloc[1:]
    i_chr = _find(heads, raw, body, _CHROM, numeric=False)
    i_id = _find(heads, raw, body, _ID, numeric=False)
    i_s = _find(heads, raw, body, _START, ("start", "begin", "left"))
    i_e = _find(heads, raw, body, _END, ("end", "stop", "right"))
    if i_s is not None and i_s == i_e:
        i_e = None

    def chrom_of_row(row):
        if i_chr is not None and _cell(row, i_chr):
            return _cell(row, i_chr)
        return chromosome_of(_cell(row, i_id)) if i_id is not None else ""

    out = []
    # a QTL table (QTL Mapping): positions are in cM, but the flanking marker names carry the physical position
    i_lm = next((i for i, h in enumerate(heads) if "marker" in h and "left" in h), None)
    i_rm = next((i for i, h in enumerate(heads) if "marker" in h and "right" in h), None)
    unitless = i_s is None or i_e is None or all(_unit(raw[i]) == 1.0 and not heads[i].endswith("bp") for i in (i_s, i_e))
    if i_lm is not None and i_rm is not None and i_chr is not None and unitless:
        cols = {"chr": i_chr, "lm": i_lm, "rm": i_rm}
        for key, head in (("pm", "peakmarker"), ("lp", "leftposition"), ("pp", "peakposition"), ("rp", "rightposition"),
                          ("cs", "cistartcm"), ("ce", "ciendcm")):
            cols[key] = next((i for i, h in enumerate(heads) if h == head), None)
        skipped = []
        for n, (_, row) in enumerate(body.iterrows(), 1):
            chrom, start, stop, why = _qtl_bp_interval(row, cols, positions, mode, state)
            rid = _cell(row, i_id) if i_id is not None else ""
            if start is None:
                skipped.append(f"{rid or 'row ' + str(n)}: {why}")
                continue
            out.append({"Region_ID": rid or f"{chrom}:{int(start)}-{int(stop)}", "Chromosome": chrom, "Start": start,
                        "Stop": stop, "Left_Marker": _cell(row, i_lm), "Right_Marker": _cell(row, i_rm), "Interval_Basis": why})
        if out and state.get("name_chr"):
            notes.append("For some QTL the chromosome of the QTL table differs from the chromosome written in the names of their markers "
                         "(for example a QTL listed on chromosome 4 with markers named chr01_...). The chromosome and position of the "
                         "marker names were used: check that the QTL table and the marker file come from the same data.")
        if out:
            notes.append("The positions of this table are in cM, so the base-pair interval was worked out from the "
                         + ("base-pair positions of its markers in the marker position file (and in the marker names where "
                            "the file lacks them)." if positions is not None else "physical positions in the names of its markers."))
            if skipped:
                notes.append(f"{len(skipped)} QTL(s) could not be converted to base pairs and were left out: "
                             + "; ".join(skipped[:5]) + (" …" if len(skipped) > 5 else ""))
            return out, "interval"
        if positions is not None:
            raise GeneExtractionError(
                "None of the QTL could be converted to base pairs with the marker position file ("
                + "; ".join(skipped[:3]) + f"). {len(positions)} of the markers named in the QTL table were found in the file. "
                "The marker names of the file (the ID column of a VCF) must be the same as the Left_Marker / Right_Marker "
                "names of the QTL table.")
        raise GeneExtractionError("None of the QTL of this table could be converted to base pairs.")
    if i_s is not None and i_e is not None and (i_chr is not None or i_id is not None):
        state["interval_table"] = True                         # Chromosome + Start + End: an interval file, not a QTL table
    return None


def read_regions(path, positions_path, interval="ci"):
    """Regions of the QTL table of QTL Mapping (positions in cM). The base-pair position of its left / peak / right markers comes from
    `positions_path`, which is REQUIRED (a Marker / Chromosome / bp table, or a VCF). interval: 'ci' or 'flank' (see _qtl_bp_interval).
    Returns (regions, notes)."""
    if not positions_path:
        raise GeneExtractionError("Please select the marker position file (Marker, Chromosome, bp) or, if there is none, a VCF file: "
                                  "it gives the base-pair position of the markers of the QTL table.")
    notes, found, sheet = [], None, ""
    tables = _read_tables(path)
    # a workbook of a QTL run: look first at the sheet that holds the intervals
    tables.sort(key=lambda t: 0 if re.search(r"qtl|interval|region|peak|sig", t[0].lower()) else 1)
    state = {"used": False}
    wanted = {_canon(v) for _n, d in tables for v in d.to_numpy().ravel() if isinstance(v, str)}
    positions = read_marker_positions(positions_path, wanted)
    for name, df in tables:
        found = _parse_table(df, notes, positions, interval, state)
        if found:
            sheet = name
            break
    if found is None and all(df.shape[0] < 2 for _, df in tables):
        raise GeneExtractionError(
            "The input file has no rows, only a header: that QTL Mapping run found no QTL, so there is nothing to extract. "
            "Use the QTL table (qtl_significant_*.csv) of a run that found QTL.")
    if found is None and state.get("interval_table"):
        raise GeneExtractionError(
            "This file is a table of intervals (Chromosome, Start, End), which Genes Extraction does not take. Use the QTL table of "
            "QTL Mapping (qtl_significant_*.csv, positions in cM) together with the marker position file or a VCF.")
    if found is None:
        raise GeneExtractionError(
            "No QTL could be read from the input file. Genes Extraction takes the QTL table of QTL Mapping (qtl_significant_*.csv) "
            "with its Left_Marker / Right_Marker columns. A table of single SNP positions, such as a GWAS result, is not a QTL result.")
    regions, _kind = found
    if max(r["Stop"] for r in regions) < 1000:
        raise GeneExtractionError(
            "The positions in the input file are all below 1,000, so they cannot be base pairs (they look like cM "
            "or Mb). Give positions in bp, name the column with its unit (for example Left_Mb, Right_Mb), or use the "
            "QTL table of QTL Mapping together with a marker position file.")
    if sheet and len(tables) > 1:
        notes.append(f"The regions were read from the sheet '{sheet}'.")
    if len(regions) > MAX_REGIONS:
        raise GeneExtractionError(f"The file holds {len(regions):,} regions. Genes Extraction is meant for the QTL of a run.")
    return regions, notes


# ------------------------------------------------------------
# EXTRACTION
# ------------------------------------------------------------

def _wrap(seq, n=60):
    return "\n".join(seq[i:i + n] for i in range(0, len(seq), n))


def _number(config, key, default=None):
    v = config.get(key)
    if v is None or str(v).strip() == "":
        return default
    try:
        return float(v)
    except ValueError:
        raise GeneExtractionError(f"'{key}' must be a number, not '{v}'.")


def run_gff_extraction(config, results_dir):
    """config: marker_path (the QTL table of QTL Mapping, positions in cM), positions_path (REQUIRED marker position file or VCF),
    gff_path, optional fasta_path, interval ('ci' or 'flank'),
    match ('inside' = gene completely inside the region, default; 'overlap' = any overlap)."""
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    interval = "flank" if config.get("interval") == "flank" else "ci"
    regions, notes = read_regions(config["marker_path"], config.get("positions_path"), interval)
    overlap = config.get("match") == "overlap"

    genes = read_gff_genes(config["gff_path"])
    by_chr = {}
    for g in genes:
        by_chr.setdefault(_chr_key(g["Chromosome"]), []).append(g)

    if not any(_chr_key(r["Chromosome"]) in by_chr for r in regions):
        raise GeneExtractionError(
            "None of the chromosome names in the annotation ("
            + ", ".join(list(dict.fromkeys(g["Chromosome"] for g in genes))[:4])
            + ") matches the chromosomes of the input file ("
            + ", ".join(list(dict.fromkeys(str(r["Chromosome"]) for r in regions))[:4]) + ").")
    n_off = sum(_chr_key(r["Chromosome"]) not in by_chr for r in regions)
    warnings = list(notes)
    if n_off:
        warnings.append(f"{n_off} region(s) are on chromosomes that are not in the annotation and were skipped.")

    fasta = None
    if config.get("fasta_path"):
        fasta = Fasta(config["fasta_path"])
        if not any(fasta.has(r["Chromosome"]) for r in regions):
            fasta.close()
            raise GeneExtractionError("None of the chromosome names in the reference FASTA matches the "
                                      "chromosomes of the input file.")

    rows, summary, attr_keys = [], [], {}
    try:
        for reg in regions:
            cand = by_chr.get(_chr_key(reg["Chromosome"]), [])
            if overlap:
                hit = [g for g in cand if g["Start"] <= reg["Stop"] and g["Stop"] >= reg["Start"]]
            else:
                hit = [g for g in cand if g["Start"] >= reg["Start"] and g["Stop"] <= reg["Stop"]]
            hit.sort(key=lambda g: (g["Start"], g["Stop"]))
            info = {"Region_ID": reg["Region_ID"], "Chromosome": reg["Chromosome"]}
            info.update({"Region_Start": int(reg["Start"]), "Region_End": int(reg["Stop"])})
            info.update({k: reg[k] for k in ("Left_Marker", "Right_Marker", "Interval_Basis") if k in reg})
            info["Genes"] = len(hit)
            summary.append(info)
            for g in hit:
                row = {k: g[k] for k in CORE}
                row.update({"Region_ID": reg["Region_ID"], "Region_Start": int(reg["Start"]),
                            "Region_End": int(reg["Stop"])})
                row.update({k: g[k] for k in DETAIL})
                row["Length_bp"] = g["Stop"] - g["Start"] + 1
                for k, v in g["attrs"].items():
                    if k in attr_keys or len(attr_keys) < MAX_ATTR_COLS:
                        attr_keys[k] = True
                        row[k] = v
                if fasta is not None:
                    seq = fasta.fetch(g["Chromosome"], g["Start"], g["Stop"])
                    row["Sequence"] = (revcomp(seq) if g["Strand"] == "-" else seq) if seq else ""
                rows.append(row)
    finally:
        if fasta is not None:
            fasta.close()

    base = CORE + ["Region_ID", "Region_Start", "Region_End"] + DETAIL + ["Length_bp"]
    cols = base + list(attr_keys) + (["Sequence"] if fasta is not None else [])
    df = pd.DataFrame(rows, columns=cols).fillna("")

    # ---- gene ID list: every gene once, with the regions it was found in
    id_rows = []
    if len(df):
        for gid, grp in df.groupby("Gene_id", sort=False):
            r0 = grp.iloc[0]
            id_rows.append({"Gene_id": gid, "Chromosome": r0["Chromosome"], "Start": r0["Start"], "Stop": r0["Stop"],
                            "Strand": r0["Strand"], "Regions": "; ".join(dict.fromkeys(grp["Region_ID"].astype(str)))})
    id_df = pd.DataFrame(id_rows, columns=["Gene_id", "Chromosome", "Start", "Stop", "Strand", "Regions"])
    list_file = results_dir / "Gene_ID_list.txt"
    list_file.write_text("".join(f"{g}\n" for g in id_df["Gene_id"]), encoding="utf-8", newline="\n")

    fasta_file = None
    if fasta is not None:
        no_seq = int((df["Sequence"] == "").sum())
        if no_seq:
            warnings.append(f"{no_seq} gene(s) have no sequence: their chromosome is missing from the reference FASTA.")
        fasta_file = results_dir / "Gene_sequences.fasta"
        seen = set()
        with open(fasta_file, "w", encoding="utf-8", newline="\n") as fh:
            for _, r in df.iterrows():
                key = (r["Gene_id"], r["Start"], r["Stop"])
                if not r["Sequence"] or key in seen:
                    continue
                seen.add(key)
                fh.write(f">{r['Gene_id']} {r['Chromosome']}:{r['Start']}-{r['Stop']}({r['Strand']}) {r['Region_ID']}"
                         f" {r['Description']}".rstrip() + "\n" + _wrap(r["Sequence"]) + "\n")
        long_ = df["Sequence"].str.len() > EXCEL_CELL
        if long_.any():
            warnings.append(f"{int(long_.sum())} sequence(s) are longer than an Excel cell can hold; "
                            "they are left out of the workbook and are in Gene_sequences.fasta.")
            df.loc[long_, "Sequence"] = "(longer than an Excel cell - see Gene_sequences.fasta)"

    genes_file = results_dir / "Genes_in_gene_regions.xlsx"
    with pd.ExcelWriter(genes_file) as xw:
        df.to_excel(xw, sheet_name="Genes", index=False)
        pd.DataFrame(summary).to_excel(xw, sheet_name="Regions", index=False)
        id_df.to_excel(xw, sheet_name="Gene_ID_list", index=False)

    result = {
        "status": "success",
        "input_kind": "interval",
        "genes": df.drop(columns=["Sequence"], errors="ignore").to_dict("records"),
        "columns": [c for c in cols if c != "Sequence"],
        "genes_file": str(genes_file),
        "gene_list_file": str(list_file),
        "gene_ids": int(len(id_df)),
        "has_sequences": fasta is not None,
        "warnings": warnings,
        "region_summary": summary,
        "region_columns": list(summary[0].keys()) if summary else [],
    }
    if fasta_file is not None:
        result["fasta_file"] = str(fasta_file)
    return result
