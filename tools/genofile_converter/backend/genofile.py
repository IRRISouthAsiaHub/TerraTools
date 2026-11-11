"""
genofile.py - Genofile Converter: genotyping-company files -> TerraTools input files.

Reads the three formats genotyping companies deliver
    * HapMap   rs# alleles chrom pos strand ... then one column per sample (AG, AA, NN, or an IUPAC letter)
    * VCF      #CHROM POS ID REF ALT ... FORMAT, one column per sample (GT field: 0/0, 0/1, 1|1, ./.)   (.vcf or .vcf.gz)
    * Raw      one row per marker, Excel / CSV: Marker, Chromosome, P1, P2, offspring ... (A:G, A/G, A|G or AG)
and codes every offspring call against the two parents:
     0  homozygous for the P1 allele
     1  heterozygous P1/P2
     2  homozygous for the P2 allele
    -1  missing, or an allele combination the parents cannot explain

A marker is left out when
    - either parent call is missing/invalid        -> "Invalid parent"
    - either parent is heterozygous                -> "Non-homozygous parent"
    - both parents carry the same allele           -> "Non-segregating"
    - (optional) it lies on an unplaced contig     -> "Unplaced contig"

Outputs (one coded table, Marker | Chromosome | one 0/1/2/-1 column per individual):
    <name>_coded.csv    the coded table
    <name>_coded.xlsx   the same table as the FIRST sheet, then report sheets (Marker_Details, Excluded_Markers, Summary)
HapMap and VCF files are streamed line by line, so a large file never has to fit in memory.
"""

import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

MISSING_CALLS = {"", "-", "-:-", "N", "NA", "N:N", "NAN", "NONE", "?", "./.", "-/-", "N/N", "NN", "00", "--", "XX"}
REASONS = ("Invalid parent", "Non-homozygous parent", "Non-segregating", "Unplaced contig")
SEPARATORS = {":": ":", "/": "/", "|": "|", "none": ""}   # raw files: "" means two-character calls such as "AG"
FORMATS = ("hapmap", "vcf", "raw")
FORMAT_LABEL = {"hapmap": "HapMap", "vcf": "VCF", "raw": "Raw allele calls"}

# IUPAC single-letter heterozygous codes, as written in HapMap files by some pipelines
XLSX_MAX_CELLS = 1_500_000      # coded cells (markers x columns) above which no Excel input file is written

IUPAC = {"R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC"}


# ─────────────────────────── low-level reading ───────────────────────────

def _open_text(path):
    path = Path(path)
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return open(path, "rt", encoding="utf-8", errors="replace", newline="")


def _lower_name(path):
    n = Path(path).name.lower()
    return n[:-3] if n.endswith(".gz") else n


def detect_format(path):
    """'hapmap' | 'vcf' | 'raw' from the extension and the first lines."""
    name = _lower_name(path)
    if name.endswith((".xlsx", ".xlsm", ".xls")):
        return "raw"
    if name.endswith((".vcf", ".vcf.gz")) or name.endswith(".bcf"):
        return "vcf"
    try:
        with _open_text(path) as fh:
            for line in fh:
                if line.startswith("##fileformat=VCF"):
                    return "vcf"
                if line.startswith("##"):
                    continue
                head = re.split(r"[\t,]", line.rstrip("\r\n"))
                low = [h.strip().lower() for h in head]
                if low[:1] == ["#chrom"]:
                    return "vcf"
                if low[:1] in (["rs#"], ["snp"], ["rs"]) and "alleles" in low and any(c in low for c in ("chrom", "chr")):
                    return "hapmap"
                if "alleles" in low and "pos" in low and "chrom" in low:
                    return "hapmap"
                break
    except Exception:
        pass
    return "hapmap" if name.endswith((".hmp", ".hmp.txt", ".hapmap")) else "raw"


def _key(name):
    return re.sub(r"[\s_\-\.]+", "", str(name)).lower()


def _guess(columns, *candidates):
    keys = {_key(c): c for c in columns}
    for cand in candidates:
        if _key(cand) in keys:
            return keys[_key(cand)]
    return ""


def _guess_parents(samples):
    """Sample names that look like parent 1 / parent 2."""
    p1 = _guess(samples, "P1", "Parent1", "Parent 1", "Parent_A", "PA", "Female", "Parent-1")
    p2 = _guess(samples, "P2", "Parent2", "Parent 2", "Parent_B", "PB", "Male", "Parent-2")
    return p1, p2


# ─────────────────────────── calls -> allele pairs ───────────────────────────

def parse_raw_call(gt, sep=":"):
    """'A:G' -> ('A', 'G'); None for missing/unparseable calls."""
    if gt is None:
        return None
    gt = str(gt).strip().upper()
    if gt in MISSING_CALLS:
        return None
    if sep:
        parts = gt.split(sep)
        if len(parts) != 2:
            return None
        a1, a2 = (a.strip() for a in parts)
    else:
        if len(gt) != 2:
            return None
        a1, a2 = gt[0], gt[1]
    if not a1 or not a2 or a1 in ("-", "N") or a2 in ("-", "N"):
        return None
    return a1, a2


def parse_hapmap_call(gt):
    """HapMap call: 'AG' / 'AA' / 'A/G' / 'A' (hom.) / 'R' (IUPAC het) / 'NN'."""
    g = str(gt).strip().upper()
    if g in MISSING_CALLS or g == "0":
        return None
    if "/" in g or ":" in g or "|" in g:
        parts = re.split(r"[/:|]", g)
        if len(parts) != 2 or not parts[0] or not parts[1] or "N" in parts or parts == ["-", "-"]:
            return None
        return parts[0], parts[1]
    if len(g) == 1:
        if g in IUPAC:
            return IUPAC[g][0], IUPAC[g][1]
        return (g, g) if g in "ACGT+" else None
    if len(g) == 2 and g[0] in "ACGT+-" and g[1] in "ACGT+-":
        return g[0], g[1]
    return None


def make_vcf_parser(ref, alts):
    """Parser for the GT field of one VCF record: '0/1', '1|1', './.', '1' (haploid)."""
    alleles = [ref] + alts
    cache = {}

    def parse(gt):
        if gt in cache:
            return cache[gt]
        out = None
        s = gt.split(":", 1)[0]
        parts = re.split(r"[/|]", s)
        if len(parts) == 1 and parts[0] not in (".", ""):
            parts = parts * 2          # haploid call -> homozygous
        if len(parts) == 2 and all(p.isdigit() and int(p) < len(alleles) for p in parts):
            out = alleles[int(parts[0])], alleles[int(parts[1])]
        cache[gt] = out
        return out

    return parse


# ─────────────────────────── coding of one marker ───────────────────────────

def code_calls(p1, p2, calls, parse):
    """Returns (codes (int8 array) | None, status). `p1`, `p2`, `calls`: raw values, `parse(raw)` -> pair | None."""
    g1, g2 = parse(p1), parse(p2)
    if g1 is None or g2 is None:
        return None, "Invalid parent"
    if g1[0] != g1[1] or g2[0] != g2[1]:
        return None, "Non-homozygous parent"
    a1, a2 = g1[0], g2[0]
    if a1 == a2:
        return None, "Non-segregating"
    hom1, hom2, het = (a1, a1), (a2, a2), {a1, a2}
    cmap = {}
    for v in set(calls):
        c = parse(v)
        if c is None:
            cmap[v] = -1
        elif c == hom1:
            cmap[v] = 0
        elif c == hom2:
            cmap[v] = 2
        elif set(c) == het:
            cmap[v] = 1
        else:
            cmap[v] = -1
    return np.fromiter((cmap[v] for v in calls), dtype=np.int8, count=len(calls)), "Informative"


# ─────────────────────────── chromosome / marker naming ───────────────────────────

def clean_chrom(chrom, style="number"):
    """'Chr01' -> '1', 'chr12' -> '12'; labels that are not numbered keep their text."""
    c = str(chrom).strip()
    if style != "number":
        return c
    m = re.fullmatch(r"(?i)(?:chr(?:omosome)?|lg|linkage[_ ]?group)?[_ ]*0*(\d+)", c)
    return m.group(1) if m else c


def marker_name(orig, chrom, pos, style):
    if style == "chrpos" and chrom != "" and pos != "":
        n = chrom if not str(chrom).isdigit() else str(chrom).zfill(2)
        return f"chr{n}_{pos}"
    if orig in ("", ".", None):
        n = chrom if not str(chrom).isdigit() else str(chrom).zfill(2)
        return f"chr{n}_{pos}" if chrom != "" and pos != "" else f"marker_{pos or orig}"
    return str(orig)


# ─────────────────────────── streaming sources ───────────────────────────

def _split(line, delim):
    return line.rstrip("\r\n").split(delim)


def iter_hapmap(path, nmax=None):
    """Yields ('header', samples) then one dict per marker: id, chrom, pos, ref, alt, calls (list of str)."""
    with _open_text(path) as fh:
        header = None
        delim = "\t"
        n = 0
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            if header is None:
                delim = "\t" if "\t" in line else ","
                header = [h.strip() for h in _split(line, delim)]
                low = [h.lower() for h in header]
                try:
                    i_chrom = next(i for i, h in enumerate(low) if h in ("chrom", "chr", "chromosome"))
                    i_pos = next(i for i, h in enumerate(low) if h in ("pos", "position", "bp"))
                except StopIteration:
                    raise ValueError("This does not look like a HapMap file: the header needs 'chrom' and 'pos' columns.")
                i_al = low.index("alleles") if "alleles" in low else -1
                fixed = {"rs#", "rs", "snp", "id", "alleles", "chrom", "chr", "chromosome", "pos", "position", "bp", "strand", "assembly#", "assembly",
                         "center", "protlsid", "assaylsid", "panellsid", "qccode"}
                first = 1 + max(i for i, h in enumerate(low) if h in fixed)   # the samples follow the last fixed HapMap column
                yield "header", header[first:]
                continue
            f = _split(line, delim)
            if len(f) < first:
                continue
            al = f[i_al].strip().upper().split("/") if i_al >= 0 else []
            yield {"id": f[0].strip(), "chrom": f[i_chrom].strip(), "pos": f[i_pos].strip(),
                   "ref": al[0] if len(al) == 2 else "", "alt": al[1] if len(al) == 2 else "", "calls": f[first:]}
            n += 1
            if nmax and n >= nmax:
                return


def iter_vcf(path, nmax=None):
    """Yields ('header', samples) then one dict per marker: id, chrom, pos, parser (alleles), calls (GT strings)."""
    with _open_text(path) as fh:
        samples = None
        n = 0
        for line in fh:
            if line.startswith("##") or not line.strip():
                continue
            if line.startswith("#"):
                f = _split(line, "\t")
                if len(f) < 10:
                    f = _split(line, " ")
                samples = [s.strip() for s in f[9:]]
                if not samples:
                    raise ValueError("The VCF has no sample columns.")
                yield "header", samples
                continue
            if samples is None:
                raise ValueError("The VCF has no '#CHROM' header line.")
            f = _split(line, "\t")
            if len(f) < 10:
                continue
            fmt = f[8].split(":")
            if "GT" not in fmt:
                continue
            gi = fmt.index("GT")
            calls = f[9:] if gi == 0 else [c.split(":")[gi] if c.count(":") >= gi else "." for c in f[9:]]
            alts = [a for a in f[4].split(",") if a not in (".", "")]
            yield {"id": f[2], "chrom": f[0], "pos": f[1], "ref": f[3], "alt": ",".join(alts),
                   "parse": make_vcf_parser(f[3], alts), "calls": calls}
            n += 1
            if nmax and n >= nmax:
                return


def _iter_for(fmt, path, nmax=None):
    if fmt == "hapmap":
        return iter_hapmap(path, nmax)
    if fmt == "vcf":
        return iter_vcf(path, nmax)
    raise ValueError(f"Unknown format {fmt!r}.")


# ─────────────────────────── raw (Excel / CSV) ───────────────────────────

def _read_raw(path, sheet=None, nrows=None):
    path = Path(path)
    if path.suffix.lower() in (".csv", ".txt", ".tsv"):
        sep = "\t" if path.suffix.lower() == ".tsv" else ","
        return pd.read_csv(path, sep=sep, nrows=nrows, dtype=str, keep_default_na=False)
    return pd.read_excel(path, sheet_name=sheet if sheet not in (None, "") else 0, nrows=nrows, dtype=str, keep_default_na=False)


def _sheet_names(path):
    if Path(path).suffix.lower() not in (".xlsx", ".xlsm", ".xls"):
        return []
    try:
        return list(pd.ExcelFile(path).sheet_names)
    except Exception:
        return []


def _guess_roles(columns):
    return {
        "marker": _guess(columns, "Marker", "Markers", "SNP", "SNP_ID", "ID", "Locus"),
        "chromosome": _guess(columns, "Chromosome", "Chr", "Chrom", "LG", "Linkage group"),
        "p1": _guess(columns, "P1", "Parent1", "Parent 1", "Parent_A", "PA", "Female"),
        "p2": _guess(columns, "P2", "Parent2", "Parent 2", "Parent_B", "PB", "Male"),
    }


def _best_sheet(path, sheets):
    """The sheet that looks most like raw genotype data (workbooks often start with a notes sheet)."""
    best, best_score = sheets[0], (-1, -1)
    for name in sheets:
        try:
            cols = [str(c) for c in _read_raw(path, name, nrows=3).columns]
        except Exception:
            continue
        score = (sum(bool(v) for v in _guess_roles(cols).values()), len(cols))
        if score > best_score:
            best, best_score = name, score
    return best


# ─────────────────────────── inspect ───────────────────────────

def inspect_file(path, fmt=None, sheet=None):
    """Format, samples / columns, guesses and a short preview, so the page can ask for the parents."""
    fmt = fmt if fmt in FORMATS else detect_format(path)
    out = {"status": "success", "format": fmt, "format_label": FORMAT_LABEL[fmt], "file": Path(path).name}
    if fmt == "raw":
        sheets = _sheet_names(path)
        if sheet in (None, "") and sheets:
            sheet = _best_sheet(path, sheets)
        df = _read_raw(path, sheet, nrows=200)
        columns = [str(c) for c in df.columns]
        guesses = _guess_roles(columns)
        if not guesses["marker"] and columns:
            guesses["marker"] = columns[0]
        sample = [v for c in columns[2:] for v in df[c].head(50).tolist() if str(v).strip().upper() not in MISSING_CALLS][:2000]
        best, best_hits = ":", -1
        for label, sep in SEPARATORS.items():
            hits = sum(parse_raw_call(v, sep) is not None for v in sample)
            if hits > best_hits:
                best, best_hits = label, hits
        out.update(sheets=sheets, sheet=sheet or "", columns=columns, guesses=guesses, separator=best,
                   preview_rows=df.head(8).to_dict(orient="records"), preview_columns=columns[:30])
        return out

    it = _iter_for(fmt, path, nmax=8)
    first = next(it)
    if first[0] != "header":
        raise ValueError("No header line found.")
    samples = first[1]
    rows = []
    chroms = []
    for rec in it:
        chroms.append(rec["chrom"])
        rows.append({"Marker": rec["id"], "Chromosome": rec["chrom"], "Position": rec["pos"],
                     "REF": rec["ref"], "ALT": rec["alt"], **{s: rec["calls"][i] for i, s in enumerate(samples[:6])}})
    p1, p2 = _guess_parents(samples)
    out.update(samples=samples, n_samples=len(samples), guesses={"p1": p1, "p2": p2},
               preview_rows=rows, preview_columns=["Marker", "Chromosome", "Position", "REF", "ALT", *samples[:6]],
               example_chromosomes=sorted(set(chroms))[:6])
    if not rows:
        raise ValueError("No marker rows were found in the file.")
    return out


# ─────────────────────────── convert ───────────────────────────

def convert(config):
    """`config` keys:
        input_file, format ('hapmap'|'vcf'|'raw'|'' = detect), output_dir
        HapMap / VCF: p1_sample, p2_sample, exclude_samples (list), marker_naming ('original'|'chrpos'),
                      chrom_style ('number'|'keep'), skip_unplaced (bool)
        Raw:          sheet, marker_col, chromosome_col, p1_col, p2_col, separator (':' '/' '|' 'none')
    Returns a JSON-serialisable summary with the output file paths."""
    src = Path(config["input_file"])
    fmt = config.get("format") or detect_format(src)
    if fmt not in FORMATS:
        raise ValueError(f"Unknown format {fmt!r}.")
    if fmt == "raw":
        res = _convert_raw(src, config)
    else:
        res = _convert_stream(src, fmt, config)
    return _write_outputs(src, fmt, config, res)


def _convert_stream(src, fmt, config):
    p1s, p2s = config.get("p1_sample", ""), config.get("p2_sample", "")
    exclude = set(config.get("exclude_samples") or [])
    naming = config.get("marker_naming", "original")
    style = config.get("chrom_style", "number")
    skip_unplaced = bool(config.get("skip_unplaced", False)) and style == "number"

    it = _iter_for(fmt, src)
    head = next(it)
    samples = head[1]
    if not p1s or not p2s:
        raise ValueError("Choose the two parents (Parent 1 and Parent 2) from the samples of the file.")
    if p1s == p2s:
        raise ValueError("Parent 1 and Parent 2 must be different samples.")
    for s in (p1s, p2s):
        if s not in samples:
            raise ValueError(f"Parent sample '{s}' is not in the file.")
    i1, i2 = samples.index(p1s), samples.index(p2s)
    off_idx = [i for i, s in enumerate(samples) if i not in (i1, i2) and s not in exclude]
    offspring = [samples[i] for i in off_idx]
    if not offspring:
        raise ValueError("No offspring samples are left after taking out the two parents.")

    coded, detail, excluded = [], [], []
    totals = np.zeros(4, dtype=np.int64)       # codes 0, 1, 2, -1
    total_markers = 0
    names_seen = {}
    for rec in it:
        total_markers += 1
        chrom = clean_chrom(rec["chrom"], style)
        pos = rec["pos"]
        name = marker_name(rec["id"], chrom, pos, naming)
        calls = rec["calls"]
        if len(calls) < len(samples):
            calls = calls + [""] * (len(samples) - len(calls))
        pa, pb = calls[i1], calls[i2]
        if fmt == "vcf":
            parse = rec["parse"]
        else:
            parse = parse_hapmap_call
        if skip_unplaced and not chrom.isdigit():
            excluded.append({"Marker": name, "Chromosome": chrom, "P1": pa, "P2": pb, "Reason": "Unplaced contig"})
            continue
        sel = [calls[i] for i in off_idx]
        codes, status = code_calls(pa, pb, sel, parse)
        if codes is None:
            excluded.append({"Marker": name, "Chromosome": chrom, "P1": pa, "P2": pb, "Reason": status})
            continue
        totals += np.array([(codes == 0).sum(), (codes == 1).sum(), (codes == 2).sum(), (codes == -1).sum()])
        names_seen[name] = names_seen.get(name, 0) + 1
        coded.append((name, chrom, codes))
        detail.append({"Marker": name, "Chromosome": chrom, "Position": pos, "P1": pa, "P2": pb,
                       "Missing_or_unexpected": int((codes == -1).sum()), "Status": status})
    return {"offspring": offspring, "coded": coded, "detail": detail, "excluded": excluded, "totals": totals,
            "total_markers": total_markers, "skipped": [], "duplicates": sum(1 for n in names_seen.values() if n > 1),
            "parents": (p1s, p2s), "sheet": ""}


def _convert_raw(src, config):
    sep = SEPARATORS.get(config.get("separator", ":"), ":")
    marker_col = config["marker_col"]
    chrom_col = config.get("chromosome_col") or ""
    p1_col, p2_col = config["p1_col"], config["p2_col"]
    style = config.get("chrom_style", "keep")
    df = _read_raw(src, config.get("sheet"))
    df.columns = [str(c) for c in df.columns]
    for role, col in (("Marker", marker_col), ("P1", p1_col), ("P2", p2_col)):
        if col not in df.columns:
            raise ValueError(f"{role} column '{col}' not found in the file.")
    if p1_col == p2_col:
        raise ValueError("P1 and P2 must be different columns.")
    if chrom_col and chrom_col not in df.columns:
        raise ValueError(f"Chromosome column '{chrom_col}' not found in the file.")
    meta = [c for c in (marker_col, chrom_col, p1_col, p2_col) if c]
    exclude = set(config.get("exclude_samples") or [])
    candidates = [c for c in df.columns if c not in meta and c not in exclude]

    def parse(v):
        return parse_raw_call(v, sep)

    offspring = [c for c in candidates if any(parse(v) is not None for v in df[c])]
    skipped = [c for c in candidates if c not in offspring]
    if not offspring:
        raise ValueError("No offspring genotype columns could be read. Check the allele separator (e.g. 'A:G' uses ':', 'AG' uses 'none').")

    coded, detail, excluded = [], [], []
    totals = np.zeros(4, dtype=np.int64)
    names_seen = {}
    for rec in df[meta + offspring].itertuples(index=False, name=None):
        row = dict(zip(meta + offspring, rec))
        name = row[marker_col]
        chrom = clean_chrom(row[chrom_col], style) if chrom_col else ""
        codes, status = code_calls(row[p1_col], row[p2_col], [row[c] for c in offspring], parse)
        if codes is None:
            excluded.append({"Marker": name, "Chromosome": chrom, "P1": row[p1_col], "P2": row[p2_col], "Reason": status})
            continue
        totals += np.array([(codes == 0).sum(), (codes == 1).sum(), (codes == 2).sum(), (codes == -1).sum()])
        names_seen[name] = names_seen.get(name, 0) + 1
        coded.append((name, chrom, codes))
        detail.append({"Marker": name, "Chromosome": chrom, "Position": "", "P1": row[p1_col], "P2": row[p2_col],
                       "Missing_or_unexpected": int((codes == -1).sum()), "Status": status})
    return {"offspring": offspring, "coded": coded, "detail": detail, "excluded": excluded, "totals": totals,
            "total_markers": len(df), "skipped": skipped, "duplicates": sum(1 for n in names_seen.values() if n > 1),
            "parents": (p1_col, p2_col), "sheet": config.get("sheet") or "(first)"}


def _write_outputs(src, fmt, config, r):
    offspring, coded = r["offspring"], r["coded"]
    cols = ["Marker", "Chromosome", *offspring]
    table = pd.DataFrame({"Marker": [c[0] for c in coded], "Chromosome": [c[1] for c in coded]})
    if coded:
        mat = np.vstack([c[2] for c in coded])
        table = pd.concat([table, pd.DataFrame(mat, columns=offspring)], axis=1)
    else:
        table = pd.DataFrame(columns=cols)
    # Linkage Mapping / Mendelian read a numeric chromosome column best when every label is a number
    if len(table) and table["Chromosome"].astype(str).str.fullmatch(r"\d+").all():
        table["Chromosome"] = table["Chromosome"].astype(int)
    excluded = pd.DataFrame(r["excluded"], columns=["Marker", "Chromosome", "P1", "P2", "Reason"])
    details = pd.DataFrame(r["detail"], columns=["Marker", "Chromosome", "Position", "P1", "P2", "Missing_or_unexpected", "Status"])
    reasons = {k: int((excluded["Reason"] == k).sum()) for k in REASONS}
    totals = r["totals"]
    summary = pd.DataFrame([
        ("Input file", src.name), ("Input format", FORMAT_LABEL[fmt]), ("Sheet", r["sheet"] or "-"),
        ("Parent 1", r["parents"][0]), ("Parent 2", r["parents"][1]),
        ("Markers in file", r["total_markers"]), ("Informative markers retained", len(table)), ("Markers left out", len(excluded)),
        *((f"  left out: {k}", n) for k, n in reasons.items() if n),
        ("Individuals (offspring)", len(offspring)),
        ("Columns skipped (no genotype calls)", ", ".join(r["skipped"]) or "none"),
        ("Coding", "0 = hom. P1, 1 = heterozygous, 2 = hom. P2, -1 = missing / unexpected"),
    ], columns=["Item", "Value"])

    out_dir = Path(config["output_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = re.sub(r"[^\w\-]+", "_", re.sub(r"(?i)\.(vcf|hmp)(\.txt)?(\.gz)?$|\.[^.]+$", "", src.name)).strip("_") or "genotypes"
    csv = out_dir / f"{stem}_coded.csv"
    xlsx = out_dir / f"{stem}_coded.xlsx"
    table.to_csv(csv, index=False)
    # An Excel sheet of millions of cells takes minutes to write and is slow to open: very large tables are given as CSV
    # (Linkage Mapping reads a CSV just as well) with a small report workbook beside it.
    big = table.shape[0] * max(table.shape[1], 1) > XLSX_MAX_CELLS
    if big:
        xlsx = out_dir / f"{stem}_report.xlsx"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
        if not big:
            table.to_excel(w, sheet_name="Coded_Markers", index=False)     # first sheet = the tool input
        details.head(100000).to_excel(w, sheet_name="Marker_Details", index=False)
        excluded.head(100000).to_excel(w, sheet_name="Excluded_Markers", index=False)
        summary.to_excel(w, sheet_name="Summary", index=False)
    coded_file = csv if big else xlsx

    n_calls = max(int(totals.sum()), 1)
    chrom_counts = table["Chromosome"].astype(str).value_counts().to_dict() if len(table) else {}
    return {
        "status": "success",
        "format": fmt,
        "coded_csv": str(csv),
        "coded_file": str(coded_file),
        "report_xlsx": str(xlsx),
        "big": bool(big),
        "total_markers": int(r["total_markers"]),
        "informative": int(len(table)),
        "excluded": int(len(excluded)),
        "reasons": reasons,
        "n_offspring": len(offspring),
        "skipped_columns": r["skipped"],
        "duplicate_names": int(r["duplicates"]),
        "code_counts": {"0": int(totals[0]), "1": int(totals[1]), "2": int(totals[2]), "-1": int(totals[3])},
        "missing_pct": round(100 * int(totals[3]) / n_calls, 2),
        "n_chromosomes": len(chrom_counts),
        "parents": list(r["parents"]),
        "preview_rows": table.head(15).to_dict(orient="records"),
        "preview_columns": ["Marker", "Chromosome", *offspring[:40]],
        "excluded_rows": excluded.head(500).to_dict(orient="records"),
    }


if __name__ == "__main__":
    # Command-line use:  python genofile.py <file> [parent1 parent2]
    import sys
    import tempfile
    info = inspect_file(sys.argv[1])
    cfg = {"input_file": sys.argv[1], "format": info["format"], "output_dir": tempfile.mkdtemp()}
    if info["format"] == "raw":
        g = info["guesses"]
        cfg.update(sheet=info["sheet"], marker_col=g["marker"], chromosome_col=g["chromosome"], p1_col=g["p1"], p2_col=g["p2"], separator=info["separator"])
    else:
        cfg.update(p1_sample=sys.argv[2] if len(sys.argv) > 2 else info["guesses"]["p1"],
                   p2_sample=sys.argv[3] if len(sys.argv) > 3 else info["guesses"]["p2"])
    res = convert(cfg)
    print(json.dumps({k: v for k, v in res.items() if not k.endswith("rows")}, indent=2))
