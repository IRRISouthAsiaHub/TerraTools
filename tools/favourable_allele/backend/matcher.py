"""
matcher.py — Favourable Allele Analysis: match a user's rice 4K haplopanel genotype file against the built-in list of
trait-gene markers (data/trait_genes_4k.json: 327 SNPs / indels in 166 trait genes with REF, ALT and the favourable allele).

The annotation never leaves the backend. The user uploads a genotype file; this module finds which of its marker IDs are in the
annotation and sends back only those markers (with their gene, trait and alleles) and the genotypes of the user's lines at them.

Accepted files
  * CSV / TSV / TXT: one row per marker, one column per line. The marker-ID column is found by matching (it need not be the
    first column); columns named chrom / pos / ref / alt are read when present and are not taken as lines.
  * VCF (also .vcf.gz): the ID column, or chrom_pos when the ID is ".", and the GT field of every sample.
  Genotype values: 0 / 1 / 2 (REF/REF, REF/ALT, ALT/ALT) with -1 or NA for no call; 0/0, 0/1, 1/1, ./. (also with |);
  nucleotides (A/G, AG, a single IUPAC letter).

Marker IDs are matched in this order: the complete ID (chr1_7268879_LOC_Os01g13060), chromosome_position (chr1_7268879),
chromosome and the position of the annotation.
"""

import csv
import gzip
import io
import json
import re
from functools import lru_cache
from pathlib import Path

DATA = Path(__file__).resolve().parent / "data" / "trait_genes_4k.json"

ID_NAMES = {"marker", "marker_id", "markerid", "marker_name", "snp", "snp_id", "snpid", "id", "rs", "rsid", "variant", "variant_id", "marker_in_8k", "marker_in_4k", "locus"}
CHROM_NAMES = {"chrom", "chr", "chromosome", "#chrom", "contig", "chromosome_no"}
POS_NAMES = {"pos", "position", "bp", "pos_bp", "physical_position", "position_bp", "bp_position"}
REF_NAMES = {"ref", "reference", "ref_allele"}
ALT_NAMES = {"alt", "alternate", "alt_allele", "alternative"}
MISSING = {"", "na", "nan", "n/a", "none", ".", "./.", ".|.", "nn", "n", "-", "--", "-1", "./", "?", "??", "-1/-1", "null", "missing"}
COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}
IUPAC = {"R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC"}


def _revcomp(s: str) -> str:
    return "".join(COMP.get(c, "N") for c in reversed(s.upper()))


@lru_cache(maxsize=1)
def annotation() -> dict:
    raw = json.loads(DATA.read_text(encoding="utf-8"))
    snps = []
    for i, s in enumerate(raw["snps"]):
        s = dict(s)
        s["i"] = i
        s["label"] = s.get("name") or s.get("gene")
        snps.append(s)
    return {"snps": snps, "lost": raw.get("lostGenes", []), "name": raw.get("name", "")}


def info() -> dict:
    a = annotation()
    genes = {s["label"] for s in a["snps"]}
    lost = set(a["lost"])
    return {"panel": "Rice 4K haplopanel", "n_snps": len(a["snps"]), "n_genes": len(genes), "lost_genes": sorted(lost)}


def chrom_number(tok) -> str | None:
    t = re.sub(r"^(chromosome|chrom|chr)[_\-\s]*", "", str(tok).strip(), flags=re.I)
    if re.fullmatch(r"\d+", t):
        return str(int(t))
    return t.upper() or None


def key_chrpos(marker: str):
    """('1', '7268879') from chr1_7268879, chr01_7268879_LOC_Os01g13060 or 1_7268879; None when the ID is not of that form."""
    m = re.match(r"^(?:chromosome|chrom|chr)?[_\-\s]*0*(\d+)[_:\-.](\d+)(?:[_:\-.].*)?$", marker.strip(), flags=re.I)
    return (m.group(1), m.group(2)) if m else None


@lru_cache(maxsize=1)
def _indexes():
    exact, chrpos, position = {}, {}, {}
    for s in annotation()["snps"]:
        exact.setdefault(s["marker"].lower(), []).append(s["i"])
        k = key_chrpos(s["marker"])
        if k:
            chrpos.setdefault(k, []).append(s["i"])
        position.setdefault((chrom_number(s["chrom"]), str(s["pos"])), []).append(s["i"])
    return exact, chrpos, position


def match_id(marker: str, chrom=None, pos=None):
    """([annotation indexes], method) for one file marker; ([], None) when it is not in the annotation."""
    exact, chrpos, position = _indexes()
    m = (marker or "").strip()
    if m.lower() in exact:
        return exact[m.lower()], "complete ID"
    k = key_chrpos(m)
    if k and k in chrpos:
        return chrpos[k], "chromosome_position"
    if chrom not in (None, "") and pos not in (None, ""):
        p = re.sub(r"[^\d]", "", str(pos))
        key = (chrom_number(chrom), str(int(p)) if p else "")
        if key in position:
            return position[key], "chromosome and position"
    if k and (k[0], k[1]) in position:
        return position[(k[0], k[1])], "chromosome and position"
    return [], None


# ───────────────────────── reading the file ─────────────────────────
def _open_text(path: str):
    p = Path(path)
    raw = p.read_bytes()
    if p.suffix.lower() == ".gz" or raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    for enc in ("utf-8-sig", "latin-1"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode("utf-8", errors="replace")


def read_table(path: str):
    """(header, rows) of a delimited file or a VCF. A VCF is turned into the columns ID, #CHROM, POS, REF, ALT and the samples (GT only)."""
    text = _open_text(path)
    lines = text.splitlines()
    vcf = any(l.startswith("##fileformat=VCF") for l in lines[:5]) or str(path).lower().endswith((".vcf", ".vcf.gz"))
    if vcf:
        header, rows = None, []
        for l in lines:
            if l.startswith("##") or not l.strip():
                continue
            parts = l.rstrip("\n").split("\t")
            if l.startswith("#CHROM"):
                header = parts
                samples = parts[9:]
                continue
            if header is None or len(parts) < 10:
                continue
            ident = parts[2] if parts[2] not in (".", "") else f"{parts[0]}_{parts[1]}"
            gts = [p.split(":")[0] for p in parts[9:]]
            rows.append([ident, parts[0], parts[1], parts[3], parts[4]] + gts)
        if header is None:
            raise ValueError("No #CHROM header line was found in the VCF.")
        return ["ID", "#CHROM", "POS", "REF", "ALT"] + header[9:], rows
    sample = "\n".join(lines[:20])
    delim = max([",", "\t", ";"], key=lambda d: sample.count(d))
    if sample.count(delim) == 0:
        delim = None
    rdr = csv.reader(io.StringIO(text), delimiter=delim or ",")
    rows = [r for r in rdr if any(c.strip() for c in r)]
    if len(rows) < 2:
        raise ValueError("The file has no data rows.")
    return [c.strip() for c in rows[0]], rows[1:]


def decode(v: str, ref: str, alt: str) -> str:
    """One genotype as '0', '1', '2' or 'x' (no call), read against the annotation's REF and ALT."""
    s = (v or "").strip()
    if s.lower() in MISSING:
        return "x"
    if re.fullmatch(r"-?\d+(\.0+)?", s):
        n = int(float(s))
        return str(n) if n in (0, 1, 2) else "x"
    parts = re.split(r"[/|]", s)
    if len(parts) == 1:
        u = s.upper()
        if u in IUPAC:
            parts = list(IUPAC[u])
        elif len(u) == 2 and u.isalpha():
            parts = [u[0], u[1]]
        elif len(u) == 1 and u.isalpha():
            parts = [u, u]
        else:
            return "x"
    if len(parts) != 2:
        return "x"
    if all(p in ("0", "1") for p in parts):
        return str(parts[0].count("1") + parts[1].count("1"))
    if any(p in (".", "") for p in parts):
        return "x"
    r, a = (ref or "").upper(), (alt or "").upper()

    def which(p, flip):
        p = p.upper()
        rr, aa = (_revcomp(r), _revcomp(a)) if flip and r not in ("-", "") and a not in ("-", "") else (r, a)
        return 0 if p == rr else 1 if p == aa else None

    for flip in (False, True):
        w = [which(p, flip) for p in parts]
        if None not in w:
            return str(sum(w))
    return "x"


def _allele_relation(f_ref, f_alt, ref, alt):
    """'same', 'swapped' or None for the file's REF/ALT against the annotation's (also on the opposite strand)."""
    fr, fa, r, a = (x.strip().upper() for x in (f_ref, f_alt, ref, alt))
    if (fr, fa) == (r, a):
        return "same"
    if (fr, fa) == (a, r):
        return "swapped"
    rr, ar = _revcomp(r), _revcomp(a)
    if (fr, fa) == (rr, ar):
        return "same"
    if (fr, fa) == (ar, rr):
        return "swapped"
    return None


def load(path: str) -> dict:
    """Read the user's file, match it against the annotation and return the data the page works on."""
    try:
        if not path or not Path(path).is_file():
            return {"status": "error", "message": "Please choose the genotype file."}
        header, rows = read_table(path)
        ann = annotation()
        snps = ann["snps"]
        names = [h.strip().lower() for h in header]
        ncol = len(header)
        rows = [r + [""] * (ncol - len(r)) if len(r) < ncol else r[:ncol] for r in rows]

        # marker-ID column: the one among the first columns (or named like one) whose values match the annotation best
        cand = list(range(min(5, ncol)))
        cand += [i for i, n in enumerate(names) if n in ID_NAMES and i not in cand]
        best, best_n = None, -1
        for c in cand:
            n = sum(1 for r in rows if match_id(r[c])[0])
            if n > best_n or (n == best_n and best is not None and names[c] in ID_NAMES and names[best] not in ID_NAMES):
                best, best_n = c, n
        id_col = best if best is not None else 0

        def col_of(options):
            for i, n in enumerate(names):
                if i != id_col and n in options:
                    return i
            return None

        c_chrom, c_pos, c_ref, c_alt = col_of(CHROM_NAMES), col_of(POS_NAMES), col_of(REF_NAMES), col_of(ALT_NAMES)
        meta = {id_col, c_chrom, c_pos, c_ref, c_alt} | set(range(id_col))
        sample_cols = [i for i in range(ncol) if i not in meta and names[i] not in {"qual", "filter", "info", "format", "gene", "gene_id", "trait", "cm", "gen_pos"}]
        if not sample_cols:
            return {"status": "error", "message": "No genotype columns (lines) were found next to the marker-ID column."}
        if best_n <= 0:
            return {"status": "error", "message": ("None of the marker IDs in the file matched the trait-gene markers of the 4K haplopanel. "
                                                    f"The IDs were read from the column “{header[id_col]}” (for example {rows[0][id_col]!r}). "
                                                    "The IDs must look like chr1_7268879 (chromosome_position), optionally followed by the gene (chr1_7268879_LOC_Os01g13060); "
                                                    "alternatively give chromosome and position columns.")}

        matched = {}          # annotation index -> (file row, method)
        report, mism = [], []
        file_dups, matched_rows = set(), set()
        for ri, r in enumerate(rows):
            idxs, method = match_id(r[id_col], r[c_chrom] if c_chrom is not None else None, r[c_pos] if c_pos is not None else None)
            for ai in idxs:
                s = snps[ai]
                rel = "same"
                if c_ref is not None and c_alt is not None and r[c_ref].strip() and r[c_alt].strip():
                    rel = _allele_relation(r[c_ref], r[c_alt], s["ref"] or "", s["alt"] or "")
                if rel is None:
                    mism.append({"file_id": r[id_col], "marker": s["marker"], "gene": s["label"],
                                 "file_alleles": f"{r[c_ref]}/{r[c_alt]}", "panel_alleles": f"{s['ref']}/{s['alt']}"})
                    continue
                if ai in matched:
                    file_dups.add(r[id_col])
                    continue
                matched[ai] = (ri, rel)
                matched_rows.add(ri)
                report.append({"file_id": r[id_col], "marker": s["marker"], "gene": s["label"], "trait": s["trait"] or "Unspecified", "method": method,
                               "ref": s["ref"], "alt": s["alt"], "favourable": s["fav"]})
        if not matched:
            return {"status": "error", "message": "The marker IDs matched, but the REF/ALT alleles in the file disagree with the annotation for every one of them."}

        order = sorted(matched)
        used = [snps[i] for i in order]
        lines = [{"id": f"L{j + 1}", "name": header[c].strip() or f"line {j + 1}", "panel": "", "code": ""} for j, c in enumerate(sample_cols)]
        calls = []
        for c in sample_cols:
            chars = []
            for ai in order:
                ri, rel = matched[ai]
                s = snps[ai]
                g = decode(rows[ri][c], s["ref"] or "", s["alt"] or "")
                if rel == "swapped" and g in "02":
                    g = "2" if g == "0" else "0"
                chars.append(g)
            calls.append("".join(chars))

        genes_all = {}
        for s in snps:
            genes_all.setdefault(s["label"], {"n": 0, "found": 0})["n"] += 1
        for s in used:
            genes_all[s["label"]]["found"] += 1
        no_gene = sorted(g for g, v in genes_all.items() if v["found"] == 0)
        lost = set(ann["lost"])
        called = sum(1 for cs in calls for ch in cs if ch != "x")
        return {
            "status": "success",
            "file": Path(path).name,
            "snps": [{k: s[k] for k in ("gene", "name", "trait", "sub", "marker", "chrom", "pos", "ref", "alt", "fav", "unfav")} for s in used],
            "varieties": lines,
            "calls": calls,
            "lostGenes": sorted(lost),
            "match": {
                "file_markers": len(rows), "id_column": header[id_col], "matched_markers": len(matched), "panel_markers": len(snps),
                "unmatched_file_markers": len(rows) - len(matched_rows),
                "genes_total": len(genes_all), "genes_found": len(genes_all) - len(no_gene), "genes_missing": no_gene,
                "genes_partial": sorted((g, v["found"], v["n"]) for g, v in genes_all.items() if 0 < v["found"] < v["n"]),
                "lines": len(lines), "calls_total": len(lines) * len(matched), "calls_made": called,
                "by_method": {m: sum(1 for r in report if r["method"] == m) for m in sorted({r["method"] for r in report})},
                "allele_mismatch": mism[:50], "n_allele_mismatch": len(mism), "duplicates": sorted(file_dups)[:20],
            },
            "report": report,
        }
    except Exception as e:  # noqa: BLE001
        import traceback
        return {"status": "error", "message": f"The file could not be read: {e}\n{traceback.format_exc()[-400:]}"}
