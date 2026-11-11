"""
VCF input: reading, validation, and the SNP table of the two bulks.

Bulks are pooled DNA, so their GT field (if any) is NOT an ordinary diploid genotype and is never used: the allele
depths (AD) are the data. GT / GQ are used only for the parents, which are inbred individuals.
"""

from __future__ import annotations

import gzip
import io
from pathlib import Path

import numpy as np
import pandas as pd

from tools.qtl_seq.backend.calling import homozygous_call

BASES = ["A", "C", "G", "T"]
_COLS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
_EMPTY = ["Chromosome", "Position", "Ref", "Alt", "Qual", "Filter", "AD_Ref", "AD_Alt", "DP", "GT", "GQ"]


class VcfError(Exception):
    pass


def _open_text(path):
    p = str(path)
    if p.lower().endswith((".gz", ".bgz")):
        return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8", errors="replace")
    return open(p, "r", encoding="utf-8", errors="replace")


def read_header(path):
    """{'samples': [...], 'contigs': {name: length or None}, 'n_header_lines': int}"""
    contigs, samples, n = {}, [], 0
    try:
        with _open_text(path) as fh:
            for line in fh:
                n += 1
                if line.startswith("##contig"):
                    body = line[line.find("<") + 1:line.rfind(">")]
                    d = dict(kv.split("=", 1) for kv in body.split(",") if "=" in kv)
                    if "ID" in d:
                        contigs[d["ID"]] = int(d["length"]) if d.get("length", "").isdigit() else None
                elif line.startswith("#CHROM"):
                    samples = line.rstrip("\r\n").split("\t")[9:]
                    break
                elif not line.startswith("#"):
                    raise VcfError(f"{Path(path).name}: no #CHROM header line before the first record.")
    except (OSError, EOFError) as e:
        raise VcfError(f"{Path(path).name}: cannot be read as a (gzipped) VCF ({e}).") from None
    if not samples:
        raise VcfError(f"{Path(path).name}: the VCF has no sample columns, so there are no allele depths to use.")
    return {"samples": samples, "contigs": contigs, "n_header_lines": n}


def read_vcf(path, sample="", role="bulk", max_records=None, chunk=400_000):
    """Read the biallelic SNPs of one sample.

    Returns (DataFrame, info). DataFrame columns: Chromosome, Position, Ref, Alt, Qual, Filter, AD_Ref, AD_Alt, DP,
    GT, GQ (NaN / '' where the VCF does not give them). Records that cannot be used are counted in `info`:
    indels, multiallelic sites, missing allele depths, duplicate positions."""
    header = read_header(path)
    samples = header["samples"]
    if sample and sample not in samples:
        raise VcfError(f"{Path(path).name}: the sample '{sample}' ({role}) is not in the file. Samples: {', '.join(samples)}.")
    name = sample or samples[0]
    info = {"file": Path(path).name, "samples": samples, "sample": name, "contigs": header["contigs"], "records": 0,
            "snps": 0, "indels": 0, "multiallelic": 0, "missing_ad": 0, "duplicates": 0, "no_alt": 0,
            "has_ad": False, "has_gt": False, "has_gq": False, "has_dp": False}
    parts = []
    reader = pd.read_csv(path, sep="\t", header=None, names=_COLS + samples, skiprows=header["n_header_lines"],
                         dtype=str, chunksize=chunk, usecols=list(range(9 + len(samples))), na_filter=False,
                         compression="infer")
    for ch in reader:
        info["records"] += len(ch)
        ref, alt = ch["REF"].str.upper(), ch["ALT"].str.upper()
        no_alt = alt.isin([".", ""])
        info["no_alt"] += int(no_alt.sum())
        multi = alt.str.contains(",", regex=False)
        snp = (ref.str.len() == 1) & (alt.str.len() == 1) & ref.isin(BASES) & alt.isin(BASES)
        info["multiallelic"] += int((multi & (ref.str.len() == 1)).sum())
        info["indels"] += int((~snp & ~multi & ~no_alt).sum())
        ch = ch[snp]
        if ch.empty:
            continue
        out = pd.DataFrame({"Chromosome": ch["CHROM"], "Position": pd.to_numeric(ch["POS"]).astype(np.int64),
                            "Ref": ref[snp], "Alt": alt[snp], "Qual": pd.to_numeric(ch["QUAL"], errors="coerce"),
                            "Filter": ch["FILTER"]})
        for c in ("AD_Ref", "AD_Alt", "DP", "GQ"):
            out[c] = np.nan
        out["GT"] = ""
        smp = ch[name]
        for fmt, labels in ch.groupby("FORMAT").groups.items():
            keys = fmt.split(":")
            cols = smp.loc[labels].str.split(":", expand=True)
            rows = out.index.get_indexer(labels)

            def put(col, values):
                out.iloc[rows, out.columns.get_loc(col)] = values

            if "AD" in keys and keys.index("AD") < cols.shape[1]:
                ad = cols[keys.index("AD")].str.split(",", expand=True)
                if ad.shape[1] >= 2:
                    put("AD_Ref", pd.to_numeric(ad[0], errors="coerce").to_numpy())
                    put("AD_Alt", pd.to_numeric(ad[1], errors="coerce").to_numpy())
                info["has_ad"] = True
            for key, col in (("DP", "DP"), ("GQ", "GQ")):
                if key in keys and keys.index(key) < cols.shape[1]:
                    put(col, pd.to_numeric(cols[keys.index(key)], errors="coerce").to_numpy())
                    info["has_" + key.lower()] = True
            if "GT" in keys and keys.index("GT") < cols.shape[1]:
                put("GT", cols[keys.index("GT")].fillna("").to_numpy())
                info["has_gt"] = True
        parts.append(out)
        if max_records and info["records"] >= max_records:
            break
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=_EMPTY)
    dup = df.duplicated(["Chromosome", "Position"], keep="first")
    info["duplicates"] = int(dup.sum())
    df = df[~dup].reset_index(drop=True)
    info["missing_ad"] = int((df["AD_Ref"].isna() | df["AD_Alt"].isna()).sum())
    info["snps"] = int(len(df))
    return df, info


# ============================================================
# BULKS
# ============================================================

def pair_bulks(high, low):
    """Both bulks at the same positions. Returns (table, stats).

    table columns: Chromosome, Position, Ref, Alt (VCF alleles), Ref_High, Alt_High, Ref_Low, Alt_Low, Qual and the
    FILTER value of each file. A position is used only when REF and ALT agree in the two files and both give AD."""
    stats = {}
    m = high.merge(low, on=["Chromosome", "Position"], how="outer", suffixes=("_H", "_L"), indicator=True)
    stats["high_only"] = int((m["_merge"] == "left_only").sum())
    stats["low_only"] = int((m["_merge"] == "right_only").sum())
    m = m[m["_merge"] == "both"]
    bad = (m["Ref_H"] != m["Ref_L"]) | (m["Alt_H"] != m["Alt_L"])
    stats["allele_mismatch"] = int(bad.sum())
    m = m[~bad]
    miss = m[["AD_Ref_H", "AD_Alt_H", "AD_Ref_L", "AD_Alt_L"]].isna().any(axis=1)
    stats["missing_ad"] = int(miss.sum())
    m = m[~miss]
    out = pd.DataFrame({
        "Chromosome": m["Chromosome"].to_numpy(), "Position": m["Position"].to_numpy(),
        "Ref": m["Ref_H"].to_numpy(), "Alt": m["Alt_H"].to_numpy(),
        "Ref_High": m["AD_Ref_H"].to_numpy(np.int64), "Alt_High": m["AD_Alt_H"].to_numpy(np.int64),
        "Ref_Low": m["AD_Ref_L"].to_numpy(np.int64), "Alt_Low": m["AD_Alt_L"].to_numpy(np.int64),
        "Qual": np.fmin(m["Qual_H"].to_numpy(float), m["Qual_L"].to_numpy(float)),
        "Filter_H": m["Filter_H"].to_numpy(), "Filter_L": m["Filter_L"].to_numpy(),
    })
    return out.reset_index(drop=True), stats


# ============================================================
# PARENTS
# ============================================================

_HOM_REF = ["0/0", "0|0", "0"]
_HOM_ALT = ["1/1", "1|1", "1"]


def parent_genotypes(df, min_depth=7, min_gq=30, error=0.01):
    """Homozygous calls of a parent: Class 'ref' / 'alt' ('' for het or missing), Depth, GQ, Pass.

    The VCF's own GT and GQ are used when present; without a GT the genotype is called from the allele depths.
    Where GQ is missing it is computed from the depths when they exist (otherwise it is not required)."""
    out = df[["Chromosome", "Position", "Ref", "Alt"]].copy()
    ad = df[["AD_Ref", "AD_Alt"]].to_numpy(float)
    have_ad = ~np.isnan(ad).any(axis=1)
    depth = np.where(have_ad, np.nansum(ad, axis=1), df["DP"].to_numpy(float))
    gt = df["GT"].fillna("").astype(str).to_numpy()
    cls = np.where(np.isin(gt, _HOM_REF), "ref", np.where(np.isin(gt, _HOM_ALT), "alt", "")).astype(object)
    gq = df["GQ"].to_numpy(float).copy()
    nogt = np.isin(gt, ["", ".", "./.", ".|."])
    sel = nogt & have_ad
    if sel.any():
        top, _d, q = homozygous_call(ad[sel].T, error)
        cls[sel] = np.where(top == 0, "ref", "alt")
        gq[sel] = q
    fill = np.isnan(gq) & have_ad
    if fill.any():
        _t, _d, q = homozygous_call(ad[fill].T, error)
        gq[fill] = q
    out["Class"], out["Depth"], out["GQ"] = cls, depth, gq
    out["Pass"] = (cls != "") & (np.nan_to_num(depth, nan=0.0) >= min_depth) & ~(gq < min_gq)
    return out


def index_is_alt(known_parent, parent_class, index_parent):
    """True where the index allele (the allele of `index_parent`) is the VCF ALT.

    Known parent == index parent: its own allele. Otherwise: the allele the known parent does not carry."""
    own_alt = np.asarray(parent_class) == "alt"
    return own_alt if known_parent == index_parent else ~own_alt


def orient(table, is_alt):
    """Swap Ref/Alt (and their counts) where the index allele is the VCF REF, so that Alt is always the index allele."""
    t = table.copy()
    swap = ~np.asarray(is_alt, dtype=bool)
    for a, b in (("Ref", "Alt"), ("Ref_High", "Alt_High"), ("Ref_Low", "Alt_Low")):
        x, y = t[a].to_numpy().copy(), t[b].to_numpy().copy()
        t[a], t[b] = np.where(swap, y, x), np.where(swap, x, y)
    return t


def parent_samples(header, cfg):
    """(parent 1 sample, parent 2 sample) of the parental VCF; None where that parent is not in the file."""
    names = header["samples"]
    s1, s2 = cfg.get("parent1_sample", ""), cfg.get("parent2_sample", "")
    for s in (s1, s2):
        if s and s not in names:
            raise VcfError(f"The parental VCF has no sample '{s}'. Samples: {', '.join(names)}.")
    if not s1 and not s2:
        return names[0], (names[1] if len(names) > 1 else None)
    return (s1 or None), (s2 or None)


def _parents(cfg, log):
    path = cfg.get("parent_vcf") or ""
    if not path:
        return None, None
    if not Path(path).is_file():
        raise VcfError(f"The parental VCF was not found: {path}")
    header = read_header(path)
    s1, s2 = parent_samples(header, cfg)
    out = []
    for which, s in (("parent1", s1), ("parent2", s2)):
        if s is None:
            out.append(None)
            continue
        df, info = read_vcf(path, s, which)
        if not (info["has_gt"] or info["has_ad"]):
            raise VcfError("The parental VCF has neither GT nor AD, so the parents' genotypes cannot be determined.")
        log(f"Parental VCF: {which} = sample '{info['sample']}' ({info['snps']:,} SNPs)")
        out.append(df)
    return out[0], out[1]


def _lookup(table, parent, tag, cfg):
    """Merge a parent's class / pass flag onto the bulk table (same order)."""
    g = parent_genotypes(parent, cfg["min_parent_depth"], cfg["min_parent_gq"], cfg["base_error"])
    g = g[["Chromosome", "Position", "Ref", "Alt", "Class", "Pass"]].rename(
        columns={"Ref": "R" + tag, "Alt": "A" + tag, "Class": "C" + tag, "Pass": "P" + tag})
    return table[["Chromosome", "Position", "Ref", "Alt"]].merge(g, on=["Chromosome", "Position"], how="left")


def build_vcf_snps(cfg, log=print, cancelled=None):
    """VCF mode. Returns (snps, ledger, chrom_lengths, warnings, input info).

    snps has the SNP-table columns with Alt = the index allele."""
    p = cfg
    warnings, ledger = [], {}
    for key, label in (("high_vcf", "high-bulk VCF"), ("low_vcf", "low-bulk VCF")):
        if not p.get(key) or not Path(p[key]).is_file():
            raise VcfError(f"The {label} was not found: {p.get(key) or '(not selected)'}")
    log(f"Reading the high-bulk VCF {Path(p['high_vcf']).name}")
    high, ih = read_vcf(p["high_vcf"], p.get("high_sample", ""), "high bulk")
    log(f"Reading the low-bulk VCF {Path(p['low_vcf']).name}")
    low, il = read_vcf(p["low_vcf"], p.get("low_sample", ""), "low bulk")
    if cancelled and cancelled():
        raise InterruptedError("cancelled")
    if Path(p["high_vcf"]).resolve() == Path(p["low_vcf"]).resolve() and ih["sample"] == il["sample"]:
        raise VcfError("The high and low bulk are the same sample of the same file. Choose two different samples.")
    for info, label in ((ih, "high-bulk"), (il, "low-bulk")):
        if not info["has_ad"]:
            raise VcfError(f"The {label} VCF has no AD (allele depth) field. The analysis needs the number of reads of each "
                           "allele; call the variants again keeping AD (for example bcftools mpileup -a AD).")
        for k, txt in (("multiallelic", "multiallelic sites"), ("indels", "indel records"),
                       ("duplicates", "duplicate positions"), ("missing_ad", "SNPs with a missing allele depth")):
            if info[k]:
                warnings.append(f"{label} VCF: {info[k]:,} {txt} were left out.")
    ca, cb = set(high["Chromosome"]), set(low["Chromosome"])
    if not ca & cb:
        raise VcfError("The two bulk VCFs share no chromosome names (for example 'chr1' vs '1'). They must use the same reference.")
    if ca ^ cb:
        warnings.append("Chromosomes present in only one bulk VCF: " + ", ".join(sorted(ca ^ cb)[:10]) + ".")
    table, ps = pair_bulks(high, low)
    if ps["high_only"] or ps["low_only"]:
        warnings.append(f"{ps['high_only']:,} positions occur only in the high-bulk VCF and {ps['low_only']:,} only in the "
                        "low-bulk VCF. Only positions present in both are compared (a jointly called VCF keeps more).")
    if ps["allele_mismatch"]:
        warnings.append(f"{ps['allele_mismatch']:,} positions have different REF/ALT alleles in the two bulks and were removed.")
    ledger["total_snps"] = int(len(table))
    if table.empty:
        raise VcfError("The two bulk VCFs have no SNP position in common that carries allele depths.")

    keep = np.ones(len(table), dtype=bool)
    q_bad = np.zeros(len(table), dtype=bool)
    if p.get("min_qual", 0):
        q_bad |= ~(table["Qual"].to_numpy(float) >= p["min_qual"])
    if p.get("require_pass"):
        bad_f = lambda s: ~np.isin(s, ["PASS", ".", ""])
        q_bad |= bad_f(table["Filter_H"].to_numpy()) | bad_f(table["Filter_L"].to_numpy())

    par1, par2 = _parents(p, log)
    index_parent = p.get("index_parent", "parent2")
    mode = "two_parents" if (par1 is not None and par2 is not None) else (
        "one_parent" if (par1 is not None or par2 is not None) else "no_parent")
    ledger["mode"] = mode
    is_alt = np.ones(len(table), dtype=bool)
    ledger["removed_parental_quality"] = 0
    if mode == "no_parent":
        warnings.append("No parental data: every SNP of the bulk VCFs is used and the index allele is the VCF ALT allele. "
                        "Parental polymorphism cannot be checked, so sites that are not polymorphic between the parents "
                        "(and sequencing errors) can remain, and the null model assumes the ALT allele has frequency 1/2 "
                        "in the population.")
        ledger["removed_parental"] = 0
    elif mode == "one_parent":
        which, par = ("parent1", par1) if par1 is not None else ("parent2", par2)
        warnings.append(f"Only {which} was supplied: the other parent is assumed to carry the other allele, so parental "
                        "polymorphism cannot be completely validated.")
        m = _lookup(table, par, "1", p)
        found = m["C1"].notna().to_numpy()
        same = ((m["R1"] == m["Ref"]) & (m["A1"] == m["Alt"])).to_numpy()
        good = found & same & (m["P1"] == True).to_numpy()
        ledger["removed_parental"] = int((~good & ~q_bad).sum())
        keep &= good
        is_alt = index_is_alt(which, m["C1"].fillna("").to_numpy(), index_parent)
    else:
        m1, m2 = _lookup(table, par1, "1", p), _lookup(table, par2, "2", p)
        c1, c2 = m1["C1"].fillna("").to_numpy(), m2["C2"].fillna("").to_numpy()
        same = ((m1["R1"] == m1["Ref"]) & (m1["A1"] == m1["Alt"]) & (m2["R2"] == m2["Ref"]) & (m2["A2"] == m2["Alt"])).to_numpy()
        good = ((m1["P1"] == True).to_numpy() & (m2["P2"] == True).to_numpy() & (c1 != "") & (c2 != "")
                & (c1 != c2) & same)
        ledger["removed_parental"] = int((~good & ~q_bad).sum())
        keep &= good
        is_alt = (c2 if index_parent == "parent2" else c1) == "alt"
    ledger["removed_quality"] = int(q_bad.sum())
    keep &= ~q_bad
    table = orient(table, is_alt)[keep].reset_index(drop=True)
    ledger["passing_parental"] = int(len(table))
    present = set(table["Chromosome"])
    lengths = {}
    for info in (ih, il):
        for c, ln in info["contigs"].items():
            if ln and c in present:
                lengths[c] = max(lengths.get(c, 0), ln)
    keep_cols = ["Chromosome", "Position", "Ref", "Alt", "Ref_High", "Alt_High", "Ref_Low", "Alt_Low"]
    return table[keep_cols], ledger, lengths, warnings, {"high": _brief(ih), "low": _brief(il)}


def _brief(info):
    return {k: v for k, v in info.items() if k != "contigs"}


def validate_vcf_inputs(cfg, fasta=None, sample_records=200_000):
    """Quick check for the "Check inputs" button: reads the first `sample_records` records of each file.

    Returns {'errors': [...], 'warnings': [...], 'files': {...}}. Checks: readable, samples, AD present, chromosome
    names, duplicate positions, multiallelic sites, missing allele depths, REF/ALT agreement between the files,
    REF against the reference FASTA (when given), read depth."""
    errors, warnings, files = [], [], {}
    data = {}
    for key, role in (("high_vcf", "high bulk"), ("low_vcf", "low bulk")):
        path = cfg.get(key)
        if not path:
            errors.append(f"The {role} VCF is not selected.")
            continue
        try:
            df, info = read_vcf(path, cfg.get(key.replace("_vcf", "_sample"), ""), role, max_records=sample_records)
        except (VcfError, ValueError, OSError) as e:
            errors.append(str(e))
            continue
        data[key] = df
        files[role] = _brief(info)
        if not info["has_ad"]:
            errors.append(f"{role} VCF ({info['file']}): no AD field, so allele depths are not available.")
        for k, txt in (("multiallelic", "multiallelic sites"), ("indels", "indel records"),
                       ("duplicates", "duplicate positions"), ("missing_ad", "SNPs without allele depth")):
            if info[k]:
                warnings.append(f"{role} VCF: {info[k]:,} {txt} in the first {info['records']:,} records (they are left out).")
        if info["snps"]:
            depth = (df["AD_Ref"] + df["AD_Alt"]).dropna()
            if len(depth):
                med = float(depth.median())
                files[role]["median_depth"] = med
                if med < cfg.get("min_depth", 7):
                    warnings.append(f"{role} VCF: median read depth {med:.0f} is below the minimum depth "
                                    f"({cfg.get('min_depth', 7)}); most SNPs would be removed.")
        else:
            errors.append(f"{role} VCF ({info['file']}): no biallelic SNP with the required fields was found.")
    if len(data) == 2:
        a, b = data["high_vcf"], data["low_vcf"]
        if not set(a["Chromosome"]) & set(b["Chromosome"]):
            errors.append("The two bulk VCFs share no chromosome names (for example 'chr1' vs '1').")
        else:
            t, ps = pair_bulks(a, b)
            if ps["allele_mismatch"]:
                warnings.append(f"{ps['allele_mismatch']:,} shared positions have different REF/ALT in the two files.")
            if len(t) == 0:
                errors.append("The two bulk VCFs share no SNP position in their first records.")
    if fasta is not None and data:
        df = next(iter(data.values()))
        bad = n = 0
        for r in df.sample(min(len(df), 500), random_state=1).itertuples():
            base = fasta.fetch(r.Chromosome, r.Position, r.Position)
            if base:
                n += 1
                bad += base != r.Ref
        if n and bad / n > 0.05:
            warnings.append(f"{bad} of {n} checked SNPs have a REF base that differs from the reference FASTA "
                            "(wrong reference or coordinates?).")
        if n == 0:
            warnings.append("None of the checked SNPs could be found in the reference FASTA (chromosome names differ?).")
    return {"errors": errors, "warnings": warnings, "files": files}


def write_filtered_vcf(path, snps, sample_high="HIGH", sample_low="LOW", note=""):
    """filtered_informative_SNPs.vcf.gz (block-gzipped): the SNPs kept for the analysis. REF / ALT are written in the
    orientation of the analysis, ALT = the index allele; the AD of the two bulks and their SNP-index are included."""
    from tools.qtl_seq.backend.bgzf import BgzfWriter
    with BgzfWriter(path) as w:
        w.write("##fileformat=VCFv4.2\n##source=TerraTools QTL-seq\n")
        w.write("##note=REF/ALT are in the orientation of the analysis: ALT is the allele the SNP-index counts"
                + (f"; {note}" if note else "") + "\n")
        w.write('##INFO=<ID=DELTA,Number=1,Type=Float,Description="Delta SNP-index (high bulk minus low bulk)">\n')
        w.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Reads supporting REF and ALT">\n')
        w.write('##FORMAT=<ID=SI,Number=1,Type=Float,Description="SNP-index">\n')
        w.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_high + "\t" + sample_low + "\n")
        buf = []
        for r in snps.itertuples(index=False):
            buf.append(f"{r.Chromosome}\t{r.Position}\t.\t{r.Ref}\t{r.Alt}\t.\tPASS\tDELTA={r.Delta_SNP_index:.4f}\tAD:SI\t"
                       f"{r.Ref_High},{r.Alt_High}:{r.SNP_index_High:.4f}\t{r.Ref_Low},{r.Alt_Low}:{r.SNP_index_Low:.4f}\n")
            if len(buf) >= 20000:
                w.write("".join(buf))
                buf = []
        if buf:
            w.write("".join(buf))
