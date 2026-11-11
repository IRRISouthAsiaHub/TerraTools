"""
gwas_inputs.py — turn "phenotypes + a genotype source" into the GWAS-ready file the engine reads.

GWAS-ready file: first column = sample ID, second column = phenotype, then one column per SNP
named ``Chr<chromosome>_<position>`` with dosage 0 / 1 / 2 (empty = missing).  The chromosome and
position are read back from the column names, so the Manhattan and circos plots use real positions.

Sources
-------
* ``3k``   the bundled rice 3K panel of HapGUI (3,024 accessions x 52,600 SNPs, IRGSP-1.0, PLINK files).
           The user supplies only a phenotype file whose first column holds the accession IDs
           (IRIS IDs such as B001, or the accession names); chromosome and position come from the .bim file.
* ``vcf``  any VCF / VCF.gz (SNPs, GT field) plus a phenotype file whose IDs match the VCF samples.
* ``file`` an already prepared GWAS-ready file (handled by the engine itself).
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd

from tools.gwas.backend.genotype_io import InputError, parse_gt, reject_pooled

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
RICE3K_PLINK = DATA_DIR / "rice3k"                    # rice3k.bed / .bim / .fam (the HapGUI panel)
RICE3K_ACC = DATA_DIR / "rice3k_accessions.csv"

MIN_SAMPLES = 20


# ───────────────────────── phenotype file ─────────────────────────

def read_table(path) -> pd.DataFrame:
    path = str(path)
    low = path.lower()
    try:
        if low.endswith((".xlsx", ".xlsm", ".xls")):
            return pd.read_excel(path)
        return pd.read_csv(path, sep=None, engine="python")
    except Exception as e:
        raise InputError(f"Could not read the phenotype file: {e}") from None


def phenotype_columns(path) -> list[str]:
    """Names of the columns after the ID column that hold numbers (the traits)."""
    df = read_table(path)
    out = []
    for c in df.columns[1:]:
        if pd.to_numeric(df[c], errors="coerce").notna().sum() >= 2:
            out.append(str(c))
    return out


def _norm(s) -> str:
    return re.sub(r"[^0-9a-z]", "", str(s).lower())


def _load_phenotype(path, trait):
    df = read_table(path)
    if df.shape[1] < 2:
        raise InputError("The phenotype file needs an ID column followed by at least one trait column.")
    id_col = df.columns[0]
    cols = [str(c) for c in df.columns]
    if trait not in cols:
        raise InputError(f"Trait '{trait}' is not a column of the phenotype file. Columns: {', '.join(cols[1:])}")
    t = df.columns[cols.index(trait)]
    y = pd.to_numeric(df[t], errors="coerce")
    ids = df[id_col].astype(str).str.strip()
    ok = y.notna() & ids.ne("") & ids.ne("nan")
    sub = pd.DataFrame({"id": ids[ok].values, "y": y[ok].values})
    dups = int(sub.id.duplicated().sum())
    if dups:
        sub = sub.groupby("id", sort=False, as_index=False)["y"].mean()   # repeated plots / reps -> mean
    return sub, id_col, dups, int((~ok).sum())


def match_samples(pheno_ids, geno_ids, aliases=None):
    """Map each phenotype ID to a genotype ID: exact first, then ignoring case / punctuation
    (IRIS 313-9917 == iris_313_9917), then through `aliases` (a list of alternative names per
    genotype ID, e.g. the accession name). A name shared by several samples is never used.
    Returns ([(pheno_idx, geno_idx)], unmatched pheno IDs)."""
    geno_ids = list(geno_ids)
    exact = {g: i for i, g in enumerate(geno_ids)}
    loose = {}

    def add(name, i):
        k = _norm(name)
        if k:
            loose.setdefault(k, set()).add(i)

    for i, g in enumerate(geno_ids):
        add(g, i)
    for i, names in enumerate(aliases or []):
        for a in names:
            add(a, i)
    pairs, missing, used = [], [], set()
    for i, p in enumerate(pheno_ids):
        j = exact.get(p)
        if j is None:
            hit = loose.get(_norm(p))
            j = next(iter(hit)) if hit and len(hit) == 1 else None
        if j is None or j in used:
            missing.append(p)
        else:
            used.add(j)
            pairs.append((i, j))
    return pairs, missing


def _need_matches(pairs, pheno_ids, missing, geno_ids, what):
    if len(pairs) >= MIN_SAMPLES:
        return
    ex_p = ", ".join(map(str, list(pheno_ids)[:3]))
    ex_g = ", ".join(map(str, list(geno_ids)[:3]))
    raise InputError(
        f"Only {len(pairs)} of {len(pheno_ids)} phenotype IDs match a {what} sample (at least {MIN_SAMPLES} are needed).\n"
        f"Your first IDs: {ex_p}\n{what} IDs look like: {ex_g}\n"
        "The first column of the phenotype file must hold the sample IDs exactly as they are written in the genotype data.")


# ───────────────────────── writing the GWAS-ready file ─────────────────────────

def _unique_names(chrom, pos):
    names, seen = [], {}
    for c, p in zip(chrom, pos):
        n = f"Chr{c}_{p}"
        k = seen.get(n, 0)
        seen[n] = k + 1
        names.append(n if k == 0 else f"{n}_{k + 1}")
    return names


def write_gwas_ready(path, sample_ids, y, dosage, chrom, pos):
    """dosage: (samples x markers) int8 array with -1 for missing."""
    # Written by hand: a DataFrame with tens of thousands of columns takes minutes to build and to_csv.
    names = _unique_names(chrom, pos)
    dosage = np.asarray(dosage)
    text = np.array(["NA", "0", "1", "2"])[np.clip(dosage, -1, 2).astype(np.int8) + 1]          # (samples x markers)
    y = np.round(np.asarray(y, dtype=float), 6)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(",".join(["ID", "Phenotype"] + names) + "\n")
        for sid, yy, row in zip(sample_ids, y, text):
            sid = str(sid)
            if "," in sid or '"' in sid:
                sid = '"' + sid.replace('"', '""') + '"'
            fh.write(f"{sid},{float(yy)!r},{','.join(row)}\n")
    return str(path)


# ───────────────────────── rice 3K ─────────────────────────
#
# The panel is the PLINK binary set used by HapGUI: rice3k.bed / .bim / .fam (3,024 accessions x 52,600
# SNPs, IRGSP-1.0) with rice3k_accessions.csv (IRIS_ID, accession name, subpopulation, country, IRGC no.).
# Its IDs are "B001" ... and the "IRIS_313-xxxx" identifiers of the earlier TerraTools panel (which it contains).

def _have_plink_panel() -> bool:
    return all(Path(str(RICE3K_PLINK) + ext).exists() for ext in (".bed", ".bim", ".fam"))


def read_fam_ids(prefix) -> list[str]:
    with open(str(prefix) + ".fam", "r", encoding="utf-8", errors="replace") as fh:
        return [ln.split()[1] for ln in fh if ln.strip()]


def read_bim(prefix) -> pd.DataFrame:
    return pd.read_csv(str(prefix) + ".bim", sep=r"\s+", header=None, dtype={0: str},
                       names=["chrom", "id", "cm", "pos", "a1", "a2"])


def read_plink(prefix, sample_idx=None) -> np.ndarray:
    """Dosage of allele A1 (the first allele of the .bim file) as an int8 array (samples x markers),
    -1 for a missing call. `sample_idx` selects (and orders) samples by their position in the .fam file."""
    n = len(read_fam_ids(prefix))
    m = len(read_bim(prefix))
    raw = np.fromfile(str(prefix) + ".bed", dtype=np.uint8)
    if raw.size < 3 or raw[0] != 0x6C or raw[1] != 0x1B:
        raise InputError("The rice 3K panel file (.bed) is not a PLINK binary file.")
    if raw[2] != 1:
        raise InputError("The rice 3K panel is not stored SNP-major.")
    nb = (n + 3) // 4
    if raw.size - 3 != m * nb:
        raise InputError("The rice 3K panel files (.bed, .bim, .fam) do not agree in size.")
    packed = raw[3:].reshape(m, nb)
    codes = ((packed[:, :, None] >> np.array([0, 2, 4, 6], dtype=np.uint8)) & 3).reshape(m, nb * 4)[:, :n]
    if sample_idx is not None:
        codes = codes[:, np.asarray(sample_idx, dtype=int)]
    lut = np.array([2, -1, 1, 0], dtype=np.int8)          # 00 = 2 x A1, 01 = missing, 10 = A1/A2, 11 = 2 x A2
    return lut[codes].T.copy()


def _need_panel():
    if not _have_plink_panel():
        raise InputError("The bundled rice 3K panel is missing (tools/gwas/data/rice3k.bed / .bim / .fam).")


def _accession_table(fam_ids) -> pd.DataFrame:
    acc = pd.read_csv(RICE3K_ACC, dtype=str)
    return acc.drop_duplicates("IRIS_ID").set_index("IRIS_ID").reindex(fam_ids)


def rice3k_info() -> dict:
    _need_panel()
    ids = read_fam_ids(RICE3K_PLINK)
    with open(str(RICE3K_PLINK) + ".bim", "rb") as fh:
        n_markers = sum(1 for _ in fh)
    return {"n_samples": len(ids), "n_markers": n_markers, "examples": ids[:5], "id_column": "IRIS_ID", "ids": ids}


def rice3k_template(path) -> str:
    """A phenotype file to fill in: every accession of the panel with its name, subpopulation and country."""
    _need_panel()
    ids = read_fam_ids(RICE3K_PLINK)
    acc = _accession_table(ids)
    pd.DataFrame({"IRIS_ID": ids, "Accession": acc["Accession"].values, "Subpopulation": acc["Subpopulation"].values,
                  "Country": acc["Country"].values, "Trait1": ""}).to_csv(path, index=False)
    return str(path)


def build_from_3k(pheno_path, trait, out_dir, maf=0.05, max_missing=0.2, log=print):
    _need_panel()
    sub, _id_col, dups, n_bad = _load_phenotype(pheno_path, trait)
    ids = sub.id.tolist()
    gids = read_fam_ids(RICE3K_PLINK)
    acc = _accession_table(gids)
    # an accession can be named by its IRIS ID, its name, or (when unique) its IRGC number
    aliases = [[a for a in (n, g) if isinstance(a, str) and a] for n, g in zip(acc["Accession"], acc["IRGC_no"])]
    pairs, missing = match_samples(ids, gids, aliases)
    _need_matches(pairs, ids, missing, gids, "3K panel")

    pi = [p for p, _ in pairs]
    gi = [g for _, g in pairs]
    dos = read_plink(RICE3K_PLINK, gi)                                       # samples x markers, -1 = missing
    bim = read_bim(RICE3K_PLINK)
    called = dos >= 0
    n_called = called.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        freq = np.where(n_called > 0, np.where(called, dos, 0).sum(0) / (2.0 * np.maximum(n_called, 1)), 0.0)
    minor = np.minimum(freq, 1 - freq)
    keep = (n_called >= (1 - max_missing) * dos.shape[0]) & (minor >= maf) & (minor > 0)     # monomorphic SNPs carry nothing
    chrom = pd.to_numeric(bim.chrom, errors="coerce")
    keep &= chrom.notna().to_numpy()
    if not keep.any():
        raise InputError(f"No SNP of the 3K panel passes the filters for these {len(pairs)} accessions "
                         f"(missing > {int(max_missing * 100)}% or MAF < {maf}).")
    log(f"[GWAS] 3K panel: {len(pairs)} accessions matched, {int(keep.sum()):,} of {len(keep):,} SNPs kept")
    sample_ids = sub.id.to_numpy()[pi]
    path = write_gwas_ready(Path(out_dir) / "gwas_ready_input.csv", sample_ids, sub.y.to_numpy()[pi],
                            dos[:, keep], chrom.to_numpy()[keep].astype(int), bim.pos.to_numpy()[keep])
    a1, a2 = bim.a1.to_numpy()[keep], bim.a2.to_numpy()[keep]
    info = {"source": "Rice 3K panel", "trait": trait, "phenotype_ids": len(sub), "matched": len(pairs),
            "unmatched": len(missing), "unmatched_examples": missing[:8], "markers": int(keep.sum()),
            "duplicate_ids_averaged": dups, "rows_without_value": n_bad,
            "sample_ids": [str(s) for s in sample_ids],
            "subpopulation": [None if not isinstance(v, str) else v for v in acc["Subpopulation"].to_numpy()[gi]],
            "alleles": [f"{b}/{a}" for a, b in zip(a1, a2)]}          # other/counted: the dosage counts copies of A1
    return path, info


# ───────────────────────── VCF ─────────────────────────

def _chrom_number(c: str):
    c = re.sub(r"^(chromosome|chr)[_\-]?", "", str(c).strip(), flags=re.I)
    return int(c) if c.isdigit() else None


def build_from_vcf(vcf_path, pheno_path, trait, out_dir, maf=0.05, max_missing=0.2, log=print):
    """VCF -> GWAS-ready file. Genotype calls are read strictly: diploid biallelic SNPs only, phased (|) and unphased (/) calls
    alike; a call with a missing allele is missing; another ploidy or an allele number above 1 stops the run with the position
    named; samples named like pooled bulks are refused. Markers filtered here (call rate, MAF) are counted in the result."""
    sub, _id_col, dups, n_bad = _load_phenotype(pheno_path, trait)
    opener = gzip.open if str(vcf_path).lower().endswith(".gz") else open
    samples = None
    sel = pairs = missing = None
    rows, chrom, pos, allele_txt = [], [], [], []
    n_seen = n_snp = n_skip_chr = n_low_cr = n_low_maf = n_multi = n_indel = n_no_gt = n_partial = 0
    with opener(vcf_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = [s.strip() for s in f[9:]]
                if not samples:
                    raise InputError("The VCF has no sample columns.")
                pairs, missing = match_samples(sub.id.tolist(), samples)
                _need_matches(pairs, sub.id.tolist(), missing, samples, "VCF")
                reject_pooled([samples[j] for _, j in pairs])
                sel = [9 + j for _, j in pairs]
                continue
            if samples is None or len(f) < 10:
                continue
            n_seen += 1
            ref, alt = f[3], f[4]
            if "," in alt:
                n_multi += 1
                continue
            if len(ref) != 1 or len(alt) != 1 or alt == ".":          # SNPs, biallelic only
                n_indel += 1
                continue
            fmt = f[8].split(":")
            if "GT" not in fmt:
                n_no_gt += 1
                continue
            gt_i = fmt.index("GT")
            n_snp += 1
            c = _chrom_number(f[0])
            if c is None:
                n_skip_chr += 1
                continue
            where = f" at {f[0]}:{f[1]}"
            if gt_i == 0:
                vals = np.fromiter((parse_gt(f[i], where) for i in sel), dtype=float, count=len(sel))
            else:
                vals = np.fromiter((parse_gt(f[i].split(":")[gt_i], where) for i in sel), dtype=float, count=len(sel))
            nan = np.isnan(vals)
            if nan.any():                                              # a call with only one missing allele is counted separately
                for i in np.flatnonzero(nan):
                    gt = f[sel[i]].split(":")[gt_i] if gt_i else f[sel[i]].split(":", 1)[0]
                    if gt not in ("./.", ".|.", "."):
                        n_partial += 1
            g = np.where(nan, -1, vals).astype(np.int8)
            ok = g >= 0
            if ok.mean() < 1 - max_missing:
                n_low_cr += 1
                continue
            p = g[ok].mean() / 2
            if min(p, 1 - p) < maf:
                n_low_maf += 1
                continue
            rows.append(g)
            chrom.append(c)
            pos.append(int(f[1]))
            allele_txt.append(f"{ref}/{alt}")                          # other/counted: the dosage counts copies of ALT
            if n_seen % 200_000 == 0:
                log(f"[GWAS] VCF: {n_seen:,} variants read, {len(rows):,} kept")
    if not rows:
        raise InputError(
            f"No SNP passed the filters ({n_seen:,} variants read: {n_snp:,} biallelic SNPs, {n_skip_chr:,} on "
            f"non-numeric chromosomes, {n_low_cr:,} with > {int(max_missing * 100)}% missing, {n_low_maf:,} with MAF < {maf}). "
            "The VCF needs genotype calls (GT) and chromosomes named 1, 2, … or Chr1, Chr2, …")
    mat = np.stack(rows, axis=1)                                         # samples x markers
    pi = [p for p, _ in pairs]
    order = np.lexsort((np.array(pos), np.array(chrom)))
    path = write_gwas_ready(Path(out_dir) / "gwas_ready_input.csv", sub.id.to_numpy()[pi], sub.y.to_numpy()[pi],
                            mat[:, order], np.array(chrom)[order], np.array(pos)[order])
    info = {"source": "VCF", "trait": trait, "phenotype_ids": len(sub), "matched": len(pairs), "unmatched": len(missing),
            "unmatched_examples": missing[:8], "markers": len(rows), "variants_read": n_seen, "snps": n_snp,
            "skipped_non_numeric_chromosome": n_skip_chr, "dropped_missing": n_low_cr, "dropped_maf": n_low_maf,
            "skipped_multiallelic": n_multi, "skipped_indel": n_indel, "skipped_no_gt": n_no_gt,
            "partially_missing_calls": n_partial, "prefilter_maf": maf, "prefilter_missing": max_missing,
            "duplicate_ids_averaged": dups, "rows_without_value": n_bad,
            "sample_ids": [str(s) for s in sub.id.to_numpy()[pi]], "subpopulation": None,
            "alleles": [allele_txt[i] for i in order]}
    return path, info


# ───────────────────────── entry point ─────────────────────────

def prepare_input(cfg, out_dir, log=print):
    """-> (path of the GWAS-ready file, info dict or None)."""
    mode = str(getattr(cfg, "input_mode", "file") or "file").lower()
    if mode == "file":
        return cfg.input_file, None
    if not cfg.phenotype_file:
        raise InputError("Please select the phenotype file.")
    if not cfg.trait:
        raise InputError("Please choose the trait (phenotype column) to analyse.")
    if mode == "3k":
        # No MAF / call-rate filter here: the engine applies them per model (the multi-locus models, like the
        # R package and HapGUI's 3K route, see every SNP).
        return build_from_3k(cfg.phenotype_file, cfg.trait, out_dir, maf=0.0, max_missing=1.0, log=log)
    if mode == "vcf":
        if not cfg.vcf_file:
            raise InputError("Please select the VCF file.")
        return build_from_vcf(cfg.vcf_file, cfg.phenotype_file, cfg.trait, out_dir, maf=float(cfg.maf),
                              max_missing=float(getattr(cfg, "max_snp_missing", 0.2)), log=log)
    raise InputError(f"Unknown input type: {mode}")


def describe(info: dict) -> str:
    lines = [f"Input: {info['source']} — trait '{info['trait']}'",
             f"Phenotype IDs matched to genotypes: {info['matched']} of {info['phenotype_ids']}"]
    if info["unmatched"]:
        lines.append(f"  not found ({info['unmatched']}): " + ", ".join(map(str, info["unmatched_examples"]))
                     + (" …" if info["unmatched"] > len(info["unmatched_examples"]) else ""))
    if info.get("duplicate_ids_averaged"):
        lines.append(f"  {info['duplicate_ids_averaged']} repeated IDs were averaged")
    if info["source"] == "VCF":
        lines.append(f"VCF: {info['variants_read']:,} variants read, {info['snps']:,} biallelic SNPs with a GT field, {info['markers']:,} kept "
                     f"({info['dropped_missing']:,} too many missing, {info['dropped_maf']:,} rare, "
                     f"{info['skipped_non_numeric_chromosome']:,} on non-numeric chromosomes)")
        lines.append(f"  not used: {info['skipped_multiallelic']:,} multi-allelic and {info['skipped_indel']:,} indel records, "
                     f"{info['skipped_no_gt']:,} without GT; {info['partially_missing_calls']:,} calls with one missing allele were "
                     "treated as missing")
    else:
        lines.append(f"SNPs in the 3K panel with data: {info['markers']:,}")
    return "\n".join(lines)
