"""Builds the example files of the Genofile Converter from raw allele sheets, as the three formats of a genotyping company:
    Demo_Input_Files/19_Genofile_Converter/input   from the 345 markers x 176 RILs of the Linkage Mapping demo (HapMap, VCF, VCF.GZ, raw CSV)
    docs/sample_data/genofile                      from the small 60 x 80 example of the manual (HapMap, VCF, raw Excel)
Run once:   python docs/make_genofile_demo.py"""
import gzip
import re
import shutil
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / "Demo_Input_Files"
OUT = DEMO / "19_Genofile_Converter" / "input"
SAMPLE = ROOT / "docs" / "sample_data"


def alleles_of(row, samples):
    seen = []
    for v in row[samples]:
        for a in v.split(":"):
            if a and a in "ACGT" and a not in seen:
                seen.append(a)
    if len(seen) == 1:
        seen.append({"A": "G", "C": "T", "G": "A", "T": "C"}[seen[0]])
    return seen[0], seen[1]


def export(df, folder, stem, chrom_label, pos_of, gz=False):
    """df: Marker, Chromosome, P1, P2, then the offspring, calls like A:G -> HapMap and VCF files."""
    folder.mkdir(parents=True, exist_ok=True)
    samples = list(df.columns[2:])
    hap = ["\t".join(["rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center", "protLSID", "assayLSID", "panelLSID", "QCcode", *samples])]
    vcf = ["##fileformat=VCFv4.2", '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
           "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples])]
    for i, row in df.iterrows():
        ref, alt = alleles_of(row, samples)
        chrom, pos = chrom_label(row["Chromosome"]), pos_of(row["Marker"], i)
        h, g = [], []
        for v in row[samples]:
            p = v.split(":")
            if len(p) != 2 or "-" in p or "N" in p:
                h.append("NN"); g.append("./.")
                continue
            h.append("".join(p))
            g.append("/".join(str(0 if a == ref else 1) for a in p))
        hap.append("\t".join([row["Marker"], f"{ref}/{alt}", chrom, str(pos), "+", "NA", "NA", "NA", "NA", "NA", "NA", *h]))
        vcf.append("\t".join([chrom, str(pos), row["Marker"], ref, alt, ".", ".", ".", "GT", *g]))
    (folder / f"{stem}_hapmap.hmp.txt").write_text("\n".join(hap) + "\n", encoding="utf-8")
    (folder / f"{stem}.vcf").write_text("\n".join(vcf) + "\n", encoding="utf-8")
    if gz:
        with gzip.open(folder / f"{stem}.vcf.gz", "wt", encoding="utf-8") as fh:
            fh.write("\n".join(vcf) + "\n")


# 1. the demo: the raw sheet of the Linkage Mapping demo (it now starts here) in the three formats
src = DEMO / "03_Linkage_Mapping" / "input" / "1_raw_allele_calls.csv"
if src.exists():
    OUT.mkdir(parents=True, exist_ok=True)
    shutil.move(str(src), str(OUT / "raw_allele_calls.csv"))
    key = DEMO / "03_Linkage_Mapping" / "input" / "answer_key_faulty_markers.csv"
    if key.exists():
        shutil.move(str(key), str(OUT / key.name))
df = pd.read_csv(OUT / "raw_allele_calls.csv", dtype=str, keep_default_na=False)


def pos_demo(marker, i):
    m = re.search(r"_(\d+)$", marker)
    return int(m.group(1)) if m else 100000 + 1000 * i


export(df, OUT, "genotypes", lambda c: "chr%02d" % int(c), pos_demo, gz=True)

# 2. the manual's small example: the raw Excel example of the Linkage chapter, moved here
raw_x = SAMPLE / "linkage" / "linkage_raw_alleles_example.xlsx"
dst = SAMPLE / "genofile"
dst.mkdir(parents=True, exist_ok=True)
if raw_x.exists():
    shutil.move(str(raw_x), str(dst / "genofile_raw_alleles_example.xlsx"))
small = pd.read_excel(dst / "genofile_raw_alleles_example.xlsx", dtype=str, keep_default_na=False)
export(small, dst, "genofile", lambda c: str(c), lambda m, i: 100000 + 25000 * (int(re.findall(r"\d+", m)[-1])), gz=False)
for old, new in (("genofile_hapmap.hmp.txt", "genofile_hapmap_example.hmp.txt"), ("genofile.vcf", "genofile_vcf_example.vcf")):
    (dst / old).replace(dst / new)
print(sorted(p.name for p in OUT.iterdir()), sorted(p.name for p in dst.iterdir()))
