"""
Generates the small example input files that accompany the TerraTools user manual
(docs/sample_data/<tool>/...). The data is simulated — it shows the exact file layout each
tool expects and runs quickly; it is not real biology.

    python docs/make_sample_data.py
"""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent / "sample_data"
rng = np.random.default_rng(2026)


def out(tool, name):
    d = OUT / tool
    d.mkdir(parents=True, exist_ok=True)
    return d / name


def ids(prefix, n, width=3):
    return [f"{prefix}{i + 1:0{width}d}" for i in range(n)]


def sim_f2(n, positions_by_chr, het=True):
    """F2 / RIL genotypes along chromosomes (Haldane); 0 / 1 / 2 (RIL: 0 / 2)."""
    cols = []
    for pos in positions_by_chr:
        m = len(pos)
        G = np.zeros((n, m))
        for i in range(n):
            haps = []
            for _ in range(2 if het else 1):
                h = np.zeros(m)
                h[0] = rng.integers(0, 2)
                for j in range(1, m):
                    r = 0.5 * (1 - np.exp(-2 * (pos[j] - pos[j - 1]) / 100))
                    h[j] = 1 - h[j - 1] if rng.random() < r else h[j - 1]
                haps.append(h)
            G[i] = sum(haps) if het else 2 * haps[0]
        cols.append(G)
    return np.hstack(cols)


# ───────────────────────── 1. Mendelian Segregation ─────────────────────────
n, per, nchr = 60, 30, 5
chrom = np.repeat(np.arange(1, nchr + 1), per)
G = rng.choice([0, 2], size=(per * nchr, n)).astype(object)
for r in (7, 8, 61, 95, 130):                       # a few distorted markers
    G[r] = rng.choice([0, 2], size=n, p=[0.82, 0.18])
G[rng.random(G.shape) < 0.03] = "NA"
df = pd.DataFrame(G, columns=ids("RIL", n))
df.insert(0, "Chromosome", chrom)
df.insert(0, "Marker", [f"M{c}_{i + 1:02d}" for c in range(1, nchr + 1) for i in range(per)])
df.to_csv(out("mendelian", "mendelian_RIL_example.csv"), index=False)

# ───────────────────────── 2. Linkage Mapping ─────────────────────────
n, per = 80, 20
pos = [np.arange(per) * 4.0 for _ in range(3)]
G = sim_f2(n, pos, het=True).T.astype(int)
G[rng.random(G.shape) < 0.04] = -1
marker = [f"SNP{c}_{i + 1:02d}" for c in range(1, 4) for i in range(per)]
chrom = np.repeat([1, 2, 3], per)
df = pd.DataFrame(G, columns=ids("F2_", n))
df.insert(0, "Chromosome", chrom)
df.insert(0, "Marker", marker)
df.to_excel(out("linkage", "linkage_F2_coded_example.xlsx"), index=False)

alleles = [("A", "G"), ("C", "T"), ("A", "C"), ("G", "T")]
raw = pd.DataFrame({"Marker": marker, "Chromosome": chrom})
a1, a2 = zip(*[alleles[i % 4] for i in range(len(marker))])
raw["P1"] = [f"{a}:{a}" for a in a1]
raw["P2"] = [f"{b}:{b}" for b in a2]
sym = {0: lambda a, b: f"{a}:{a}", 1: lambda a, b: f"{a}:{b}", 2: lambda a, b: f"{b}:{b}", -1: lambda a, b: "-"}
for j, col in enumerate(ids("F2_", n)):
    raw[col] = [sym[int(G[i, j])](a1[i], a2[i]) for i in range(len(marker))]
raw.to_excel(out("linkage", "linkage_raw_alleles_example.xlsx"), index=False)

# ───────────────────────── 3. QTL Mapping ─────────────────────────
n = 120
pos = [np.arange(0, 101, 5.0) for _ in range(5)]
G = sim_f2(n, pos, het=True)
mk = [f"c{c}m{j + 1:02d}" for c in range(1, 6) for j in range(len(pos[0]))]
chr_row = np.repeat(np.arange(1, 6), len(pos[0]))
cm_row = np.concatenate(pos)
sid = ids("L", n)
y = (1.0 * (G[:, 21 + 8] - 1) + 0.8 * (G[:, 63 + 4] - 1) - 0.6 * (G[:, 63 + 15] - 1) + rng.normal(0, 1, n))
with open(out("qtl", "qtl_genotypes_example.csv"), "w") as f:
    f.write(",".join(["Marker"] + mk) + "\n")
    f.write(",".join(["Chromosome"] + [str(c) for c in chr_row]) + "\n")
    f.write(",".join(["Position"] + [f"{p:g}" for p in cm_row]) + "\n")
    for i in range(n):
        f.write(",".join([sid[i]] + [str(int(v)) for v in G[i]]) + "\n")
pd.DataFrame({"ID": sid, "Yield": y.round(3), "Height": (100 + 8 * rng.normal(size=n)).round(2)}).to_csv(
    out("qtl", "qtl_phenotypes_example.csv"), index=False)

# ───────────────────────── 4. Genes Extraction ─────────────────────────
pd.DataFrame({"QTL_ID": ["QTL_3.1", "QTL_3.2", "QTL_5.1"],
              "QTL_Start": [4_000_000, 10_000_000, 22_500_000],
              "QTL_End": [4_600_000, 11_500_000, 23_300_000]}).to_csv(
    out("gene_extractor", "qtl_intervals_example.csv"), index=False)
grng = np.random.default_rng(4)      # own generator, so the other example files stay as they were
_gene_pos = {3: [], 5: []}
with open(out("gene_extractor", "annotation_example.gff3"), "w") as f:
    f.write("##gff-version 3\n")
    for c in (3, 5):
        for k in range(60):
            s = int(grng.integers(500_000, 30_000_000))
            e = s + int(grng.integers(800, 6000))
            _gene_pos[c].append((s, e))
            f.write(f"{c}\texample\tgene\t{s}\t{e}\t.\t{grng.choice(['+', '-'])}\t.\t"
                    f"ID=gene:Os{c:02d}g{k + 1:07d};Name=OsEx{c}_{k + 1};biotype=protein_coding;"
                    f"description=hypothetical protein {k + 1}\n")

# the QTL table of QTL Mapping (positions in cM) whose markers are called M001 ... and so carry no position in their name,
# and the marker position file (Marker, Chromosome, bp) that gives their base-pair positions
_qtl = [("QTL_3.1", 3, "M001", 18.0, "M002", 22.5, "M003", 27.0, 19.5, 25.5, 3_700_000, 4_250_000, 4_800_000),
        ("QTL_5.1", 5, "M004", 40.0, "M005", 44.0, "M006", 49.0, 41.5, 47.0, 22_000_000, 22_800_000, 23_600_000)]
_lod = (9.4, 6.1)
pd.DataFrame([{"QTL_ID": q[0], "Chromosome": q[1], "Left_Marker": q[2], "Left_Position": q[3], "Peak_Marker": q[4],
               "Peak_Position": q[5], "Right_Marker": q[6], "Right_Position": q[7], "CI_Start_cM": q[8], "CI_End_cM": q[9],
               "LOD": lod} for q, lod in zip(_qtl, _lod)]).to_csv(out("gene_extractor", "qtl_table_example.csv"), index=False)
pd.DataFrame({"Marker": [q[k] for q in _qtl for k in (2, 4, 6)],
              "Chromosome": [f"chr{q[1]:02d}" for q in _qtl for _ in range(3)],
              "bp": [q[k] for q in _qtl for k in (10, 11, 12)]}).to_csv(out("gene_extractor", "marker_positions_example.csv"), index=False)

# ───────────────────────── 4b. Gene Prioritization ─────────────────────────
# Real rice (IRGSP-1.0 / RAP-DB) gene IDs, so the databases can answer: four well-known drought genes (OsbZIP23, SNAC1-like
# NAC22, DREB2A, OsNCED1) and two genes with no drought annotation. 
out("gene_prioritization", "candidate_genes_rice.txt").write_text(
    "Os02g0766700\nOs03g0133000\nOs01g0165000\nOs03g0645900\nOs09g0439200\nOs05g0426100\n", encoding="utf-8")

# ───────────────────────── 5. PYBSA G-statistics ─────────────────────────
_pb_positions = {}


def gatk_table(path, samples, n_sv=1500):
    rows = []
    for c in range(1, 4):
        drawn = np.sort(rng.integers(10_000, 20_000_000, n_sv // 3))
        pos = _pb_positions.setdefault(c, drawn)        # parents and bulks describe the same variants
        for p in pos:
            rows.append([str(c), int(p), "A", "G"])
    df = pd.DataFrame(rows, columns=["CHROM", "POS", "REF", "ALT"])
    df["QUAL"] = 50.0
    for s, kind in samples:
        gt, ad, gq = [], [], []
        for c, p in zip(df.CHROM, df.POS):
            depth = int(rng.poisson(30)) + 5
            if kind == "ref":
                rng.binomial(depth, 0.01)    # (draw kept so the other example files stay unchanged)
                alt = 0                      # parents are fully homozygous: no reads of the other allele
            elif kind == "alt":
                rng.binomial(depth, 0.99)
                alt = depth
            else:
                shift = 0.3 * np.exp(-abs(p - 8_000_000) / 2_000_000) if c == "2" else 0
                alt = rng.binomial(depth, 0.5 + (shift if kind == "bulkA" else -shift))
            ref = depth - alt
            gt.append("A/A" if alt / depth < 0.2 else "G/G" if alt / depth > 0.8 else "A/G")
            ad.append(f"{ref},{alt}")
            gq.append(int(rng.integers(30, 99)))
        df[f"{s}.GT"], df[f"{s}.AD"], df[f"{s}.GQ"] = gt, ad, gq
    df.to_csv(path, sep="\t", index=False)


gatk_table(out("pybsaseq", "parents_example.tsv"), [("P1", "ref"), ("P2", "alt")])
gatk_table(out("pybsaseq", "bulks_example.tsv"), [("HighBulk", "bulkA"), ("LowBulk", "bulkB")])

# ───────────────────────── 6. QTL-Seq (VCF) ─────────────────────────
with gzip.open(out("qtl_seq", "qtlseq_example.vcf.gz"), "wt") as f:
    f.write("##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCultivar\tHighBulk\tLowBulk\n")
    for c in range(1, 4):
        for p in np.sort(rng.integers(10_000, 6_000_000, 2000)):
            shift = 0.38 * np.exp(-abs(p - 3_000_000) / 700_000) if c == 2 else 0
            def reads(freq):
                d = int(rng.poisson(30)) + 4
                alt = int(rng.binomial(d, min(max(freq, 0.0), 1.0)))
                return f"./.:{d - alt},{alt}"
            f.write(f"Chr{c}\t{int(p)}\t.\tA\tG\t50\tPASS\t.\tGT:AD\t{reads(0.003)}\t"
                    f"{reads(rng.binomial(60, 0.5 + shift) / 60)}\t{reads(rng.binomial(60, 0.5 - shift) / 60)}\n")

# ───────────────────────── 7. Diversity Analysis ─────────────────────────
varieties = ["Aus", "Indica", "Japonica", "Basmati"]
rows = []
for vi, v in enumerate(varieties):
    base = rng.integers(0, 4, 40)
    for k in range(10):
        blk = np.where(rng.random(40) < 0.25, rng.integers(0, 4, 40), base)
        blk = blk.astype(object)
        blk[rng.random(40) < 0.05] = -1
        rows.append([v, f"{v[:3].upper()}_{k + 1:02d}"] + list(blk))
pd.DataFrame(rows, columns=["Variety", "Sample"] + [f"Block{i + 1:02d}" for i in range(40)]).to_csv(
    out("diversity", "haplotype_blocks_example.csv"), index=False)

# ───────────────────────── 8. GWAS ─────────────────────────
n = 150
nm_chr, nchr = 40, 5
pos_by = [np.sort(rng.integers(100_000, 25_000_000, nm_chr)) for _ in range(nchr)]
G = rng.choice([0, 1, 2], size=(n, nm_chr * nchr), p=[0.45, 0.1, 0.45]).astype(float)
y = 1.1 * G[:, 25] + 0.7 * G[:, 100] + rng.normal(0, 1, n)
G[rng.random(G.shape) < 0.02] = np.nan
names = [f"Chr{c + 1}_{p}" for c in range(nchr) for p in pos_by[c]]
df = pd.DataFrame(G, columns=names)
df.insert(0, "Phenotype", y.round(3))
df.insert(0, "ID", ids("G", n))
df.to_csv(out("gwas", "gwas_ready_example.csv"), index=False, na_rep="NA")
# 3K panel phenotype example (real accession IDs of the bundled HapGUI panel, B001 ...)
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tools.gwas.backend import gwas_inputs as _g3k
_fam = _g3k.read_fam_ids(_g3k.RICE3K_PLINK)
_pick = np.sort(np.random.default_rng(3).choice(len(_fam), 250, replace=False))   # own generator: the shared stream is untouched
_dos = _g3k.read_plink(_g3k.RICE3K_PLINK, _pick).astype(float)
_dos[_dos < 0] = 1.0
_frq = _dos.mean(0) / 2
_cand = np.where((_frq > 0.25) & (_frq < 0.75))[0]
m1, m2 = _cand[1500], _cand[900]                    # two simulated QTL; the rest of the trait is noise
pheno = (5 + 0.8 * _dos[:, m1] + 0.5 * _dos[:, m2] + rng.normal(0, 0.6, 250)).round(2)
pd.DataFrame({"SampleID": [_fam[i] for i in _pick], "Yield": pheno, "Height": rng.normal(110, 12, 250).round(1)}).to_csv(
    out("gwas", "rice3k_phenotype_example.csv"), index=False)
# VCF + matching phenotype
with gzip.open(out("gwas", "gwas_example.vcf.gz"), "wt") as f:
    f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids("G", n)) + "\n")
    gt = {0: "0/0", 1: "0/1", 2: "1/1"}
    for j, nme in enumerate(names):
        c, p = nme[3:].split("_")
        f.write(f"{c}\t{p}\t.\tA\tG\t50\tPASS\t.\tGT\t" + "\t".join(gt.get(v, "./.") if v == v else "./." for v in G[:, j]) + "\n")
pd.DataFrame({"ID": ids("G", n), "Yield": y.round(3), "Height": rng.normal(100, 9, n).round(1)}).to_csv(
    out("gwas", "gwas_vcf_phenotype_example.csv"), index=False)

# ───────────────────────── 9. Clustering & Genomic Selection ─────────────────────────
n_per, m = 30, 100
rows = []
for g in range(3):
    f = rng.random(m) * 0.8 + 0.1
    for k in range(n_per):
        rows.append(rng.binomial(2, np.clip(f + rng.normal(0, 0.05, m), 0.02, 0.98)))
G = np.array(rows).astype(int)
G[rng.random(G.shape) < 0.02] = -1
df = pd.DataFrame(G, columns=[f"SNP{i + 1:03d}" for i in range(m)])
df.insert(0, "SampleID", ids("Acc", 3 * n_per))
df.to_csv(out("clustering", "clustering_genotypes_example.csv"), index=False)

# ───────────────────────── 10. Additive Dominance ─────────────────────────
n, m = 40, 200
G = rng.choice([0, 1, 2], size=(n, m), p=[0.4, 0.2, 0.4]).astype(object)
G[rng.random(G.shape) < 0.02] = ""
df = pd.DataFrame(G, columns=[f"SNP{i + 1:03d}" for i in range(m)])
df.insert(0, "SampleID", ids("Line", n))
df.to_csv(out("additive_dominance", "snp_dosage_example.csv"), index=False)

# ───────────────────────── 11. BreedSight Genomic Prediction ─────────────────────────
rng.choice([0, 1, 2], size=(150, 150), p=[0.4, 0.2, 0.4])        # (draws kept so later example files stay unchanged)
rng.normal(0, 0.4, 150); rng.random(150); rng.normal(0, 0.6, 150)
r3 = np.random.default_rng(88)
ntr, nte, m = 320, 40, 150
Gt = r3.choice([0, 1, 2], size=(ntr + nte, m), p=[0.4, 0.2, 0.4])
eff = r3.normal(0, 0.5, m) * (r3.random(m) < 0.12)
yy = (Gt @ eff + r3.normal(0, 0.8, ntr + nte)).round(3) + 5
snp = [f"SNP{i + 1:03d}" for i in range(m)]
sid = ids("Line", ntr + nte)
pd.DataFrame({"sample_id": sid[:ntr], "phenotypes": yy[:ntr]}).to_csv(out("genomic_prediction", "train_phenotypes_example.csv"), index=False)
pd.DataFrame(Gt[:ntr], columns=snp).assign(sample_id=sid[:ntr])[["sample_id"] + snp].to_csv(out("genomic_prediction", "train_genotypes_example.csv"), index=False)
pd.DataFrame({"sample_id": sid[ntr:], "phenotypes": yy[ntr:]}).to_csv(out("genomic_prediction", "test_phenotypes_example.csv"), index=False)
pd.DataFrame(Gt[ntr:], columns=snp).assign(sample_id=sid[ntr:])[["sample_id"] + snp].to_csv(out("genomic_prediction", "test_genotypes_example.csv"), index=False)

# ───────────────────────── 12. Optimal Cross Selection ─────────────────────────
n, m = 40, 300
rng.choice([0, 1, 2], size=(n, m), p=[0.45, 0.1, 0.45])          # (draw kept so later example files stay unchanged)
crit = rng.normal(10, 2, n).round(3)
# four families built from overlapping founders -> a realistic spread of kinship between lines
r2 = np.random.default_rng(77)
founders = r2.choice([0, 2], size=(6, m))
rows = []
for fam in range(4):
    trio = [founders[(fam + k) % 6] for k in range(3)]
    for _ in range(10):
        pick = r2.integers(0, 3, m)
        rows.append([trio[pick[j]][j] for j in range(m)])
Gc = np.array(rows, dtype=float)
p = Gc.mean(0) / 2
geno = ids("Line", n)

# marker data for the usefulness mode: the same 40 inbred lines as a 0 / 2 dosage matrix, additive marker effects and a genetic map
r5 = np.random.default_rng(78)
mk = [f"SNP{i + 1:03d}" for i in range(m)]
pd.DataFrame(Gc.astype(int), columns=mk).assign(Genotype=geno)[["Genotype"] + mk].to_csv(out("cross_selection", "dosage_matrix_example.csv"), index=False)
eff_m = (r5.normal(0, 0.25, m) * (r5.random(m) < 0.25)).round(4)
pd.DataFrame({"Marker": mk, "Effect": eff_m}).to_csv(out("cross_selection", "marker_effects_example.csv"), index=False)
nchr = 5
per = m // nchr
gm_rows = []
for c in range(nchr):
    cm = np.cumsum(r5.uniform(0.4, 2.4, per)).round(2)
    gm_rows += [(mk[c * per + k], f"Chr{c + 1}", cm[k]) for k in range(per)]
pd.DataFrame(gm_rows, columns=["Marker", "Chromosome", "cM"]).to_csv(out("cross_selection", "genetic_map_example.csv"), index=False)
# the same markers as a VCF (physical positions, 4 cM per Mb), for the case without a genetic map
with open(out("cross_selection", "markers_example.vcf"), "w", newline="\n") as fh:
    fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for name, chrom, cm in gm_rows:
        fh.write(f"{chrom}\t{int(round(cm / 4.0 * 1e6))}\t{name}\tA\tG\t.\t.\t.\n")

# ───────────────────────── 13. Heterosis ─────────────────────────
ntr, npar, m = 100, 12, 120
snp = [f"SNP{i + 1:03d}" for i in range(m)]
Gt = rng.choice([0, 1, 2], size=(ntr, m)).astype(float)
Gp = rng.choice([0, 2], size=(npar, m)).astype(float)
Gte = rng.choice([0, 1, 2], size=(10, m)).astype(float)
eff = rng.normal(0, 0.3, m) * (rng.random(m) < 0.2)
yy = (Gt @ eff + 6 + rng.normal(0, 0.5, ntr)).round(3)      # a trait on a typical scale (e.g. grain yield, t/ha)
pd.DataFrame({"Sample": ids("Train", ntr), "Yield": yy}).to_csv(out("heterosis", "train_phenotype_example.csv"), index=False)
for name, G_, idx in (("train_genotype_example.csv", Gt, ids("Train", ntr)), ("parent_genotype_example.csv", Gp, ids("Parent", npar, 2)),
                      ("test_genotype_example.csv", Gte, ids("Test", 10, 2))):
    d = pd.DataFrame(G_.astype(int), columns=snp, index=idx)
    d.index.name = "Sample"
    d.to_csv(out("heterosis", name))

# ───────────────────────── 14. Phenotypic Data Analysis ─────────────────────────
G_names = ids("Geno", 20, 2)
env = ["Env_2023_Wet", "Env_2023_Dry", "Env_2024_Wet"]
gen_eff = rng.normal(0, 1.2, 20)
env_eff = {e: v for e, v in zip(env, (0, -1.2, 0.5))}
rows = []
for e in env:
    for r in (1, 2, 3):
        for gi, g in enumerate(G_names):
            blk = 1 + (gi // 5)
            yv = 6 + gen_eff[gi] + env_eff[e] + 0.4 * gen_eff[gi] * (env.index(e) - 1) + rng.normal(0, 0.5)
            rows.append([g, e, r, blk, 1 + gi // 5, 1 + gi % 5, f"{r}{gi + 1:02d}", int(e.split("_")[1]), e.split("_")[2],
                         round(yv, 2), round(100 + 6 * gen_eff[gi] + rng.normal(0, 3), 1), round(85 + 2 * gen_eff[gi] + rng.normal(0, 2), 1)])
pd.DataFrame(rows, columns=["Genotype", "Environment", "Rep", "Block", "Row", "Column", "Plot", "Year", "Season",
                            "Yield", "PlantHeight", "DaysToFlower"]).to_csv(out("phenotypic", "multi_environment_trial_example.csv"), index=False)

# The same file layout for the three experimental designs of the Phenotypic Data Analysis tool
PH_COLS = ["Genotype", "Environment", "Rep", "Block", "Row", "Column", "Plot", "Year", "Season",
           "Yield", "PlantHeight", "DaysToFlower"]


def ph_row(g, e, rep, blk, row, col, gen_eff, ge_slope, blk_eff, plot):
    ei = env.index(e)
    yv = 6 + gen_eff + env_eff[e] + ge_slope * gen_eff * (ei - 1) + blk_eff + rng.normal(0, 0.45)
    return [g, e, rep, blk, row, col, plot, int(e.split("_")[1]), e.split("_")[2], round(yv, 2),
            round(100 + 6 * gen_eff + rng.normal(0, 3), 1), round(85 + 2 * gen_eff + rng.normal(0, 2), 1)]


# RCBD: every genotype once in each complete block (replicate); Block repeats Rep
g24 = ids("Geno", 24, 2)
eff24 = dict(zip(g24, rng.normal(0, 1.2, 24)))
rows = []
for e in env:
    for r in (1, 2, 3):
        rep_eff = rng.normal(0, 0.3)
        for k, g in enumerate(rng.permutation(g24)):
            rows.append(ph_row(g, e, r, r, 1 + k // 6, 1 + k % 6, eff24[g], 0.4, rep_eff, f"{r}{k + 1:02d}"))
pd.DataFrame(rows, columns=PH_COLS).to_csv(out("phenotypic", "rcbd_example.csv"), index=False)

# Alpha-lattice: 24 genotypes in 6 incomplete blocks of 4 per replicate, re-randomised in each replicate;
# Block is numbered 1-6 inside every replicate
rows = []
for e in env:
    for r in (1, 2, 3):
        rep_eff = rng.normal(0, 0.3)
        blk_eff = rng.normal(0, 0.5, 6)
        for k, g in enumerate(rng.permutation(g24)):
            b = k // 4
            rows.append(ph_row(g, e, r, b + 1, 1 + b, 1 + k % 4, eff24[g], 0.4, rep_eff + blk_eff[b], f"{r}{k + 1:02d}"))
pd.DataFrame(rows, columns=PH_COLS).to_csv(out("phenotypic", "alpha_lattice_example.csv"), index=False)

# Augmented RCBD: 4 checks in every one of 6 blocks, 36 new entries once each (6 per block); one Rep
checks = ["Check1", "Check2", "Check3", "Check4"]
new_entries = ids("New", 36, 2)
eff_aug = {g: v for g, v in zip(checks + new_entries, rng.normal(0, 1.2, 40))}
rows = []
for e in env:
    order = rng.permutation(new_entries)
    for b in range(6):
        blk_eff = rng.normal(0, 0.6)
        plots = checks + list(order[b * 6:(b + 1) * 6])
        for k, g in enumerate(rng.permutation(plots)):
            rows.append(ph_row(g, e, 1, b + 1, 1 + b, 1 + k, eff_aug[g], 0.4, blk_eff, f"{b + 1}{k + 1:02d}"))
pd.DataFrame(rows, columns=PH_COLS).to_csv(out("phenotypic", "augmented_rcbd_example.csv"), index=False)

# ── Genotype Viewer: two parents and eight RILs of rice (simple A / B / H / N calls) and the same data as a VCF ──
import json as _json
_gv_rng = np.random.default_rng(17)
_lengths = _json.load(open(Path(__file__).resolve().parent.parent / "tools" / "genotype_viewer" / "backend" / "genosee" / "chromosome_length_database.json"))["Oryza_sativa"]
_gv_lines = ["Parent1", "Parent2"] + ids("RIL", 8, 2)
_gv_table = []
for _chrom, _len in _lengths.items():
    _pos = np.sort(_gv_rng.choice(np.arange(200_000, _len - 200_000), size=max(4, int(_len / 4_000_000)), replace=False))
    _calls = np.empty((len(_pos), len(_gv_lines)), dtype="<U1")
    _calls[:, 0], _calls[:, 1] = "A", "B"
    for _j in range(2, len(_gv_lines)):       # a RIL starts with one parent and switches parent at every crossover
        _cuts = np.sort(_gv_rng.uniform(0, _len, _gv_rng.poisson(_len / 25_000_000)))
        _start = int(_gv_rng.integers(0, 2))
        _calls[:, _j] = ["AB"[(_start + int((_cuts < q).sum())) % 2] for q in _pos]
    _noise = _gv_rng.random(_calls.shape)
    _calls[2:][_noise[2:] < 0.03] = "H"        # a few heterozygous and missing calls
    _calls[2:][_noise[2:] > 0.98] = "N"
    for _i, _q in enumerate(_pos):
        _gv_table.append([_chrom, f"{_chrom}_{_q}", int(_q)] + list(_calls[_i]))
_gv = pd.DataFrame(_gv_table, columns=["chr", "marker_name", "pos"] + _gv_lines)
_gv.to_csv(out("genotype_viewer", "rice_RILs_simple_example.csv"), index=False)
_gt = {"A": "0/0", "B": "1/1", "H": "0/1", "N": "./."}
with open(out("genotype_viewer", "rice_RILs_example.vcf"), "w", newline="\n") as f:
    f.write("##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(_gv_lines) + "\n")
    for _r in _gv_table:
        f.write("\t".join([_r[0], str(_r[2]), _r[1], "A", "G", ".", ".", ".", "GT"] + [_gt[c] for c in _r[3:]]) + "\n")

print({p.name: len(list(p.iterdir())) for p in sorted(OUT.iterdir())})
