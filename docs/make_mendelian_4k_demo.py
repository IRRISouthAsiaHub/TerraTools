"""Simulated RIL demo for Mendelian Segregation: about 4,500 markers x 150 lines on 12 rice-sized chromosomes.

A proper RIL (F6-F8): every line carries one parent's allele, each chromosome is a mosaic with ~2 breakpoints,
heterozygous calls are rare (residual heterozygosity in short blocks), and about 4 % of the calls are missing (-1).
Most markers segregate 1:1; three segregation-distortion regions (chr 3, chr 6, chr 9) make a real, localised
Mendelian result to find. Run:  python docs/make_mendelian_4k_demo.py
"""
import csv
import random
from pathlib import Path

SEED = 20241119
N_LINES = 150
N_MARKERS = 4500
CHR_MB = [43.3, 35.9, 36.4, 35.5, 29.9, 31.2, 29.7, 28.4, 23.0, 23.2, 29.0, 27.5]   # rice (MSU7), Mb
BREAKS_PER_MB = 0.055          # ~1.6 breakpoints per 30 Mb in an advanced RIL
MISSING = 0.04
HET_BLOCK_PER_CHR = 0.20       # chance that a line carries a short heterozygous block on a chromosome
MISCALL = 0.002
# (chromosome, centre Mb, half-width Mb, weight of the parent-2 allele); parent 1 has weight 1
SDR = [(3, 24.0, 4.0, 0.45), (6, 10.0, 3.5, 0.50), (9, 15.0, 3.0, 1.60)]

OUT = Path(__file__).resolve().parent.parent / "Demo_Input_Files" / "02_Mendelian_Segregation" / "input" / "coded_genotypes_RIL_4k_simulated.csv"


def main():
    rnd = random.Random(SEED)
    per_chr = [round(N_MARKERS * L / sum(CHR_MB)) for L in CHR_MB]
    markers = []                                   # (name, chr, pos_bp)
    for c, (L, n) in enumerate(zip(CHR_MB, per_chr), start=1):
        pos = sorted(rnd.randint(100_000, int(L * 1e6) - 100_000) for _ in range(n))
        markers += [(f"chr{c}_{p}", c, p) for p in pos]
    lines = [f"RIL{i + 1:03d}" for i in range(N_LINES)]
    geno = {m[0]: [] for m in markers}

    for _ in lines:
        calls = {}
        for c, L in enumerate(CHR_MB, start=1):
            cm = [m for m in markers if m[1] == c]
            while True:                                                  # rejection step = selection on the chromosome
                t = 0.0
                breaks = []
                while True:
                    t += rnd.expovariate(BREAKS_PER_MB)
                    if t >= L:
                        break
                    breaks.append(t)
                allele = rnd.choice([0, 2])
                seg = []
                last = 0.0
                for b in breaks + [L]:
                    seg.append((last, b, allele))
                    allele = 2 - allele
                    last = b
                def at(pos_mb):
                    for a, b, al in seg:
                        if a <= pos_mb < b or b == L:
                            return al
                # the allele at each distortion centre decides whether the chromosome is kept
                keep = 1.0
                for sc, mid, half, w2 in SDR:
                    if sc == c:
                        al = at(mid)
                        keep *= (w2 if al == 2 else 1.0) / max(1.0, w2)
                if rnd.random() <= keep:
                    break
            het = None
            if rnd.random() < HET_BLOCK_PER_CHR:
                a = rnd.uniform(0, L - 3)
                het = (a, a + rnd.uniform(0.5, 3.0))
            for name, _c, p in cm:
                mb = p / 1e6
                g = at(mb)
                if het and het[0] <= mb <= het[1]:
                    g = 1
                elif rnd.random() < MISCALL:
                    g = 2 - g
                calls[name] = g
        for name, _c, _p in markers:
            g = calls[name]
            geno[name].append(-1 if rnd.random() < MISSING else g)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Marker", "Chromosome"] + lines)
        for name, c, _p in markers:
            w.writerow([name, c] + geno[name])
    print("wrote", OUT, len(markers), "markers x", N_LINES, "lines")


if __name__ == "__main__":
    main()
