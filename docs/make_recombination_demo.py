"""Writes Demo_Input_Files/20_Recombination_Analysis/input/mapping_population.csv: a simulated backcross-derived
population with known crossovers (the truth goes to the scratch file `truth` of simulate()).

Chromosomes 1-3, 300 markers each (every gene has two markers, named chr<N>_<pos>_LOC_Os0<N>g<gene>), a recurrent parent
and 80 lines that carry donor segments, plus 2 % genotyping errors, 4 % missing and 1.5 % heterozygous calls.
Run:  python docs/make_recombination_demo.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent.parent / "Demo_Input_Files" / "20_Recombination_Analysis" / "input"
LENGTHS_MB = {1: 42.0, 2: 36.0, 3: 34.0}
N_MARKERS = 300
N_LINES = 80


def simulate(seed=7, err=0.02, miss=0.04, het=0.015):
    rng = np.random.default_rng(seed)
    rows, truth = [], {}
    for c, L in LENGTHS_MB.items():
        genes = np.sort(rng.uniform(0.2, L - 0.2, N_MARKERS // 2))
        pos = np.sort(np.concatenate([genes, genes + rng.uniform(0.001, 0.02, len(genes))]))
        gid = np.repeat(np.arange(1, len(genes) + 1), 2)[np.argsort(np.argsort(np.concatenate([genes, genes + 1e-9])))]
        for p, g in zip(pos, np.sort(np.repeat(np.arange(1, len(genes) + 1), 2))):
            rows.append((c, int(p * 1e6), f"chr{c}_{int(p * 1e6)}_LOC_Os{c:02d}g{g * 10:05d}"))
    markers = pd.DataFrame(rows, columns=["CHROM", "POS", "ID"])
    G = np.zeros((len(markers), N_LINES), dtype=int)
    n_co = np.zeros(N_LINES, dtype=int)
    for j in range(N_LINES):
        for c, L in LENGTHS_MB.items():
            ix = np.flatnonzero(markers["CHROM"].to_numpy() == c)
            p = markers["POS"].to_numpy()[ix] / 1e6
            state = np.zeros(len(ix), dtype=int)
            for _ in range(rng.poisson(1.4)):
                s = rng.uniform(0, L)
                e = s + rng.exponential(5.0)
                state[(p >= s) & (p <= e)] = 2
            n_co[j] += int((np.diff(state) != 0).sum())
            G[ix, j] = state
    truth = n_co.copy()
    flip = rng.random(G.shape) < err
    G = np.where(flip & (G == 0), 2, np.where(flip & (G == 2), 0, G))
    G = np.where(rng.random(G.shape) < het, 1, G)
    G = np.where(rng.random(G.shape) < miss, -1, G)
    return markers, G, truth


def main():
    markers, G, _ = simulate()
    df = markers.copy()
    df.insert(3, "REF", "A")
    df.insert(4, "ALT", "G")
    df["IRRI-PARENT_GT"] = 0
    # the donor is not genotyped; the parent column is homozygous at every marker
    for j in range(G.shape[1]):
        df[f"RIL-{j + 1:03d}_GT"] = G[:, j]
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "mapping_population.csv", index=False)
    print("wrote", OUT / "mapping_population.csv", df.shape)


if __name__ == "__main__":
    main()
