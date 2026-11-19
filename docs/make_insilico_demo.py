"""Simulated F2 population for the In silico BSA demo and tests (known answers).

300 F2 plants, 3 chromosomes (100 cM, 20 Mb each), 300 markers per chromosome, 1 % genotyping errors, 3 % missing calls.
Three QTL, each explaining about 8 % of the trait variance:
  chr1 at 35 cM: additive-dominant (d = 0.5 a)      - found by every BSA statistic
  chr2 at 50 cM: FULLY over-dominant (a = 0, d != 0) - equal allele frequencies in the pools: only Gw can see it
  chr3 at 70 cM: completely dominant (d = a)         - found by every BSA statistic
Run:  python docs/make_insilico_demo.py   -> Demo_Input_Files/21_InSilico_BSA/input/
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
N, NCHR, NM, LEN_CM, LEN_BP = 300, 3, 300, 100.0, 20_000_000
QTL = [(1, 35.0, 1.0, 0.5, "additive-dominant"), (2, 50.0, 0.0, 1.4, "fully over-dominant"), (3, 70.0, 1.0, 1.0, "completely dominant")]


def simulate(seed=2022, n=N):
    rng = np.random.default_rng(seed)
    cm = np.linspace(0, LEN_CM, NM)
    pos = np.round(np.linspace(100_000, LEN_BP - 100_000, NM)).astype(int)
    d = np.diff(cm) / 100.0
    r = 0.5 * (1 - np.exp(-2 * d))
    markers, cols = [], []
    for c in range(1, NCHR + 1):
        def gamete():
            h = np.empty((n, NM), dtype=np.int8)
            h[:, 0] = rng.integers(0, 2, n)
            flips = rng.random((n, NM - 1)) < r
            h[:, 1:] = (h[:, [0]] + np.cumsum(flips, axis=1)) % 2
            return h
        cols.append(gamete() + gamete())                              # 0, 1, 2 = copies of the parent-2 allele
        markers.append(pd.DataFrame({"CHROM": c, "POS": pos, "ID": [f"chr{c}_{p}" for p in pos]}))
    G = np.hstack(cols)                                                # lines x markers
    markers = pd.concat(markers, ignore_index=True)
    y = np.zeros(n)
    truth = []
    for c, qcm, a, dd, kind in QTL:
        j = (c - 1) * NM + int(np.argmin(np.abs(cm - qcm)))
        g = G[:, j]
        y += np.where(g == 0, -a, np.where(g == 1, dd, a))
        truth.append({"chromosome": c, "position_bp": int(pos[j - (c - 1) * NM]), "type": kind})
    var_q = y.var()
    y = y + rng.normal(0, np.sqrt(var_q * 1.6), n)                      # three QTL together explain ~ 38 % of the variance
    y = 100 + 10 * (y - y.mean()) / y.std()
    obs = G.copy()
    err = rng.random(G.shape) < 0.01
    obs[err] = rng.integers(0, 3, err.sum())
    obs[rng.random(G.shape) < 0.03] = -1
    return markers, obs.T.astype(np.int8), y, truth


def main():
    markers, G, y, truth = simulate()
    out = ROOT / "Demo_Input_Files" / "21_InSilico_BSA" / "input"
    out.mkdir(parents=True, exist_ok=True)
    names = [f"F2-{i + 1:03d}" for i in range(G.shape[1])]
    df = pd.concat([markers, pd.DataFrame(G, columns=names)], axis=1)
    df.to_csv(out / "genotypes_F2.csv", index=False)
    pd.DataFrame({"Sample": names, "PlantHeight": np.round(y, 2)}).to_csv(out / "phenotypes.csv", index=False)
    print("written", out, truth)


if __name__ == "__main__":
    main()
