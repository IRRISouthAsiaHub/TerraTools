"""
Write one simulated data set in the formats the reference programs read, plus this program's own p-values, so every program sees
the same markers, samples, phenotype, relationship matrix and covariates.

    python -m tools.gwas.benchmark.make_inputs --scenario causal_3 --seed 7 --out benchmark_run

Files written to --out
  numeric_genotypes.csv   samples x markers, 0/1/2, NA missing (QC-passed, mean-imputed copy used for the scan)
  map.csv                 SNP, Chromosome, Position
  phenotype.csv           Taxa, Trait
  kinship.csv             VanRaden relationship matrix (same samples, same order)
  pcs.csv                 first five principal components
  truth.json              causal markers and effects
  terratools_<model>.csv  SNP, P (this program's p-values for OLS, MLM, EMMAX, GEMMA, FarmCPU, BLINK)
"""

import argparse
import json
import os

import numpy as np
import pandas as pd

from tools.gwas.backend import analysis as an
from tools.gwas.config import GWASConfig
from tools.gwas.validation.simulate import SCENARIOS, simulate


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", default="causal_3", choices=None)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--pcs", type=int, default=3)
    ap.add_argument("--out", default="benchmark_run")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    sim = simulate(SCENARIOS[a.scenario], a.seed)
    cfg = GWASConfig(models=["OLS", "MLM", "EMMAX", "GEMMA", "FarmCPU", "BLINK"], permutations=0, ols_pcs=a.pcs, mlm_pcs=a.pcs,
                     emmax_pcs=a.pcs, gemma_pcs=a.pcs, compare_pcs="", n_pcs=5)
    A = an.analyse(cfg, sim.data, labels=sim.labels, log=lambda m: None)
    pd.DataFrame(A.G, index=A.ids, columns=A.snp).to_csv(os.path.join(a.out, "numeric_genotypes.csv"), index_label="Taxa")
    pd.DataFrame({"SNP": A.snp, "Chromosome": A.chrom, "Position": A.pos}).to_csv(os.path.join(a.out, "map.csv"), index=False)
    pd.DataFrame({"Taxa": A.ids, "Trait": A.y}).to_csv(os.path.join(a.out, "phenotype.csv"), index=False)
    pd.DataFrame(A.K, index=A.ids, columns=A.ids).to_csv(os.path.join(a.out, "kinship.csv"), index_label="Taxa")
    pcs = np.asarray(A.pc["scores"])[:, :5]
    pd.DataFrame(pcs, index=A.ids, columns=[f"PC{i + 1}" for i in range(pcs.shape[1])]).to_csv(os.path.join(a.out, "pcs.csv"), index_label="Taxa")
    json.dump({"scenario": a.scenario, "seed": a.seed, "causal": sim.truth["causal_names"], "effects": sim.causal_beta.tolist(),
               "pcs_used": a.pcs}, open(os.path.join(a.out, "truth.json"), "w"), indent=1)
    for name, r in A.results.items():
        pd.DataFrame({"SNP": A.snp, "P": r.p}).to_csv(os.path.join(a.out, f"terratools_{name}.csv"), index=False)
    print(f"Wrote inputs and this program's results to {a.out}")


if __name__ == "__main__":
    main()
