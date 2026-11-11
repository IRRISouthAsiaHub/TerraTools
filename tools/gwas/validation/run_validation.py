"""
Simulation study for the GWAS models: detection, false-positive control, p-value calibration, effect accuracy and stability.

    python -m tools.gwas.validation.run_validation --reps 10 --out validation_report
    python -m tools.gwas.validation.run_validation --quick

Every replicate uses a different seed derived from --seed, so a run is reproducible. For each scenario and model the report gives
  lambda            median and spread over replicates (reported, never used to rescale p-values)
  FWER              share of replicates with at least one Bonferroni hit in a null scenario (target <= alpha)
  type I at 0.05    share of tested markers with raw p <= 0.05 outside the causal markers (target 0.05)
  power             share of causal markers reaching the Bonferroni threshold (single-locus) or selected as a QTN (multi-locus)
  false hits        Bonferroni / QTN hits that are not causal markers, per replicate
  effect bias       mean (estimated - true) / true effect of the causal markers that were detected (allele-count scale)
  stability         standard deviation of power and of lambda across replicates
Markers are independent apart from structure and siblings: linkage disequilibrium inside chromosomes is not simulated.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import time

import numpy as np

from tools.gwas.backend import analysis as an
from tools.gwas.config import GWASConfig
from tools.gwas.validation.simulate import SCENARIOS, simulate

SINGLE = ("OLS", "MLM", "EMMAX", "GEMMA")


def run_one(sim, models, maf, alpha, pcs):
    cfg = GWASConfig(models=list(models), maf=maf, alpha=alpha, permutations=0, ols_pcs=pcs.get("OLS", 0), mlm_pcs=pcs.get("MLM", 0),
                     emmax_pcs=pcs.get("EMMAX", 0), gemma_pcs=pcs.get("GEMMA", 0), compare_pcs="", n_pcs=5, top_n=5)
    A = an.analyse(cfg, sim.data, labels=sim.labels, log=lambda m: None)
    names = {str(s): i for i, s in enumerate(A.snp)}
    truth = [names.get(n) for n in sim.truth["causal_names"]]
    rows = []
    for model, r in A.results.items():
        sig = r.info["sig"]
        p = np.asarray(r.p, float)
        ok = np.isfinite(p) & (p > 0)
        hit = sig["genome_wide"]
        is_c = np.zeros(len(p), dtype=bool)
        for t in truth:
            if t is not None:
                is_c[t] = True
        found = [t for t in truth if t is not None and hit[t]]
        bias = []
        for k, t in enumerate(truth):
            if t is not None and hit[t] and np.isfinite(r.beta[t]):
                true = sim.causal_beta[k]
                bias.append((r.beta[t] - true) / abs(true))
        rows.append({"model": model, "lambda": r.info["diag"]["lambda"], "tested": int(ok.sum()),
                     "type1_05": float((p[ok & ~is_c] <= 0.05).mean()), "false_hits": int((hit & ~is_c).sum()),
                     "any_false_hit": bool((hit & ~is_c).any()), "causal": len(truth), "found": len(found),
                     "bias": float(np.mean(bias)) if bias else float("nan"), "failed": False})
    for model, msg in A.errors.items():
        rows.append({"model": model, "failed": True, "error": msg})
    return rows, A.qc.ledger[-1]["snps_left"], A.structure


def summarise(rows):
    out = []
    for model in sorted({r["model"] for r in rows}):
        rr = [r for r in rows if r["model"] == model]
        good = [r for r in rr if not r["failed"]]
        if not good:
            out.append({"model": model, "replicates": len(rr), "failed": len(rr), "error": rr[0].get("error", "")})
            continue
        lam = np.array([r["lambda"] for r in good], float)
        c = sum(r["causal"] for r in good)
        f = sum(r["found"] for r in good)
        per_rep_power = np.array([r["found"] / r["causal"] for r in good if r["causal"]], float)
        bias = np.array([r["bias"] for r in good], float)
        out.append({"model": model, "replicates": len(rr), "failed": len(rr) - len(good), "lambda_median": float(np.nanmedian(lam)),
                    "lambda_sd": float(np.nanstd(lam)), "type1_at_0.05": float(np.mean([r["type1_05"] for r in good])),
                    "FWER": float(np.mean([r["any_false_hit"] for r in good])),
                    "false_hits_per_rep": float(np.mean([r["false_hits"] for r in good])),
                    "power": (f / c) if c else float("nan"), "power_sd": float(per_rep_power.std()) if len(per_rep_power) else float("nan"),
                    "effect_bias": float(np.nanmean(bias)) if np.isfinite(bias).any() else float("nan")})
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--seed", type=int, default=2024)
    ap.add_argument("--models", default="OLS,MLM,EMMAX,GEMMA,FarmCPU,BLINK")
    ap.add_argument("--scenarios", default=",".join(SCENARIOS))
    ap.add_argument("--maf", default="0.01,0.05", help="MAF thresholds to compare")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--out", default="gwas_validation")
    ap.add_argument("--quick", action="store_true", help="3 replicates, fewer models, smaller panel")
    a = ap.parse_args(argv)
    models = a.models.split(",")
    reps = a.reps
    if a.quick:
        reps, models = 3, ["OLS", "MLM", "EMMAX"]
    os.makedirs(a.out, exist_ok=True)
    allrows, table = [], []
    for sname in a.scenarios.split(","):
        sc = SCENARIOS[sname]
        if a.quick:
            sc = type(sc)(**{**sc.__dict__, "m": 1200, "n_hidden": 300})
        for maf in [float(x) for x in a.maf.split(",")]:
            for pc_label, pcs in (("no PCs", {}), ("3 PCs", {k: 3 for k in SINGLE})):
                rows = []
                t0 = time.time()
                for rep in range(reps):
                    sim = simulate(sc, a.seed + 1000 * rep + 17)
                    r, kept, _ = run_one(sim, models, maf, a.alpha, pcs)
                    for x in r:
                        x.update(scenario=sname, maf=maf, rep=rep, pcs=pc_label)
                    rows += r
                for s in summarise(rows):
                    s.update(scenario=sname, maf_threshold=maf, covariates=pc_label)
                    table.append(s)
                allrows += rows
                print(f"{sname:24s} maf>={maf:<5g} {pc_label:7s} {time.time() - t0:6.1f}s", flush=True)
    with open(os.path.join(a.out, "replicates.json"), "w") as fh:
        json.dump(allrows, fh, indent=1, default=float)
    cols = ["scenario", "maf_threshold", "covariates", "model", "replicates", "failed", "lambda_median", "lambda_sd", "type1_at_0.05",
            "FWER", "false_hits_per_rep", "power", "power_sd", "effect_bias", "error"]
    with open(os.path.join(a.out, "summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(table)
    print(f"\nWrote {a.out}/summary.csv")
    return table


if __name__ == "__main__":
    main()
