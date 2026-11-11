"""
Compare this program's p-values with a reference program's, marker by marker, against stated tolerances.

    python -m tools.gwas.benchmark.compare benchmark_run --reference GAPIT_MLM=gapit_mlm.csv:SNP:P.value --ours MLM

Each --reference is label=file:snp_column:p_column. Programs that use different model details (for example PCs, kinship centring or
variance-component estimation) are not expected to agree exactly; the tolerances below are for the same model on identical inputs:

  Spearman correlation of -log10 p   >= 0.99
  largest |difference in -log10 p|   <= 0.5
  overlap of the top 1% (min 5) markers >= 0.8

A reference that is not installed is reported as "not run"; nothing is substituted for it.
"""

import argparse
import json
import os
import sys

import pandas as pd

from tools.gwas.backend.gemma_external import TOLERANCE, compare


def load(path, snp_col, p_col):
    t = pd.read_csv(path)
    return dict(zip(t[snp_col].astype(str), pd.to_numeric(t[p_col], errors="coerce")))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--ours", required=True, help="model name of this program's results to use, e.g. MLM")
    ap.add_argument("--reference", action="append", required=True, help="label=file:snp_column:p_column")
    a = ap.parse_args(argv)
    mine = load(os.path.join(a.folder, f"terratools_{a.ours}.csv"), "SNP", "P")
    rows = []
    for spec in a.reference:
        label, rest = spec.split("=", 1)
        parts = rest.rsplit(":", 2)
        path, sc, pc = parts if len(parts) == 3 else (rest, "SNP", "P")
        if not os.path.isfile(path):
            rows.append({"reference": label, "status": "not run", "message": f"{path} does not exist; run the reference program first (see README.md)"})
            continue
        ref = load(path, sc, pc)
        common = sorted(set(mine) & set(ref))
        rep = compare([mine[s] for s in common], [ref[s] for s in common])
        rep.update(reference=label, ours=a.ours, status="compared", tolerances=TOLERANCE)
        rows.append(rep)
    print(json.dumps(rows, indent=1))
    return 0 if all(r.get("agrees") for r in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
