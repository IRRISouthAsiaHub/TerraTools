"""
genomic_prediction_cli.py — subprocess entrypoint for the BreedSight
genomic prediction pipeline.

CRITICAL: this file (and therefore TensorFlow) must NEVER be imported
directly into the main GUI process. TensorFlow bundles its own copy of
LLVM; Qt WebEngine (used by every TerraTools window) bundles a different
one. Both register the same global LLVM command-line flags at import
time, and the second registration triggers a hard, uncatchable
LLVM_ERROR abort() that kills the entire application — not just this
tool. Running the pipeline in a separate `python` process (this script,
launched via subprocess by genomic_prediction_runner.py) keeps
TensorFlow's native library loading completely isolated from Qt's.

Prints each log line to stdout as it happens (flushed immediately) so the
parent process can stream progress live, then writes the full JSON result
to <outdir>/summary.json once the pipeline finishes.
"""

import argparse
import json
import sys
from pathlib import Path

# Importing the engine (and therefore TensorFlow) only happens here, inside
# this standalone subprocess — never inside the GUI process.
from tools.genomic_prediction.backend.genomic_prediction_engine import run_breedsight_pipeline


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train_pheno", required=True)
    ap.add_argument("--train_geno", required=True)
    ap.add_argument("--test_pheno", required=True)
    ap.add_argument("--test_geno", required=True)
    ap.add_argument("--outer_n_splits", type=int, required=True)
    ap.add_argument("--epochs", type=int, required=True)
    ap.add_argument("--batch_size", type=int, required=True)
    ap.add_argument("--learning_rate", type=float, required=True)
    ap.add_argument("--l2_reg", type=float, required=True)
    ap.add_argument("--dropout_rate", type=float, required=True)
    ap.add_argument("--rf_n_estimators", type=int, required=True)
    ap.add_argument("--rf_max_depth", type=int, required=True)
    ap.add_argument("--alpha", type=float, required=True)
    ap.add_argument("--heritability", type=float, required=True)
    ap.add_argument("--use_pca", type=str, required=True)
    ap.add_argument("--n_components", type=int, required=True)
    ap.add_argument("--feature_selection", type=str, required=True)
    ap.add_argument("--rfe_n_features", type=int, required=True)
    ap.add_argument("--use_raw_genotypes", type=str, required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    summary_path = outdir / "summary.json"

    def as_bool(s):
        return str(s).strip().lower() in ("true", "1", "yes")

    def on_log(msg):
        print(msg, flush=True)

    try:
        result = run_breedsight_pipeline(
            train_pheno_file=args.train_pheno,
            train_geno_file=args.train_geno,
            test_pheno_file=args.test_pheno,
            test_geno_file=args.test_geno,
            outer_n_splits=args.outer_n_splits,
            epochs=args.epochs,
            batch_size=args.batch_size,
            learning_rate=args.learning_rate,
            l2_reg=args.l2_reg,
            dropout_rate=args.dropout_rate,
            rf_n_estimators=args.rf_n_estimators,
            rf_max_depth=args.rf_max_depth,
            alpha=args.alpha,
            heritability=args.heritability,
            use_pca=as_bool(args.use_pca),
            n_components=args.n_components,
            feature_selection=as_bool(args.feature_selection),
            rfe_n_features=args.rfe_n_features,
            use_raw_genotypes=as_bool(args.use_raw_genotypes),
            outdir=outdir,
            on_log=on_log,
        )
    except Exception as e:
        import traceback
        result = {"status": "error", "message": f"{e}\n\n{traceback.format_exc()}"}

    with open(summary_path, "w") as f:
        json.dump(result, f)

    sys.exit(0 if result.get("status") == "success" else 1)


if __name__ == "__main__":
    main()
