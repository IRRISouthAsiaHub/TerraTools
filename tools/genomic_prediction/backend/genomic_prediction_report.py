"""
genomic_prediction_report.py — runs the BreedSight engine (BreedSight_all_models_final_v2.py) and turns what
it wrote into the result the TerraTools screen shows.

The engine trains BreedSight and, only when asked (baselines=[...]), Lasso / RRBLUP_1 / GBLUP_1 for comparison.
This module only calls run_complete_analysis(), then reads
its CSV files and diagnostic plots, adds two overview plots and the download list. Import it only inside the prediction subprocess (it pulls in TensorFlow).
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd


def _records(df, n=None):
    if n:
        df = df.head(n)
    return df.astype(object).where(df.notna(), None).to_dict("records")


# ── two overview plots (the engine writes the per-fold regression plots and the R² comparison) ──────

def _plot_fold_r2(outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.5, 5), dpi=120)
    found = False
    for f in sorted(outdir.glob("*_detailed_results.csv")):
        d = pd.read_csv(f)
        if {"Fold", "Val_R2"} <= set(d.columns):
            ax.plot(d["Fold"], d["Val_R2"], "-o", ms=4, lw=1.4, label=f.name.replace("_detailed_results.csv", ""))
            found = True
    if not found:
        plt.close(fig)
        return None
    ax.set_xlabel("Cross-validation fold")
    ax.set_ylabel("Validation R²")
    ax.set_title("Validation R² per fold, all models")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    path = outdir / "plot_folds.png"
    fig.savefig(path, facecolor="white")
    plt.close(fig)
    return str(path)


def _plot_test_grid(test_df, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if "True_Phenotype" not in test_df.columns:
        return None
    cols = [c for c in test_df.columns if c.startswith("Predicted_Phenotype_")]
    if not cols:
        return None
    n = len(cols)
    ncol = min(4, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 3.4 * nrow), dpi=110, squeeze=False)
    y = test_df["True_Phenotype"].to_numpy(float)
    for ax, c in zip(axes.ravel(), cols):
        p = test_df[c].to_numpy(float)
        ax.scatter(y, p, s=18, alpha=0.7, color="#2563EB", edgecolor="none")
        lo, hi = float(min(y.min(), p.min())), float(max(y.max(), p.max()))
        ax.plot([lo, hi], [lo, hi], "r--", lw=1)
        r = np.corrcoef(y, p)[0, 1] if len(y) > 2 else float("nan")
        ax.set_title(f"{c.replace('Predicted_Phenotype_', '')}  (r = {r:.2f})", fontsize=9)
        ax.set_xlabel("True", fontsize=8)
        ax.set_ylabel("Predicted", fontsize=8)
        ax.grid(alpha=0.3)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    fig.suptitle("Test set: predicted versus true phenotype", fontsize=11)
    fig.tight_layout()
    path = outdir / "plot_test_grid.png"
    fig.savefig(path, facecolor="white")
    plt.close(fig)
    return str(path)


# ── input check ────────────────────────────────────────────────────────────

_HOW_TO_BUILD = ("Build them with the Clustering tool (training / testing files, with a phenotype file) and the "
                 "Additive Dominance tool (Training + Testing mode).")


def _header(path):
    return [str(c) for c in pd.read_csv(path, nrows=0).columns]


def _ids(path):
    return pd.read_csv(path, usecols=[0]).iloc[:, 0].astype(str).str.strip().tolist()


def _same_lines(name_a, a, name_b, b):
    """BreedSight pairs the rows of two files by position, so they must list the same lines in the same order."""
    if len(a) != len(b):
        raise ValueError(f"The {name_a} has {len(a)} lines but the {name_b} has {len(b)}. "
                         "They must hold the same lines, in the same order.")
    if a != b:
        i = next(k for k, (x, y) in enumerate(zip(a, b)) if x != y)
        if set(a) == set(b):
            raise ValueError(f"The {name_a} and the {name_b} hold the same lines in a different order "
                             f"(row {i + 1}: {a[i]} / {b[i]}). BreedSight matches them by row, so sort them the same way.")
        raise ValueError(f"The {name_a} and the {name_b} do not hold the same lines "
                         f"(row {i + 1}: {a[i]} / {b[i]}). {_HOW_TO_BUILD}")


def check_inputs(train_pheno_file, train_geno_file, test_pheno_file, test_geno_file):
    """Plain-language checks of the four files before the (slow) engine starts. The engine pairs rows by
    position and only compares marker names as sets, so a mismatch would otherwise fail late or, worse,
    run on misaligned data."""
    if "phenotypes" not in _header(train_pheno_file):
        raise ValueError("The training phenotypic file has no column called 'phenotypes'. It needs sample_id and "
                         "phenotypes (the Training file of the Clustering tool has it when a phenotype file was chosen there).")

    ids = {"tp": _ids(train_pheno_file), "tg": _ids(train_geno_file), "sp": _ids(test_pheno_file), "sg": _ids(test_geno_file)}
    _same_lines("training phenotypic file", ids["tp"], "training genotypic/additive file", ids["tg"])
    _same_lines("testing phenotypic file", ids["sp"], "testing genotypic/additive file", ids["sg"])

    tr_cols, te_cols = _header(train_geno_file)[1:], _header(test_geno_file)[1:]
    if tr_cols != te_cols:
        common = len(set(tr_cols) & set(te_cols))
        msg = (f"The testing genotypic/additive file must have exactly the same columns as the training one, in the same "
               f"order. Training has {len(tr_cols)} columns, testing has {len(te_cols)}, {common} in common. ")
        if common == len(tr_cols) == len(te_cols):
            msg += "The columns are the same but in a different order."
        elif te_cols and set(te_cols) == set(ids["sg"]):
            msg += ("The columns of the testing file are the testing lines themselves: that matrix was computed on the testing "
                    "lines alone. A relationship matrix needs testing lines as rows and the TRAINING lines as columns, taken "
                    "from one matrix computed over all lines. ")
        elif te_cols and common == len(te_cols) < len(tr_cols):
            msg += (f"The testing file has only {len(te_cols)} of the {len(tr_cols)} training columns, so the other "
                    f"{len(tr_cols) - len(te_cols)} are missing. Every testing line needs a value for every training column "
                    "(for a relationship matrix: one column per training line). The two files were probably not made in the "
                    "same run, or the testing file was cut down. ")
        raise ValueError(msg.rstrip() + " " + _HOW_TO_BUILD)

    both = set(ids["tg"]) & set(ids["sg"])
    if both:
        print(f"WARNING: {len(both)} line(s) are in both the training and the testing file (e.g. {sorted(both)[0]}). "
              "Their accuracy on the test set will look better than it is.", flush=True)
    print(f"Input check passed: {len(ids['tg'])} training lines, {len(ids['sg'])} testing lines, {len(tr_cols)} columns.", flush=True)


# ── entry point ────────────────────────────────────────────────────────────

def run_breedsight(train_pheno_file, train_geno_file, test_pheno_file, test_geno_file, outdir, *,
                   outer_n_splits, epochs, batch_size, learning_rate, l2_reg, dropout_rate, rf_n_estimators,
                   rf_max_depth, alpha, heritability, use_pca, n_components, feature_selection, rfe_n_features,
                   use_raw_genotypes, baselines=()):
    check_inputs(train_pheno_file, train_geno_file, test_pheno_file, test_geno_file)
    outdir = Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    # the engine writes its diagnostic plots to the relative folder "output1/diagnostic_plots"
    os.chdir(outdir)

    from tools.genomic_prediction.backend.BreedSight_all_models_final_v2 import run_complete_analysis

    run_complete_analysis(
        training_file=train_pheno_file, training_additive_file=train_geno_file,
        testing_file=test_pheno_file, testing_additive_file=test_geno_file,
        output_dir=str(outdir),
        outer_n_splits=int(outer_n_splits), epochs=int(epochs), batch_size=int(batch_size),
        learning_rate=float(learning_rate), l2_reg=float(l2_reg), dropout_rate=float(dropout_rate),
        rf_n_estimators=int(rf_n_estimators), rf_max_depth=None if int(rf_max_depth) <= 0 else int(rf_max_depth),
        alpha=float(alpha), heritability=float(heritability), use_pca=bool(use_pca),
        n_components=int(n_components), feature_selection=bool(feature_selection),
        rfe_n_features=int(rfe_n_features), use_raw_genotypes=bool(use_raw_genotypes), generate_plots=True,
        baselines=list(baselines or []))

    metrics = pd.read_csv(outdir / "model_metrics.csv")
    test_df = pd.read_csv(outdir / "testing_predictions.csv")
    best = metrics.loc[metrics["Avg_Val_R2"].idxmax()]

    plots = {}
    diag = outdir / "output1" / "diagnostic_plots"
    bs = pd.read_csv(outdir / "BreedSight_detailed_results.csv")
    best_fold = int(bs.loc[bs["Val_R2"].idxmax(), "Fold"])
    reg = diag / f"validation_BreedSight_fold_{best_fold}_regression.png"
    plots["regression"] = str(reg) if reg.exists() else None
    cmp_ = diag / "model_comparison_summary.png"
    plots["comparison"] = str(cmp_) if cmp_.exists() else None
    for key, fn in (("folds", lambda: _plot_fold_r2(outdir)), ("test_grid", lambda: _plot_test_grid(test_df, outdir))):
        try:
            plots[key] = fn()
        except Exception as e:                       # a missing overview plot must not lose the analysis
            print(f"Could not draw the {key} plot: {e}", flush=True)
            plots[key] = None

    downloads = {k: str(outdir / f) for k, f in (
        ("metrics_csv", "model_metrics.csv"), ("test_csv", "testing_predictions.csv"),
        ("val_csv", "validation_predictions.csv"), ("test_metrics_csv", "test_metrics.csv"),
        ("summary_txt", "analysis_summary.txt")) if (outdir / f).exists()}

    test_metrics = None
    if (outdir / "test_metrics.csv").exists():
        test_metrics = _records(pd.read_csv(outdir / "test_metrics.csv"))

    return {
        "status": "success",
        "metrics_table": _records(metrics),
        "test_metrics_table": test_metrics,
        "test_predictions_table": _records(test_df, 100),
        "best_model": str(best["Model"]),
        "best_val_r2": float(best["Avg_Val_R2"]),
        "plots": plots,
        "downloads": downloads,
    }
