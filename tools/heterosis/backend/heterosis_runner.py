# -*- coding: utf-8 -*-
"""
heterosis_runner.py

Backend computation for the Heterosis Prediction app.

Ported from the original Gradio implementation (Sequential Keras DNN ->
GEBV prediction -> mid-parent / better-parent heterosis screening) into
the same "runner" interface used by qtl_runner.py / linkage_runner.py,
so it can be driven from a PySide6 QWebChannel bridge instead of Gradio.

Differences from the original Gradio app:
  * No `gr.File` / `gr.Progress` objects — plain file paths in, a
    `progress_callback(fraction, message)` out.
  * Every DataFrame returned to the UI is converted to
    `{"rows": [...], "columns": [...]}` (JSON-safe: no numpy scalars).
  * All artifacts (plot PNG, result CSVs) are written under the
    `results_dir` the bridge passes in, and their paths are returned
    in the result dict rather than as Gradio component values.
"""

import os
import random
import tempfile

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import pearsonr

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, BatchNormalization, LeakyReLU
from tensorflow.keras import regularizers
from tensorflow.keras.optimizers import Adam
import tensorflow as tf


# ============================================================
# REPRODUCIBILITY
# ============================================================

SEED = 42
os.environ["PYTHONHASHSEED"] = str(SEED)
os.environ["TF_DETERMINISTIC_OPS"] = "1"
random.seed(SEED)
np.random.seed(SEED)
tf.random.set_seed(SEED)


# ============================================================
# MODEL
# ============================================================

def build_model(input_dim, l1_reg=0.001, l2_reg=0.001, dropout_rate=0.2):

    model = Sequential()

    for units in [512, 512, 256, 128, 64, 32, 16, 8]:

        model.add(Dense(
            units,
            input_shape=(input_dim,) if units == 512 else [],
            kernel_initializer="he_normal",
            kernel_regularizer=regularizers.l1_l2(l1=l1_reg, l2=l2_reg)
        ))

        model.add(BatchNormalization())
        model.add(Dropout(dropout_rate))
        model.add(LeakyReLU(alpha=0.1))

    model.add(Dense(1, activation="relu"))

    model.compile(
        loss="mse",
        optimizer=Adam(),
        metrics=["mae"]
    )

    return model


# ============================================================
# HELPERS
# ============================================================

def _df_to_records(df, decimals=4):
    """Convert a DataFrame into JSON-safe {rows, columns}."""

    out = df.copy()

    for col in out.select_dtypes(include=[np.floating]).columns:
        out[col] = out[col].round(decimals)

    out = out.replace({np.nan: None})

    return {
        "columns": list(out.columns),
        "rows": out.to_dict(orient="records")
    }


def _to_tmp_csv(df, name, results_dir):

    path = os.path.join(results_dir, name)
    df.to_csv(path, index=False)

    return path


def load_preview(path, n_rows=10):
    """Quick preview of an uploaded CSV — used for all four file inputs."""

    try:

        df = pd.read_csv(path)

        preview = df.head(n_rows).replace({np.nan: None})

        return {
            "status": "success",
            "filename": os.path.basename(path),
            "n_rows": int(df.shape[0]),
            "n_cols": int(df.shape[1]),
            "preview_columns": list(df.columns),
            "preview_rows": preview.to_dict(orient="records")
        }

    except Exception as e:

        return {
            "status": "error",
            "message": str(e)
        }


# ============================================================
# MAIN ENTRY POINT
# ============================================================

def run_heterosis(config, results_dir, progress_callback=None):
    """
    config keys:
      train_pheno_file, train_geno_file, parent_geno_file  (paths, required)
      test_geno_file                                       (path, optional)
      num_F1, num_P                                        (int)
      l1_reg, l2_reg, dropout_rate                          (float)
      k_folds, epochs, batch_size                            (int)

    progress_callback(fraction: float, message: str) is called
    periodically so the UI can show live status / a progress bar.
    """

    def notify(fraction, message):
        if progress_callback:
            try:
                progress_callback(fraction, message)
            except Exception:
                pass

    logs = []

    def log(msg):
        logs.append(msg)
        notify(None, msg)
        return "\n".join(logs)

    os.makedirs(results_dir, exist_ok=True)

    # --------------------------------------------------------
    # Load data
    # --------------------------------------------------------

    try:

        train_pheno = pd.read_csv(config["train_pheno_file"])
        train_geno = pd.read_csv(config["train_geno_file"], index_col=0)
        parent_geno_df = pd.read_csv(config["parent_geno_file"], index_col=0)

    except Exception as e:

        return {
            "status": "error",
            "message": f"File loading error: {e}",
            "logs": "\n".join(logs)
        }

    parent_index = parent_geno_df.index

    train_Y = np.array(pd.to_numeric(train_pheno.iloc[:, 1], errors="coerce"))
    train_X = train_geno.to_numpy(dtype=np.float32)
    parent_X_raw = parent_geno_df.to_numpy(dtype=np.float32)

    scaler = StandardScaler()
    train_X = scaler.fit_transform(train_X)
    parent_X = scaler.transform(parent_X_raw)

    # --------------------------------------------------------
    # K-Fold cross-validation
    # --------------------------------------------------------

    k_folds = int(config.get("k_folds", 5))
    epochs = int(config.get("epochs", 100))
    batch_size = int(config.get("batch_size", 16))
    l1_reg = float(config.get("l1_reg", 0.001))
    l2_reg = float(config.get("l2_reg", 0.001))
    dropout_rate = float(config.get("dropout_rate", 0.2))
    num_F1 = int(config.get("num_F1", 5))
    num_P = int(config.get("num_P", 5))

    kfold = KFold(n_splits=k_folds, shuffle=True, random_state=60)
    fold_metrics = []
    fold_records = []

    log(f"Starting {k_folds}-fold cross-validation (epochs={epochs}, batch={batch_size})")

    for fold, (train_idx, val_idx) in enumerate(kfold.split(train_X)):

        notify(fold / k_folds, f"Fold {fold + 1}/{k_folds}")
        log(f"  Fold {fold + 1}/{k_folds} ...")

        model = build_model(train_X.shape[1], l1_reg, l2_reg, dropout_rate)

        model.fit(
            train_X[train_idx], train_Y[train_idx],
            validation_data=(train_X[val_idx], train_Y[val_idx]),
            epochs=epochs, batch_size=batch_size, verbose=0
        )

        val_pred = model.predict(train_X[val_idx], verbose=0).flatten()

        mse = mean_squared_error(train_Y[val_idx], val_pred)
        rmse = float(np.sqrt(mse))
        r2 = r2_score(train_Y[val_idx], val_pred)
        pearson_corr, _ = pearsonr(train_Y[val_idx], val_pred)

        fold_metrics.append((mse, rmse, r2, pearson_corr))

        fold_records.append({
            "Fold": fold + 1,
            "MSE": round(float(mse), 4),
            "RMSE": round(rmse, 4),
            "R2": round(float(r2), 4),
            "Pearson_r": round(float(pearson_corr), 4)
        })

        log(f"    MSE={mse:.4f}  RMSE={rmse:.4f}  R2={r2:.4f}  Pearson r={pearson_corr:.4f}")

    avg_mse = float(np.mean([m[0] for m in fold_metrics]))
    avg_rmse = float(np.mean([m[1] for m in fold_metrics]))
    avg_r2 = float(np.mean([m[2] for m in fold_metrics]))
    avg_pearson = float(np.mean([m[3] for m in fold_metrics]))

    log(f"Average CV -> MSE={avg_mse:.4f}  RMSE={avg_rmse:.4f}  R2={avg_r2:.4f}  Pearson r={avg_pearson:.4f}")

    # --------------------------------------------------------
    # Final model on full training data
    # --------------------------------------------------------

    notify(0.85, "Training final model...")
    log("Training final model on full dataset...")

    final_model = build_model(train_X.shape[1], l1_reg, l2_reg, dropout_rate)
    final_model.fit(train_X, train_Y, epochs=epochs, batch_size=batch_size, verbose=0)

    # --------------------------------------------------------
    # Parent GEBVs
    # --------------------------------------------------------

    parent_GEBV_values = final_model.predict(parent_X, verbose=0).flatten()
    parent_GEBV_dict = dict(zip(parent_index, parent_GEBV_values))

    # --------------------------------------------------------
    # Hybrid predictions (all pairwise combinations)
    # --------------------------------------------------------

    hybrid_geno, F1_names, SCA_list = [], [], []
    Parent1_GEBV_list, Parent2_GEBV_list = [], []

    for i in range(len(parent_X) - 1):
        for j in range(i + 1, len(parent_X)):

            hybrid = (parent_X[i] + parent_X[j]) / 2
            hybrid_geno.append(hybrid)
            F1_names.append(f"{parent_index[i]} : {parent_index[j]}")

            P1 = parent_GEBV_dict[parent_index[i]]
            P2 = parent_GEBV_dict[parent_index[j]]

            Parent1_GEBV_list.append(P1)
            Parent2_GEBV_list.append(P2)
            SCA_list.append((P1 + P2) / 2)

    hybrid_geno = np.array(hybrid_geno)

    notify(0.95, f"Scoring {len(hybrid_geno)} hybrid combinations...")

    GEBV_pred = final_model.predict(hybrid_geno, verbose=0).flatten()

    parent_mean = np.mean(train_Y)
    better_parent = np.max(train_Y)
    MPH = ((GEBV_pred - parent_mean) / parent_mean) * 100
    BPH = ((GEBV_pred - better_parent) / better_parent) * 100

    hybrid_results = pd.DataFrame({
        "Hybrid": F1_names,
        "GEBV": GEBV_pred,
        "MPH": MPH,
        "BPH": BPH,
        "Parent1_GEBV": Parent1_GEBV_list,
        "Parent2_GEBV": Parent2_GEBV_list
    }).sort_values(by="GEBV", ascending=False)

    F1_out = hybrid_results.head(num_F1)

    parent_GEBV_df = (
        hybrid_results
        .assign(
            Parent1=[x.split(":")[0].strip() for x in hybrid_results["Hybrid"]],
            Parent2=[x.split(":")[1].strip() for x in hybrid_results["Hybrid"]]
        )
        .melt(id_vars=["GEBV"], value_vars=["Parent1", "Parent2"],
              var_name="Type", value_name="Parent")
        .groupby("Parent")
        .agg(GEBV=("GEBV", "mean"))
        .reset_index()
        .sort_values(by="GEBV", ascending=False)
    )

    P_out = parent_GEBV_df.head(num_P)

    # --------------------------------------------------------
    # Optional test-set predictions
    # --------------------------------------------------------

    test_out = None
    test_csv = None

    if config.get("test_geno_file"):

        test_df = pd.read_csv(config["test_geno_file"], index_col=0)
        test_X = scaler.transform(test_df.to_numpy(dtype=np.float32))
        test_GEBV = final_model.predict(test_X, verbose=0).flatten()

        test_out = pd.DataFrame({
            "Test_Genotype": test_df.index,
            "GEBV": test_GEBV
        })

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------

    notify(0.97, "Generating plots...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.patch.set_facecolor("#0f1117")

    top_data = hybrid_results.head(min(20, len(hybrid_results)))
    colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(top_data)))

    axes[0].barh(range(len(top_data)), top_data["GEBV"].values, color=colors)
    axes[0].set_yticks(range(len(top_data)))
    axes[0].set_yticklabels(top_data["Hybrid"].values, fontsize=7, color="white")
    axes[0].set_xlabel("Predicted GEBV", color="white")
    axes[0].set_title("Top Hybrid Combinations", color="white", fontweight="bold")
    axes[0].set_facecolor("#1a1d2e")
    axes[0].tick_params(colors="white")
    for spine in axes[0].spines.values():
        spine.set_edgecolor("#444")

    axes[1].hist(hybrid_results["GEBV"].values, bins=30,
                 color="#7c5cbf", edgecolor="#b388ff", alpha=0.85)
    axes[1].set_xlabel("Predicted GEBV", color="white")
    axes[1].set_ylabel("Count", color="white")
    axes[1].set_title("GEBV Distribution Across All Hybrids", color="white", fontweight="bold")
    axes[1].set_facecolor("#1a1d2e")
    axes[1].tick_params(colors="white")
    for spine in axes[1].spines.values():
        spine.set_edgecolor("#444")

    cv_text = (f"CV Results | MSE: {avg_mse:.4f}  "
               f"RMSE: {avg_rmse:.4f}  R2: {avg_r2:.4f}  Pearson r: {avg_pearson:.4f}")
    fig.suptitle(cv_text, color="#aab4ff", fontsize=9, y=1.01)
    plt.tight_layout()

    plot_path = os.path.join(results_dir, "heterosis_plot.png")
    plt.savefig(plot_path, bbox_inches="tight", facecolor=fig.get_facecolor(), dpi=150)
    plt.close(fig)

    # --------------------------------------------------------
    # CSV exports
    # --------------------------------------------------------

    f1_csv = _to_tmp_csv(F1_out, "Top_F1_Hybrids.csv", results_dir)
    p_csv = _to_tmp_csv(P_out, "Top_Parents.csv", results_dir)
    all_csv = _to_tmp_csv(hybrid_results, "All_Hybrid_Results.csv", results_dir)

    if test_out is not None:
        test_csv = _to_tmp_csv(test_out, "Test_GEBV.csv", results_dir)

    log(f"Done. {len(hybrid_results)} hybrid combinations evaluated.")
    log(f"Top {num_F1} F1 hybrids and top {num_P} parents identified.")

    notify(1.0, "Complete.")

    # --------------------------------------------------------
    # Result payload
    # --------------------------------------------------------

    result = {
        "status": "success",
        "logs": "\n".join(logs),

        "cv_summary": {
            "mse": round(avg_mse, 4),
            "rmse": round(avg_rmse, 4),
            "r2": round(avg_r2, 4),
            "pearson_r": round(avg_pearson, 4)
        },
        "fold_metrics": fold_records,

        "f1_hybrids": _df_to_records(F1_out)["rows"],
        "f1_hybrids_columns": _df_to_records(F1_out)["columns"],

        "top_parents": _df_to_records(P_out)["rows"],
        "top_parents_columns": _df_to_records(P_out)["columns"],

        "test_predictions": _df_to_records(test_out)["rows"] if test_out is not None else [],
        "test_predictions_columns": _df_to_records(test_out)["columns"] if test_out is not None else [],

        "n_hybrids_evaluated": int(len(hybrid_results)),
        "n_parents": int(len(parent_index)),

        "plot_path": plot_path,

        "f1_csv": f1_csv,
        "parents_csv": p_csv,
        "all_csv": all_csv,
        "test_csv": test_csv
    }

    return result
