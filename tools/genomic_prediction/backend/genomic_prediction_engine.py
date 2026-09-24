# -*- coding: utf-8 -*-
"""
genomic_prediction_engine.py — BreedSight Genomic Prediction engine.

Ported from the Gradio app's run_pipeline() and all its supporting
functions (metrics, GRM computation, feature selection, the BreedSight
DNN+RF ensemble model, and every plotting function). No modeling or
statistical logic has been changed — only the orchestration layer
(run_pipeline -> run_breedsight_pipeline) was adapted to take plain file
paths and a log callback instead of Gradio's gr.File objects and
gr.Progress().

The asyncio event-loop patching from the original file (needed to keep a
long-running Gradio server process from warning on shutdown) was removed —
irrelevant outside a Gradio server context.
"""

import os
import sys
import warnings

warnings.filterwarnings("ignore", category=ResourceWarning)
warnings.filterwarnings("ignore", category=UserWarning)

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import io
import traceback
from pathlib import Path

try:
    import tensorflow as tf
    from tensorflow.keras.layers import Dense, BatchNormalization, Dropout, LeakyReLU
    from tensorflow.keras.optimizers import Adam
    from tensorflow.keras import regularizers
    from tensorflow.keras.callbacks import EarlyStopping, ReduceLROnPlateau
    TF_AVAILABLE = True
except ImportError:
    TF_AVAILABLE = False

try:
    from sklearn.preprocessing import StandardScaler
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.feature_selection import RFE
    from sklearn.metrics import mean_squared_error, r2_score
    from sklearn.decomposition import PCA
    from sklearn.model_selection import KFold
    from sklearn.linear_model import Lasso, Ridge
    from sklearn.kernel_ridge import KernelRidge
    from scipy.stats import pearsonr
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False

RANDOM_STATE = 90
np.random.seed(RANDOM_STATE)
if TF_AVAILABLE:
    tf.random.set_seed(RANDOM_STATE)

DEFAULT_PARAMS = {
    'epochs': 300,
    'batch_size': 32,
    'learning_rate': 0.001,
    'l2_reg': 0.01,
    'dropout_rate': 0.3,
    'rf_n_estimators': 200,
    'rf_max_depth': 15,
    'alpha': 0.5,
    'fixed_epochs_final': 300,
    'outer_n_splits': 5,
    'feature_selection': True,
    'heritability': 0.82,
    'use_raw_genotypes': False,
    'use_pca': True,
    'n_components': 50,
    'rfe_n_features': 100,
}

# ═══════════════════════════════════════════════════════════════════════════
#  PALETTE
# ═══════════════════════════════════════════════════════════════════════════
PALETTE = {
    'bg':       '#FFFFFF',
    'surface':  '#F8FAFC',
    'surface2': '#EFF6FF',
    'border':   '#E2E8F0',
    'accent':   '#2563EB',
    'accent2':  '#0EA5E9',
    'green':    '#059669',
    'orange':   '#D97706',
    'red':      '#DC2626',
    'purple':   '#7C3AED',
    'text':     '#1E293B',
    'subtext':  '#64748B',
    'gridline': '#E2E8F0',
}

MODEL_COLORS = {
    'BreedSight': '#2563EB',
    'Lasso':      '#059669',
    'RBLUP':      '#D97706',
    'GBLUP':      '#7C3AED',
}

MODEL_LIGHT = {
    'BreedSight': '#DBEAFE',
    'Lasso':      '#D1FAE5',
    'RBLUP':      '#FEF3C7',
    'GBLUP':      '#EDE9FE',
}

# ═══════════════════════════════════════════════════════════════════════════
#  GLOBAL MATPLOTLIB SETTINGS (applied once at module level)
# ═══════════════════════════════════════════════════════════════════════════
_PLOT_DPI       = 150
MAX_SCATTER_PTS = 3000   # max scatter points per panel

plt.rcParams.update({
    'font.family':         'DejaVu Sans',
    'axes.spines.top':     False,
    'axes.spines.right':   False,
    'axes.grid':           True,
    'grid.color':          '#E2E8F0',
    'grid.linewidth':      0.7,
    'grid.linestyle':      '--',
    'axes.facecolor':      '#FAFBFD',
    'figure.facecolor':    '#FFFFFF',
    'text.color':          '#1E293B',
    'axes.labelcolor':     '#334155',
    'xtick.color':         '#64748B',
    'ytick.color':         '#64748B',
    'axes.edgecolor':      '#CBD5E1',
    'axes.titleweight':    'bold',
    'axes.titlesize':      11,
    'axes.labelsize':      9,
    'legend.frameon':      False,
    'legend.fontsize':     8,
    'figure.dpi':          _PLOT_DPI,
    'savefig.dpi':         _PLOT_DPI,
    'figure.autolayout':   False,   # we control layout explicitly
})


# ═══════════════════════════════════════════════════════════════════════════
#  CORE FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════

def calculate_metrics_with_ci(true_vals, pred_vals, n_bootstraps=300):
    """Compute regression metrics with bootstrap 95% CIs."""
    if len(true_vals) == 0 or len(pred_vals) == 0:
        return {k: np.nan for k in ['mse', 'rmse', 'pearson_r', 'pearson_p', 'r2',
                                     'r2_ci', 'pearson_ci']}
    true_vals = np.asarray(true_vals, dtype=float)
    pred_vals = np.asarray(pred_vals, dtype=float)
    mse  = mean_squared_error(true_vals, pred_vals)
    rmse = np.sqrt(mse)
    if len(np.unique(pred_vals)) == 1 or len(np.unique(true_vals)) == 1:
        corr, p_value, r2 = 0.0, 1.0, 0.0
    else:
        corr, p_value = pearsonr(true_vals, pred_vals)
        r2 = r2_score(true_vals, pred_vals)
    n = len(true_vals)
    r2_boot, corr_boot = [], []
    rng = np.random.default_rng(RANDOM_STATE)
    for _ in range(n_bootstraps):
        idx = rng.integers(0, n, n)
        try:
            if len(np.unique(true_vals[idx])) > 1 and len(np.unique(pred_vals[idx])) > 1:
                r2_boot.append(r2_score(true_vals[idx], pred_vals[idx]))
                corr_boot.append(pearsonr(true_vals[idx], pred_vals[idx])[0])
        except Exception:
            pass
    r2_ci   = list(np.percentile(r2_boot,   [2.5, 97.5])) if len(r2_boot) > 5   else [r2,   r2  ]
    corr_ci = list(np.percentile(corr_boot, [2.5, 97.5])) if len(corr_boot) > 5 else [corr, corr]
    return dict(mse=mse, rmse=rmse, pearson_r=float(corr), pearson_p=float(p_value),
                r2=float(r2), r2_ci=r2_ci, pearson_ci=corr_ci)


def compute_grm_efficient(X_train, X_val=None, ref_features=None, is_train=True,
                          max_markers=10000, use_raw_genotypes=False,
                          use_pca=False, n_components=50):
    if is_train:
        if X_train.shape[1] > max_markers:
            np.random.seed(RANDOM_STATE)
            sel = np.random.choice(X_train.shape[1], max_markers, replace=False)
            X_train = X_train[:, sel]
            if X_val is not None:
                X_val = X_val[:, sel]
        else:
            sel = None

        marker_means = np.mean(X_train, axis=0)
        Xc = X_train - marker_means
        p  = Xc.shape[1]
        G  = np.dot(Xc, Xc.T) / p
        I  = G * G
        md = np.mean(np.diag(I))
        In = I / md if md != 0 else I
        X_tr = np.concatenate([G, In], axis=1)

        X_va = None
        if X_val is not None:
            Vc  = X_val - marker_means
            Gv  = np.dot(Vc, Xc.T) / p
            Iv  = Gv * Gv
            In2 = Iv / md if md != 0 else Iv
            X_va = np.concatenate([Gv, In2], axis=1)

        ref = dict(marker_means=marker_means, mean_diag=md, n_markers=p,
                   selected_markers=sel, X_train_centered_original=Xc)

        if use_pca:
            max_nc = min(X_tr.shape[0], X_tr.shape[1]) - 1
            nc     = max(1, min(int(n_components), max_nc))
            pca    = PCA(n_components=nc)
            X_tr   = pca.fit_transform(X_tr)
            if X_va is not None:
                X_va = pca.transform(X_va)
            ref['pca'] = pca

        if use_raw_genotypes:
            ref['centering_only'] = True
            return Xc, (X_val - marker_means if X_val is not None else None), ref

        return X_tr, X_va, ref

    else:
        if ref_features is None:
            raise ValueError("ref_features must be provided for test phase")
        sel     = ref_features.get('selected_markers')
        if sel is not None:
            X_train = X_train[:, sel]
        Xc_orig = ref_features['X_train_centered_original']
        mm      = ref_features['marker_means']
        p       = ref_features['n_markers']
        md      = ref_features['mean_diag']
        Tc  = X_train - mm
        G   = np.dot(Tc, Xc_orig.T) / p
        I   = G * G
        In  = I / md if md != 0 else I
        X_te = np.concatenate([G, In], axis=1)
        if use_pca and 'pca' in ref_features:
            X_te = ref_features['pca'].transform(X_te)
        return X_te, None, ref_features


def safe_feature_selection(X_train, y_train, X_val, n_features=100):
    if X_train.shape[1] <= n_features:
        return X_train, X_val, None
    rf  = RandomForestRegressor(n_estimators=50, random_state=RANDOM_STATE, n_jobs=-1)
    n_sel = min(n_features, X_train.shape[1])
    sel = RFE(estimator=rf, n_features_to_select=n_sel)
    sel.fit(X_train, y_train)
    Xtr = sel.transform(X_train)
    Xva = sel.transform(X_val) if X_val is not None else None
    return Xtr, Xva, sel


def BreedSight(trainX, trainy, valX=None, valy=None, testX=None,
               epochs=300, batch_size=32, learning_rate=0.001,
               l2_reg=0.01, dropout_rate=0.3,
               rf_n_estimators=200, rf_max_depth=15,
               alpha=0.5, verbose=0, use_early_stopping=True):

    feature_scaler = StandardScaler()
    Xtr = feature_scaler.fit_transform(trainX)
    Xva = feature_scaler.transform(valX)  if valX  is not None else None
    Xte = feature_scaler.transform(testX) if testX is not None else None

    def build_dnn(dim):
        inp = tf.keras.Input(shape=(dim,))
        x   = Dense(256, kernel_initializer='he_normal',
                    kernel_regularizer=regularizers.l2(l2_reg))(inp)
        x   = BatchNormalization()(x); x = Dropout(dropout_rate)(x); x = LeakyReLU(0.1)(x)
        x   = Dense(128, kernel_initializer='he_normal',
                    kernel_regularizer=regularizers.l2(l2_reg))(x)
        x   = BatchNormalization()(x); x = Dropout(dropout_rate)(x); x = LeakyReLU(0.1)(x)
        x   = Dense(64,  kernel_initializer='he_normal',
                    kernel_regularizer=regularizers.l2(l2_reg))(x)
        x   = BatchNormalization()(x); x = Dropout(dropout_rate)(x); x = LeakyReLU(0.1)(x)
        out = Dense(1, activation='linear')(x)
        m   = tf.keras.Model(inp, out)
        m.compile(loss=tf.keras.losses.Huber(delta=1.0),
                  optimizer=Adam(learning_rate=learning_rate, clipvalue=0.5),
                  metrics=['mse'])
        return m

    dnn       = build_dnn(trainX.shape[1])
    callbacks = []
    val_data  = (Xva, valy) if (Xva is not None and valy is not None) else None
    if use_early_stopping and val_data:
        callbacks = [
            EarlyStopping(monitor='val_loss', patience=15,
                          restore_best_weights=True, verbose=0),
            ReduceLROnPlateau(monitor='val_loss', factor=0.5,
                              patience=10, min_lr=1e-7, verbose=0),
        ]
    dnn.fit(Xtr, trainy, epochs=epochs, batch_size=batch_size,
            validation_data=val_data, verbose=verbose,
            callbacks=callbacks, shuffle=True)

    rf = RandomForestRegressor(n_estimators=rf_n_estimators, max_depth=rf_max_depth,
                               random_state=RANDOM_STATE, n_jobs=-1,
                               min_samples_split=5, min_samples_leaf=2)
    rf.fit(trainX, trainy.ravel())

    pd_tr = alpha * dnn.predict(Xtr, verbose=0).flatten() + (1 - alpha) * rf.predict(trainX)
    pd_va = (alpha * dnn.predict(Xva, verbose=0).flatten() + (1 - alpha) * rf.predict(valX)
             if valX is not None else None)
    pd_te = (alpha * dnn.predict(Xte, verbose=0).flatten() + (1 - alpha) * rf.predict(testX)
             if testX is not None else None)

    arts = dict(feature_scaler=feature_scaler, dnn_model=dnn, rf_model=rf)
    return pd_tr, pd_va, pd_te, None, rf, arts


# ═══════════════════════════════════════════════════════════════════════════
#  INTERPRETATION HELPERS
# ═══════════════════════════════════════════════════════════════════════════

def interpret_r2(r2):
    if r2 >= 0.80:
        return ("🟢 Excellent",
                "The model explains ≥80% of phenotypic variance. "
                "This level of accuracy is suitable for genomic selection in breeding programs.")
    elif r2 >= 0.60:
        return ("🟡 Good",
                "The model explains 60–80% of phenotypic variance. "
                "Useful for ranking candidates but may need more markers or larger training sets.")
    elif r2 >= 0.40:
        return ("🟠 Moderate",
                "The model explains 40–60% of variance. "
                "Consider increasing marker density, training population size, or tuning hyperparameters.")
    elif r2 >= 0.20:
        return ("🔴 Poor",
                "The model explains only 20–40% of variance. "
                "Review data quality, check for population stratification, and consider trait heritability.")
    else:
        return ("⛔ Very Poor",
                "The model explains <20% of variance. "
                "Check data integrity, ensure correct file assignment, and verify marker encoding (0/1/2).")


def interpret_pearson(r):
    if r >= 0.90:
        return "Near-perfect prediction accuracy — strong genomic signal captured."
    elif r >= 0.75:
        return "High prediction accuracy — good candidate for genomic selection."
    elif r >= 0.55:
        return "Moderate accuracy — acceptable for low-heritability traits."
    else:
        return "Low accuracy — model may need more data or feature engineering."


def compare_models_text(metrics_df):
    best  = metrics_df.loc[metrics_df['Avg_Val_R2'].idxmax()]
    worst = metrics_df.loc[metrics_df['Avg_Val_R2'].idxmin()]
    lines = []
    lines.append(f"**Best model:** {best['Model']} "
                 f"(Val R² = {best['Avg_Val_R2']:.4f} ± {best['Std_Val_R2']:.4f})")
    lines.append(f"**Weakest model:** {worst['Model']} "
                 f"(Val R² = {worst['Avg_Val_R2']:.4f})")
    diff = best['Avg_Val_R2'] - worst['Avg_Val_R2']
    lines.append(f"**Performance gap:** {diff:.4f} R² units between best and worst")
    bs_row = metrics_df[metrics_df['Model'] == 'BreedSight']
    if not bs_row.empty:
        bs_r2 = bs_row['Avg_Val_R2'].values[0]
        label, msg = interpret_r2(bs_r2)
        lines.append(f"\n**BreedSight Accuracy Rating:** {label}")
        lines.append(f"{msg}")
        lines.append(f"\n**Pearson r interpretation:** "
                     f"{interpret_pearson(bs_row['Avg_Val_r'].values[0])}")
    return "\n\n".join(lines)


def heritability_check(metrics_df, h2):
    lines = ["### Heritability-Based Accuracy Assessment"]
    lines.append(f"Trait heritability (h²) set to **{h2:.2f}**")
    theoretical_max = np.sqrt(h2)
    lines.append(f"Theoretical maximum Pearson r = √h² = **{theoretical_max:.4f}**")
    lines.append("")
    for _, row in metrics_df.iterrows():
        r   = row['Avg_Val_r']
        pct = (r / theoretical_max * 100) if theoretical_max > 0 else 0
        status = "✅" if pct >= 70 else ("⚠️" if pct >= 40 else "❌")
        lines.append(f"{status} **{row['Model']}**: r = {r:.4f} ({pct:.1f}% of theoretical max)")
    lines.append("")
    lines.append("*A model achieving ≥70% of the theoretical maximum is considered "
                 "highly accurate for GS.*")
    return "\n\n".join(lines)


# ═══════════════════════════════════════════════════════════════════════════
#  PLOTTING UTILITIES
# ═══════════════════════════════════════════════════════════════════════════

def _style_ax(ax):
    """Apply consistent styling to a single axes."""
    ax.set_facecolor('#FAFBFD')
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    ax.spines['left'].set_color('#CBD5E1')
    ax.spines['bottom'].set_color('#CBD5E1')
    ax.tick_params(colors='#64748B', length=3)
    ax.grid(True, color='#E2E8F0', linewidth=0.7, linestyle='--', alpha=0.8)


def _downsample(arr1, arr2, max_pts=MAX_SCATTER_PTS):
    """Return downsampled pair of arrays for scatter plots."""
    n = len(arr1)
    if n <= max_pts:
        return arr1, arr2
    rng = np.random.default_rng(RANDOM_STATE)
    idx = rng.choice(n, max_pts, replace=False)
    return arr1[idx], arr2[idx]


def _safe_range(arr):
    """Return (min, max) with a guard against zero-range arrays."""
    lo, hi = float(np.nanmin(arr)), float(np.nanmax(arr))
    if np.isclose(lo, hi):
        lo -= 1.0; hi += 1.0
    return lo, hi


def _save_fig(fig, filepath):
    """Save figure → PNG and close it. Returns filepath on success, None on error."""
    try:
        fig.savefig(
            filepath,
            format='png',
            dpi=_PLOT_DPI,
            bbox_inches='tight',
            facecolor='white',
            edgecolor='none',
            pad_inches=0.2,
        )
        return filepath
    except Exception as e:
        print(f"[_save_fig] Error saving {filepath}: {e}")
        return None
    finally:
        plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────
# PLOT 1 – Regression + Residuals + Distribution
# ─────────────────────────────────────────────────────────────────────────
def make_regression_plot(true_vals, pred_vals, title="", model_name="BreedSight",
                         save_path=None):
    """Three-panel regression analysis: scatter | residuals | distribution."""
    fig = None
    try:
        true_arr = np.asarray(true_vals, dtype=float)
        pred_arr = np.asarray(pred_vals, dtype=float)

        if len(true_arr) < 2:
            raise ValueError("Need at least 2 samples for regression plot.")

        metrics = calculate_metrics_with_ci(true_arr, pred_arr)
        color   = MODEL_COLORS.get(model_name, PALETTE['accent'])
        light   = MODEL_LIGHT.get(model_name,  '#DBEAFE')

        # ── Figure / layout ──────────────────────────────────────────────
        fig = plt.figure(figsize=(15, 5))
        fig.patch.set_facecolor('white')
        gs  = gridspec.GridSpec(1, 3, figure=fig, width_ratios=[2.5, 1.5, 1.5],
                                wspace=0.32, left=0.07, right=0.97,
                                top=0.84, bottom=0.14)
        ax_main = fig.add_subplot(gs[0])
        ax_res  = fig.add_subplot(gs[1])
        ax_dist = fig.add_subplot(gs[2])

        # ── Panel 1: Scatter ─────────────────────────────────────────────
        _style_ax(ax_main)
        lo_t, hi_t = _safe_range(true_arr)
        lo_p, hi_p = _safe_range(pred_arr)
        lo = min(lo_t, lo_p)
        hi = max(hi_t, hi_p)
        pad = (hi - lo) * 0.07

        t_plot, p_plot = _downsample(true_arr, pred_arr)
        ax_main.scatter(t_plot, p_plot, alpha=0.60, s=40,
                        color=color, edgecolors='white', linewidth=0.4, zorder=3)
        lim = [lo - pad, hi + pad]
        ax_main.plot(lim, lim, '--', color='#94A3B8', lw=1.5,
                     label='Perfect fit', zorder=2)

        # Regression line + CI
        z  = np.polyfit(true_arr, pred_arr, 1)
        xs = np.linspace(lo - pad, hi + pad, 300)
        ys = np.poly1d(z)(xs)
        ax_main.plot(xs, ys, '-', color=PALETTE['red'], lw=2,
                     label='Regression fit', zorder=4)
        resid_se = np.std(pred_arr - np.poly1d(z)(true_arr))
        ax_main.fill_between(xs, ys - 1.96 * resid_se, ys + 1.96 * resid_se,
                             color=color, alpha=0.08, zorder=1, label='95% CI band')

        stats_txt = (
            f"R²   = {metrics['r2']:.4f}  [{metrics['r2_ci'][0]:.3f}, {metrics['r2_ci'][1]:.3f}]\n"
            f"r    = {metrics['pearson_r']:.4f}  [{metrics['pearson_ci'][0]:.3f}, {metrics['pearson_ci'][1]:.3f}]\n"
            f"RMSE = {metrics['rmse']:.4f}\n"
            f"n    = {len(true_arr)}"
        )
        ax_main.text(0.04, 0.97, stats_txt, transform=ax_main.transAxes,
                     fontsize=8, va='top',
                     bbox=dict(boxstyle='round,pad=0.45', facecolor=light,
                               edgecolor=color, alpha=0.92, linewidth=1.2),
                     color=PALETTE['text'], family='monospace')
        ax_main.set_xlabel('True Phenotype')
        ax_main.set_ylabel('Predicted Phenotype')
        ax_main.set_title(title or 'Predicted vs True', pad=8)
        ax_main.legend(fontsize=7.5, loc='lower right')
        ax_main.set_xlim(lim); ax_main.set_ylim(lim)

        # ── Panel 2: Residuals ───────────────────────────────────────────
        _style_ax(ax_res)
        residuals = pred_arr - true_arr
        p_ds, r_ds = _downsample(pred_arr, residuals)
        ax_res.scatter(p_ds, r_ds, alpha=0.55, s=28,
                       color=color, edgecolors='white', linewidth=0.3, zorder=3)
        ax_res.axhline(0,               color='#94A3B8', lw=1.5, ls='--', zorder=2)
        ax_res.axhline(residuals.mean(), color=PALETTE['red'], lw=1.2,
                       label=f'Mean={residuals.mean():.3f}', zorder=4)
        lo_p2, hi_p2 = _safe_range(pred_arr)
        ax_res.fill_between(
            [lo_p2, hi_p2],
            residuals.mean() - 1.96 * residuals.std(),
            residuals.mean() + 1.96 * residuals.std(),
            color=color, alpha=0.08)
        ax_res.set_xlabel('Predicted')
        ax_res.set_ylabel('Residual (Pred − True)')
        ax_res.set_title('Residuals', pad=8)
        ax_res.legend(fontsize=7)

        # ── Panel 3: Distribution ────────────────────────────────────────
        _style_ax(ax_dist)
        n_bins = max(8, min(25, len(true_arr) // 5))
        ax_dist.hist(true_arr, bins=n_bins, alpha=0.55, color='#64748B',
                     label='True', density=True, edgecolor='white', linewidth=0.4)
        ax_dist.hist(pred_arr, bins=n_bins, alpha=0.60, color=color,
                     label='Predicted', density=True, edgecolor='white', linewidth=0.4)
        ax_dist.set_xlabel('Phenotype Value')
        ax_dist.set_ylabel('Density')
        ax_dist.set_title('Distribution', pad=8)
        ax_dist.legend(fontsize=7.5)

        fig.suptitle(f'Prediction Analysis — {model_name}', fontsize=12,
                     fontweight='bold', color=PALETTE['text'], y=0.98)

        return _save_fig(fig, save_path) if save_path else (_save_fig(fig, None) or None)

    except Exception as e:
        print(f"[make_regression_plot] Error: {e}\n{traceback.format_exc()}")
        if fig is not None:
            plt.close(fig)
        return None


# ─────────────────────────────────────────────────────────────────────────
# PLOT 2 – Model Comparison Bar
# ─────────────────────────────────────────────────────────────────────────
def make_comparison_bar(metrics_df, save_path=None):
    """Grouped bar chart: Train vs Val R² and Pearson r for each model."""
    fig = None
    try:
        models = metrics_df['Model'].tolist()
        tr_r2  = [float(v) for v in metrics_df['Avg_Train_R2']]
        va_r2  = [float(v) for v in metrics_df['Avg_Val_R2']]
        va_err = [float(v) for v in metrics_df['Std_Val_R2']]
        tr_err = [float(v) for v in metrics_df['Std_Train_R2']]
        tr_r   = [float(v) for v in metrics_df['Avg_Train_r']]
        va_r   = [float(v) for v in metrics_df['Avg_Val_r']]

        x = np.arange(len(models))
        w = 0.36

        fig, axes = plt.subplots(1, 2, figsize=(13, 5))
        fig.patch.set_facecolor('white')
        fig.subplots_adjust(left=0.08, right=0.97, top=0.88, bottom=0.12,
                            wspace=0.28)

        colors = [MODEL_COLORS.get(m, PALETTE['accent']) for m in models]
        lights = [MODEL_LIGHT.get(m, '#DBEAFE') for m in models]

        for ax_idx, (metric_tr, metric_va, err_tr, err_va, ylabel, plot_title) in enumerate([
            (tr_r2, va_r2, tr_err, va_err, 'R²',       'Train vs Validation R²'),
            (tr_r,  va_r,  None,   None,   'Pearson r', 'Train vs Validation Pearson r'),
        ]):
            ax = axes[ax_idx]
            _style_ax(ax)

            kw_tr = dict(capsize=4, error_kw=dict(ecolor='#94A3B8', lw=1.3)) if err_tr else {}
            kw_va = dict(capsize=4, error_kw=dict(ecolor='#475569', lw=1.3)) if err_va else {}

            bars_tr = ax.bar(x - w / 2, metric_tr, w,
                             yerr=err_tr if err_tr else None,
                             color=lights, edgecolor=colors, linewidth=1.4,
                             label='Train', **kw_tr)
            bars_va = ax.bar(x + w / 2, metric_va, w,
                             yerr=err_va if err_va else None,
                             color=colors, edgecolor=colors, linewidth=1.4,
                             label='Validation', alpha=0.88, **kw_va)

            # Value labels above bars
            max_err_tr = max(err_tr) if err_tr else 0
            max_err_va = max(err_va) if err_va else 0
            for bar, val in zip(bars_tr, metric_tr):
                h = bar.get_height()
                ax.text(bar.get_x() + bar.get_width() / 2,
                        h + max_err_tr + 0.012,
                        f'{val:.3f}', ha='center', va='bottom',
                        fontsize=7.5, color=PALETTE['subtext'])
            for bar, val in zip(bars_va, metric_va):
                h = bar.get_height()
                ax.text(bar.get_x() + bar.get_width() / 2,
                        h + max_err_va + 0.012,
                        f'{val:.3f}', ha='center', va='bottom',
                        fontsize=8, color=PALETTE['text'], fontweight='bold')

            ax.set_xticks(x)
            ax.set_xticklabels(models, fontsize=9)
            ax.set_ylabel(ylabel)
            ax.set_title(plot_title, pad=8)
            ax.legend(fontsize=8)

            upper = max(max(metric_tr, default=0), max(metric_va, default=0))
            ax.set_ylim(0, max(upper * 1.30 + 0.10, 0.15))
            if ylabel == 'R²':
                ax.axhline(0.5, color='#94A3B8', lw=0.8, ls=':', alpha=0.7)
                ax.text(len(models) - 0.5, 0.51, 'R²=0.5', fontsize=6.5, color='#94A3B8')

        fig.suptitle('Model Performance Comparison', fontsize=13,
                     fontweight='bold', color=PALETTE['text'])
        return _save_fig(fig, save_path) if save_path else (_save_fig(fig, None) or None)

    except Exception as e:
        print(f"[make_comparison_bar] Error: {e}\n{traceback.format_exc()}")
        if fig is not None:
            plt.close(fig)
        return None


# ─────────────────────────────────────────────────────────────────────────
# PLOT 3 – Fold-by-Fold Lines
# ─────────────────────────────────────────────────────────────────────────
def make_fold_r2_lines(results_df_dict, save_path=None):
    """Fold-by-fold R² and Pearson r with 95% CI shading."""
    fig = None
    try:
        fig, axes = plt.subplots(1, 2, figsize=(13, 5))
        fig.patch.set_facecolor('white')
        fig.subplots_adjust(left=0.08, right=0.97, top=0.88, bottom=0.12, wspace=0.28)

        spec = [
            (axes[0], 'Val_R2',        'Val_R2_CI_lower',        'Val_R2_CI_upper',
             'Validation R²',           'R² per Fold'),
            (axes[1], 'Val_Pearson_r', 'Val_Pearson_CI_lower', 'Val_Pearson_CI_upper',
             'Validation Pearson r', 'Pearson r per Fold'),
        ]

        for ax, val_col, lo_col, hi_col, ylabel, ptitle in spec:
            _style_ax(ax)
            for model, df in results_df_dict.items():
                c     = MODEL_COLORS.get(model, PALETTE['accent'])
                folds = df['Fold']
                vals  = df[val_col].astype(float)
                lo    = df[lo_col].astype(float).fillna(vals)
                hi    = df[hi_col].astype(float).fillna(vals)

                ax.plot(folds, vals, marker='o', lw=2.0, color=c,
                        label=model, alpha=0.95, markersize=6,
                        markeredgecolor='white', markeredgewidth=1.3, zorder=3)
                ax.fill_between(folds, lo, hi, color=c, alpha=0.10, zorder=1)

                # Annotate last point
                ax.annotate(f'{vals.iloc[-1]:.3f}',
                            xy=(folds.iloc[-1], vals.iloc[-1]),
                            xytext=(5, 2), textcoords='offset points',
                            fontsize=7, color=c, va='bottom')

            ax.set_xlabel('Fold')
            ax.set_ylabel(ylabel)
            ax.set_title(ptitle, pad=8)
            ax.legend(fontsize=7.5)
            ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

        fig.suptitle('Cross-Validation Performance per Fold (95% CI shaded)',
                     fontsize=12, fontweight='bold', color=PALETTE['text'])
        return _save_fig(fig, save_path) if save_path else (_save_fig(fig, None) or None)

    except Exception as e:
        print(f"[make_fold_r2_lines] Error: {e}\n{traceback.format_exc()}")
        if fig is not None:
            plt.close(fig)
        return None


# ─────────────────────────────────────────────────────────────────────────
# PLOT 4 – Test Prediction Scatter Grid
# ─────────────────────────────────────────────────────────────────────────
def make_prediction_scatter_grid(test_pred_df, has_true, save_path=None):
    """2×N scatter grid for test predictions across all models."""
    fig = None
    try:
        if not has_true:
            return None

        models = [c.replace('Predicted_Phenotype_', '')
                  for c in test_pred_df.columns
                  if c.startswith('Predicted_Phenotype_')]
        n = len(models)
        if n == 0:
            return None

        nc      = min(2, n)
        nr      = (n + nc - 1) // nc
        cell_w  = 6.0
        cell_h  = 5.5
        fig, axes = plt.subplots(
            nr, nc,
            figsize=(cell_w * nc, cell_h * nr),
            squeeze=False,
        )
        fig.patch.set_facecolor('white')
        fig.subplots_adjust(
            left=0.08, right=0.97,
            top=1.0 - 0.06 / nr,   # leave room for suptitle
            bottom=0.10,
            wspace=0.32, hspace=0.42,
        )

        true_arr  = test_pred_df['True_Phenotype'].values.astype(float)
        axes_flat = axes.flatten()

        for i, model in enumerate(models):
            ax    = axes_flat[i]
            _style_ax(ax)
            pred  = test_pred_df[f'Predicted_Phenotype_{model}'].values.astype(float)
            color = MODEL_COLORS.get(model, PALETTE['accent'])
            light = MODEL_LIGHT.get(model, '#DBEAFE')

            lo = min(float(np.nanmin(true_arr)), float(np.nanmin(pred)))
            hi = max(float(np.nanmax(true_arr)), float(np.nanmax(pred)))
            if np.isclose(lo, hi):
                lo -= 1.0; hi += 1.0
            pad = (hi - lo) * 0.07
            lim = [lo - pad, hi + pad]

            t_ds, p_ds = _downsample(true_arr, pred)
            ax.scatter(t_ds, p_ds, alpha=0.55, s=36, color=color,
                       edgecolors='white', linewidth=0.4, zorder=3)
            ax.plot(lim, lim, '--', color='#94A3B8', lw=1.5, zorder=2)

            if len(np.unique(true_arr)) > 1 and len(np.unique(pred)) > 1:
                z  = np.polyfit(true_arr, pred, 1)
                xs = np.linspace(lim[0], lim[1], 200)
                ax.plot(xs, np.poly1d(z)(xs), '-', color=color, lw=2, zorder=4)

            r2_val, _ = (r2_score(true_arr, pred), None)
            r_val, _  = pearsonr(true_arr, pred)
            rmse_val  = np.sqrt(mean_squared_error(true_arr, pred))
            label, _  = interpret_r2(r2_val)

            ax.text(0.04, 0.97,
                    f"R² = {r2_val:.4f}\nr  = {r_val:.4f}\nRMSE = {rmse_val:.4f}\n{label}",
                    transform=ax.transAxes, fontsize=8, va='top',
                    bbox=dict(boxstyle='round,pad=0.40', facecolor=light,
                              edgecolor=color, alpha=0.92, linewidth=1.2),
                    color=PALETTE['text'], family='monospace')

            ax.set_title(model, fontsize=10, fontweight='bold', color=color, pad=7)
            ax.set_xlabel('True Phenotype')
            ax.set_ylabel('Predicted Phenotype')
            ax.set_xlim(lim); ax.set_ylim(lim)

        # Hide unused subplots
        for j in range(n, len(axes_flat)):
            axes_flat[j].set_visible(False)

        fig.suptitle('Test Set Predictions — All Models', fontsize=13,
                     fontweight='bold', color=PALETTE['text'])
        return _save_fig(fig, save_path) if save_path else (_save_fig(fig, None) or None)

    except Exception as e:
        print(f"[make_prediction_scatter_grid] Error: {e}\n{traceback.format_exc()}")
        if fig is not None:
            plt.close(fig)
        return None


# ─────────────────────────────────────────────────────────────────────────
# PLOT 5 – RMSE Heatmap
# ─────────────────────────────────────────────────────────────────────────
def make_rmse_heatmap(results_df_dict, save_path=None):
    """Heatmap: models × folds for RMSE and R²."""
    fig = None
    try:
        rows = []
        for model, df in results_df_dict.items():
            for _, row in df.iterrows():
                rows.append({'Model':    model,
                             'Fold':     int(row['Fold']),
                             'Val_RMSE': float(row['Val_RMSE']),
                             'Val_R2':   float(row['Val_R2'])})
        df_all = pd.DataFrame(rows)

        fig, axes = plt.subplots(1, 2, figsize=(13, max(4.0, len(results_df_dict) * 1.0 + 2.0)))
        fig.patch.set_facecolor('white')
        fig.subplots_adjust(left=0.14, right=0.97, top=0.88, bottom=0.14, wspace=0.36)

        for ax, metric, cmap_name, fmt, ptitle in [
            (axes[0], 'Val_RMSE', 'YlOrRd', '.4f', 'Validation RMSE (lower = better)'),
            (axes[1], 'Val_R2',   'Blues',  '.4f', 'Validation R² (higher = better)'),
        ]:
            pivot = df_all.pivot(index='Model', columns='Fold', values=metric)
            cmap  = matplotlib.colormaps[cmap_name]
            im    = ax.imshow(pivot.values, cmap=cmap, aspect='auto')
            fig.colorbar(im, ax=ax, shrink=0.80, pad=0.03)

            ax.set_xticks(range(pivot.shape[1]))
            ax.set_xticklabels([f'F{c}' for c in pivot.columns], fontsize=8)
            ax.set_yticks(range(pivot.shape[0]))
            ax.set_yticklabels(pivot.index, fontsize=8)
            ax.set_title(ptitle, pad=8, fontsize=10, fontweight='bold')
            for sp in ax.spines.values():
                sp.set_visible(False)
            ax.tick_params(colors='#64748B', length=0)

            vmin, vmax = pivot.values.min(), pivot.values.max()
            vrange = vmax - vmin if not np.isclose(vmin, vmax) else 1.0
            for (r_i, c_i), val in np.ndenumerate(pivot.values):
                norm_v    = (val - vmin) / vrange
                txt_color = 'white' if norm_v > 0.60 else PALETTE['text']
                ax.text(c_i, r_i, f'{val:{fmt}}',
                        ha='center', va='center', fontsize=8, color=txt_color)

        fig.suptitle('Performance Heatmap Across Folds', fontsize=12,
                     fontweight='bold', color=PALETTE['text'])
        return _save_fig(fig, save_path) if save_path else (_save_fig(fig, None) or None)

    except Exception as e:
        print(f"[make_rmse_heatmap] Error: {e}\n{traceback.format_exc()}")
        if fig is not None:
            plt.close(fig)
        return None


# ─────────────────────────────────────────────────────────────────────────
# PLOT 6 – Prediction Distribution (Violin + Box)
# ─────────────────────────────────────────────────────────────────────────
def make_prediction_distribution(test_pred_df, has_true, save_path=None):
    """Violin + box plot of predicted distributions per model."""
    fig = None
    try:
        models = [c.replace('Predicted_Phenotype_', '')
                  for c in test_pred_df.columns
                  if c.startswith('Predicted_Phenotype_')]

        fig, ax = plt.subplots(figsize=(11, 5))
        fig.patch.set_facecolor('white')
        fig.subplots_adjust(left=0.10, right=0.97, top=0.88, bottom=0.14)
        _style_ax(ax)

        data_list, labels, colors_v = [], [], []
        if has_true and 'True_Phenotype' in test_pred_df.columns:
            data_list.append(test_pred_df['True_Phenotype'].values.astype(float))
            labels.append('True')
            colors_v.append('#64748B')
        for m in models:
            col = f'Predicted_Phenotype_{m}'
            if col in test_pred_df.columns:
                data_list.append(test_pred_df[col].values.astype(float))
                labels.append(m)
                colors_v.append(MODEL_COLORS.get(m, PALETTE['accent']))

        if not data_list:
            plt.close(fig)
            return None

        positions = list(range(len(data_list)))

        # Violin only for arrays with non-zero std
        viol_data, viol_pos, viol_cols = [], [], []
        for idx, d in enumerate(data_list):
            if np.std(d) > 1e-9 and len(d) >= 4:
                viol_data.append(d)
                viol_pos.append(positions[idx])
                viol_cols.append(colors_v[idx])

        if viol_data:
            parts = ax.violinplot(viol_data, positions=viol_pos,
                                  showmedians=True, showextrema=True,
                                  widths=0.55)
            for pc, c in zip(parts['bodies'], viol_cols):
                pc.set_facecolor(c); pc.set_alpha(0.30); pc.set_edgecolor(c)
            parts['cmedians'].set_colors(viol_cols)
            for key in ('cbars', 'cmaxes', 'cmins'):
                parts[key].set_colors('#CBD5E1')

        # Box plots on top
        bp = ax.boxplot(data_list, positions=positions,
                        widths=0.12, patch_artist=True,
                        boxprops=dict(facecolor='white', linewidth=1.2),
                        medianprops=dict(linewidth=2),
                        whiskerprops=dict(linewidth=1.0, color='#94A3B8'),
                        capprops=dict(linewidth=1.0, color='#94A3B8'),
                        flierprops=dict(marker='o', ms=3, alpha=0.35, linestyle='none'))
        for patch, c in zip(bp['boxes'],   colors_v): patch.set_edgecolor(c)
        for med,   c in zip(bp['medians'], colors_v): med.set_color(c)

        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel('Phenotype Value')
        ax.set_title('Distribution of True and Predicted Phenotypes', pad=8)

        patches = [mpatches.Patch(facecolor=c, label=l, alpha=0.70)
                   for c, l in zip(colors_v, labels)]
        ax.legend(handles=patches, loc='upper right', fontsize=8)

        fig.suptitle('Phenotype Distribution Analysis', fontsize=12,
                     fontweight='bold', color=PALETTE['text'])
        return _save_fig(fig, save_path) if save_path else (_save_fig(fig, None) or None)

    except Exception as e:
        print(f"[make_prediction_distribution] Error: {e}\n{traceback.format_exc()}")
        if fig is not None:
            plt.close(fig)
        return None


# ═══════════════════════════════════════════════════════════════════════════
#  MAIN PIPELINE
# ═══════════════════════════════════════════════════════════════════════════



# ═══════════════════════════════════════════════════════════════════════════
#  MAIN PIPELINE (adapted from the Gradio app's run_pipeline())
# ═══════════════════════════════════════════════════════════════════════════
# Changes from the original:
#   - Takes plain file path strings instead of gr.File objects (.name attr)
#   - Takes an `outdir` Path instead of a hardcoded "/tmp/breedSight_results"
#   - Takes an optional `on_log` callback instead of appending to a Gradio
#     log textbox; if provided, it's called with each new log line so the
#     desktop UI can stream progress in real time
#   - Returns a JSON-serializable dict instead of an 18-tuple of Gradio
#     component values
# No modeling/statistical logic has been changed.

def run_breedsight_pipeline(
    train_pheno_file, train_geno_file,
    test_pheno_file, test_geno_file,
    outer_n_splits, epochs, batch_size, learning_rate,
    l2_reg, dropout_rate, rf_n_estimators, rf_max_depth,
    alpha, heritability, use_pca, n_components,
    feature_selection, rfe_n_features, use_raw_genotypes,
    outdir,
    on_log=None,
):
    logs = []
    def log(msg):
        logs.append(msg)
        if on_log:
            try:
                on_log(msg)
            except Exception:
                pass
        return "\n".join(logs)

    if not SKLEARN_AVAILABLE:
        return {"status": "error", "message": "scikit-learn is not installed.", "log": log("scikit-learn not installed")}
    if not TF_AVAILABLE:
        return {"status": "error", "message": "TensorFlow is not installed.", "log": log("TensorFlow not installed")}

    try:
        log("Loading data...")

        training_data     = pd.read_csv(train_pheno_file)
        training_additive = pd.read_csv(train_geno_file)
        testing_data       = pd.read_csv(test_pheno_file)
        testing_additive   = pd.read_csv(test_geno_file)

        log(f"  Training phenotypes : {training_data.shape}")
        log(f"  Training genotypes  : {training_additive.shape}")
        log(f"  Testing  phenotypes : {testing_data.shape}")
        log(f"  Testing  genotypes  : {testing_additive.shape}")

        if training_data.shape[1] > training_additive.shape[1]:
            log("WARNING: Training phenotypic file has MORE columns than genotypic file.")
            log("  -> Did you swap the phenotypic and genotypic uploads?")

        if 'phenotypes' not in training_data.columns:
            return {"status": "error", "message": "'phenotypes' column missing from training file.", "log": "\n".join(logs)}

        if training_data.shape[0] != training_additive.shape[0]:
            return {"status": "error", "message": f"Row count mismatch: training pheno={training_data.shape[0]}, geno={training_additive.shape[0]}", "log": "\n".join(logs)}

        if testing_data.shape[0] != testing_additive.shape[0]:
            return {"status": "error", "message": f"Row count mismatch: testing pheno={testing_data.shape[0]}, geno={testing_additive.shape[0]}", "log": "\n".join(logs)}

        log("Validating data...")
        train_add_raw = training_additive.iloc[:, 1:].values.astype(float)
        test_add_raw  = testing_additive.iloc[:,  1:].values.astype(float)
        y_train       = training_data['phenotypes'].values.astype(float)

        has_test_pheno = 'phenotypes' in testing_data.columns
        y_test   = testing_data['phenotypes'].values.astype(float) if has_test_pheno else None
        test_ids = testing_data.iloc[:, 0].values

        if np.isnan(train_add_raw).any() or np.isnan(y_train).any():
            return {"status": "error", "message": "Missing values in training data.", "log": "\n".join(logs)}
        if np.isnan(test_add_raw).any():
            return {"status": "error", "message": "Missing values in testing genotypes.", "log": "\n".join(logs)}
        log("  No missing values detected")

        log("Data Summary:")
        log(f"  Training samples   : {len(y_train)}")
        log(f"  Testing  samples   : {len(test_ids)}")
        log(f"  Markers (SNPs)     : {train_add_raw.shape[1]}")
        log(f"  Pheno mean +/- std : {y_train.mean():.4f} +/- {y_train.std():.4f}")
        log(f"  Pheno range        : [{y_train.min():.4f}, {y_train.max():.4f}]")

        log(f"Starting {outer_n_splits}-Fold Cross-Validation...")

        results_dict  = {m: [] for m in ['BreedSight', 'Lasso', 'RBLUP', 'GBLUP']}
        val_pred_list = {m: [] for m in results_dict}
        kf = KFold(n_splits=int(outer_n_splits), shuffle=True, random_state=RANDOM_STATE)
        rf_depth = None if int(rf_max_depth) <= 0 else int(rf_max_depth)

        for fold, (tr_idx, va_idx) in enumerate(kf.split(train_add_raw), 1):
            log(f"-- Fold {fold}/{outer_n_splits} --")

            Xtr_raw = train_add_raw[tr_idx]
            Xva_raw = train_add_raw[va_idx]
            ytr     = y_train[tr_idx]
            yva     = y_train[va_idx]

            Xtr_g, Xva_g, ref = compute_grm_efficient(
                Xtr_raw, Xva_raw, is_train=True,
                use_raw_genotypes=bool(use_raw_genotypes),
                use_pca=bool(use_pca), n_components=int(n_components))

            if feature_selection:
                Xtr_f, Xva_f, _ = safe_feature_selection(
                    Xtr_g, ytr, Xva_g, int(rfe_n_features))
            else:
                Xtr_f, Xva_f = Xtr_g, Xva_g

            # BreedSight
            log("    BreedSight...")
            pd_tr, pd_va, _, _, _, _ = BreedSight(
                Xtr_f, ytr, Xva_f, yva, None,
                epochs=int(epochs), batch_size=int(batch_size),
                learning_rate=float(learning_rate), l2_reg=float(l2_reg),
                dropout_rate=float(dropout_rate), rf_n_estimators=int(rf_n_estimators),
                rf_max_depth=rf_depth, alpha=float(alpha),
                verbose=0, use_early_stopping=True)

            for model_name, (pdt, pdv) in [('BreedSight', (pd_tr, pd_va))]:
                tm = calculate_metrics_with_ci(ytr, pdt)
                vm = calculate_metrics_with_ci(yva, pdv)
                results_dict[model_name].append({
                    'Fold': fold,
                    'Train_R2': tm['r2'], 'Train_R2_CI_lower': tm['r2_ci'][0],
                    'Train_R2_CI_upper': tm['r2_ci'][1],
                    'Val_R2': vm['r2'],   'Val_R2_CI_lower': vm['r2_ci'][0],
                    'Val_R2_CI_upper': vm['r2_ci'][1],
                    'Train_MSE': tm['mse'],   'Val_MSE': vm['mse'],
                    'Train_RMSE': tm['rmse'], 'Val_RMSE': vm['rmse'],
                    'Train_Pearson_r': tm['pearson_r'],
                    'Train_Pearson_CI_lower': tm['pearson_ci'][0],
                    'Train_Pearson_CI_upper': tm['pearson_ci'][1],
                    'Val_Pearson_r': vm['pearson_r'],
                    'Val_Pearson_CI_lower': vm['pearson_ci'][0],
                    'Val_Pearson_CI_upper': vm['pearson_ci'][1],
                })
                val_pred_list[model_name].append(
                    pd.DataFrame({'Sample_ID': training_data.iloc[va_idx, 0].values,
                                  f'Predicted_{model_name}': pdv,
                                  'True_Phenotype': yva, 'Fold': fold}))
                log(f"      BreedSight  Train R2={tm['r2']:.4f}  Val R2={vm['r2']:.4f}  "
                    f"Val r={vm['pearson_r']:.4f}")

            # Traditional models
            for mname, mdl in [
                ('Lasso', Lasso(alpha=0.1, random_state=RANDOM_STATE, max_iter=10000)),
                ('RBLUP', Ridge(alpha=1.0, random_state=RANDOM_STATE)),
                ('GBLUP', KernelRidge(alpha=1.0)),
            ]:
                log(f"    {mname}...")
                mdl.fit(Xtr_f, ytr)
                pdt2 = mdl.predict(Xtr_f)
                pdv2 = mdl.predict(Xva_f)
                tm2  = calculate_metrics_with_ci(ytr, pdt2)
                vm2  = calculate_metrics_with_ci(yva, pdv2)
                results_dict[mname].append({
                    'Fold': fold,
                    'Train_R2': tm2['r2'], 'Train_R2_CI_lower': tm2['r2_ci'][0],
                    'Train_R2_CI_upper': tm2['r2_ci'][1],
                    'Val_R2': vm2['r2'],   'Val_R2_CI_lower': vm2['r2_ci'][0],
                    'Val_R2_CI_upper': vm2['r2_ci'][1],
                    'Train_MSE': tm2['mse'],   'Val_MSE': vm2['mse'],
                    'Train_RMSE': tm2['rmse'], 'Val_RMSE': vm2['rmse'],
                    'Train_Pearson_r': tm2['pearson_r'],
                    'Train_Pearson_CI_lower': tm2['pearson_ci'][0],
                    'Train_Pearson_CI_upper': tm2['pearson_ci'][1],
                    'Val_Pearson_r': vm2['pearson_r'],
                    'Val_Pearson_CI_lower': vm2['pearson_ci'][0],
                    'Val_Pearson_CI_upper': vm2['pearson_ci'][1],
                })
                val_pred_list[mname].append(
                    pd.DataFrame({'Sample_ID': training_data.iloc[va_idx, 0].values,
                                  f'Predicted_{mname}': pdv2,
                                  'True_Phenotype': yva, 'Fold': fold}))
                log(f"      {mname:12} Train R2={tm2['r2']:.4f}  Val R2={vm2['r2']:.4f}  "
                    f"Val r={vm2['pearson_r']:.4f}")

        # ── Final model on full training set ─────────────────────────────
        log("Training final models on full training set...")

        Xtr_g_final, _, ref_final = compute_grm_efficient(
            train_add_raw, is_train=True,
            use_raw_genotypes=bool(use_raw_genotypes),
            use_pca=bool(use_pca), n_components=int(n_components))
        Xte_g_final, _, _ = compute_grm_efficient(
            test_add_raw, ref_features=ref_final, is_train=False,
            use_raw_genotypes=bool(use_raw_genotypes),
            use_pca=bool(use_pca), n_components=int(n_components))

        if feature_selection:
            Xtr_f2, Xte_f2, _ = safe_feature_selection(
                Xtr_g_final, y_train, Xte_g_final, int(rfe_n_features))
        else:
            Xtr_f2, Xte_f2 = Xtr_g_final, Xte_g_final

        final_preds = {}
        _, _, pd_te_bs, _, _, _ = BreedSight(
            Xtr_f2, y_train, None, None, Xte_f2,
            epochs=int(epochs), batch_size=int(batch_size),
            learning_rate=float(learning_rate), l2_reg=float(l2_reg),
            dropout_rate=float(dropout_rate), rf_n_estimators=int(rf_n_estimators),
            rf_max_depth=rf_depth, alpha=float(alpha),
            verbose=0, use_early_stopping=False)
        final_preds['BreedSight'] = pd_te_bs

        for mname, mdl in [
            ('Lasso', Lasso(alpha=0.1, random_state=RANDOM_STATE, max_iter=10000)),
            ('RBLUP', Ridge(alpha=1.0, random_state=RANDOM_STATE)),
            ('GBLUP', KernelRidge(alpha=1.0)),
        ]:
            mdl.fit(Xtr_f2, y_train)
            final_preds[mname] = mdl.predict(Xte_f2)
        log("  Final models trained")

        # ── Compile results ──────────────────────────────────────────────
        results_df_dict = {m: pd.DataFrame(v) for m, v in results_dict.items()}

        metrics_summary = {}
        for mname, df in results_df_dict.items():
            metrics_summary[mname] = {
                'Model':          mname,
                'Avg_Train_R2':   float(df['Train_R2'].mean()),
                'Std_Train_R2':   float(df['Train_R2'].std()),
                'Avg_Val_R2':     float(df['Val_R2'].mean()),
                'Std_Val_R2':     float(df['Val_R2'].std()),
                'Avg_Train_RMSE': float(df['Train_RMSE'].mean()),
                'Avg_Val_RMSE':   float(df['Val_RMSE'].mean()),
                'Avg_Train_r':    float(df['Train_Pearson_r'].mean()),
                'Avg_Val_r':      float(df['Val_Pearson_r'].mean()),
                'N_Folds':        len(df),
            }
        metrics_df = pd.DataFrame(list(metrics_summary.values()))

        test_pred_df = pd.DataFrame({'Sample_ID': test_ids})
        for mname, preds in final_preds.items():
            test_pred_df[f'Predicted_Phenotype_{mname}'] = preds
        if has_test_pheno:
            test_pred_df['True_Phenotype'] = y_test

        # ── Save CSVs ────────────────────────────────────────────────────
        outdir.mkdir(parents=True, exist_ok=True)
        metrics_csv  = outdir / "model_metrics.csv"
        test_csv     = outdir / "test_predictions.csv"
        val_all      = pd.concat(
            [pd.concat(v, ignore_index=True) for v in val_pred_list.values()],
            ignore_index=True)
        val_csv = outdir / "validation_predictions.csv"
        metrics_df.to_csv(metrics_csv, index=False)
        test_pred_df.to_csv(test_csv,  index=False)
        val_all.to_csv(val_csv,        index=False)

        # ── Generate plots ───────────────────────────────────────────────
        log("Generating plots...")

        bs_df   = results_df_dict['BreedSight']
        best_f  = int(bs_df.loc[bs_df['Val_R2'].idxmax(), 'Fold'])
        bs_fold_frames = val_pred_list['BreedSight']
        bs_fold_val    = next(
            (d for d in bs_fold_frames if d['Fold'].iloc[0] == best_f), None)

        reg_path  = str(outdir / "plot_regression.png")
        cmp_path  = str(outdir / "plot_comparison.png")
        fold_path = str(outdir / "plot_folds.png")
        heat_path = str(outdir / "plot_heatmap.png")
        dist_path = str(outdir / "plot_distribution.png")
        test_path_plot = str(outdir / "plot_test_grid.png") if has_test_pheno else None

        if bs_fold_val is not None:
            r = make_regression_plot(
                bs_fold_val['True_Phenotype'].values,
                bs_fold_val['Predicted_BreedSight'].values,
                f'Best Validation Fold {best_f}', 'BreedSight',
                save_path=reg_path)
            reg_path = r if r else None
            log(f"  {'OK' if r else 'FAILED'}: Regression plot")
        else:
            reg_path = None

        r = make_comparison_bar(metrics_df, save_path=cmp_path)
        cmp_path = r if r else None
        log(f"  {'OK' if r else 'FAILED'}: Comparison bar chart")

        r = make_fold_r2_lines(results_df_dict, save_path=fold_path)
        fold_path = r if r else None
        log(f"  {'OK' if r else 'FAILED'}: Fold R2 lines")

        r = make_rmse_heatmap(results_df_dict, save_path=heat_path)
        heat_path = r if r else None
        log(f"  {'OK' if r else 'FAILED'}: Performance heatmap")

        r = make_prediction_distribution(test_pred_df, has_test_pheno, save_path=dist_path)
        dist_path = r if r else None
        log(f"  {'OK' if r else 'FAILED'}: Distribution plot")

        if has_test_pheno:
            r = make_prediction_scatter_grid(test_pred_df, has_test_pheno,
                                             save_path=test_path_plot)
            test_path_plot = r if r else None
            log(f"  {'OK' if r else 'FAILED'}: Test grid scatter")
        else:
            test_path_plot = None

        # ── Interpretation texts ─────────────────────────────────────────
        interp_model = compare_models_text(metrics_df)
        interp_h2    = heritability_check(metrics_df, float(heritability))

        best_model  = metrics_df.loc[metrics_df['Avg_Val_R2'].idxmax(), 'Model']
        best_val_r2 = metrics_df.loc[metrics_df['Avg_Val_R2'].idxmax(), 'Avg_Val_R2']
        log(f"Best model: {best_model}  (Val R2 = {best_val_r2:.4f})")
        log("Pipeline complete!")

        return {
            "status": "success",
            "log": "\n".join(logs),
            "metrics_table": metrics_df.to_dict(orient="records"),
            "test_predictions_table": test_pred_df.head(100).to_dict(orient="records"),
            "interpretation_model": interp_model,
            "interpretation_heritability": interp_h2,
            "best_model": str(best_model),
            "best_val_r2": float(best_val_r2),
            "plots": {
                "regression": reg_path,
                "comparison": cmp_path,
                "folds": fold_path,
                "test_grid": test_path_plot,
                "heatmap": heat_path,
                "distribution": dist_path,
            },
            "downloads": {
                "metrics_csv": str(metrics_csv),
                "test_csv": str(test_csv),
                "val_csv": str(val_csv),
            },
        }

    except Exception as e:
        tb = traceback.format_exc()
        log(f"Error: {e}\n{tb}")
        return {"status": "error", "message": str(e), "log": "\n".join(logs)}
