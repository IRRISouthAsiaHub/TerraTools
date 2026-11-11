"""
Linear and mixed-model association scans, and the null-model significance threshold.

Model            y = X0 b0 + g b + u + e        u ~ N(0, s2 K)   e ~ N(0, s2 d I)       (V = s2 (K + d I))
K = U diag(w) U' is decomposed once. In the rotated space y* = U'y, X* = U'X the covariance is diagonal, s2 (w_i + d), so
every quantity below is a weighted least-squares problem with weights 1/(w_i + d).

REML / ML      -2 l = (n - p) log s2 + sum log(w_i + d) + log|X' H X| (REML), with H = diag(1/(w_i + d)); the full constants are
               kept so the value can be compared with a dense implementation. d is searched on a log grid and refined by a
               bounded one-dimensional search; a solution on the edge of the grid (no polygenic variance / no residual
               variance) is reported, not hidden.

Models (what is estimated once, what is re-estimated for every marker):
  ols    no random effect (d -> infinity). Covariates optional. t test, n - p0 - 1 degrees of freedom.
  mlm    variance ratio d from the null model (REML), fixed for all markers; s2 re-estimated per marker; t test.
  emmax  variance ratio AND s2 from the null model, both fixed for all markers; the Wald statistic is then normal (chi-square, 1 df).
  gemma  d re-estimated per marker by REML; Wald test (t, n - p0 - 1 df) or likelihood-ratio test (maximum likelihood fits of the
         null and alternative models, chi-square 1 df). This is a per-marker REML scan written for this program; it is only
         called "equivalent" to an external program if a comparison run says so (see gemma_external.py).
Leave-one-chromosome-out: the null model and the eigendecomposition are redone with the K that excludes the marker's chromosome.

Significance by simulation under the fitted null model (null_threshold): replicate phenotypes keep the covariate effects and the
covariance structure and carry no marker effect.
  residual_permutation   the whitened null residuals v = H^(1/2)(y* - X* b) are exchangeable under the null model; permuting them
                         and transforming back gives y_r = X b + U diag(sqrt(w + d)) v_r. Covariates are kept; structure is kept
                         because the permuting happens after whitening.
  parametric_bootstrap   v_r is drawn from N(0, s2 I) instead.
For every replicate all markers are tested with the same statistic as the real scan (variance components fixed at the null
estimates), and the maximum |statistic| is recorded. The threshold at level a is the (1 - a) quantile of these maxima, reported
with a distribution-free interval. It controls the family-wise error rate for THAT model and THAT data set; thresholds must not be
carried across models. It is an approximation for the per-marker REML model (the variance ratio is held at the null estimate in the
replicates) and is not available for iterative multi-locus procedures.
"""

from __future__ import annotations

import math

import numpy as np
from scipy import stats
from scipy.optimize import minimize_scalar

from tools.gwas.backend.population import eig_psd

LOG_LO, LOG_HI = -6.0, 6.0
GRID_POINTS = 49
METHODS = ("ols", "mlm", "emmax", "gemma")
STAT_NAMES = {"ols": "t", "mlm": "t", "emmax": "z", "gemma": "t"}


class ModelError(ValueError):
    """The requested model cannot be fitted to these data."""


# ============================================================================================ design

def design_matrix(n, covariates=None, names=None):
    """[1, covariates] with dependent columns removed.

    Returns (X0, used_names, dropped) where dropped lists (name, reason). A column is dropped when it is (numerically) a linear
    combination of the columns before it, so a covariate that repeats the intercept or an earlier covariate never enters twice.
    """
    cols, used, dropped = [np.ones(n)], ["Intercept"], []
    basis = [np.ones(n) / math.sqrt(n)]
    if covariates is not None:
        C = np.asarray(covariates, dtype=float)
        C = C.reshape(n, -1) if C.size else C.reshape(n, 0)
        names = list(names) if names is not None else [f"cov{i + 1}" for i in range(C.shape[1])]
        for j in range(C.shape[1]):
            c = C[:, j]
            if not np.all(np.isfinite(c)):
                raise ModelError(f"Covariate '{names[j]}' has missing values; remove those samples or fill them in first.")
            r = c.copy()
            for _ in range(2):
                for b in basis:
                    r = r - b * (b @ r)
            nc = float(np.linalg.norm(c))
            if nc == 0 or float(np.linalg.norm(r)) < 1e-6 * nc:
                dropped.append((names[j], "constant or a linear combination of earlier columns"))
                continue
            basis.append(r / np.linalg.norm(r))
            cols.append(c)
            used.append(names[j])
    return np.column_stack(cols), used, dropped


# ============================================================================================ null model

def _null_eval(ys, X0s, w, delta, mode):
    n, p = X0s.shape
    h = 1.0 / (w + delta)
    Xh = X0s * h[:, None]
    A = Xh.T @ X0s
    beta = np.linalg.solve(A, Xh.T @ ys)
    res = ys - X0s @ beta
    rss = float(res @ (h * res))
    sumlog = float(np.log(w + delta).sum())
    logdetA = np.linalg.slogdet(A)[1]
    if mode == "REML":
        s2 = rss / (n - p)
        ll = -0.5 * ((n - p) * math.log(2 * math.pi * s2) + sumlog + logdetA + (n - p)) + 0.5 * np.linalg.slogdet(X0s.T @ X0s)[1]
    else:
        s2 = rss / n
        ll = -0.5 * (n * math.log(2 * math.pi * s2) + sumlog + n)
    return ll, beta, s2, rss


def null_fit(ys, X0s, w, mode="REML"):
    """Variance ratio d = s2_e / s2_g by REML (or ML) in the rotated space, with the fit's diagnostics."""
    n, p = X0s.shape
    grid = np.linspace(LOG_LO, LOG_HI, GRID_POINTS)
    lls = np.array([_null_eval(ys, X0s, w, 10 ** t, mode)[0] for t in grid])
    i = int(np.argmax(lls))
    lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
    best_t, best_ll = grid[i], lls[i]
    try:
        r = minimize_scalar(lambda t: -_null_eval(ys, X0s, w, 10 ** t, mode)[0], bounds=(lo, hi), method="bounded",
                            options={"xatol": 1e-7})
        ok = bool(r.success)
        if -r.fun >= best_ll:
            best_t, best_ll = float(r.x), float(-r.fun)
    except Exception:
        ok = False
    delta = 10 ** best_t
    ll, beta, s2, rss = _null_eval(ys, X0s, w, delta, mode)
    edge = "lower" if i == 0 else "upper" if i == len(grid) - 1 else None
    # independent check: a much finer grid must not find a better likelihood than the one reported
    fine = np.linspace(LOG_LO, LOG_HI, 481)
    ll_fine = max(_null_eval(ys, X0s, w, 10 ** t, mode)[0] for t in fine)
    return {"mode": mode, "delta": float(delta), "sigma2_g": float(s2), "sigma2_e": float(s2 * delta),
            "h2": float(1.0 / (1.0 + delta)), "loglik": float(ll), "beta": beta, "n": int(n), "n_fixed": int(p),
            "converged": ok and edge is None, "at_bound": edge, "check_ok": bool(ll >= ll_fine - 1e-6),
            "check_best_grid_loglik": float(ll_fine)}


# ============================================================================================ per-marker REML / ML

def _batch_eval(ys, X0s, Gb, w, delta, mode):
    """Log-likelihood, effect and SE factor for each marker at its own variance ratio. Returns arrays over markers."""
    n, p0 = X0s.shape
    m = Gb.shape[1]
    p = p0 + 1
    h = 1.0 / (w[:, None] + delta[None, :])
    XX = (X0s[:, :, None] * X0s[:, None, :]).reshape(n, p0 * p0)
    A = (h.T @ XX).reshape(m, p0, p0)
    hg = h * Gb
    b = hg.T @ X0s
    c = (hg * Gb).sum(0)
    r0 = h.T @ (X0s * ys[:, None])
    s = hg.T @ ys
    yy = h.T @ (ys * ys)
    M = np.zeros((m, p, p))
    M[:, :p0, :p0] = A
    M[:, :p0, p0] = b
    M[:, p0, :p0] = b
    M[:, p0, p0] = c
    rhs = np.concatenate([r0, s[:, None]], axis=1)
    try:
        Minv = np.linalg.inv(M)
    except np.linalg.LinAlgError:
        Minv = np.linalg.pinv(M)
    beta = np.einsum("mij,mj->mi", Minv, rhs)
    rss = np.maximum(yy - (beta * rhs).sum(1), 1e-300)
    sumlog = np.log(w[:, None] + delta[None, :]).sum(0)
    with np.errstate(all="ignore"):
        logdetM = np.linalg.slogdet(M)[1]
        if mode == "REML":
            s2 = rss / (n - p)
            ll = -0.5 * ((n - p) * np.log(s2) + sumlog + logdetM)
        else:
            s2 = rss / n
            ll = -0.5 * (n * np.log(s2) + sumlog)
    ll = np.where(np.isfinite(ll), ll, -np.inf)
    return ll, beta[:, p0], Minv[:, p0, p0], s2


def _opt_delta(ys, X0s, Gb, w, mode):
    """Per-marker variance ratio: coarse log grid, then four rounds of local refinement. Returns (delta, at_edge)."""
    m = Gb.shape[1]
    grid = np.linspace(LOG_LO, LOG_HI, GRID_POINTS)
    best_ll = np.full(m, -np.inf)
    best_t = np.full(m, grid[0])
    for t in grid:
        ll = _batch_eval(ys, X0s, Gb, w, np.full(m, 10 ** t), mode)[0]
        better = ll > best_ll
        best_ll[better] = ll[better]
        best_t[better] = t
    edge = (best_t <= grid[0] + 1e-9) | (best_t >= grid[-1] - 1e-9)
    width = grid[1] - grid[0]
    for _ in range(4):
        centre = best_t.copy()
        for off in np.linspace(-width, width, 9):
            t = np.clip(centre + off, LOG_LO, LOG_HI)
            ll = _batch_eval(ys, X0s, Gb, w, 10 ** t, mode)[0]
            better = ll > best_ll
            best_ll[better] = ll[better]
            best_t[better] = t[better]
        width /= 4.0
    return 10 ** best_t, edge


# ============================================================================================ scans

def _project(Xw):
    Q, _ = np.linalg.qr(Xw)
    return Q


def _t_test(gy, gg, yy, df, sigma2_fixed=None):
    """Effect, SE, statistic and p from the projected cross-products (all arrays over markers)."""
    with np.errstate(all="ignore"):
        beta = gy / gg
        if sigma2_fixed is None:
            s2 = np.maximum(yy - beta * gy, 0.0) / df
            se = np.sqrt(s2 / gg)
            stat = beta / se
            p = 2.0 * stats.t.sf(np.abs(stat), df)
        else:
            se = np.sqrt(sigma2_fixed / gg)
            stat = beta / se
            p = 2.0 * stats.norm.sf(np.abs(stat))
    return beta, se, stat, p


def _blocks(m, chunk):
    for s in range(0, m, chunk):
        yield slice(s, min(s + chunk, m))


def ols_scan(y, X0, G, chunk=1000):
    """Ordinary least squares for every marker, adjusting for the columns of X0 (intercept + any covariates)."""
    y = np.asarray(y, float)
    n, p0 = X0.shape
    df = n - p0 - 1
    if df < 3:
        raise ModelError("Too few samples for the number of covariates.")
    Q = _project(X0)
    yr = y - Q @ (Q.T @ y)
    m = G.shape[1]
    beta = np.full(m, np.nan); se = beta.copy(); stat = beta.copy(); p = beta.copy()
    for sl in _blocks(m, chunk):
        Gb = G[:, sl]
        Gr = Gb - Q @ (Q.T @ Gb)
        gg = (Gr ** 2).sum(0)
        ok = gg > 1e-10 * (Gb ** 2).sum(0)
        gy = Gr.T @ yr
        gg_safe = np.where(ok, gg, np.nan)
        b, s, t, pv = _t_test(gy, gg_safe, float(yr @ yr), df)
        beta[sl], se[sl], stat[sl], p[sl] = b, s, t, pv
    return _pack(beta, se, stat, p, "ols", "wald", df=df, extra={"covariates": p0 - 1})


def _pack(beta, se, stat, p, method, test, df=None, delta=None, null=None, eig=None, extra=None):
    invalid = ~np.isfinite(p)
    info = {"method": method, "test": test, "stat_name": ("chi2" if test == "lrt" else STAT_NAMES[method]), "df": df,
            "n_tested": int((~invalid).sum()), "n_invalid": int(invalid.sum())}
    if null is not None:
        info["null"] = {k: v for k, v in null.items() if k != "beta"}
    if eig is not None:
        info["kinship_eigen"] = eig
    info.update(extra or {})
    return {"beta": beta, "se": se, "stat": stat, "p": p, "delta": delta if delta is not None else np.full(len(p), np.nan),
            "invalid": invalid, "info": info}


def mixed_scan(y, X0, G, K, method="mlm", test="wald", chunk=500, eig=None):
    """Mixed-model scan of every column of G with the relationship matrix K (see the module description)."""
    if method not in ("mlm", "emmax", "gemma"):
        raise ModelError(f"Unknown mixed model '{method}'.")
    if test not in ("wald", "lrt"):
        raise ModelError("test must be 'wald' or 'lrt'.")
    if test == "lrt" and method != "gemma":
        raise ModelError("The likelihood-ratio test belongs to the per-marker REML model only.")
    y = np.asarray(y, float)
    n, p0 = X0.shape
    df = n - p0 - 1
    if df < 3:
        raise ModelError("Too few samples for the number of covariates.")
    w, U, einfo = eig if eig is not None else eig_psd(K)
    ys, X0s = U.T @ y, U.T @ X0
    nf = null_fit(ys, X0s, w, "REML")
    m = G.shape[1]
    beta = np.full(m, np.nan); se = beta.copy(); stat = beta.copy(); p = beta.copy(); delta = np.full(m, np.nan)
    extra = {}
    if method in ("mlm", "emmax"):
        sw = 1.0 / np.sqrt(w + nf["delta"])
        Q = _project(X0s * sw[:, None])
        yw = ys * sw
        yr = yw - Q @ (Q.T @ yw)
        yy = float(yr @ yr)
        for sl in _blocks(m, chunk):
            Gw = (U.T @ G[:, sl]) * sw[:, None]
            Gr = Gw - Q @ (Q.T @ Gw)
            gg = (Gr ** 2).sum(0)
            ok = gg > 1e-10 * (Gw ** 2).sum(0)
            gy = Gr.T @ yr
            b, s, t, pv = _t_test(gy, np.where(ok, gg, np.nan), yy, df, nf["sigma2_g"] if method == "emmax" else None)
            beta[sl], se[sl], stat[sl], p[sl] = b, s, t, pv
            delta[sl] = nf["delta"]
    else:
        n_edge = 0
        nf0 = null_fit(ys, X0s, w, "ML") if test == "lrt" else None
        for sl in _blocks(m, chunk):
            Gs = U.T @ G[:, sl]
            if test == "wald":
                d_j, edge = _opt_delta(ys, X0s, Gs, w, "REML")
                ll, b, sef, s2 = _batch_eval(ys, X0s, Gs, w, d_j, "REML")
                sej = np.sqrt(s2 * sef)
                t = b / sej
                beta[sl], se[sl], stat[sl] = b, sej, t
                p[sl] = 2.0 * stats.t.sf(np.abs(t), df)
            else:
                d_j, edge = _opt_delta(ys, X0s, Gs, w, "ML")
                ll1, b, sef, s2 = _batch_eval(ys, X0s, Gs, w, d_j, "ML")
                lr = np.maximum(2.0 * (ll1 - nf0["loglik"]), 0.0)
                beta[sl], se[sl], stat[sl] = b, np.sqrt(s2 * sef), lr
                p[sl] = stats.chi2.sf(lr, 1)
            delta[sl] = d_j
            n_edge += int(edge.sum())
        extra = {"n_delta_at_grid_edge": n_edge,
                 "note": "Variance ratio re-estimated for every marker; Wald test uses REML fits, the likelihood-ratio test uses "
                         "maximum-likelihood fits of the null and the alternative model."}
        if nf0 is not None:
            extra["null_ml"] = {k: v for k, v in nf0.items() if k != "beta"}
    return _pack(beta, se, stat, p, method, test, df=(1 if test == "lrt" else df), delta=delta, null=nf, eig=einfo, extra=extra)


def mixed_scan_loco(y, X0, G, chrom, loco_iter, method="mlm", test="wald", chunk=500):
    """Leave-one-chromosome-out: markers of chromosome c are tested with the kinship built without c."""
    chrom = np.asarray(chrom)
    m = G.shape[1]
    out = None
    per_chr = {}
    for c, Kc, kinfo in loco_iter:
        idx = np.flatnonzero(chrom == c)
        if len(idx) == 0:
            continue
        r = mixed_scan(y, X0, G[:, idx], Kc, method, test, chunk)
        if out is None:
            out = {k: np.full(m, np.nan) for k in ("beta", "se", "stat", "p", "delta")}
        for k in ("beta", "se", "stat", "p", "delta"):
            out[k][idx] = r[k]
        per_chr[int(c)] = {"n_snps_in_kinship": kinfo["n_snps"], "delta": r["info"]["null"]["delta"], "h2": r["info"]["null"]["h2"],
                           "n_tested": r["info"]["n_tested"]}
        proto = r
    if out is None:
        raise ModelError("Leave-one-chromosome-out needs at least two chromosomes with markers.")
    res = _pack(out["beta"], out["se"], out["stat"], out["p"], method, test, df=proto["info"]["df"], delta=out["delta"],
                extra={"loco": True, "per_chromosome": per_chr})
    return res


# ============================================================================================ simulation-based threshold

NULL_METHODS = ("residual_permutation", "parametric_bootstrap")
LIMITATIONS = (
    "The threshold is specific to this model, this data set and this covariate set; do not apply it to another model.",
    "It is a Monte Carlo estimate. With few replicates the (1 - a) quantile of the maxima is imprecise; the interval shows by how much.",
    "Residual permutation assumes the fitted null model (kinship + covariates) describes the correlation among samples; "
    "an unmodelled structure remains in the replicates as it is in the data.",
)


def null_threshold(y, X0, G, K=None, method="mlm", n_perm=1000, alphas=(0.05,), null="residual_permutation", seed=1,
                   chunk=500, eig=None, progress=None):
    """Genome-wide threshold on the p-value scale from replicate phenotypes generated under the fitted null model."""
    if method not in METHODS:
        raise ModelError(f"No significance simulation for '{method}'.")
    if null not in NULL_METHODS:
        raise ModelError(f"Unknown null model '{null}'.")
    y = np.asarray(y, float)
    n, p0 = X0.shape
    df = n - p0 - 1
    R = int(n_perm)
    if R < 20:
        raise ModelError("At least 20 replicates are needed for a threshold.")
    rng = np.random.default_rng(seed)
    if method == "ols":
        U, sw = None, np.ones(n)
        ys, Xs = y, X0
        Q = _project(Xs)
        v = ys - Q @ (Q.T @ ys)
        sigma2 = float(v @ v) / (n - p0)
        nf = None
    else:
        if K is None and eig is None:
            raise ModelError("A relationship matrix is required.")
        w, U, _ = eig if eig is not None else eig_psd(K)
        ys, Xs = U.T @ y, U.T @ X0
        nf = null_fit(ys, Xs, w, "REML")
        sw = 1.0 / np.sqrt(w + nf["delta"])
        Xw = Xs * sw[:, None]
        Q = _project(Xw)
        v = ys * sw - Xw @ np.linalg.lstsq(Xw, ys * sw, rcond=None)[0]
        sigma2 = nf["sigma2_g"]
    if null == "residual_permutation":
        V = np.empty((n, R))
        for r in range(R):
            V[:, r] = v[rng.permutation(n)]
    else:
        V = math.sqrt(sigma2) * rng.standard_normal((n, R))
    Yr = V - Q @ (Q.T @ V)
    yy = (Yr ** 2).sum(0)
    fixed_var = method == "emmax"
    maxabs = np.zeros(R)
    m = G.shape[1]
    done = 0
    for sl in _blocks(m, chunk):
        Gb = G[:, sl]
        Gw = (Gb if U is None else U.T @ Gb) * sw[:, None]
        Gr = Gw - Q @ (Q.T @ Gw)
        gg = (Gr ** 2).sum(0)
        ok = gg > 1e-10 * (Gw ** 2).sum(0)
        GY = Gr.T @ Yr
        with np.errstate(all="ignore"):
            if fixed_var:
                stat = GY / (np.sqrt(gg)[:, None] * math.sqrt(sigma2))
            else:
                rss = np.maximum(yy[None, :] - GY ** 2 / gg[:, None], 1e-300)
                stat = GY / np.sqrt(gg[:, None] * rss / df)
        stat[~ok, :] = 0.0
        maxabs = np.maximum(maxabs, np.nan_to_num(np.abs(stat), nan=0.0).max(0))
        done += sl.stop - sl.start
        if progress:
            progress(done / m)
    good = np.isfinite(maxabs) & (maxabs > 0)
    mx = np.sort(maxabs[good])[::-1]
    Rg = len(mx)
    if Rg < 20:
        raise ModelError("Too few valid replicates to set a threshold.")

    def to_p(s):
        return float(2.0 * (stats.norm.sf(s) if fixed_var else stats.t.sf(s, df)))

    out = {}
    notes = []
    for a in alphas:
        k = max(int(math.ceil(a * Rg)), 1)
        s_thr = float(mx[k - 1])
        lo = int(max(stats.binom.ppf(0.025, Rg, a), 0))
        hi = int(min(stats.binom.ppf(0.975, Rg, a) + 1, Rg - 1))
        s_hi, s_lo = float(mx[lo]), float(mx[hi])
        out[f"{a:g}"] = {"alpha": float(a), "p": to_p(s_thr), "p_stricter": to_p(s_hi), "p_looser": to_p(s_lo),
                          "statistic": s_thr, "neg_log10_p": -math.log10(max(to_p(s_thr), 1e-300))}
        if a * Rg < 10:
            notes.append(f"Only {a * Rg:.0f} replicates lie beyond the {a:g} threshold; it is imprecise (use more replicates).")
    return {"method": method, "null": null, "n_requested": R, "n_successful": int(Rg), "seed": int(seed), "thresholds": out,
            "statistic_name": STAT_NAMES[method], "df": None if fixed_var else df,
            "max_stat_quantiles": {q: float(np.quantile(maxabs[good], q)) for q in (0.5, 0.9, 0.95, 0.99)},
            "limitations": list(LIMITATIONS) + notes
            + (["For the per-marker REML model the variance ratio is held at the null-model estimate in the replicates "
                "(approximation)."] if method == "gemma" else [])}
