"""Mixed-model and linear-model scans checked against dense brute-force computations."""
import numpy as np
import pytest
from scipy import stats
from scipy.optimize import minimize_scalar

from tools.gwas.backend import mixed, population as pop


def simulate(n=60, m=40, seed=0, h2=0.5, effect=0.0):
    rng = np.random.default_rng(seed)
    G = rng.binomial(2, rng.uniform(0.15, 0.5, m), size=(n, m)).astype(float)
    K, _ = pop.vanraden(G)
    L = np.linalg.cholesky(K + 1e-6 * np.eye(n))
    u = L @ rng.standard_normal(n) * np.sqrt(h2)
    y = u + rng.standard_normal(n) * np.sqrt(1 - h2) + effect * G[:, 0] + 2.0
    return y, G, K


def dense_reml(y, X, K, delta):
    n, p = X.shape
    V = K + delta * np.eye(n)
    Vi = np.linalg.inv(V)
    A = X.T @ Vi @ X
    b = np.linalg.solve(A, X.T @ Vi @ y)
    r = y - X @ b
    s2 = float(r @ Vi @ r) / (n - p)
    S = s2 * V
    Si = np.linalg.inv(S)
    ll = (-0.5 * (n - p) * np.log(2 * np.pi) - 0.5 * np.linalg.slogdet(S)[1] - 0.5 * np.linalg.slogdet(X.T @ Si @ X)[1]
          + 0.5 * np.linalg.slogdet(X.T @ X)[1] - 0.5 * float(r @ Si @ r))
    return ll, b, s2


def dense_gls_scan(y, X0, K, delta, g, s2_fixed=None):
    n = len(y)
    X = np.column_stack([X0, g])
    V = K + delta * np.eye(n)
    Vi = np.linalg.inv(V)
    A = X.T @ Vi @ X
    b = np.linalg.solve(A, X.T @ Vi @ y)
    r = y - X @ b
    s2 = float(r @ Vi @ r) / (n - X.shape[1]) if s2_fixed is None else s2_fixed
    return b[-1], np.sqrt(s2 * np.linalg.inv(A)[-1, -1])


def test_null_reml_matches_dense_and_grid():
    y, G, K = simulate()
    X0 = np.ones((len(y), 1))
    w, U, _ = pop.eig_psd(K)
    nf = mixed.null_fit(U.T @ y, U.T @ X0, w)
    assert nf["check_ok"]
    ll_dense, _, s2 = dense_reml(y, X0, K, nf["delta"])
    assert nf["loglik"] == pytest.approx(ll_dense, abs=1e-6)
    assert nf["sigma2_g"] == pytest.approx(s2, rel=1e-6)
    best = minimize_scalar(lambda t: -dense_reml(y, X0, K, 10 ** t)[0], bounds=(-6, 6), method="bounded", options={"xatol": 1e-8})
    assert np.log10(nf["delta"]) == pytest.approx(best.x, abs=2e-3)


def test_mlm_and_emmax_match_dense_gls():
    y, G, K = simulate(seed=3, effect=0.8)
    X0 = np.column_stack([np.ones(len(y)), np.random.default_rng(1).standard_normal(len(y))])
    for method in ("mlm", "emmax"):
        r = mixed.mixed_scan(y, X0, G, K, method)
        d = r["info"]["null"]["delta"]
        s2 = r["info"]["null"]["sigma2_g"]
        for j in (0, 5, 17):
            b, se = dense_gls_scan(y, X0, K, d, G[:, j], s2 if method == "emmax" else None)
            assert r["beta"][j] == pytest.approx(b, rel=1e-6)
            assert r["se"][j] == pytest.approx(se, rel=1e-6)
            stat = b / se
            p = 2 * (stats.norm.sf(abs(stat)) if method == "emmax" else stats.t.sf(abs(stat), len(y) - 3))
            assert r["p"][j] == pytest.approx(p, rel=1e-5)


def test_emmax_and_mlm_differ_only_in_variance_handling():
    y, G, K = simulate(seed=4, effect=0.5)
    X0 = np.ones((len(y), 1))
    a, b = mixed.mixed_scan(y, X0, G, K, "mlm"), mixed.mixed_scan(y, X0, G, K, "emmax")
    assert np.allclose(a["beta"], b["beta"])           # same fixed-effect estimates
    assert not np.allclose(a["se"], b["se"])           # different SE (s2 re-estimated per marker vs fixed)
    assert a["info"]["stat_name"] == "t" and b["info"]["stat_name"] == "z"


def test_gemma_style_wald_matches_dense_per_marker_reml():
    y, G, K = simulate(n=50, m=12, seed=5, effect=0.7)
    X0 = np.ones((len(y), 1))
    r = mixed.mixed_scan(y, X0, G, K, "gemma", "wald")
    for j in (0, 3, 9):
        X = np.column_stack([X0, G[:, j]])
        best = minimize_scalar(lambda t: -dense_reml(y, X, K, 10 ** t)[0], bounds=(-5, 5), method="bounded", options={"xatol": 1e-9})
        b, se = dense_gls_scan(y, X0, K, 10 ** best.x, G[:, j])
        assert r["delta"][j] == pytest.approx(10 ** best.x, rel=0.03)
        assert r["beta"][j] == pytest.approx(b, rel=0.02)
        assert r["se"][j] == pytest.approx(se, rel=0.02)


def test_gemma_lrt_is_nonnegative():
    y, G, K = simulate(n=60, m=60, seed=6)
    r = mixed.mixed_scan(y, np.ones((60, 1)), G, K, "gemma", "lrt")
    assert (r["stat"] >= 0).all() and r["info"]["stat_name"] == "chi2"
    assert np.isfinite(r["p"]).all()


def test_lrt_only_for_gemma():
    y, G, K = simulate()
    with pytest.raises(mixed.ModelError):
        mixed.mixed_scan(y, np.ones((len(y), 1)), G, K, "mlm", "lrt")


def test_ols_matches_statsmodels_with_covariates():
    import statsmodels.api as sm
    y, G, _ = simulate(n=80, m=10, seed=7, effect=0.6)
    cov = np.random.default_rng(2).standard_normal((80, 2))
    X0, names, dropped = mixed.design_matrix(80, cov, ["a", "b"])
    r = mixed.ols_scan(y, X0, G)
    for j in range(10):
        fit = sm.OLS(y, np.column_stack([X0, G[:, j]])).fit()
        assert r["beta"][j] == pytest.approx(fit.params[-1], rel=1e-8)
        assert r["se"][j] == pytest.approx(fit.bse[-1], rel=1e-8)
        assert r["p"][j] == pytest.approx(fit.pvalues[-1], rel=1e-6)


def test_redundant_covariates_are_dropped():
    n = 30
    a = np.random.default_rng(0).standard_normal(n)
    X0, used, dropped = mixed.design_matrix(n, np.column_stack([a, 2 * a + 1, np.ones(n)]), ["a", "dup", "const"])
    assert used == ["Intercept", "a"] and {d[0] for d in dropped} == {"dup", "const"}
    with pytest.raises(mixed.ModelError):
        mixed.design_matrix(n, np.array([np.nan] + [1.0] * (n - 1)), ["bad"])


def test_loco_uses_kinship_without_the_chromosome():
    y, G, K = simulate(n=60, m=60, seed=8)
    chrom = np.repeat([1, 2, 3], 20)
    r = mixed.mixed_scan_loco(y, np.ones((60, 1)), G, chrom, pop.loco_kinship(G, chrom), "mlm")
    assert r["info"]["loco"] and np.isfinite(r["p"]).all()
    Kc, _ = pop.vanraden(G[:, chrom != 2])
    direct = mixed.mixed_scan(y, np.ones((60, 1)), G[:, chrom == 2], Kc, "mlm")
    assert np.allclose(r["p"][chrom == 2], direct["p"], rtol=1e-6)


def test_threshold_is_calibrated_under_the_null():
    """About 5% of null data sets should have a marker beyond the 5% genome-wide threshold."""
    hits, trials = 0, 40
    for s in range(trials):
        y, G, K = simulate(n=60, m=60, seed=100 + s, h2=0.5)
        X0 = np.ones((60, 1))
        thr = mixed.null_threshold(y, X0, G, K, "mlm", n_perm=300, seed=s)["thresholds"]["0.05"]["p"]
        p = mixed.mixed_scan(y, X0, G, K, "mlm")["p"]
        hits += int(p.min() <= thr)
    assert hits / trials <= 0.2          # wide band: 40 trials around 5%


def test_threshold_reports_replicates_and_limits():
    y, G, K = simulate()
    out = mixed.null_threshold(y, np.ones((len(y), 1)), G, K, "emmax", n_perm=100, alphas=(0.05, 0.1), null="parametric_bootstrap")
    assert out["n_requested"] == 100 and out["n_successful"] <= 100 and set(out["thresholds"]) == {"0.05", "0.1"}
    assert out["limitations"] and out["thresholds"]["0.05"]["p"] < out["thresholds"]["0.1"]["p"]
    assert out["thresholds"]["0.05"]["p_stricter"] <= out["thresholds"]["0.05"]["p_looser"]
    with pytest.raises(mixed.ModelError):
        mixed.null_threshold(y, np.ones((len(y), 1)), G, K, "mlm", n_perm=5)


def test_kinship_matches_brute_force():
    rng = np.random.default_rng(0)
    G = rng.binomial(2, 0.3, size=(25, 40)).astype(float)
    G[rng.random(G.shape) < 0.05] = np.nan
    K, info = pop.vanraden(G)
    p = np.nanmean(G, 0) / 2
    n, m = G.shape
    brute = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            s = 0.0
            for k in range(m):
                zi = 0.0 if np.isnan(G[i, k]) else G[i, k] - 2 * p[k]
                zj = 0.0 if np.isnan(G[j, k]) else G[j, k] - 2 * p[k]
                s += zi * zj
            brute[i, j] = s / (2 * np.sum(p * (1 - p)))
    assert np.allclose(K, brute) and np.allclose(K, K.T)
    w, _, e = pop.eig_psd(K)
    assert e["psd"] and (w >= 0).all() and info["n_snps"] == m
