"""Multiple testing, significance categories, diagnostics and PVE."""
import numpy as np
import pytest
from scipy import stats

from tools.gwas.backend import inference as inf


def res(p, kind="single", **kw):
    p = np.asarray(p, float)
    return inf.ModelResult(name="M", p=p, beta=np.zeros(len(p)), se=np.ones(len(p)), kind=kind, **kw)


def test_bonferroni_and_bh_match_definitions():
    p = np.array([0.001, 0.01, 0.02, 0.04, 0.5, np.nan])
    bonf, bh = inf.adjust(p)
    assert np.allclose(bonf[:5], np.minimum(p[:5] * 5, 1)) and np.isnan(bonf[5]) and np.isnan(bh[5])
    # BH by hand: p*m/rank then running minimum from the end
    manual = np.minimum.accumulate((p[:5] * 5 / np.arange(1, 6))[::-1])[::-1]
    assert np.allclose(bh[:5], manual)


def test_bh_is_monotone_and_never_below_raw_p():
    rng = np.random.default_rng(1)
    p = rng.uniform(size=500)
    _, bh = inf.adjust(p)
    assert (bh >= p - 1e-15).all() and (bh <= 1).all()
    o = np.argsort(p)
    assert (np.diff(bh[o]) >= -1e-15).all()


def test_categories_follow_the_chosen_rule_only():
    p = np.concatenate([[1e-9, 2e-4, 3e-4], np.random.default_rng(0).uniform(0.01, 1, 997)])
    r = res(p)
    s = inf.evaluate(r, alpha=0.05, rule="bonferroni")
    assert list(s["category"][:3]) == ["Genome-wide significant", "Suggestive", "Suggestive"]
    assert s["genome_wide"][0] and not s["genome_wide"][1]
    assert s["thresholds"]["bonferroni_p"] == pytest.approx(0.05 / 1000)
    s_bh = inf.evaluate(r, alpha=0.05, rule="bh")
    assert s_bh["rule_used"] == "bh" and s_bh["flag_bonferroni"].sum() == s["flag_bonferroni"].sum()


def test_permutation_rule_without_threshold_falls_back_and_says_so():
    s = inf.evaluate(res(np.linspace(1e-8, 1, 100)), rule="permutation")
    assert s["rule_used"] == "bonferroni" and s["notes"]


def test_top_ranked_markers_are_never_significant():
    p = np.random.default_rng(2).uniform(0.05, 1, 2000)            # nothing near 2.5e-5
    r = res(p)
    s = inf.evaluate(r)
    assert s["genome_wide"].sum() == 0
    ex = inf.exploratory_idx(r, 10)
    assert len(ex) == 10 and not s["genome_wide"][ex].any()
    assert np.array_equal(np.sort(p[ex]), np.sort(p)[:10])
    assert inf.diagnostics(r)["n_genome_wide"] == 0
    assert "No genome-wide significant" in inf.NO_HITS


def test_multilocus_keeps_its_own_qtn_rule_and_no_p_is_invented():
    p = np.random.default_rng(3).uniform(0.01, 1, 500)
    r = res(p, kind="multilocus", qtn_idx=np.array([7, 40]), qtn_lod=np.array([4.1, 3.5]))
    s = inf.evaluate(r)
    assert set(np.flatnonzero(s["genome_wide"])) == {7, 40} and s["category"][7] == "QTN"
    assert np.allclose(r.p, p)                                      # LOD never converted back to p


def test_lambda_is_one_for_uniform_p_and_is_not_applied():
    p = np.random.default_rng(4).uniform(size=20000)
    assert inf.genomic_lambda(p) == pytest.approx(1.0, abs=0.04)
    chi = stats.chi2.isf(p, 1) * 2.0
    assert inf.genomic_lambda(stats.chi2.sf(chi, 1)) == pytest.approx(2.0, abs=0.1)
    r = res(p)
    before = r.p.copy()
    inf.diagnostics(r)
    assert np.array_equal(r.p, before)                              # no genomic control


def test_invalid_statistics_are_counted_not_tested():
    r = res([0.5, np.nan, 0.0, 1.5, 0.01])
    d = inf.diagnostics(r)
    assert d["snps_tested"] == 2 and d["invalid_statistics"] == 3


def test_pve_ols_marginal_equals_r2_gain():
    rng = np.random.default_rng(5)
    n = 200
    X0 = np.column_stack([np.ones(n), rng.normal(size=n)])
    G = rng.binomial(2, 0.3, size=(n, 4)).astype(float)
    y = X0 @ [1, 2] + 0.8 * G[:, 1] + rng.normal(size=n)
    out = inf.pve_ols_marginal(y, X0, G)
    for j in range(4):
        def rss(X):
            b = np.linalg.lstsq(X, y, rcond=None)[0]
            return ((y - X @ b) ** 2).sum()
        expect = (rss(X0) - rss(np.column_stack([X0, G[:, j]]))) / ((y - y.mean()) ** 2).sum()
        assert out[j] == pytest.approx(expect, rel=1e-8, abs=1e-12)
    assert out[1] == out.max()


def test_pve_is_unavailable_for_a_constant_marker_or_phenotype():
    rng = np.random.default_rng(6)
    n = 50
    X0 = np.ones((n, 1))
    G = np.column_stack([np.ones(n), rng.binomial(2, 0.4, n)]).astype(float)
    out = inf.pve_ols_marginal(rng.normal(size=n), X0, G)
    assert np.isnan(out[0]) and np.isfinite(out[1])
    assert np.isnan(inf.pve_ols_marginal(np.ones(n), X0, G)).all()


def test_pve_conditional_is_additive_for_orthogonal_markers_and_bounded():
    rng = np.random.default_rng(7)
    n = 400
    G = rng.binomial(2, 0.4, size=(n, 3)).astype(float)
    y = 0.6 * G[:, 0] + 0.4 * G[:, 2] + rng.normal(size=n)
    per, joint = inf.pve_conditional(y, np.ones((n, 1)), G, [0, 2])
    assert (per >= 0).all() and 0 <= joint <= 1
    assert per.sum() <= joint + 0.02
    assert np.isnan(inf.pve_conditional(y, np.ones((n, 1)), G, [])[0]).size == 0
