"""The simulation, the external-GEMMA comparison helpers and the benchmark inputs."""
import numpy as np

from tools.gwas.backend import analysis as an, gemma_external as ge
from tools.gwas.config import GWASConfig
from tools.gwas.validation.simulate import SCENARIOS, Scenario, simulate


def small(**kw):
    return Scenario("t", n=120, m=400, n_hidden=100, chromosomes=4, **kw)


def test_simulation_is_reproducible_and_has_the_requested_features():
    a, b = simulate(small(n_causal=2), 5), simulate(small(n_causal=2), 5)
    assert np.array_equal(a.data.G, b.data.G, equal_nan=True) and np.array_equal(a.data.y, b.data.y)
    assert not np.array_equal(a.data.y, simulate(small(n_causal=2), 6).data.y)
    assert np.isnan(a.data.G).mean() > 0.02                      # missing calls present
    assert len(a.causal) == 2 and a.data.marker_map.valid
    maf = np.minimum(np.nanmean(a.data.G, 0) / 2, 1 - np.nanmean(a.data.G, 0) / 2)
    assert maf.min() < 0.05 and maf.max() > 0.3                  # a spread of allele frequencies
    for g in ("indica", "japonica", "aus"):
        assert (a.labels == g).sum() > 10


def test_structure_inflates_ols_and_mixed_model_stays_calibrated():
    sim = simulate(small(), 3)
    cfg = GWASConfig(models=["OLS", "MLM"], permutations=0, compare_pcs="", n_pcs=3)
    A = an.analyse(cfg, sim.data, labels=sim.labels, log=lambda m: None)
    assert A.results["OLS"].info["diag"]["lambda"] > 3
    assert 0.8 < A.results["MLM"].info["diag"]["lambda"] < 1.25


def test_planted_marker_is_found_by_the_mixed_model():
    sim = simulate(small(n_causal=1, causal_pve=0.3), 12)
    cfg = GWASConfig(models=["MLM"], permutations=0, compare_pcs="", n_pcs=3)
    A = an.analyse(cfg, sim.data, labels=sim.labels, log=lambda m: None)
    j = list(A.snp).index(sim.truth["causal_names"][0])
    assert A.results["MLM"].info["sig"]["genome_wide"][j]


def test_external_comparison_numbers():
    rng = np.random.default_rng(0)
    p = rng.uniform(size=500) ** 3
    same = ge.compare(p, p)
    assert same["agrees"] and same["spearman"] == 1.0 and same["max_abs_delta_neglog10p"] == 0
    shuffled = ge.compare(p, rng.permutation(p))
    assert not shuffled["agrees"]


def test_missing_gemma_binary_is_reported_not_faked(monkeypatch):
    monkeypatch.delenv("GEMMA", raising=False)
    monkeypatch.setattr(ge.shutil, "which", lambda _: None)
    r = ge.run_gemma(np.zeros((3, 3)), np.zeros(3), ["a", "b", "c"], np.ones(3))
    assert r["status"] == "not run" and "Install" in r["message"] or "install" in r["message"].lower()


def test_benchmark_inputs_are_written(tmp_path):
    from tools.gwas.benchmark import make_inputs
    SCENARIOS["tiny"] = small(n_causal=1)
    try:
        make_inputs.main(["--scenario", "tiny", "--out", str(tmp_path)])
    finally:
        del SCENARIOS["tiny"]
    for f in ("numeric_genotypes.csv", "map.csv", "phenotype.csv", "kinship.csv", "pcs.csv", "truth.json", "terratools_MLM.csv"):
        assert (tmp_path / f).exists(), f
