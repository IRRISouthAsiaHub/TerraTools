"""mrMLM and FASTmrMLM come from the mrMLM R package; these tests check the files TerraTools writes and, when the package is
installed in the bundled R, one real run."""
import numpy as np
import pandas as pd
import pytest

from tools.gwas.backend import mrmlm_package as mp


def test_input_files_follow_the_package_format(tmp_path):
    rng = np.random.default_rng(0)
    G = rng.binomial(2, 0.4, size=(6, 4)).astype(float)
    gen, phe, ps = mp._write_inputs(str(tmp_path), G, [f"s{i}" for i in range(6)], ["Chr1_10", "Chr1_20", "Chr2_5", "Chr2_9"],
                                    [1, 1, 2, 2], [10, 20, 5, 9], rng.normal(size=6), ["T", "T", "G", "A"], rng.normal(size=(6, 2)))
    g = pd.read_csv(gen)
    assert list(g.columns[:4]) == ["rs#", "chrom", "pos", "genotype for code 1"] and list(g.columns[4:]) == [f"s{i}" for i in range(6)]
    assert np.array_equal(g.iloc[:, 4:].to_numpy(), G.T - 1.0)                    # dosage 0/1/2 written as -1/0/1
    lines = open(phe).read().splitlines()
    assert lines[0] == "<Phenotype>,trait1" and len(lines) == 7
    ps_lines = open(ps).read().splitlines()
    assert ps_lines[0].startswith("<PopStr>") and ps_lines[1] == "<ID>,PC1,PC2" and ps_lines[2].startswith("s0,")


def test_no_structure_file_without_covariates(tmp_path):
    G = np.zeros((4, 2))
    _, _, ps = mp._write_inputs(str(tmp_path), G, list("abcd"), ["a", "b"], [1, 1], [1, 2], np.arange(4.0), ["A", "A"], None)
    assert ps == "none"


def test_duplicate_marker_names_are_refused():
    with pytest.raises(mp.MrmlmError):
        mp.run(np.zeros((3, 2)), list("abc"), ["x", "x"], [1, 1], [1, 2], np.arange(3.0))


@pytest.mark.skipif(mp.package_version() is None, reason="mrMLM R package not installed")
def test_real_package_run_finds_a_planted_marker():
    from tools.gwas.validation.simulate import Scenario, simulate
    sim = simulate(Scenario("t", n=150, m=300, n_hidden=100, chromosomes=3, n_causal=1, causal_pve=0.35), 4)
    d = sim.data
    G = np.where(np.isnan(d.G), np.nanmean(d.G, 0), d.G)
    out = mp.run(G, d.ids, d.snp_ids, d.chrom, d.pos, d.y, methods=("FASTmrMLM",))
    assert out["version"]
    nlp = out["FASTmrMLM"]["neglog10p"]
    j = list(d.snp_ids).index(sim.truth["causal_names"][0])
    assert nlp[j] >= np.nanpercentile(nlp, 99)
