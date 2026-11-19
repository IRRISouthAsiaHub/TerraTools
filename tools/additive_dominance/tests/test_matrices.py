"""Relationship matrices: the additive matrix is VanRaden (2008), the split mode uses the training lines only,
and the dominance matrix is Nishio & Satoh (2014) or Vitezica et al. (2013)."""
import json
import subprocess

import numpy as np
import pandas as pd
import pytest

from tools.common.r_locator import find_rscript

CLI = __import__("pathlib").Path(__file__).resolve().parents[1] / "backend" / "additive_dominance_cli.R"
RSCRIPT = find_rscript()
pytestmark = pytest.mark.skipif(RSCRIPT is None, reason="R is not available")


def _write(path, X, ids):
    df = pd.DataFrame(X, columns=[f"M{i}" for i in range(X.shape[1])])
    df.insert(0, "ID", ids)
    df.to_csv(path, index=False)
    return str(path)


def _run(tmp_path, name, cfg):
    c, o = tmp_path / f"{name}.json", tmp_path / name
    c.write_text(json.dumps(cfg))
    proc = subprocess.run([str(RSCRIPT), str(CLI), str(c), str(o)], capture_output=True, text=True)
    res = o / "additive_dominance_results.json"
    if not res.exists():                       # a refused configuration answers on stdout only
        return json.loads(proc.stdout.strip().splitlines()[-1]), o
    return json.loads(res.read_text()), o


def _data(seed=3, n=100, m=600):
    rng = np.random.default_rng(seed)
    p = rng.uniform(0.1, 0.5, m)
    return rng.binomial(2, p, size=(n, m)).astype(float)


def test_split_additive_has_the_vanraden_scale_and_ignores_the_test_lines(tmp_path):
    X = _data()
    tr, te = np.arange(80), np.arange(80, 100)
    ftr = _write(tmp_path / "tr.csv", X[tr], [f"L{i}" for i in tr])
    fte = _write(tmp_path / "te.csv", X[te], [f"L{i}" for i in te])
    cfg = dict(train_file=ftr, test_file=fte, impute_method="mean", matrix_types=["additive", "dominance"])
    j, o = _run(tmp_path, "a", cfg)
    assert j["status"] == "success"
    A = pd.read_csv(o / "training_additive.csv", index_col=0).values
    p = X[tr].mean(0) / 2
    vr = (X[tr] - 2 * p) @ (X[tr] - 2 * p).T / (2 * (p * (1 - p)).sum())
    assert np.allclose(A, vr, atol=1e-5)                       # VanRaden with the training frequencies
    assert 0.9 < np.mean(np.diag(A)) < 1.1
    X2 = X.copy()
    X2[te] = np.random.default_rng(9).binomial(2, 0.5, size=(20, X.shape[1]))
    fte2 = _write(tmp_path / "te2.csv", X2[te], [f"L{i}" for i in te])
    j2, o2 = _run(tmp_path, "b", dict(cfg, test_file=fte2))
    for f in ("training_additive.csv", "training_dominance.csv"):
        assert np.array_equal(pd.read_csv(o / f).values, pd.read_csv(o2 / f).values), f


def test_vitezica_dominance_is_orthogonal_to_the_additive_effect_and_has_unit_scale(tmp_path):
    X = _data(seed=4, n=400, m=300)
    f = _write(tmp_path / "x.csv", X, [f"L{i}" for i in range(len(X))])
    j, o = _run(tmp_path, "v", dict(snp_file=f, impute_method="mean", matrix_types=["dominance"], dominance_method="vitezica"))
    assert j["status"] == "success" and "Vitezica" in j["summary"]["dominance_method"]
    D = pd.read_csv(o / "Dominance_Matrix.csv", index_col=0).values
    p = X.mean(0) / 2
    W = -2 * p**2 + (2 * p + 1) * X - X**2                   # dominance deviations
    assert np.allclose(D, W @ W.T / ((2 * p * (1 - p)) ** 2).sum(), atol=1e-5)
    assert 0.9 < np.mean(np.diag(D)) < 1.1
    r = [abs(np.corrcoef(W[:, k], X[:, k])[0, 1]) for k in range(X.shape[1])]
    assert np.median(r) < 0.05


def test_dominance_of_inbred_lines_gets_a_warning_and_an_unknown_method_is_refused(tmp_path):
    X = np.random.default_rng(5).binomial(1, 0.5, size=(60, 300)) * 2.0          # no heterozygotes
    f = _write(tmp_path / "ril.csv", X, [f"L{i}" for i in range(60)])
    j, _ = _run(tmp_path, "r", dict(snp_file=f, impute_method="mean", matrix_types=["additive", "dominance"]))
    assert j["status"] == "success" and any("heterozygous" in w for w in j["summary"]["warnings"])
    j, _ = _run(tmp_path, "bad", dict(snp_file=f, impute_method="mean", matrix_types=["dominance"], dominance_method="x"))
    assert j["status"] == "error"


def test_summary_has_no_empty_objects_for_an_additive_only_run(tmp_path):
    """jsonlite writes a NULL inside a list as {}, which the page then tried to render as text (React error 31)."""
    f = _write(tmp_path / "x.csv", _data(n=40, m=200), [f"L{i}" for i in range(40)])

    def walk(v, path=""):
        if isinstance(v, dict):
            assert v, f"empty object at {path}"
            for k, x in v.items():
                walk(x, f"{path}.{k}")
        elif isinstance(v, list):
            for i, x in enumerate(v):
                walk(x, f"{path}[{i}]")

    for types in (["additive"], ["dominance"], ["additive", "dominance"]):
        j, _ = _run(tmp_path, "a" + "".join(t[0] for t in types), dict(snp_file=f, impute_method="mean", matrix_types=types))
        assert j["status"] == "success"
        walk(j)
        assert isinstance(j["summary"]["warnings"], list)
        assert ("dominance_method" in j["summary"]) == ("dominance" in types)
