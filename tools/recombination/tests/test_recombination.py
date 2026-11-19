import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "docs"))

from tools.recombination.backend import recombination as R    # noqa: E402


def _runs(v):
    return [(int(a), n) for a, n in R._runs(list(v))]


def test_legacy_smoothing_first_short_run_takes_the_previous_value():
    ab = np.array([1] * 10 + [2] * 2 + [1] * 10 + [2] * 10, dtype=np.int8)
    assert _runs(R.smooth_legacy(ab, 5)) == [(1, 22), (2, 10)]
    # the first run takes the value of the run after it
    assert _runs(R.smooth_legacy(np.array([1] * 2 + [2] * 10, dtype=np.int8), 5)) == [(2, 12)]


def test_robust_smoothing_keeps_long_blocks_and_merges_the_shortest_first():
    ab = np.array([1] * 20 + [2] * 3 + [1] * 2 + [2] * 20, dtype=np.int8)
    out = R.smooth_robust(ab, 5)
    assert _runs(out) == [(1, 20), (2, 25)] or _runs(out) == [(1, 25), (2, 20)]
    assert all(n >= 5 for _, n in _runs(out)) and len(out) == len(ab)
    assert np.array_equal(R.smooth_robust(np.array([1, 1, 1, 1, 1, 1], dtype=np.int8), 3), np.ones(6, dtype=np.int8))


def _sim_file(tmp_path, coded=False):
    import make_recombination_demo as M
    markers, G, truth = M.simulate()
    df = markers.copy()
    if coded:
        out = pd.DataFrame({"Marker": markers["ID"], "Chromosome": markers["CHROM"]})
        for j in range(G.shape[1]):
            out[f"RIL-{j + 1:03d}"] = G[:, j]
        f = tmp_path / "coded.csv"
        out.to_csv(f, index=False)
        return f, truth
    df["REF"], df["ALT"] = "A", "G"
    df["IRRI-PARENT_GT"] = 0
    for j in range(G.shape[1]):
        df[f"RIL-{j + 1:03d}_GT"] = G[:, j]
    f = tmp_path / "pop.csv"
    df.to_csv(f, index=False)
    return f, truth


def test_crossovers_are_recovered_from_a_simulated_population(tmp_path):
    f, truth = _sim_file(tmp_path)
    a = R.analyse({"input_file": str(f), "parent": "IRRI-PARENT", "smoothing": "legacy", "min_block_pct": 2}, log=lambda m: None)
    per = a["per_ril"].set_index("RIL")["Crossovers"].reindex([f"RIL-{j + 1:03d}" for j in range(len(truth))]).to_numpy()
    assert np.corrcoef(per, truth)[0, 1] > 0.85
    assert 0.7 * truth.sum() < per.sum() < 1.2 * truth.sum()
    # without smoothing every genotyping error is a crossover
    raw = R.analyse({"input_file": str(f), "parent": "IRRI-PARENT", "smoothing": "none"}, log=lambda m: None)
    assert len(raw["co"]) > 3 * len(a["co"])
    assert a["has_genes"] and set(a["co"]["CO_type"]) <= {"intragenic", "intergenic"}


def test_a_coded_table_gives_the_same_result_as_the_parent_column(tmp_path):
    f, _ = _sim_file(tmp_path)
    fc, _ = _sim_file(tmp_path, coded=True)
    a = R.analyse({"input_file": str(f), "parent": "IRRI-PARENT", "min_block_pct": 2}, log=lambda m: None)
    b = R.analyse({"input_file": str(fc), "parent": "", "min_block_pct": 2}, log=lambda m: None)
    assert len(a["co"]) == len(b["co"]) and a["n_poly"] == b["n_poly"]
    assert np.allclose(a["per_ril"]["Donor_pct"], b["per_ril"]["Donor_pct"])


def test_input_problems_are_reported_in_words(tmp_path):
    f, _ = _sim_file(tmp_path)
    for cfg, msg in (({"parent": "NOPE"}, "not a sample"), ({"parent": "IRRI-PARENT", "donor": "IRRI-PARENT"}, "different"),
                     ({"parent": "IRRI-PARENT", "smoothing": "x"}, "Unknown smoothing"), ({"parent": "IRRI-PARENT", "bin_mb": 0}, "bin size")):
        with pytest.raises(R.RecombError, match=msg):
            R.analyse({"input_file": str(f), **cfg})
    bad = tmp_path / "bad.csv"
    bad.write_text("Marker,a,b,c\nm1,0,1,2\nm2,0,1,2\nm3,0,1,2\n")
    with pytest.raises(R.RecombError, match="chromosome"):
        R.read_genotypes(bad)


def test_run_writes_the_files_and_figures(tmp_path):
    f, _ = _sim_file(tmp_path)
    r = R.run({"input_file": str(f), "parent": "IRRI-PARENT", "min_block_pct": 2, "output_dir": str(tmp_path / "out")})
    assert r["status"] == "success" and len(r["plots"]) == 7
    for k in ("workbook", "zip", "crossovers", "per_ril", "per_chr", "bins"):
        assert Path(r["files"][k]).is_file()
    assert all(Path(p["path"]).is_file() for p in r["plots"])
    assert R.run({"input_file": str(tmp_path / "missing.csv"), "output_dir": str(tmp_path / "o2")})["status"] == "error"
