"""Clustering: K-means on the principal components finds the populations, the K-selection aid points at the right K,
and the split is stratified by cluster."""
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from tools.clustering.backend.clustering_engine import cluster_split_data


def _three_populations(tmp_path):
    rng = np.random.default_rng(7)
    m, fst = 1500, 0.1
    anc = rng.uniform(0.1, 0.9, m)
    X, lab = [], []
    for k in range(3):
        pk = rng.beta(anc * (1 - fst) / fst, (1 - anc) * (1 - fst) / fst)
        X.append(rng.binomial(2, pk, size=(60, m)))
        lab += [k] * 60
    df = pd.DataFrame(np.vstack(X), columns=[f"M{i}" for i in range(m)])
    df.insert(0, "ID", [f"L{i}" for i in range(180)])
    f = tmp_path / "g.csv"
    df.to_csv(f, index=False)
    return str(f), np.array(lab)


def test_clusters_recover_the_populations_with_all_or_few_components(tmp_path):
    f, lab = _three_populations(tmp_path)
    for npc in (0, 5):
        r = cluster_split_data(f, 3, 0.2, tmp_path / f"o{npc}", n_pcs=npc)
        cl = pd.read_csv(r["cluster_table_csv"]).set_index("ID").loc[[f"L{i}" for i in range(180)], "Cluster"].values
        assert adjusted_rand_score(lab, cl) > 0.99
        assert r["n_pcs_used"] == (r["n_pcs_available"] if npc == 0 else 5)


def test_k_selection_suggests_the_true_number_of_populations(tmp_path):
    f, _ = _three_populations(tmp_path)
    r = cluster_split_data(f, 5, 0.2, tmp_path / "o")
    assert r["k_suggested"] == 3
    assert [row["K"] for row in r["k_selection"]] == list(range(2, 11))
    assert r["k_selection_png"].endswith("k_selection.png")


def test_split_is_stratified_and_tiny_files_still_run(tmp_path):
    f, lab = _three_populations(tmp_path)
    r = cluster_split_data(f, 3, 0.2, tmp_path / "o")
    t = pd.read_csv(r["cluster_table_csv"])
    assert (t.groupby("Cluster")["Set"].apply(lambda s: (s == "Test").sum()) == 12).all()
    small = pd.read_csv(f).iloc[:6, :40]
    s = tmp_path / "s.csv"
    small.to_csv(s, index=False)
    assert cluster_split_data(str(s), 2, 0.3, tmp_path / "o2")["status"] == "success"
