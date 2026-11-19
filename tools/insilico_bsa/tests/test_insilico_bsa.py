import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "docs"))

from tools.insilico_bsa.backend import insilico_bsa as B    # noqa: E402


def _files(tmp_path, seed=2022):
    import make_insilico_demo as M
    markers, G, y, truth = M.simulate(seed=seed)
    names = [f"F2-{i + 1:03d}" for i in range(G.shape[1])]
    g = tmp_path / "geno.csv"
    pd.concat([markers, pd.DataFrame(G, columns=names)], axis=1).to_csv(g, index=False)
    p = tmp_path / "pheno.csv"
    pd.DataFrame({"Sample": names, "Height": y}).to_csv(p, index=False)
    return g, p, truth, names, y, markers, G


def _cfg(g, p, out, **kw):
    return {"input_file": str(g), "pheno_file": str(p), "trait": "Height", "output_dir": str(out), "permutations": 200, "population": "f2", **kw}


def test_g_statistic_matches_the_definition():
    # one marker: high pool 30 QQ / 20 Qq / 10 qq, low pool 10 / 20 / 30
    O = np.array([[[30.0], [20.0], [10.0]], [[10.0], [20.0], [30.0]]])
    g = float(B._g_table(O)[0])
    E = np.array([[20, 20, 20], [20, 20, 20]], dtype=float)
    ref = sum(o * np.log(o / e) for o, e in zip(O[..., 0].ravel(), E.ravel()))
    assert g == pytest.approx(ref)
    assert B._g_table(np.zeros((2, 3, 1)))[0] == 0                       # empty marker is not an error


def test_a_fully_over_dominant_locus_is_seen_by_gw_but_not_by_allele_statistics(tmp_path):
    g, p, truth, *_ = _files(tmp_path)
    a = B.analyse(_cfg(g, p, tmp_path / "o", permutations=300), log=lambda m: None)
    reg = a["regions"]
    chrs = set(reg["Chromosome"])
    assert chrs == {1, 2, 3}                                             # all three planted loci
    od = reg[reg["Chromosome"] == 2].iloc[0]
    q = [t for t in truth if t["chromosome"] == 2][0]["position_bp"]
    assert od["Start_bp"] <= q <= od["End_bp"]
    assert not od["Also_significant_by_G1"] and not od["Also_significant_by_dAF"]
    assert "over-dominant" in od["Pattern"]
    add = reg[reg["Chromosome"] == 1].iloc[0]
    assert add["Also_significant_by_G1"]                                 # additive / dominant loci are seen by both


def test_inbred_lines_make_gw_equal_g1(tmp_path):
    g, p, *_ , names, y, markers, G = _files(tmp_path)
    Gi = np.where(G == 1, -1, G)                                         # RIL-like: no heterozygous calls
    gi = tmp_path / "ril.csv"
    pd.concat([markers, pd.DataFrame(Gi, columns=names)], axis=1).to_csv(gi, index=False)
    a = B.analyse(_cfg(gi, p, tmp_path / "o", permutations=0, population="ril"), log=lambda m: None)
    mk = a["markers"]
    assert np.allclose(mk["Gw"], mk["G1"], atol=1e-9)
    assert not a["diploid"]


def test_the_null_gives_few_false_regions(tmp_path):
    g, p, _, names, y, markers, G = _files(tmp_path)
    rng = np.random.default_rng(5)
    pn = tmp_path / "null.csv"
    pd.DataFrame({"Sample": names, "Height": rng.normal(size=len(names))}).to_csv(pn, index=False)
    hits = 0
    for s in range(4):
        a = B.analyse(_cfg(g, pn, tmp_path / "o", permutations=200, seed=s), log=lambda m: None)
        hits += int(len(a["regions"]) > 0)
    assert hits <= 1                                                     # alpha 0.05 genome-wide: at most the odd false positive


def test_permutation_is_reproducible_and_pools_are_the_extremes(tmp_path):
    g, p, *_ = _files(tmp_path)
    a1 = B.analyse(_cfg(g, p, tmp_path / "o", permutations=100, seed=7), log=lambda m: None)
    a2 = B.analyse(_cfg(g, p, tmp_path / "o", permutations=100, seed=7), log=lambda m: None)
    assert a1["thr_use"] == a2["thr_use"]
    y = a1["y"]
    assert y[a1["lo_i"]].max() < np.median(y) < y[a1["hi_i"]].min()
    assert len(a1["lo_i"]) == len(a1["hi_i"]) == a1["pool_size"] == 60


def test_gwas_ready_layout_needs_no_phenotype_file(tmp_path):
    g, p, truth, names, y, markers, G = _files(tmp_path)
    df = pd.DataFrame(np.where(G == -1, np.nan, G).T, columns=[f"chr{c}_{pp}" for c, pp in zip(markers["CHROM"], markers["POS"])])
    df.insert(0, "Pheno", y)
    df.insert(0, "ID", names)
    f = tmp_path / "gwas_ready.csv"
    df.to_csv(f, index=False)
    info = B.inspect_inputs(str(f))
    assert info["traits"] == ["Pheno"] and info["suggest_population"] == "f2"
    a = B.analyse({"input_file": str(f), "output_dir": str(tmp_path / "o"), "permutations": 100}, log=lambda m: None)
    assert a["trait"] == "Pheno" and len(a["regions"]) >= 2


def test_problems_are_reported_in_words(tmp_path):
    g, p, *_ = _files(tmp_path)
    for cfg, msg in (({"trait": "Nope"}, "not a numeric column"), ({"pool_fraction": 0.9}, "between 5"), ({"permutations": 10 ** 6}, "permutations must be")):
        with pytest.raises(B.BsaError, match=msg):
            B.analyse(_cfg(g, p, tmp_path / "o", **cfg), log=lambda m: None)
    bad = tmp_path / "bad.csv"
    pd.DataFrame({"Sample": ["x1", "x2", "x3", "x4", "x5", "x6"], "Height": [1, 2, 3, 4, 5, 6]}).to_csv(bad, index=False)
    with pytest.raises(B.BsaError, match="match"):
        B.inspect_inputs(str(g), str(bad))
    assert B.run({"input_file": str(tmp_path / "missing.csv"), "output_dir": str(tmp_path / "o2")})["status"] == "error"


def test_run_writes_the_files_and_figures(tmp_path):
    g, p, truth, *_ = _files(tmp_path)
    r = B.run(_cfg(g, p, tmp_path / "out"), log=lambda m: None)
    assert r["status"] == "success" and r["counts"]["regions"] >= 3 and len(r["plots"]) == 7
    loci = pd.DataFrame(r["loci"])
    for t in truth:                                                  # every planted locus lies in a locus that all three window sizes found
        hit = loci[(loci["Chromosome"] == t["chromosome"]) & (loci["Start_bp"] <= t["position_bp"]) & (loci["End_bp"] >= t["position_bp"])]
        assert len(hit) == 1 and hit.iloc[0]["Windows_that_found_it"] == "3 of 3"
    assert [c["Window"] for c in r["comparison"]] == ["250kb", "500kb", "1Mb"]
    for k in ("workbook", "zip", "markers", "regions", "pools", "comparison", "loci", "windows"):
        assert Path(r["files"][k]).is_file()
    assert all(Path(p_["path"]).is_file() for p_ in r["plots"])
    sheets = pd.ExcelFile(r["files"]["workbook"]).sheet_names
    assert sheets == ["Summary", "Settings", "Window comparison", "Loci", "Regions", "SNP-index regions", "SNP-index track", "Thresholds", "Windows", "Pools", "Markers"]
    assert r["plots"][0]["key"] == "snpindex" and Path(r["files"]["snp_track"]).is_file() and Path(r["files"]["snp_regions"]).is_file()


# ───────────── explicit population, valid-marker smoothing, distance merging, GFF annotation ─────────────

def test_population_is_chosen_explicitly_and_checked_against_the_data(tmp_path):
    g, p, *_, names, y, markers, G = _files(tmp_path)
    a = B.analyse(_cfg(g, p, tmp_path / "o", permutations=0), log=lambda m: None)
    assert a["population"] == "f2" and not any("guessed" in n for n in a["notes"])
    auto = B.analyse(_cfg(g, p, tmp_path / "o", permutations=0, population="auto"), log=lambda m: None)
    assert any("Choose the population type yourself" in n for n in auto["notes"])      # the guess is announced, never silent
    ril = B.analyse(_cfg(g, p, tmp_path / "o", permutations=0, population="ril"), log=lambda m: None)
    assert ril["population"] == "ril" and not ril["diploid"] and any("treated as missing" in n for n in ril["notes"])
    bc = B.analyse(_cfg(g, p, tmp_path / "o", permutations=0, population="bc1"), log=lambda m: None)
    assert any("looks like an F2" in n for n in bc["notes"])                              # an F2 called a backcross is questioned
    with pytest.raises(B.BsaError, match="Unknown population"):
        B.analyse(_cfg(g, p, tmp_path / "o", population="tetraploid"), log=lambda m: None)


def test_backcross_data_are_recognised_and_analysed(tmp_path):
    g, p, truth, names, y, markers, G = _files(tmp_path)
    Gb = np.where(G == 2, 1, G)                                           # BC1: only 0 and 1 (heterozygous) remain
    gb = tmp_path / "bc.csv"
    pd.concat([markers, pd.DataFrame(Gb, columns=names)], axis=1).to_csv(gb, index=False)
    assert B.inspect_inputs(str(gb), str(p))["suggest_population"] == "bc1"
    a = B.analyse(_cfg(gb, p, tmp_path / "o", permutations=0, population="bc1"), log=lambda m: None)
    assert a["diploid"] and not any("looks like" in n for n in a["notes"])


def test_window_statistics_follow_the_maximum_statistic_formula_for_each_window_size(tmp_path):
    g, p, truth, names, y, markers, G = _files(tmp_path)
    a = B.analyse(_cfg(g, p, tmp_path / "o", permutations=200), log=lambda m: None)
    assert a["sizes"] == [250.0, 500.0, 1000.0] and list(a["windows"]) == ["250kb", "500kb", "1Mb"]
    assert not any("Smooth" in c for c in a["markers"].columns) and not hasattr(B, "smoothing_matrix") and not hasattr(B, "find_regions")
    mk = a["markers"]
    for lab, w in a["windows"].items():
        tab, mx = w["table"], w["maxima"]["Gw"]
        assert len(mx) == 200 and w["thr_use"]["Gw"] == pytest.approx(np.nanquantile(mx, 0.95))
        # the window value is the mean Gw of the scored markers it holds
        size = w["kb"] * 1000
        r0 = tab.iloc[len(tab) // 2]
        inw = mk[(mk["Chromosome"] == r0.Chromosome) & (mk["Position"] >= r0.Window_start) & (mk["Position"] < r0.Window_end) & mk["Scored"]]
        assert r0.Gw == pytest.approx(inw["Gw"].mean()) and r0.Scored_markers == len(inw)
        assert tab["Window_end"].sub(tab["Window_start"]).eq(size).all()
        # P_adj = (1 + #{T_b >= Gw}) / (B + 1)
        ok = tab[tab["Gw"].notna()]
        for i in (ok["Gw"].idxmax(), ok["Gw"].idxmin()):
            assert tab.loc[i, "P_adj_Gw"] == pytest.approx((1 + (mx >= tab.loc[i, "Gw"]).sum()) / 201)
        assert set(tab.index[tab["Sig_Gw"]]) == set(tab.index[tab["P_adj_Gw"] <= 0.05 + 1e-12])      # significant = adjusted p <= alpha


def test_windows_with_too_few_scored_markers_are_not_tested():
    pos = np.array([10_000, 20_000, 30_000, 600_000])                      # the second 250 kb window holds one marker only
    A, tab = B.build_windows(np.ones(4, dtype=int), pos, 250_000)
    assert list(tab["Window_start"]) == [0, 500_000] and A.shape == (2, 4)
    st = {"ok": np.array([[1], [1], [1], [1]], dtype=bool), "Gw": np.array([[6.0], [6.0], [0.0], [9.0]]), "G1": np.zeros((4, 1)), "dAF": np.zeros((4, 1))}
    st["Gw"][2] = 0.0
    means, valid = B.window_means(A, st, 3)
    assert valid[:, 0].tolist() == [True, False] and means["Gw"][0, 0] == pytest.approx(4.0) and np.isnan(means["Gw"][1, 0])
    st["ok"][2, 0] = False                                                  # an unscored marker is left out of the mean, it is not a zero
    means2, valid2 = B.window_means(A, st, 2)
    assert means2["Gw"][0, 0] == pytest.approx(6.0)


def test_markers_with_too_few_calls_get_no_statistic_and_the_locus_is_still_found(tmp_path):
    g, p, truth, names, y, markers, G = _files(tmp_path)
    Gm = G.copy()
    q = [t for t in truth if t["chromosome"] == 1][0]["position_bp"]
    near = np.flatnonzero((markers["CHROM"] == 1) & (abs(markers["POS"] - q) < 600_000))
    Gm[near[::2], :] = np.where(np.arange(Gm.shape[1]) % 12 == 0, Gm[near[::2], :], -1)     # every other marker around the QTL nearly empty
    gm = tmp_path / "gaps.csv"
    pd.concat([markers, pd.DataFrame(Gm, columns=names)], axis=1).to_csv(gm, index=False)
    a = B.analyse(_cfg(gm, p, tmp_path / "o", permutations=150, max_missing=0.97), log=lambda m: None)
    mk = a["markers"]
    assert (~mk["Scored"]).sum() > 0 and mk.loc[~mk["Scored"], "Gw"].isna().all()
    assert 1 in set(a["regions"]["Chromosome"])                          # the locus is still found
    assert any("no statistic" in n for n in a["notes"])


def test_significant_windows_form_regions_and_loci_are_compared_across_sizes():
    def tab(size, starts, sig):
        t = pd.DataFrame({"Chromosome": 1, "Window_start": starts, "Window_end": [x + size for x in starts]})
        t["Markers"] = 5
        t["Scored_markers"] = 5
        t["Gw"] = [9.0 if s_ else 1.0 for s_ in sig]
        t["G1"], t["abs_dAF"] = t["Gw"], 0.3
        for st_ in ("Gw", "G1", "dAF"):
            t["Sig_" + st_] = sig
        t["P_adj_Gw"] = 0.01
        return t
    pos = np.arange(0, 3_000_000, 50_000)
    mk = pd.DataFrame({"Marker": [f"m{i}" for i in range(len(pos))], "Chromosome": 1, "Position": pos, "Scored": True, "Gw": 5.0, "dAF": 0.2,
                       "AF2_low": .4, "AF2_high": .6, "Het_low": .5, "Het_high": .5})
    t5 = tab(500_000, [0, 500_000, 1_000_000, 2_000_000], [True, True, False, True])
    near = B.window_regions(mk, t5, "500kb", 500, merge_gap_bp=0)
    assert len(near) == 2 and sorted(near["Length_kb"]) == [500.0, 1000.0]          # two touching windows are one region
    merged = B.window_regions(mk, t5, "500kb", 500, merge_gap_bp=1_000_000)
    assert len(merged) == 1 and merged.iloc[0]["Start_bp"] == 0 and merged.iloc[0]["End_bp"] == 2_500_000
    t25 = tab(250_000, [0, 250_000, 750_000], [True, True, False])
    wins = {"250kb": {"kb": 250.0, "table": t25, "regions": B.window_regions(mk, t25, "250kb", 250), "thr_use": {"Gw": 5.0}},
            "500kb": {"kb": 500.0, "table": t5, "regions": near, "thr_use": {"Gw": 5.0}}}
    comp, loci = B.compare_windows(wins, [250.0, 500.0])
    assert list(comp["Regions"]) == [1, 2]
    assert len(loci) == 2 and loci.iloc[0]["Windows_that_found_it"] == "2 of 2" and loci.iloc[1]["Found_by"] == "500kb"     # chr1 0-1 Mb by both, 2-2.5 Mb only at 500 kb


def test_positions_come_from_the_pos_column_not_from_the_marker_name(tmp_path):
    rng = np.random.default_rng(3)
    n = 40
    rows = [[1, 1_000_000 + i * 1000, f"chr1_{5 + i * 7}_LOC_Os01g{i:05d}"] + list(rng.integers(0, 3, n)) for i in range(60)]
    df = pd.DataFrame(rows, columns=["CHROM", "POS", "Marker"] + [f"IRRI-{400 + i}" for i in range(n)])
    f = tmp_path / "tab.tsv"
    df.to_csv(f, sep="\t", index=False)
    d = B.read_genotypes(str(f))
    assert d["pos_source"] == "column 'POS'" and d["markers"]["POS"].iloc[0] == 1_000_000
    assert any("differs from the column 'POS'" in n for n in d["notes"])


def test_candidate_genes_from_a_gff_are_listed_per_region(tmp_path):
    from tools.insilico_bsa.backend import annotation as AN
    gff = tmp_path / "a.gff3"
    gff.write_text("##gff-version 3\n"
                   "Chr1\tMSU\tgene\t1000\t3000\t.\t+\t.\tID=LOC_Os01g00010;Name=OsA;Note=expressed%20protein\n"
                   "Chr1\tMSU\tmRNA\t1000\t3000\t.\t+\t.\tID=LOC_Os01g00010.1;Parent=LOC_Os01g00010\n"
                   "Chr1\tMSU\tgene\t9000\t9500\t.\t-\t.\tID=LOC_Os01g00020;Note=kinase\n"
                   "Chr1\tMSU\tgene\t50000\t52000\t.\t+\t.\tID=LOC_Os01g00030\n"
                   "Chr2\tMSU\tgene\t1000\t3000\t.\t+\t.\tID=LOC_Os02g00010\n")
    reg = pd.DataFrame({"Region": ["R1"], "Chromosome": [1], "Start_bp": [2000], "End_bp": [9100], "Peak_bp": [9300], "Peak_marker": ["m"]})
    cg = AN.candidate_genes(reg, AN.read_genes(str(gff)), flank_kb=0)
    assert list(cg["Gene_ID"]) == ["LOC_Os01g00020", "LOC_Os01g00010"]        # nearest to the peak first; the mRNA and chr2 are not genes here
    assert list(cg["Location"]) == ["inside the region"] * 2 and cg.iloc[1]["Description"] == "expressed protein"
    wide = AN.candidate_genes(reg, AN.read_genes(str(gff)), flank_kb=45)
    assert "LOC_Os01g00030" in set(wide["Gene_ID"]) and wide.set_index("Gene_ID").loc["LOC_Os01g00030", "Location"] == "in the flank"
    out = AN.run({"regions_file": _write(tmp_path / "r.csv", reg), "gff_file": str(gff), "flank_kb": 0, "output_dir": str(tmp_path / "ann")})
    assert out["status"] == "success" and out["genes_total"] == 2 and Path(out["files"]["xlsx"]).is_file()
    reg2 = pd.concat([reg.assign(Window="250kb"), reg.assign(Window="500kb", Region="R2")], ignore_index=True)
    rf = _write(tmp_path / "r2.csv", reg2)
    one = AN.run({"regions_file": rf, "window": "500kb", "gff_file": str(gff), "flank_kb": 0, "output_dir": str(tmp_path / "ann2")})
    assert one["status"] == "success" and one["genes_total"] == 2 and set(x["Region"] for x in one["genes"]) == {"R2"}
    assert AN.run({"regions_file": rf, "window": "1Mb", "gff_file": str(gff), "output_dir": str(tmp_path / "ann3")})["status"] == "error"
    assert AN.run({"regions_file": str(tmp_path / "none.csv"), "gff_file": str(gff), "output_dir": str(tmp_path)})["status"] == "error"


def _write(path, df):
    df.to_csv(path, index=False)
    return str(path)


def test_the_result_is_valid_json_without_nan(tmp_path):
    import json
    g, p, *_ = _files(tmp_path)
    r = B.run(_cfg(g, p, tmp_path / "out", merge_gap_kb=0), log=lambda m: None)
    assert r["status"] == "success"
    txt = json.dumps(B.json_safe(r), allow_nan=False)          # raises on NaN / inf
    assert "NaN" not in txt and "Infinity" not in txt
    assert any(c["Median_region_length_kb"] is None or isinstance(c["Median_region_length_kb"], float) for c in json.loads(txt)["comparison"])
    assert B.json_safe({"a": float("nan"), "b": [np.float64("inf"), np.int64(3), np.bool_(True)]}) == {"a": None, "b": [None, 3, True]}
    df = pd.DataFrame({"x": [1.5, np.nan], "y": ["a", None]})
    assert json.dumps(B._records(df), allow_nan=False) == '[{"x": 1.5, "y": "a"}, {"x": null, "y": null}]'


def test_a_window_size_without_significant_regions_lists_the_closest_candidates(tmp_path):
    g, p, truth, names, y, markers, G = _files(tmp_path)
    p0 = tmp_path / "null.csv"
    pd.DataFrame({"Line": names, "Null": np.random.default_rng(1).normal(size=len(names))}).to_csv(p0, index=False)    # no QTL at all
    a = B.analyse(_cfg(g, p0, tmp_path / "o", trait="Null", permutations=300), log=lambda m: None)
    assert len(a["regions"]) == 0 and len(a["near"]) > 0
    n = a["near"]
    need = {"Left_marker", "Left_bp", "Peak_marker", "Peak_bp", "Right_marker", "Right_bp", "P_adj_Gw_at_peak", "LOD", "Delta_SNP_index", "Fraction_of_threshold"}
    assert need <= set(n.columns) and (n["Fraction_of_threshold"] < 1.0).all() and (n["P_adj_Gw_at_peak"] > 0.05).all()
    assert set(n["Window"]) == {"250kb", "500kb", "1Mb"} and n["Region"].str.contains("-N").all()
    r = B.run(_cfg(g, p0, tmp_path / "o2", trait="Null", permutations=300), log=lambda m: None)
    assert r["status"] == "success" and len(r["near_regions"]) > 0 and Path(r["files"]["near"]).is_file()
    assert "Closest to significant" in pd.ExcelFile(r["files"]["workbook"]).sheet_names


# ───────────── SNP-index view: signed Delta, the bulk each value comes from, the allele of the chosen parent, permutation intervals ─────────────

def test_snp_index_of_each_pool_and_the_signed_delta_follow_the_chosen_parent(tmp_path):
    g, p, truth, *_ = _files(tmp_path)
    a2 = B.analyse(_cfg(g, p, tmp_path / "o2", permutations=60), log=lambda m: None)
    a1 = B.analyse(_cfg(g, p, tmp_path / "o1", permutations=60, snp_index_parent=1), log=lambda m: None)
    m2, m1 = a2["markers"], a1["markers"]
    assert np.allclose(m2["Delta_SNP_index"], m2["SNP_index_high"] - m2["SNP_index_low"], equal_nan=True)     # the difference is high pool minus low pool
    assert np.allclose(m2["Delta_SNP_index"], m2["AF2_high"] - m2["AF2_low"], equal_nan=True)                # parent 2: the allele frequency of code 2
    assert np.allclose(m1["SNP_index_high"], 1 - m2["SNP_index_high"], equal_nan=True)                        # parent 1 is the complement
    assert np.allclose(m1["Delta_SNP_index"], -m2["Delta_SNP_index"], equal_nan=True)
    r2, r1 = a2["snp_regions"], a1["snp_regions"]
    assert len(r2) >= 1 and len(r1) == len(r2)
    assert (np.sign(r1["Delta_SNP_index_peak"].to_numpy()) == -np.sign(r2["Delta_SNP_index_peak"].to_numpy())).all()
    with pytest.raises(B.BsaError):
        B.analyse(_cfg(g, p, tmp_path / "o3", snp_index_parent=3), log=lambda m: None)


def test_delta_snp_index_intervals_come_from_the_permutations_and_regions_lie_outside_them(tmp_path):
    g, p, truth, *_ = _files(tmp_path)
    a = B.analyse(_cfg(g, p, tmp_path / "o", permutations=300), log=lambda m: None)
    tr, regs = a["snp_track"], a["snp_regions"]
    ok = tr["CI99_upper"].notna()
    assert ok.any() and a["snp_null_perms"] == 300
    assert (tr.loc[ok, "CI99_upper"] >= tr.loc[ok, "CI95_upper"]).all() and (tr.loc[ok, "CI99_lower"] <= tr.loc[ok, "CI95_lower"]).all()
    assert (tr.loc[ok, "CI95_upper"] > 0).all() and (tr.loc[ok, "CI95_lower"] < 0).all()                       # the null is centred on zero
    assert len(regs) >= 2
    for t in truth:                                                                                           # the additive and dominant loci are peak regions
        if t["chromosome"] in (1, 3):
            assert ((regs["Chromosome"] == t["chromosome"]) & (regs["Left_bp"] <= t["position_bp"]) & (regs["Right_bp"] >= t["position_bp"])).any()
    for r in regs.itertuples():                                                                              # a peak region lies outside the 99 % interval at its peak
        assert abs(r.Delta_SNP_index_peak) > abs(r.CI99_at_peak) - 1e-9
    od = [t for t in truth if t["chromosome"] == 2][0]                                                       # the over-dominant locus has no allele-frequency peak
    assert not ((regs["Chromosome"] == 2) & (regs["Left_bp"] <= od["position_bp"]) & (regs["Right_bp"] >= od["position_bp"])).any()


def test_without_permutations_there_are_no_intervals_and_no_peak_regions(tmp_path):
    g, p, *_ = _files(tmp_path)
    a = B.analyse(_cfg(g, p, tmp_path / "o", permutations=0), log=lambda m: None)
    assert a["snp_track"]["CI99_upper"].isna().all() and len(a["snp_regions"]) == 0


def test_a_sparse_map_gets_a_wider_snp_index_window_and_a_dense_map_keeps_the_primary_one(tmp_path):
    g, p, truth, names, y, markers, G = _files(tmp_path)
    dense = B.analyse(_cfg(g, p, tmp_path / "o1", permutations=40, primary_window_kb=500), log=lambda m: None)
    assert dense["snp_window_auto"] is False and dense["snp_window_kb"] == 500
    keep = np.arange(0, len(markers), 6)                      # a thin, gene-like map
    g2 = tmp_path / "sparse.csv"
    pd.concat([markers.iloc[keep].reset_index(drop=True), pd.DataFrame(G[keep], columns=names)], axis=1).to_csv(g2, index=False)
    sparse = B.analyse(_cfg(g2, p, tmp_path / "o2", permutations=40, primary_window_kb=500, min_maf=0.0), log=lambda m: None)
    assert sparse["snp_window_auto"] is True and sparse["snp_window_kb"] > 500
    assert any("sparse" in n for n in sparse["notes"])
