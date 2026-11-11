import numpy as np
import pandas as pd
import pytest

from tools.qtl_seq.backend import stats
from tools.qtl_seq.backend.calling import homozygous_call


def snp_table(rows):
    return pd.DataFrame(rows, columns=["Chromosome", "Position", "Ref", "Alt", "Ref_High", "Alt_High", "Ref_Low", "Alt_Low"])


# ---------------------------------------------------------------- 1. SNP-index
def test_snp_index_is_alt_over_total_reads():
    assert stats.snp_index(3, 7) == pytest.approx(0.3)
    assert stats.snp_index(0, 10) == 0.0
    assert stats.snp_index(10, 0) == 1.0
    assert np.isnan(stats.snp_index(0, 0))
    got = stats.snp_index(np.array([1, 2, 0]), np.array([3, 2, 0]))
    assert got[:2] == pytest.approx([0.25, 0.5]) and np.isnan(got[2])


# ---------------------------------------------------------------- 2. Delta SNP-index
def test_delta_sign_convention_and_range():
    t = stats.add_indices(snp_table([["c", 1, "A", "G", 2, 18, 16, 4], ["c", 2, "A", "G", 20, 0, 0, 20], ["c", 3, "A", "G", 0, 20, 20, 0]]))
    assert t["SNP_index_High"][0] == pytest.approx(0.9) and t["SNP_index_Low"][0] == pytest.approx(0.2)
    assert t["Delta_SNP_index"][0] == pytest.approx(0.7)          # index allele enriched in the HIGH bulk -> positive
    assert t["Delta_SNP_index"][1] == pytest.approx(-1.0)         # enriched in the LOW bulk -> negative
    assert t["Delta_SNP_index"][2] == pytest.approx(1.0)
    assert t["Delta_SNP_index"].between(-1, 1).all()
    assert stats.delta_snp_index(0.4, 0.4) == 0


# ---------------------------------------------------------------- 3. parental genotype filtering
def test_homozygous_call_confidence_depends_on_counts():
    counts = np.array([[20, 18, 10, 0], [0, 2, 10, 20], [0, 0, 0, 0], [0, 0, 0, 0]])
    top, depth, gq = homozygous_call(counts)
    assert list(top[:2]) == [0, 0] and top[3] == 1
    assert gq[0] >= 30 and gq[3] >= 30                  # clean homozygotes
    assert gq[1] < 30                                   # 2 reads of a second allele: not confident
    assert gq[2] < 5                                    # 10/10: heterozygous


def test_parent_genotypes_use_gt_gq_and_depth():
    from tools.qtl_seq.backend.vcfio import parent_genotypes, index_is_alt
    df = pd.DataFrame({"Chromosome": "c", "Position": [1, 2, 3, 4, 5], "Ref": "A", "Alt": "G",
                       "AD_Ref": [20, 0, 10, 3, np.nan], "AD_Alt": [0, 20, 10, 3, np.nan], "DP": [20, 20, 20, 6, 20],
                       "GT": ["0/0", "1/1", "0/1", "1/1", "1/1"], "GQ": [60, 60, 50, 60, 10]})
    g = parent_genotypes(df, min_depth=7, min_gq=30)
    assert list(g["Class"]) == ["ref", "alt", "", "alt", "alt"]
    assert list(g["Pass"]) == [True, True, False, False, False]      # het, depth 6 and GQ 10 are removed
    assert list(index_is_alt("parent2", g["Class"], "parent2")) == [False, True, False, True, True]
    assert list(index_is_alt("parent1", g["Class"], "parent2")) == [True, False, True, False, False]


# ---------------------------------------------------------------- 4. read-depth and SNP-index filtering
def test_filters_follow_the_documented_rules_and_are_counted():
    t = snp_table([
        ["c", 1, "A", "G", 3, 3, 10, 10],        # depth 6 in the high bulk -> removed (low depth)
        ["c", 2, "A", "G", 10, 10, 10, 10],      # fine
        ["c", 3, "A", "G", 20, 1, 20, 1],        # index < 0.3 in both bulks -> removed
        ["c", 4, "A", "G", 14, 6, 20, 0],        # high 0.3 -> kept, low index 0 is kept as it is
        ["c", 5, "A", "G", 70, 70, 10, 10],      # depth 140 -> removed only when max_depth is set
    ])
    kept, led = stats.apply_filters(t, min_depth=7, max_depth=0, min_snp_index=0.3)
    assert list(kept["Position"]) == [2, 4, 5]
    assert led == {"input": 5, "removed_low_depth": 1, "removed_high_depth": 0, "removed_snp_index": 1, "retained": 3}
    assert kept[kept["Position"] == 4]["SNP_index_Low"].iloc[0] == 0.0
    kept, led = stats.apply_filters(t, min_depth=7, max_depth=100, min_snp_index=0.3)
    assert list(kept["Position"]) == [2, 4] and led["removed_high_depth"] == 1
    kept, led = stats.apply_filters(t, min_depth=7, use_snp_index=False)
    assert 3 in list(kept["Position"])


# ---------------------------------------------------------------- 5. sliding windows
def test_sliding_window_means_counts_and_minimum():
    pos = np.array([100, 200, 300, 1500, 1600])
    w = stats.sliding_windows(pos, {"v": np.array([1.0, 2.0, 3.0, 10.0, 20.0])}, chrom_length=2000, size=1000, step=500, min_snps=2)
    # windows start 1, 501, 1001 (1001 + 999 = 2000 fits); the one at 501 holds no SNP of 501..1500 except 1500
    assert list(w["Window_Start"]) == [1, 1001]
    assert list(w["Window_End"]) == [1000, 2000]
    assert list(w["SNP_Count"]) == [3, 2]
    assert w["v"] == pytest.approx([2.0, 15.0])


def test_window_grid_and_short_chromosomes():
    assert list(stats.window_starts(2_000_000, 2_000_000, 10_000)) == [1]
    assert list(stats.window_starts(1_000_000, 2_000_000, 10_000)) == [1]
    s = stats.window_starts(2_030_000, 2_000_000, 10_000)
    assert s[0] == 1 and s[-1] == 30_001 and len(s) == 4


# ---------------------------------------------------------------- 6. confidence intervals
def test_limits_are_reproducible_ordered_and_symmetric():
    a = stats.simulate_limits([12], [12], "RIL", 20, 20, 10000, seed=7).frame
    b = stats.simulate_limits([12], [12], "RIL", 20, 20, 10000, seed=7).frame
    c = stats.simulate_limits([12], [12], "RIL", 20, 20, 10000, seed=8).frame
    pd.testing.assert_frame_equal(a, b)
    assert not a.equals(c)
    r = a.iloc[0]
    assert r["CI99_Upper"] > r["CI95_Upper"] > r["CI90_Upper"] > 0
    assert r["CI99_Lower"] < r["CI95_Lower"] < r["CI90_Lower"] < 0
    assert r["CI95_Upper"] == pytest.approx(-r["CI95_Lower"], abs=0.08)


def test_limits_narrow_with_depth_and_plants_and_follow_the_population():
    f = stats.simulate_limits([7, 30, 100], [7, 30, 100], "RIL", 20, 20, 10000, seed=1).frame
    up = list(f.sort_values("Depth_High")["CI95_Upper"])
    assert up[0] > up[1] > up[2]
    ril = stats.simulate_limits([30], [30], "RIL", 20, 20, 10000, seed=1).frame["CI95_Upper"][0]
    f2 = stats.simulate_limits([30], [30], "F2", 20, 20, 10000, seed=1).frame["CI95_Upper"][0]
    big = stats.simulate_limits([30], [30], "RIL", 100, 100, 10000, seed=1).frame["CI95_Upper"][0]
    assert f2 < ril                    # an F2 plant contributes two alleles
    assert big < ril


def test_null_matches_the_analytic_variance():
    rng = np.random.default_rng(3)
    sim = stats.simulate_null(30, 30, "RIL", 20, 20, 200_000, rng)
    # Var = 2 * [ p(1-p)/n + p(1-p)/d * (1 - 1/n) ]  with p = 1/2
    n, d = 20, 30
    expected = 2 * (0.25 / n + 0.25 / d * (1 - 1 / n))
    assert sim.var() == pytest.approx(expected, rel=0.03)
    assert sim.mean() == pytest.approx(0, abs=0.005)


def test_backcross_null_is_centred_on_a_quarter():
    rng = np.random.default_rng(1)
    assert stats.bulk_frequency("BC", 40, 50000, rng).mean() == pytest.approx(0.25, abs=0.005)
    assert stats.bulk_frequency("BC", 40, 50000, rng, "index").mean() == pytest.approx(0.75, abs=0.005)


def test_unknown_population_is_rejected():
    with pytest.raises(ValueError):
        stats.bulk_frequency("DH", 20, 10, np.random.default_rng(0))


# ---------------------------------------------------------------- 7. significant windows
def make_windows(deltas, up95=0.5, level_up=(0.4, 0.5, 0.7), start=1, step=10_000, size=100_000, chrom="c"):
    n = len(deltas)
    w = pd.DataFrame({"Chromosome": chrom, "Window_Start": start + step * np.arange(n), "Window_End": start + step * np.arange(n) + size - 1,
                      "SNP_Count": 20, "Mean_SNP_index_High": 0.5, "Mean_SNP_index_Low": 0.5, "Mean_Delta_SNP_index": deltas})
    for lv, u in zip(stats.LEVELS, level_up):
        w[f"CI{lv}_Upper"], w[f"CI{lv}_Lower"] = u, -u
    return stats.flag_windows(w)


def test_windows_are_flagged_at_the_right_level_in_both_directions():
    w = make_windows([0.0, 0.45, 0.55, 0.75, -0.45, -0.6, -0.8])
    assert list(w["Significant_90"]) == [False, True, True, True, True, True, True]
    assert list(w["Significant_95"]) == [False, False, True, True, False, True, True]
    assert list(w["Significant_99"]) == [False, False, False, True, False, False, True]
    sw = stats.significant_windows(w, 95)
    assert list(sw["Level"]) == [95, 99, 95, 99] and list(sw["Sign"]) == [1, 1, -1, -1]


def test_a_delta_of_one_is_not_significant_by_itself():
    w = make_windows([0.99, 0.8], level_up=(1.2, 1.3, 1.5))
    assert not w[["Significant_90", "Significant_95", "Significant_99"]].to_numpy().any()


# ---------------------------------------------------------------- 8. QTL interval merging
def test_adjacent_windows_merge_into_one_interval_with_peak_inside():
    d = [0] * 3 + [0.45, 0.6, 0.8, 0.65, 0.5] + [0] * 3
    w = make_windows(d)
    q, sig = stats.call_qtls(w, {"c": np.arange(1, 200_000, 500)}, window_step=10_000, min_level=90, merge_gap=0, basis="window_extent")
    assert len(q) == 1
    r = q.iloc[0]
    assert r["Significance_Level"] == "99%" and r["Peak_Delta_SNP_index"] == pytest.approx(0.8)
    assert r["Start"] == int(w["Window_Start"][3]) and r["End"] == int(w["Window_End"][7])
    assert r["Start"] < r["Peak_Position"] < r["End"]
    assert r["Window_Count"] == 5 and r["Confidence_Threshold"] == pytest.approx(0.7)
    assert r["Minimum_Delta_SNP_index"] == pytest.approx(0.45) and r["Maximum_Delta_SNP_index"] == pytest.approx(0.8)
    assert r["SNP_Count"] == int(((np.arange(1, 200_000, 500) >= r["Start"]) & (np.arange(1, 200_000, 500) <= r["End"])).sum())
    assert list(sig["QTL_ID"]) == ["QTL1"] * 5
    # the strictest-limit basis narrows the interval to the 99 % stretch
    q2, _ = stats.call_qtls(w, {"c": np.arange(1, 200_000, 500)}, 10_000, 90, 0, "top_level")
    assert q2.iloc[0]["Start"] == int(w["Window_Start"][5]) and q2.iloc[0]["End"] == int(w["Window_End"][5])
    q3, _ = stats.call_qtls(w, {"c": np.arange(1, 200_000, 500)}, 10_000, 90, 0, "window_midpoints")
    assert q3.iloc[0]["Start"] == int(w["Window_Mid"][3]) and q3.iloc[0]["End"] == int(w["Window_Mid"][7])


def test_gaps_signs_and_chromosomes_split_intervals():
    d = [0.6, 0.6, 0.0, 0.0, 0.6, 0.6, -0.6, -0.6]
    w = make_windows(d)
    q, _ = stats.call_qtls(w, {}, 10_000, 90, merge_gap=0, basis="window_extent")
    assert len(q) == 3 and list(q["Peak_Delta_SNP_index"] > 0) == [True, True, False]
    assert "negative" in q.iloc[2]["Direction"] and "positive" in q.iloc[0]["Direction"]
    # a gap of two missing windows (20 kb) is bridged by merge_gap = 30 kb but not by 0; the sign change always splits
    q, _ = stats.call_qtls(w, {}, 10_000, 90, merge_gap=30_000, basis="window_extent")
    assert len(q) == 2
    other = make_windows([0.6, 0.6], chrom="d")
    q, _ = stats.call_qtls(pd.concat([make_windows([0.6, 0.6]), other], ignore_index=True), {}, 10_000, 90, 0)
    assert list(q["Chromosome"]) == ["c", "d"]


def test_no_significant_window_gives_no_qtl_and_level_counts():
    q, sig = stats.call_qtls(make_windows([0.1, 0.0, -0.2]), {}, 10_000, 90, 0)
    assert q.empty and sig.empty
    w = make_windows([0.8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.55, 0.55])
    q, _ = stats.call_qtls(w, {}, 10_000, 90, 0)
    assert stats.level_counts(q) == {"90": 2, "95": 2, "99": 1}
    q, _ = stats.call_qtls(w, {}, 10_000, 99, 0)
    assert len(q) == 1


# ---------------------------------------------------------------- genome-wide windows with limits
def test_build_windows_end_to_end_on_a_planted_signal():
    rng = np.random.default_rng(5)
    pos = np.sort(rng.choice(np.arange(1, 4_000_000), 3000, replace=False))
    p_h = np.where((pos > 1_500_000) & (pos < 2_500_000), 0.9, 0.5)
    ah, al = rng.binomial(20, p_h), rng.binomial(20, np.where(p_h > 0.5, 0.1, 0.5))
    t = snp_table([["c", int(x), "A", "G", 20 - a, a, 20 - b, b] for x, a, b in zip(pos, ah, al)])
    kept, _ = stats.apply_filters(t, 7, 0, 0.3)
    lim = stats.simulate_limits(kept["Depth_High"], kept["Depth_Low"], "RIL", 20, 20, 5000, seed=1)
    for c, a in lim.limits_for(kept["Depth_High"], kept["Depth_Low"]).items():
        kept[c] = a
    p = {"window_size": 1_000_000, "window_step": 50_000, "min_snps_per_window": 10, "ci_method": "snp_average"}
    win = stats.build_windows(kept, {"c": 4_000_000}, p, lim)
    q, _ = stats.call_qtls(win, {"c": kept["Position"].to_numpy()}, 50_000, 95, 100_000)
    assert len(q) == 1 and 1_700_000 < q.iloc[0]["Peak_Position"] < 2_300_000
    p["ci_method"] = "window"
    p.update(population="RIL", n_high=20, n_low=20, replicates=5000, seed=1)
    win2 = stats.build_windows(kept, {"c": 4_000_000}, p, None)
    assert (win2["CI95_Upper"] < win["CI95_Upper"]).all()          # more SNPs narrow the window-level limits
    assert win2["Significant_95"].any()
