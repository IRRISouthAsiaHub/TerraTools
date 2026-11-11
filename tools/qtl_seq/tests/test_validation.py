import pytest

from tools.qtl_seq.backend import validation


def test_reference_values_are_expectations_not_results():
    ref = validation.REFERENCE
    q = ref["blast_ril"]["qtls"][0]
    assert (q["chromosome"], q["start"], q["end"], q["max_delta"]) == ("chr6", 2_390_000, 4_390_000, 0.79)
    assert ref["blast_ril"]["population"] == "RIL" and ref["blast_ril"]["n_high"] == 20 and ref["blast_ril"]["total_lines"] == 241
    f2 = ref["vigour_f2"]
    assert f2["population"] == "F2" and f2["n_high"] == 50 and f2["total_lines"] == 531
    assert [(x["chromosome"], x["max_delta"], x["level"]) for x in f2["qtls"]] == [("chr3", 0.61, 99), ("chr1", 0.67, 95)]


@pytest.mark.parametrize("key", ["blast_ril", "vigour_f2"])
def test_simulated_experiment_recovers_the_planted_qtls(key):
    r = validation.run_validation(key, seed=1)
    assert r["status"] == "success"
    expected = [c for c in r["checks"] if not c.get("info")]
    assert expected and all(c["found"] for c in expected), r["checks"]
    assert r["passed"], r["checks"]
    # the planted QTLs are positive: the index allele is enriched in the high bulk
    assert all(q["Peak_Delta_SNP_index"] > 0 for q in r["qtls"] if q["Chromosome"] in {e["chromosome"] for e in validation.REFERENCE[key]["qtls"]})


def test_comparison_flags_a_missing_or_weak_qtl():
    import pandas as pd
    empty = pd.DataFrame(columns=["QTL_ID", "Chromosome", "Start", "End", "Peak_Position", "Peak_Delta_SNP_index", "Significance_Level"])
    checks = validation.compare(empty, "blast_ril")
    assert not checks[0]["found"] and not checks[0]["passed"]
    weak = pd.DataFrame([{"QTL_ID": "QTL1", "Chromosome": "chr6", "Start": 2_000_000, "End": 5_000_000, "Peak_Position": 3_000_000,
                          "Peak_Delta_SNP_index": 0.5, "Significance_Level": "95%"}])
    c = validation.compare(weak, "blast_ril")[0]
    assert c["found"] and not c["passed"] and "95%" in c["note"] and "differs" in c["note"]
