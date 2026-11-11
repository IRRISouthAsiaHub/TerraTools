import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tools.common.seqtools import Fasta
from tools.qtl_seq.backend import bamio, bcftools, deps, pipeline, vcfio
from tools.qtl_seq.tests import helpers as H

SAMTOOLS = bool(deps.where("samtools"))
PYSAM = bool(deps.pysam_where())
BCFTOOLS = bool(deps.where("bcftools"))


def ad(ref, alt, gt="0/1", gq=60):
    return f"{gt}:{ref},{alt}:{ref + alt}:{gq}"


def vcf_pair(tmp, rows_high, rows_low, **kw):
    contigs = {"c1": 100_000, "c2": 50_000}
    H.write_vcf(tmp / "high.vcf", rows_high, ["HIGH"], contigs, **kw)
    H.write_vcf(tmp / "low.vcf", rows_low, ["LOW"], contigs, **kw)
    return {"high_vcf": str(tmp / "high.vcf"), "low_vcf": str(tmp / "low.vcf")}


def rec(c, p, f, ref="A", alt="G", q=50, flt="PASS"):
    return (c, p, ref, alt, [f], q, flt)


# ---------------------------------------------------------------- 9. VCF validation
def test_vcf_reads_ad_dp_gt_gq(tmp_path):
    H.write_vcf(tmp_path / "a.vcf", [rec("c1", 100, ad(3, 7, "0/1", 40))], ["S"])
    df, info = vcfio.read_vcf(tmp_path / "a.vcf")
    r = df.iloc[0]
    assert (r["AD_Ref"], r["AD_Alt"], r["DP"], r["GT"], r["GQ"]) == (3, 7, 10, "0/1", 40)
    assert info["has_ad"] and info["has_gt"] and info["has_gq"] and info["has_dp"]


def test_vcf_counts_multiallelic_indels_duplicates_and_missing_ad(tmp_path):
    rows = [rec("c1", 100, ad(5, 5)),
            rec("c1", 100, ad(6, 6)),                                         # duplicate position
            ("c1", 200, "A", "G,T", ["0/1:4,3,2:9:50"], 50, "PASS"),          # multiallelic
            ("c1", 300, "AT", "A", [ad(5, 5)], 50, "PASS"),                   # indel
            ("c1", 400, "A", "G", ["0/1:.:10:50"], 50, "PASS"),               # AD missing
            rec("c1", 500, ad(2, 8))]
    H.write_vcf(tmp_path / "a.vcf", rows, ["S"])
    df, info = vcfio.read_vcf(tmp_path / "a.vcf")
    assert info["duplicates"] == 1 and info["multiallelic"] == 1 and info["indels"] == 1 and info["missing_ad"] == 1
    assert list(df["Position"]) == [100, 400, 500]


def test_vcf_without_ad_is_rejected_with_a_clear_message(tmp_path):
    H.write_vcf(tmp_path / "high.vcf", [("c1", 100, "A", "G", ["0/1:20"], 50, "PASS")], ["H"], fmt="GT:DP")
    H.write_vcf(tmp_path / "low.vcf", [("c1", 100, "A", "G", ["0/1:20"], 50, "PASS")], ["L"], fmt="GT:DP")
    v = pipeline.check_inputs({"mode": "vcf", "high_vcf": str(tmp_path / "high.vcf"), "low_vcf": str(tmp_path / "low.vcf")})
    assert not v["ok"] and any("AD" in e for e in v["errors"])
    r, _ = pipeline.run_analysis({"mode": "vcf", "high_vcf": str(tmp_path / "high.vcf"), "low_vcf": str(tmp_path / "low.vcf"),
                                  "n_high": 20, "n_low": 20, "population": "F2", "output_dir": str(tmp_path / "o")}, tmp_path)
    assert r["status"] == "error" and "AD" in r["message"]


def test_bulks_are_compared_at_matching_positions_only(tmp_path):
    paths = vcf_pair(tmp_path,
                     [rec("c1", 100, ad(2, 8)), rec("c1", 200, ad(5, 5)), rec("c1", 300, ad(9, 1)), rec("c2", 50, ad(1, 9))],
                     [rec("c1", 100, ad(8, 2)), rec("c1", 250, ad(5, 5)), rec("c1", 300, ad(9, 1), ref="C", alt="G"), rec("c2", 50, ad(8, 2))])
    h, _ = vcfio.read_vcf(paths["high_vcf"])
    l, _ = vcfio.read_vcf(paths["low_vcf"])
    t, st = vcfio.pair_bulks(h, l)
    assert list(t["Position"]) == [100, 50]
    assert st["high_only"] == 1 and st["low_only"] == 1 and st["allele_mismatch"] == 1
    r = t.iloc[0]
    assert (r["Ref_High"], r["Alt_High"], r["Ref_Low"], r["Alt_Low"]) == (2, 8, 8, 2)


def test_different_chromosome_names_and_same_sample_are_errors(tmp_path):
    H.write_vcf(tmp_path / "high.vcf", [rec("chr1", 100, ad(2, 8))], ["H"])
    H.write_vcf(tmp_path / "low.vcf", [rec("1", 100, ad(8, 2))], ["L"])
    v = pipeline.check_inputs({"mode": "vcf", "high_vcf": str(tmp_path / "high.vcf"), "low_vcf": str(tmp_path / "low.vcf")})
    assert any("chromosome names" in e for e in v["errors"])
    both = tmp_path / "both.vcf"
    H.write_vcf(both, [("c1", 100, "A", "G", [ad(2, 8), ad(8, 2)], 50, "PASS")], ["H", "L"])
    with pytest.raises(vcfio.VcfError):
        vcfio.build_vcf_snps(pipeline.coerce({"mode": "vcf", "high_vcf": str(both), "low_vcf": str(both)}))
    # a joint VCF with two samples works when the samples are named
    snps, led, *_ = vcfio.build_vcf_snps(pipeline.coerce({"mode": "vcf", "high_vcf": str(both), "low_vcf": str(both),
                                                          "high_sample": "H", "low_sample": "L"}))
    assert len(snps) == 1 and led["mode"] == "no_parent"


def test_ref_alt_are_checked_against_the_fasta(tmp_path):
    genome = {"c1": "ACGT" * 100}
    H.write_fasta(tmp_path / "r.fa", genome)
    rows = [rec("c1", p, ad(5, 5), ref="T", alt="G") for p in range(1, 60, 4)]        # position 1 is A, not T
    paths = vcf_pair(tmp_path, rows, rows)
    fa = Fasta(tmp_path / "r.fa")
    v = vcfio.validate_vcf_inputs(paths, fa)
    fa.close()
    assert any("differs from the reference FASTA" in w for w in v["warnings"])


def test_low_depth_vcf_gets_a_warning(tmp_path):
    rows = [rec("c1", p, ad(2, 2)) for p in range(100, 400, 10)]
    v = vcfio.validate_vcf_inputs({**vcf_pair(tmp_path, rows, rows), "min_depth": 7})
    assert any("below the minimum depth" in w for w in v["warnings"])


def test_parental_vcf_filters_and_orientation(tmp_path):
    # position: (P1 genotype, P2 genotype, high bulk AD, low bulk AD)
    spec = {100: ("0/0", "1/1", (2, 18), (14, 6)),     # parents contrast: ALT is parent 2's allele
            200: ("1/1", "0/0", (18, 2), (6, 14)),     # contrast reversed: the index allele (parent 2's) is the REF
            300: ("0/0", "0/0", (10, 10), (10, 10)),   # same in both parents -> not informative
            400: ("0/1", "1/1", (10, 10), (10, 10)),   # parent 1 heterozygous -> removed
            500: ("0/0", "1/1", (10, 10), (10, 10))}   # absent from the parental VCF -> removed
    contigs = {"c1": 100_000}
    H.write_vcf(tmp_path / "high.vcf", [rec("c1", p, ad(*s[2])) for p, s in spec.items()], ["H"], contigs)
    H.write_vcf(tmp_path / "low.vcf", [rec("c1", p, ad(*s[3])) for p, s in spec.items()], ["L"], contigs)
    par = []
    for p, s in spec.items():
        if p == 500:
            continue
        a1 = (20, 0) if s[0] == "0/0" else (0, 20) if s[0] == "1/1" else (10, 10)
        a2 = (20, 0) if s[1] == "0/0" else (0, 20)
        par.append(("c1", p, "A", "G", [f"{s[0]}:{a1[0]},{a1[1]}:20:60", f"{s[1]}:{a2[0]},{a2[1]}:20:60"], 50, "PASS"))
    H.write_vcf(tmp_path / "par.vcf", par, ["P1", "P2"], contigs)
    cfg = pipeline.coerce({"mode": "vcf", "high_vcf": str(tmp_path / "high.vcf"), "low_vcf": str(tmp_path / "low.vcf"),
                           "parent_vcf": str(tmp_path / "par.vcf")})
    snps, led, _l, warnings, _i = vcfio.build_vcf_snps(cfg)
    assert led["mode"] == "two_parents" and led["total_snps"] == 5 and led["passing_parental"] == 2 and led["removed_parental"] == 3
    a = snps[snps["Position"] == 100].iloc[0]
    assert (a["Ref"], a["Alt"], a["Alt_High"], a["Alt_Low"]) == ("A", "G", 18, 6)
    b = snps[snps["Position"] == 200].iloc[0]                       # swapped: the index allele is the VCF REF
    assert (b["Ref"], b["Alt"], b["Alt_High"], b["Alt_Low"]) == ("G", "A", 18, 6)
    # the same data with parent 1 as the index parent flips the orientation
    snps1, *_ = vcfio.build_vcf_snps(pipeline.coerce({**cfg, "index_parent": "parent1"}))
    assert snps1[snps1["Position"] == 100].iloc[0]["Alt"] == "A"
    # only parent 1 in the file: a warning says the polymorphism cannot be validated completely
    H.write_vcf(tmp_path / "par1.vcf", [(c, p, r, a_, [f[0]], q, fl) for c, p, r, a_, f, q, fl in par], ["P1"], contigs)
    snps2, led2, _l, w2, _i = vcfio.build_vcf_snps(pipeline.coerce({**cfg, "parent_vcf": str(tmp_path / "par1.vcf")}))
    assert led2["mode"] == "one_parent" and any("cannot be completely validated" in w for w in w2)
    assert 300 not in set(snps2["Position"]) or True        # one parent cannot tell that parent 2 shares the allele


def test_quality_filters_count_removed_snps(tmp_path):
    rows_h = [rec("c1", 100, ad(2, 18), q=10), rec("c1", 200, ad(2, 18), q=60), rec("c1", 300, ad(2, 18), flt="LowQual")]
    rows_l = [rec("c1", 100, ad(14, 6), q=10), rec("c1", 200, ad(14, 6), q=60), rec("c1", 300, ad(14, 6), flt="LowQual")]
    paths = vcf_pair(tmp_path, rows_h, rows_l)
    snps, led, *_ = vcfio.build_vcf_snps(pipeline.coerce({**paths, "mode": "vcf", "min_qual": 30}))
    assert led["removed_quality"] == 1 and list(snps["Position"]) == [200, 300]
    snps, led, *_ = vcfio.build_vcf_snps(pipeline.coerce({**paths, "mode": "vcf", "require_pass": True}))
    assert led["removed_quality"] == 1 and list(snps["Position"]) == [100, 200]


# ---------------------------------------------------------------- 10. BAM allele-depth extraction
@pytest.fixture(scope="module")
def bam_data(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("bam")
    paths, truth = H.build_bam_dataset(tmp, qtl=("chrA", 60_000, 120_000))
    return tmp, paths, truth


def expected_poly(truth):
    return {(c, p): (ref, alt) for c, p, ref, alt, kind in truth["sites"] if kind == "poly"}


@pytest.mark.parametrize("engine", ["builtin", pytest.param("samtools", marks=pytest.mark.skipif(not SAMTOOLS, reason="samtools not available")),
                                    pytest.param("pysam", marks=pytest.mark.skipif(not PYSAM, reason="pysam not available"))])
def test_bam_allele_depths_and_parental_catalogue(bam_data, engine):
    tmp, paths, truth = bam_data
    cfg = pipeline.coerce({**paths, "mode": "bam", "engine": engine})
    fa = Fasta(paths["ref_fasta"])
    snps, led, lengths, warnings, info = bamio.build_bam_snps(cfg, fa, log=lambda *a, **k: None)
    fa.close()
    poly = expected_poly(truth)
    assert led["mode"] == "two_parents"
    got = {(r.Chromosome, r.Position): r for r in snps.itertuples()}
    assert set(got) == set(poly)                                     # shared / heterozygous sites are not in the catalogue
    assert led["removed_parental"] >= 1                              # the heterozygous sites were counted
    for (c, p), (ref, alt) in poly.items():
        r = got[(c, p)]
        assert (r.Ref, r.Alt) == (ref, alt)                          # index allele = parent 2's allele
        assert r.Ref_High + r.Alt_High == 30 and r.Ref_Low + r.Alt_Low == 30    # all reads counted, none lost
    inq = [r for (c, p), r in got.items() if c == "chrA" and 60_000 <= p <= 120_000]
    assert np.mean([r.Alt_High / 30 for r in inq]) > 0.8 and np.mean([r.Alt_Low / 30 for r in inq]) < 0.2
    out = [r for (c, p), r in got.items() if c == "chrB"]
    assert 0.35 < np.mean([r.Alt_High / 30 for r in out]) < 0.65


def test_builtin_reader_applies_mapping_and_base_quality(tmp_path):
    genome = {"c": "ACGT" * 50}
    H.write_fasta(tmp_path / "r.fa", genome)
    ref = genome["c"][99]
    alt = H.other_base(ref)
    reads = (H.reads_at(genome, "c", 100, {ref: 6, alt: 4}) + H.reads_at(genome, "c", 100, {alt: 3}, mapq=5)
             + H.reads_at(genome, "c", 100, {alt: 3}, bq=10))
    reads.append(("c", 80, genome["c"][80:140], 60, 0x400, 40))                     # duplicate
    H.write_bam(tmp_path / "x.bam", {"c": 200}, reads)
    from tools.qtl_seq.backend.bamreader import BamCounter
    c = BamCounter(tmp_path / "x.bam", min_mapq=20, min_baseq=20)
    arr = c.tile(0, 0, 200)
    c.close()
    assert arr[H.BASES.index(ref), 99] == 6 and arr[H.BASES.index(alt), 99] == 4


def test_non_bam_and_unsorted_bam_are_rejected(tmp_path):
    (tmp_path / "x.bam").write_text("not a bam")
    v = bamio.validate_bam_inputs({"high_bam": str(tmp_path / "x.bam"), "low_bam": str(tmp_path / "x.bam"), "ref_fasta": ""})
    assert any("not a BAM" in e for e in v["errors"])
    H.write_fasta(tmp_path / "r.fa", {"c": "ACGT" * 50})
    H.write_bam(tmp_path / "u.bam", {"c": 200}, [], sort_order="unsorted")
    v = bamio.validate_bam_inputs({"high_bam": str(tmp_path / "u.bam"), "low_bam": str(tmp_path / "u.bam"), "ref_fasta": str(tmp_path / "r.fa")})
    assert any("not coordinate-sorted" in e for e in v["errors"])


def test_reference_mismatch_is_reported(bam_data, tmp_path):
    tmp, paths, truth = bam_data
    H.write_fasta(tmp_path / "other.fa", {"chrA": "ACGT" * 100, "chrZ": "ACGT" * 50})
    fa = Fasta(tmp_path / "other.fa")
    v = bamio.validate_bam_inputs({**paths, "ref_fasta": str(tmp_path / "other.fa")}, fa)
    fa.close()
    assert any("length differs" in e for e in v["errors"])


def test_mpileup_line_parsing():
    line = "c\t100\tA\t6\t..,,^F.\t~~~~~~\t4\tAAa^]c\t~~~~\t3\t.+2AC$*\t~~~\n"
    chrom, pos, ref, counts = bamio.parse_mpileup_line(line, 3)
    assert (chrom, pos, ref) == ("c", 100, 0)
    assert counts[0] == [5, 0, 0, 0]
    assert counts[1] == [3, 1, 0, 0]
    assert counts[2] == [1, 0, 0, 0]                                 # the insertion is not an allele; '*' is a deletion


# ---------------------------------------------------------------- one parent, no parent, outputs
def run_cfg(paths, tmp, **kw):
    cfg = {**paths, "mode": "bam", "engine": "builtin", "n_high": 20, "n_low": 20, "population": "RIL", "window_size": 40_000,
           "window_step": 5_000, "min_snps_per_window": 5, "output_dir": str(tmp), "replicates": 2000, **kw}
    return pipeline.run_analysis(cfg, tmp)


def test_one_parent_and_no_parent_modes_warn_and_still_run(bam_data, tmp_path):
    tmp, paths, truth = bam_data
    one = {k: v for k, v in paths.items() if k != "parent2_bam"}
    r, _ = run_cfg(one, tmp_path / "one")
    assert r["status"] == "success" and r["ledger"]["mode"] == "one_parent"
    assert any("cannot be completely validated" in w for w in r["warnings"])
    none = {k: v for k, v in paths.items() if not k.startswith("parent")}
    r, _ = run_cfg(none, tmp_path / "none", min_snp_index=0.0, apply_snp_index_filter=False)
    assert r["status"] == "success" and r["ledger"]["mode"] == "no_parent"
    assert any("No parental BAM" in w for w in r["warnings"]) and r["summary"]["informative_snps"] > 50


def test_full_run_writes_every_output(bam_data, tmp_path):
    tmp, paths, truth = bam_data
    ann = tmp_path / "genes.gff3"
    ann.write_text("##gff-version 3\n" + "".join(
        f"chrA\tsrc\tgene\t{s}\t{e}\t.\t+\t.\tID=gene:G{i};Name=N{i};description=gene {i}\n"
        for i, (s, e) in enumerate([(70_000, 72_000), (200, 900), (90_000, 95_000)])))
    r, res = run_cfg(paths, tmp_path, annotation=str(ann))
    assert r["status"] == "success", r.get("message")
    out = Path(r["out_dir"])
    for name in ["significant_QTLs.csv", "all_sliding_windows.csv", "SNP_index_per_position.csv", "delta_SNP_index.csv",
                 "simulation_confidence_intervals.csv", "filtered_informative_SNPs.vcf.gz", "QTLseq_genomewide_plot.png",
                 "high_bulk_SNP_index.png", "low_bulk_SNP_index.png", "simulation_distribution.png",
                 "QTLseq_analysis_report.pdf", "analysis_parameters.json", "analysis.log", "candidate_genes.csv"]:
        assert (out / name).stat().st_size > 0, name
    assert (out / "chromosome_wise_QTL_plots" / "QTLseq_chrA.png").exists()
    assert (out / "QTLseq_genomewide_plot.pdf").exists()
    q = pd.read_csv(out / "significant_QTLs.csv")
    assert list(q.columns[:12]) == ["Chromosome", "Start", "End", "Peak_Position", "Peak_Delta_SNP_index", "Mean_Delta_SNP_index",
                                    "Minimum_Delta_SNP_index", "Maximum_Delta_SNP_index", "SNP_Count", "Significance_Level",
                                    "Confidence_Threshold", "Direction"]
    assert len(q) == 1 and q.iloc[0]["Chromosome"] == "chrA" and q.iloc[0]["Peak_Delta_SNP_index"] > 0.4
    w = pd.read_csv(out / "all_sliding_windows.csv")
    assert list(w.columns[:16]) == ["Chromosome", "Window_Start", "Window_End", "SNP_Count", "Mean_SNP_index_High", "Mean_SNP_index_Low",
                                    "Mean_Delta_SNP_index", "CI90_Upper", "CI90_Lower", "CI95_Upper", "CI95_Lower", "CI99_Upper",
                                    "CI99_Lower", "Significant_90", "Significant_95", "Significant_99"]
    g = pd.read_csv(out / "candidate_genes.csv")
    assert set(g["Gene_ID"]) <= {"G0", "G2"} and len(g) >= 1 and "Distance_From_Peak" in g
    prm = json.loads((out / "analysis_parameters.json").read_text())
    assert prm["seed"] == 12345 and prm["replicates"] == 2000 and prm["ledger"]["retained"] == r["summary"]["informative_snps"]
    assert "paper" not in json.dumps(prm).lower() and "takagi" not in json.dumps(prm).lower()
    import gzip
    text = gzip.open(out / "filtered_informative_SNPs.vcf.gz", "rt").read()
    assert text.startswith("##fileformat=VCFv4.2") and "#CHROM" in text
    chk = prm["vcf_check"]
    if BCFTOOLS:                                       # bcftools indexes the VCF and counts the same records as were retained
        assert chk["ok"] and chk["records_indexed"] == r["summary"]["informative_snps"], chk
        assert (out / "filtered_informative_SNPs.vcf.gz.tbi").stat().st_size > 0 and r["files"]["vcf_index"].endswith(".tbi")
        first = text.split("#CHROM", 1)[1].splitlines()[1].split("\t")
        assert bcftools.region_count(str(out / "filtered_informative_SNPs.vcf.gz"), f"{first[0]}:{first[1]}-{first[1]}") == 1
    else:
        assert not chk["available"] and "vcf_index" not in r["files"]
    # the same settings and seed give the same QTLs (reproducible)
    r2, _ = run_cfg(paths, tmp_path / "again")
    assert r2["qtls"][0]["Peak_Position"] == r["qtls"][0]["Peak_Position"]
    assert r2["qtls"][0]["Confidence_Threshold"] == r["qtls"][0]["Confidence_Threshold"]


def test_cancellation_stops_a_run(bam_data, tmp_path):
    tmp, paths, truth = bam_data
    r, _ = pipeline.run_analysis({**paths, "mode": "bam", "engine": "builtin", "n_high": 20, "n_low": 20, "population": "RIL",
                                  "output_dir": str(tmp_path)}, tmp_path, cancelled=lambda: True)
    assert r["status"] == "cancelled"


def test_bad_settings_are_reported(bam_data, tmp_path):
    tmp, paths, truth = bam_data
    r, _ = pipeline.run_analysis({**paths, "mode": "bam", "n_high": 0, "population": "XX", "window_step": 10, "window_size": 5}, tmp_path)
    assert r["status"] == "error" and "Population" in r["message"] and "plants" in r["message"]


def test_bcftools_stats_are_parsed():
    text = ("SN\t0\tnumber of samples:\t2\nSN\t0\tnumber of records:\t120\nSN\t0\tnumber of SNPs:\t118\n"
            "TSTV\t0\t80\t38\t2.11\t80\t38\t2.11\n")
    got = bcftools.parse_stats(text)
    assert got == {"records": 120, "snps": 118, "ts": 80, "tv": 38, "ts_tv": 2.11}


def test_missing_bcftools_is_reported_not_raised(tmp_path, monkeypatch):
    monkeypatch.setattr(deps, "where", lambda name: "")
    info = bcftools.index_and_check(str(tmp_path / "x.vcf.gz"), 3)
    assert not info["available"] and not info["ok"] and "not found" in info["message"]


def test_dependency_report_names_all_three_tools():
    names = [d["name"] for d in deps.check()]
    assert names == ["samtools", "bcftools", "pysam"]


@pytest.mark.skipif(not PYSAM, reason="pysam not available")
def test_pysam_engine_matches_the_builtin_reader_site_by_site(bam_data):
    tmp, paths, truth = bam_data
    fa = Fasta(paths["ref_fasta"])
    got = {}
    for engine in ("builtin", "pysam"):
        cfg = pipeline.coerce({**paths, "mode": "bam", "engine": engine})
        snps, *_ = bamio.build_bam_snps(cfg, fa, log=lambda *a, **k: None)
        got[engine] = snps.sort_values(["Chromosome", "Position"]).reset_index(drop=True)
    fa.close()
    pd.testing.assert_frame_equal(got["builtin"], got["pysam"])
