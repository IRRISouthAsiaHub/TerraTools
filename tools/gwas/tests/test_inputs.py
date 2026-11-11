"""Genotype encodings, VCF genotype parsing and quality control."""
import numpy as np
import pytest

from tools.gwas.backend import genotype_io as gio, qc as qcmod


# ------------------------------------------------------------------ encodings
def test_dosage_is_kept():
    G = np.array([[0, 1, 2], [2, 1, np.nan]])
    out, rep = gio.to_dosage(G, "auto")
    assert rep["encoding"] == "dosage" and not rep["converted"]
    assert np.array_equal(out, G, equal_nan=True)


def test_centred_is_converted_and_logged():
    G = np.array([[-1, 0, 1], [1, 0, np.nan]])
    out, rep = gio.to_dosage(G, "auto")
    assert rep["encoding"] == "centred" and rep["converted"]
    assert np.array_equal(out, np.array([[0, 1, 2], [2, 1, np.nan]]), equal_nan=True)
    assert "Converted" in rep["message"]


@pytest.mark.parametrize("G", [np.array([[0, 1], [1, 0], [0, 0]], float), np.array([[0.0, 1.0, np.nan], [1.0, 1.0, 0.0]])])
def test_zero_one_matrix_is_never_guessed(G):
    with pytest.raises(gio.InputError, match="never guesses"):
        gio.to_dosage(G, "auto")


def test_zero_one_matrix_with_explicit_choice():
    G = np.array([[0, 1], [1, 0]], float)
    d, _ = gio.to_dosage(G, "dosage")
    c, rep = gio.to_dosage(G, "centred")
    assert np.array_equal(d, G) and np.array_equal(c, G + 1) and rep["converted"]


def test_values_must_fit_the_chosen_encoding():
    with pytest.raises(gio.InputError, match="not valid"):
        gio.to_dosage(np.array([[0, 1, 2]], float), "centred")
    with pytest.raises(gio.InputError, match="not valid"):
        gio.to_dosage(np.array([[-1, 0, 1]], float), "dosage")


def test_mixed_negative_and_large_values_are_rejected():
    with pytest.raises(gio.InputError):
        gio.to_dosage(np.array([[-1, 0, 2]], float), "auto")


def test_missing_codes_are_not_genotypes():
    G = np.array([[0, 1, 2, -9], [2, -99, 1, 0]], float)
    out, rep = gio.to_dosage(G, "auto")
    assert np.isnan(out[0, 3]) and np.isnan(out[1, 1]) and rep["missing_calls"] == 2
    # -1 is a missing code only if the user says so under dosage coding
    G2 = np.array([[0, 1, 2, -1]], float)
    out2, _ = gio.to_dosage(G2, "dosage", extra_missing=[-1])
    assert np.isnan(out2[0, 3])


def test_read_file_preserves_ids_and_map(tmp_path):
    f = tmp_path / "g.csv"
    f.write_text("ID,Phenotype,Chr1_100,Chr1_200,Chr2_50\nA,1.5,0,1,2\nB,2.5,2,NA,0\nC,3.5,1,1,1\n")
    d = gio.read_gwas_file(path=str(f), encoding="dosage")
    assert list(d.ids) == ["A", "B", "C"] and list(d.snp_ids) == ["Chr1_100", "Chr1_200", "Chr2_50"]
    assert list(d.chrom) == [1, 1, 2] and list(d.pos) == [100, 200, 50] and d.marker_map.valid
    assert np.isnan(d.G[1, 1]) and d.G[0, 2] == 2


def test_markers_without_positions_are_marked_invalid(tmp_path):
    f = tmp_path / "g.csv"
    f.write_text("ID,Phenotype,m1,m2,m3\nA,1,0,1,2\nB,2,2,1,0\n")
    d = gio.read_gwas_file(path=str(f), encoding="dosage")
    assert not d.marker_map.valid


def test_duplicate_ids_and_markers_are_errors(tmp_path):
    f = tmp_path / "g.csv"
    f.write_text("ID,Phenotype,Chr1_1,Chr1_2,Chr1_3\nA,1,0,1,2\nA,2,2,1,0\n")
    with pytest.raises(gio.InputError, match="more than once"):
        gio.read_gwas_file(path=str(f), encoding="dosage")
    f.write_text("ID,Phenotype,Chr1_1,Chr1_1,Chr1_3\nA,1,0,1,2\nB,2,2,1,0\n")
    with pytest.raises(gio.InputError, match="more than once"):
        gio.read_gwas_file(path=str(f), encoding="dosage")


# ------------------------------------------------------------------ VCF genotypes
@pytest.mark.parametrize("gt,expected", [("0/0", 0), ("0|0", 0), ("0/1", 1), ("1|0", 1), ("1/1", 2), ("1|1", 2), ("0/1:35,12:47", 1)])
def test_parse_gt_diploid(gt, expected):
    assert gio.parse_gt(gt) == expected


@pytest.mark.parametrize("gt", ["./.", ".|.", "./1", "1|.", "0/."])
def test_parse_gt_missing_and_partially_missing(gt):
    assert np.isnan(gio.parse_gt(gt))


@pytest.mark.parametrize("gt", ["1", "0", "0/0/1", "0/1/1/1", "0/2", "2|1"])
def test_parse_gt_rejects_other_ploidy_and_alleles(gt):
    with pytest.raises(gio.VcfGenotypeError):
        gio.parse_gt(gt, " at Chr1:100")


def test_pooled_bulk_samples_are_refused():
    gio.reject_pooled(["IR64", "Nipponbare"])
    with pytest.raises(gio.InputError, match="pooled"):
        gio.reject_pooled(["IR64", "HIGH_BULK"])
    with pytest.raises(gio.InputError, match="pooled"):
        gio.reject_pooled(["pool1"])


def test_vcf_reader_end_to_end(tmp_path):
    from tools.gwas.backend import gwas_inputs as gi
    samples = [f"S{i}" for i in range(24)]
    hdr = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n"
    rng = np.random.default_rng(0)
    rows = []
    for pos in range(100, 160, 10):
        gts = [rng.choice(["0/0", "0|1", "1/1", "./."]) for _ in samples]
        gts[0] = "./1"                                    # one missing allele
        rows.append(f"Chr1\t{pos}\t.\tA\tG\t50\tPASS\t.\tGT:DP\t" + "\t".join(g + ":10" for g in gts))
    rows.append("Chr1\t200\t.\tA\tG,T\t50\tPASS\t.\tGT\t" + "\t".join("0/1" for _ in samples))      # multi-allelic: skipped
    rows.append("Chr1\t210\t.\tAT\tA\t50\tPASS\t.\tGT\t" + "\t".join("0/1" for _ in samples))        # indel: skipped
    v = tmp_path / "a.vcf"
    v.write_text(hdr + "\n".join(rows) + "\n")
    ph = tmp_path / "p.csv"
    ph.write_text("ID,Trait\n" + "\n".join(f"{s},{i}" for i, s in enumerate(samples)) + "\n")
    path, info = gi.build_from_vcf(str(v), str(ph), "Trait", str(tmp_path), maf=0.0, max_missing=1.0, log=lambda m: None)
    assert info["skipped_multiallelic"] == 1 and info["skipped_indel"] == 1 and info["partially_missing_calls"] == 6
    d = gio.read_gwas_file(path=path, encoding="dosage")
    assert d.G.shape == (24, 6) and np.isnan(d.G[0]).all()
    # a haploid call stops the run
    v.write_text(hdr + "Chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t" + "\t".join("1" for _ in samples) + "\n")
    with pytest.raises(gio.InputError, match="ploidy"):
        gi.build_from_vcf(str(v), str(ph), "Trait", str(tmp_path), log=lambda m: None)
    # a pooled bulk sample stops the run
    v.write_text(hdr.replace("S3", "high_bulk") + "Chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t" + "\t".join("0/1" for _ in samples) + "\n")
    ph.write_text("ID,Trait\n" + "\n".join(f"{s},{i}" for i, s in enumerate(["high_bulk" if s == "S3" else s for s in samples])) + "\n")
    with pytest.raises(gio.InputError, match="pooled"):
        gi.build_from_vcf(str(v), str(ph), "Trait", str(tmp_path), log=lambda m: None)


# ------------------------------------------------------------------ quality control
def make(n=40, m=30, seed=0):
    rng = np.random.default_rng(seed)
    G = rng.binomial(2, rng.uniform(0.2, 0.5, m), size=(n, m)).astype(float)
    chrom = np.repeat([1, 2], m // 2)
    pos = np.tile(np.arange(m // 2) * 1000 + 1, 2)
    return G, [f"s{i}" for i in range(m)], chrom, pos, [f"a{i}" for i in range(n)]


def test_ledger_reports_every_stage():
    G, snp, ch, pos, ids = make()
    G[:, 0] = 1                                       # monomorphic
    G[:30, 1] = np.nan                                # SNP call rate
    G[5, :] = np.nan                                  # sample call rate
    G[:, 2] = (np.arange(40) == 3).astype(float)      # one carrier: MAF 1/80 (sample 5 is gone)
    pos2 = pos.copy(); pos2[4] = pos2[3]; ch2 = ch.copy(); ch2[4] = ch2[3]       # duplicate position
    r = qcmod.run_qc(G, snp, ch2, pos2, ids, qcmod.QCSettings(min_maf=0.05, imputation="mean"))
    st = {s["stage"]: s for s in r.ledger}
    assert st["Duplicate positions"]["snps_removed"] == 1
    assert st["SNP call rate"]["snps_removed"] >= 1
    assert st["Sample call rate"]["samples_removed"] == 1
    assert st["Monomorphic SNPs"]["snps_removed"] >= 1
    assert st["Minor allele frequency"]["snps_removed"] >= 1
    assert r.G.shape == (len(r.sample_idx), len(r.snp_idx)) and not np.isnan(r.G).any()
    assert r.snp_table["Status"][0] != "kept" and r.sample_table["Status"][5] == "sample call rate"


def test_fewer_than_two_snps_is_a_clear_error_never_a_bypass():
    G, snp, ch, pos, ids = make()
    with pytest.raises(qcmod.QCError, match="Only 0 SNP"):
        qcmod.run_qc(G, snp, ch, pos, ids, qcmod.QCSettings(min_maf=0.51))     # no marker can have MAF above 0.5
    G2 = np.zeros_like(G); G2[:, 0] = np.arange(40) % 3
    with pytest.raises(qcmod.QCError, match="at least 2"):
        qcmod.run_qc(G2, snp, ch, pos, ids, qcmod.QCSettings())


def test_qc_never_alters_the_input_and_keeps_observed_frequencies():
    G, snp, ch, pos, ids = make()
    G[2, 7] = np.nan
    before = G.copy()
    r = qcmod.run_qc(G, snp, ch, pos, ids, qcmod.QCSettings())
    assert np.array_equal(G, before, equal_nan=True)
    j = list(r.snp_idx).index(7)
    obs = G[:, 7][~np.isnan(G[:, 7])]
    assert r.maf[j] == pytest.approx(min(obs.mean() / 2, 1 - obs.mean() / 2))
    assert np.isnan(r.G_observed[2, j]) and r.G[2, j] == pytest.approx(obs.mean())     # only the missing cell changed
    assert r.imputed_calls == 1 and r.notes and "not haplotype" in r.notes[0]


def test_imputation_options():
    G, snp, ch, pos, ids = make()
    G[1, 3] = np.nan
    removed = qcmod.run_qc(G, snp, ch, pos, ids, qcmod.QCSettings(imputation="remove"))
    assert 3 not in removed.snp_idx and not np.isnan(removed.G).any()
    with pytest.raises(qcmod.QCError, match="imputation is switched off"):
        qcmod.run_qc(G, snp, ch, pos, ids, qcmod.QCSettings(imputation="none"))
    ok = qcmod.run_qc(make()[0], snp, ch, pos, ids, qcmod.QCSettings(imputation="none"))
    assert ok.imputed_calls == 0
