"""Gene Extractor: the QTL table of QTL Mapping (cM) plus a REQUIRED marker position file or VCF."""
from pathlib import Path

from tools.gene_extractor.backend.gene_extractor_runner import run_extraction

D = Path(__file__).resolve().parents[3] / "Demo_Input_Files" / "10_Genes_Extraction" / "input"
QTL = D / "qtl_table_cM_markers_without_position.csv"
MAP, VCF, GFF = D / "marker_positions.csv", D / "marker_positions.vcf", D / "annotation_IRGSP1_regions.gff3"


def _run(tmp_path, **kw):
    cfg = {"marker_path": str(QTL), "gff_path": str(GFF), "positions_path": str(MAP), "match": "overlap"}
    cfg.update(kw)
    return run_extraction(cfg, tmp_path / "out")


def test_qtl_table_with_a_marker_table_or_a_vcf_gives_the_same_genes(tmp_path):
    a = _run(tmp_path)
    b = _run(tmp_path, positions_path=str(VCF))
    assert a["status"] == "success" and b["status"] == "success"
    assert a["gene_ids"] == b["gene_ids"] == 22 and len(a["region_summary"]) == 2


def test_the_marker_position_file_is_required(tmp_path):
    r = _run(tmp_path, positions_path="")
    assert r["status"] == "error" and "marker position file" in r["message"]
    cfg = {"marker_path": str(QTL), "gff_path": str(GFF)}
    assert run_extraction(cfg, tmp_path / "o2")["status"] == "error"


def test_a_table_of_intervals_in_bp_is_no_longer_an_input(tmp_path):
    files = {
        "headed.csv": "QTL_ID,QTL_Start,QTL_End\nQTL_1.1,20510000,20600000\nQTL_5.1,7230000,7410000\n",
        "plain.csv": "Chromosome,Start,End\n1,20510000,20600000\n5,7230000,7410000\n",
        "noheader.csv": "QTL_1.1,20510000,20600000\nQTL_5.1,7230000,7410000\n",
    }
    for name, text in files.items():
        f = tmp_path / name
        f.write_text(text)
        r = _run(tmp_path, marker_path=str(f))
        assert r["status"] == "error" and "table of intervals" in r["message"] and "QTL table" in r["message"], name


def test_marker_names_match_across_chr_prefix_and_the_name_chromosome_is_used(tmp_path):
    import pandas as pd
    from tools.gene_extractor.backend.gff_genes import read_regions
    q = pd.DataFrame({"QTL_ID": ["QTL_4.1"], "Chromosome": [4], "Left_Marker": ["chr01_18102213"], "Peak_Marker": ["chr01_20562598"],
                      "Right_Marker": ["chr01_37920337"], "Left_Position": [20], "Peak_Position": [25], "Right_Position": [30],
                      "CI_Start_cM": [20], "CI_End_cM": [30]})
    m = pd.DataFrame({"Marker": ["01_18102213", "01_37920337"], "Chromosome": [1, 1], "bp": [18102213, 37920337]})
    qp, mp = tmp_path / "q.csv", tmp_path / "m.csv"
    q.to_csv(qp, index=False)
    m.to_csv(mp, index=False)
    regions, notes = read_regions(str(qp), str(mp))          # 01_... in the file matches chr01_... in the QTL table
    assert [(r["Chromosome"], int(r["Start"]), int(r["Stop"])) for r in regions] == [("1", 18102213, 37920337)]
    regions, notes = read_regions(str(qp), str(tmp_path / "m.csv").replace("m.csv", "m.csv"))
    assert regions
    # no marker file entry at all: the name carries chromosome 1 although the QTL table says 4 -> the name wins, with a note
    m.iloc[:0].to_csv(mp, index=False)
    m2 = pd.DataFrame({"Marker": ["other"], "Chromosome": [1], "bp": [5]})
    m2.to_csv(mp, index=False)
    regions, notes = read_regions(str(qp), str(mp))
    assert regions[0]["Chromosome"] == "1" and any("differs" in n for n in notes)
