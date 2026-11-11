"""
Captures every screenshot used in the TerraTools user manual from the real app, using the
example files in docs/sample_data.   python docs/make_screenshots.py [scenario ...]
Output: docs/manual/shots/*.jpg + manifest.json (callout captions, read by build_manual.py).
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
os.environ.setdefault("TERRATOOLS_RSCRIPT", str(HERE.parent / "r_runtime_build" / "bin" / "Rscript.exe"))
from shotlib import App   # noqa: E402

S = HERE / "sample_data"
OUT = HERE / "manual" / "shots"


def card(title, note, **kw):
    return dict(text=title, sel="h3", closest=".card", note=note, **kw)


def label(text, note, **kw):
    return dict(text=text, sel="label", closest=".field", note=note, **kw)


# ───────────────────────── general ─────────────────────────
def sc_home(a):
    yield from a.route("", height=1000)
    yield 1200
    yield from a.shot("home", title="The TerraTools home screen", callouts=[
        dict(sel=".sidebar", global_=True, note="Sidebar: every tool, grouped as Trait Mapping, Genomic Selection and Data Analysis. Click a name to open it."),
        dict(sel=".quick-input", global_=True, note="Quick open (Ctrl+K): type part of a tool's name and press Enter."),
    ])


def sc_tour(a):
    yield from a.route("qtl", height=1000)
    yield from a.shot("interface_tour", title="Anatomy of a tool page", callouts=[
        dict(sel=".back-btn", global_=True, note="Back: returns to the screen you came from (Alt+← also works)."),
        dict(sel=".crumbs", global_=True, note="Breadcrumbs: where you are. Click Home to return to the dashboard."),
        dict(sel=".tab-btn", nth=0, parent=1, note="Tabs: most tools have an input tab first and a Results tab that opens automatically when a run finishes."),
        dict(sel=".file-drop", nth=0, note="File box: click Browse to choose a file. The path appears here once it is loaded."),
        dict(sel=".run-bar", note="Run bar: starts the analysis. The line underneath shows progress and any error."),
        dict(sel=".collapse-btn", global_=True, note="Collapse the sidebar to gain screen width."),
    ])


# ───────────────────────── 0. Genofile Converter ─────────────────────────
def sc_genofile(a):
    G = HERE.parent / "Demo_Input_Files" / "19_Genofile_Converter" / "input"
    yield from a.route("genofile", height=1250)
    yield from a.browse(0, G / "genotypes.vcf.gz", wait=3000)
    yield from a.shot("genofile_input", title="Genofile Converter — a VCF, its parents and the output settings", height=1900, callouts=[
        dict(sel=".file-drop", nth=0, note="The company file: HapMap, VCF / VCF.GZ or an Excel / CSV sheet of allele calls. The format is detected; change it with the buttons if it is wrong."),
        card("How calls are coded", "How an individual is coded against the two parents, and when a marker is left out."),
        card("Parents and samples", "Choose Parent 1 and Parent 2 among the samples of the file (P1 and P2 were found automatically); samples that do not belong to the population can be left out."),
        card("Output settings", "Marker names (as in the file or chromosome_position), chromosome labels (Chr01 becomes 1) and whether markers on scaffolds are left out."),
    ])
    yield from a.run_and_wait(120, settle=2500)
    yield from a.shot("genofile_result", title="Genofile Converter — the converted files", height=1900, callouts=[
        dict(sel=".summary-grid", note="Markers in the file, coded (kept), left out, individuals and the share of calls that are missing or unexpected."),
        card("Your converted file", "What the file contains (a row per marker, a column per individual, coded 0 / 1 / 2 / −1), the CSV and Excel downloads, and shortcuts to continue in another tool."),
        card("Why markers were left out", "Parent missing, parent heterozygous, or both parents with the same allele; the genotype codes of all kept markers are counted on the right."),
    ])


# ───────────────────────── 1. Mendelian ─────────────────────────
def sc_mendelian(a):
    yield from a.route("mendelian", height=1000)
    yield from a.browse(0, S / "mendelian" / "mendelian_RIL_example.csv")
    yield from a.shot("mendelian_input", title="Mendelian Segregation — input and settings", callouts=[
        dict(sel=".file-drop", note="Genotype CSV: Marker, Chromosome, then one column per individual (0 / 1 / 2)."),
        card("Population Type", "Population type: RIL tests a 1:1 ratio, F2 tests 1:2:1."),
        card("Significance", "Threshold (alpha) and the multiple-testing correction applied to the p-values (FDR by default)."),
        card("Missing Data", "Values to treat as missing, separated by commas."),
    ])
    yield from a.run_and_wait(120)
    yield from a.shot("mendelian_results", title="Mendelian Segregation — results", height=1250, callouts=[
        dict(sel=".result-card", nth=0, note="−log10(p) for every marker: high points are markers whose genotype ratio departs from the expectation."),
        dict(sel=".result-card", nth=1, note="Share of distorted markers on each chromosome."),
        dict(text="Filtered Markers", sel="button", note="Downloads: markers that pass the test, the full table and the LM-ready file for Linkage Mapping."),
        dict(sel=".summary-box", note="Summary: markers tested, flagged distorted, passing the quality check, and low-power markers."),
    ])


# ───────────────────────── 2. Linkage ─────────────────────────
def sc_linkage(a):
    yield from a.route("linkage", height=1250)
    yield from a.browse(0, S / "linkage" / "linkage_F2_coded_example.xlsx", wait=3000)
    yield from a.shot("linkage_upload", title="Linkage Mapping — Upload tab", callouts=[
        dict(sel=".file-drop", nth=0, note="Genotype Excel file: Marker, Chromosome, then one column per individual coded 0 / 1 / 2 / −1."),
        label("Population Type", "Population type: F2 or RIL."),
        dict(sel=".convert-hint", nth=0, note="Genotypes from a genotyping company (HapMap, VCF, A:G calls)? Convert them with the Genofile Converter first."),
        card("Data Preview", "Preview of the first rows — check that columns were read as expected."),
        card("Chromosome Assignments", "Markers found per chromosome."),
    ])
    yield from a.tab("Parameters")
    yield from a.shot("linkage_params", title="Linkage Mapping — Parameters tab", height=1100, callouts=[
        card("Analysis Parameters", "Missing-data limit, maximum recombination fraction, optional segregation-distortion filter and map scaling."),
        card("Parameter Guide", "What each parameter does."),
    ])
    yield from a.run_and_wait(180)
    yield from a.shot("linkage_results", title="Linkage Mapping — Results tab", height=1300, callouts=[
        dict(text="File for QTL", sel="h3", closest=".card", note="Ready-to-use input file for QTL Mapping: the map positions with the same markers, transposed."),
        card("Downloads", "Excel exports: full results, all markers with positions, retained markers and recombination table."),
    ])
    yield from a.tab("Map Plot", wait=1800)
    yield from a.shot("linkage_map", title="Linkage Mapping — Map Plot tab", height=1300, callouts=[
        card("Plot Options", "Choose chromosomes, theme, colours and sizes, then click Update Plot."),
        card("Genetic Linkage Map", "The linkage map: each bar is a linkage group, ticks are markers, positions in cM (Kosambi)."),
    ])


# ───────────────────────── 3. QTL ─────────────────────────
def sc_qtl(a):
    yield from a.route("qtl", height=1150)
    yield from a.browse(0, S / "qtl" / "qtl_genotypes_example.csv")
    yield from a.browse(1, S / "qtl" / "qtl_phenotypes_example.csv", wait=2200)
    yield from a.shot("qtl_input", title="QTL Mapping — input files and settings", callouts=[
        dict(sel=".file-drop", nth=0, note="Genotype CSV: three header rows (marker names, chromosome, position in cM), then one row per individual."),
        dict(sel=".file-drop", nth=1, note="Phenotype CSV: first column = individual ID, one column per trait."),
        label("Trait / Phenotype", "Trait: choose the phenotype column to map."),
        label("Mapping Method", "Mapping method: SIM, CIM and MQM are regression methods; RF, SVR and XGBoost are machine-learning scores."),
        label("Significance threshold", "Threshold: a fixed LOD, or a genome-wide permutation test."),
        label("Support interval", "Support-interval LOD drop: 1.5 gives roughly a 95% interval around each peak."),
    ])
    yield from a.run_and_wait(240)
    yield from a.shot("qtl_results", title="QTL Mapping — results", height=1350, callouts=[
        dict(text="Chromosome-wise Manhattan Plot", sel=".plot-title", note="LOD curve of every chromosome; shaded bands are the QTL regions."),
    ])
    yield from a.shot("qtl_table", title="QTL Mapping — the QTL table", height=1100, focus=dict(text="QTL regions", sel=".plot-title"), callouts=[
        dict(text="QTL regions", sel=".plot-title", note="One row per QTL: left, peak and right flanking marker, interval size, LOD, p-value, effects."),
        dict(sel="table", nth=0, note="Scroll the table sideways for more columns (effect size, PVE, confidence intervals)."),
    ])


# ───────────────────────── 4. Genes extraction ─────────────────────────
def sc_gene(a):
    yield from a.route("gene_extractor", height=950)
    yield from a.browse(0, S / "gene_extractor" / "qtl_table_example.csv")
    yield from a.browse(1, S / "gene_extractor" / "marker_positions_example.csv")
    yield from a.browse(2, S / "gene_extractor" / "annotation_example.gff3")
    yield from a.shot("gene_input", title="Genes Extraction — input files", height=1850, callouts=[
        card("QTL Results (input)", "The QTL Intervals file or the QTL table of QTL Mapping (Chromosome, Start, End) — or click “Use QTL Mapping output” to take the last QTL Mapping run."),
        card("Genetic Map (required)", "Marker, Chromosome and bp columns, which place the QTL (in cM) on the genome. With no map, choose “No map: use a VCF”."),
        card("Annotation File (GFF3)", "Gene annotation; the genes inside each interval are read from its gene lines."),
        card("Reference Genome FASTA (optional)", "Adds the sequence of every gene."),
        card("Options", "For a QTL table: confidence interval or flanking markers; and whether a gene must lie completely inside the region."),
    ])
    yield from a.run_and_wait(120)
    yield from a.shot("gene_results", title="Genes Extraction — genes in the region", height=1100, callouts=[
        dict(sel=".table-wrap", note="Each QTL region (with the basis of its interval) and the number of genes in it."),
    ])
    yield from a.tab("Genes in Gene Regions")
    yield from a.shot("gene_results_genes", title="Genes Extraction — genes in gene regions", height=1100, callouts=[
        dict(sel=".table-wrap", note="The genes with their annotation details; the buttons underneath download the gene ID list, the table and the sequences."),
    ])


# ───────────────────────── 4b. Gene prioritization ─────────────────────────
def sc_geneprio(a):
    yield from a.route("gene_prioritization", height=1500)
    yield from a.text_input("Trait", "Drought tolerance")
    yield from a.browse(0, S / "gene_prioritization" / "candidate_genes_rice.txt")
    yield from a.shot("geneprio_input", title="Gene Prioritization — species, trait and inputs", height=1800, callouts=[
        card("Species and trait", "Species (rice, wheat or maize) and the trait in words; related terms are added automatically."),
        card("Gene ID list", "The genes to rank: paste IDs or choose a TXT / CSV / Excel file, for example the Gene_ID_list.txt of Genes Extraction."),
        card("Options", "An option to ignore the cached answers and ask the databases again."),
        card("Database connections", "Ensembl Plants, Gramene, Plant Reactome, NCBI and (rice) Oryzabase: press Check connections to test them."),
    ])
    yield from a.run_and_wait(900, settle=3500)
    yield from a.shot("geneprio_summary", title="Gene Prioritization — summary", height=1500, callouts=[
        dict(sel=".summary-grid", note="How many genes were read, annotated and linked to the trait."),
        dict(text="Top genes", sel="h3", note="The top genes with their final score from the databases."),
        dict(text="Download Ranked List (CSV)", sel="button", note="Downloads: ranked list (CSV), workbook with all evidence (Excel) and report (HTML)."),
    ])
    yield from a.tab("Ranked Genes")
    yield from a.shot("geneprio_ranking", title="Gene Prioritization — ranked genes", height=1200, callouts=[
        dict(sel=".table-tools", note="Search by gene, symbol or description."),
        dict(sel=".table-wrap", note="Click a column to sort, a row to open the evidence. A dash means no evidence of that kind could be retrieved; it is not a score of zero."),
    ])
    yield from a.click(sel="tbody tr", text="Os02g0766700")
    yield 1500
    yield from a.shot("geneprio_evidence", title="Gene Prioritization — evidence behind one gene", height=1700, callouts=[
        dict(sel=".card h3", note="The gene, its rank, position and source links."),
        dict(sel=".table-wrap", note="Score of every evidence category and how many points it added to the final score."),
    ])


# ───────────────────────── 5. PYBSA ─────────────────────────
def sc_pybsa(a):
    yield from a.route("pybsaseq", height=1350)
    yield from a.browse(0, S / "pybsaseq" / "parents_example.tsv")
    yield from a.browse(1, S / "pybsaseq" / "bulks_example.tsv")
    yield from a.text_input("Bulk sizes", "30,30")
    yield from a.text_input("Number of replications", "1000")
    yield from a.text_input("Sliding window", "2000000,100000")
    yield from a.shot("pybsa_input", title="PYBSA G-statistics — input and settings", callouts=[
        card("Parents File", "Parents file: GATK4 table (tab or comma separated) with CHROM, POS, REF, ALT and GT / AD / GQ columns for the two parents."),
        card("Bulks File", "Bulks file: the same layout for the two bulks (columns named <bulk>.GT, <bulk>.AD, <bulk>.GQ)."),
        card("Population & Bulk Sizes", "Number of individuals in each bulk, population type (F2, RIL, BC), p-value cut-offs and simulation replications."),
        card("Windows & Smoothing", "Sliding-window size and step (bp), smoothing, read length."),
        card("Region & Misc", "Optional region of interest, GQ cut-off and other filters."),
    ])
    status = yield from a.run_and_wait(900)
    yield from a.shot("pybsa_results", title="PYBSA G-statistics — results", height=1300, callouts=[
        dict(sel=".subtab-btn", nth=0, parent=1, note="Sub-tabs: the HTML report, plots, every output file and the log."),
    ])


# ───────────────────────── 6. QTL-seq ─────────────────────────
def sc_qtlseq(a):
    yield from a.route("qtl_seq", height=1500)
    yield from a.shot("qtlseq_reads", title="QTL-Seq — starting from BAM files", callouts=[
        card("Reference and bulks", "Reference genome (FASTA) and the BAM files of the high and the low bulk."),
        card("Parents (optional)", "Two parents give the strictest filter; one parent or none gives a weaker analysis, and the page says so."),
        card("Experiment", "Population type, plants per bulk and the allele the SNP-index counts."),
    ])
    yield from a.click(text="VCF files", sel=".chip")
    yield from a.browse(0, S / "qtl_seq" / "qtlseq_example.vcf.gz", wait=1200)
    yield from a.browse(1, S / "qtl_seq" / "qtlseq_example.vcf.gz", wait=1200)
    yield from a.tab("Settings")
    yield from a.click(sel=".switch", nth=1)
    yield from a.text_input("High bulk sample", "HighBulk")
    yield from a.text_input("Low bulk sample", "LowBulk")
    yield from a.shot("qtlseq_vcf", title="QTL-Seq — settings for a VCF", height=2100, callouts=[
        card("SNP filters", "Minimum read depth (7) and minimum SNP-index (0.3) of at least one bulk."),
        card("Sliding window", "Window size (2 Mb), step (10 kb) and the minimum number of SNPs in a window (10)."),
        card("Confidence limits (simulation)", "Replicates, random seed and how the limits of a window are formed."),
        card("VCF sample names", "Needed when one VCF holds both bulks: the sample column of each bulk."),
    ])
    yield from a.run_and_wait(240)
    yield from a.shot("qtlseq_results", title="QTL-Seq — results", height=1900, callouts=[
        dict(sel=".summary-grid", note="Informative SNPs, QTLs and how many reach the 99 %, 95 % and 90 % limits."),
        dict(text="Genome-wide plots", sel="h3", note="Δ SNP-index along the genome with the sliding mean, the three limits and the QTL intervals."),
        dict(text="Significant QTLs", sel="h3", note="One row per QTL. A click opens its chromosome with the genes that overlap it."),
    ])


# ───────────────────────── 7. Diversity ─────────────────────────
def sc_diversity(a):
    yield from a.route("diversity", height=1150)
    yield from a.browse(0, S / "diversity" / "haplotype_blocks_example.csv")
    yield from a.shot("diversity_input", title="Diversity Analysis — input and settings", callouts=[
        card("Input Matrix", "CSV: column 1 = Variety, column 2 = Sample ID, columns 3+ = haplotype block codes (−1 or blank = missing)."),
        card("Number of Clusters", "Groups for hierarchical clustering."),
        card("Top Variable Blocks", "How many of the most informative blocks to show."),
        card("Radar Chart", "Plot the variety radar on a 0–1 scale or with real mean values."),
    ])
    yield from a.run_and_wait(180)
    yield from a.shot("diversity_results", title="Diversity Analysis — results summary", height=1300, callouts=[
        dict(sel=".summary-box", note="Summary of the dataset and of every diversity metric (mean ± standard deviation, range), per variety and per cluster."),
    ])
    yield from a.click(text="PCA", sel=".subtab-btn")
    yield 3500
    yield from a.shot("diversity_plots", title="Diversity Analysis — charts", height=1300,
                      focus=dict(sel=".subtab-btn", nth=0, parent=1), callouts=[
        dict(sel=".subtab-btn", nth=0, parent=1, note="Sub-tabs: block analysis, metrics heatmap, variety comparison, PCA & clustering, NJ tree and the results table."),
        dict(sel=".result-card", nth=0, note="Every chart is a PNG on a white background; click it to open it in a pop-up, or use Download PNG. Here: PCA coloured by variety, and by cluster."),
    ])
    yield from a.click(text="NJ Tree", sel=".subtab-btn")
    yield 2500
    yield from a.shot("diversity_tree", title="Diversity Analysis — circular neighbour-joining tree", height=1500,
                      focus=dict(sel=".subtab-btn", nth=0, parent=1), callouts=[
        dict(sel=".result-card", nth=0, note="Circular tree: pale sectors mark the varieties, concentric circles give the genetic distance, red triangles are branches supported by at least 70% of bootstrap trees. Click to enlarge."),
    ])
    yield from a.click(sel=".result-card .img-wrap")
    yield 1500
    yield from a.shot("diversity_popup", title="Diversity Analysis — the tree opened in a pop-up", height=1500, callouts=[])


# ───────────────────────── 8. GWAS ─────────────────────────
def sc_gwas(a):
    yield from a.route("gwas", height=1450)
    yield from a.browse(0, S / "gwas" / "rice3k_phenotype_example.csv", wait=2200)
    yield from a.shot("gwas_3k", title="GWAS — rice 3K panel (phenotype file only)", callouts=[
        label("Use", "Genotype source: the built-in rice 3K panel, your own VCF, or a ready-made GWAS file."),
        dict(sel=".file-drop", nth=0, note="Phenotype file: first column = accession ID (B001 …) or accession name, other columns = traits."),
        label("Trait to analyse", "Trait: pick the column to analyse."),
        dict(text="Download ID template", sel="button", note="Downloads all 3,024 accession IDs (with name, subpopulation, country) and an empty trait column to fill in."),
    ])
    yield from a.select("Use", "vcf")
    yield from a.browse(0, S / "gwas" / "gwas_example.vcf.gz", wait=900, sel=".file-drop button")
    yield from a.browse(1, S / "gwas" / "gwas_vcf_phenotype_example.csv", wait=2200)
    yield from a.shot("gwas_vcf", title="GWAS — your own genotypes as a VCF", callouts=[
        dict(sel=".file-drop", nth=0, note="VCF (.vcf / .vcf.gz) with genotype calls (GT). Only biallelic SNPs on numbered chromosomes are used."),
        dict(sel=".file-drop", nth=1, note="Phenotype file: first column = sample ID exactly as in the VCF header."),
        label("Trait to analyse", "Trait to analyse."),
    ])
    yield from a.select("Use", "3k")
    yield from a.browse(0, S / "gwas" / "rice3k_phenotype_example.csv", wait=2200)
    yield from a.slider("Permutations", 30)
    yield from a.run_and_wait(600)
    yield from a.shot("gwas_results", title="GWAS — results tables", height=1450, callouts=[
        dict(text="GWAS-ready input", sel="button", note="Downloads: Excel report, MTA table, the GWAS-ready file built from your phenotype data, and all plots."),
        dict(text="Marker–trait associations", sel=".plot-title", note="MTA table: every marker that passes the permutation LOD threshold, per model."),
        dict(text="Model summary", sel=".plot-title", note="One row per model, with genomic inflation (lambda)."),
    ])
    yield from a.shot("gwas_plots", title="GWAS — plots", height=1500, focus=dict(text="Plots", sel=".plot-title"), callouts=[
        dict(sel=".result-card", nth=0, note="Manhattan plot: LOD per marker; diamonds are MTAs."),
        dict(sel=".result-card", nth=1, note="QQ plot: points on the diagonal mean a well-calibrated test."),
    ])


# ───────────────────────── 9. Clustering ─────────────────────────
def sc_clustering(a):
    yield from a.route("clustering", height=950)
    yield from a.browse(0, S / "clustering" / "clustering_genotypes_example.csv")
    yield from a.shot("clustering_input", title="Clustering & Genomic Selection — input", callouts=[
        card("Genotype File", "CSV or TSV: sample IDs in the first column, one column per SNP (0 / 1 / 2, −1 = missing)."),
        card("Number of Clusters", "K for K-means clustering."),
        card("Test Size", "Fraction of each cluster held out as the testing set."),
    ])
    yield from a.run_and_wait(120)
    yield from a.shot("clustering_results", title="Clustering & Genomic Selection — results", height=1350, callouts=[
        dict(sel=".summary-grid", note="Samples, markers, clusters, train / test sizes and variance explained by PC1 and PC2."),
        dict(sel=".result-card", nth=0, note="PCA plot coloured by K-means cluster."),
        dict(text="Training Data", sel="button", note="Downloads: training and testing genotype files (cluster-stratified), cluster membership and the plot."),
    ])


# ───────────────────────── 10. Additive Dominance ─────────────────────────
def sc_additive(a):
    yield from a.route("additive_dominance", height=950)
    yield from a.browse(0, S / "additive_dominance" / "snp_dosage_example.csv")
    yield from a.shot("additive_input", title="Additive Dominance — input", callouts=[
        card("SNP File", "CSV or TSV: sample IDs in the first column, SNP markers after (0 / 1 / 2). Blank cells are missing."),
        card("Imputation", "How missing genotypes are filled."),
        card("Matrix Type", "Additive, dominance or both relationship matrices."),
    ])
    yield from a.run_and_wait(120)
    yield from a.shot("additive_results", title="Additive Dominance — results", height=1300, callouts=[
        dict(sel=".summary-grid", note="Samples, SNPs, missing data and markers with MAF > 0.05."),
        dict(sel=".result-card", nth=0, note="Heatmap of the additive relationship matrix (clusters of related lines show as blocks)."),
        dict(sel=".subtab-btn", nth=0, parent=1, note="Sub-tabs: Heatmaps, Histogram, Data Table and Downloads (CSV of each matrix)."),
    ])


# ───────────────────────── 11. BreedSight ─────────────────────────
def sc_breedsight(a):
    yield from a.route("genomic_prediction", height=1650)
    for i, f in enumerate(["train_phenotypes_example.csv", "train_genotypes_example.csv", "test_phenotypes_example.csv", "test_genotypes_example.csv"]):
        yield from a.browse(i, S / "genomic_prediction" / f, wait=900)
    yield from a.slider("Number of CV folds", 2)
    yield from a.slider("Training epochs", 50)
    yield from a.shot("breedsight_input", title="BreedSight Genomic Prediction — input and settings", callouts=[
        card("Training Phenotypic Data", "Training phenotypes: sample_id, phenotypes."),
        card("Training Genotypic Data", "Training genotypes: sample_id, then SNP columns (0 / 1 / 2), no missing values."),
        card("Testing Phenotypic Data", "Testing phenotypes: sample_id and an optional phenotypes column."),
        card("Testing Genotypic Data", "Testing genotypes: same SNP columns as the training file."),
        card("Cross-Validation & Heritability", "Number of folds and the trait heritability used in the interpretation."),
    ])
    yield from a.run_and_wait(1200)
    yield from a.shot("breedsight_results", title="BreedSight Genomic Prediction — results", height=1350, callouts=[
        dict(sel="table", nth=0, note="Accuracy of every model (BreedSight DNN+RF, Lasso, GBLUP_1, RRBLUP_1) on training and validation data."),
        dict(sel=".subtab-btn", nth=0, parent=1, note="Sub-tabs: metrics and interpretation, test predictions, plots, downloads."),
    ])


# ───────────────────────── 12. Cross selection ─────────────────────────
def sc_cross(a):
    yield from a.route("cross_selection", height=1700)
    yield from a.browse(0, S / "cross_selection" / "dosage_matrix_example.csv", wait=800)
    yield from a.browse(1, S / "cross_selection" / "marker_effects_example.csv", wait=800)
    yield from a.browse(2, S / "cross_selection" / "genetic_map_example.csv", wait=800)
    yield from a.click(text="Preview & Validate", sel="button")
    yield 2500
    yield from a.shot("cross_markers_input", title="Optimal Cross Selection — input files", height=1700, callouts=[
        dict(sel=".file-drop", nth=0, note="Allele dosage matrix: genotypes in the rows, markers in the columns, values 0 / 1 / 2. The relationship matrix is calculated from it."),
        dict(sel=".file-drop", nth=1, note="Marker effects: the columns Marker and Effect."),
        dict(sel=".file-drop", nth=2, note="Genetic map: Marker, Chromosome and cM, to predict the variance of the progeny."),
        card("Optimal Cross Selection", "Number of crosses, penalty degree, and the population type, generation, selected proportion and candidate pool."),
    ])
    yield from a.run_and_wait(300)
    yield from a.shot("cross_markers_results", title="Optimal Cross Selection — usefulness of the crosses", height=1600, callouts=[
        dict(sel=".summary-grid", note="Crosses, parents, the expected usefulness of the plan, group coancestry and progeny inbreeding."),
        dict(text="Crosses to make", sel="h3", closest=".card", note="Each cross with its usefulness, mean and progeny standard deviation."),
    ])

    # no genetic map: a VCF and a recombination rate instead
    yield from a.route("cross_selection", height=1500)
    yield from a.tab("Configure")
    yield from a.click(sel=".chip", text="No map")
    yield 500
    yield from a.browse(2, S / "cross_selection" / "markers_example.vcf", wait=800)
    yield from a.text_input("Recombination rate (cM per Mb)", "4")
    yield from a.click(text="Preview & Validate", sel="button")
    yield 2500
    yield from a.shot("cross_markers_vcf", title="Optimal Cross Selection — a VCF instead of a genetic map", height=1700, callouts=[
        dict(text="No map: use a VCF", sel=".chip", note="Choose this when you have no genetic map."),
        dict(sel=".file-drop", nth=2, note="The VCF: its CHROM, POS and ID columns give the position of every marker."),
        dict(text="Recombination rate (cM per Mb)", sel="label", closest=".field", note="The rate that turns base pairs into cM, the same for the whole genome."),
    ])
    yield from a.run_and_wait(300)
    yield from a.shot("cross_markers_vcf_results", title="Optimal Cross Selection — plan from a VCF", height=1400, callouts=[
        dict(sel=".summary-grid", note="The plan is calculated in the same way; the note above says the positions came from the VCF."),
    ])


# ───────────────────────── 13. Heterosis ─────────────────────────
def sc_heterosis(a):
    yield from a.route("heterosis", height=1250)
    for i, f in enumerate(["train_phenotype_example.csv", "train_genotype_example.csv", "parent_genotype_example.csv", "test_genotype_example.csv"]):
        yield from a.browse(i, S / "heterosis" / f, wait=1200)
    yield from a.slider("K-Fold CV splits", 3)
    yield from a.shot("heterosis_input", title="Heterosis — input files and parameters", callouts=[
        dict(sel=".file-drop", nth=0, note="Training phenotype CSV: column 0 = sample ID, column 1 = numeric trait value."),
        dict(sel=".file-drop", nth=1, note="Training genotype CSV: first column = sample ID, then numeric marker columns. Row order must match the phenotype file."),
        dict(sel=".file-drop", nth=2, note="Parent genotype CSV: first column = parent ID, same marker columns as the training file."),
        dict(sel=".file-drop", nth=3, note="Test genotype CSV (optional): same marker columns, for out-of-sample predictions."),
        card("Parameters", "Number of hybrids / parents to report, cross-validation folds and network settings."),
    ])
    yield from a.run_and_wait(900)
    yield from a.shot("heterosis_results", title="Heterosis — results", height=1700, callouts=[
        dict(text="Downloads", sel="h3", closest=".card", note="Top F1 hybrids, top parents, all hybrid results and test predictions as CSV."),
        dict(text="Cross-Validation Summary", sel="h3", closest=".card", note="Prediction accuracy from cross-validation (MSE, RMSE, R², Pearson r)."),
        dict(text="Top F1 Hybrid", sel="h3", closest=".card", note="Best hybrid combinations with GEBV, mid-parent (MPH) and better-parent (BPH) heterosis."),
    ])


# ───────────────────────── 14. Phenotypic ─────────────────────────
def sc_phenotypic(a):
    yield from a.route("phenotypic", height=1300)
    yield from a.browse(0, S / "phenotypic" / "multi_environment_trial_example.csv", wait=3000)
    for _ in range(60):                      # the file is read by R in the background: wait for the summary cards
        ready = yield from a.call("!!document.querySelector('.view.active .summary-card')")
        if ready:
            break
        yield 1000
    yield 1500
    yield from a.select("Design of the trial", "rcbd")     # this trial has every genotype in every replicate
    yield 1500
    for _ in range(60):                      # the file is read again for the new design
        busy = yield from a.call("document.body.innerText.includes('Reading the data file')")
        if not busy:
            break
        yield 1000
    yield 1500
    yield from a.shot("phenotypic_input", title="Phenotypic Data Analysis — data and design", callouts=[
        card("Phenotypic Data File", "CSV or Excel, one row per plot. Genotype, Environment, Rep (and Block) columns are found by name; every other numeric column is a trait. The file is read once and shared by every analysis."),
        card("Experimental Design", "Augmented RCBD, RCBD or alpha-lattice: decides the blocking terms of the models."),
        card("Environments", "Analyse all environments or only chosen ones; 0 can count as missing."),
        card("Design Check", "Layout of every environment: plots, accessions, replications, blocks and checks."),
    ])
    yield from a.tab("Analyses")
    yield from a.run_and_wait(300)
    yield from a.shot("phenotypic_results", title="Phenotypic Data Analysis — an analysis tab", height=1700, callouts=[
        dict(sel=".subtab-btn", nth=0, parent=1, note="One tab per analysis: descriptives, environments, BLUEs, ANOVA, heritability, correlation, PCA, clustering and stability."),
        dict(text="Traits to analyse", sel="h3", closest=".card", note="Each analysis has its own traits and options, and its own Run button."),
    ])
    yield from a.click(text="ANOVA", sel=".subtab-btn")
    yield 800
    yield from a.run_and_wait(400)
    yield from a.shot("phenotypic_anova", title="Phenotypic Data Analysis — ANOVA", height=1500, callouts=[
        dict(sel=".anova-tbl", nth=0, note="One table per trait: every environment and the pooled analysis side by side (DF, SS, MS, F, P, significance), with plots, mean and CV underneath."),
    ])
    yield from a.click(text="Heritability", sel=".subtab-btn")
    yield 800
    yield from a.run_and_wait(300)
    yield from a.shot("phenotypic_heritability", title="Phenotypic Data Analysis — heritability", height=1500, callouts=[
        dict(sel=".plot-title", nth=0, note="Variance components and heritability (entry-mean and plot basis) per trait and environment."),
    ])


# ───────────────────────── Genotype Viewer ─────────────────────────
def sc_genotype(a):
    yield from a.route("genotype_viewer", height=2000)
    yield from a.browse(0, S / "genotype_viewer" / "rice_RILs_simple_example.csv", wait=1200)
    yield from a.click(sel=".chip", text="3-color")
    yield from a.click(sel=".chip", text="Red-Blue-Yellow-Grey")
    yield from a.text_input("Sample names, separated by commas", "Parent1, Parent2, RIL01, RIL02")
    yield from a.text_input("REF (first colour)", "Parent 1")
    yield from a.text_input("ALT (second colour)", "Parent 2")
    yield from a.text_input("HET (heterozygous)", "Heterozygous")
    yield from a.shot("genotype_input", title="Genotype Viewer — file, species, drawing mode, colours and image quality", height=2200, callouts=[
        dict(sel=".file-drop", nth=0, note="Genotype file: CSV or VCF. The first column is the chromosome, then the marker name, the position (bp) and one column per sample; the column names do not matter."),
        card("Species", "Species whose chromosome lengths are used. The chromosome names of the file are matched to its names (1, chr1 and chr01 are the same)."),
        card("Samples to draw (optional)", "Name some samples to draw only those; empty draws every sample."),
        card("Drawing mode", "Normal: one image per sample. Compare: one image per chromosome with all samples side by side. Zoomed: one region of one chromosome for all samples."),
        card("Colouring", "2 or 3 colours, fill between markers, marker names (and, in compare / zoomed mode, the dashed line from each name)."),
        card("Colour palette", "GenoSee's palettes (bright, Aqua, grays, reds, blues) and a mixed red-blue-yellow-grey palette."),
        card("Legend names (optional)", "Words that replace REF / ALT / HET in the legend, for example the names of the two parents."),
        card("Image quality", "PNG resolution (200, 300 or 600 dpi) and an extra vector PDF for publication."),
    ])
    yield from a.run_and_wait(300)
    yield from a.shot("genotype_results", title="Genotype Viewer — the graphical genotypes", height=1900, callouts=[
        dict(sel=".summary-grid", note="Images drawn, samples, markers and chromosomes."),
        dict(sel=".results-grid", nth=0, note="One image per sample (normal mode): each chromosome is a bar coloured by the genotype; the legend is at the right."),
    ])


# ───────────────────────── Favourable Allele Analysis ─────────────────────────
def sc_favallele(a):
    yield from a.route("favourable_allele", height=1200, wait=2500)
    yield from a.shot("favallele_upload_empty", title="Favourable Allele Analysis — scope and the genotype file", height=1200, callouts=[
        dict(sel=".note", nth=0, note="Scope: the tool is for the rice 4K haplopanel only, for now. The trait-gene markers are kept inside TerraTools; it is not scaled to other panels or crops yet."),
        card("Genotype file", "Your 4K haplopanel genotypes: one row per marker, one column per line; 0 / 1 / 2 / -1, VCF genotypes or nucleotides."),
    ])
    yield from a.browse(0, HERE.parent / "Demo_Input_Files" / "18_Favourable_Allele" / "input" / "rice_4K_haplopanel_demo.csv", wait=2500)
    yield from a.shot("favallele_setup", title="Favourable Allele Analysis — matched marker IDs, filters and lines", height=2600, callouts=[
        dict(sel=".summary-grid", note="Markers in your file, how many of the trait-gene markers were found by ID, the trait genes covered and the number of lines."),
        card("Matched marker IDs", "Every ID of your file that matched, with the gene, the trait and how it was matched. Download matched IDs saves the list as CSV."),
        card("Filters", "Restrict the traits and genes (the list follows the traits) and choose the status to show."),
        card("Select lines to compare", "Search and tick the lines of your file. The coloured names are the selected ones: click a name to make it the line shown in detail, × removes it."),
    ])
    yield from a.tab("Overview", wait=2200)
    yield from a.shot("favallele_overview", title="Favourable Allele Analysis — overview of the selected line", height=1650, callouts=[
        dict(sel=".fa-bar", note="The selected lines (click one to show it in detail) and the filters in force. Change takes you back to the file tab."),
        dict(sel=".fa-stats", note="SNP and gene counts of the line shown: favourable, unfavourable, heterozygous and missing."),
        card("Line profile", "One dot for every gene, placed on the plant part that the trait belongs to. Hover for the gene, click for all its SNPs; colour the dots by trait or by allele status."),
        card("SNPs by trait", "Number of SNPs of the chosen status for every trait, one bar per selected line."),
    ])
    yield from a.js("document.querySelectorAll('.fa-dot')[5].dispatchEvent(new MouseEvent('click', {bubbles: true})); 1")
    yield 1200
    yield from a.shot("favallele_gene", title="Favourable Allele Analysis — the details of one gene", height=1100, callouts=[
        {"sel": ".fa-modal", "global": True, "note": "Trait, position and the status of every SNP of the gene in every selected line (REF/ALT written out as the two nucleotides)."},
    ])
    yield from a.js("document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true})); 1")
    yield 800
    yield from a.tab("Genes", wait=2200)
    yield from a.shot("favallele_genes", title="Favourable Allele Analysis — genes and the allele map", height=2300, callouts=[
        card("Allele status by gene", "A bar per gene for the line shown: favourable, heterozygous, unfavourable and missing SNPs."),
        card("Allele map of the selected lines", "Rows are the selected lines, columns the SNPs of each gene, coloured by status — to see where the lines differ."),
    ])
    yield from a.tab("Heatmap", wait=1500)
    for t in ("DSR trait", "Floral", "Hybrid fertility"):
        yield from a.click(sel=".fa-quickbtn", text=t)
    yield from a.click(sel=".fa-quick-row .btn", text="All lines")
    yield from a.shot("favallele_heatmap_setup", title="Favourable Allele Analysis — Heatmap tab: genes and lines", height=1700, callouts=[
        card("Genes", "Pick genes from the list (search by name), add all the genes of a trait with one click, or take the genes of the filters."),
        card("Lines", "Pick lines from the list, or take all the lines of your file."),
        card("Figure", "Resolution of the PNG."),
    ])


ORDER = [sc_home, sc_tour, sc_genofile, sc_mendelian, sc_linkage, sc_genotype, sc_qtl, sc_gene, sc_geneprio, sc_pybsa, sc_qtlseq, sc_diversity, sc_gwas,
         sc_clustering, sc_additive, sc_breedsight, sc_cross, sc_heterosis, sc_favallele, sc_phenotypic]
SCENARIOS = {f.__name__[3:]: f for f in ORDER}

if __name__ == "__main__":
    pick = sys.argv[1:] or list(SCENARIOS)
    app = App(OUT)
    app.run([SCENARIOS[p] for p in pick])
