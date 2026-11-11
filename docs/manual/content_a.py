"""Manual content, part A: getting started and the trait-mapping tools."""
from mdkit import box, example, figure, ol, sample_files, table, ul

TOOLS = [
    ("genofile", "Genofile Converter", "Data Analysis", "The start of the chain: turns the HapMap, VCF or raw allele files of a genotyping company into the 0 / 1 / 2 input files of Mendelian Segregation and Linkage Mapping."),
    ("phenotypic", "Phenotypic Data Analysis", "Data Analysis", "Multi-environment trial statistics: ANOVA, mixed models, BLUEs, heritability, G×E stability."),
    ("mendelian", "Mendelian Segregation", "Trait Mapping", "Tests every marker for the expected segregation ratio (1:1 or 1:2:1) to find distorted markers."),
    ("linkage", "Linkage Mapping", "Trait Mapping", "Builds genetic linkage maps (Kosambi cM) from coded marker data, with a file ready for QTL Mapping."),
    ("genotype", "Genotype Viewer", "Trait Mapping", "Draws graphical genotypes of every chromosome for each sample (CSV or VCF), with bright or publication-ready palettes."),
    ("qtl", "QTL Mapping", "Trait Mapping", "Locates QTL as regions with left / peak / right markers using six mapping methods; compares methods."),
    ("qtlseq", "QTL-Seq", "Trait Mapping", "Bulk segregant analysis with the Δ(SNP-index) method from reads, BAM files or a VCF."),
    ("diversity", "Diversity Analysis", "Trait Mapping", "Haplotype diversity metrics, clustering, PCA and a neighbour-joining tree for varieties."),
    ("gwas", "GWAS", "Trait Mapping", "Genome-wide association with eight models; built-in rice 3K panel or your own VCF."),
    ("pybsa", "PYBSA G-statistics", "Trait Mapping", "Bulk segregant analysis from GATK variant tables using the G-statistic and Δ(allele frequency)."),
    ("gene", "Genes Extraction", "Trait Mapping", "Turns the QTL of QTL Mapping into base-pair intervals and lists the genes inside them."),
    ("prioritize", "Gene Prioritization", "Trait Mapping", "Ranks candidate genes for a trait from Ensembl Plants, Gramene, Plant Reactome, PubMed and (rice) Oryzabase."),
    ("clustering", "Clustering & Genomic Selection", "Genomic Selection", "PCA and K-means clustering with a cluster-stratified training / testing split."),
    ("additive", "Additive Dominance", "Genomic Selection", "Additive and dominance genomic relationship matrices from SNP dosages."),
    ("breedsight", "BreedSight GUI", "Genomic Selection", "Deep-neural-network + random-forest genomic prediction benchmarked against classical models."),
    ("cross", "Optimal Cross Selection", "Genomic Selection", "Plans the crosses with the highest usefulness while keeping diversity."),
    ("heterosis", "Heterosis", "Genomic Selection", "Deep-learning GEBV prediction and screening of hybrid combinations."),
    ("favallele", "Favourable Allele Analysis", "Genomic Selection", "Matches your rice 4K haplopanel genotype file with 166 known trait genes and shows which favourable alleles each line carries (rice 4K only, for now)."),
]


def front_welcome():
    rows = [[f"{i + 1}", f"<b>{n}</b>", c, d] for i, (_k, n, c, d) in enumerate(TOOLS)]
    return f"""
<p>TerraTools is a desktop application from IRRI-SAH that brings nineteen plant-breeding and genetics analyses into one window.
You open a tool from the sidebar, load your files, press <b>Run</b>, and read the results as plots and tables that you can download.
Everything runs on your own computer, and R and Python are already built in; only Gene Prioritization needs an internet connection, because it asks public gene databases.</p>
<h3>The tools</h3>
{table(["#", "Tool", "Group", "What it does"], rows)}
<h3>What you need</h3>
{ul(["A Windows 10 / 11 (64-bit) PC, or a Mac with Apple Silicon (M1 or later) running macOS 12 or newer.",
     "At least 8 GB of memory; 16 GB is more comfortable for deep-learning tools (BreedSight, Heterosis) and for large VCF files.",
     "Your data as CSV or Excel files in the layouts described in this manual.",
     "About 3 GB of free disk space for the program, plus room for your results."])}
{box("tip", "Every chapter lists the exact file layout the tool expects and points to a small example file in the <code>sample_data</code> folder that accompanies this manual. Run an example first to see how a tool behaves, then replace it with your own data.")}
"""


def front_install():
    return f"""
<h3>Windows</h3>
{ol(["Run the installer (<code>TerraTools_Setup.exe</code>) and follow the prompts; tick <i>Create a desktop shortcut</i> if you want one. It installs TerraTools (by default in <i>Program Files</i>) together with its own R and Python runtimes — nothing else needs to be installed.",
     "Start TerraTools from the Start menu, the desktop shortcut or the <i>Launch TerraTools</i> option on the last installer page. The first start takes a little longer while the program prepares itself.",
     "To remove it later use <b>Settings → Apps → TerraTools → Uninstall</b>."])}
<h3>macOS</h3>
{ol(["Open <code>TerraTools.dmg</code> and drag <b>TerraTools</b> onto the <b>Applications</b> folder.",
     "Open TerraTools from Applications. The very first launch unpacks the built-in R runtime, which takes about a minute; later launches are quick.",
     "If macOS says the app cannot be verified, the build you have was not notarized by Apple. Open <b>System Settings → Privacy &amp; Security</b>, scroll down to the message about TerraTools and click <b>Open Anyway</b> (on macOS 14 and older you can instead right-click the app and choose <b>Open</b>). You only need to do this once."])}
{box("note", "A notarized build opens without any prompt. Which kind you have depends on how the release was published; the disk image contains a <code>READ ME FIRST.txt</code> only when approval is needed.")}
<h3>Where are my results?</h3>
<p>While a tool runs, it writes its working files to a temporary folder. Nothing is lost when you switch tools, but the temporary copies are not meant to be kept:
use the <b>Download</b> buttons (Excel, CSV, PNG, ZIP) to save results where you want them. Each button opens an ordinary <b>Save As</b> window.</p>
"""


def front_tour():
    return f"""
<p>All nineteen tools share the same layout, so once you know one you know them all.</p>
{figure("home")}
<p>The <b>sidebar</b> on the left lists the tools by group. A green dot beside a tool means a run finished there while you were looking at another page;
a spinning dot means it is still running — you can leave a tool while it works and come back later, and a notification appears when it is done.</p>
{figure("interface_tour")}
<h3>Working with plots and tables</h3>
{ul(["<b>Plots:</b> every plot is an image on a white background, shown whole and never squeezed; the plots of a tab are arranged in two columns (one column in a narrow window); click it to open it full size in a pop-up (click outside it or press Esc to close). Every plot has its own <b>Download PNG</b> button, and pages with several plots offer <b>Download all as PNG</b>.",
     "<b>Tables:</b> use the <i>Filter rows…</i> box to search, scroll sideways for more columns, and use the download button next to the table to get the complete data as CSV or Excel.",
     "<b>Status line:</b> under the Run button. Errors are written there in plain words — read them first, they usually say exactly which file or setting to fix."])}
"""


def front_data():
    coding = table(["Tool", "Genotype coding", "Missing values"], [
        ["Genofile Converter (input)", "alleles as sent by the company: HapMap <code>AG</code>, VCF <code>0/1</code>, raw <code>A:G</code>; the output is coded 0 / 1 / 2 / −1", "<code>NN</code>, <code>./.</code>, <code>-</code>, <code>N:N</code> … (becomes <code>-1</code>)"],
        ["Mendelian Segregation", "0 / 1 / 2 (RIL: 0 and 2)", "Codes you list in the settings (default <code>NA</code>, <code>.</code>, <code>-9</code>); add <code>-1</code> if you use it"],
        ["Linkage Mapping", "0 = parent 1, 1 = heterozygous, 2 = parent 2", "<code>-1</code> (or blank)"],
        ["Genotype Viewer", "<code>A</code> / <code>H</code> / <code>B</code>, or <code>0|0</code> / <code>0|1</code> / <code>1|1</code>, or <code>0/0</code> / <code>0/1</code> / <code>1/1</code> (or a VCF)", "<code>N</code>, <code>.|.</code> or <code>./.</code>"],
        ["QTL Mapping", "0 / 1 / 2 (RIL: 0 and 2)", "any negative value (<code>-1</code>) or blank; filled automatically from neighbouring markers"],
        ["GWAS", "0 / 1 / 2 dosage or −1 / 0 / 1 centred, as you choose (a file with only 0 and 1 must be told which); VCF: diploid genotype calls", "<code>NA</code> or blank, <code>-9</code>, <code>-99</code>, <code>-999</code>; other codes can be listed"],
        ["Clustering & Genomic Selection", "0 / 1 / 2", "<code>-1</code>"],
        ["Additive Dominance", "0 / 1 / 2", "blank or <code>NA</code> (filled by mean, median or zero)"],
        ["Diversity Analysis", "haplotype block numbers (0 is a real code)", "<code>-1</code> or blank"],
        ["BreedSight, Heterosis", "numeric marker values (0 / 1 / 2)", "none allowed — impute first"],
    ])
    return f"""
<p>Most problems people meet come from the input files, not from the analyses. These rules apply to every tool.</p>
{ul(["<b>File types.</b> CSV (comma-separated, UTF-8) or Excel <code>.xlsx</code>, as stated for each tool. Save from Excel with <i>File → Save As → CSV UTF-8</i>.",
     "<b>One header row</b> with column names, unless the tool says otherwise (QTL Mapping's genotype file has three header rows).",
     "<b>IDs must match across files, exactly.</b> <code>Line_01</code>, <code>line_01</code> and <code>Line 01</code> are three different names. Avoid leading or trailing spaces. Chromosomes must be written the same way everywhere (<code>3</code> in one file and <code>Chr3</code> in another will not match).",
     "<b>No merged cells, no notes or totals rows, no blank rows</b> inside the data. Put comments in a separate sheet or file.",
     "<b>Decimals use a dot</b> (<code>5.25</code>, not <code>5,25</code>). Do not use thousands separators in positions.",
     "<b>Keep names simple:</b> letters, digits, underscore and hyphen. Very long names are fine but make tables hard to read."])}
<h3>Genotype coding at a glance</h3>
{coding}
{box("warn", "Different tools treat the <i>missing</i> code differently. A <code>-1</code> that is read as a real genotype (for example in Mendelian Segregation if you forget to list it) silently changes the result. Check the table above before loading a file.")}
<h3>Example data</h3>
<p>The <code>sample_data</code> folder that comes with this manual holds one small, simulated example for every tool, in exactly the layout the tool expects.
The screenshots in this manual were made with these files. They are quick to run and a good template: open one in Excel, replace the values with your own, and keep the header and the column order.</p>
"""


# ───────────────────────── tools 1-7 ─────────────────────────
def ch_mendelian():
    return f"""
<p class="lead">Checks whether each marker segregates in the ratio Mendel predicts, and flags the markers that do not (segregation distortion).
Run it before building a linkage map or mapping QTL, so that badly behaved markers can be removed.</p>
<h3>What you need</h3>
<p>One CSV file: a marker per row, the first column the marker name, the second the chromosome, and one column for every individual in the population.</p>
<p>Genotypes as they come from a genotyping company (HapMap, VCF, <code>A:G</code> calls) are converted to this layout by the <b>Genofile Converter</b> (chapter 5), which has a <b>Mendelian Segregation</b> shortcut under <i>Your converted file</i>.</p>
{example("mendelian/mendelian_RIL_example.csv", rows=5, cols=7)}
{ul(["Genotypes are <code>0</code>, <code>1</code> and <code>2</code> (homozygous parent A, heterozygous, homozygous parent B). In a RIL population most calls are 0 or 2.",
     "Every marker with at least one non-missing call is tested. When an <i>expected</i> class count is below 5 the χ² is unreliable, so the marker is flagged <i>low power</i> and counted in the summary.",
     "Missing values are whatever you list in <i>Missing Data</i> (see below).",
     "The file that Linkage Mapping writes for QTL Mapping (markers across, one row per individual, with a <code>Chromosome</code> and a <code>Position</code> row) is accepted too; it is turned into one row per marker automatically."])}
{sample_files("mendelian")}
<h3>Step by step</h3>
{ol(["Open <b>Mendelian Segregation</b> and click <b>Browse</b> to load the genotype CSV.",
     "Choose the population type, the significance threshold (α) and the multiple-testing correction, and check the missing-data codes.",
     "Click <b>Run Segregation Test</b>. The Results tab opens when it finishes."])}
{figure("mendelian_input")}
{table(["Setting", "What it does"], [
    ["Population type", "<b>RIL</b>: expects 1:1 between the two homozygous classes (χ², 1 d.f.). <b>F2</b>: expects 1 : 2 : 1 (χ², 2 d.f.)."],
    ["Significance threshold (α)", "A marker is called <i>distorted</i> when its <b>corrected</b> p-value (<code>p_adj</code>) is below α."],
    ["Multiple-testing correction", "<b>FDR (Benjamini-Hochberg)</b>, the default, <b>Bonferroni</b>, or <b>None</b> (raw p). Testing hundreds of markers at a raw α of 0.05 flags many by chance alone, so a correction is recommended. It is applied over all markers tested, before unplaced markers are dropped."],
    ["Missing data", "Comma-separated codes counted as missing. The default is <code>NA,.,-9,-1</code>."],
    ["Keep only markers with a known chromosome", "Drops markers without a chromosome (unplaced markers) from the results."]])}
<h3>Reading the results</h3>
{figure("mendelian_results")}
{ul(["<b>−log10(p) by marker:</b> every point is a marker; red points are distorted after the correction. With Bonferroni or no correction a dashed line marks the cut-off; FDR has no single line, so the colour carries the call. A cluster of red points in one place usually means a genuinely distorted region, not several bad markers.",
     "<b>Distortion rate by chromosome:</b> the share of distorted markers on each chromosome.",
     "<b>Summary:</b> markers tested, flagged distorted, passing the quality check, and the number with an expected count below 5 (read their p-values with caution).",
     "<b>Full Results (CSV):</b> every marker with its counts (<code>c0, c1, c2</code>), χ², p, <code>p_adj</code>, <code>distorted</code>, <code>low_power</code>, heterozygosity, −log10 p and the genotype calls.",
     "<b>Filtered Markers (CSV):</b> the markers that are <i>not</i> distorted, i.e. that pass the quality check, with their genotypes.",
     "<b>LM-ready File (XLSX):</b> the same markers as just Marker, Chromosome and the genotype calls, ready to load in Linkage Mapping."])}
{box("tip", "A distorted region is not always a problem: it can point to a gametic selection factor or a sterility locus. Compare with the map position before deleting a whole stretch of markers.")}
"""


def ch_linkage():
    return f"""
<p class="lead">Builds a genetic linkage map the way JoinMap and Lep-MAP3 do: it checks which markers belong together (LOD), finds the marker order that the genotypes support, and estimates all distances together (multipoint maximum likelihood) in centiMorgans.
It writes a file that QTL Mapping can read directly. Genotypes from a genotyping company are converted to the 0 / 1 / 2 codes first, by the <b>Genofile Converter</b> (chapter 5).</p>
<h3>What you need</h3>
<p>An Excel file (<code>.xlsx</code>) with one row per marker: <b>Marker</b>, <b>Chromosome</b>, then one column per individual, coded <code>0</code> (homozygous for parent 1),
<code>1</code> (heterozygous), <code>2</code> (homozygous for parent 2) or <code>-1</code> (missing).</p>
{example("linkage/linkage_F2_coded_example.xlsx", rows=5, cols=7)}
{box("key", "<b>How the order is decided.</b> Markers with identical genotypes share one position. The order is searched from the genotypes (smallest recombination between neighbours, then improved with the multipoint likelihood). If the marker names carry physical positions (<code>chr02_16919033</code>) that order is kept unless the genotypes favour another order by more than 3 LOD; markers without a position in their name are placed where they fit best. Markers on chromosome 0 (unplaced) are left out.")}
{box("note", "Genotypes as they come from a genotyping company (HapMap, VCF, or <code>A:G</code> calls) are not read by this tool. Convert them with the <b>Genofile Converter</b> (chapter 5) first: its <b>Linkage Mapping</b> shortcut under <i>Your converted file</i> loads the converted workbook straight into the Upload tab. If the marker names carry the physical position (<code>chr02_16919033</code>), choose <i>Marker names: chromosome_position</i> there.")}
{sample_files("linkage")}
<h3>Step 1: upload and set parameters</h3>
{ol(["On the <b>Upload</b> tab, browse to the coded Excel or CSV file (the file written for QTL Mapping can be used again here) and choose the population type (<b>F2</b> or <b>RIL</b>). A preview and the number of markers per chromosome appear.",
     "On the <b>Parameters</b> tab set the filters (see the table) and click <b>Run Analysis</b>."])}
{figure("linkage_upload")}
{figure("linkage_params")}
{table(["Parameter", "What it does", "Typical value"], [
    ["Missing data threshold", "Markers with a larger share of missing calls are removed.", "0.30"],
    ["Max recombination fraction (r)", "Neighbouring markers with r at or above this value are not considered linked and get a 50 cM gap (the QC report warns when this happens).", "0.45"],
    ["Marker order", "The starting / reference order: automatic (re-sorts a chromosome by the physical position in the marker names only if the file lists it alphabetically), by physical position, or as in the file.", "automatic"],
    ["Map method", "<b>Find the best marker order</b> (recommended) or <b>keep my marker order</b>. Both estimate Kosambi distances with all markers of a chromosome together and allow for 1 % genotyping error, so a single wrong call does not inflate the map. Markers linked more strongly to another chromosome, not linked to their own, or that add more than 10 cM to the map are removed and listed in the Removal Log.", "find best order"],
    ["Segregation distortion filter", "Optional χ² filter that removes distorted markers (with optional Bonferroni correction). Monomorphic markers are always removed.", "off"],
    ])}
<h3>Reading the results</h3>
{ul(["<b>QC Report:</b> how many markers were removed and why (missing data, monomorphic, distortion, chromosome 0); the full list can be downloaded.",
     "<b>Results:</b> the linkage-map table (marker, chromosome, position in cM), retained markers, all markers with their status, a chromosome summary (length, spacing, largest gap, log-likelihood, which order was used) and the Marker Fit table (suspect genotypes per map position — calls that are improbable given their neighbours).",
     "<b>File for QTL Mapping:</b> the linkage map joined to your genotypes by marker ID and transposed — one row per individual, with marker, chromosome and position as the first three rows. Load it directly as the genotype file in <b>QTL Mapping</b>.",
     "<b>Map Plot, Density &amp; RF, Cytogenetic:</b> the map figure, marker density per bin, the recombination heat map of adjacent markers, and an interactive ideogram that shows the allele composition along each chromosome."])}
{figure("linkage_results")}
{figure("linkage_map")}
{box("tip", "Choose chromosomes, theme and sizes in <b>Plot Options</b> and click <b>Update Plot</b> — changes apply to the finished analysis without re-running it.")}
"""


def ch_qtl():
    return f"""
<p class="lead">Finds quantitative trait loci. Each QTL is reported as a <b>region</b>: a left flanking marker, a peak marker and a right flanking marker, with the
LOD, p-value and effect at the peak — following the usual interval-mapping practice (Lander &amp; Botstein, Haley &amp; Knott, Zeng, Churchill &amp; Doerge, Broman &amp; Sen).</p>
<h3>What you need</h3>
<p><b>1. Genotype CSV</b> — three header rows, then one row per individual. Row 1 holds the marker names, row 2 the chromosome of each marker and row 3 its position in cM.
(The linkage-mapping tool writes exactly this file for you.)</p>
{example("qtl/qtl_genotypes_example.csv", rows=6, cols=7, header_rows=1, caption="Genotype file (first column and first markers shown). Rows 1–3 are marker, chromosome and position; each following row is one individual.")}
<p><b>2. Phenotype CSV</b> — one row per individual: the first column is the ID (the same IDs as in the genotype file), then one column per trait.</p>
{example("qtl/qtl_phenotypes_example.csv", rows=5, cols=4)}
{ul(["At least 8 individuals must appear in both files; individuals present in only one file are ignored.",
     "Genotypes are <code>0</code>, <code>1</code>, <code>2</code> (RIL: 0 and 2); negative values and blanks are missing, and are filled in from the neighbouring markers.",
     "Markers must be sorted by position inside each chromosome."])}
{sample_files("qtl")}
<h3>Step by step</h3>
{ol(["Open <b>QTL Mapping</b>. Load the genotype CSV, then the phenotype CSV. The <b>Trait</b> list fills with the phenotype columns — choose one.",
     "Pick the mapping method and the threshold (see below). Leave the other settings at their defaults for a first run.",
     "Click <b>Run QTL Analysis</b>. The Results tab shows the Manhattan plot and the QTL table."])}
{figure("qtl_input")}
{table(["Setting", "What it does"], [
    ["Mapping method", "<b>SIM</b> — one marker at a time. <b>CIM</b> — adds a few background markers (cofactors) to remove the effect of other QTL; the usual first choice. <b>MQM</b> — a larger background model. <b>RF, SVR, XGBoost</b> — machine-learning importance scores on a LOD-like scale (not true LODs)."],
    ["Significance threshold", "<b>Fixed LOD</b> (default 2.5), or a <b>genome-wide permutation test</b> (SIM / CIM / MQM): the phenotype is shuffled many times and the 95th percentile of the genome-wide maximum LOD becomes the threshold. Permutation is the safer choice for publications."],
    ["Peak p-value cut-off", "A peak must also have a nominal p-value below this (default 0.05)."],
    ["Peak separation (LOD drop)", "Two peaks on one chromosome count as separate QTL only if the LOD curve dips by at least this much between them (default 2)."],
    ["Support interval — LOD drop", "The QTL region extends left and right of the peak until the LOD has fallen by this amount; 1.5 gives roughly a 95 % interval (default 1.5). The first marker outside on each side is the flanking marker."],
    ["CIM cofactors / window", "Maximum number of background markers and the distance (cM) around the tested position in which cofactors are ignored."],
    ["Hold-out fraction", "Fraction of individuals set aside before mapping. Leave at 0 to use everyone."]])}
<h3>Reading the results</h3>
{figure("qtl_results")}
<p>Each panel is one chromosome: the curve is the LOD along the chromosome, the dashed line the threshold, the red diamond the QTL peak, and the shaded band the QTL region.</p>
{figure("qtl_table")}
{table(["Column", "Meaning"], [
    ["QTL_ID", "Name of the QTL (chromosome and order along it)."],
    ["Left_Marker, Peak_Marker, Right_Marker", "The flanking markers and the marker with the highest LOD. The QTL lies between the flanking markers."],
    ["Interval_Size_cM, CI_Start/End", "Width of the region in cM, and the interpolated positions where the LOD curve crosses the support-interval level."],
    ["LOD, P_Value, Genomewide_P", "Strength of the evidence at the peak. <i>Genomewide_P</i> appears when the permutation threshold is used."],
    ["PVE_Marker_Percent", "Percent of the trait variance explained by the peak marker."],
    ["Additive_Effect, SE_Additive", "Half the difference between the two homozygous classes, in trait units. Positive means the higher-coded allele (genotype 2) raises the trait. For CIM / MQM it is estimated with the background markers in the model."],
    ["Dominance_Effect, Dominance_Class", "Departure of the heterozygote from the mid-point (F2 only) and a label: additive, partial dominance, dominance, overdominance."]])}
<h3>QTL interval file</h3>
<p>Besides the QTL table, every run offers <b>QTL Intervals CSV</b>: one row per QTL with three columns — <code>QTL_ID</code>, <code>QTL_Start</code> and <code>QTL_End</code>.
The interval is the <i>CI_Start / CI_End</i> support interval. Because the map is in cM and a genome annotation in base pairs, the two ends are converted to bp by linear interpolation between the neighbouring markers, using the physical position written in the marker names (e.g. <code>chr03_5000000</code>). If you give the <b>Linkage map</b> (optional, under Input Files) — the <i>Retained Markers</i> or <i>Linkage Map</i> Excel file saved from Linkage Mapping — its marker positions are used for this conversion. It must be the map the genotype file was built from: if its cM positions differ from the genotype file's, it is ignored and a note says so. A QTL on a chromosome with fewer than two such markers cannot be converted and is listed under the download buttons instead.
This file is the input of <b>Genes Extraction</b> (<i>Use QTL Mapping output</i>).</p>
<h3>Compare Models and Common QTL</h3>
<p>The <b>Compare Models</b> tab runs several methods on the same data and lines their QTL regions up; two methods found “the same QTL” when their regions overlap.
<b>Common QTL</b> then keeps the regions that at least <i>n</i> methods agree on. Both tables can be downloaded as CSV.</p>
{box("tip", "The quickest way to a convincing result: run CIM with the permutation threshold, then compare CIM, SIM and MQM. QTL found by all three are the ones to trust.")}
"""


def ch_gene():
    return f"""
<p class="lead">Takes the QTL found by <b>QTL Mapping</b> (or any table of QTL intervals), finds the genes of a GFF3 annotation in them, and lists each gene with the details of the annotation (name, function) and, if you give the reference genome, its sequence. The result can be downloaded as a table, as a <b>gene ID list</b> and as a FASTA file.</p>
<h3>What you need</h3>
<p><b>1. QTL results</b> — a CSV, TSV or Excel table. The columns are found by their <b>name</b>, in any order. Genes Extraction takes QTL results only: a table of single SNP positions, such as a GWAS result, is refused.</p>
{table(["Your file has", "It is read as", "Example columns"], [
    ["A chromosome, a start and an end", "<b>Intervals</b>, used exactly as given.", "<code>Chromosome, Start, End</code> · <code>QTL_ID, QTL_Start, QTL_End</code> · <code>Chromosome, Left, Right</code> · <code>Chromosome, Left_Mb, Right_Mb</code> (Mb and kb names are converted to bp)"],
    ["The QTL table of QTL Mapping (<code>qtl_significant_*.csv</code>), with Left_Marker and Right_Marker", "<b>Intervals</b>. The table's positions are in cM, so the tool needs the base-pair position of its markers: from the <b>genetic map</b> (below) or, if you have no map, from a <b>VCF</b>. The confidence interval is placed between the QTL's left, peak and right markers.", "<code>QTL_ID, Chromosome, Left_Marker, Right_Marker, Left_Position …</code>"]])}
{example("gene_extractor/qtl_table_example.csv", rows=2, cols=9, caption="The QTL table of QTL Mapping (positions in cM; the marker names carry no position).")}
{example("gene_extractor/marker_positions_example.csv", rows=6, cols=3, caption="Genetic map: the base-pair position of those markers (a VCF can be given instead).")}
{example("gene_extractor/qtl_intervals_example.csv", rows=3, cols=3, caption="QTL intervals.")}
{ul(["Positions are in <b>base pairs</b> (a column name that ends in <code>Mb</code> or <code>kb</code> is converted). The chromosome of an interval table without a Chromosome column is read from the QTL ID (<code>QTL_3.1</code> is chromosome 3, <code>QTL_Chr03.2</code> is <code>Chr03</code>).",
     "An Excel workbook is read sheet by sheet: a sheet named like <i>QTL</i>, <i>Interval</i> or <i>Regions</i> is used first.",
     "A header row is expected. A file with no header and three columns — ID, start, end — is also understood.",
     "A position in <b>cM</b> is never used as base pairs: a column named <code>…_cM</code> is skipped, and a table whose positions are all below 1,000 is refused with a message. The QTL Intervals file of a QTL Mapping run that found no QTL has only a header; the tool says so, and <b>Use QTL Mapping output</b> takes the newest run that did find QTL."])}
<p>You do not have to write the QTL file yourself: <b>QTL Mapping</b> writes it (<b>QTL Intervals CSV</b>, built from the support interval of every QTL), and <b>Use QTL Mapping output</b> takes the file of your last QTL Mapping run.</p>
<p><b>2. Genetic map</b> (required) — the base-pair position of the markers: a table with <b>Marker</b>, <b>Chromosome</b> and <b>bp</b> columns (the names must be the same as <code>Left_Marker</code> / <code>Right_Marker</code> of the QTL table). The <i>Physical Map</i> sheet of a Linkage Mapping workbook works too. <b>If you have no genetic map, choose <i>No map: use a VCF</i></b> and select a VCF: its <code>CHROM</code> and <code>POS</code> give the position and its <code>ID</code> column the marker name (a VCF without IDs cannot be matched). The file needs only the markers of the QTL. It is used for a QTL table in cM; when the QTL results already have positions in bp it is not used, and the run says so. Marker positions come from the array manifest or the VCF of the markers, or from aligning the marker sequences to the genome.</p>
<p><b>3. Annotation (GFF3)</b> — the genome annotation, plain or <code>.gz</code>. The genes are read from its <code>gene</code> lines.</p>
<p><b>4. Reference genome FASTA</b> (optional) — plain <code>.fa</code> / <code>.fasta</code>, not zipped.</p>
{box("warn", "Chromosome names are matched ignoring style (<code>Chr01</code>, <code>chr1</code> and <code>1</code> are the same chromosome), but the FASTA and the GFF3 must be of the same genome version as the positions in the QTL results (and in the marker position file).")}
{sample_files("gene_extractor")}
<h3>Options</h3>
{table(["Option", "What it does"], [
    ["QTL interval", "Only for a QTL table in cM. <b>Confidence interval</b> (default) places the support interval of the QTL on the genome, as the QTL Intervals file does; <b>Left to right flanking marker</b> takes the whole span between the flanking markers (wider, more genes)."],
    ["A gene is taken when it lies", "<b>Completely inside</b> the region (default): start and stop both within it. <b>Overlapping</b> also takes the genes cut by the border of the region."]])}
<h3>Step by step</h3>
{ol(["Open <b>Genes Extraction</b>. Browse to the QTL results (or click <b>Use QTL Mapping output</b> for the QTL Intervals file, <b>Use QTL table (cM)</b> for the QTL table) and to the genetic map (or choose <b>No map: use a VCF</b> and select the VCF) and to the annotation; add the FASTA if you want the sequences.", "Set the options if needed, then click <b>Extract Data</b>.",
     "Open <b>Regions</b> to see how many genes each region holds, and <b>Genes in Gene Regions</b> for the genes. There, <b>Download Gene ID List</b> saves the gene IDs, <b>Download Genes Data</b> the table, and <b>Download Gene Sequences (FASTA)</b> the sequences."])}
{figure("gene_input")}
{figure("gene_results")}
{figure("gene_results_genes")}
<h3>Reading the results</h3>
{ul(["<b>Regions:</b> every region read from the file and its number of genes. A region with 0 genes usually means a very short interval, or that the chromosome name or the genome version does not match; regions on chromosomes that are not in the annotation are skipped and a note says how many.",
     "<b>Genes in Gene Regions:</b> one row per gene and region; a gene found in two regions (overlapping QTL) appears once for each. Columns: Chromosome, Start, Stop, Strand, Gene_id, Description (the <code>description</code> or <code>Note</code> of the annotation — many GFF3 files have none), Region_ID, Region_Start, Region_End, then Name, Feature, Source, Length_bp and every other attribute of the GFF3 line (biotype …).",
     "<b>Gene ID list</b> (<code>Gene_ID_list.txt</code>): one gene ID per line, every gene once, in genome order — ready to paste into an enrichment or ortholog tool. The Excel file has the same list on its <i>Gene_ID_list</i> sheet, with the regions each gene was found in, next to the <i>Genes</i> and <i>Regions</i> sheets.",
     "<b>Sequence</b> (with a FASTA): the genomic sequence of the gene from its start to its stop, reverse-complemented for genes on the minus strand. The sequences are also in <b>Gene Sequences (FASTA)</b>; a sequence too long for an Excel cell is only in that file."])}
"""


def ch_prioritize():
    return f"""
<p class="lead">Ranks candidate genes for a trait. It reads what public databases know about each gene — a trait annotation, its function, the pathways it is in, its orthologs in other species, and the papers that name it together with the trait — turns that into a score from 0 to 100. Give it the gene list of the QTL (the Gene_ID_list.txt written by Genes Extraction). <b>It needs an internet connection</b>, unlike the other tools.</p>
{box("warn", "The score ranks genes by the evidence found in the databases. It is <b>not</b> the probability that a gene is causal, and a low score does not mean a gene is unimportant: many genes simply have little annotation. Where a database has nothing, the tool says so; it never fills a gap with an estimate.")}
<h3>What you need</h3>
{ul(["<b>Species</b> — rice (<i>Oryza sativa</i>, IRGSP-1.0), wheat (<i>Triticum aestivum</i>) or maize (<i>Zea mays</i>). Chickpea is listed but cannot be used yet: it is not in Ensembl Plants or Gramene, so there is no verified source for its genes.",
     "<b>Trait</b> — in words, for example <i>Drought tolerance</i>, <i>Salt tolerance</i>, <i>Blast resistance</i>, <i>Grain yield</i>. Related terms are added for you (drought: water deficit, osmotic stress, abscisic acid …); the terms used are listed on the Summary tab. Built-in vocabularies cover abiotic stress (drought, salt, heat, cold, submergence, oxidative, metal toxicity, alkalinity, nitrogen, phosphorus, potassium, micronutrients), biotic stress (blast, bacterial blight, sheath blight, rusts, Fusarium / head blight, mildews, viruses, nematodes, insects, weeds, general disease), yield (grain size, grain number, tillering, harvest index, grain filling, heterosis), phenology and architecture (flowering, height, lodging, leaf angle, roots, seedling vigor, senescence, photosynthesis), seeds and reproduction (dormancy and pre-harvest sprouting, shattering, male sterility, apomixis) and quality (amylose / cooking quality, glycemic index and resistant starch, protein, oil, aroma, chalkiness, pigments and vitamins, anti-nutrients, sugar, fibre). A trait outside these groups is still searched by its own words, and the Summary tab says so; try the wording the databases use.",
     "<b>The genes to rank</b> (below)."])}
<p><b>A gene ID list</b> — paste IDs (one per line, or separated by commas, spaces or tabs) or choose a TXT / CSV / TSV / Excel file; the first column is read. For the genes of a QTL, use the <code>Gene_ID_list.txt</code> that <b>Genes Extraction</b> writes. Use the IDs of the Ensembl Plants / Gramene annotation: <code>Os02g0766700</code> (rice, RAP-DB), <code>TraesCS1A02G000100</code> (wheat), <code>Zm00001eb000010</code> (maize). Gene symbols and UniProt accessions are also looked up. Rice <b>MSU / LOC_Os IDs</b> are translated to RAP-DB IDs only for the genes that Oryzabase lists; the others are reported as not matched, so export the RAP-DB IDs instead.</p>
{sample_files("gene_prioritization")}
<h3>How the score is made</h3>
<p>Each gene receives up to five evidence scores between 0 and 1, which are combined with fixed weights:</p>
{table(["Evidence", "Weight", "Where it comes from", "What counts"], [
    ["Known gene–trait association", "30 %", "Gramene (Trait Ontology annotation of the gene; QTL covering the gene counts less); Ensembl Plants phenotypes; for rice, the curated trait classes and Trait Ontology terms of Oryzabase", "An annotation naming the trait or a related term. A direct trait term scores higher than a related one."],
    ["Functional relevance to the trait", "25 %", "Gene description, GO terms and protein domains (Gramene, Ensembl Plants); for rice, the curated gene name of Oryzabase", "Terms that match the trait words, weighted by the kind of GO evidence."],
    ["Pathway relevance", "20 %", "Plant Reactome, Gramene", "The gene is a member of a pathway whose name or summary matches the trait."],
    ["Ortholog evidence", "15 %", "Ensembl Plants homology, Gramene gene tree", "An ortholog (for example in <i>Arabidopsis</i>) that has a matching phenotype, description or trait term; one-to-one orthologs count more."],
    ["Literature support", "10 %", "PubMed (NCBI E-utilities)", "Papers in which the gene and the trait are named together. A title mention counts more than an abstract mention. This shows co-mention, not proof."]])}
<p>The <b>final score</b> is S = [Σ (weight × evidence) ÷ Σ weight] × 100, summed over the categories for which <i>data could be retrieved</i>. A category that is unavailable (for example the PubMed search failed) is left out and the other weights are scaled up — it is not treated as zero evidence. A category that was retrieved but found nothing relevant scores 0. Each piece of evidence is counted once, in one category only.</p>
<p>The <b>evidence level</b> summarises the final score: <i>Direct association</i> (a gene-level annotation names the trait or a close term), <i>Multiple lines</i> (three or more categories with evidence), <i>Two lines</i>, <i>Single line</i>, <i>No trait evidence</i>.</p>
<h3>Step by step</h3>
{ol(["Open <b>Gene Prioritization</b>. Choose the species and type the trait.", "Paste the gene IDs or browse to the gene list (for example the Gene_ID_list.txt of Genes Extraction).",
     "Click <b>Prioritize Genes</b>. The progress bar shows the step; <b>Cancel</b> stops the run. The four database connections are listed and <b>Check connections</b> tests them before you start.",
     "Read the <b>Summary</b>, the <b>Ranked Genes</b> table (click a column to sort, a row to open the gene) and the <b>Gene Evidence</b> panel. Download the ranked list (CSV), the workbook (Excel) or the report (HTML) from the Summary tab."])}
{figure("geneprio_input")}
{figure("geneprio_summary")}
{figure("geneprio_ranking")}
{figure("geneprio_evidence")}
<h3>Reading the results</h3>
{ul(["<b>Ranked Genes:</b> Rank, Gene ID, Symbol, Description, Final Score and Evidence Level. Genes are ranked by the final score. A dash means that no evidence of that kind could be retrieved; it is not a score of zero.",
     "<b>Gene Evidence:</b> the description, the position, the score of every category, each piece of evidence with its source and a link (which opens in your browser), the GO terms, domains, pathways, orthologs and the papers found. “Co-mention” means the gene and the trait are named in the same paper — read the paper before drawing a conclusion.",
     "<b>Charts:</b> the score of the top genes, what the final score of the best genes is made of (one colour per evidence category), and which of the best-ranked genes share a pathway.",
     "<b>Excel workbook:</b> sheets <i>Ranked genes</i>, <i>Evidence details</i> (one row per piece of evidence with its source, link and the date it was retrieved), <i>Publications</i>, <i>Pathways</i>, <i>GO terms</i> and <i>Run information</i> (settings, weights, service status, warnings, the identifiers that were not matched)."])}
<h3>Good to know</h3>
{ul(["<b>Time.</b> Several seconds per gene the first time (a run with about 150 candidate genes took about ten minutes in testing): every gene needs several requests to four services, which the tool spaces out to respect their limits. Start with a short list to try it out; the progress bar shows how far it is and <b>Cancel</b> stops it.",
     "<b>Cache.</b> Answers are kept on your computer for 30 days (missing answers for 7), so repeating a run is fast and works for the genes already seen. Switch on <i>Refresh cached data</i> to ask the databases again.",
     "<b>PubMed speed.</b> NCBI allows three requests per second without a key and ten with one. To use a key, set the environment variables <code>NCBI_API_KEY</code> (and <code>NCBI_EMAIL</code>) before starting TerraTools.",
     "<b>When a service is down.</b> The run continues with the services that answer; a service that fails eight times in a row is skipped for the rest of the run, so the run never waits on it for hours. The evidence it could not get is shown as unavailable and listed in the Run information. If none can be reached the tool stops and says so; nothing is invented.",
     "<b>Gaps in the databases.</b> Many plant genes have no description or GO annotation in Ensembl Plants / Gramene, and phenotype records are rare. Such a gene gets a low or empty score because the databases know little about it, not because it was judged unimportant.",
     "<b>Rice extras.</b> For rice the curated gene name, trait classes and Trait Ontology terms of Oryzabase are used (its gene list is downloaded the first time and kept in the cache; nothing is added to the program), and each gene has a link to RAP-DB."])}
"""


def ch_pybsa():
    return f"""
<p class="lead">Bulk segregant analysis (BSA-seq) from variant tables. For every variant it compares the allele depths of two bulks with Fisher's exact test, smooths the G-statistic and the
Δ(allele frequency) in sliding windows, and reports regions that exceed thresholds estimated by simulation.</p>
<h3>What you need</h3>
<p>Two tab- or comma-separated tables produced by GATK4 (<code>VariantsToTable</code>), one for the two parents and one for the two bulks. Each has the columns
<code>CHROM, POS, REF, ALT, QUAL</code> followed, for every sample, by <code>&lt;sample&gt;.GT</code>, <code>&lt;sample&gt;.AD</code> and <code>&lt;sample&gt;.GQ</code>. The parents file lists the two parents (first parent first),
the bulks file the two bulks; both must describe the <i>same</i> variants.</p>
{example("pybsaseq/bulks_example.tsv", rows=4, cols=8, caption="Bulks file (excerpt): sample names (here HighBulk and LowBulk) are taken from the column headers.")}
<p>A GATK command that writes such a table from a joint VCF:</p>
<pre class="code">gatk VariantsToTable -V bulks.vcf.gz -F CHROM -F POS -F REF -F ALT -F QUAL -GF GT -GF AD -GF GQ -O bulks.tsv</pre>
{box("warn", "<code>QUAL</code> and the allele-depth field (<code>.AD</code>) are required — without them the tool stops. The file extension must be <code>.tsv</code> (tab-separated) or <code>.csv</code>.")}
{box("key", "<b>Only variants where the two parents are completely homozygous for different alleles are used</b> (one parent has reads for the reference allele only, the other for the alternate allele only), and the same variants must be in the bulks file. Parents with a few contaminating reads at a site drop that variant.")}
{sample_files("pybsaseq")}
<h3>Step by step</h3>
{ol(["Open <b>PYBSA G-statistics</b>. Load the parents file and the bulks file.",
     "Enter the number of individuals in each bulk (<i>Bulk sizes</i>) and the population structure.",
     "Leave the other settings at their defaults for a first run, then click <b>Run Analysis</b>. Large genomes and 10 000 replications can take a long time; the Results tab opens when it is finished."])}
{figure("pybsa_input")}
{table(["Setting", "What it does"], [
    ["Bulk sizes", "Number of individuals in the first and in the second bulk, as <code>first,second</code>."],
    ["Population structure", "F2, RIL or BC; used to simulate the expected allele-frequency spread."],
    ["P-value cutoffs", "<code>real_data,simulation</code> significance levels for the observed data and for the simulated thresholds."],
    ["Number of replications", "Simulation replicates for the thresholds (more is slower but smoother)."],
    ["Sliding window", "<code>size,step</code> in bp, e.g. <code>2000000,10000</code> = 2 Mb windows moved in 10 kb steps."],
    ["Smoothing", "<code>window_len,polyorder</code> of the Savitzky–Golay curve drawn through the windows."],
    ["Region of interest", "<code>-1</code> = whole genome, or repeating triplets <code>chrom,start,end</code>."],
    ["Misc parameters", "GQ cut-off, minimum variants per window and the mirror index of Δ(allele frequency)."]])}
<h3>Reading the results</h3>
{figure("pybsa_results")}
{ul(["<b>HTML report:</b> an interactive summary with the plots and the tables; download it to share.",
     "<b>Plots:</b> number of variants, G-statistic and Δ(allele frequency) along each chromosome with their simulated thresholds. Regions where the smoothed curve crosses the threshold are the candidate QTL.",
     "<b>All Output Files:</b> every CSV the analysis wrote, for instance the list of variants, the sliding-window data and the candidate regions, each downloadable.",
     "<b>Log:</b> the program's own messages, useful when a file was rejected."])}
"""


def ch_qtlseq():
    return f"""
<p class="lead">Bulk segregant analysis from sequencing data. Two bulks of plants, one with a high and one with a low trait value, are sequenced. For every SNP the <b>SNP-index</b> of a bulk is the share of its reads that carry the <i>index allele</i>; the <b>Δ SNP-index</b> is the high bulk minus the low bulk, so it lies between −1 and 1. A QTL shows as a stretch of the genome where the sliding mean of Δ rises above (or falls below) what chance alone would give.</p>
<h3>What you need</h3>
{ul(["<b>BAM files</b> (coordinate-sorted, with an index) of the high and the low bulk, and the <b>reference genome</b> (FASTA) they were aligned to. Parents are optional: with two parents TerraTools keeps only the SNPs at which they are confidently homozygous for different alleles; with one parent that parent must be homozygous; with none, the two bulks are compared on their own and the page warns that sites which are not polymorphic between the parents can remain.",
     "<b>Or VCF files</b> of the high bulk and the low bulk (a parental VCF is optional) with an <code>AD</code> (allele depth) field. The genotype (<code>GT</code>) of a pooled bulk is never used as if it were a diploid genotype: only the read counts count. A single VCF that holds both bulks can be entered twice, with the two sample columns named under <i>Advanced settings</i>.",
     "Optionally a <b>gene annotation</b> (GFF3 or GTF) to list the genes that overlap each QTL."])}
<p>Also give the <b>population type</b> (F2, RIL or backcross) and the <b>number of plants in each bulk</b>: they set the spread of the Δ SNP-index expected by chance, which the confidence limits are built from.</p>
{example("qtl_seq/qtlseq_example.vcf.gz", rows=4, cols=11, caption="VCF example: the sample columns hold <code>GT:AD</code> values. The <code>AD</code> part (reference reads, alternate reads) is what counts. In this file the <code>HighBulk</code> and <code>LowBulk</code> columns are used.")}
{box("note", "SNP-index = (reads of the index allele) ÷ (all reads). The index allele is the allele of the parent chosen under <i>SNP-index counts the allele of</i> (parent 2 unless changed). A positive Δ means that allele is enriched in the high bulk. With no parent the index allele is the VCF ALT allele. A SNP is kept when both bulks have at least 7 reads and at least one bulk has a SNP-index of 0.3 or more; the other bulk's value is kept as it is.")}
{sample_files("qtl_seq")}
<h3>Step by step</h3>
{ol(["Open <b>QTL-Seq</b>. On the <b>Input</b> tab choose <i>BAM files</i> or <i>VCF files</i>, select the files, and fill in the population type and the number of plants in each bulk.",
     "Click <b>Check inputs</b> to see a report of the files (missing indexes, different chromosome names, missing depth fields) before the run.",
     "On <b>Settings</b> change the filters, the window, or the simulation if needed; the defaults are a good start. Click <b>Run QTL-seq</b>; the progress bar and the log follow the run, and <b>Cancel</b> stops it.",
     "Read the <b>Results</b> tab, open a chromosome by clicking a QTL, and download the tables and figures."])}
{figure("qtlseq_reads")}
{figure("qtlseq_vcf")}
{table(["Setting", "What it does", "Default"], [
    ["Minimum read depth", "SNPs covered by fewer reads in either bulk are removed.", "7"],
    ["Minimum SNP-index", "At least one bulk must reach this value; removes sequencing errors and sites that do not segregate.", "0.3"],
    ["Window size / step", "The Δ SNP-index is averaged over windows that slide along each chromosome.", "2 Mb / 10 kb"],
    ["Minimum SNPs per window", "Windows with fewer SNPs are left out.", "10"],
    ["Simulation replicates", "Random bulks drawn for every read-depth pair to form the confidence limits.", "10,000"],
    ["Random seed", "The same seed and settings always give the same limits.", "12345"],
    ["Limits of a window", "Average of the limits of its SNPs, or a simulation of the whole window.", "Average of its SNPs"],
    ["Call a QTL at", "The lowest limit (90, 95 or 99 %) a window must pass; adjacent significant windows are merged.", "90 %"]])}
<h3>Reading the results</h3>
{figure("qtlseq_results")}
{ul(["<b>Summary</b>: the informative SNPs, the number of QTLs and how many reach the 99 %, 95 % and 90 % limits (a QTL that reaches 99 % also counts at 95 % and 90 %), and a table of how many SNPs each filter removed.",
     "<b>Genome-wide plot</b>: the Δ SNP-index of every SNP, the sliding-window mean, the 90 / 95 / 99 % limits (positive and negative) and the QTL intervals. The tabs next to it show the SNP-index of each bulk. A Δ of ±1 is not a threshold; only the simulated limits decide.",
     "<b>Significant QTLs</b>: chromosome, start, end, peak position, peak / mean / minimum / maximum Δ, SNP count, the level reached, its limit and the direction (which bulk the index allele is enriched in). A click on a row opens the chromosome with the genes that overlap the interval. They are candidates by position only.",
     "<b>Simulation</b>: the distribution of the Δ SNP-index under no QTL, the limits by read depth, and the assumptions used.",
     "<b>Validation</b>: runs the pipeline on simulated data for two published designs and compares the QTLs it finds with the reported ones. The reported values are expectations, not results the program returns; the data are simulated, so a pass shows that a QTL of that size and place is recovered, not that a published analysis is reproduced.",
     "<b>Downloads</b>: the QTL table, all windows, SNP-index and Δ per position, the filtered SNPs (VCF), the figures (PNG and PDF) and a PDF report. Every run writes a folder <code>QTLseq_results</code> that also holds the parameters (JSON) and the log."])}
{box("tip", "Reads are counted with <i>samtools</i> when it is found (the fastest), otherwise with <i>pysam</i>, otherwise with a slower built-in reader that suits small genomes only. On Windows both programs are used from WSL. When <i>bcftools</i> is found, the filtered VCF is also indexed (<code>.tbi</code>) and its records and Ts/Tv ratio are checked; the page shows which programs were found. The window limits are not corrected for the number of windows tested, so a few additional short intervals can appear by chance; the 99 % level and the SNP count of an interval help to judge them.")}
"""


def ch_diversity():
    return f"""
<p class="lead">Measures haplotype diversity of each variety and sample, groups them, and draws the relationships: block statistics, a metrics heat map, comparisons between varieties,
PCA with clusters and a neighbour-joining tree.</p>
<h3>What you need</h3>
<p>A CSV with one row per sample: <b>column 1</b> the variety (cultivar group), <b>column 2</b> the sample ID, and from <b>column 3</b> the haplotype block codes — one numeric code per block per sample.</p>
{example("diversity/haplotype_blocks_example.csv", rows=5, cols=8)}
{ul(["Block values are whole numbers that name a haplotype; <code>0</code> is a real haplotype, not a missing value.",
     "Write missing blocks as <code>-1</code> or leave them empty; they are left out of each sample's metrics.",
     "Every variety needs at least a few samples for the per-variety comparisons to be meaningful."])}
{sample_files("diversity")}
<h3>Step by step</h3>
{ol(["Open <b>Diversity Analysis</b>, browse to the input matrix (CSV) and set the number of clusters, the number of top variable blocks to display and the radar-chart scale.", "Click <b>Run Full Analysis</b>."])}
{figure("diversity_input")}
<h3>Reading the results</h3>
{figure("diversity_results")}
{table(["Metric (per sample)", "Meaning"], [
    ["Richness", "Number of blocks with a call (higher = more allelic variants)."],
    ["Shannon / Allelic entropy", "Evenness of the haplotype codes: 0 = no diversity; higher = more, and more evenly spread."],
    ["Simpson / Gini-Simpson", "Probability that two random block calls differ (0 – 1)."],
    ["Heterozygosity, Gene diversity, PIC", "Expected heterozygosity and polymorphism information content across the blocks."],
    ["Effective alleles", "The number of equally frequent alleles that would give the same diversity."],
    ["Nucleotide diversity, Load", "Pairwise divergence between the sample's block calls, and the cumulative haplotype burden."]])}
{figure("diversity_plots")}
{figure("diversity_tree")}
{figure("diversity_popup")}
{ul(["<b>Block Analysis</b> ranks blocks by how much they vary; <b>Metrics Heatmap</b> shows all metrics for all samples; <b>Variety Comparison</b> has violin plots of Shannon diversity per variety and a radar chart; <b>PCA &amp; Clustering</b> show relationships; the <b>NJ Tree</b> is a circular neighbour-joining tree in the style of an iTOL export: the root is at the centre, pale coloured sectors behind the sample names mark the varieties, concentric circles and a scale bar give the genetic distance, and red triangles mark branches found in at least 70% of 100 bootstrap trees (blocks resampled with replacement; 50 trees for more than 150 samples, none above 400). Click it to enlarge it, or download it as one picture for a paper.",
     "<b>Results Table (CSV)</b> holds every metric for every sample — download it for your own plots."])}
"""
