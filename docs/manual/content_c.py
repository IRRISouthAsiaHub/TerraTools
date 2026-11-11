"""Manual content, part C: Genotype Viewer and Favourable Allele Analysis."""
from mdkit import box, example, figure, ol, sample_files, table, ul


def ch_genotype():
    return f"""
<p class="lead">The page is called <b>Graphical Genotype Visualisation</b> and is <b>used for parent recovery</b>: every chromosome of a line is drawn as a bar, coloured along its length by the genotype of the markers,
so that you can see which parent each stretch of chromosome came from and how much of the recurrent parent a line has recovered. It is built on GenoSee (Hashimoto 2024), with a choice of colour palettes,
parent names in the legend, rounded and shaded chromosomes, and PNG or vector-PDF images that are ready for a paper.</p>
<h3>What you need</h3>
<p>One genotype file, as <b>CSV</b> (or tab- / semicolon-separated text) or as a <b>VCF</b> (<code>.vcf</code> or <code>.vcf.gz</code>), with one row per marker:</p>
{example("genotype_viewer/rice_RILs_simple_example.csv", rows=5, cols=8, caption="Genotype file: chromosome, marker name, position (bp), then one column per sample (here two parents and eight RILs).")}
{ul(["<b>The columns are taken by position, not by name:</b> the first column is the chromosome, the second the marker name, the third the position in bp, and every other column is a sample. The header names do not matter.",
     "<b>Chromosome names</b> are matched to the species you choose (<code>1</code>, <code>chr1</code> and <code>chr01</code> all mean the same chromosome). Markers on chromosomes that the species does not have are not drawn, and the Results tab says how many were left out.",
     "<b>Genotype calls</b> can be written in one of three ways (use one way in the whole file):"])}
{table(["Style", "Calls", "Missing"], [
    ["Simple", "<code>A</code> (parent 1), <code>B</code> (parent 2), <code>H</code> (heterozygous)", "<code>N</code>"],
    ["Phased", "<code>0|0</code>, <code>0|1</code>, <code>1|0</code>, <code>1|1</code>", "<code>.|.</code>"],
    ["Non-phased", "<code>0/0</code>, <code>0/1</code>, <code>1/1</code>", "<code>./.</code>"]])}
{box("warn", "Multi-allelic sites (for example <code>0/2</code>) and polyploid calls (<code>0/0/0/1</code>) cannot be drawn. Remove them first, for example with <code>bcftools view --max-alleles 2</code>.")}
<p><b>A VCF</b> is converted with GenoSee's own converter: each sample's <code>GT</code> field becomes the call. The option <i>Name the markers CHROM_POS</i> names the markers
<code>chr01_1735605</code>; untick it to use the ID column of the VCF instead.</p>
{sample_files("genotype_viewer")}
<h3>Species</h3>
<p>The chromosome lengths come from GenoSee's table of species: rice (<i>Oryza sativa</i>), <i>Arabidopsis thaliana</i>, wheat, sorghum, tomato, maize, barley and soybean.
For another species or genome version, use <b>Add a species</b>: choose its genome FASTA file (<code>.fa</code>, <code>.fasta</code>, also <code>.gz</code>) and type a name; the chromosome lengths are read from the file and kept for the next time.</p>
<h3>Step by step</h3>
{ol(["Open <b>Genotype Viewer</b> (Trait Mapping, after Linkage Mapping) and click <b>Browse</b> to choose the genotype file. The number of samples and their first names appear under the box.",
     "Choose the <b>species</b>. To draw only some samples, type their names in <b>Samples to draw</b>, separated by commas.",
     "Choose the <b>drawing mode</b>, the <b>colouring</b>, the <b>colour palette</b> and the <b>image quality</b> (see the table below).",
     "Optional: type what REF and ALT stand for in <b>Legend names</b>, for example <i>Parent 1</i> and <i>Parent 2</i>.",
     "Click <b>Draw Graphical Genotypes</b>. The Results tab opens with the images."])}
{figure("genotype_input")}
{table(["Setting", "What it does"], [
    ["Drawing mode — Normal", "One image per sample, with all chromosomes side by side. A file of 200 samples gives 200 images: name the samples you need."],
    ["Drawing mode — Compare", "One image per chromosome, with all samples side by side, to compare them along that chromosome."],
    ["Drawing mode — Zoomed", "One region of one chromosome for all samples. Give the chromosome (for example <code>chr01</code>) and the start and end in bp; the end cannot be beyond the chromosome."],
    ["Colouring", "<b>2-color</b>: REF and ALT (heterozygous calls are drawn in the colour of the nearer allele call). <b>3-color</b>: REF, HET and ALT."],
    ["Graphic style", "On (default): rounded, softly shaded chromosomes with a thin outline and a light shadow, no black lines inside them, a pale grid and bold labels. Off: the plain GenoSee drawing."],
    ["Fill the area between two markers", "On: the stretch between two markers is coloured (the usual graphical genotype). Off: only the markers are drawn."],
    ["Show marker names in the figure", "Off by default (chromosomes only, cleaner for a paper). On: the marker names are written beside the chromosomes."],
    ["Dashed line from each marker name", "Only with marker names in the Compare and Zoomed modes: the dashed line that runs from each name across all samples. Switch it off to keep the names only."],
    ["Colour palette", "See the next table."],
    ["Legend names", "The words that replace REF, ALT and HET in the legend. Empty boxes keep GenoSee's words."],
    ["PNG resolution", "200 dpi (screen), 300 dpi (print) or 600 dpi (publication)."],
    ["Also save a vector PDF", "Writes a vector PDF next to every PNG; each result has a download button for it."]])}
<h3>Colour palettes</h3>
{table(["Palette", "Colours (REF / HET / ALT)", "Use"], [
    ["normal — bright", "yellow / green / blue (missing: grey)", "GenoSee's default, bright and easy to tell apart."],
    ["Aqua", "pale green / blue / dark blue", "A calm blue-green series."],
    ["grays", "black / grey / white", "For black-and-white printing."],
    ["reds", "red / light red / pale", "One-colour series."],
    ["blues", "blue / light blue / pale", "One-colour series."],
    ["Red-Blue-Yellow-Grey", "red / yellow / blue (missing: grey)", "A mixed palette in which the two parents and the heterozygotes stand out clearly."]])}
<p>The first five palettes are GenoSee's own; the sixth is added by TerraTools in the same format.</p>
<h3>Reading the results</h3>
{figure("genotype_results")}
{ul(["<b>The bars:</b> each bar is a chromosome, drawn to scale in Mb. Each stretch is coloured by the genotype of the marker: with the default coding yellow is REF, blue is ALT, green is heterozygous and grey is missing. The colour changes where the genotype changes between two markers.",
     "<b>The legend</b> (at the right in Normal and Compare mode, below the samples in Zoomed mode) names the colours. Type your own names in <i>Legend names</i> to show the parents.",
     "<b>Download:</b> every image has a <b>Download PNG</b> button, <b>Download all as PNG</b> saves them together, and with <i>Also save a vector PDF</i> each image also has a PDF."])}
{box("key", "<b>Which colour is which parent?</b> The viewer only knows allele types: <code>A</code> / <code>0|0</code> is REF (the first colour) and <code>B</code> / <code>1|1</code> is ALT (the second colour). Code the file so that REF is parent 1 and ALT is parent 2, and put the parents into the legend with <i>Legend names</i>. Adding the parents as sample columns is a good check: their bars are one colour all along.")}
{box("tip", "For a figure in a paper: leave the marker names off, choose 600 dpi and tick <i>Also save a vector PDF</i>; use the PDF for the final layout.")}
{box("note", "GenoSee is described in Hashimoto S. (2024) <i>GenoSee: a novel visualization tool for graphical genotypes</i>, Breeding Science 74(5): 454–461 (MIT licence). Please cite it when you use these figures. GGT 2.0 is a closed program without source code, so it is not part of TerraTools.")}
"""


def ch_favallele():
    return f"""
<p class="lead">Shows which <b>favourable and unfavourable alleles</b> of known trait genes each of your lines carries. You upload the genotypes of your lines (a <b>rice 4K haplopanel</b> file);
TerraTools matches the marker IDs of your file with its own list of trait-gene markers — 327 SNPs and indels in <b>166 known trait genes</b> for traits such as grain quality, yield, flowering,
biotic and abiotic stress and nutrient use — and shows the alleles of your lines at the markers it finds. It is the genetic signature of a line: the list of trait alleles it has, that a breeder can use to choose parents and lines.</p>
{box("key", "The tool is for the <b>rice 4K haplopanel only, for now</b>. The 166 trait genes are built into TerraTools and are not shown or exported (only the markers that match your file are shown). It has <b>not been scaled yet</b> to other marker panels, other crops or other genes; markers of your file that are not among the trait-gene markers are ignored.", title="Scope")}
<h3>What you need</h3>
<p>One genotype file of your lines: a CSV, TXT or tab-separated file (or a VCF, also gzipped) with <b>one row per marker and one column per line</b>; the first row holds the line names.</p>
{table(["Part of the file", "What it must look like"], [
    ["Marker ID", "<code>chr1_7268879</code> (chromosome_position). The ID may be followed by the gene, <code>chr1_7268879_LOC_Os01g13060</code>. If the ID does not match, columns named <code>chrom</code> and <code>pos</code> are used (matched with the position of the trait-gene marker). The ID column need not be the first one; TerraTools finds it."],
    ["Genotypes", "<code>0</code> = REF/REF, <code>1</code> = heterozygous, <code>2</code> = ALT/ALT, <code>-1</code> (or empty, NA) = missing. VCF genotypes (<code>0/0</code>, <code>0/1</code>, <code>1/1</code>, <code>./.</code>) and nucleotides (<code>A/G</code>, <code>AG</code>) are converted to these codes."],
    ["REF and ALT (optional)", "Columns <code>REF</code> and <code>ALT</code>, or a VCF, are compared with the alleles of the trait-gene marker: if they are swapped the genotypes are flipped, and a marker whose alleles do not fit at all is left out and reported."]])}
<h3>Step by step</h3>
{ol(["Open <b>Favourable Allele Analysis</b> (Genomic Selection). On <b>Your 4K File</b> click <b>Browse</b> and choose your genotype file. It is matched at once.",
     "Read the match report: the markers in your file, how many of the trait-gene markers were found (<b>Matched IDs</b>), the trait genes covered and the number of lines. The <b>Matched marker IDs</b> table lists every matched ID with its gene, trait and how it was matched; <b>Download matched IDs</b> saves it as CSV.",
     "Optional: restrict the genes with <b>Traits</b> and <b>Genes</b> (the list follows the traits; type to search) and choose a <b>View status</b> to look at one kind of call only. Nothing selected means all.",
     "In <b>Select lines to compare</b> search and tick one or more lines. The selected ones appear as coloured names: click a name to make it the line shown in detail, or <b>×</b> to remove it.",
     "Click <b>Show results</b>. The <b>Overview</b>, <b>Genes</b> and <b>SNP Table</b> tabs show the first selected line in detail; click another name in the bar at the top of any of them to switch. The <b>Heatmap</b> tab draws the balloon heatmap (see below)."])}
{figure("favallele_upload_empty")}
{figure("favallele_setup")}
<h3>How a SNP is scored</h3>
<p>Every trait-gene SNP has a reference allele (REF), an alternative allele (ALT) and a <b>favourable</b> and an <b>unfavourable</b> allele from the literature. For each line the call at the SNP gives one of four statuses:</p>
{table(["Status", "Meaning"], [
    ["Favourable", "The line is homozygous for the favourable allele."],
    ["Unfavourable", "The line is homozygous for the other allele."],
    ["Heterozygous", "The line carries both alleles."],
    ["Missing", "No call for the line, or the favourable allele matches neither REF nor ALT."]])}
<p>The favourable allele is compared with REF and ALT; when it matches only the reverse complement of one of them, the annotation is taken to be on the opposite strand and is flipped. A deletion, written <code>-</code>, counts as an allele like any other.</p>
{box("note", "For <b>A/T and C/G SNPs</b> the strand cannot be told from the alleles, so the favourable allele is used as written.")}
<h3>Reading the results</h3>
{figure("favallele_overview")}
{ul(["<b>The bar at the top</b> (on every result tab) lists the selected lines and the filters in force; the highlighted name is the one shown in detail.",
     "<b>Summary cards:</b> <i>Fav SNPs</i> and <i>Unfav SNPs</i> (homozygous favourable / unfavourable), <i>Het</i>, <i>Missing</i>, and <i>Fav genes</i> / <i>Unfav genes</i> — the genes with at least one favourable / unfavourable SNP. A gene with both kinds of SNP is counted in both.",
     "<b>Line profile:</b> every gene is a dot on the rice plant, placed on the plant part its trait belongs to (panicle, upper stem, leaves, lower stem, stem or roots). Hover for the gene, click for all its SNPs. Colour the dots by trait or by allele status (favourable, unfavourable, mixed, heterozygous, missing); click a trait in the legend to highlight it. <b>Download PNG</b> saves the picture.",
     "<b>SNPs by trait:</b> the number of SNPs of the chosen status for every trait, one bar per selected line."])}
{figure("favallele_gene")}
<p>The <b>gene window</b> gives the trait, sub-category and position of the gene and the call of every selected line at each SNP, written as the two nucleotides (for example <code>TG</code>); the colour shows the status.</p>
{figure("favallele_genes")}
{ul(["<b>Genes tab, allele status by gene:</b> a bar for each gene of the line shown, sorted by the number of favourable SNPs; click a bar for the gene window.",
     "<b>Allele map:</b> lines as rows and the SNPs of the best-typed genes as columns, coloured by status, to compare the signatures side by side. Click a cell or a gene name for the gene window.",
     "<b>SNP Table:</b> every SNP of the line shown with REF, ALT, favourable allele, stored code, call and status; filter the rows with the box above the table and <b>Download CSV</b> to keep it."])}
<h3>Balloon heatmap</h3>
<p>The <b>Heatmap</b> tab compares many genes and many lines in one figure. It is drawn by R (the script of the <i>Heatmap_Final</i> folder), on the genes and lines you choose from your file.</p>
{ol(["Open the <b>Heatmap</b> tab. In <b>Genes</b> tick the genes you want (type to search). <b>Add all the genes of a trait</b> adds a whole trait with one click, and <b>Genes of the filters</b> takes the genes chosen on the first tab. The order in which you pick the genes is the order on the figure.",
     "In <b>Lines</b> tick the lines (type to search), or take <b>All lines</b> or the lines selected on the first tab.",
     "Choose the <b>resolution</b> (150, 300 or 600 dpi) and click <b>Draw Balloon Heatmap</b>. The figure appears below the settings; click it to enlarge it and use the download button to save the PNG. The gene frequency matrix, the marker QC report and the favourable counts per line can be saved as CSV files."])}
{figure("favallele_heatmap_setup")}
{figure("favallele_heatmap", caption="Favourable Allele Analysis — balloon heatmap of 12 genes in 24 lines (drawn by the R script)", width="100%")}
{ul(["<b>Balloons:</b> one per line and gene. Size and blue depth show the <i>favourable allele fraction</i>: the share of the gene's markers at which the line is homozygous for the favourable allele (a heterozygous call counts as not favourable). A grey cross means no call at any marker of the gene.",
     "<b>Coloured strip on top:</b> the trait category of each gene. <b>Bars on the right:</b> the number of favourable alleles of each line, out of the markers kept.",
     "<b>Quality control:</b> before drawing, markers whose favourable allele matches neither REF nor ALT, and markers missing in 80 % or more of the <i>chosen</i> lines, are removed. The page lists them, so choosing very few lines can remove markers that are well typed in the whole file."])}
{box("note", "17 further trait genes have no scorable marker (every marker is an indel whose favourable allele cannot be matched to REF or ALT) and are not part of the analysis; the <b>About</b> tab lists them. The figure grows with the number of genes and lines; for a poster choose 600 dpi.")}
{box("tip", "Select two or three lines to compare their signatures — for example a donor of a favourable allele and the line you want to improve — and use the allele map to see at a glance where they differ.")}
"""


def ch_genofile():
    return f"""
<p class="lead">The <b>start of the chain</b>. Genotyping companies deliver alleles — <code>A/G</code>, <code>AG</code>, <code>0/1</code> — in HapMap, VCF or spreadsheet files, but the analyses of TerraTools read genotypes <b>coded against the two parents</b>
(0, 1, 2). The Genofile Converter reads the company file, codes every individual against the parents you name, and writes <b>one coded genotype table</b> — the starting point for the other analyses.</p>
<h3>What you need</h3>
<p>One genotype file from the company, in one of three formats; the format is detected from the file.</p>
{table(["Format", "What it looks like", "Calls"], [
    ["HapMap (<code>.hmp.txt</code>)", "Tab-separated, one row per marker: <code>rs#</code>, <code>alleles</code>, <code>chrom</code>, <code>pos</code>, <code>strand</code> … <code>QCcode</code>, then one column per sample.", "<code>AA</code>, <code>AG</code>, <code>GG</code>; a single letter is a homozygote (<code>A</code>) or an IUPAC code (<code>R</code> = A/G); <code>NN</code> is missing."],
    ["VCF (<code>.vcf</code>, <code>.vcf.gz</code>)", "<code>#CHROM POS ID REF ALT … FORMAT</code>, then one column per sample. The GT field is read.", "<code>0/0</code>, <code>0/1</code>, <code>1/1</code> (or with <code>|</code>); <code>./.</code> is missing. More than two alleles are handled."],
    ["Raw allele calls (Excel, CSV)", "One row per marker: <code>Marker</code>, <code>Chromosome</code>, the parents <code>P1</code> and <code>P2</code>, then one column per individual.", "<code>A:G</code>, <code>A/G</code>, <code>A|G</code> or <code>AG</code>; <code>-</code> or <code>N</code> is missing."]])}
{example("genofile/genofile_raw_alleles_example.xlsx", rows=5, cols=8, caption="Raw allele sheet (Excel or CSV): Marker, Chromosome, the parents P1 and P2, then one column per individual.")}
{sample_files("genofile")}
{box("key", "<b>The parents define the codes.</b> A marker is coded only when both parents are called, each is homozygous, and they carry different alleles. An individual is then <code>0</code> when it has Parent 1's allele twice, <code>2</code> for Parent 2's allele twice, <code>1</code> with one of each, and <code>-1</code> when it is missing or carries an allele neither parent has.")}
<h3>Step by step</h3>
{ol(["Open <b>Genofile Converter</b> (Data Analysis) and click <b>Browse</b>. Choose the file; the page shows the format it found (change it with the three buttons if it is wrong) and a preview.",
     "<b>HapMap and VCF:</b> in <b>Parents and samples</b> choose Parent 1 and Parent 2 from the samples of the file (names such as <code>P1</code> and <code>P2</code> are found automatically). Samples that are not part of the population — controls, other lines — can be listed under <b>Samples to leave out</b>. "
     "<b>Raw sheet:</b> check the columns found for Marker, Chromosome, P1 and P2 and the allele separator; every other column is an individual.",
     "HapMap and VCF: set the <b>Output settings</b> — the marker names (as in the file, or <code>chromosome_position</code>, for example <code>chr01_1044946</code>, which carries the physical position of the marker, read by analyses such as Linkage Mapping), the chromosome labels (<code>Chr01</code> becomes <code>1</code>) and whether markers on scaffolds are left out.",
     "Click <b>Convert Genotype File</b>. The Results tab opens.",
     "Download the coded file (CSV, or Excel with the reports), or continue with it in another tool with the buttons under <b>Your converted file</b>."])}
{figure("genofile_input")}
<h3>Reading the results</h3>
{figure("genofile_result")}
{ul(["<b>Summary cards:</b> markers in the file, markers coded (kept), left out, individuals and the share of calls that are missing or unexpected.",
     "<b>Your converted file:</b> one table with a row per marker and the columns <i>Marker</i>, <i>Chromosome</i> and one column per individual, filled with <code>0</code> (Parent 1 allele), <code>1</code> (heterozygous), <code>2</code> (Parent 2 allele) or <code>-1</code> (missing / unexpected). The parents are not in the table. It is offered as a <b>CSV</b> (the table) and as an <b>Excel workbook</b> (the same table on the first sheet, then <i>Marker_Details</i> — the missing and unexpected calls of each marker, <i>Excluded_Markers</i> — every marker left out and why, and <i>Summary</i> — the settings and counts).",
     "<b>Why markers were left out:</b> a parent call missing, a parent heterozygous, both parents with the same allele, or (when chosen) an unplaced contig. <b>Genotype codes</b> shows how many calls were 0, 1, 2 and −1.",
     "<b>Coded markers</b> and <b>Markers left out</b> list the first rows, with the reason for each marker that was left out."])}
{box("note", "A very large file (a table of more than 1.5 million coded cells) is not written as an Excel sheet, which would take minutes to write and open; the table is then given as the CSV only, and a small workbook holds the reports. HapMap and VCF files are read line by line, so even a large file is converted without filling the memory.")}
{box("tip", "Convert once, then use the files in every analysis of the same population. If you used <i>chromosome_position</i> marker names, Linkage Mapping can read the physical positions from them.")}
"""
