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
     "Download the coded file (CSV, or Excel with the reports), or continue with it in another tool with the buttons under <b>Your converted file</b>: <b>Mendelian Segregation</b>, <b>Linkage Mapping</b> and <b>Recombination Analysis</b> open with the converted file already loaded."])}
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


def ch_recombination():
    return f"""
<p class="lead">Recombination Analysis answers three questions about a mapping population of RILs or backcross-derived lines: <b>where did the crossovers happen</b>, <b>how much of the donor genome does each line carry</b>, and <b>where along the chromosomes is recombination frequent or rare</b>.
It replaces the set of R scripts that were run one after the other (crossover detection, rate along the genome, circos plot, intragenic crossovers) with one tool.</p>
<h3>What you need</h3>
<p>One table with <b>a row per marker</b>, giving the chromosome, the <b>position in base pairs</b> and the genotype of every line, coded <b>0 / 1 / 2</b> with <b>−1</b> for missing. Two kinds of file are read:</p>
{ul(["<b>A table with the recurrent parent in it</b> (the layout of a VCF-style table): <code>CHROM</code>, <code>POS</code>, <code>ID</code>, optional <code>REF</code>, <code>ALT</code> …, then one column per sample (<code>IRRI-370_GT</code>, <code>RIL-001_GT</code> …). The parent's allele is called <b>A</b>, the other homozygous allele <b>B</b> (the donor).",
     "<b>A coded table</b> from the Genofile Converter (<code>Marker</code>, <code>Chromosome</code>, then the individuals): there is no parent column, <code>0</code> is the allele of Parent 1 (A) and <code>2</code> the allele of Parent 2 (B). The position is read from the marker name (<code>chr01_1044946</code>). The Genofile Converter has a button that opens its result here."])}
{example("recombination/mapping_population_example.csv", rows=5, cols=9, caption="Genotype table: CHROM, POS, ID, REF, ALT, then one _GT column per sample.")}
{sample_files("recombination")}
{box("key", "If the marker names carry the gene (<code>chr1_134521_LOC_Os01g01280</code>) the crossovers are also split into <b>intragenic</b> (both flanking markers lie in the same gene) and <b>intergenic</b>, and the genes with the most intragenic crossovers are listed with their names.")}
<h3>Step by step</h3>
{ol(["Open <b>Recombination Analysis</b> (Trait Mapping) and click <b>Browse</b>. The page reads the file and shows the number of markers, chromosomes and samples, the share of heterozygous and missing calls and a preview.",
     "Under <b>Parent and lines</b> choose the <b>recurrent parent</b> (the page guesses it). Name the <b>donor</b> only if it is genotyped in the file. List samples that are not part of the population — the donor, controls, a line that failed — under <b>Samples to leave out</b>. Choose <i>none</i> for a coded table.",
     "Set the <b>smoothing</b> (below), the bin size of the rate plot and whether the rice centromeres are shaded.",
     "Click <b>Run Analysis</b>. The Results tab opens with the figures, the downloads and the tables."])}
{figure("recomb_input")}
<h3>Smoothing: the setting that matters most</h3>
<p>A genotyping error inside a block of the other parent's allele would look like two crossovers. TerraTools therefore overwrites every run of A or B calls <b>shorter than a minimum block</b>, given as a percentage of the markers of the chromosome (at least 3 markers). Three methods are offered:</p>
{table(["Method", "What it does"], [
    ["Original script", "The first run shorter than the minimum takes the value of the run before it (the first run takes the value of the next one); repeated up to 50 times. Reproduces the results of the earlier R scripts exactly."],
    ["Shortest run first", "The shortest run is merged into its longer neighbour, again and again, until every run has the minimum length. No limit on the number of passes."],
    ["None", "Every isolated call counts as two crossovers. Useful only to see how much noise there is."]])}
{box("warn", "<b>The block size changes the numbers a great deal.</b> In a simulation with 489 known crossovers in 80 lines, a block of 2 % recovered 405 of them, 3 % about 350 and 5 % only 264 to 290, depending on the method,, because real short donor segments are shorter than the block and are removed with the errors. On the lab's 4K data the same file gave 3,432 crossovers with the original script and 5,686 with the shortest-run-first method at 5 %. Choose the block from what you know about the shortest introgression you expect, report it with your results, and compare two settings before you rely on a count.")}
<h3>Reading the results</h3>
{figure("recomb_result")}
{ul(["<b>Summary cards:</b> lines, informative markers (the parent homozygous and at least one line with the other allele; the others cannot show a crossover), crossovers, mean per line, mean donor % and, when the marker names carry genes, the share of intragenic crossovers.",
     "<b>Figures:</b> crossovers per line and per chromosome; the <b>genome-wide crossover rate</b> (crossovers per line in each bin, with the centromeres shaded); donor content per line and per chromosome; the <b>haplotype mosaic</b> (one row per line, sorted by donor %; yellow = parent, red = donor, blue = heterozygous, grey = missing); a <b>circos</b> version of it; and, with genes, intragenic versus intergenic crossovers and the genes with the most intragenic crossovers.",
     "<b>Please note:</b> lines with more than half of their calls missing (they cannot be used), very heterozygous lines, and a reminder that a crossover is only located to the interval between its two flanking markers.",
     "<b>Tables and downloads:</b> the <b>Excel workbook</b> has the sheets <i>Summary, Settings, Per_chromosome, Per_RIL, Crossovers, CO_rate_by_bin</i> and <i>Genes</i>; the <b>ZIP</b> has the same tables as CSV files and every figure as a PNG. The crossover table gives, for each crossover, the line, chromosome, the two flanking markers and their positions (<i>BP_start</i>, <i>BP_end</i>), the middle of the interval, its size, the number of markers inside it, the direction (A→B or B→A) and, with genes, the type and the gene."])}
{box("note", "<b>Donor % per line</b> is the average of the chromosomes, each the share of B among its A and B calls; the column <i>Donor_pct_all_markers</i> pools all markers instead. A line with no locatable crossover is counted in the mean per line with 0.")}
{box("tip", "In backcross-derived lines (mean donor around 20 %) crossovers can only be seen inside donor segments, so the crossovers per chromosome say as much about where the donor was kept as about the recombination rate. Compare chromosomes with that in mind.")}
"""


def ch_insilico():
    return f"""
<p class="lead"><b>In silico BSA</b> is bulked segregant analysis without any pooled sequencing. You already have the genotype of every line of a mapping population and its trait value, so TerraTools sorts the lines by the trait, takes the lowest and the highest fraction as the <b>low and the high pool</b>, and compares the genotypes of the two pools at every marker.
It uses the <b>Gw statistic</b> of dQTG-seq (Li et al. 2022), which counts alleles <i>and</i> genotypes. Because it also looks at the heterozygotes, it finds <b>over-dominant</b> loci that allele-frequency methods (G′, ED, Δ SNP-index) cannot see.</p>
{box("note", "This tool is for a population in which every line is genotyped and phenotyped (F2, F3, IMF2, BC1, RIL, DH). If you have <i>sequenced pools</i> (a high bulk and a low bulk), use <b>QTL-Seq</b> or <b>PYBSA G-statistics</b> instead.")}
<h3>What you need</h3>
<p>A genotype file and a phenotype file. Both are CSV or Excel files.</p>
{ul(["<b>Genotype file</b> — one row per marker with a <code>CHROM</code> and a <code>POS</code> column (the position in base pairs), optionally an <code>ID</code> or <code>Marker</code> name, then one column per line, coded <b>0 / 1 / 2</b> with <b>−1</b> for missing: 0 and 2 are the two parents and 1 is heterozygous. The positions are read from <code>POS</code>, never from the marker name; only when there is no <code>POS</code> column are names such as <code>chr1_134521</code> used for the position. The coded files of the <b>Genofile Converter</b> work.",
     "<b>Phenotype file</b> — the line names in the first column, written exactly as in the genotype file, then one numeric column per trait. Lines without a trait value are left out. A GWAS-ready file (one row per line with a trait column) can be given as the genotype file and then needs no phenotype file.",
     "At least <b>20 lines</b> with both a genotype and a trait value."])}
{example("insilico_bsa/insilico_genotypes_example.csv", rows=4, cols=8, caption="Genotype file: CHROM, POS, ID, then one column per line (0 / 1 / 2, −1 = missing).")}
{example("insilico_bsa/insilico_phenotypes_example.csv", rows=4, cols=2, caption="Phenotype file: the line name and the trait.")}
{sample_files("insilico_bsa")}
<h3>Step by step</h3>
{ol(["Open <b>In silico BSA</b> (Trait Mapping), click <b>Browse</b> and choose the genotype file, then the phenotype file. The page reports the markers, chromosomes and lines, the share of heterozygous and missing calls, and how many line names matched.",
     "Choose the <b>trait</b> (when the file has several) and the <b>population type</b>. The page tells you what the heterozygosity of your data looks like, but the choice is yours.",
     "Set the share of lines in each <b>pool</b> (20 % by default) and check the windows and the threshold (below).",
     "Click <b>Run Analysis</b>. The Results tab opens when the run ends; with 900 markers and 1,000 permutations it takes a few seconds."])}
{figure("insilico_input")}
<h3>Population type</h3>
<p>The type decides how heterozygous calls are used, and the page checks the heterozygosity of the data against it.</p>
{table(["Population", "Heterozygous calls", "Statistic"], [
    ["<b>F2 / F3 / IMF2</b> (1 : 2 : 1)", "Used. Three genotype classes.", "Gw (alleles and genotypes)"],
    ["<b>BC1</b> (1 : 1)", "Used. Two genotype classes.", "Gw"],
    ["<b>RIL</b>, <b>DH</b>", "Treated as missing: every line carries one allele.", "Gw equals G1 (alleles only), so the tool cannot find an over-dominant locus in these populations."]])}
<h3>Settings</h3>
{table(["Setting", "Default", "What it does"], [
    ["Lines in each pool", "20 %", "The lowest and the highest fraction of the lines by trait value (5 to 50 %). A larger pool is more stable but dilutes the contrast."],
    ["Primary window", "500 kb", "The window size used for the headline threshold, the genotype figure and the comparison with G1 and ΔAF. All three sizes (250 kb, 500 kb, 1 Mb) are analysed and reported in every run."],
    ["Minimum scored markers per window", "0 = automatic", "Automatic is 3, or fewer on a sparse map (a 300-marker map has about one marker per window, so a fixed 3 would test almost nothing). Windows below the minimum are not tested."],
    ["Merge windows closer than", "1000 kb", "Significant windows of a chromosome form one region when the gap between them is at most this. Choose it from the marker density and the decay of linkage disequilibrium; the Settings sheet lists the median marker spacing."],
    ["Minimum called lines per pool", "5", "A marker with fewer calls in a pool is not scored."],
    ["Minimum minor allele frequency", "0.05", "Rare markers are left out before the pools are compared."],
    ["Permutations", "1,000", "200 to 5,000. More permutations give a smoother threshold and finer adjusted p-values."],
    ["Significance level", "0.05", "0.10 (a looser level, suggested by Li et al. for small samples), 0.05 (genome-wide) or 0.01 (stringent)."]])}
<h3>How the result is obtained</h3>
{ul(["<b>Counts.</b> At every marker the two pools give a pool × allele table (2 × 2) and a pool × genotype table (2 × 3). A marker that cannot be scored gets <i>no</i> statistic (never a zero) and does not enter any window.",
     "<b>Statistics.</b> <b>G1</b> is the G statistic of the allele table, <b>G2</b> that of the genotype table, and <b>Gw = G1/(G1+G2) · G1 + G2/(G1+G2) · G2</b>. A fully over-dominant locus has the same allele frequency in both pools (G1 near 0) but many more heterozygotes in one pool (a large G2), so only Gw sees it. The allele-frequency difference (the Δ SNP-index of QTL-seq) is reported too.",
     "<b>Windows.</b> The genome is cut into non-overlapping windows of 250 kb, 500 kb and 1 Mb. The value of a window is the mean of the marker statistic over the scored markers it holds.",
     "<b>Threshold.</b> Gw has no known distribution, so two pools of the same size are drawn at random from all lines many times and the largest window value of the genome is kept. The threshold of each window size is the (1 − α) quantile of those maxima, and the adjusted p-value of a window is (1 + the number of permutations that reach it) / (permutations + 1). This controls the error across all windows of that size."])}
<h3>Reading the results</h3>
{figure("insilico_result")}
{ul(["<b>Please note:</b> warnings about the data, for example a heterozygosity that does not fit the chosen population or lines that were left out.",
     "<b>Loci across the window sizes:</b> regions of the three sizes that overlap form one <i>locus</i>. <i>Windows_that_found_it</i> says how many sizes reach the threshold there (a locus found by all three does not depend on the window you chose), and <i>Peak_spread_kb</i> is how far apart the peak markers of the sizes lie.",
     "<b>Significant regions of every window size:</b> one row per region and window size, with the left, peak and right marker and their positions, the adjusted p-value at the peak, the LOD (Gw / ln 10) and the Δ SNP-index (the allele of parent 2: positive means it is more frequent in the high pool). The Excel report adds the genotype classes at the peak. The <i>Pattern</i> column reads them: an <b>allele-frequency shift</b> (additive or dominant), or <b>heterozygotes enriched</b> in one pool without an allele shift (over-dominant). <i>Also_significant_by_G1</i> and <i>…_by_dAF</i> say whether the allele-based statistics find the same region.",
     "<b>Closest to significant:</b> shown for any window size without a significant region: the regions nearest to the threshold, with the fraction of the threshold they reach. They are <b>not significant</b>; treat them as leads.",
     "<b>Candidate genes (a separate step):</b> choose a window size, the GFF3 or GTF annotation of the <i>same genome assembly</i> as your marker positions and a flank (100 kb by default). The tool lists the genes that overlap each region, nearest to the peak marker first, and saves them as Excel and CSV. They are positional candidates, not the causal gene.",
     "<b>Figures:</b> Gw along the genome for the three window sizes; the loci found by each size; Gw compared with G1 and the allele-frequency difference; the genotype classes of the two pools at each peak; the trait distribution with the two pools; and the permutation null distributions.",
     "<b>Downloads:</b> the Excel workbook (<i>Summary, Settings, Window comparison, Loci, Regions, Thresholds, Windows, Closest to significant, Pools, Markers</i>), a ZIP with the same tables as CSV files and all figures, and the main tables as single CSV files."])}
{box("key", "<b>The example:</b> the simulated F2 population of 300 plants has three QTL of about 8 % of the variance each. With 20 % pools, the F2 type and 1,000 permutations, all three window sizes find all three loci. The additive and the dominant QTL (chromosomes 3 and 1) are also found by G1. The over-dominant QTL on chromosome 2 (peak near 11.1 Mb) is found <b>only through the genotype counts</b>: 71 % of the high pool is heterozygous against 22 % of the low pool, with almost the same allele frequency in both.")}
{box("warn", "The permutation test treats the lines as exchangeable, which fits an F2, IMF2, backcross or RIL population from <b>one cross</b>. In a diversity panel with family structure, related lines in the same pool give false signals everywhere; use <b>GWAS</b> with a kinship matrix there. A very polygenic trait gives many regions: compare them with <b>QTL Mapping</b> before you choose candidates.")}
<p>Li P, Li G, Zhang Y-W, Zuo J-F, Liu J-Y, Zhang Y-M (2022). A combinatorial strategy to identify various types of QTLs for quantitative traits using extreme phenotype individuals in an F2 population. <i>Plant Communications</i> 3:100319. Magwene PM, Willis JH, Kelly JK (2011). The statistics of bulk segregant analysis using next generation sequencing. <i>PLoS Computational Biology</i> 7:e1002255.</p>
"""
