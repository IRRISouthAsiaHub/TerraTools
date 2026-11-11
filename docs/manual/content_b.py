"""Manual content, part B: GWAS, genomic-selection tools, phenotypic analysis, appendices."""
from mdkit import box, example, figure, ol, sample_files, table, ul


def ch_gwas():
    return f"""
<p class="lead">Genome-wide association study: tests every marker against a trait with up to eight models and reports the marker–trait associations (MTAs) that pass a stated significance criterion, together with the quality-control steps, diagnostics and limits of each model.
You can use the <b>built-in rice 3K panel</b> and supply only phenotypes, or bring your own genotypes as a VCF or as a ready-made file.</p>
<h3>Choose where the genotypes come from</h3>
{table(["Genotype source", "You provide", "Use it when"], [
    ["<b>Rice 3K panel — phenotype file only</b>", "a phenotype file with the 3K accession IDs", "your lines are rice 3K accessions (the same panel as HapGUI: 3,024 accessions and 52,600 SNPs are built in)"],
    ["<b>My own genotypes — VCF + phenotype file</b>", "a VCF (.vcf / .vcf.gz) and a phenotype file", "you have genotype calls for your own panel"],
    ["<b>My own genotypes — ready-made GWAS file</b>", "one file: ID, phenotype, SNP columns", "you already prepared the exact GWAS input"]])}
<h3>Option 1 — the rice 3K panel</h3>
<p>Prepare one phenotype file (CSV, TSV or Excel). The <b>first column holds the accession ID</b> as in the panel (for example <code>B001</code>) <b>or the accession name</b> (for example <code>HEIBIAO</code>); the other columns are traits.
Extra spaces, capital letters and punctuation do not matter (<code>b 001</code> still matches), unmatched IDs are listed in the results, and repeated IDs are averaged.
SNPs with more than 20 % missing calls or a minor-allele frequency below the filter are dropped before the analysis, so the number of SNPs analysed depends on the accessions you supply.
Files that still use the older <code>IRIS_313-…</code> identifiers are recognised too.</p>
{example("gwas/rice3k_phenotype_example.csv", rows=5, cols=3)}
{box("tip", "Not sure how the IDs are spelled? Click <b>Download ID template</b> on the page: it saves a CSV with all 3,024 accession IDs (with name, subpopulation and country) and an empty trait column — paste your phenotypes into it.")}
{ol(["Open <b>GWAS</b> and leave <b>Genotype source</b> on <i>Rice 3K panel</i>.", "Browse to your phenotype file and choose the trait to analyse.",
     "Choose the models (below) and click <b>Run GWAS</b>."])}
{figure("gwas_3k")}
<p>Behind the scenes TerraTools matches your IDs to the panel, keeps the accessions you have phenotypes for, and writes a <b>GWAS-ready file</b> (ID, phenotype, one column per SNP named like <code>Chr5_27558337</code>).
You can download that file from the Results tab and reuse it. SNP positions are derived from the marker IDs, which are positions along the 12 chromosomes laid end to end (IRGSP-1.0).</p>
<h3>Option 2 — your own VCF</h3>
{example("gwas/gwas_example.vcf.gz", rows=4, cols=12, caption="VCF (excerpt). The sample names in the header row must match the IDs in the phenotype file.")}
{example("gwas/gwas_vcf_phenotype_example.csv", rows=4, cols=3, caption="Phenotype file for the VCF: first column = sample ID (as in the VCF header).")}
{ul(["Only <b>biallelic SNPs on numbered chromosomes</b> are used (<code>1, 2, …</code> or <code>Chr1, Chr2, …</code>); indels, multi-allelic sites and scaffolds are skipped and counted.",
     "Genotype calls must be <b>diploid</b> (<code>0/0</code>, <code>0|1</code>, <code>1/1</code>; phased or unphased). A call with a missing allele (<code>./1</code>) counts as missing. Haploid or polyploid calls, and alleles above 1, stop the run with the marker named; they are not converted. Samples that are pooled bulks (names such as <i>high_bulk</i>) are refused: a pool is not a diploid individual.",
     "SNPs with too many missing calls, or a minor-allele frequency below the filter, are dropped when the VCF is read and again in the quality-control step; both are listed in the QC table.",
     "A VCF with millions of SNPs is slow: filter it first (for instance keep MAF ≥ 5 %)."])}
{figure("gwas_vcf")}
<h3>Option 3 — a ready-made file</h3>
<p>One CSV: <b>ID</b>, <b>Phenotype</b>, then one column per SNP named <code>Chr&lt;chromosome&gt;_&lt;position&gt;</code> with the genotype values (blank or <code>NA</code> for missing) in the coding you choose under <b>Genotype coding</b>.</p>
{ul(["<b>0 / 1 / 2 dosage</b>: the number of copies of the counted allele.", "<b>−1 / 0 / 1 centred</b>: converted to 0 / 1 / 2 by adding 1.",
     "A file that holds only 0 and 1 is <b>ambiguous</b> (a 0/1 dosage, or a centred matrix with no homozygous alternate). TerraTools never guesses: you must choose the coding, and the choice, the conversion and the number of missing calls are written to the log and the report.",
     "Values outside the chosen coding stop the run. <code>-9</code>, <code>-99</code>, <code>-999</code> and <code>NA</code> are read as missing; add other codes under <b>Extra values to read as missing</b>. Sample IDs, marker names, chromosomes and positions are kept as they are in the file."])}
{example("gwas/gwas_ready_example.csv", rows=5, cols=6, tail_cols=0)}
{sample_files("gwas")}
<h3>Quality control</h3>
<p>TerraTools expects genotypes that are already filtered (pruned for linkage, minor-allele frequency applied). It adds no further filtering: only markers without any variation, or without a valid chromosome and position, are dropped, and each missing call is filled with the mean dosage of its marker (not haplotype imputation: it ignores linkage disequilibrium and shrinks the variance of markers with many missing calls). Your genotype matrix itself is never edited. If fewer than 2 SNPs remain, the analysis stops with a message. The QC stage table is in the CSV and Excel report. The relationship matrix is VanRaden: K = ZZ′ / (2 Σ p(1−p)) with Z the genotypes centred by twice the allele frequency, computed from the observed calls before imputation. Under <b>Advanced</b> it can use LD-pruned markers or leave-one-chromosome-out matrices. If K is not positive semi-definite, negative eigenvalues are set to 0 and nothing else is changed. Principal components are taken from K and only used when you give a model a number of PCs (0 by default).</p>
<h3>The models</h3>
{table(["Model", "What it is"], [
    ["OLS", "Marker-by-marker linear regression, optionally with PC covariates. With 0 PCs it is labelled <b>uncorrected</b>; with PCs, <b>structure-corrected</b>. It does not model relatedness."],
    ["MLM", "y = Xβ + Zu + e, u ~ N(0, σ²ₐK), fitted by REML. The variance ratio is estimated once on the null model and the variance components are re-estimated for each marker; t test. The REML value is checked against a dense computation; the report says whether the check passed and whether the variance ratio reached the edge of its search range."],
    ["EMMAX", "The same null model, with the variance ratio <i>and</i> the variance components fixed for all markers; a faster approximation."],
    ["GEMMA-style", "The variance ratio is re-estimated for every marker (REML with a Wald test, or maximum likelihood with a likelihood-ratio test). It follows the same model as GEMMA but is not GEMMA: agreement with GEMMA is claimed only if you run the comparison (Advanced tab, needs the GEMMA program) and the numbers are within the stated tolerances."],
    ["FarmCPU, BLINK", "Iterative multi-locus models, run by R packages that come with TerraTools: <b>FarmCPU by rMVP</b> and <b>BLINK by GAPIT</b>. TerraTools writes your genotypes, map and phenotype, calls the package with its own default settings (FarmCPU has one setting here, the number of iterations), and reads back each marker's p-value and effect. They use no kinship matrix, need physical marker positions (marker names such as Chr3_123456), and the packages give no standard error for BLINK. The percent variance explained is computed by TerraTools (marginal), because the packages do not return the markers they selected."],
    ["mrMLM, FASTmrMLM", "Run by the <b>mrMLM R package</b>, which is installed in the R that comes with TerraTools: TerraTools writes your genotypes (as −1 / 0 / 1), phenotype and the structure PCs in the package's file formats, calls the package, and reads its tables back. Which markers are QTNs, their effects and their LOD are the package's own results, and nothing is recomputed. The settings are the package's search radius (kb), the critical LOD (default 3) and the number of structure PCs. The package gives a first-step p-value for every marker and no standard error; the QTN list is judged by LOD, so no p-value is derived from it. The percent variance explained shown is computed by TerraTools from the joint model of the QTNs."]])}
<h3>Significance</h3>
<p>There is nothing to set. Every run reports Bonferroni (α = 0.05 over the markers tested) and Benjamini–Hochberg separately, and, for OLS, MLM, EMMAX and GEMMA-style, a threshold simulated under the model's own null fit (1000 replicates; the report gives its interval). <i>Genome-wide significant</i> follows Bonferroni. A marker with p ≤ 1/m that is not genome-wide is <i>Suggestive</i>. The 15 top-ranked markers of each model that pass nothing are listed separately as <i>Exploratory</i> and are never called significant. If no marker passes, the result says <b>No genome-wide significant associations detected</b>. Multi-locus models keep their own QTN criterion.</p>
<h3>Reading the results</h3>
{figure("gwas_results")}
{table(["Output", "Meaning"], [
    ["Model summary", "Per model: markers tested, genome-wide and suggestive counts, the three thresholds, λ, the smallest p, and a note when a model did not run."],
    ["Associations", "One row per marker in three categories (genome-wide / QTN, suggestive, exploratory): model, SNP, chromosome, position, alleles, MAF, effect, SE, raw p, Bonferroni p, BH p, LOD, PVE and the criterion that applied. A multi-locus QTN keeps the LOD of its own model; no p-value is derived from that LOD."],
    ["Diagnostics", "Per model: SNPs tested, invalid statistics, λ, smallest p, the p-value histogram, convergence and the REML check. λ is reported only. TerraTools does not call it good or bad and never rescales p-values with it (no genomic control)."],
    ["Excel report", "Summary, associations, QC stages, diagnostics, one sheet per model with every SNP, and a Notes sheet with the definitions."],
    ["MTA table (CSV), GWAS-ready input, All plots (ZIP)", "The association rows, the file TerraTools built from your phenotypes, and the plots as PNG files."]])}
{figure("gwas_plots")}
<h3>Percent variance explained (PVE)</h3>
{table(["Model", "Definition (fraction of the phenotypic variance)"], [
    ["OLS", "(RSS without the marker − RSS with it) / total sum of squares, covariates in both fits."],
    ["MLM, EMMAX, GEMMA-style", "b²v / (b²v + σ²ₐ·mean(diag K) + σ²ₑ), b the marker effect, v its sample variance, the variances from the null REML fit. A marginal share relative to the variance model; it is not a joint estimate."],
    ["Multi-locus models", "Conditional: what the marker adds to the joint linear model of all selected markers and covariates."]])}
<p>Where an estimate cannot be defined (singular design, invalid statistic) it is shown as <b>unavailable</b>, not as 0. PVE from different models is not on a common basis and is not claimed to equal the value another program reports.</p>
<p>The <b>Manhattan plot</b> shows −log₁₀ p along the genome with the Bonferroni and suggestive lines (genome-wide hits filled, exploratory markers hollow), the <b>QQ plot</b> the p-value distribution and λ, and the <b>circos plot</b> the genome-wide markers on a ring of chromosomes.</p>
{box("note", "On a strongly structured panel such as rice (indica, japonica, aus) OLS without PCs is expected to give a λ far above 1 and many false associations; the mixed models are meant to remove them. This is checked by simulation: see <i>tools/gwas/validation</i> in the program folder, which reports false-positive rate, power and p-value calibration of every model on rice-like simulated data.")}
"""


def ch_clustering():
    return f"""
<p class="lead">Groups genotypes with principal-component analysis and K-means, then splits them into training and testing sets so that every cluster is represented in both — the standard way to
set up a fair genomic-selection experiment.</p>
<h3>What you need</h3>
<p>A CSV or TSV: the first column the sample ID, then one column per SNP (0 / 1 / 2). Missing values may be written as <code>-1</code>.</p>
{example("clustering/clustering_genotypes_example.csv", rows=5, cols=8)}
{sample_files("clustering")}
<h3>Step by step</h3>
{ol(["Open <b>Clustering &amp; Genomic Selection</b> and load the genotype file.", "Set the number of clusters <i>K</i> and the fraction of each cluster to keep for testing.", "Click <b>Run Clustering &amp; Split</b>."])}
{figure("clustering_input")}
<h3>Reading the results</h3>
{figure("clustering_results")}
{ul(["<b>The plot</b> shows the samples on the first two principal components, coloured by cluster; the axis labels give the variance each component explains.",
     "<b>Training Data</b> and <b>Testing Data</b> (CSV) are the two genotype files, ready to load in BreedSight. The split is <i>stratified</i> by cluster; if a cluster is too small for that, the tool says so and falls back to a random split.",
     "<b>Cluster membership</b> lists the cluster of every sample."])}
{box("note", "K is capped automatically so that it never exceeds the number of samples minus one.")}
"""


def ch_additive():
    return f"""
<p class="lead">Builds genomic relationship matrices from SNP dosages: the <b>additive</b> matrix (rrBLUP <code>A.mat</code>) and the <b>dominance</b> matrix (sommer <code>D.mat</code>).
They are the building blocks of GBLUP and of kinship-aware selection.</p>
<h3>What you need</h3>
<p>One CSV or TSV: <b>rows are samples, columns are SNP markers</b>, the first column holds the sample ID. Markers are coded <code>0</code> (AA), <code>1</code> (Aa), <code>2</code> (aa). Blank cells are treated as missing.</p>
{example("additive_dominance/snp_dosage_example.csv", rows=5, cols=8)}
{sample_files("additive_dominance")}
<h3>Step by step</h3>
{ol(["Open <b>Additive Dominance</b> and load the SNP file.", "Choose how missing values are filled (mean, median or zero) and tick the matrices you want.", "Click <b>Compute Matrices</b>."])}
{figure("additive_input")}
{ul(["Markers with no variation are removed automatically; at least two samples are needed.", "Missing values are filled <i>per marker</i> with the mean (default), the median or zero."])}
<h3>Reading the results</h3>
{figure("additive_results")}
{ul(["<b>Heatmaps:</b> sample-by-sample relationship; blocks along the diagonal are groups of related lines.",
     "<b>Histogram:</b> the distribution of the off-diagonal values — a peak near zero means a mostly unrelated set.",
     "<b>Data Table / Downloads:</b> the matrices as CSV (square, sample names in the first column and the header)."])}
"""


def ch_breedsight():
    return f"""
<p class="lead">Genomic prediction with <b>BreedSight</b>, a hybrid of a deep neural network and a random forest, benchmarked against Lasso, GBLUP_1 and RRBLUP_1 with cross-validation.
It tells you how well the genotypes predict the trait and predicts the breeding value of the test set.</p>
<h3>What you need</h3>
<p>Four CSV files. The sample IDs in a phenotype file and its genotype file must be the same.</p>
{example("genomic_prediction/train_phenotypes_example.csv", rows=4, cols=2, caption="Training phenotypes: <code>sample_id</code>, <code>phenotypes</code>.")}
{example("genomic_prediction/train_genotypes_example.csv", rows=4, cols=7, caption="Training genotypes: <code>sample_id</code>, then one column per SNP (0 / 1 / 2).")}
{example("genomic_prediction/test_phenotypes_example.csv", rows=4, cols=2, caption="Testing phenotypes: the same layout; the <code>phenotypes</code> column is optional for new, unphenotyped lines.")}
{example("genomic_prediction/test_genotypes_example.csv", rows=3, cols=7, caption="Testing genotypes: the same SNP columns as the training file.")}
{box("warn", "<b>No missing values.</b> Impute the genotypes before loading them (the Additive Dominance tool's imputation is one option). The test file must have exactly the same SNP columns as the training file.")}
{sample_files("genomic_prediction")}
<h3>Step by step</h3>
{ol(["Open <b>BreedSight GUI</b> and load the four files.", "Set the cross-validation folds and the trait's heritability h² (used only to put the accuracy in context).",
     "Keep the network settings at their defaults for a first run. Click <b>Run BreedSight Analysis</b>. The live log shows each fold; a run trains the network once per fold plus a final time, so it can take many minutes on larger data."])}
{figure("breedsight_input")}
{table(["Setting", "What it does"], [
    ["Number of CV folds / h²", "Cross-validation folds; trait heritability (the ceiling for prediction accuracy is about √h²)."],
    ["Epochs, batch size, learning rate", "Training schedule of the neural network."],
    ["DNN weight α", "The ensemble prediction is α × DNN + (1 − α) × random forest."],
    ["L2 regularization, dropout", "Guard against over-fitting."],
    ["RF trees / max depth", "Size of the random forest (0 = unlimited depth)."],
    ["PCA / RFE / raw genotypes", "Feature engineering: PCA of the relationship matrix, recursive feature elimination, or the raw genotypes."]])}
<h3>Reading the results</h3>
{figure("breedsight_results")}
{box("note", "The screenshots use the small example data with 2 folds and 50 epochs so that the run takes minutes, which is why the accuracy is modest and the program rates it “Very poor”. Read the page the same way with your own data: compare the models with each other, look at the validation r (not the training r — a high training r with a low validation r means over-fitting), and compare the best r with the theoretical maximum √h².")}
{ul(["<b>Metrics &amp; Interpretation:</b> R², Pearson r and RMSE for every model on the training and validation folds, and a plain-language reading of how close the best model is to the theoretical maximum for your heritability.",
     "<b>Test Predictions:</b> the predicted value for every test sample (and the observed value if you supplied it).",
     "<b>Plots:</b> regression of predicted against observed, model comparison, per-fold results, a performance heat map and the phenotype distribution.",
     "<b>Downloads:</b> model metrics, test predictions and validation predictions as CSV."])}
"""


def ch_cross():
    return f"""
<p class="lead">Decides which parents to cross. From the marker data of your candidates it predicts the <b>usefulness</b> of every possible cross — how good the best progeny of the cross are expected to be — and finds the plan of crosses that gives the most while keeping diversity. It compares the plan with simply mating the best parents.</p>
<h3>What you need</h3>
<p>Three files. The relationship matrix between the candidates is <b>calculated from the dosage matrix</b>, so there is nothing else to upload.</p>
{ul(["<b>Allele dosage matrix</b> — genotypes in the rows, markers in the columns, values 0 / 1 / 2. The parents are taken as inbred lines (0 and 2); a heterozygous call (1) is read as missing, and the number of such calls is reported.",
     "<b>Marker effects</b> — the columns <code>Marker</code> and <code>Effect</code>: the additive effect of every marker (for example from a genomic prediction model). The genomic value of each candidate is its dosage × these effects.",
     "<b>Genetic map</b> — the columns <code>Marker</code>, <code>Chromosome</code> and <code>cM</code>. The distance between markers sets how much they recombine, and therefore how much the progeny of a cross can vary (the Haldane map function is applied). <b>No map? Use a VCF instead:</b> choose <i>No map: use a VCF</i> and give the VCF and a recombination rate (below)."])}
{example("cross_selection/dosage_matrix_example.csv", rows=4, cols=6, caption="<b>Allele dosage matrix</b>: the genotype name in the first column, then one column per marker.")}
{example("cross_selection/marker_effects_example.csv", rows=4, cols=2, caption="<b>Marker effects</b>.")}
{example("cross_selection/genetic_map_example.csv", rows=4, cols=3, caption="<b>Genetic map</b>.")}
{box("warn", "The marker names must be identical in the three files; markers that are not in all three are left out and the number is reported. Preview &amp; Validate shows how many genotypes and markers remain.")}
<h4>No genetic map: a VCF and a recombination rate</h4>
<p>A VCF has the position of every marker on the genome in base pairs (<code>CHROM</code> and <code>POS</code>), but the progeny variance needs positions in cM. Choose <b>No map: use a VCF</b>, select the VCF (plain or <code>.gz</code>) and enter the <b>recombination rate in cM per Mb</b>: each position is turned into cM as <i>bp / 1,000,000 × rate</i>. A marker is found in the VCF by its <code>ID</code>, or by <code>CHROM_POS</code> / <code>CHROM:POS</code> when the ID is empty, and must carry the same name as in the dosage matrix and the marker effects. The VCF is used only for the positions; the genotypes still come from the dosage matrix.</p>
{box("warn", "One rate for the whole genome is an approximation (roughly 4 cM per Mb in rice, 1 in maize, 0.2 in wheat; use the value of your species). Real recombination is very low around the centromere and high toward the chromosome ends, so the progeny variance, and with it the usefulness, is only approximate. A genetic map from <b>Linkage Mapping</b> is more accurate; the result states which one was used.")}
{figure("cross_markers_vcf")}
{sample_files("cross_selection")}
<h3>Settings</h3>
{ul(["<b>Number of crosses</b> you can make, the <b>maximum crosses per parent</b> and the <b>penalty degree</b> (see below).",
     "<b>Population made from each cross</b> — doubled haploids (DH) or recombinant inbred lines (RIL) — and its <b>generation</b> (1 = DH made from the F1, or RILs taken from the F2; 11 or more for fully inbred RILs).",
     "<b>Selected proportion</b> — the share of the progeny you will keep (0.05 = the best 5 %); it sets the selection intensity in the usefulness.",
     "<b>Candidate pool</b> — the number of best candidates, in order of genomic value, whose crosses are all evaluated (60 candidates mean 1,770 crosses). A larger pool takes longer."])}
<h3>Step by step</h3>
{ol(["Open <b>Optimal Cross Selection</b>, load the dosage matrix, the marker effects and the genetic map (or the VCF and its rate) and click <b>Preview &amp; Validate</b>.", "Set the number of crosses, the maximum crosses per parent, the penalty degree and the population settings.", "Click <b>Run Cross Selection</b>. The <b>Optimal Plan</b> tab opens when it is done."])}
{figure("cross_markers_input")}
<h3>How the optimal plan is found</h3>
<p>For every cross among the candidates the tool calculates, with the R package SimpleMating (Lehermeier et al. 2017), the <b>mean</b> of the progeny from the marker effects, its <b>standard deviation</b> from the effects and the recombination along the map, and the <b>usefulness = mean + selection intensity × standard deviation</b>. A cross of two average parents that differ at many large-effect markers can have a higher usefulness than a cross of two good parents that are alike.</p>
<p>A plan is a list of crosses, each parent used in at most <i>max crosses per parent</i> of them. A parent in <i>n</i> of the <i>N</i> crosses contributes <i>x = n / 2N</i> to the next generation. A plan is scored by</p>
{ul(["<b>Expected usefulness</b> — the mean usefulness of the crosses of the plan.", "<b>Group coancestry</b> — <i>x′Cx</i>, the relatedness of the next generation, from the relationship matrix calculated from the dosage matrix (VanRaden 2008; half of it is the coancestry); the lower, the more diversity is kept.", "<b>Progeny inbreeding</b> — the mean coancestry between the two parents of each cross."])}
<p>The <b>penalty degree</b> sets the balance: 0° maximises usefulness only, 90° minimises coancestry only. Usefulness and coancestry are rescaled to a common range and the plan with the highest <i>cos(angle) × usefulness − sin(angle) × coancestry</i> is searched by simulated annealing, which changes one parent of one cross at a time and keeps the plan valid (no selfing, no repeated cross, no parent over its limit). The search is repeated at 0°, 1°, 5°, 10° … 90° to draw the trade-off, and the longer search is made at your angle. The best angle is a choice: read the trade-off, then pick.</p>
{box("note", "The marker effects carry no intercept, so the values are on a relative scale: only the differences between crosses matter. Negate the effects for a trait where lower is better. Negative relationships are set to 0.")}
<h3>Reading the results</h3>
{figure("cross_markers_results")}
{ul(["<b>Optimal Plan:</b> the crosses to make (with each cross's <i>Usefulness</i>, <i>Cross Mean</i>, <i>Progeny SD</i> and <i>Kinship Coefficient</i>), each parent's genomic value, number of crosses and contribution, and the usefulness / coancestry of the plan next to <i>truncation selection</i> (the best parents mated at random).",
     "<b>Usefulness versus diversity:</b> the chart and table show the trade-off at every penalty degree.",
     "<b>Downloads:</b> crossing plan, parent contributions, trade-off table (CSV) and the chart, the usefulness of all crosses, and the relationship matrix that was calculated."])}
{box("tip", "A cross between two top parents is not always best: if they are close relatives the offspring lose diversity. Compare your plan with truncation selection and move the penalty degree until the usefulness you give up buys enough diversity.")}
"""


def ch_heterosis():
    return f"""
<p class="lead">Trains a deep-learning model on a training population, predicts the genomic estimated breeding value (GEBV) of every possible hybrid between your parents, and screens them for
mid-parent (MPH) and better-parent (BPH) heterosis.</p>
<h3>What you need</h3>
<p>Three CSV files, plus one optional:</p>
{example("heterosis/train_phenotype_example.csv", rows=4, cols=2, caption="<b>Training phenotype</b>: column 0 = sample ID, column 1 = the numeric trait value (header row required).")}
{example("heterosis/train_genotype_example.csv", rows=4, cols=7, caption="<b>Training genotype</b>: first column = sample ID (the row label), then numeric marker columns.")}
{example("heterosis/parent_genotype_example.csv", rows=4, cols=7, caption="<b>Parent genotype</b>: the candidate parents, with exactly the same marker columns as the training genotype.")}
<p>The optional <b>test genotype</b> file has the same marker columns and gives predictions for extra lines.</p>
{box("warn", "<b>Row order matters.</b> The training phenotype and training genotype files are paired by <i>row position</i>, not by ID: row 5 of one must be the same individual as row 5 of the other. Sort both files the same way, and keep the marker columns identical in all genotype files. Missing values are not allowed.")}
{sample_files("heterosis")}
<h3>Step by step</h3>
{ol(["Open <b>Heterosis</b> and load the files — a short note under each box shows its size.", "Choose how many top hybrids and parents to report and the cross-validation settings; epochs and folds mostly control the run time.", "Click <b>Run Prediction</b>. The log shows each fold."])}
{figure("heterosis_input")}
<h3>Reading the results</h3>
{figure("heterosis_results")}
{box("warn", "<b>Check the cross-validation numbers before trusting the hybrid ranking.</b> The network learns the trait on its raw scale, so a trait with a large mean (for example 2,000 kg/ha) needs many more epochs than one with a small mean, and a short run gives predictions that are systematically too low or too high — which also distorts MPH and BPH. Work with a trait on a moderate scale (for example in t/ha, or centred), keep at least the default 100 epochs, and look for a clearly positive Pearson r and an R² near or above zero. The <i>ranking</i> of hybrids can be useful even when the absolute values are off.")}
{ul(["<b>Cross-validation summary:</b> MSE, RMSE, R² and Pearson r averaged over the folds, and per fold. A Pearson r near zero means the markers do not predict this trait.",
     "<b>Top F1 hybrids:</b> hybrids ranked by predicted GEBV, with <i>MPH</i> and <i>BPH</i>. Both are computed against the <b>training population</b> — MPH against its mean, BPH against its best value — so read them as 'better than the training population by x %'.",
     "<b>Top parents:</b> parents ranked by their mean GEBV over all hybrids they take part in.",
     "<b>Test genotype predictions</b> (when a test file was loaded), and the hybrid-performance plot."])}
"""


def ch_phenotypic():
    return f"""
<p class="lead">Analyses multi-environment field trials laid out as an augmented RCBD, an RCBD or an alpha-lattice: descriptive statistics, ANOVA, mixed models with BLUEs, variance components and heritability,
correlations, clustering, PCA and genotype-by-environment stability (including AMMI).</p>
<h3>What you need</h3>
<p>A CSV or Excel file with <b>one row per plot</b>. Columns are found by their <b>name</b>, in any order (capital letters and spaces around a name do not matter). There are two kinds of columns.</p>
{table(["Column name", "Needed for", "What it holds"], [
    ["<b>Genotype</b>", "all designs", "Name of the variety or line. Checks and new entries both go here; a genotype that appears more than once in an environment is treated as replicated."],
    ["<b>Environment</b>", "all designs", "Location, season or year label. At least two different environments are needed."],
    ["<b>Rep</b>", "all designs", "Replicate number. Augmented RCBD: use 1 in every row if the trial is not replicated. RCBD and alpha-lattice: each genotype appears once in every Rep."],
    ["<b>Block</b>", "augmented RCBD, alpha-lattice", "Block number. Augmented RCBD: every block holds all the checks and some new entries. Alpha-lattice: the incomplete block, numbered inside its Rep (1, 2, 3 … again in every Rep). An RCBD does not need it."],
    ["<b>Any other column</b>, for example Yield, PlantHeight, DaysToFlower", "—", "A <b>trait</b>. Every other numeric column is a trait. Leave a cell <b>empty</b> (or write <code>NA</code>) when a trait was not measured. A <code>0</code> is also read as <i>not measured</i> (a note lists where this happened); switch this off with <b>Treat 0 as missing</b> if zeros are real values."],
    ["Plot, Row, Column, Range, Year, Season, Entry, Line.type, Remarks", "—", "Layout and label columns. They are not treated as traits."]])}
{example("phenotypic/multi_environment_trial_example.csv", rows=5, cols=12)}
{ul(["At least <b>3 genotypes</b> and <b>2 environments</b> are needed.", "Trait columns must be numeric. The tool reads the file once and lists the traits it found (those with at least one value); every analysis tab then lets you tick the traits it should use.",
     "The file can hold other columns in any position (a plot number, a field map, comments). If a column the design needs is missing, a warning under the file name names it.",
     "Each genotype should appear in each environment (ideally repeated in reps); unbalanced data still run but give less precise estimates."])}
<p>The first rows of each design look like this (the Plot, Yield and PlantHeight columns are examples; you can have as many trait columns as you like):</p>
{table(["Design", "Genotype", "Environment", "Rep", "Block", "Plot", "Yield"], [
    ["Augmented RCBD", "Check1", "Env1", "1", "1", "101", "5.4"], ["", "New07", "Env1", "1", "1", "102", "6.1"], ["", "Check2", "Env1", "1", "1", "103", "5.3"], ["", "New21", "Env1", "1", "2", "201", "5.8"],
    ["RCBD", "Geno05", "Env1", "1", "(none)", "101", "6.8"], ["", "Geno03", "Env1", "1", "(none)", "102", "7.8"], ["", "Geno05", "Env1", "2", "(none)", "201", "6.5"],
    ["Alpha-lattice", "Geno05", "Env1", "1", "1", "101", "6.8"], ["", "Geno03", "Env1", "1", "1", "102", "7.8"], ["", "Geno12", "Env1", "1", "2", "105", "5.9"], ["", "Geno05", "Env1", "2", "3", "201", "6.5"]])}
{sample_files("phenotypic")}
<h3>Experimental design</h3>
<p>Choose the layout of your trial under <b>Experimental Design</b>. It decides which blocking terms the ANOVA, the mixed models (BLUEs) and the heritability use.</p>
{table(["Design", "Columns needed", "How the file should look", "Blocking terms"], [
    ["Augmented RCBD", "Genotype, Environment, Rep, Block", "Checks repeated in every block, new entries once. Use one Rep (all 1) if the trial has no replicates.", "Block, fixed. Rep is used only if every Rep + Block group holds a check"],
    ["RCBD", "Genotype, Environment, Rep", "Every genotype once in each replicate; the replicate is the block. The Block column is not used.", "Rep"],
    ["Alpha-lattice", "Genotype, Environment, Rep, Block", "Incomplete blocks inside each replicate. Number the blocks inside each Rep (1, 2, 3 … in every Rep); block 1 of Rep 1 and block 1 of Rep 2 are treated as different blocks.", "Rep, and Block within Rep"]])}
{ul(["The <b>augmented RCBD</b> follows the classical augmented-design analysis. The checks are the genotypes that occur more than once in the environment. A block is the value of the Block column; if the file has several Reps and every Rep + Block group holds a check, the blocks are the Rep + Block groups instead (a Rep column that only numbers the repeats of the checks is ignored, and the Design Check tab says so). Block effects are estimated from the checks alone, every genotype's adjusted mean (BLUE) is its value corrected for its block, and the error term has (checks − 1) × (blocks − 1) degrees of freedom. The ANOVA of each environment has two parts — treatments adjusted for blocks, and treatments ignoring blocks, each split into checks, new entries and checks versus new entries — and the <b>Design Check</b> tab adds the standard errors of the differences between two checks, two new entries in the same block, two new entries in different blocks, and a new entry and a check. If the checks are not repeated across the blocks (or every block holds the same genotypes) the genotypes cannot be separated from the blocks; the analysis then treats the blocks as random and says so in the Design Check tab.",
     "With <b>RCBD</b> and <b>alpha-lattice</b> the replicates (and blocks) are treated as nested in their environment, so “Rep 1” of one environment is not the same as “Rep 1” of another. The pooled mixed model (BLUEs across environments) and the heritability use random blocks for all designs. The pooled ANOVA table of the RCBD and alpha-lattice is a true analysis-of-variance table (degrees of freedom, sums of squares, F and p); in the alpha-lattice the genotype line is adjusted for the incomplete blocks.",
     "The <b>Design Check</b> tab lists, for each environment, the number of plots, genotypes, replicates and blocks, the block size and how often each genotype occurs. It also warns when the file does not look like the chosen design — for example genotypes missing from a replicate in an RCBD, a single block per replicate in an alpha-lattice, or no repeated checks in an augmented design."])}
<p>The example files <code>augmented_rcbd_example.csv</code>, <code>rcbd_example.csv</code> and <code>alpha_lattice_example.csv</code> show the layout of each design; <code>multi_environment_trial_example.csv</code> is a three-replicate trial in which every genotype is in every replicate; analyse it as an RCBD (under the augmented design it falls back to random blocks and says so).</p>
{box("note", "Pick the design that matches how the trial was planted. The same file can be analysed with any of them, but a wrong choice gives wrong error terms: an RCBD model on an augmented trial ignores the block effects, and an alpha-lattice model needs Block numbers that really are incomplete blocks.")}
<h3>Step by step</h3>
{ol(["Open <b>Phenotypic Data Analysis</b> and load the file in <b>Data &amp; Design</b>. A summary shows the design, plots, accessions, environments and traits, with a table of the traits and the design check.", "Choose the <b>experimental design</b>, whether to analyse all environments or only some, whether 0 counts as missing, and the colour scheme.",
     "Open <b>Analyses</b>, pick an analysis (Descriptives, Environments, BLUEs, ANOVA, Heritability, Correlation, PCA, Clustering or Stability &amp; G×E), choose its traits and options and click its <b>Run</b> button. Every analysis runs on its own, so you can run only what you need; results stay until you change the data settings.",
     "Use <b>Downloads</b> for the Excel workbook of everything run so far, the ANOVA workbook, all figures (300-dpi PNG and vector PDF) or everything in one zip."])}
{figure("phenotypic_input")}
<h3>Reading the results</h3>
{figure("phenotypic_results")}
{table(["Analysis tab", "What you get"], [
    ["Descriptives", "Mean, SD, CV, range, quartiles and shape of every trait (all environments, then each), box and violin plots per environment and the distribution of one trait."],
    ["Environments", "Accession means and ranks in every environment, how the ranks change, and the correlation between environments."],
    ["BLUEs", "Adjusted means (BLUEs) of every accession within each environment and over all environments from a mixed model; the standard errors of differences are added to the Design Check."],
    ["ANOVA", "One analysis-of-variance table per trait with every environment and the pooled analysis side by side (DF, SS, MS, F, P, CV), a plot of where the variation comes from, and residual and normal Q-Q plots of one trait."],
    ["Heritability", "Genotypic, G×E, block and error variance components and heritability on a plot and on an entry-mean basis, per trait, environment and pooled."],
    ["Correlation, PCA, Clustering", "Relationships between traits (Pearson, Spearman or Kendall with P values), a PCA biplot and variance explained, and clusters of accessions (K-means or Ward hierarchical) with their trait means."],
    ["Stability &amp; G×E", "Wricke ecovalence, Finlay-Wilkinson slope, AMMI scores and ASV of every accession, with the stability, G×E and AMMI1 biplots."]])}
<p>Every figure can be saved as a 300-dpi PNG or as a vector PDF, and every table as a CSV. A message box above the results lists notes, warnings and errors of that analysis (for example when an environment has too little data).</p>
{figure("phenotypic_anova")}
{figure("phenotypic_heritability")}
{box("tip", "Heritability is a property of the trial, not only of the trait: it grows with replication and with the number of environments. Compare it between trials with that in mind.")}
"""


def appendix_workflows():
    return f"""
<h3>A. From raw marker calls to QTL</h3>
{ol(["<b>Genofile Converter</b> — turn the HapMap, VCF or A:G file of the genotyping company into the 0 / 1 / 2 files (name the two parents).", "<b>Mendelian Segregation</b> — screen the converted markers for distortion (optional).",
     "<b>Linkage Mapping → Parameters → Run</b> — build the map; check the QC report.", "<b>Linkage Mapping → Results → File for QTL Mapping</b> — download <i>QTL Mapping input (CSV)</i>.",
     "<b>QTL Mapping</b> — load that file as the genotype CSV and your phenotype CSV; run CIM with the permutation threshold.", "<b>Genes Extraction</b> — use the QTL Mapping output, the genetic map (or a VCF) with the base-pair position of its markers and your GFF3 annotation to list the candidate genes of each QTL.", "<b>Gene Prioritization</b> — give the gene ID list that Genes Extraction wrote (rice, wheat or maize), and the trait; read which candidates the databases support most."])}
<h3>B. GWAS on the rice 3K panel</h3>
{ol(["Download the <b>ID template</b> on the GWAS page and paste in your phenotypes.", "Choose the trait and the models (MLM and a multi-locus model such as FarmCPU or BLINK).", "Run; read the MTA table, then the Manhattan and circos plots; keep the Excel report."])}
<h3>C. Genomic selection</h3>
{ol(["<b>Clustering &amp; Genomic Selection</b> — split the lines into training and testing sets.", "<b>BreedSight GUI</b> — train and validate; check that the accuracy is useful.",
     "<b>Additive Dominance</b> — compute the kinship matrix of the candidate parents.", "<b>Optimal Cross Selection</b> — with the marker effects, the dosage matrix and a genetic map, choose the crosses."])}
<h3>D. Bulk segregant analysis</h3>
{ol(["Call variants for the parent and the two bulks.", "<b>QTL-Seq</b> with the two bulks (and the parents, when you have them) as sorted BAM files or VCF files: Δ SNP-index, simulated confidence limits and QTL intervals.", "Alternative: <b>PYBSA G-statistics</b> with the two GATK tables, for a G-statistic analysis with a report."])}
<h3>E. A figure of the population's chromosomes</h3>
{ol(["Prepare a file with the chromosome, the marker name, the position (bp) and one column per sample, parents included (or use a VCF).", "<b>Genotype Viewer</b> — choose the species, <i>Compare</i> mode for one chromosome or <i>Normal</i> for one image per sample, a bright or mixed palette and, for a paper, 600 dpi plus the vector PDF.", "Type the parents' names into <i>Legend names</i> and download the images."])}
<h3>F. Choosing parents by their trait alleles</h3>
{ol(["<b>Favourable Allele Analysis</b> — choose the panel and the traits you breed for.", "Select the candidate parents and read the cards and the heat map: which favourable alleles each carries and where they differ.", "Take the pair that complements best into <b>Optimal Cross Selection</b>, or cross them directly."])}
"""


def appendix_troubleshooting():
    return table(["What you see", "Likely cause and what to do"], [
        ["“Please select …” right after clicking Run", "A required file or choice is missing. Check every file box and drop-down on the input tab."],
        ["Result tables are empty, or “0 matching samples”", "The IDs in your files are not written identically (case, spaces, <code>Chr3</code> vs <code>3</code>). Open both files side by side and compare the first IDs."],
        ["“Phenotype column not found”, “Trait not found”", "The phenotype file has a different header than you expect, or the trait column is not numeric. Re-select the file after fixing it."],
        ["QTL Mapping: “Only n matching samples”", "The first column of the phenotype file and the first column of the genotype rows do not share IDs, or fewer than 8 do."],
        ["Genofile Converter: every marker is left out, or almost none are coded", "Parent 1 and Parent 2 are not the two parents (or are the same line): a marker is only coded when both parents are called, homozygous and carry different alleles. Pick the right samples; the <i>Markers left out</i> table gives the reason of every marker."],
        ["Linkage Mapping: a marker is missing from the map", "Look at the QC Report: it lists every removed marker with the reason (too much missing data, monomorphic, distorted, chromosome 0)."],
        ["Linkage map is far longer than expected", "Markers are not in map order in the input file, or markers from different regions are on one chromosome. Fix the order and the chromosome column."],
        ["GWAS: “Only 0 of n phenotype IDs match”", "For the 3K panel use the IDs from the <b>Download ID template</b> button. For a VCF use the sample names from the VCF header row."],
        ["GWAS with OLS finds hundreds of MTAs", "OLS ignores population structure. Switch to MLM or a multi-locus model and compare the genomic inflation λ."],
        ["GWAS: “genotype coding is ambiguous”", "The file holds only 0 and 1. Choose <i>0/1/2 dosage</i> or <i>−1/0/1 centred</i> under Genotype coding; TerraTools does not guess."],
        ["GWAS: “Only n SNP(s) pass quality control”", "Fewer than 2 variable markers with a valid chromosome and position were left. Check the genotype file and its coding."],
        ["QTL-Seq: “samtools was not found”", "Reading BAM files is much faster with samtools (Linux / macOS, or WSL on Windows). Without it TerraTools uses pysam if it is installed, otherwise a slower built-in reader; for a whole genome start from a VCF or install samtools."],
        ["QTL-Seq: no QTL is reported", "No window passed the limit of the chosen level. Check the number of plants per bulk and the population type (they set the limits), lower the level to 90 %, or check that the bulks were not swapped."],
        ["Gene Prioritization: “No internet connection …” or a service shown as not reachable", "The tool asks Ensembl Plants, Gramene, Plant Reactome, NCBI and (rice) Oryzabase. Check the connection (a proxy or firewall can block them) and press <b>Check connections</b>. With only some services down the run still works and says which evidence is missing."],
        ["Gene Prioritization: IDs “not matched”", "Use the IDs of the Ensembl Plants / Gramene annotation (rice RAP-DB <code>Os02g0766700</code>, wheat <code>TraesCS…</code>, maize <code>Zm00001eb…</code>). Rice MSU IDs (<code>LOC_Os…</code>) are translated only for the genes Oryzabase lists."],
        ["Gene Prioritization: most genes have a score of 0 or “—”", "The databases hold little for those genes (no trait annotation, description or pathway). The score is then low or unavailable; it is not a judgement on the gene."],
        ["Genotype Viewer: “no markers on the chromosomes of this species”", "The chromosome names or the species do not fit. Names such as <code>1</code>, <code>chr1</code> and <code>chr01</code> are accepted; check that the species is the one the markers were called on, or add the species from its FASTA file."],
        ["Genotype Viewer: some markers are not drawn", "They lie on chromosomes or scaffolds the species table does not have (the Results tab gives the number). This is normal for unplaced scaffolds."],
        ["Genotype Viewer: the VCF takes a long time", "Every sample and every site is converted to a table first. Filter the VCF to the markers you need (for example with <code>bcftools view -R</code>) and keep biallelic sites only."],
        ["Genotype Viewer: the colours do not match the parents", "The viewer only knows REF and ALT. Code the file so that parent 1 is <code>A</code> / <code>0|0</code> and parent 2 is <code>B</code> / <code>1|1</code>, and name them in <i>Legend names</i>."],
        ["Favourable Allele Analysis: the page is blank", "The dashboard is built into the program and needs no files. Close the tool page and open it again; if it stays blank, restart TerraTools."],
        ["BreedSight / Heterosis is very slow", "They train neural networks. Lower the epochs and the folds for a trial run, and reduce the number of markers."],
        ["A plot looks empty or a table says “Not computed”", "That analysis was not applicable to the data (for example dominance without heterozygous calls). The message on the page says which."],
        ["The page looks stuck", "Look at the status line below the Run button — long analyses show their progress there. You can switch to another tool; a notification appears when the run ends."],
        ["macOS: “cannot be opened because the developer cannot be verified”", "System Settings → Privacy &amp; Security → <b>Open Anyway</b> (see Chapter 2)."]])


def appendix_glossary():
    return table(["Term", "Meaning"], [
        ["BLUE", "Best linear unbiased estimate: a genotype's mean after removing the effects of environment, rep and block."],
        ["cM (centiMorgan)", "Unit of genetic distance; 1 cM ≈ 1 % recombination between two markers for short distances. TerraTools uses the Kosambi mapping function."],
        ["Δ(SNP-index)", "Difference in the share of reads carrying the other parent's allele between two bulks; around zero where no QTL is near."],
        ["Favourable / unfavourable allele", "Favourable Allele Analysis: the allele of a trait gene that improves the trait (from the literature), and the opposite allele. A variety is favourable at a SNP when it is homozygous for the favourable allele."],
        ["F2 / RIL", "F2: second filial generation (1 : 2 : 1 genotypes). RIL: recombinant inbred line (mostly homozygous, 1 : 1)."],
        ["Evidence level", "Gene Prioritization: Direct association, Multiple lines, Two lines, Single line or No trait evidence, from how many evidence categories support the gene."],
        ["Final score", "Gene Prioritization: the final score (0–100) weighs what the databases know about a gene and the trait. It is not the probability that the gene is causal."],
        ["GEBV", "Genomic estimated breeding value: a line's predicted genetic value from its markers."],
        ["GenoSee", "The graphical-genotype program (Hashimoto 2024, MIT licence) on which the Genotype Viewer is based."],
        ["Heritability (H²)", "Share of the phenotypic variance that is genetic."],
        ["Kinship", "Genetic relatedness of two lines; 0 = unrelated, 1 = identical."],
        ["LOD", "Logarithm of odds. In QTL Mapping the likelihood-ratio LOD of a marker; in GWAS the 1-df χ² statistic divided by 2 ln 10. Higher = stronger evidence."],
        ["Lambda (λ)", "Genomic inflation factor of a GWAS: about 1 when p-values are well calibrated."],
        ["MAF", "Minor-allele frequency: the frequency of the rarer allele at a marker."],
        ["MPH / BPH", "Mid-parent / better-parent heterosis: how far a hybrid exceeds the parents' mean / the better parent, in percent."],
        ["MTA", "Marker–trait association: a marker that passes the significance threshold in GWAS."],
        ["Permutation threshold", "A significance level obtained by shuffling the phenotype many times and taking the 95th percentile of the highest genome-wide statistic."],
        ["PVE", "Percent of the phenotypic variance explained by a marker or QTL."],
        ["QTL", "Quantitative trait locus: a genome region that influences a quantitative trait."],
        ["SNP-index", "Share of a bulk's reads that carry the allele of the non-baseline parent."],
        ["Segregation distortion", "A marker whose genotype ratio departs from Mendel's expectation."]])


def appendix_formats():
    return table(["Tool", "Files", "Layout of the first rows / columns", "Example file"], [
        ["Genofile Converter", "1 HapMap, VCF / VCF.GZ or raw Excel / CSV", "HapMap: rs#, alleles, chrom, pos … then samples; VCF: #CHROM … FORMAT then samples (GT); raw: Marker, Chromosome, P1, P2, offspring (A:G …); the two parents are named in the tool", "genofile/*"],
        ["Phenotypic Data Analysis", "1 CSV/Excel", "Columns by name: Genotype, Environment, Rep, Block (not for an RCBD); every other numeric column is a trait", "phenotypic/*"],
        ["Mendelian Segregation", "1 CSV", "Marker, Chromosome, one column per individual (0/1/2); the QTL Mapping layout is also accepted", "mendelian/mendelian_RIL_example.csv"],
        ["Linkage Mapping", "1 Excel or CSV", "Marker, Chromosome, one column per individual (0/1/2/−1) — made by the Genofile Converter from company files", "linkage/linkage_F2_coded_example.xlsx"],
        ["Genotype Viewer", "1 CSV or VCF", "By position: chromosome, marker, position (bp), then one column per sample (A/H/B/N, 0|1 or 0/1); or a VCF / VCF.GZ", "genotype_viewer/*"],
        ["QTL Mapping", "2 CSV", "Genotype: 3 header rows (marker, chromosome, cM) then one row per individual; Phenotype: ID, traits", "qtl/qtl_genotypes_example.csv, qtl_phenotypes_example.csv"],
        ["QTL-Seq", "2 BAM + FASTA, or 2 VCF", "Bulk VCFs with GT:AD per sample (a parental VCF is optional)", "qtl_seq/qtlseq_example.vcf.gz"],
        ["Diversity Analysis", "1 CSV", "Variety, Sample ID, haplotype block codes (−1 = missing)", "diversity/haplotype_blocks_example.csv"],
        ["GWAS", "phenotype + (3K | VCF | ready file)", "Phenotype: ID, traits · ready file: ID, Phenotype, Chr1_12345 …", "gwas/*"],
        ["PYBSA G-statistics", "2 tables (TSV/CSV)", "CHROM, POS, REF, ALT, then &lt;sample&gt;.GT / .AD / .GQ", "pybsaseq/parents_example.tsv, bulks_example.tsv"],
        ["Genes Extraction", "CSV / Excel + GFF3", "QTL intervals (Chromosome, Start, End / QTL_ID, QTL_Start, QTL_End) or the QTL table of QTL Mapping · genetic map (Marker, Chromosome, bp) or, without a map, a VCF · GFF3 annotation · optional FASTA", "gene_extractor/*"],
        ["Gene Prioritization", "Gene IDs", "Gene IDs (TXT / CSV / Excel, first column, or pasted; the Gene_ID_list.txt of Genes Extraction) · internet connection", "gene_prioritization/*"],
        ["Clustering & Genomic Selection", "1 CSV/TSV", "SampleID, SNP columns (0/1/2, −1 missing)", "clustering/clustering_genotypes_example.csv"],
        ["Additive Dominance", "1 CSV/TSV", "SampleID, SNP columns (0/1/2, blank missing)", "additive_dominance/snp_dosage_example.csv"],
        ["BreedSight GUI", "4 CSV", "sample_id + phenotypes · sample_id + SNPs (no missing), for train and test", "genomic_prediction/*"],
        ["Optimal Cross Selection", "3 files", "Allele dosage matrix (genotypes × markers, 0/1/2) · Marker, Effect · Marker, Chromosome, cM — or a VCF and a recombination rate", "cross_selection/*"],
        ["Heterosis", "3–4 CSV", "pheno: ID, trait · genotypes: ID, marker columns (same columns, same row order)", "heterosis/*"],
        ["Favourable Allele Analysis", "1 file", "Rice 4K haplopanel genotypes: marker IDs (chr1_7268879) × lines, codes 0 / 1 / 2 / -1 — or a VCF", "favourable_allele/*"]])
