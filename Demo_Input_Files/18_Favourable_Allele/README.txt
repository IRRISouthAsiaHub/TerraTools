FAVOURABLE ALLELE ANALYSIS - favourable alleles of 166 known trait genes in your rice 4K haplopanel lines   (run time: the matching takes about 1 s, the heatmap about 5 s)

SCOPE   Rice 4K haplopanel only, for now. The 327 SNPs / indels of the 166 trait genes are kept inside TerraTools (not shown, not exported); the user only supplies a genotype
        file, and TerraTools matches its marker IDs with them. Not scaled yet to other marker panels, crops or genes.

INPUT   input/rice_4K_haplopanel_demo.csv   1,027 markers x 24 rice lines. The first column is the marker ID (chromosome_position, e.g. chr1_7268879), then one column per line.
        Genotypes are the codes 0 = REF/REF, 1 = heterozygous, 2 = ALT/ALT, -1 = missing.  327 of the markers are trait-gene markers, 700 are other markers of the panel
        (random positions and genotypes, to show that markers that are not trait-gene markers are ignored).  The 24 lines are real genotypes of improved varieties and landraces.
        The marker IDs may also carry the gene (chr1_7268879_LOC_Os01g13060), a VCF (GT field) or nucleotide genotypes (A/G) can be given instead, and chromosome / position /
        REF / ALT columns are used when the ID alone does not match.

STEPS   Genomic Selection > Favourable Allele Analysis > tab "Your 4K File" > Browse: rice_4K_haplopanel_demo.csv.
        1. Read the match report: 1,027 markers in the file, 327 / 327 trait-gene markers matched (all by chromosome_position), 166 / 166 trait genes found, 24 lines.
           "Download matched IDs (CSV)" saves the list of matched IDs.
        2. Tick lines to compare (the first three are selected), click "Show results", then look at Overview, Genes and SNP Table.
        3. Tab "Heatmap": add the genes of DSR trait, Floral and Hybrid fertility (or pick genes by name), "All lines", then "Draw Balloon Heatmap".

EXPECT  1. 327 / 327 matched IDs, 700 markers of the file ignored, 5,534 of 7,848 genotype calls present (71 %).
        3. A balloon heatmap drawn by R (the script of the Heatmap_Final folder): one balloon per line and gene, the trait strip on top, the favourable counts per line on the right.
           expected_output/gene_balloon.png (12 genes x 24 lines, 26 markers kept),  gene_freq_matrix.csv,  qc_marker_filter_report.csv,  locus_counts.csv
