RECOMBINATION ANALYSIS - crossovers, donor content and recombination along the chromosomes   (run time: about 5 s)

INPUT   input/mapping_population.csv   900 markers on 3 chromosomes (CHROM, POS, ID, REF, ALT, then one _GT column per sample), the recurrent parent IRRI-PARENT and 80 simulated
        backcross-derived lines with 2 % genotyping errors, 4 % missing and 1.5 % heterozygous calls. Every gene has two markers (chr1_12345678_LOC_Os01g00010).
        SIMULATED with known answers: the 80 lines carry 489 true crossovers (mean 6.1 per line). Made by docs/make_recombination_demo.py.

STEPS   Trait Mapping > Recombination Analysis > Browse: mapping_population.csv.  Recurrent parent: IRRI-PARENT (found automatically), no donor.
        Settings: smoothing "Original script", minimum block 2 %, bin 0.5 Mb.  > Run Analysis.

EXPECT  900 of 900 markers informative; 405 crossovers in 80 lines (mean 5.06 per line) = 83 % of the 489 planted; the rest are donor segments shorter than the 2 % block
        (about 6 markers) or crossovers hidden by missing calls. With a 5 % block only about 260 are found: the block size matters (see the Method Guide).
        Mean donor content 15.3 %; 37 intragenic crossovers (9.1 %).
        expected_output/crossover_breakpoints.csv      every crossover: line, chromosome, flanking markers and positions, size, direction, intragenic / intergenic, gene
        expected_output/per_RIL_summary.csv, per_chromosome_summary.csv, recombination_summary.csv, crossover_rate_by_bin.csv, genes_with_intragenic_crossovers.csv
        expected_output/Recombination_Analysis.xlsx    all tables in one workbook; Recombination_Analysis_all_files.zip holds everything
        expected_output/*.png                           crossovers per line, genome-wide rate, donor content, haplotype mosaic, circos plot, intragenic vs intergenic, top genes
