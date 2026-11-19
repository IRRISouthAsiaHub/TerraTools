IN SILICO BSA - bulked segregant analysis from the genotype and trait of every line   (run time: about 5 s)

INPUT   input/genotypes_F2.csv   900 markers on 3 chromosomes (CHROM, POS, ID, then one column per line), 300 simulated F2 plants coded 0 / 1 / 2 / -1
        (0 = parent 1, 1 = heterozygous, 2 = parent 2; 1 % genotyping errors, 3 % missing calls).
        input/phenotypes.csv     Sample, PlantHeight (the first column holds the line names exactly as in the genotype file).
        SIMULATED with known answers (docs/make_insilico_demo.py): three QTL of about 8 % of the variance each:
        chr1 at 7.05 Mb additive-dominant (d = 0.5 a), chr2 at 9.97 Mb FULLY over-dominant (a = 0), chr3 at 13.94 Mb completely dominant (d = a).

STEPS   Trait Mapping > In silico BSA > Browse: genotypes_F2.csv, then the phenotype file phenotypes.csv (trait PlantHeight).
        Pools 20 % (60 + 60 lines), population type: F2 / IMF2, primary window 500 kb, merge regions closer than 1000 kb, 1000 permutations, level 0.05.  > Run Analysis.

EXPECT  3 loci, each found by all three window sizes (250 kb, 500 kb and 1 Mb); 3 significant regions at every size. Gw threshold of the 500 kb windows 5.92
        (250 kb: 6.30, 1 Mb: 5.46), G1 4.78, |dAF| 0.20. At the 500 kb size:
        R1 chr3 peak 13.7 Mb (adjusted p 0.001, LOD 10.0), R2 chr1 peak 6.7 Mb (LOD 7.8) - allele-frequency shifts, also found by G1 and the delta allele frequency;
        R3 chr2 peak 11.1 Mb (LOD 5.6) - the over-dominant locus: heterozygotes are 71 % of the high pool against 22 % of the low pool, with almost the same allele
        frequency in both pools. Only Gw finds it ("Also_significant_by_G1 = False"; see methods_compared.png).
        expected_output/loci_comparison.csv         one row per locus: span, peak, how many window sizes found it, length and p-value at each size
        expected_output/significant_regions.csv     one row per region and window size: left / peak / right marker, genotype classes at the peak, pattern, adjusted p-value
        expected_output/window_comparison.csv       the three window sizes side by side (threshold, windows tested, regions)
        expected_output/marker_statistics.csv       G1, G2, Gw and the window values of every marker
        expected_output/InSilico_BSA.xlsx           Summary, Settings, Window comparison, Loci, Regions, Thresholds, Windows, Pools, Markers; InSilico_BSA_all_files.zip holds everything
        expected_output/*.png                       gw_windows_compared, loci_compared, methods_compared, genotypes_at_peaks, pools, permutation_null

NOTE    The permutations use a fixed seed, so the thresholds repeat. Real data differ: a very polygenic trait gives many regions, and in a RIL population
        (no heterozygous lines) Gw equals G1. Method: Li et al. 2022, Plant Communications 3:100319 (dQTG-seq).
