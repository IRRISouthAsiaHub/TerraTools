MENDELIAN SEGREGATION - find markers that do not follow the expected ratio  (run time about 10 s)

INPUT   input/coded_genotypes_RIL.csv   334 markers x 176 rice RILs (BSP2 population), coded 0/1/2: Marker, Chromosome, then one column per line. REAL genotypes.

STEPS   Trait Mapping > Mendelian Segregation > Browse: coded_genotypes_RIL.csv
        Population type: RIL (classes 0 and 2 in 1:1, exact binomial test).  Alpha 0.05.  Correction: FDR (Benjamini-Hochberg, the default).  Missing codes: NA, ., -9, -1.  > Run Segregation Test.

EXPECT  334 markers tested, none inconclusive (every marker has far more than 6 calls). The 2069 heterozygous calls of the RILs are not part of the test and are listed in the summary.
        FDR (default):   37 flagged distorted, 297 pass the quality check (results_filtered.csv, markers_LMready.xlsx).
                         Distorted per chromosome: chr1 16, chr6 10, chr3 7, chr8 3, chr7 1.
        Bonferroni:       8 distorted, 326 pass (chr1 5, chr3 2, chr6 1).
        None (raw p):   95 distorted, 239 pass (chr1 25, chr6 23, chr3 14, chr12 12, chr7 10, chr8 5, chr10 5, chr5 1) - uncorrected testing flags far more markers by chance.
        Three plots: manhattan.png (segregation plot): -log10(p) along the markers, one colour per chromosome, dashed line = cut-off of the correction (FDR: raw p 0.0055, markers above it are distorted); heterozygosity.png: % heterozygous calls per marker, same colours (3-9 % here, as expected in RILs); rateplot.png (distortion summary): % distorted markers per chromosome.
        markers_LMready.xlsx: the passing markers (Marker, Chromosome, genotypes), ready to load in Linkage Mapping.

OTHER POPULATIONS  Choose BC1P1 / BC1P2 (1:1), F2 (1:2:1), F3 (3:2:3) or F4 (7:2:7) for those generations: each has its own model, F1 is refused, and a marker with too few calls to reach alpha is reported as inconclusive (column call), neither distorted nor passed.

SECOND DEMO (a dense map)  input/coded_genotypes_RIL_4k_simulated.csv   4,499 markers x 150 RILs on 12 rice-sized chromosomes, coded 0/1/2 with -1 missing. SIMULATED (docs/make_mendelian_4k_demo.py, fixed seed): a proper advanced-generation RIL, about 4 % missing, about 1 % heterozygous, with three segregation-distortion regions placed on purpose.
        Same settings: Population type RIL, Alpha 0.05, FDR (BH), Missing codes NA, ., -9, -1.
EXPECT  4,499 markers tested, none inconclusive; 7,075 heterozygous calls are not part of the test.
        FDR (default):   156 flagged distorted, 4,343 pass. Distorted per chromosome: chr3 70 (about 21 to 29 Mb), chr9 50 (13 to 19 Mb), chr6 35 (10 to 15 Mb), chr2 1.
        Bonferroni:       27 distorted, all on chr3 (the strongest region).
        None (raw p):    571 distorted, spread over nine chromosomes: uncorrected testing flags far more markers by chance.
        manhattan.png: flat near 0 on most chromosomes with three clear peaks, on chr3, chr6 and chr9, above the dashed cut-off (FDR: raw p 0.0017); heterozygosity.png: 0 to 5 % per marker (mean 1.1 %); rateplot.png: distortion concentrated on chr3, chr6 and chr9. Outputs are in expected_output/RIL_4k_simulated/.
        Use this file to see what a normal RIL looks like. A population that is mostly one parent on every chromosome (for example a backcross-derived set run as RIL) shows almost every marker distorted; choose a type that fits the generation (see below) or treat the test only as a quality check.
