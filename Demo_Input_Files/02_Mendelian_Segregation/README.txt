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
