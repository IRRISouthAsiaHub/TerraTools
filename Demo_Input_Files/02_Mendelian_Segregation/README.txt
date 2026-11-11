MENDELIAN SEGREGATION - find markers that do not follow the expected ratio  (run time about 10 s)

INPUT   input/coded_genotypes_RIL.csv   334 markers x 176 rice RILs (BSP2 population), coded 0/1/2: Marker, Chromosome, then one column per line. REAL genotypes.

STEPS   Trait Mapping > Mendelian Segregation > Browse: coded_genotypes_RIL.csv
        Population type: RIL (1:1).  Alpha 0.05.  Correction: FDR (Benjamini-Hochberg, the default).  Missing codes: NA, ., -9, -1.  > Run Segregation Test.

EXPECT  334 markers tested, none with a low-power (expected count < 5) warning.
        FDR (default):   39 flagged distorted, 295 pass the quality check (results_filtered.csv, markers_LMready.xlsx).
                         Distorted per chromosome: chr1 17, chr6 11, chr3 7, chr8 3, chr7 1.
        Bonferroni:       8 distorted, 326 pass (chr1 5, chr3 2, chr6 1).
        None (raw p):   104 distorted, 230 pass (chr1 25, chr6 23, chr3 17, chr7 12, chr12 12, chr8 6, chr10 6, chr5 2, chr2 1) - uncorrected testing flags far more markers by chance.
        manhattan.png: -log10(p) along the markers, red = distorted after the correction; rateplot.png: share of distorted markers per chromosome.
        markers_LMready.xlsx: the passing markers (Marker, Chromosome, genotypes), ready to load in Linkage Mapping.
