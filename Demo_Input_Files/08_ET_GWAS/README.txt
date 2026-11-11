ET-GWAS - extreme-trait GWAS, good for rare alleles  (run time about 10 s)

INPUT   input/etgwas_panel.csv   300 unrelated lines x 480 SNPs (Chr_Position) + Phenotype: a small diversity panel with weak sub-population structure. SIMULATED.
        answer_key_planted_loci.csv: a common SNP on chr3 (MAF 0.37) and a RARE SNP on chr8 (MAF 0.09) control the trait.

STEPS   Trait Mapping > ET-GWAS > Data: ready-made file > etgwas_panel.csv.   Bulks 20 %, MAF 0.05, FDR 0.05, merge 100 kb, genomic control ON.  > Run ET-GWAS.

EXPECT  Bulks of 60 lines (low / random / high), 437 SNPs tested, inflation factor lambda = 2.0 (corrected by genomic control).
        Exactly two MTAs at FDR < 0.05, both Bonferroni significant:  Chr8_24014817 (rare allele, MAF 0.087, PVE 14.9 %, q = 6e-07)
        and Chr3_9825966 (MAF 0.37, PVE 24.5 %, q = 8e-07) - the two planted loci.  Both alleles are enriched in the high bulk.
