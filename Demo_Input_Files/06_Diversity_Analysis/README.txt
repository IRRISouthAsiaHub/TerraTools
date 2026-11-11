DIVERSITY ANALYSIS - diversity indices of varieties from haplotype blocks  (run time about 10 s)

INPUT   input/diversity_haplotype_blocks.csv   120 samples (5 varieties x 24) x 150 haplotype blocks coded 0/1/2, 3 % missing (-1). SIMULATED:
        Variety_1 is a narrow variety, Variety_5 a diverse pool.

STEPS   Trait Mapping > Diversity Analysis > Browse: diversity_haplotype_blocks.csv.   Clusters: 5, top blocks 25, normalised radar on.  > Run Full Analysis.

EXPECT  Mean Shannon diversity per variety rises steadily: Variety_1 0.36, Variety_2 0.77, Variety_3 0.98, Variety_4 1.04, Variety_5 1.09
        (Effective alleles 1.2, 1.8, 2.4, 2.7, 2.9).  The PCA and the 5 hierarchical clusters reproduce the 5 varieties exactly (24 samples each).
        The Variety Comparison tab shows the gradient; the neighbour-joining tree separates five groups.  results.csv has every metric for every sample.
