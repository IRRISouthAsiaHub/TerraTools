PHENOTYPIC DATA ANALYSIS - field trial in three environments  (run time about 20 s)

INPUT   input/phenotypic_trial.csv   180 plots: Genotype, Environment, Rep, Block, Yield, Plant_Height, Days_to_Flowering.
        40 new entries (G001-G040) + 4 checks, augmented design with 5 blocks, 3 environments (Dry_2024, Wet_2024, Dry_2025). SIMULATED.

STEPS   Data Analysis > Phenotypic Data Analysis > Data & Design > Browse: phenotypic_trial.csv  (design: Augmented RCBD is picked from the Block column).
        Analyses tab: pick an analysis, keep all three traits ticked, press its Run button. Each analysis runs on its own (about 5-30 s each);
        run Descriptives, Heritability, BLUPs, ANOVA, PCA, Clustering (3 clusters, k-means) and Stability & GxE (trait Yield) to fill the tabs.

EXPECT  Descriptive statistics: mean yield Dry_2024 42.2, Dry_2025 47.3, Wet_2024 51.1 - the dry environments yield less.
        Heritability (pooled over environments), on a plot basis / on an entry-mean basis:
              Yield 0.49 / 0.74     Plant_Height 0.52 / 0.76     Days_to_Flowering 0.17 / 0.38
        BLUPs (Yield, pooled over environments): 40 new entries, H2 Cullis 0.73, mean reliability 0.71; best entry G011 (BLUP +4.77, SE 1.56, predicted yield 51.6).
              Downloads > BLUPs (.xlsx) gives the BLUP table and the variance/heritability table. The 4 checks are fixed effects and have no BLUP.
        ANOVA: genotypes (adjusted for blocks) are significant for every trait in the pooled analysis (Yield P < 0.001, Plant_Height P < 0.001, Days_to_Flowering P < 0.05).
        Design Check: 5 blocks, 4 repeated checks, 40 single entries per environment.  PCA, Clustering and Stability tabs fill with plots.
        expected_output/ holds every table (CSV), the main figures (PNG) and the two Excel workbooks (Phenotypic_Analysis.xlsx, ANOVA_tables.xlsx).
        (The values come from simulated data; yours should match to the printed digits.)
