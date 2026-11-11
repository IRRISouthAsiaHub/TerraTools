HETEROSIS - predict hybrid performance from the parents with a neural network  (run time about 90 s with the demo settings)

INPUT   input/het_train_phenotype.csv (ID, Yield_gain) and het_train_genotype.csv: 150 training lines x 288 SNPs.   input/het_parent_genotype.csv: 10 candidate parents.
        input/het_test_genotype.csv: 10 further lines to predict.  REAL genotypes, SIMULATED trait (8 QTL, h2 0.9), a gain relative to a check (mean about 1.5).
        The trait is kept near 0-3 on purpose: the network is not given a rescaled trait (see BreedSight).

STEPS   Genomic Selection > Heterosis > the four files.   F1 hybrids to report 5, parents 5, folds 3, epochs 100, batch 16, dropout 0.2, L1 0.001, L2 0.001.  > Run Prediction.
        (Defaults are 5 folds; the demo uses 3 so that it takes about a minute and a half.)

EXPECT  Cross-validation: Pearson r 0.63, R2 0.28, RMSE 0.70 - the network ranks lines well.  45 hybrids (10 parents) are scored.
        Top F1 hybrid BSP2-121 x BSP2-80, predicted 1.78, mid-parent heterosis +18.6 %.  Mid-parent and best-parent heterosis are measured against the mean and the maximum of the
        TRAINING lines, so best-parent heterosis is negative (the best training line is better than any hybrid).  Top parents and test-line predictions are listed as well.
