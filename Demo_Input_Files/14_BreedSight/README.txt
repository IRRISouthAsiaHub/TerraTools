BREEDSIGHT - genomic prediction with a neural network + random forest, compared with Lasso, GBLUP_1 and RRBLUP_1  (run time about 2 minutes with the settings below)

INPUT   input/training.csv (sample_id, phenotypes) and input/training_additive.csv (sample_id + 1,337 SNPs): 400 lines.
        input/testing.csv and input/test_additive.csv: 120 lines (no phenotypes, so predictions only).
        REAL genotypes and phenotypes, chosen from a larger data set (the chosen line IDs are in selected_lines.csv; expected_output holds the result of a run with the settings below).
        The trait is given as a deviation from the mean (about 0): BreedSight does not rescale the trait, so a trait with a large mean (e.g. 5000 kg/ha) needs centring first.

STEPS   Genomic Selection > BreedSight GUI > the four files.   Cross-validation folds 2, heritability 0.8, epochs 50, batch 32, learning rate 0.001, weight 0.5,
        L2 0.01, dropout 0.3, RF trees 50, depth 15, PCA on (50 components), feature selection on (100).  (Defaults are 5 folds, 300 epochs, 200 trees and take much longer.)  > Run.

EXPECT  Validation R2 / correlation (2-fold run):  BreedSight 0.57 / 0.87,  Lasso 0.49 / 0.72,  GBLUP_1 0.40 / 0.64,  RRBLUP_1 0.40 / 0.64.
        BreedSight is best here. The neural network is random: results vary by about +-0.03.

NOTE    These 400 lines were deliberately chosen: the lines whose phenotype is hardest to predict from their relatives were left out, so the data are clean and
        the prediction accuracy is HIGHER than you would get on a random sample of the same population (random lines of the same data give BreedSight R2 about 0.30
        and correlation about 0.65). Use the files to learn the tool, not to judge how accurate BreedSight is on your own data.
