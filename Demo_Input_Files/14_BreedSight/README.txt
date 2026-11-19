BREEDSIGHT - RICE DIVERSITY PANEL, SEED WIDTH  (real public data; run time about 1-2 minutes with the settings below)

DATA    377 of the 413 rice accessions of the Rice Diversity Panel 1 with a seed-width record: 317 training lines (input/training.csv)
        and 60 testing lines (input/testing.csv; their phenotypes are included, so BreedSight also reports an independent test result).
        The lines were split at random (seed 7). No line was chosen or left out because of how well it is predicted.

SOURCE  Rice Diversity Panel 1 (RDP1), 44K SNP set, Cornell University / USDA-ARS (McCouch lab), http://www.ricediversity.org/data/index.cfm
        Paper:      Zhao K, Tung C-W, Eizenga GC, et al. (2011) Genome-wide association mapping reveals a rich genetic architecture of complex
                    traits in Oryza sativa. Nature Communications 2:467. https://doi.org/10.1038/ncomms1467
        Genotypes:  http://www.ricediversity.org/data/sets/44kgwas/RiceDiversity.44K.MSU6.Genotypes_PLINK.zip   (sativas413.ped / .map: 413 accessions, 36,901 SNPs)
        Phenotypes: http://www.ricediversity.org/data/sets/44kgwas/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt   (column "Seed width", mm; 377 accessions have a value)
        The sample_id is the HybID of those files (RiceDiversity.44K.germplasm.csv gives the accession names).

WHAT WAS DONE TO THE RAW FILES
        1. SNPs with more than 10% missing calls or a minor allele frequency below 0.05 were removed: 29,564 of 36,901 SNPs kept.
        2. Genotypes were coded 0/1/2 (copies of the minor allele); the few missing calls were set to the rounded column mean.
        3. The additive matrix is VanRaden (2008): G = (M - 2p)(M - 2p)' / 2 sum p(1-p), with the allele frequencies p of the 317 TRAINING lines only
           (the same matrix as Additive Dominance, split mode). training_additive.csv is training x training (columns 1..317 = the training lines in the order
           of training.csv); test_additive.csv is testing x training. The testing lines therefore never influence the training matrix.
        4. Seed width was turned into a z-score with the mean and SD of the TRAINING lines (mean 3.118 mm, SD 0.393 mm; range 2.20 to 3.99 mm),
           because BreedSight does not rescale the trait. To get mm back: mm = 3.118 + 0.393 x phenotypes.

STEPS   Genomic Selection > BreedSight GUI > the four files (input/training.csv, input/training_additive.csv, input/testing.csv, input/test_additive.csv).
        Cross-validation folds 2, heritability 0.8, epochs 50; every other setting at its default (batch 32, learning rate 0.001, DNN weight 0.4, L2 0.001,
        dropout 0.2, RF trees 300, RF depth 0 = None, PCA on with 50 components, feature selection on (100), "Use raw genotypes" OFF); compare with Lasso,
        RBLUP and GBLUP.  > Run.  (The default of 5 folds and 300 epochs takes longer.)

EXPECT  Validation R2 / correlation of the 2-fold cross-validation (expected_output holds the run), and the same on the 60 independent testing lines:
                          2-fold cross-validation        60 testing lines
        BreedSight           0.57 / 0.76                    0.48 / 0.69   (95% CI of R2 0.28 - 0.63, of r 0.57 - 0.80)
        GBLUP                0.41 / 0.65                    0.28 / 0.55
        RBLUP                0.38 / 0.63                    0.40 / 0.64
        Lasso                0.38 / 0.62                    0.33 / 0.58
        BreedSight is best. The neural network is random: results vary by about +-0.03 between runs. The test R2 is lower than the cross-validation R2
        and, with only 60 lines, uncertain (see the interval); the cross-validation figure rests on 317 lines.

HONEST NOTES
        - Ten rice traits were run through BreedSight with the same settings and the same kind of random split (seed length, plant height, flowering time,
          amylose content, alkali spreading value, protein content, seed length/width ratio, seed surface area, brown-rice seed width and seed width).
          Seed width was kept because it gave the highest cross-validation R2 (above 0.5); brown-rice seed width (0.55) was similar. Because the trait was
          picked after seeing these results, the figures above are somewhat optimistic for a new trait. Seed dimensions are controlled by a few large-effect
          genes (GS3, GW5, GW2), which is why they are easier to predict than plant height or protein content.
        - A panel of 413 diverse accessions with strong structure (indica, japonica, aus, ...) is easier to predict than a narrow breeding population.
        - Please cite Zhao et al. (2011) when you use these files.
