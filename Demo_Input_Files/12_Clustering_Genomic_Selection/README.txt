CLUSTERING & GENOMIC SELECTION - PCA + k-means and a training/testing split  (run time about 5 s)

INPUT   input/snp_matrix.csv   176 rice RILs x 288 SNPs (Sample, then SNPs coded 0/1/2). REAL genotypes.

STEPS   Genomic Selection > Clustering & Genomic Selection > Browse: snp_matrix.csv.   Clusters 5, test size 0.20.  > Run Clustering & Split.

EXPECT  176 samples, 288 markers, 5 clusters (sizes 54 / 42 / 37 / 34 / 9), a stratified split: 140 training and 36 testing samples.
        PC1 explains 10.1 % and PC2 5.2 % of the variance.  A second figure (K selection) shows the silhouette width and the elbow for K = 2 to 10; "Principal components used" (default all) can be set to the first 5 to 20.  Download the training and testing sets: they are the input of the next steps.

NEXT    With a phenotype file (sample_id, phenotypes) chosen as well, the training file also gets the phenotypes column.  The "Relationship matrices" button then opens
        Additive Dominance with the training and testing files already filled in, to compute the matrices of those two files.
