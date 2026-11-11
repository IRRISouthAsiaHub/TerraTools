ADDITIVE DOMINANCE - genomic relationship matrices (rrBLUP A.mat, sommer D.mat)  (run time about 10 s)

INPUT   input/snp_matrix.csv   176 rice RILs x 288 SNPs, coded 0/1/2. REAL genotypes.

STEPS   Genomic Selection > Additive Dominance > Browse: snp_matrix.csv.   Missing values: mean.   Matrices: Additive + Dominance.  > Compute Matrices.

EXPECT  176 samples x 288 markers, no missing data, 288 markers above the MAF filter, none monomorphic.
        Additive matrix: diagonal about 1.93 (RILs are nearly homozygous), mean off-diagonal about 0.  Heat maps show blocks of related lines along the diagonal.
        Downloads: Additive_Matrix.csv and Dominance_Matrix.csv.

NEXT    For a training/testing split use the Training + Testing mode instead: choose the Training and Testing files of the Clustering tool.  Markers are centred on the training lines and the
        additive matrix is X Xt / number of markers, written as training_additive.csv (training x training) and test_additive.csv (testing x training), with the training lines as the columns of both.
