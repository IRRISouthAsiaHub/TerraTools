ADDITIVE DOMINANCE - genomic relationship matrices (VanRaden / rrBLUP A.mat; dominance after Nishio & Satoh 2014 or Vitezica et al. 2013)  (run time about 10 s)

INPUT   input/snp_matrix.csv   176 rice RILs x 288 SNPs, coded 0/1/2. REAL genotypes.

STEPS   Genomic Selection > Additive Dominance > Browse: snp_matrix.csv.   Missing values: mean.   Matrices: Additive + Dominance.  > Compute Matrices.

EXPECT  176 samples x 288 markers, no missing data, 288 markers above the MAF filter, none monomorphic.
        Additive matrix: diagonal about 1.93 (RILs are nearly homozygous), mean off-diagonal about 0.  Heat maps show blocks of related lines along the diagonal.
        Downloads: Additive_Matrix.csv and Dominance_Matrix.csv.
        A warning says that only a few percent of the calls are heterozygous (RILs): the dominance matrix then carries almost no information. Use dominance for F2, BC1, F3 or hybrid material.

NEXT    For a training/testing split use the Training + Testing mode instead: choose the Training and Testing files of the Clustering tool.  Allele frequencies, centring and imputation come from the training lines only (additive and dominance), and the
        additive matrix is VanRaden's X_c X_c' / 2 sum p(1-p) (diagonal about 1), written as training_additive.csv (training x training) and test_additive.csv (testing x training), with the training lines as the columns of both.
