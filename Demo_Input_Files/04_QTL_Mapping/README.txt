QTL MAPPING - SIM and CIM on three traits  (run time about 5 s per run)

INPUT   input/genotypes_QTL_format.csv + linkage_map.xlsx   output of Linkage Mapping.   input/phenotypes.csv   Yield, Plant_Height, No_QTL_Trait (SIMULATED).
        answer_key_planted_QTL.csv: Yield QTL on chr1 (+2.0) and chr5 (+1.2); Plant_Height on chr9 (+4.0) and chr5 (+2.5); No_QTL_Trait none.

STEPS   Trait Mapping > QTL Mapping > load the three files.
        Run 1: trait Yield, SIM, fixed LOD 3.     Run 2: trait Yield, CIM, permutation test (1000 permutations, 5 %).
        Run 3: Plant_Height, CIM, permutation.     Run 4: No_QTL_Trait, CIM, permutation (the negative control).

EXPECT  SIM at LOD 3 finds 4 QTL: the two planted ones and two "ghost" peaks on chr4 and chr7.
        CIM with the permutation threshold (LOD 3.54): only chr1 (LOD about 62) and chr5 (about 32) remain.  Plant_Height: chr9 and chr5.  No_QTL_Trait: 0 QTL.
