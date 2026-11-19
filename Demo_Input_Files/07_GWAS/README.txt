GWAS - genome-wide association with three models  (run time about 12 s)

INPUT   input/gwas_ready.csv   ID, Phenotype and 288 SNP columns (Chr_Position) of 176 RILs - the ready-made format.
        input/genotypes.vcf.gz + phenotypes.csv   the same data as VCF + trait table (second input route).  Genotypes REAL, trait SIMULATED (QTL on chr1 and chr5).

        input/gwas_ready_position_rows.csv   the same genotypes in the layout with a Position row and a Chromosome row under the header
        (marker names are plain m1, m2, ...; positions come from the rows - the layout Linkage Mapping writes).

STEPS   Trait Mapping > GWAS > Data & Models: ready-made GWAS file > gwas_ready.csv.  Models: OLS, MLM, FarmCPU.  > Run GWAS.

EXPECT  Manhattan plot with one panel per model plus the QQ plot.  Genome-wide significant SNPs: OLS 31 (in blocks, linkage disequilibrium),
        MLM 4, FarmCPU 2 - one on chr1 and one on chr5, the two planted regions.  Inflation (lambda): OLS 3.27, MLM 0.58, FarmCPU 1.08.
