QTL-SEQ - bulked segregant analysis from a VCF  (run time about 10 s with 1000 simulation replicates)

INPUT   input/qtlseq_bulks.vcf   one joint VCF with four samples (Parent1, Parent2, HighBulk, LowBulk): 3 chromosomes x 2500 SNPs. SIMULATED, a QTL planted on chr2 at 12 Mb.

STEPS   Trait Mapping > QTL-Seq > Input: VCF files.  Choose qtlseq_bulks.vcf for high bulk, low bulk and parental VCF (the same file three times).
        Plants per bulk: 20 and 20, population F2, index parent: parent 2.  Settings tab: simulation replicates 1000 (default 10000).
        Turn on "advanced settings" and type the sample names: HighBulk, LowBulk, Parent1, Parent2.  > Run QTL-seq.

EXPECT  7,500 informative SNPs retained; one significant QTL at the 99 % level: chr2, 8.77-15.22 Mb, peak at 11.85 Mb, peak delta SNP-index 0.59,
        enriched in the high bulk.  Chromosomes 1 and 3 stay flat inside the confidence limits.  (Planted at 12 Mb.)
