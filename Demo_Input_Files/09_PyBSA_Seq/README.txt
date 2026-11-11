PyBSA-Seq - bulked segregant analysis with the G statistic  (run time about 50 s with 500 simulation replicates)

INPUT   input/pybsa_parents.tsv  and  input/pybsa_bulks.tsv   GATK4 VariantsToTable tables (CHROM, POS, REF, ALT, QUAL and <sample>.GT/.AD/.GQ):
        parents Parent1 / Parent2 (homozygous, different) and bulks HighBulk / LowBulk; 3 chromosomes x 1800 variants. SIMULATED, a QTL planted on chr2.

STEPS   Trait Mapping > PYBSA G-statistics > Parents file: pybsa_parents.tsv, Bulks file: pybsa_bulks.tsv.
        Bulk sizes 30,30   Population F2   p-values 0.01,0.01   Replication 500 (default 10000)   Sliding window 2000000,10000   Smoothing 51,3
        Read length 100   Region -1 (all)   Misc 20,5,1   Parent reference: no.  > Run Analysis.

EXPECT  An HTML report and figures.  In the four-panel figure the G statistic and the delta allele frequency leave their simulated 99 % thresholds only on chromosome 2:
        one clear peak at about 10 Mb (the largest significant region on chr2 spans 4.66-15.06 Mb, 1041 sliding windows).
        Chromosomes 1 and 3 stay inside the thresholds.  (The sSV/totalSV panel is sensitive and marks a few weak regions elsewhere; use G and delta AF.)
