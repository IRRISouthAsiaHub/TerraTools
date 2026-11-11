LINKAGE MAPPING - coded genotypes to a genetic map  (run time: map about 25 s)

INPUT   input/2_coded_genotypes_RIL.csv   345 markers x 176 RILs coded 0/1/2 (Marker, Chromosome, then one column per RIL), 11 deliberately faulty markers.
        Raw allele calls from a genotyping company (HapMap, VCF, A:G sheets) are converted first by the Genofile Converter: see 19_Genofile_Converter, which writes exactly this file.

STEPS   Trait Mapping > Linkage Mapping > Upload tab > Browse: 2_coded_genotypes_RIL.csv (or "Open in Linkage Mapping" from the Genofile Converter).
        Upload tab: population RIL.  Parameters: defaults (missing 0.30, max recombination 0.45, best marker order).  > Run Analysis.

EXPECT  QC Report: 334 markers kept.  Map: 12 chromosomes, about 1360 cM in total.  The "QTL Mapping input (CSV)" button writes the genotypes with the map positions.
