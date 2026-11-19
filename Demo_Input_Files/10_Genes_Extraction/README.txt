GENES EXTRACTION - genes inside the QTL of QTL Mapping  (run time 1-2 s)

The QTL table AND the marker position file are both required (the file gives the base-pair position of the QTL markers). If you have no such table, give a VCF instead.

INPUT   input/qtl_table_cM_markers_without_position.csv  the QTL table of QTL Mapping (positions in cM) whose markers are called M001..M006, so the names say nothing about the genome.
        input/marker_positions.csv  the genetic map with the bp of those markers (Marker, Chromosome, bp).
        input/marker_positions.vcf  the same markers as a VCF (CHROM, POS, ID), for the case WITHOUT a genetic map.
        input/annotation_IRGSP1_regions.gff3   real rice genes (IRGSP-1.0, Ensembl Plants) of the regions of interest.

STEPS   Trait Mapping > Genes Extraction > QTL results: the QTL table > Marker position file: marker_positions.csv > annotation: the GFF3 > overlap > Run.
        With no marker table: choose "Use a VCF" and select marker_positions.vcf.

EXPECT  Two QTL regions (chr1 20.51-20.60 Mb, chr5 7.23-7.41 Mb) and 22 genes (position, strand, biotype, description) and a gene ID list for Gene Prioritization. The confidence
        interval in cM is placed on the genome between the markers of each QTL. Genes Extraction takes the QTL table of QTL Mapping only: a GWAS table or a table of intervals in bp is refused.

