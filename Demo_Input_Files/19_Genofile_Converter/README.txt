GENOFILE CONVERTER - the start of the chain: genotyping-company files -> one coded 0/1/2 genotype table for the other analyses   (run time: about 1 s)

INPUT   The same 345 markers x 176 RILs (+ the two parents P1 and P2), delivered the three ways a genotyping company sends them. Any one of them gives the same result:
        input/genotypes_hapmap.hmp.txt   HapMap: rs#, alleles, chrom, pos, strand ... QCcode, then one column per sample (calls AA, AG, GG, NN)
        input/genotypes.vcf              VCF: #CHROM POS ID REF ALT ... FORMAT, GT of every sample (0/0, 0/1, 1/1, ./.)
        input/genotypes.vcf.gz           the same VCF, gzipped
        input/raw_allele_calls.csv       raw sheet: Marker, Chromosome, P1, P2, then one column per RIL, calls like C:G
        input/answer_key_faulty_markers.csv   the 11 deliberately faulty markers (heterozygous parent, same allele in both parents, missing parent) and why each is left out

STEPS   Data Analysis > Genofile Converter > Browse: genotypes.vcf.gz (the format is detected).
        Parents and samples: Parent 1 = P1, Parent 2 = P2 (found automatically).  Output settings: defaults (marker names as in the file, chromosome numbers).  > Convert Genotype File.
        Then download the coded CSV / Excel file, or continue with it in another tool with the buttons under "Your converted file".
        For the raw CSV: Browse raw_allele_calls.csv, check the column roles (Marker, Chromosome, P1, P2) and the separator (colon), Convert.

EXPECT  345 markers read, 334 coded, 11 left out (3 missing parent, 4 heterozygous parent, 4 same allele in both parents); 176 individuals; 0.14 % of the calls missing / unexpected;
        code counts 0: 26,304   1: 2,069   2: 30,326   -1: 85.  12 chromosomes.
        expected_output/genotypes_coded.csv    the coded table: Marker, Chromosome, then one column per individual (0 = P1 allele, 1 = heterozygous, 2 = P2 allele, -1 = missing / unexpected);
                                               identical to 02_Mendelian_Segregation/input/coded_genotypes_RIL.csv
        expected_output/genotypes_coded.xlsx   the same table on the first sheet, then Marker_Details, Excluded_Markers (the 11 markers and why) and Summary
        Marker names such as chr01_1044946 carry the physical position, which analyses such as Linkage Mapping use to check the marker order.
