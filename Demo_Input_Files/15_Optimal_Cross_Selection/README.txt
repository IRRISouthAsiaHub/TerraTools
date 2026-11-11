OPTIMAL CROSS SELECTION - which crosses to make: the highest usefulness, balanced against diversity  (run time about 12 s)

INPUT   input/dosage_matrix.csv    the 60 lines x 288 SNP allele dosage matrix (0 / 1 / 2). REAL genotypes (rice, BSP2 lines).
        input/marker_effects.csv   Marker, Effect: additive effect of each SNP. SIMULATED: ridge regression of a simulated breeding value on the dosages.
        input/genetic_map.csv      Marker, Chromosome, cM. SIMULATED: 4 cM per Mb of the physical position (a stand-in for a real linkage map;
                                   a real map from Linkage Mapping is used the same way).
        input/markers.vcf          the same SNPs as a VCF (their real physical positions), for the case WITHOUT a genetic map.
        There is no relationship matrix to prepare: it is calculated from the dosage matrix (VanRaden).

STEPS   Genomic Selection > Optimal Cross Selection > the three files > Preview & Validate (60 genotypes x 288 markers on 12 chromosomes;
        618 heterozygous calls are read as missing).
        Population: Doubled haploids, generation 1, selected proportion 0.05, candidate pool 60.  Crosses 20, max per parent 4, penalty 30.  > Run.
        No genetic map? Choose "No map: use a VCF", select input/markers.vcf and enter a recombination rate of 4 cM per Mb: this gives the same plan
        of crosses (the usefulness values differ only by rounding, because the demo map was made at that rate). For real data use the rate of your
        species; one rate for the whole genome is an approximation, so a genetic map from Linkage Mapping is more accurate.

EXPECT  For every one of the 1,770 crosses among the 60 lines the tool calculates the mean of the progeny, its standard deviation (from the marker
        effects and the recombination along the map) and the usefulness = mean + selection intensity x standard deviation. The plan of 20 crosses uses
        27 parents and maximises the mean usefulness (-3.73) at a group coancestry of 0.096 and a progeny inbreeding of 0; truncation selection
        gives -2.83 at a coancestry of 0.18 and an inbreeding of 0.079. Move the penalty toward 0 degrees for a higher usefulness and fewer parents.
        The values are on a relative scale (marker effects carry no intercept): only the differences between crosses matter.
        Outputs: the crossing plan, the contribution of every parent, the usefulness / diversity table and plot for every penalty degree, the
        usefulness of all crosses (ocs_all_crosses.csv) and the relationship matrix that was calculated (ocs_relationship_matrix.csv).
