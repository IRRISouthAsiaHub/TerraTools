GENOTYPE VIEWER - graphical genotypes of every chromosome for each sample, drawn with GenoSee  (run time about 5 s)

INPUT   input/genotypes_12_RILs.csv   308 markers x 12 rice RILs (BSP2 population): chromosome, marker name, position (bp), then one column per line.
        Calls are A / B / H and N for missing.  REAL genotypes of the Linkage Mapping demo; the markers whose name carries no position were left out, the
        positions are the ones in the marker names, and the codes 0 / 2 / 1 / -1 became A / B / H / N (so A and B are the two parents, REF and ALT in the picture).
        input/genotypes_12_RILs.vcf   the same data as a VCF (GT field: 0/0, 1/1, 0/1, ./.), to try the VCF input.
        The column names do not matter: the viewer takes the first column as the chromosome, the second as the marker, the third as the position.

STEPS   Trait Mapping > Genotype Viewer > Browse: genotypes_12_RILs.csv.   Species: Oryza_sativa.   (the chromosome names chr01 ... chr12 match the species)
        1. Drawing mode Normal, 3-color, palette normal (bright), "Show marker names in the figure" off (the default), 300 dpi.  Legend names: REF = Parent 1, ALT = Parent 2, HET = Heterozygous.  Samples to draw: BSP2-79.  > Draw Graphical Genotypes.
        2. Drawing mode Compare (one image per chromosome with all 12 lines side by side), same colouring.
        3. Drawing mode Zoomed: Chromosome chr01, Start 0, End 15000000, palette Aqua.
        4. VCF: browse genotypes_12_RILs.vcf, leave "Name the markers CHROM_POS" ticked.  Same markers, same picture.
        For a figure: choose 600 dpi and tick "Also save a vector PDF".

EXPECT  1. 1 image: the 12 chromosomes of BSP2-79, 308 markers, 3 colours (yellow Parent 1, blue Parent 2, green Heterozygous, grey MISSING).   expected_output/normal_3-color_bright.png
        2. 12 images (chr01 ... chr12), 12 lines each.                                                                      expected_output/compare_chr01.png
        3. 1 image of the first 15 Mb of chr01 for the 12 lines.                                                            expected_output/zoomed_chr01_0-15Mb_Aqua.png
        4. 308 markers, 12 samples, as for the CSV.
        Palettes (GenoSee's own): normal = bright, Aqua, grays, reds, blues.  Chromosome lengths come from GenoSee's species table (Oryza_sativa).
        Cite: Hashimoto S. (2024), GenoSee, https://github.com/hashimotoshumpei/GenoSee
