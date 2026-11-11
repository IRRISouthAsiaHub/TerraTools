GENE PRIORITIZATION - rank candidate genes with public databases  (needs the internet; about 3 minutes the first time)

INPUT   input/candidate_genes.txt   8 rice genes (RAP-DB IDs): six well-known yield genes - GS3 (Os03g0407400), GW2 (Os02g0244100), DEP1 (Os09g0441900),
        IPA1/WFP (Os08g0509600), GW8 (Os08g0531600), Ghd7 (Os07g0261200) - and two controls with no yield role (Os01g0545100, Os05g0228400).

STEPS   Trait Mapping > Gene Prioritization > Species: rice.  Trait: grain yield.  Gene list: candidate_genes.txt. > Run.

EXPECT  The six known yield genes rank on top with the evidence level "Direct association" (scores about 49-62 of 100) and the two controls at 0 "No trait evidence".
        Evidence comes from Gramene trait annotations, Oryzabase curated gene names and trait classes (downloaded once and cached, nothing is added to the app), Plant Reactome pathways and PubMed co-mentions; each gene also gets a RAP-DB link; the score ranks the evidence, it is not a probability.
NOTE    The time depends on how fast Ensembl Plants, Gramene, Plant Reactome, NCBI and Oryzabase answer: about 3 minutes (the Oryzabase list, about 11 MB, is fetched once) for the first run, seconds for a repeat (the app caches answers).
        Scores can change a little over time. See input/slow_example_from_genes_extraction/ for the longer chain from Genes Extraction (41 genes, 10-14 minutes).
