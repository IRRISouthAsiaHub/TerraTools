import {
  Dna, Link2, Waypoints, ScanSearch, Activity, Network, Globe, Group, Scale,
  Eye, Crosshair, StarHalf, ChartColumn, MapPinned, LayoutGrid, Layers, ListOrdered, Rows3, Fingerprint, FileCog, Shuffle, FlaskConical,
} from 'lucide-react'

export const CATEGORIES = [
  { id: 'analysis', label: 'Data Analysis', icon: ChartColumn },
  { id: 'trait', label: 'Trait Mapping', icon: MapPinned },
  { id: 'selection', label: 'Genomic Selection', icon: Dna },
]
export const ALL_TAB = { id: 'all', label: 'All Applications', icon: LayoutGrid }

// `id` is what ToolManager.open_tool() expects on the Python side.
export const TOOLS = [
  { id: 'genofile', n: 19, cat: 'analysis', name: 'Genofile Converter', icon: FileCog,
    desc: 'Start here: convert the HapMap, VCF or raw allele files of a genotyping company into the 0 / 1 / 2 input files of Mendelian Segregation and Linkage Mapping.' },
  { id: 'phenotypic', n: 1, cat: 'analysis', name: 'Phenotypic Data Analysis', icon: ChartColumn,
    desc: 'Analyze and visualize phenotypic data with a range of statistical methods, mixed models, and publication-ready plotting options.' },

  { id: 'mendelian', n: 2, cat: 'trait', name: 'Mendelian Segregation', icon: Dna,
    desc: 'Analyze genetic inheritance patterns and test for Mendelian segregation ratios in plant populations.' },
  { id: 'recombination', n: 20, cat: 'trait', name: 'Recombination Analysis', icon: Shuffle,
    desc: 'Locate the crossovers of RILs or backcross-derived lines, measure the donor genome of every line, and map recombination hot and cold regions, including intragenic crossovers.' },
  { id: 'linkage', n: 3, cat: 'trait', name: 'Linkage Mapping', icon: Link2,
    desc: 'Construct genetic linkage maps using molecular marker data and recombination frequencies.' },
  { id: 'genotype_viewer', n: 17, cat: 'trait', name: 'Genotype Viewer', icon: Rows3,
    desc: 'Graphical genotype visualisation, used for parent recovery: draw every chromosome of each line coloured by the parent it came from, in bright, publication-ready figures (GenoSee).' },
  { id: 'qtl', n: 4, cat: 'trait', name: 'QTL Mapping', icon: Waypoints,
    desc: 'Identify quantitative trait loci (QTL) that affect complex traits using genetic marker data.' },
  { id: 'qtl_seq', n: 5, cat: 'trait', name: 'QTL-Seq', icon: Dna,
    desc: 'Map QTLs from bulked sequencing data (BAM or VCF) with the Δ SNP-index, simulated confidence limits and candidate genes.' },
  { id: 'insilico_bsa', n: 21, cat: 'trait', name: 'In silico BSA', icon: FlaskConical,
    desc: 'Bulked segregant analysis without sequencing pools: pick the extreme lines of an F2 or mapping population by trait value and compare their genotypes, with the Gw statistic that also finds over-dominant loci.' },
  { id: 'diversity', n: 6, cat: 'trait', name: 'Diversity Analysis', icon: Network,
    desc: 'Assess genetic diversity and population structure across germplasm collections using marker data.' },
  { id: 'gwas', n: 7, cat: 'trait', name: 'GWAS', icon: Globe,
    desc: 'Conduct genome-wide association studies to identify marker-trait associations across diverse panels.' },
  { id: 'et_gwas', n: 8, cat: 'trait', name: 'ET-GWAS (Extreme Trait)', icon: Layers,
    desc: 'Compare allele frequencies of low, random and high phenotype bulks to uncover associations, including rare variants.' },
  { id: 'pybsaseq', n: 9, cat: 'trait', name: 'PYBSA G-statistics', icon: Activity,
    desc: 'Perform bulk segregant analysis using G-statistics to detect genomic regions associated with traits.' },
  { id: 'gene_extractor', n: 10, cat: 'trait', name: 'Genes Extraction', icon: ScanSearch,
    desc: 'Extract the genes, their annotation and sequences from the QTL found by QTL Mapping.' },

  { id: 'gene_prioritization', n: 11, cat: 'trait', name: 'Gene Prioritization', icon: ListOrdered,
    desc: 'Rank candidate genes for a trait by combining Ensembl Plants, Gramene, Plant Reactome, PubMed and (rice) Oryzabase evidence; give it the gene list of Genes Extraction.' },

  { id: 'clustering', n: 12, cat: 'selection', name: 'Clustering & Genomic Selection', icon: Group,
    desc: 'Cluster lines from SNP genotypes (PCA + K-Means) and split them into stratified training and testing sets; add a phenotype file to get a phenotyped training set.' },
  { id: 'additive_dominance', n: 13, cat: 'selection', name: 'Additive Dominance', icon: Scale,
    desc: 'Compute additive (rrBLUP) and dominance (sommer) genomic relationship matrices, or the training and testing matrices of the two files made by the Clustering tool.' },
  { id: 'genomic_prediction', n: 14, cat: 'selection', name: 'BreedSight GUI', icon: Eye,
    desc: 'Genomic prediction with a neural network + random forest ensemble, benchmarked against Lasso, GBLUP and RRBLUP, from the training and testing files built by Clustering and Additive Dominance.' },
  { id: 'cross_selection', n: 15, cat: 'selection', name: 'Optimal Cross Selection (OCS)', icon: Crosshair,
    desc: 'Identify the crosses with the highest usefulness, predicted from marker data and a genetic map, while maintaining diversity.' },
  { id: 'heterosis', n: 16, cat: 'selection', name: 'Heterosis', icon: StarHalf,
    desc: 'Quantify and analyze heterosis effects in hybrid combinations to guide parental line selection.' },
  { id: 'favourable_allele', n: 18, cat: 'selection', name: 'Favourable Allele Analysis', icon: Fingerprint,
    desc: 'Upload your rice 4K haplopanel genotypes and see which favourable alleles of 166 known trait genes each line carries (rice 4K only, for now).' }
]
