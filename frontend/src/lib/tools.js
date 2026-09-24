import {
  Dna, Link2, Waypoints, ScanSearch, Activity, Network, Globe, Group, Scale,
  Eye, Crosshair, StarHalf, ChartColumn, MapPinned, LayoutGrid,
} from 'lucide-react'

export const CATEGORIES = [
  { id: 'trait', label: 'Trait Mapping', icon: MapPinned },
  { id: 'selection', label: 'Genomic Selection', icon: Dna },
  { id: 'analysis', label: 'Data Analysis', icon: ChartColumn },
]
export const ALL_TAB = { id: 'all', label: 'All Applications', icon: LayoutGrid }

// `id` is what ToolManager.open_tool() expects on the Python side.
export const TOOLS = [
  { id: 'mendelian', n: 1, cat: 'trait', name: 'Mendelian Segregation', icon: Dna,
    desc: 'Analyze genetic inheritance patterns and test for Mendelian segregation ratios in plant populations.' },
  { id: 'linkage', n: 2, cat: 'trait', name: 'Linkage Mapping', icon: Link2,
    desc: 'Construct genetic linkage maps using molecular marker data and recombination frequencies.' },
  { id: 'qtl', n: 3, cat: 'trait', name: 'QTL Mapping', icon: Waypoints,
    desc: 'Identify quantitative trait loci (QTL) that affect complex traits using genetic marker data.' },
  { id: 'gene_extractor', n: 4, cat: 'trait', name: 'Genes Extraction', icon: ScanSearch,
    desc: 'Extract and analyze gene sequences from genomic regions of interest identified during mapping.' },
  { id: 'pybsaseq', n: 5, cat: 'trait', name: 'PYBSA G-statistics', icon: Activity,
    desc: 'Perform bulk segregant analysis using G-statistics to detect genomic regions associated with traits.' },
  { id: 'qtl_seq', n: 6, cat: 'trait', name: 'QTL-Seq', icon: Dna,
    desc: 'Perform QTL analysis using sequencing-based bulk segregant approaches for rapid gene mapping.' },
  { id: 'diversity', n: 8, cat: 'trait', name: 'Diversity Analysis', icon: Network,
    desc: 'Assess genetic diversity and population structure across germplasm collections using marker data.' },
  { id: 'gwas', n: 9, cat: 'trait', name: 'GWAS', icon: Globe,
    desc: 'Conduct genome-wide association studies to identify marker-trait associations across diverse panels.' },

  { id: 'clustering', n: 10, cat: 'selection', name: 'Clustering & Genomic Selection', icon: Group,
    desc: 'Cluster germplasm and apply genomic selection models to predict breeding values across populations.' },
  { id: 'additive_dominance', n: 11, cat: 'selection', name: 'Additive Dominance', icon: Scale,
    desc: 'Decompose genetic variance into additive and dominance components to improve genomic prediction accuracy.' },
  { id: 'genomic_prediction', n: 12, cat: 'selection', name: 'BreedSight GUI', icon: Eye,
    desc: 'Visual interface for managing breeding programs, tracking selections, and summarizing trial results.' },
  { id: 'cross_selection', n: 13, cat: 'selection', name: 'Optimal Cross Selection (OCS)', icon: Crosshair,
    desc: 'Identify optimal crosses that maximize genetic gain while maintaining diversity using OCS algorithms.' },
  { id: 'heterosis', n: 14, cat: 'selection', name: 'Heterosis', icon: StarHalf,
    desc: 'Quantify and analyze heterosis effects in hybrid combinations to guide parental line selection.' },

  { id: 'phenotypic', n: 16, cat: 'analysis', name: 'Phenotypic Data Analysis', icon: ChartColumn,
    desc: 'Analyze and visualize phenotypic data with a range of statistical methods, mixed models, and publication-ready plotting options.' },
]
