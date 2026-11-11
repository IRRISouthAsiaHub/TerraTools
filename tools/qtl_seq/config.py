"""
QTL-seq (SNP-index) — configuration and default analysis parameters.

Every threshold the analysis uses is listed in DEFAULTS and can be changed in the "Advanced settings" of the page.
The values written to analysis_parameters.json are the ones actually used.
"""

from pathlib import Path

from tools.common.paths import results_dir

TOOL_DIR = Path(__file__).resolve().parent
RESULT_DIR = results_dir("qtl_seq", TOOL_DIR / "ui" / "results")

LEVELS = (90, 95, 99)                       # confidence levels of the simulated limits (percent)
POPULATIONS = ("F2", "RIL", "BC")           # F2, recombinant inbred lines, backcross
MAX_SIM_DEPTH = 200                         # depths above this are simulated as 200 (the limit no longer changes)

DEFAULTS = {
    # --- input
    "mode": "bam",                          # "bam" or "vcf"
    "output_dir": "",
    # --- design
    "high_name": "High bulk",
    "low_name": "Low bulk",
    "parent1_name": "Parent 1",
    "parent2_name": "Parent 2",
    "index_parent": "parent2",              # whose allele the SNP-index counts: "parent1" or "parent2"
    "population": "F2",
    "n_high": 20,
    "n_low": 20,
    "bc_recurrent": "non_index",            # backcross only: recurrent parent is the one WITHOUT the index allele
    # --- SNP filters
    "min_depth": 7,                         # reads covering the SNP in EACH bulk
    "max_depth": 0,                         # 0 = no upper limit
    "min_snp_index": 0.3,                   # at least one bulk must reach this SNP-index
    "apply_snp_index_filter": True,
    "min_parent_depth": 7,
    "min_parent_gq": 30,                    # Phred-scaled confidence of a homozygous parent call
    "min_qual": 0,                          # VCF QUAL (0 = not used)
    "require_pass": False,                  # VCF FILTER must be PASS or "."
    "min_mapq": 20,                         # BAM counting
    "min_baseq": 20,
    "min_alt_reads": 2,                     # no/one parent: reads of the other allele needed (both bulks together)
    "min_alt_fraction": 0.05,               # ... and its share of the pooled bulk reads
    "base_error": 0.01,                     # per-read error rate of the genotype caller
    # --- sliding window
    "window_size": 2_000_000,
    "window_step": 10_000,
    "min_snps_per_window": 10,
    # --- simulation
    "replicates": 10_000,
    "seed": 12345,
    "depth_mode": "auto",                   # "auto" = depth of every SNP, "fixed" = one depth for all
    "fixed_depth": 20,
    "ci_method": "snp_average",             # "snp_average" or "window"  (see stats.py)
    # --- QTL calling
    "qtl_min_level": 90,                    # lowest confidence level that makes a window significant
    "merge_gap": 100_000,                   # significant windows closer than this are one QTL
    "interval_basis": "top_level",          # "top_level", "window_extent" or "window_midpoints"
    # --- BAM engine
    "engine": "auto",                       # "auto", "samtools" or "builtin"
    "chromosomes": "",                      # optional comma separated list; empty = every chromosome
    # --- VCF samples (empty = first sample of the file)
    "high_sample": "",
    "low_sample": "",
    "parent1_sample": "",
    "parent2_sample": "",
}


def merged(cfg):
    out = dict(DEFAULTS)
    out.update({k: v for k, v in (cfg or {}).items() if v is not None})
    return out
