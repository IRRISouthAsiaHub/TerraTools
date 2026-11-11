from dataclasses import dataclass, field
from typing import List


@dataclass
class GWASConfig:
    """Configuration for one GWAS pipeline run. Built from whatever the
    frontend (HTML form or native Qt widgets) collects, then handed to
    GWASEngine.run()."""

    # Where the genotypes come from:
    #   "3k"   bundled rice 3K panel + a phenotype file (IDs as in the panel)
    #   "vcf"  a VCF + a phenotype file
    #   "file" an already prepared GWAS-ready file (ID, phenotype, SNP columns)
    input_mode: str = "file"
    input_file: str = ""
    phenotype_file: str = ""
    trait: str = ""
    vcf_file: str = ""

    models: List[str] = field(default_factory=lambda: ["OLS", "MLM"])

    # Genotype encoding of a prepared file: "auto" (decided only when the values make it unambiguous), "dosage" (0/1/2) or
    # "centred" (-1/0/1). A file that holds only 0 and 1 is never guessed: the run stops and asks for the choice.
    genotype_encoding: str = "auto"
    missing_codes: str = ""             # extra values to read as missing, comma separated (-9, -99, -999 and NA are always missing)

    # Quality control (each stage is reported with the numbers it removed)
    maf: float = 0.0                    # genotypes are expected to be filtered already; only unusable markers are dropped
    max_snp_missing: float = 1.0
    max_sample_missing: float = 1.0
    drop_duplicate_positions: bool = False
    imputation: str = "mean"            # "mean" | "remove" | "none" (data already imputed externally)

    # Relationship matrix and population structure
    kinship_method: str = "vanraden"    # "vanraden" (all markers) | "ld_pruned"
    ld_prune_r2: float = 0.2
    loco: bool = False                  # leave-one-chromosome-out kinship for the mixed models
    n_pcs: int = 3
    pca_ld_prune: bool = False
    compare_pcs: str = ""

    # Significance
    significance_rule: str = "bonferroni"   # category rule: "bonferroni" | "bh" | "permutation"
    alpha: float = 0.05
    permutations: int = 1000            # replicates of the model-specific null simulation (0 = none)
    permutation_null: str = "residual_permutation"   # or "parametric_bootstrap"
    permutation_seed: int = 1
    top_n: int = 15                     # exploratory markers listed per model (never called significant)
    minimum_pve: float = 0.0

    # Covariates per model (0 = none): OLS with 0 is the uncorrected baseline; the mixed models rely on K alone unless asked
    ols_pcs: int = 0
    mlm_pcs: int = 0
    emmax_pcs: int = 0
    gemma_pcs: int = 0
    gemma_test: str = "wald"            # "wald" | "lrt"

    # Iterative models
    farmcpu_iterations: int = 10

    # External reference program (optional; skipped with an explanation when absent)
    gemma_binary: str = ""
    compare_gemma: bool = False

    # Plot toggles
    show_bonferroni: bool = True
    show_suggestive: bool = True
    show_lod: bool = True
    show_pve: bool = False

    # mrMLM / FASTmrMLM (run by the mrMLM R package): search radius in kb, critical LOD, structure principal components
    mrmlm_radius_kb: float = 20.0
    mrmlm_lod: float = 3.0
    mrmlm_n_pca: int = 3

    def __repr__(self):
        return (
            f"GWASConfig(input_file={self.input_file!r}, models={self.models}, encoding={self.genotype_encoding}, "
            f"maf={self.maf}, imputation={self.imputation}, kinship={self.kinship_method}, loco={self.loco}, "
            f"permutations={self.permutations}, rule={self.significance_rule}, top_n={self.top_n})"
        )

    @classmethod
    def from_dict(cls, d: dict) -> "GWASConfig":
        """Build a GWASConfig from the JSON dict the HTML frontend sends
        over the Qt WebChannel bridge."""
        d = dict(d)
        valid_fields = {f for f in cls.__dataclass_fields__}
        clean = {k: v for k, v in d.items() if k in valid_fields}
        return cls(**clean)

    def missing_code_list(self):
        out = []
        for tok in str(self.missing_codes or "").replace(";", ",").split(","):
            tok = tok.strip()
            if tok:
                try:
                    out.append(float(tok))
                except ValueError:
                    pass
        return out

    def pc_compare_list(self):
        out = []
        for tok in str(self.compare_pcs or "").replace(";", ",").split(","):
            tok = tok.strip()
            if tok.isdigit():
                out.append(int(tok))
        return sorted(set(out))
