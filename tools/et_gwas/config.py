from dataclasses import dataclass


@dataclass
class EtGwasConfig:
    """Settings of one Et-GWAS run (Gnanapragasam et al. 2024, doi:10.26508/lsa.202302352)."""

    # Where the genotypes come from (same three sources as the GWAS tool):
    #   "3k"   the bundled rice 3K panel + a phenotype file
    #   "vcf"  a VCF + a phenotype file
    #   "file" a ready-made file (ID, phenotype, then one column per SNP named Chr<chr>_<pos>, dosage 0/1/2)
    input_mode: str = "3k"
    input_file: str = ""
    phenotype_file: str = ""
    trait: str = ""
    vcf_file: str = ""

    bulk_pct: float = 20.0          # % of the lines in each extreme bulk (the paper compared 10, 15 and 20; 20 was best)
    maf: float = 0.05               # minimum minor-allele frequency over the panel
    min_alleles: int = 0            # minimum called alleles per bulk (0 = the bulk size, as in the paper)
    use_gc: bool = True             # divide the statistics by the inflation factor
    fdr: float = 0.05
    suggestive: float = 5e-5
    call_by: str = "fdr"            # "fdr" or "suggestive": what counts as a marker-trait association
    merge_kb: float = 100.0         # MTAs closer than this form one region (the paper bins at 100 kb)
    seed: int = 1                   # random bulk

    @classmethod
    def from_dict(cls, d: dict) -> "EtGwasConfig":
        valid = set(cls.__dataclass_fields__)
        return cls(**{k: v for k, v in d.items() if k in valid})
