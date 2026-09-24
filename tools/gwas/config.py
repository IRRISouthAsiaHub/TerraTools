from dataclasses import dataclass, field
from typing import List


@dataclass
class GWASConfig:
    """Configuration for one GWAS pipeline run. Built from whatever the
    frontend (HTML form or native Qt widgets) collects, then handed to
    GWASEngine.run()."""

    input_file: str

    models: List[str] = field(default_factory=lambda: ["OLS", "MLM"])

    maf: float = 0.05
    permutations: int = 100
    top_n: int = 15
    minimum_pve: float = 0.0

    # Advanced / per-model settings
    mlm_pcs: int = 3
    gemma_pcs: int = 3
    farmcpu_iterations: int = 10
    farmcpu_bin_size: int = 10          # scaled internally to bp (x10,000)
    blink_ld_r2: float = 0.7

    # Plot toggles
    show_bonferroni: bool = True
    show_suggestive: bool = True
    show_lod: bool = True
    show_pve: bool = False

    # mrMLM / FASTmrMLM — not exposed in the UI yet, sensible defaults
    mrmlm_screen_thresh: float = 0.01
    mrmlm_max_candidates: int = 20
    mrmlm_drop_thresh: float = 0.05
    mrmlm_n_pca: int = 3
    fastmrmlm_n_pca: int = 3

    def __repr__(self):
        return (
            f"GWASConfig(input_file={self.input_file!r}, models={self.models}, "
            f"maf={self.maf}, permutations={self.permutations}, top_n={self.top_n}, "
            f"minimum_pve={self.minimum_pve}, mlm_pcs={self.mlm_pcs}, "
            f"gemma_pcs={self.gemma_pcs}, farmcpu_iterations={self.farmcpu_iterations}, "
            f"farmcpu_bin_size={self.farmcpu_bin_size}, blink_ld_r2={self.blink_ld_r2}, "
            f"show_bonferroni={self.show_bonferroni}, show_suggestive={self.show_suggestive}, "
            f"show_lod={self.show_lod}, show_pve={self.show_pve})"
        )

    @classmethod
    def from_dict(cls, d: dict) -> "GWASConfig":
        """Build a GWASConfig from the JSON dict the HTML frontend sends
        over the Qt WebChannel bridge."""
        valid_fields = {f for f in cls.__dataclass_fields__}
        clean = {k: v for k, v in d.items() if k in valid_fields}
        return cls(**clean)
