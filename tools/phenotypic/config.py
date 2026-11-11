from dataclasses import dataclass, field
from typing import List

MODULES = ("descriptive", "environments", "blues", "anova", "herit", "correlation", "pca", "cluster", "stability")
DESIGNS = ("augmented", "rcbd", "alpha")
SCHEMES = ("default", "blue", "viridis", "plasma", "magma")


@dataclass
class PhenotypicConfig:
    """Data settings (shared by every analysis) plus the settings of the one analysis being run."""
    input_file: str
    design: str = "augmented"           # augmented | rcbd | alpha
    analysis_mode: str = "all"          # "all" or "specific"
    selected_environments: List[str] = field(default_factory=list)
    color_scheme: str = "default"       # default | blue | viridis | plasma | magma
    zero_as_na: bool = True             # a 0 in a trait column means "not measured"
    module: str = ""                    # the analysis to run (one of MODULES)
    settings: dict = field(default_factory=dict)   # that analysis's settings: traits / trait, k, method, ...
    outdir: str = ""                    # folder of the earlier runs on the same data settings ("" = start a new one)

    @staticmethod
    def from_dict(d: dict) -> "PhenotypicConfig":
        return PhenotypicConfig(
            input_file=d.get("input_file", ""),
            design=d.get("design") if d.get("design") in DESIGNS else "augmented",
            analysis_mode="specific" if d.get("analysis_mode") == "specific" else "all",
            selected_environments=[str(e) for e in (d.get("selected_environments") or [])],
            color_scheme=d.get("color_scheme") if d.get("color_scheme") in SCHEMES else "default",
            zero_as_na=bool(d.get("zero_as_na", True)),
            module=d.get("module", ""),
            settings=d.get("settings") or {},
            outdir=d.get("outdir", "") or "",
        )

    def params(self) -> dict:
        """What phenotypic_cli.R reads (params.json)."""
        return {
            "input_file": self.input_file.replace(chr(92), "/"),
            "design": self.design,
            "analysis_mode": self.analysis_mode,
            "selected_environments": self.selected_environments,
            "color_scheme": self.color_scheme,
            "zero_as_na": self.zero_as_na,
            "module": self.module,
            "settings": self.settings,
        }
