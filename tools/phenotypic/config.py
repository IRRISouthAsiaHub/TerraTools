from dataclasses import dataclass, field
from typing import List


@dataclass
class PhenotypicConfig:
    input_file: str
    trait: str = ""
    include_interaction: bool = False
    analysis_mode: str = "all"          # "all" or "specific"
    selected_environments: List[str] = field(default_factory=list)
    n_clusters: int = 3
    clustering_method: str = "kmeans"   # "kmeans" or "hierarchical"
    color_scheme: str = "default"       # default | blue | viridis | plasma | magma

    @staticmethod
    def from_dict(d: dict) -> "PhenotypicConfig":
        return PhenotypicConfig(
            input_file=d.get("input_file", ""),
            trait=d.get("trait", ""),
            include_interaction=bool(d.get("include_interaction", False)),
            analysis_mode=d.get("analysis_mode", "all"),
            selected_environments=d.get("selected_environments", []) or [],
            n_clusters=int(d.get("n_clusters", 3)),
            clustering_method=d.get("clustering_method", "kmeans"),
            color_scheme=d.get("color_scheme", "default"),
        )
