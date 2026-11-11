from dataclasses import dataclass


@dataclass
class DiversityConfig:
    input_file: str
    n_clusters: int = 5      # matches Gradio slider: min 2, max 20, default 5
    top_blocks: int = 25     # matches Gradio slider: min 5, max 50, step 5, default 25
    radar_normalized: bool = True   # radar chart: 0-1 scale (True) or actual metric values (False)

    @staticmethod
    def from_dict(d: dict) -> "DiversityConfig":
        return DiversityConfig(
            input_file=d.get("input_file", ""),
            n_clusters=int(d.get("n_clusters", 5)),
            top_blocks=int(d.get("top_blocks", 25)),
            radar_normalized=bool(d.get("radar_normalized", True)),
        )
