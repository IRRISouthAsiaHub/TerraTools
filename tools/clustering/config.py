from dataclasses import dataclass


@dataclass
class ClusteringConfig:
    input_file: str
    n_clusters: int = 5
    test_size: float = 0.2  # 0.1 - 0.8, matches the Gradio slider range
    pheno_file: str = ""            # optional: sample IDs + trait; adds a `phenotypes` column to the training file
    keep_test_pheno: bool = False   # also put the phenotypes in the test file (to measure test accuracy)

    @staticmethod
    def from_dict(d: dict) -> "ClusteringConfig":
        return ClusteringConfig(
            input_file=d.get("input_file", ""),
            n_clusters=int(d.get("n_clusters", 5)),
            test_size=float(d.get("test_size", 0.2)),
            pheno_file=d.get("pheno_file", "") or "",
            keep_test_pheno=bool(d.get("keep_test_pheno", False)),
        )
