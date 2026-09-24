from dataclasses import dataclass


@dataclass
class OCSConfig:
    criterion_file: str
    kmat_file: str
    threshold: float = 0.15         # kinship threshold, matches app.R slider (0.05-0.5)
    max_per_cluster: int = 5        # matches app.R slider (1-20)

    @staticmethod
    def from_dict(d: dict) -> "OCSConfig":
        return OCSConfig(
            criterion_file=d.get("criterion_file", ""),
            kmat_file=d.get("kmat_file", ""),
            threshold=float(d.get("threshold", 0.15)),
            max_per_cluster=int(d.get("max_per_cluster", 5)),
        )