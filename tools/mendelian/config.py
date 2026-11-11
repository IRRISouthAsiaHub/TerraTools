from dataclasses import dataclass


@dataclass
class MendelianConfig:
    input_file: str
    pop: str = "RIL"                 # "RIL" or "F2" — matches app.R's radio choices
    pthr: float = 0.05               # alpha, applied to the corrected p-value
    mtc_method: str = "fdr"          # multiple-testing correction: none | bonferroni | fdr (Benjamini-Hochberg)
    missing_codes: str = "NA,.,-9,-1"   # comma-separated
    drop_unplaced: bool = True

    @staticmethod
    def from_dict(d: dict) -> "MendelianConfig":
        return MendelianConfig(
            input_file=d.get("input_file", ""),
            pop=d.get("pop", "RIL"),
            pthr=float(d.get("pthr", 0.05)),
            mtc_method=d.get("mtc_method", "fdr") if d.get("mtc_method") in ("none", "bonferroni", "fdr") else "fdr",
            missing_codes=d.get("missing_codes", "NA,.,-9,-1"),
            drop_unplaced=bool(d.get("drop_unplaced", True)),
        )
