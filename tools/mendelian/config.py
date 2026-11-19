from dataclasses import dataclass

# Population types with an explicit segregation model (see mendelian_engine.R). Anything else is rejected.
SUPPORTED_POPULATIONS = ("RIL", "DH", "BC1P1", "BC1P2", "F2", "F3", "F4")


def population_error(pop) -> str:
    """Message for a population type that has no segregation model, or "" when it is supported."""
    key = str(pop or "").strip().upper()
    if key in SUPPORTED_POPULATIONS:
        return ""
    why = ""
    if key == "F1":
        why = " F1 plants are all heterozygous, so there is no segregation to test."
    elif key == "BC1":
        why = " Say which parent the backcross is to: BC1P1 or BC1P2."
    return f"Unsupported population type '{pop}'. Supported types: {', '.join(SUPPORTED_POPULATIONS)}.{why}"


@dataclass
class MendelianConfig:
    input_file: str
    pop: str = "RIL"                 # one of SUPPORTED_POPULATIONS
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
