from dataclasses import dataclass


@dataclass
class MendelianConfig:
    input_file: str
    pop: str = "RIL"                 # "RIL" or "F2" — matches app.R's radio choices
    pthr: float = 0.05
    missing_codes: str = "NA,.,-9"   # comma-separated, matches app.R's textInput default
    drop_unplaced: bool = True

    @staticmethod
    def from_dict(d: dict) -> "MendelianConfig":
        return MendelianConfig(
            input_file=d.get("input_file", ""),
            pop=d.get("pop", "RIL"),
            pthr=float(d.get("pthr", 0.05)),
            missing_codes=d.get("missing_codes", "NA,.,-9"),
            drop_unplaced=bool(d.get("drop_unplaced", True)),
        )
