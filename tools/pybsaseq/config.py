from dataclasses import dataclass


@dataclass
class PyBSASeqConfig:
    parents_file: str
    bulks_file: str
    bulksizes: str = "430,385"
    popstruct: str = "F2"                    # F2 | RIL | BC
    pvalues: str = "0.01,0.01"
    replication: int = 10000
    slidingwindow: str = "2000000,10000"
    smoothing: str = "51,3"
    chromosome_order: str = ""               # optional, comma-separated
    read_length: int = 100
    region: str = "-1"
    misc_params: str = "20,5,1"
    adjust_gap: bool = False
    parent_ref: str = "no"

    @staticmethod
    def from_dict(d: dict) -> "PyBSASeqConfig":
        return PyBSASeqConfig(
            parents_file=d.get("parents_file", ""),
            bulks_file=d.get("bulks_file", ""),
            bulksizes=d.get("bulksizes", "430,385"),
            popstruct=d.get("popstruct", "F2"),
            pvalues=d.get("pvalues", "0.01,0.01"),
            replication=int(d.get("replication", 10000)),
            slidingwindow=d.get("slidingwindow", "2000000,10000"),
            smoothing=d.get("smoothing", "51,3"),
            chromosome_order=d.get("chromosome_order", ""),
            read_length=int(d.get("read_length", 100)),
            region=d.get("region", "-1"),
            misc_params=d.get("misc_params", "20,5,1"),
            adjust_gap=bool(d.get("adjust_gap", False)),
            parent_ref=d.get("parent_ref", "no"),
        )
