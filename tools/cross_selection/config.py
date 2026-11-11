from dataclasses import dataclass


@dataclass
class OCSConfig:
    geno_file: str = ""             # allele dosage matrix (genotypes x markers, 0 / 1 / 2)
    effects_file: str = ""          # Marker, Effect
    map_file: str = ""              # genetic map (Marker, Chromosome, cM) or a VCF
    map_kind: str = "map"           # "map" | "vcf"
    cm_per_mb: float = 0.0          # recombination rate, used only with a VCF
    cross_type: str = "DH"          # DH | RIL
    generation: int = 1
    prop_sel: float = 0.05          # proportion of the progeny that is selected (selection intensity)
    pool_size: int = 60             # the best candidates whose crosses are evaluated
    n_crosses: int = 20             # crosses in the plan
    max_crosses_per_parent: int = 4
    penalty_degree: float = 30.0    # 0 = usefulness only, 90 = diversity only
    seed: int = 1

    @staticmethod
    def from_dict(d: dict) -> "OCSConfig":
        return OCSConfig(
            geno_file=d.get("geno_file", ""),
            effects_file=d.get("effects_file", ""),
            map_file=d.get("map_file", ""),
            map_kind=d.get("map_kind", "map"),
            cm_per_mb=float(d.get("cm_per_mb") or 0),
            cross_type=d.get("cross_type", "DH"),
            generation=int(d.get("generation", 1)),
            prop_sel=float(d.get("prop_sel", 0.05)),
            pool_size=int(d.get("pool_size", 60)),
            n_crosses=int(d.get("n_crosses", 20)),
            max_crosses_per_parent=int(d.get("max_crosses_per_parent", 4)),
            penalty_degree=float(d.get("penalty_degree", 30)),
            seed=int(d.get("seed", 1)),
        )
