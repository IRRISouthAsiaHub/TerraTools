from dataclasses import dataclass, field


@dataclass
class GenomicPredictionConfig:
    train_pheno_file: str
    train_geno_file: str
    test_pheno_file: str
    test_geno_file: str

    outer_n_splits: int = 5
    heritability: float = 0.82

    epochs: int = 300
    batch_size: int = 32
    learning_rate: float = 0.001
    l2_reg: float = 0.001
    dropout_rate: float = 0.2
    rf_n_estimators: int = 300
    rf_max_depth: int = 0    # 0 means None (unlimited) — matches the Gradio slider convention
    alpha: float = 0.4       # DNN weight; RF weight = 1 - alpha

    use_pca: bool = True
    feature_selection: bool = True
    use_raw_genotypes: bool = False
    n_components: int = 50
    rfe_n_features: int = 100

    # optional models to compare BreedSight with ('Lasso', 'RRBLUP_1', 'GBLUP_1'); none = BreedSight only
    baselines: list = field(default_factory=list)

    @staticmethod
    def from_dict(d: dict) -> "GenomicPredictionConfig":
        return GenomicPredictionConfig(
            train_pheno_file=d.get("train_pheno_file", ""),
            train_geno_file=d.get("train_geno_file", ""),
            test_pheno_file=d.get("test_pheno_file", ""),
            test_geno_file=d.get("test_geno_file", ""),
            outer_n_splits=int(d.get("outer_n_splits", 5)),
            heritability=float(d.get("heritability", 0.82)),
            epochs=int(d.get("epochs", 300)),
            batch_size=int(d.get("batch_size", 32)),
            learning_rate=float(d.get("learning_rate", 0.001)),
            l2_reg=float(d.get("l2_reg", 0.001)),
            dropout_rate=float(d.get("dropout_rate", 0.2)),
            rf_n_estimators=int(d.get("rf_n_estimators", 300)),
            rf_max_depth=int(d.get("rf_max_depth", 0)),
            alpha=float(d.get("alpha", 0.4)),
            use_pca=bool(d.get("use_pca", True)),
            feature_selection=bool(d.get("feature_selection", True)),
            use_raw_genotypes=bool(d.get("use_raw_genotypes", False)),
            n_components=int(d.get("n_components", 50)),
            rfe_n_features=int(d.get("rfe_n_features", 100)),
            baselines=[str(b) for b in (d.get("baselines") or [])],
        )
