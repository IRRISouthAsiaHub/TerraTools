"""
Reproducible rice-like validation data for the GWAS models.

Population model
  Three subpopulations (labelled indica, japonica, aus) split from one ancestral population. The allele frequency of each marker in
  each subpopulation is drawn from a Balding-Nichols model with a chosen Fst (rice varietal groups are strongly differentiated, so
  the default is 0.25). Ancestral frequencies come from a skewed Beta distribution, so most markers are rare-to-intermediate and a
  few are common. A share of the samples are close relatives (full sibs of another sample of the same subpopulation). Genotype calls
  are missing at random with a per-marker rate, plus a few markers and samples with high missingness.
Phenotype model
  y = subpopulation mean shift + polygenic background + causal markers + noise
  The polygenic background is the sum of many small effects at HIDDEN markers that are simulated with the panel but never tested, so
  it creates the confounding that kinship and structure correction must absorb without being a true association in the tested set.
  Causal markers are tested markers (MAF >= causal_min_maf) with known effects, so detection, effect accuracy and false positives are
  all measured against a known truth. Markers are independent apart from the structure and the siblings; linkage disequilibrium
  within chromosomes is not simulated, which is stated in the validation report.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from tools.gwas.backend import genotype_io

GROUPS = ("indica", "japonica", "aus")


@dataclass
class Scenario:
    name: str
    n: int = 250
    m: int = 3000
    n_hidden: int = 600
    chromosomes: int = 12
    fst: float = 0.25
    group_sizes: tuple = (0.45, 0.35, 0.20)
    sib_fraction: float = 0.25
    missing_rate: float = 0.03
    maf_shape: tuple = (0.6, 1.2)            # Beta(a, b) of the ancestral frequency
    h2_polygenic: float = 0.3
    group_shift: float = 1.0                 # sd of the phenotype carried by subpopulation membership
    n_causal: int = 0
    causal_pve: float = 0.06                 # each causal marker, as a fraction of the phenotypic variance
    causal_min_maf: float = 0.10
    distribution: str = "normal"             # "normal" | "skewed" | "heavy_tailed"
    note: str = ""


@dataclass
class Simulated:
    data: genotype_io.GenotypeData
    labels: np.ndarray
    causal: np.ndarray                       # indices into data.snp_ids
    causal_beta: np.ndarray                  # per allele-count effect used to build y
    scenario: Scenario
    seed: int
    truth: dict = field(default_factory=dict)


def _frequencies(rng, sc, total):
    a, b = sc.maf_shape
    p0 = np.clip(rng.beta(a, b, total), 0.02, 0.98)
    f = sc.fst
    k = (1 - f) / f
    pops = np.vstack([rng.beta(p0 * k, (1 - p0) * k) for _ in GROUPS])          # Balding-Nichols
    return np.clip(pops, 0.002, 0.998)


def _group_of(rng, n, sizes):
    cut = np.round(np.cumsum(sizes) * n).astype(int)
    g = np.zeros(n, dtype=int)
    g[cut[0]:cut[1]] = 1
    g[cut[1]:] = 2
    return g


def _genotypes(rng, sc, freq, group):
    n = sc.n
    total = freq.shape[1]
    G = np.empty((n, total))
    for i in range(n):
        G[i] = rng.binomial(2, freq[group[i]])
    # full sibs: a share of samples are re-drawn from two parents of the same subpopulation
    n_sib = int(round(sc.sib_fraction * n))
    if n_sib:
        pick = rng.choice(n, size=n_sib, replace=False)
        for i in pick:
            same = np.flatnonzero(group == group[i])
            same = same[same != i]
            a, b = rng.choice(same, size=2, replace=False)
            # gametes from two sample genotypes treated as parents
            gam1 = rng.binomial(1, G[a] / 2.0)
            gam2 = rng.binomial(1, G[b] / 2.0)
            sib_of = rng.choice(same)
            G[i] = gam1 + gam2
            # a second child of the same parents so that sibling pairs exist
            if sib_of != i:
                G[sib_of] = rng.binomial(1, G[a] / 2.0) + rng.binomial(1, G[b] / 2.0)
    return G


def _noise(rng, n, kind):
    if kind == "normal":
        return rng.normal(size=n)
    if kind == "skewed":
        e = rng.lognormal(0.0, 0.9, n)
        return (e - e.mean()) / e.std()
    if kind == "heavy_tailed":
        e = rng.standard_t(4, n)
        return e / e.std()
    raise ValueError(f"unknown distribution {kind!r}")


def simulate(sc: Scenario, seed: int) -> Simulated:
    rng = np.random.default_rng(seed)
    total = sc.m + sc.n_hidden
    group = _group_of(rng, sc.n, sc.group_sizes)
    freq = _frequencies(rng, sc, total)
    G_all = _genotypes(rng, sc, freq, group)
    G, H = G_all[:, :sc.m], G_all[:, sc.m:]

    def standardise(M):
        s = M.std(0)
        s[s == 0] = 1
        return (M - M.mean(0)) / s

    poly = standardise(H) @ rng.normal(size=H.shape[1])
    poly = (poly - poly.mean()) / poly.std()
    shift = np.array([0.0, 1.0, -0.8])[group]
    shift = (shift - shift.mean()) / shift.std()

    resid_share = max(1 - sc.h2_polygenic, 0.05)
    other = sc.group_shift * shift + np.sqrt(sc.h2_polygenic) * poly + np.sqrt(resid_share) * _noise(rng, sc.n, sc.distribution)
    causal = np.empty(0, dtype=int)
    cbeta = np.empty(0)
    y = other
    if sc.n_causal:
        maf = np.minimum(G.mean(0) / 2, 1 - G.mean(0) / 2)
        pool = np.flatnonzero(maf >= sc.causal_min_maf)
        causal = rng.choice(pool, size=sc.n_causal, replace=False)
        total_var = other.var() / (1 - sc.n_causal * sc.causal_pve)               # so that each marker explains causal_pve of the final variance
        cbeta = np.array([np.sqrt(sc.causal_pve * total_var / G[:, j].var()) * rng.choice([-1, 1]) for j in causal])
        y = other + (G[:, causal] - G[:, causal].mean(0)) @ cbeta

    miss = rng.random(G.shape) < sc.missing_rate
    bad_snps = rng.choice(sc.m, size=max(sc.m // 100, 2), replace=False)       # a few markers with very high missingness
    miss[:, bad_snps] |= rng.random((sc.n, len(bad_snps))) < 0.6
    bad_samples = rng.choice(sc.n, size=max(sc.n // 60, 1), replace=False)
    miss[bad_samples] |= rng.random((len(bad_samples), sc.m)) < 0.5
    Gm = G.copy()
    Gm[miss] = np.nan

    per = sc.m // sc.chromosomes
    chrom = np.minimum(np.arange(sc.m) // per + 1, sc.chromosomes)
    pos = np.array([int(1000 + (k - np.flatnonzero(chrom == c)[0]) * 30000) for k, c in enumerate(chrom)])
    names = np.array([f"Chr{c}_{p}" for c, p in zip(chrom, pos)])
    enc = {"encoding": "dosage", "converted": False, "missing_calls": int(np.isnan(Gm).sum()),
           "message": "Simulated 0/1/2 alternate-allele dosage."}
    mm = genotype_io.MarkerMap(chrom=chrom, pos=pos, chrom_source="marker names", pos_source="marker names", valid=True)
    data = genotype_io.GenotypeData(ids=np.array([f"S{i + 1}" for i in range(sc.n)]), y=y, snp_ids=names, chrom=chrom, pos=pos,
                                    G=Gm, encoding=enc, marker_map=mm)
    return Simulated(data=data, labels=np.array(GROUPS, dtype=object)[group], causal=causal, causal_beta=cbeta, scenario=sc, seed=seed,
                     truth={"phenotype_sd": float(np.std(y)), "causal_names": [str(names[j]) for j in causal],
                            "bad_snps": int(len(bad_snps)), "bad_samples": int(len(bad_samples))})


SCENARIOS = {
    "null_structured": Scenario("null_structured", note="Strong structure and polygenic background, no causal marker: every hit is false."),
    "null_unstructured": Scenario("null_unstructured", fst=0.02, group_shift=0.0, sib_fraction=0.0, h2_polygenic=0.0,
                                  note="No structure, no relatedness, no background: the plain calibration check."),
    "causal_1": Scenario("causal_1", n_causal=1, causal_pve=0.10, note="One common causal marker, 10% of the variance."),
    "causal_3": Scenario("causal_3", n_causal=3, causal_pve=0.05, note="Three causal markers, 5% each."),
    "causal_skewed": Scenario("causal_skewed", n_causal=2, causal_pve=0.06, distribution="skewed", note="Skewed residuals."),
    "causal_heavy_tail": Scenario("causal_heavy_tail", n_causal=2, causal_pve=0.06, distribution="heavy_tailed",
                                  note="Heavy-tailed residuals."),
    "causal_weak_structure": Scenario("causal_weak_structure", n_causal=2, causal_pve=0.06, fst=0.05, group_shift=0.3,
                                      note="Weak structure, for comparison with the strongly structured rice-like case."),
}
