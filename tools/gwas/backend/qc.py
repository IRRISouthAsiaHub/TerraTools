"""
Genotype quality control for GWAS.

Stages, in this order (every stage reports how many samples and SNPs it removed):
  1. duplicate marker positions     of two markers at the same chromosome and position the one with the better call rate stays
  2. SNP call rate                  markers with more than `max_snp_missing` missing calls
  3. sample call rate               samples with more than `max_sample_missing` missing calls (on the markers that remain)
  4. monomorphic markers            markers without variation among the remaining samples
  5. minor allele frequency         markers with MAF below `min_maf`
The raw matrix is never changed: QC returns index sets, per-marker and per-sample tables and, only if the user asked for it,
a separate imputed copy. The analysis stops with an explanation when fewer than two markers (or ten samples) remain; QC is
never skipped to keep a run going.

Missing genotypes
  "mean"    each missing call is replaced by the marker's mean dosage over the observed calls (2 x allele frequency). This is a
            simple, explicit option. It is NOT haplotype-based imputation: it ignores linkage disequilibrium, shrinks the
            variance of markers with many missing calls and biases results when missingness is high or related to the trait or
            to subpopulation. Prefer externally imputed genotypes (e.g. a haplotype-based tool) when much is missing.
  "remove"  markers with any missing call are removed.
  "none"    the data are taken as complete (externally imputed); the run stops if a missing call is left.
Allele frequencies are always computed from the observed calls only, so imputation does not move them.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from tools.gwas.backend.genotype_io import InputError

MIN_SNPS = 2
MIN_SAMPLES = 10
IMPUTATIONS = ("mean", "remove", "none")


class QCError(InputError):
    """The data cannot be analysed after quality control (the message says why)."""


@dataclass
class QCSettings:
    min_maf: float = 0.05
    max_snp_missing: float = 0.2
    max_sample_missing: float = 0.2
    drop_duplicate_positions: bool = True
    imputation: str = "mean"


@dataclass
class QCResult:
    snp_idx: np.ndarray                 # indices (into the raw marker list) that pass
    sample_idx: np.ndarray              # indices (into the raw sample list) that pass
    G: np.ndarray                       # dosage of the retained samples and markers; imputed only if requested
    G_observed: np.ndarray              # the same without imputation (NaN = missing)
    maf: np.ndarray                     # MAF of the retained markers
    call_rate: np.ndarray               # call rate of the retained markers (before imputation)
    ledger: list = field(default_factory=list)
    snp_table: dict = field(default_factory=dict)
    sample_table: dict = field(default_factory=dict)
    imputed_calls: int = 0
    notes: list = field(default_factory=list)


def allele_freq(G):
    """Frequency of the counted allele per marker from the observed calls only (NaN where nothing is observed)."""
    with np.errstate(invalid="ignore"):
        return np.nanmean(G, axis=0) / 2.0


def maf_of(G):
    p = allele_freq(G)
    return np.minimum(p, 1.0 - p)


def call_rate_of(G, axis=0):
    return 1.0 - np.isnan(G).mean(axis=axis)


def impute_mean(G):
    """Copy of G with each missing call replaced by the marker mean (the only values changed are the missing ones)."""
    out = np.array(G, dtype=float, copy=True)
    miss = np.isnan(out)
    if miss.any():
        mean = np.nanmean(out, axis=0)
        mean = np.where(np.isnan(mean), 0.0, mean)
        out[miss] = np.broadcast_to(mean, out.shape)[miss]
    return out


def _stage(ledger, name, rule, n_samples, n_snps, removed_samples=0, removed_snps=0):
    ledger.append({"stage": name, "rule": rule, "samples_removed": int(removed_samples), "snps_removed": int(removed_snps),
                   "samples_left": int(n_samples), "snps_left": int(n_snps)})


def run_qc(G, snp_ids, chrom, pos, sample_ids, settings: QCSettings):
    G = np.asarray(G, dtype=float)
    if settings.imputation not in IMPUTATIONS:
        raise InputError(f"Unknown missing-genotype option '{settings.imputation}' (use mean, remove or none).")
    n0, m0 = G.shape
    ledger = []
    alive_s = np.ones(n0, dtype=bool)
    alive_m = np.ones(m0, dtype=bool)
    reason_m = np.array([""] * m0, dtype=object)
    reason_s = np.array([""] * n0, dtype=object)
    cr0 = call_rate_of(G)
    maf0 = maf_of(G)
    _stage(ledger, "Input", f"{n0:,} samples, {m0:,} SNPs read", n0, m0)

    # 1 duplicate positions
    if settings.drop_duplicate_positions:
        key = chrom.astype(np.int64) * (1 << 40) + pos.astype(np.int64)
        order = np.lexsort((np.arange(m0), -cr0, key))           # within a position: best call rate first, then file order
        first = np.ones(m0, dtype=bool)
        first[1:] = key[order][1:] != key[order][:-1]
        dup = np.ones(m0, dtype=bool)
        dup[order[first]] = False
        reason_m[dup] = "duplicate position"
        alive_m &= ~dup
        _stage(ledger, "Duplicate positions", "one SNP per chromosome and position (best call rate kept)", n0, alive_m.sum(),
               0, dup.sum())
    # 2 SNP call rate
    low = alive_m & (cr0 < 1.0 - settings.max_snp_missing)
    reason_m[low] = "SNP call rate"
    alive_m &= ~low
    _stage(ledger, "SNP call rate", f"missing calls <= {settings.max_snp_missing:.0%}", n0, alive_m.sum(), 0, low.sum())
    # 3 sample call rate (on the SNPs that remain)
    if alive_m.any():
        cr_s = call_rate_of(G[:, alive_m], axis=1)
    else:
        cr_s = np.zeros(n0)
    bad_s = cr_s < 1.0 - settings.max_sample_missing
    reason_s[bad_s] = "sample call rate"
    alive_s &= ~bad_s
    _stage(ledger, "Sample call rate", f"missing calls <= {settings.max_sample_missing:.0%}", alive_s.sum(), alive_m.sum(),
           bad_s.sum(), 0)
    # 4 monomorphic, 5 MAF (computed on the remaining samples)
    sub = G[np.ix_(alive_s, np.arange(m0))] if alive_s.any() else G[:0]
    with np.errstate(invalid="ignore"):
        maf_s = maf_of(sub) if len(sub) else np.full(m0, np.nan)
        var_s = np.nanvar(sub, axis=0) if len(sub) else np.zeros(m0)
    mono = alive_m & ~(np.nan_to_num(var_s) > 1e-12)
    reason_m[mono] = "monomorphic"
    alive_m &= ~mono
    _stage(ledger, "Monomorphic SNPs", "no variation among the remaining samples", alive_s.sum(), alive_m.sum(), 0, mono.sum())
    rare = alive_m & (np.nan_to_num(maf_s, nan=0.0) < settings.min_maf)
    reason_m[rare] = "MAF"
    alive_m &= ~rare
    _stage(ledger, "Minor allele frequency", f"MAF >= {settings.min_maf:g}", alive_s.sum(), alive_m.sum(), 0, rare.sum())

    # 6 missing calls that are left
    Gq = G[np.ix_(alive_s, alive_m)] if (alive_s.any() and alive_m.any()) else G[:0, :0]
    notes = []
    imputed = 0
    if Gq.size and settings.imputation == "remove":
        has_missing = np.isnan(Gq).any(axis=0)
        idx = np.flatnonzero(alive_m)
        reason_m[idx[has_missing]] = "missing calls (remove option)"
        alive_m[idx[has_missing]] = False
        Gq = Gq[:, ~has_missing]
        _stage(ledger, "Missing calls", "SNPs with any missing call removed", alive_s.sum(), alive_m.sum(), 0, has_missing.sum())
    n, m = alive_s.sum(), alive_m.sum()
    if m < MIN_SNPS or n < MIN_SAMPLES:
        raise QCError(_failure_message(ledger, n, m, settings))
    Gobs = Gq
    if np.isnan(Gq).any():
        imputed = int(np.isnan(Gq).sum())
        if settings.imputation == "none":
            raise QCError(f"{imputed:,} missing genotype calls are left after quality control and imputation is switched off. "
                          "Provide externally imputed genotypes, choose mean imputation, or remove the SNPs with missing calls.")
        Gq = impute_mean(Gq)
        notes.append(f"Mean imputation replaced {imputed:,} missing calls ({imputed / Gq.size:.2%} of the genotype matrix). "
                     "This is not haplotype-based imputation: it ignores linkage disequilibrium and shrinks the variance of "
                     "markers with many missing calls. Allele frequencies come from the observed calls.")
        _stage(ledger, "Imputation", "mean dosage of the marker replaces each missing call", n, m, 0, 0)
    snp_idx, sample_idx = np.flatnonzero(alive_m), np.flatnonzero(alive_s)
    status_m = np.where(alive_m, "kept", reason_m)
    status_s = np.where(alive_s, "kept", reason_s)
    snp_table = {"SNP": np.asarray(snp_ids).astype(str), "Chromosome": chrom, "Position": pos, "Call_rate": cr0, "MAF_raw": maf0,
                 "Status": status_m}
    sample_table = {"Sample": np.asarray(sample_ids).astype(str), "Call_rate": cr_s, "Status": status_s}
    return QCResult(snp_idx=snp_idx, sample_idx=sample_idx, G=Gq, G_observed=Gobs, maf=maf_of(Gobs), call_rate=call_rate_of(Gobs),
                    ledger=ledger, snp_table=snp_table, sample_table=sample_table, imputed_calls=imputed, notes=notes)


def _failure_message(ledger, n, m, settings):
    lines = []
    for st in ledger[1:]:
        if st["samples_removed"] or st["snps_removed"]:
            lines.append(f"  {st['stage']}: removed {st['samples_removed']:,} samples and {st['snps_removed']:,} SNPs ({st['rule']})")
    why = "\n".join(lines) or "  (nothing was removed: the input itself is too small)"
    if m < MIN_SNPS:
        head = f"Only {m} SNP(s) pass quality control; at least {MIN_SNPS} are needed, so the analysis cannot proceed."
    else:
        head = f"Only {n} sample(s) pass quality control; at least {MIN_SAMPLES} are needed, so the analysis cannot proceed."
    first = ledger[0]["rule"]
    return (f"{head}\nStart: {first}.\n{why}\n"
            f"Thresholds used: MAF >= {settings.min_maf:g}, SNP missing <= {settings.max_snp_missing:.0%}, sample missing <= "
            f"{settings.max_sample_missing:.0%}. Relax the threshold that removed most (for example the MAF filter) or check the input.")
