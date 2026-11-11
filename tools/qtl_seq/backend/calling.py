"""Genotype call of an inbred parent from allele read counts (the same caller serves BAM and VCF input)."""

from __future__ import annotations

import numpy as np


def homozygous_call(counts, error=0.01):
    """Is a parent homozygous, and for which allele?

    counts: array (alleles, sites) of read counts. Two hypotheses are compared with a binomial likelihood and equal
    priors: homozygous for the most frequent allele (every other read is an error with probability `error`/3
    each) against heterozygous (every read is either of two alleles, 1/2 each). Returns
    (allele index, depth, gq) where gq is the Phred-scaled confidence of the homozygous call (capped at 99)."""
    c = np.asarray(counts, dtype=float)
    depth = c.sum(axis=0)
    top = c.argmax(axis=0)
    major = c.max(axis=0)
    minor = depth - major
    e = float(error)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        ll_hom = major * np.log(1 - e) + minor * np.log(e / 3.0)
        ll_het = depth * np.log(0.5)
        post = 1.0 / (1.0 + np.exp(np.clip(ll_het - ll_hom, -700, 700)))
        gq = -10.0 * np.log10(np.maximum(1.0 - post, 1e-10))
    gq = np.where(depth > 0, np.minimum(gq, 99.0), 0.0)
    return top, depth.astype(np.int64), gq
