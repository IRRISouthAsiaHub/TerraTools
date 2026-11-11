"""
Validation mode.

REFERENCE holds the design and the findings reported for two rice experiments. They are expectations to compare
with, never results that the program returns.

The original sequencing data are not part of the program, so validation runs on SIMULATED data: SNPs, read depths
and bulks are generated for the same design (population type, plants per bulk, QTL position and size) and sent
through the same filters, simulation, windows and QTL calling as real data. This shows that the pipeline recovers a
QTL of the reported size at the reported place with the reported kind of significance. It is not a replication of the
published numbers: that needs the original reads, aligned and filtered the same way.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from tools.qtl_seq import config as cfgmod
from tools.qtl_seq.backend import pipeline, stats

# rice chromosome lengths (IRGSP-1.0)
RICE = {"chr1": 43_270_923, "chr2": 35_937_250, "chr3": 36_413_819, "chr4": 35_502_694, "chr5": 29_958_434,
        "chr6": 31_248_787, "chr7": 29_697_621, "chr8": 28_443_022, "chr9": 23_012_720, "chr10": 23_207_287,
        "chr11": 29_021_106, "chr12": 27_531_856}

REFERENCE = {
    "blast_ril": {
        "title": "Experiment 1 — rice blast resistance (RIL)",
        "cross": "Nortai × Hitomebore", "population": "RIL", "generation": "F7", "total_lines": 241,
        "n_high": 20, "n_low": 20, "high_name": "Resistant bulk", "low_name": "Susceptible bulk",
        "qtls": [{"name": "qPi-nor1(t)", "chromosome": "chr6", "start": 2_390_000, "end": 4_390_000, "max_delta": 0.79,
                  "level": 99, "p": "P < 0.01"}],
        "linkage_lod": 7.38,
    },
    "vigour_f2": {
        "title": "Experiment 2 — seedling vigour (F2)",
        "cross": "Dunghan Shali × Hitomebore", "population": "F2", "generation": "F2", "total_lines": 531,
        "n_high": 50, "n_low": 50, "high_name": "High-vigour bulk", "low_name": "Low-vigour bulk",
        "qtls": [{"name": "Chromosome 3 interval", "chromosome": "chr3", "start": 36_210_000, "end": 37_310_000,
                  "max_delta": 0.61, "level": 99, "p": "P < 0.01"},
                 {"name": "Chromosome 1 interval", "chromosome": "chr1", "start": 39_080_000, "end": 41_080_000,
                  "max_delta": 0.67, "level": 95, "p": "P < 0.05"}],
    },
}


def reference_summary():
    return [{"key": k, **{kk: vv for kk, vv in v.items()}} for k, v in REFERENCE.items()]


def simulate_experiment(key, seed=1, snp_spacing=4000, mean_depth=12, cm_bp=250_000, block=1_000_000):
    """Simulated SNP table (read counts of two bulks) for a reference design.

    Allele share of each bulk: 1/2 away from a QTL, pulled towards the index allele (high bulk) or away from it (low bulk)
    near a QTL, decaying with the recombination fraction (1 cM per `cm_bp`). The share of the n plants is a binomial draw,
    changing every `block` bp; reads are binomial draws from that share with a negative-binomial depth."""
    ref = REFERENCE[key]
    rng = np.random.default_rng(seed)
    pop, nh, nl = ref["population"], ref["n_high"], ref["n_low"]
    rows = []
    # A reported interval can lie beyond the end of the IRGSP-1.0 chromosome (chromosome 3: 37.31 Mb reported, 36.41 Mb in
    # that assembly, because another reference was used), so the simulated chromosome is lengthened to hold it.
    lengths = {c: max(L, max([q["end"] for q in ref["qtls"] if q["chromosome"] == c], default=0) + 2_000_000)
               for c, L in RICE.items()}
    for chrom, length in lengths.items():
        n = int(length / snp_spacing)
        pos = np.sort(rng.choice(np.arange(1, length), n, replace=False))
        nb = int(length // block) + 1
        qtl = [q for q in ref["qtls"] if q["chromosome"] == chrom]
        bid = pos // block
        mid = bid * block + block / 2

        def share(bulk_n, sign):
            p = np.full(nb, 0.5)
            centres = np.arange(nb) * block + block / 2
            for q in qtl:
                # the effect at the QTL is set so that the window mean reaches about the reported maximum
                a = min(0.98, q["max_delta"] / 0.93)
                d = np.abs(centres - (q["start"] + q["end"]) / 2) / cm_bp / 100.0           # Morgans
                p = p + sign * 0.5 * a * np.exp(-2 * d)
            p = np.clip(p, 0.001, 0.999)
            # realised share of the index allele among the plants of the bulk
            return (rng.binomial(2 * bulk_n, p) / (2.0 * bulk_n)) if pop == "F2" else (rng.binomial(bulk_n, p) / float(bulk_n))

        fh, fl = share(nh, +1), share(nl, -1)
        dh = rng.negative_binomial(8, 8 / (8 + mean_depth), n)
        dl = rng.negative_binomial(8, 8 / (8 + mean_depth), n)
        ah, al = rng.binomial(dh, fh[bid]), rng.binomial(dl, fl[bid])
        rows.append(pd.DataFrame({"Chromosome": chrom, "Position": pos, "Ref": "A", "Alt": "G",
                                  "Ref_High": dh - ah, "Alt_High": ah, "Ref_Low": dl - al, "Alt_Low": al}))
    return pd.concat(rows, ignore_index=True), lengths


def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def compare(qtls, key, tolerance=0.15, slack=2_000_000):
    """Compare called QTLs with the reference findings of experiment `key`. Returns a list of checks."""
    out = []
    ref = REFERENCE[key]
    used = set()
    for e in ref["qtls"]:
        hits = qtls[qtls["Chromosome"] == e["chromosome"]] if len(qtls) else qtls
        best, best_ov = None, 0
        for i, q in (hits.iterrows() if len(hits) else []):
            ov = overlap(e["start"], e["end"], q["Start"], q["End"])
            if ov > best_ov or (best is None and e["start"] <= q["Peak_Position"] <= e["end"]):
                best, best_ov = q, ov
                used.add(q["QTL_ID"])
        item = {"expected": e["name"], "chromosome": e["chromosome"], "expected_interval": f"{e['start'] / 1e6:.2f}–{e['end'] / 1e6:.2f} Mb",
                "expected_max_delta": e["max_delta"], "expected_significance": e["p"], "expected_level": e["level"]}
        if best is None:
            item.update(found=False, passed=False, note="No QTL interval overlaps the expected region.")
        else:
            lvl = int(str(best["Significance_Level"]).rstrip("%"))
            # the peak of a 2 Mb sliding mean is only located to within about one window
            peak_in = e["start"] - slack <= best["Peak_Position"] <= e["end"] + slack
            dd = abs(abs(best["Peak_Delta_SNP_index"]) - e["max_delta"])
            item.update(found=True, called_interval=f"{best['Start'] / 1e6:.2f}–{best['End'] / 1e6:.2f} Mb",
                        peak_position=int(best["Peak_Position"]), peak_in_expected=bool(peak_in),
                        peak_delta=round(float(best["Peak_Delta_SNP_index"]), 3), delta_difference=round(float(dd), 3),
                        level=lvl, level_ok=lvl >= e["level"], delta_ok=bool(dd <= tolerance),
                        passed=bool(peak_in and lvl >= e["level"] and dd <= tolerance), qtl_id=best["QTL_ID"])
            item["note"] = ("" if item["passed"] else
                            "; ".join(x for x in (None if peak_in else "peak more than one window away from the expected interval",
                                                  None if lvl >= e["level"] else f"only {lvl}% reached",
                                                  None if dd <= tolerance else f"Δ differs by {dd:.2f}") if x))
        out.append(item)
    extra = int(len(qtls) - len(used)) if len(qtls) else 0
    out.append({"expected": "Other QTL intervals", "chromosome": "", "expected_interval": "none reported", "found": extra == 0,
                "passed": True, "info": True,
                "note": "none" if extra == 0 else f"{extra} additional interval(s) elsewhere or on the same chromosome. The limits are "
                "not corrected for the many windows tested, so occasional extra intervals are expected."})
    return out


def run_validation(key, seed=1, overrides=None, progress=None):
    """Simulate experiment `key`, analyse it with the reference design, compare. Returns a JSON-ready dict."""
    if key not in REFERENCE:
        return {"status": "error", "message": f"Unknown validation experiment: {key}"}
    ref = REFERENCE[key]
    p = cfgmod.merged({"population": ref["population"], "n_high": ref["n_high"], "n_low": ref["n_low"], "seed": seed,
                       "high_name": ref["high_name"], "low_name": ref["low_name"], "qtl_min_level": 95, **(overrides or {})})
    p = pipeline.coerce(p)
    snps, lengths = simulate_experiment(key, seed)
    res = pipeline.analyze_snps(snps, lengths, p, log=(lambda m, pct=None: progress(pct or 0, m)) if progress else (lambda *a, **k: None))
    checks = compare(res["qtls"], key)
    return {"status": "success", "key": key, "title": ref["title"], "reference": ref, "checks": checks,
            "passed": all(c["passed"] for c in checks), "qtls": pipeline._records(res["qtls"]),
            "plot": pipeline.plot_payload(res), "counts": res["counts"], "ledger": res["ledger"],
            "note": ("Simulated data for this design (not the original reads): the check shows that a QTL of the reported size and "
                     "position is recovered by this pipeline. It does not reproduce the published numbers. Chromosome lengths are "
                     "those of IRGSP-1.0, lengthened where the reported interval lies beyond the end of the chromosome.")}
