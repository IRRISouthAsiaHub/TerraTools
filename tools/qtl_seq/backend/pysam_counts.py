"""
Allele counting with pysam (A, C, G, T read counts per site and BAM).

It is kept free of numpy and of other TerraTools modules so that it can also run as a stand-alone script with the Python of
WSL (pysam has no Windows build):

    python3 pysam_counts.py --chrom chr1 --mapq 20 --baseq 20 [--positions sites.txt] a.bam b.bam ...

Output: one line per covered site: the 1-based position, then A, C, G, T counts of BAM 1, of BAM 2, ... (tab separated).

Reads are filtered like `samtools mpileup`: unmapped, secondary, QC-fail and duplicate reads are skipped, as are reads below the
mapping quality and bases below the base quality. Unlike mpileup, the two mates of an overlapping pair are both counted, and
no base-alignment-quality recalculation is done (mpileup -B).
"""

from __future__ import annotations

import sys

SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x400          # unmapped, secondary, QC fail, duplicate
BATCH = 100_000


def _open(pysam, path):
    f = pysam.AlignmentFile(path, "rb")
    if not f.has_index():
        f.close()
        pysam.index(path)
        f = pysam.AlignmentFile(path, "rb")
    return f


def spans(length, positions, tile):
    """0-based half-open [start, end) pieces to count: tiles of the chromosome, or the tiles that hold requested positions."""
    if positions is None:
        return [(s, min(s + tile, length)) for s in range(0, length, tile)]
    out, cur = [], None
    for p in positions:
        t = (p - 1) // tile
        if cur is None or t != cur[0]:
            cur = [t, p - 1, p]
            out.append(cur)
        else:
            cur[2] = p
    return [(a, b) for _t, a, b in out]


def iter_counts(bams, chrom, min_mapq, min_baseq, positions=None, tile=1_000_000, should_stop=None):
    """Yields (list of 1-based positions, list of rows); a row holds 4 counts (A, C, G, T) for every BAM, in order.

    positions: sorted list of 1-based positions to report, or None for every covered site."""
    import pysam

    files = [_open(pysam, b) for b in bams]
    keep = lambda r: r.mapping_quality >= min_mapq and not (r.flag & SKIP_FLAGS)     # noqa: E731
    try:
        length = files[0].get_reference_length(chrom)
        if length is None:
            raise ValueError(f"{chrom} is not in {bams[0]}")
        wanted = None if positions is None else set(positions)
        pos, rows = [], []
        for start, end in spans(length, None if positions is None else sorted(positions), tile):
            if should_stop and should_stop():
                return
            per = []
            for f in files:
                a, c, g, t = f.count_coverage(chrom, start, end, quality_threshold=min_baseq, read_callback=keep)
                per += [a, c, g, t]
            for i, vals in enumerate(zip(*per)):
                p = start + i + 1
                if wanted is not None:
                    if p not in wanted:
                        continue
                elif not any(vals):
                    continue
                pos.append(p)
                rows.append(vals)
                if len(pos) >= BATCH:
                    yield pos, rows
                    pos, rows = [], []
        if pos:
            yield pos, rows
    finally:
        for f in files:
            f.close()


def main(argv):
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--mapq", type=int, default=20)
    ap.add_argument("--baseq", type=int, default=20)
    ap.add_argument("--positions", default="")
    ap.add_argument("bams", nargs="+")
    a = ap.parse_args(argv)
    positions = None
    if a.positions:
        with open(a.positions) as fh:
            positions = [int(x.split("\t")[-1]) for x in fh if x.strip()]
    out = sys.stdout
    for pos, rows in iter_counts(a.bams, a.chrom, a.mapq, a.baseq, positions):
        out.write("".join(str(p) + "\t" + "\t".join(map(str, r)) + "\n" for p, r in zip(pos, rows)))
    out.flush()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
