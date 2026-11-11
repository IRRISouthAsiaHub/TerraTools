"""Small data generators for the tests: a FASTA, coordinate-sorted BAM files and VCF files with known allele counts."""

from __future__ import annotations

import struct
from pathlib import Path

import numpy as np

from tools.qtl_seq.backend.bgzf import BgzfWriter

BASES = "ACGT"
CODE = {b: i for i, b in enumerate("=ACMGRSVTWYHKDBN")}


def reg2bin(beg, end):
    end -= 1
    if beg >> 14 == end >> 14:
        return ((1 << 15) - 1) // 7 + (beg >> 14)
    if beg >> 17 == end >> 17:
        return ((1 << 12) - 1) // 7 + (beg >> 17)
    if beg >> 20 == end >> 20:
        return ((1 << 9) - 1) // 7 + (beg >> 20)
    if beg >> 23 == end >> 23:
        return ((1 << 6) - 1) // 7 + (beg >> 23)
    if beg >> 26 == end >> 26:
        return ((1 << 3) - 1) // 7 + (beg >> 26)
    return 0


def random_genome(contigs, seed=0):
    rng = np.random.default_rng(seed)
    return {n: "".join(rng.choice(list(BASES), L)) for n, L in contigs.items()}


def write_fasta(path, genome, width=60):
    with open(path, "w") as f:
        for n, s in genome.items():
            f.write(f">{n}\n")
            for i in range(0, len(s), width):
                f.write(s[i:i + width] + "\n")


def write_bam(path, contigs, reads, sort_order="coordinate"):
    """reads: list of (contig name, 0-based pos, sequence, mapq, flag, base quality). All reads are M (match) only."""
    names = list(contigs)
    ids = {n: i for i, n in enumerate(names)}
    reads = sorted(reads, key=lambda r: (ids[r[0]], r[1]))
    with BgzfWriter(path) as w:
        text = f"@HD\tVN:1.6\tSO:{sort_order}\n" + "".join(f"@SQ\tSN:{n}\tLN:{L}\n" for n, L in contigs.items())
        body = b"BAM\x01" + struct.pack("<i", len(text)) + text.encode() + struct.pack("<i", len(names))
        for n, L in contigs.items():
            body += struct.pack("<i", len(n) + 1) + n.encode() + b"\x00" + struct.pack("<i", L)
        w.write(body)
        for k, (chrom, pos, seq, mapq, flag, bq) in enumerate(reads):
            name = f"r{k}".encode() + b"\x00"
            n = len(seq)
            cigar = struct.pack("<I", (n << 4) | 0)
            codes = [CODE[b] for b in seq]
            if len(codes) % 2:
                codes.append(0)
            packed = bytes((codes[i] << 4) | codes[i + 1] for i in range(0, len(codes), 2))
            rec = struct.pack("<iiBBHHHiiii", ids[chrom], pos, len(name), mapq, reg2bin(pos, pos + n), 1, flag, n, -1, -1, 0)
            rec += name + cigar + packed + bytes([bq]) * n
            w.write(struct.pack("<i", len(rec)) + rec)


def reads_at(genome, chrom, pos1, allele_counts, read_len=60, mapq=60, bq=40, rng=None):
    """Reads that all cover position `pos1` (1-based) and carry the given alleles there: {allele base: n reads}."""
    rng = rng or np.random.default_rng(0)
    out = []
    ref = genome[chrom]
    for base, n in allele_counts.items():
        for _ in range(int(n)):
            off = int(rng.integers(5, read_len - 5))
            start = max(0, pos1 - 1 - off)
            seq = list(ref[start:start + read_len])
            seq[pos1 - 1 - start] = base
            out.append((chrom, start, "".join(seq), mapq, 0, bq))
    return out


def other_base(b):
    return BASES[(BASES.index(b) + 1) % 4]


def build_bam_dataset(tmp, seed=3, depth=30, sites_per_contig=None, qtl=None, qtl_p=(0.9, 0.1)):
    """Parents (P1 = reference allele, P2 = alternative allele at every parental SNP), a high and a low bulk.

    Returns (cfg paths dict, truth dict). `truth['sites']` lists (chrom, pos, ref base, alt base, kind) where kind is
    'poly' (parents differ), 'shared' (both parents carry the alt allele: not informative) or 'het' (parent 1 heterozygous).
    `qtl` = (chrom, start, end): bulk allele share qtl_p inside, 0.5 outside."""
    tmp = Path(tmp)
    rng = np.random.default_rng(seed)
    contigs = {"chrA": 200_000, "chrB": 120_000}
    sites_per_contig = sites_per_contig or {"chrA": 150, "chrB": 90}
    genome = random_genome(contigs, seed)
    write_fasta(tmp / "ref.fa", genome)
    sites = []
    for c, n in sites_per_contig.items():
        pos = np.sort(rng.choice(np.arange(500, contigs[c] - 500, 700), n, replace=False))
        for i, p in enumerate(pos):
            ref = genome[c][p - 1]
            kind = "poly"
            if i % 17 == 5:
                kind = "shared"
            elif i % 23 == 7:
                kind = "het"
            sites.append((c, int(p), ref, other_base(ref), kind))
    reads = {"P1": [], "P2": [], "H": [], "L": []}
    for c, p, ref, alt, kind in sites:
        if kind == "het":
            reads["P1"] += reads_at(genome, c, p, {ref: depth // 2, alt: depth - depth // 2}, rng=rng)
        else:
            reads["P1"] += reads_at(genome, c, p, {ref: depth}, rng=rng) if kind != "shared" else reads_at(genome, c, p, {alt: depth}, rng=rng)
        reads["P2"] += reads_at(genome, c, p, {alt: depth}, rng=rng)
        inq = qtl is not None and c == qtl[0] and qtl[1] <= p <= qtl[2]
        for bulk, pa in (("H", qtl_p[0] if inq else 0.5), ("L", qtl_p[1] if inq else 0.5)):
            na = int(rng.binomial(depth, pa))
            reads[bulk] += reads_at(genome, c, p, {ref: depth - na, alt: na}, rng=rng)
    paths = {"ref_fasta": str(tmp / "ref.fa")}
    for key, name in (("P1", "parent1_bam"), ("P2", "parent2_bam"), ("H", "high_bam"), ("L", "low_bam")):
        write_bam(tmp / f"{key}.bam", contigs, reads[key])
        paths[name] = str(tmp / f"{key}.bam")
    return paths, {"sites": sites, "contigs": contigs, "genome": genome}


def write_vcf(path, records, samples, contigs=None, fmt="GT:AD:DP:GQ", gz=False):
    """records: list of (chrom, pos, ref, alt, [sample field strings], qual, filter)."""
    lines = ["##fileformat=VCFv4.2"]
    for c, L in (contigs or {}).items():
        lines.append(f"##contig=<ID={c},length={L}>")
    lines.append('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">')
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples))
    for chrom, pos, ref, alt, fields, qual, flt in records:
        lines.append("\t".join([chrom, str(pos), ".", ref, alt, str(qual), flt, ".", fmt] + list(fields)))
    text = "\n".join(lines) + "\n"
    if gz:
        with BgzfWriter(path) as w:
            w.write(text)
    else:
        Path(path).write_text(text)
