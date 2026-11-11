"""
Sequence and annotation helpers shared by the tools that read a reference FASTA or a GFF3: chromosome-name
normalisation, random access to an indexed FASTA, reverse complement and GFF attribute parsing.
"""

from __future__ import annotations

from pathlib import Path
from urllib.parse import unquote

_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def norm_chrom(c) -> str:
    """'chr01', 'Chr1', 'chromosome_1', '01' and '1' are the same chromosome."""
    s = str(c).strip()
    low = s.lower()
    for pre in ("chromosome", "chrom", "chr"):
        if low.startswith(pre):
            s = s[len(pre):]
            break
    s = s.lstrip("_-. ")
    return str(int(s)) if s.isdigit() else s.upper()


# ============================================================
# FASTA (random access through a .fai index)
# ============================================================

class Fasta:
    def __init__(self, path):
        self.path = str(path)
        if self.path.lower().endswith((".gz", ".bgz")):
            raise ValueError("A compressed FASTA cannot be read at random. Unzip the reference (.fa / .fasta) first.")
        self._idx = {}
        for name, length, off, lb, lw in self._load_index():
            self._idx[norm_chrom(name)] = (name, length, off, lb, lw)
        if not self._idx:
            raise ValueError(f"{self.path} has no sequences.")
        self._fh = open(self.path, "rb")

    def _load_index(self):
        fai = Path(self.path + ".fai")
        out = []
        if fai.exists():
            for line in fai.read_text().splitlines():
                f = line.split("\t")
                if len(f) >= 5:
                    out.append((f[0], int(f[1]), int(f[2]), int(f[3]), int(f[4])))
            if out:
                return out
        name = None
        offset = 0
        length = seq_off = 0
        lb = lw = None
        with open(self.path, "rb") as fh:
            for line in fh:
                n = len(line)
                if line.startswith(b">"):
                    if name is not None and lb:
                        out.append((name, length, seq_off, lb, lw))
                    name = line[1:].split()[0].decode("ascii", "replace")
                    length, seq_off, lb, lw = 0, offset + n, None, None
                else:
                    s = line.rstrip(b"\r\n")
                    if lb is None and s:
                        lb, lw = len(s), n
                    length += len(s)
                offset += n
            if name is not None and lb:
                out.append((name, length, seq_off, lb, lw))
        try:
            fai.write_text("".join(f"{a}\t{b}\t{c}\t{d}\t{e}\n" for a, b, c, d, e in out))
        except OSError:
            pass
        return out

    def contigs(self) -> dict:
        """{sequence name: length} as written in the FASTA."""
        return {name: length for name, length, *_ in self._idx.values()}

    def has(self, chrom) -> bool:
        return norm_chrom(chrom) in self._idx

    def fetch(self, chrom, start: int, end: int) -> str:
        """Bases start..end (1-based, inclusive), upper case; '' when out of range."""
        e = self._idx.get(norm_chrom(chrom))
        if e is None:
            return ""
        _name, length, off, lb, lw = e
        start, end = max(start, 1), min(end, length)
        if end < start:
            return ""
        s0, e0 = start - 1, end
        l0, l1 = s0 // lb, (e0 - 1) // lb
        b0 = off + l0 * lw + (s0 - l0 * lb)
        b1 = off + l1 * lw + ((e0 - 1) - l1 * lb) + 1
        self._fh.seek(b0)
        return self._fh.read(b1 - b0).replace(b"\n", b"").replace(b"\r", b"").decode("ascii", "replace").upper()

    def close(self):
        self._fh.close()


def _attrs(s: str) -> dict:
    d = {}
    for kv in s.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k.strip()] = unquote(v) if "%" in v else v
    return d


def _strip_prefix(s: str) -> str:
    for pre in ("gene:", "gene-", "transcript:", "rna-"):
        if s.startswith(pre):
            return s[len(pre):]
    return s
