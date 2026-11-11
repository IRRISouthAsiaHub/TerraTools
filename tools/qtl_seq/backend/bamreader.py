"""
Built-in BAM reader (numpy only): base counts per reference position, used when samtools is not available.

A BAM is BGZF = multi-member gzip, so the standard library can decompress it. A coordinate-sorted BAM is streamed
once and the reads are piled up tile by tile. Counting rules, as in `samtools mpileup -q <mapq> -Q <baseq>`:
reads that are unmapped, secondary, QC-fail, duplicates or supplementary are skipped; mapping quality >= min_mapq;
base quality >= min_baseq; insertions and soft clips are ignored; deletions and N-skips are not alleles.
Not reproduced: BAQ and the overlapping-mate correction. It is slower than samtools (thousands of reads per second),
fine for small genomes, targeted data and tests; use samtools for whole genomes.
"""

from __future__ import annotations

import gzip
import io
import struct

import numpy as np

_LUT = np.full(16, 4, dtype=np.uint8)               # BAM 4-bit base codes =ACMGRSVTWYHKDBN -> A0 C1 G2 T3, else 4
_LUT[1], _LUT[2], _LUT[4], _LUT[8] = 0, 1, 2, 3
_EXCLUDE = 0x4 | 0x100 | 0x200 | 0x400 | 0x800
_REC_HEAD = struct.Struct("<iiBBHHHiiii")


class BamError(RuntimeError):
    pass


def read_header(path):
    """(contigs [(name, length)], header text). BamError when the file is not a BAM."""
    try:
        with gzip.open(str(path), "rb") as f:
            if f.read(4) != b"BAM\x01":
                raise BamError(f"{path} is not a BAM file (it must be .bam, not SAM or CRAM).")
            (l_text,) = struct.unpack("<i", f.read(4))
            text = f.read(l_text).decode("utf-8", "replace")
            (n_ref,) = struct.unpack("<i", f.read(4))
            refs = []
            for _ in range(n_ref):
                (l_name,) = struct.unpack("<i", f.read(4))
                name = f.read(l_name).rstrip(b"\x00").decode("ascii", "replace")
                (l_ref,) = struct.unpack("<i", f.read(4))
                refs.append((name, l_ref))
        return refs, text
    except (OSError, EOFError, struct.error):
        raise BamError(f"{path} is not a BAM file (it must be a BGZF-compressed .bam, not SAM or CRAM).") from None


class BamCounter:
    """Streams one coordinate-sorted BAM and returns the A/C/G/T counts of consecutive tiles."""

    def __init__(self, path, min_mapq=20, min_baseq=20):
        self.path = str(path)
        self.min_mapq, self.min_baseq = int(min_mapq), int(min_baseq)
        self.refs, self.header_text = read_header(path)
        self._f = io.BufferedReader(gzip.GzipFile(self.path, "rb"), 1 << 22)
        self._skip_header()
        self._pending = None
        self._eof = False
        self._last = (-1, -1)
        self.n_reads = 0
        self._buf = None
        self._buf_ref = None
        self._origin = 0

    def _skip_header(self):
        f = self._f
        f.read(4)
        (l_text,) = struct.unpack("<i", f.read(4))
        f.read(l_text)
        (n_ref,) = struct.unpack("<i", f.read(4))
        for _ in range(n_ref):
            (l_name,) = struct.unpack("<i", f.read(4))
            f.read(l_name + 4)

    def _next_record(self):
        raw = self._f.read(4)
        if len(raw) < 4:
            return None
        (size,) = struct.unpack("<i", raw)
        data = self._f.read(size)
        ref_id, pos, l_name, mapq, _bin, n_cigar, flag, l_seq, *_ = _REC_HEAD.unpack_from(data, 0)
        return ref_id, pos, mapq, flag, data, n_cigar, l_name, l_seq

    def _peek(self):
        if self._pending is None and not self._eof:
            rec = self._next_record()
            if rec is None:
                self._eof = True
            else:
                if rec[0] >= 0 and (rec[0], rec[1]) < self._last:
                    raise BamError(f"{self.path} is not coordinate-sorted. Sort it first (samtools sort).")
                self._last = (rec[0], rec[1])
                self._pending = rec
        return self._pending

    def _grow(self, need):
        if need > self._buf.shape[1]:
            extra = max(need - self._buf.shape[1], self._buf.shape[1])
            self._buf = np.pad(self._buf, ((0, 0), (0, extra)))

    def _add_read(self, rec):
        _ref_id, pos, mapq, flag, data, n_cigar, l_name, l_seq = rec
        if (flag & _EXCLUDE) or mapq < self.min_mapq or l_seq == 0 or n_cigar == 0:
            return
        off = 32 + l_name
        cig = struct.unpack_from(f"<{n_cigar}I", data, off)
        off += 4 * n_cigar
        span = sum(c >> 4 for c in cig if (c & 15) in (0, 2, 3, 7, 8))
        rp = pos - self._origin
        if rp < 0:
            return
        self._grow(rp + span + 1)
        packed = np.frombuffer(data, np.uint8, (l_seq + 1) // 2, off)
        off += (l_seq + 1) // 2
        codes = np.empty(2 * packed.size, dtype=np.uint8)
        codes[0::2], codes[1::2] = packed >> 4, packed & 15
        base = _LUT[codes[:l_seq]]
        qual = np.frombuffer(data, np.uint8, l_seq, off)
        base = np.where((qual >= self.min_baseq) | (qual == 255), base, 4)
        qi, buf = 0, self._buf
        for c in cig:
            op, ln = c & 15, c >> 4
            if op in (0, 7, 8):
                seg = base[qi:qi + ln]
                for k in range(4):
                    buf[k, rp:rp + ln] += seg == k
                qi += ln
                rp += ln
            elif op in (1, 4):
                qi += ln
            elif op in (2, 3):
                rp += ln

    def tile(self, ref_id, start, end):
        """Counts (4, end-start) uint16 for reference `ref_id`, 0-based positions [start, end). Call with consecutive
        tiles of one reference, then the next reference."""
        if self._buf is None or self._buf_ref != ref_id:
            self._buf = np.zeros((4, (end - start) + 1000), dtype=np.uint16)
            self._buf_ref, self._origin = ref_id, start
        while True:                                   # reads of references that were not asked for are skipped
            rec = self._peek()
            if rec is not None and 0 <= rec[0] < ref_id:
                self._pending = None
                continue
            break
        while True:
            rec = self._peek()
            if rec is None or rec[0] != ref_id or rec[1] >= end:
                break
            self._pending = None
            self._add_read(rec)
            self.n_reads += 1
        n = end - start
        self._grow(n)
        out = self._buf[:, :n].copy()
        carry = self._buf[:, n:]
        self._buf = np.pad(carry, ((0, 0), (0, max(n + 1000 - carry.shape[1], 0))))
        self._origin = end
        return out

    def close(self):
        self._f.close()
