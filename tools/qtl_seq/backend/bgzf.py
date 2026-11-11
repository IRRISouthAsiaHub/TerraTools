"""Minimal BGZF writer: the block-gzip format that samtools, bcftools and tabix read (.vcf.gz, .bam)."""

from __future__ import annotations

import struct
import zlib

_BLOCK = 0xFF00                                     # uncompressed bytes per block
_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def _block(data: bytes) -> bytes:
    co = zlib.compressobj(6, zlib.DEFLATED, -15)
    body = co.compress(data) + co.flush()
    size = len(body) + 25                           # header 18 + trailer 8 - 1
    head = struct.pack("<BBBBIBBHBBHH", 31, 139, 8, 4, 0, 0, 255, 6, 66, 67, 2, size)
    return head + body + struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data))


class BgzfWriter:
    def __init__(self, path):
        self._f = open(path, "wb")
        self._buf = bytearray()

    def write(self, data):
        if isinstance(data, str):
            data = data.encode("utf-8")
        self._buf += data
        while len(self._buf) >= _BLOCK:
            self._f.write(_block(bytes(self._buf[:_BLOCK])))
            del self._buf[:_BLOCK]

    def close(self):
        if self._buf:
            self._f.write(_block(bytes(self._buf)))
            self._buf.clear()
        self._f.write(_EOF)
        self._f.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
