"""
The genotype layouts used across TerraTools.

Marker rows (Mendelian Segregation, Linkage Mapping):
    Marker | Chromosome | Ind_1 | Ind_2 ...        one row per marker
QTL layout (what Linkage Mapping exports for QTL Mapping):
    Marker     | m1 | m2 | m3 ...                  three header rows, then one row per individual
    Chromosome | 1  | 1  | 2  ...
    Position   | 0.0| 12.4| 0.0 ...
    Ind_1      | 0  | 1  | 2  ...

to_marker_rows() lets a file in the QTL layout be used where marker rows are expected.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

_TEXT = (".csv", ".tsv", ".txt")


def _sep(path):
    return "\t" if str(path).lower().endswith((".tsv", ".txt")) else ","


def is_qtl_layout(path) -> bool:
    """True when the file's first column starts with the rows 'Chromosome' and 'Position' (or 'cM')."""
    path = str(path)
    if not path.lower().endswith(_TEXT):
        return False
    try:
        head = pd.read_csv(path, sep=_sep(path), dtype=str, nrows=3, header=0)
    except Exception:
        return False
    if head.shape[0] < 2 or head.shape[1] < 3:
        return False
    first = [str(v).strip().lower() for v in head.iloc[:2, 0]]
    return first[0] in ("chromosome", "chr", "linkage group", "lg") and first[1] in ("position", "pos", "cm", "position (cm)")


def to_marker_rows(path, out_dir) -> str:
    """Write <name>_marker_rows.csv (Marker, Chromosome, one column per individual) and return its path."""
    path = str(path)
    df = pd.read_csv(path, sep=_sep(path), dtype=str, keep_default_na=False)
    markers = [str(c).strip() for c in df.columns[1:]]
    chrom = df.iloc[0, 1:].tolist()
    individuals = df.iloc[2:, 0].astype(str).tolist()
    body = df.iloc[2:, 1:].to_numpy().T
    out = pd.DataFrame(body, columns=individuals)
    out.insert(0, "Chromosome", chrom)
    out.insert(0, "Marker", markers)
    dest = Path(out_dir) / f"{Path(path).stem}_marker_rows.csv"
    out.to_csv(dest, index=False)
    return str(dest)
