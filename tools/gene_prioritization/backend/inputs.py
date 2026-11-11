"""
Reading the user's input: the genes to rank.

  Gene IDs       TXT / CSV / TSV / XLSX file and/or pasted text -> a list of identifiers
                 (for example the Gene_ID_list.txt written by Genes Extraction)
"""

from __future__ import annotations

import re

import pandas as pd


class InputError(Exception):
    pass


def _norm(h):
    return re.sub(r"[^a-z0-9]", "", str(h).lower())


# ===================================================================== gene lists
GENE_HEADERS = {"gene", "geneid", "genename", "locus", "id", "geneids", "gene_id", "ensemblid", "symbol", "name", "candidategene"}


def split_ids(text):
    out = []
    for tok in re.split(r"[\s,;]+", text or ""):
        tok = tok.strip().strip("\"'")
        if tok and tok.upper() not in ("NA", "N/A", "NAN", "NULL", "NONE"):
            out.append(tok)
    return out


def read_gene_ids(path="", pasted=""):
    """-> (ids, notes). Order is kept, repeats are dropped."""
    notes, ids = [], []
    if path:
        low = str(path).lower()
        if low.endswith((".xlsx", ".xls")):
            frames = list(pd.read_excel(path, header=None, dtype=str, sheet_name=None).values())
            df = frames[0] if frames else pd.DataFrame()
            ids = _ids_from_frame(df)
        elif low.endswith((".csv", ".tsv")):
            df = pd.read_csv(path, header=None, dtype=str, sep="\t" if low.endswith(".tsv") else ",",
                             engine="python", encoding="utf-8-sig", on_bad_lines="skip")
            ids = _ids_from_frame(df)
        else:
            with open(path, encoding="utf-8-sig", errors="replace") as fh:
                ids = split_ids(fh.read())
            if ids and _norm(ids[0]) in GENE_HEADERS:
                ids = ids[1:]
    if pasted and pasted.strip():
        ids += split_ids(pasted)
    seen, uniq = set(), []
    for i in ids:
        if i not in seen:
            seen.add(i)
            uniq.append(i)
    if len(uniq) < len(ids):
        notes.append(f"{len(ids) - len(uniq)} repeated gene ID(s) were counted once.")
    return uniq, notes


def _ids_from_frame(df):
    if df.empty:
        return []
    head = [_norm(c) for c in df.iloc[0]]
    col = next((i for i, h in enumerate(head) if h in GENE_HEADERS), None)
    if col is not None:
        return split_ids(" ".join(str(v) for v in df.iloc[1:, col].dropna()))
    return split_ids(" ".join(str(v) for v in df.iloc[:, 0].dropna()))
