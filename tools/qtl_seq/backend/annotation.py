"""Gene annotation (GFF3 / GTF) and the genes that overlap a QTL interval."""

from __future__ import annotations

import gzip
import io
import re

import numpy as np
import pandas as pd

from tools.common.seqtools import _attrs, _strip_prefix, norm_chrom

_GTF = re.compile(r'(\w+)\s+"([^"]*)"')
DESC_KEYS = ("description", "Description", "product", "Note", "note", "function", "gene_description", "gene_biotype", "biotype")
NAME_KEYS = ("Name", "gene_name", "symbol", "Alias")


def _open(path):
    p = str(path)
    if p.lower().endswith((".gz", ".bgz")):
        return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8", errors="replace")
    return open(p, "r", encoding="utf-8", errors="replace")


def load_genes(path):
    """Genes of a GFF3 or GTF file as a DataFrame: Chromosome, Gene_ID, Start, End, Strand, Name, Annotation.
    Genes are the `gene` features; for files without them the transcripts are grouped by their gene."""
    rows, tx = [], {}
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            kind = f[2]
            gtf = '"' in f[8] and "=" not in f[8].split(";")[0]
            a = {k: v for k, v in _GTF.findall(f[8])} if gtf else _attrs(f[8])
            try:
                s, e = int(f[3]), int(f[4])
            except ValueError:
                continue
            gid = a.get("gene_id") if gtf else a.get("ID")
            if kind == "gene" and gid:
                desc = next((a[k] for k in DESC_KEYS if a.get(k)), "")
                name = next((a[k] for k in NAME_KEYS if a.get(k)), "")
                rows.append((f[0], _strip_prefix(gid), s, e, f[6], name, desc))
            elif kind in ("mRNA", "transcript") and not gtf:
                par = a.get("Parent")
                if par:
                    g = tx.setdefault(_strip_prefix(par.split(",")[0]), [f[0], 10**12, 0, f[6], "", ""])
                    g[1], g[2] = min(g[1], s), max(g[2], e)
                    g[4] = g[4] or next((a[k] for k in NAME_KEYS if a.get(k)), "")
                    g[5] = g[5] or next((a[k] for k in DESC_KEYS if a.get(k)), "")
    if not rows and tx:
        rows = [(v[0], k, v[1], v[2], v[3], v[4], v[5]) for k, v in tx.items()]
    df = pd.DataFrame(rows, columns=["Chromosome", "Gene_ID", "Start", "End", "Strand", "Name", "Annotation"])
    df["_norm"] = df["Chromosome"].map(norm_chrom)
    return df


def candidate_genes(qtls, genes):
    """Genes that overlap each QTL interval. These are candidates by position only; nothing here says a gene is causal.
    Distance_From_Peak is 0 when the peak lies inside the gene, otherwise the bp to its nearest end."""
    cols = ["QTL_ID", "Chromosome", "Gene_ID", "Gene_Name", "Gene_Start", "Gene_End", "Strand", "Annotation",
            "Distance_From_Peak", "Gene_Contains_Peak", "QTL_Start", "QTL_End", "Peak_Position"]
    out = []
    if genes is None or genes.empty or qtls.empty:
        return pd.DataFrame(columns=cols)
    for q in qtls.itertuples(index=False):
        g = genes[(genes["_norm"] == norm_chrom(q.Chromosome)) & (genes["End"] >= q.Start) & (genes["Start"] <= q.End)]
        for r in g.itertuples(index=False):
            inside = r.Start <= q.Peak_Position <= r.End
            dist = 0 if inside else int(min(abs(q.Peak_Position - r.Start), abs(q.Peak_Position - r.End)))
            out.append((q.QTL_ID, q.Chromosome, r.Gene_ID, r.Name, r.Start, r.End, r.Strand, r.Annotation, dist,
                        bool(inside), q.Start, q.End, q.Peak_Position))
    df = pd.DataFrame(out, columns=cols)
    if not len(df):
        return df
    df["_n"] = df["QTL_ID"].str.extract(r"(\d+)", expand=False).astype(int)         # QTL2 before QTL10
    return df.sort_values(["_n", "Distance_From_Peak"], kind="stable").drop(columns="_n").reset_index(drop=True)
