"""
Candidate genes of the significant regions: a separate step run after the analysis, from a genome annotation (GFF3 / GTF).

Every gene that overlaps a region, widened by a flank on both sides, is listed with its position relative to the region and to the peak
marker and with the description found in the annotation. This lists positional candidates; it does not say which gene is causal.
"""

from __future__ import annotations

import gzip
import os
import re
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

from tools.insilico_bsa.backend.insilico_bsa import BsaError

MAX_ROWS = 400
GENE_TYPES = ("gene", "pseudogene", "ncrna_gene", "transposable_element_gene", "tRNA_gene", "rRNA_gene")
DESC_KEYS = ("Note", "description", "product", "Description", "function", "gene_biotype")
NAME_KEYS = ("Name", "gene_name", "symbol", "Alias")


def _attrs(text: str) -> dict:
    """GFF3 'key=value;key=value' or GTF 'key "value"; key "value"' -> dict (percent-encoding decoded)."""
    out = {}
    if "=" in text:
        for part in text.strip().split(";"):
            if "=" in part:
                k, v = part.split("=", 1)
                out[k.strip()] = unquote(v.strip())
    else:
        for k, v in re.findall(r'(\w+)\s+"([^"]*)"', text):
            out[k] = v
    return out


def read_genes(path: str) -> pd.DataFrame:
    """-> DataFrame[Chromosome, Start, End, Strand, Gene_ID, Gene_Name, Description] of the gene features of a GFF3 / GTF file."""
    from tools.recombination.backend import recombination as R
    path = str(path)
    if not os.path.isfile(path):
        raise BsaError("The annotation file was not found.")
    opener = gzip.open if path.lower().endswith(".gz") else open
    rows = []
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line or line[0] == "#":
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] not in GENE_TYPES:
                    continue
                try:
                    start, end = int(f[3]), int(f[4])
                except ValueError:
                    continue
                a = _attrs(f[8])
                gid = a.get("ID") or a.get("gene_id") or a.get("locus_tag") or ""
                gid = re.sub(r"^(gene|GENE)[:\-]", "", gid)
                name = next((a[k] for k in NAME_KEYS if a.get(k)), "")
                desc = next((a[k] for k in DESC_KEYS if a.get(k)), "")
                rows.append((f[0], min(start, end), max(start, end), f[6], gid, name if name != gid else "", desc))
    except OSError as e:
        raise BsaError(f"Could not read the annotation file: {e}") from None
    if not rows:
        raise BsaError("No gene features were found in the annotation file (expected a GFF3 / GTF with lines of type 'gene').")
    df = pd.DataFrame(rows, columns=["Seqid", "Start", "End", "Strand", "Gene_ID", "Gene_Name", "Description"])
    df["Chromosome"] = df["Seqid"].map(R.chromosome_number)
    df = df[df["Chromosome"].notna()].copy()
    if df.empty:
        raise BsaError("The sequence names of the annotation (e.g. '" + str(rows[0][0]) + "') could not be matched to chromosome numbers.")
    df["Chromosome"] = df["Chromosome"].astype(int)
    return df.drop(columns="Seqid")


def candidate_genes(regions: pd.DataFrame, genes: pd.DataFrame, flank_kb: float = 0.0) -> pd.DataFrame:
    """Genes overlapping each region widened by `flank_kb` on both sides, nearest to the peak marker first."""
    cols = ["Region", "Chromosome", "Gene_ID", "Gene_Name", "Gene_Start", "Gene_End", "Strand", "Location", "Distance_to_peak_bp", "Description"]
    flank = float(flank_kb) * 1000.0
    out = []
    for r in regions.itertuples():
        g = genes[(genes["Chromosome"] == int(r.Chromosome)) & (genes["End"] >= r.Start_bp - flank) & (genes["Start"] <= r.End_bp + flank)]
        for x in g.itertuples():
            inside = x.End >= r.Start_bp and x.Start <= r.End_bp
            dist = 0 if x.Start <= r.Peak_bp <= x.End else int(min(abs(x.Start - r.Peak_bp), abs(x.End - r.Peak_bp)))
            out.append({"Region": r.Region, "Chromosome": int(r.Chromosome), "Gene_ID": x.Gene_ID, "Gene_Name": x.Gene_Name, "Gene_Start": int(x.Start),
                        "Gene_End": int(x.End), "Strand": x.Strand, "Location": "inside the region" if inside else "in the flank",
                        "Distance_to_peak_bp": dist, "Description": x.Description})
    df = pd.DataFrame(out, columns=cols)
    if len(df):
        order = {r: i for i, r in enumerate(regions["Region"])}
        df = df.assign(_o=df["Region"].map(order)).sort_values(["_o", "Distance_to_peak_bp"]).drop(columns="_o").reset_index(drop=True)
    return df


def run(config: dict) -> dict:
    """config: regions_file (significant_regions.csv of a run), window (optional: one window size, e.g. "500kb"), gff_file, flank_kb, output_dir. Writes candidate_genes.csv / .xlsx."""
    try:
        rp = config.get("regions_file")
        if not rp or not os.path.isfile(str(rp)):
            raise BsaError("Run the analysis first: there are no regions to annotate.")
        reg = pd.read_csv(rp)
        win = str(config.get("window") or "")
        if win and "Window" in reg.columns:                            # regions of one window size (the file holds all sizes)
            if win not in set(reg["Window"].astype(str)):
                raise BsaError(f"The analysis found no significant region with {win} windows, so there is nothing to annotate.")
            reg = reg[reg["Window"].astype(str) == win].reset_index(drop=True)
        if reg.empty:
            raise BsaError("The analysis found no significant region, so there is nothing to annotate.")
        genes = read_genes(config.get("gff_file"))
        flank = float(config.get("flank_kb", 0) or 0)
        if flank < 0 or flank > 5000:
            raise BsaError("The flank must be between 0 and 5,000 kb.")
        cg = candidate_genes(reg, genes, flank)
        out = Path(config["output_dir"])
        out.mkdir(parents=True, exist_ok=True)
        csv, xlsx = out / "candidate_genes.csv", out / "candidate_genes.xlsx"
        cg.to_csv(csv, index=False)
        per = (cg.groupby("Region").agg(Genes=("Gene_ID", "count"), Inside=("Location", lambda s: int((s == "inside the region").sum())))
               if len(cg) else pd.DataFrame(columns=["Genes", "Inside"]))
        summ = reg[["Region", "Chromosome", "Start_bp", "End_bp", "Peak_marker", "Peak_bp"]].merge(per, left_on="Region", right_index=True, how="left").fillna({"Genes": 0, "Inside": 0})
        summ[["Genes", "Inside"]] = summ[["Genes", "Inside"]].astype(int)
        with pd.ExcelWriter(xlsx, engine="openpyxl") as xw:
            cg.to_excel(xw, sheet_name="Candidate genes", index=False)
            summ.to_excel(xw, sheet_name="Per region", index=False)
        rec = cg.head(MAX_ROWS).copy()
        return {"status": "success", "message": f"{len(cg):,} genes in {int((summ['Genes'] > 0).sum())} of {len(summ)} regions (flank {flank:g} kb)",
                "genes": rec.where(pd.notna(rec), None).to_dict("records"), "genes_total": int(len(cg)),
                "per_region": summ.to_dict("records"), "n_genes_in_annotation": int(len(genes)),
                "files": {"csv": str(csv), "xlsx": str(xlsx)}}
    except BsaError as e:
        return {"status": "error", "message": str(e)}
    except Exception as e:                                       # noqa: BLE001
        import traceback
        traceback.print_exc()
        return {"status": "error", "message": f"The annotation failed: {e}"}
