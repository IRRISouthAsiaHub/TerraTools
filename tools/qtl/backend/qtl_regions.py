"""
qtl_regions.py — compare QTL *regions* across mapping models.

Each model reports one row per QTL (left flanking marker, peak, right flanking
marker).  Two models found "the same QTL" when their regions overlap on the same
chromosome, not when they happen to name the same peak marker — peaks of
different methods routinely differ by a marker or two.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

NUM_COLS = [
    "Left_Position", "Peak_Position", "Right_Position", "LOD", "P_Value", "PVE_Marker_Percent", "Additive_Effect",
    "Interval_Size_cM", "CI_Start_cM", "CI_End_cM",
]


_CHR_TOKEN = re.compile(r"^(?:chromosome|chr|ch|lg|c|s)", re.I)


def physical_mb(marker, chromosome):
    """Physical position (Mb) read from a marker NAME: the number that follows the token naming the
    chromosome (chr02_16919033, IRGSP1_C02_34852782, MSU7_1_4797026_C-T, S2_1234567). None when the name
    carries no position. Same rule as the linkage tool's marker ordering."""
    try:
        chrom = int(float(chromosome))
    except (TypeError, ValueError):
        return None
    tok = re.split(r"[_:.\-]", str(marker))
    for i in range(len(tok) - 1):
        t = _CHR_TOKEN.sub("", tok[i], count=1)
        if t.isdigit() and int(t) == chrom and re.fullmatch(r"\d{4,}", tok[i + 1]):
            return round(int(tok[i + 1]) / 1e6, 6)
    return None


def physical_position(marker):
    """(chromosome, bp) read from a marker NAME without knowing the chromosome in advance: a chromosome
    number followed by a position of four or more digits (chr02_16919033, IRGSP1_C02_34852782,
    MSU7_9_14567436_Indel). None for names that carry no position (NAS3_2, Pi33_2)."""
    tok = re.split(r"[_:.\-]", str(marker))
    for i in range(len(tok) - 1):
        t = _CHR_TOKEN.sub("", tok[i], count=1)
        if t.isdigit() and re.fullmatch(r"\d{4,}", tok[i + 1]):
            return int(t), int(tok[i + 1])
    return None


def order_check(left_mb, peak_mb, right_mb):
    """Does the map (cM) run in the same direction as the genome (Mb) between the three markers?"""
    if None in (left_mb, peak_mb, right_mb):
        return ""
    if left_mb <= peak_mb <= right_mb:
        return "map follows physical order"
    if left_mb >= peak_mb >= right_mb:
        return "map runs against physical order"
    return "map order differs from physical order"


def chr_key(c):
    try:
        return (0, float(c), "")
    except (TypeError, ValueError):
        return (1, 0.0, str(c))


def load_regions(model_results: dict) -> pd.DataFrame:
    """All models' QTL regions in one table (adds a Model column)."""
    frames = []
    for model, res in model_results.items():
        path = res.get("qtl_csv")
        if not path or not Path(path).exists():
            continue
        df = pd.read_csv(path)
        if df.empty or "Peak_Marker" not in df.columns:
            continue
        df = df.copy()
        for c in NUM_COLS:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors="coerce")
        df["Chromosome"] = df["Chromosome"].astype(str)
        df["Model"] = model
        frames.append(df)
    if not frames:
        return pd.DataFrame(columns=["Model", "Chromosome", *NUM_COLS])
    return pd.concat(frames, ignore_index=True)


def cluster_regions(regions: pd.DataFrame) -> list[dict]:
    """Group regions from different models that overlap on a chromosome.

    The strongest region seeds a cluster; a weaker region from another model joins
    it when the two intervals overlap.  Joining is judged against the seed only
    (never a growing union), and a model contributes at most one region per
    cluster, so two distinct neighbouring QTL of one model are not blurred together.
    """
    out = []
    if regions is None or regions.empty:
        return out
    for chrom in sorted(regions["Chromosome"].unique(), key=chr_key):
        sub = regions[regions["Chromosome"] == chrom].sort_values("LOD", ascending=False)
        clusters: list[dict] = []
        for _, r in sub.iterrows():
            best, best_d = None, None
            for c in clusters:
                if r["Model"] in c["models"]:
                    continue
                if r["Left_Position"] <= c["seed_right"] and r["Right_Position"] >= c["seed_left"]:
                    d = abs(r["Peak_Position"] - c["seed_peak"])
                    if best is None or d < best_d:
                        best, best_d = c, d
            if best is None:
                clusters.append({"seed_left": r["Left_Position"], "seed_right": r["Right_Position"],
                                 "seed_peak": r["Peak_Position"], "models": {r["Model"]: r}})
            else:
                best["models"][r["Model"]] = r

        rows = []
        for c in clusters:
            members = list(c["models"].values())
            top = max(members, key=lambda m: m["LOD"])
            left = min(members, key=lambda m: m["Left_Position"])
            right = max(members, key=lambda m: m["Right_Position"])
            lods = [m["LOD"] for m in members]
            pves = [m["PVE_Marker_Percent"] for m in members if pd.notna(m.get("PVE_Marker_Percent"))]
            row = {
                "Chromosome": chrom,
                "Left_Marker": left["Left_Marker"], "Left_Position": left["Left_Position"],
                "Peak_Marker": top["Peak_Marker"], "Peak_Position": top["Peak_Position"],
                "Right_Marker": right["Right_Marker"], "Right_Position": right["Right_Position"],
                "Region_Size_cM": right["Right_Position"] - left["Left_Position"],
                "N_Models_Detected": len(members),
                "Models_List": ", ".join(sorted(c["models"])),
                "Peak_Range_cM": f"{min(m['Peak_Position'] for m in members):g}-{max(m['Peak_Position'] for m in members):g}",
                "Mean_LOD": round(sum(lods) / len(lods), 3),
                "Max_LOD": round(max(lods), 3),
                "Mean_PVE": round(sum(pves) / len(pves), 3) if pves else None,
            }
            # physical positions from the marker names, and whether the map agrees with the genome
            row["Left_Mb"] = physical_mb(left["Left_Marker"], chrom)
            row["Peak_Mb"] = physical_mb(top["Peak_Marker"], chrom)
            row["Right_Mb"] = physical_mb(right["Right_Marker"], chrom)
            row["Order_Check"] = order_check(row["Left_Mb"], row["Peak_Mb"], row["Right_Mb"])
            for model, m in c["models"].items():
                row[f"LOD_{model}"] = round(float(m["LOD"]), 3)
            rows.append(row)
        rows.sort(key=lambda x: x["Peak_Position"])
        for i, row in enumerate(rows, 1):
            row["Region_ID"] = f"R_{chrom}.{i}"
        out.extend(rows)
    return out
