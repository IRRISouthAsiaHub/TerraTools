"""
linkage_qtl_export.py — Build the genotype file QTL Mapping reads, from a
finished Linkage Mapping run.

The linkage map (Marker, Chromosome, Position) is matched by marker ID to the
same genotype file that went into Linkage Mapping, keeping the linkage order,
and the result is transposed:

    Marker      | m1   | m2   | m3  ...      <- row 1: marker names
    Chromosome  | 1    | 1    | 2   ...      <- row 2: linkage group
    Position    | 0.0  | 12.4 | 0.0 ...      <- row 3: map position (cM)
    Ind_1       | 0    | 1    | 2   ...      <- one row per individual,
    Ind_2       | 2    | 0    | -1  ...         codes 0 / 1 / 2, -1 = missing

That is exactly the layout of the QTL Mapping "Genotype CSV".

Genotypes are the ORIGINAL input calls (0/1/2/-1; or a/ab/b/- converted to
0/1/2/-1), not the filled-in values used while building the map, so missing
calls stay missing.
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd

# One sheet of an .xlsx can hold 16,384 columns; the CSV has no such limit.
XLSX_MAX_COLS = 16384
XLSX_MAX_CELLS = 4_000_000

TEXT_CODES = {"a": 0, "ab": 1, "b": 2, "-": -1, "": -1, "na": -1, "nan": -1}


def _read_genotypes(path):
    path = Path(path)
    if path.suffix.lower() in (".csv", ".txt", ".tsv"):
        df = pd.read_csv(path, sep="\t" if path.suffix.lower() == ".tsv" else ",", dtype=object)
    else:
        df = pd.read_excel(path, sheet_name=0, dtype=object)   # same sheet Linkage Mapping reads
    df.columns = [str(c).strip() for c in df.columns]
    return df


def _to_codes(col):
    """Column of raw calls -> integer codes 0/1/2, missing = -1."""
    num = pd.to_numeric(col, errors="coerce")
    text = col.astype(str).str.strip().str.lower().map(TEXT_CODES)
    codes = num.where(num.isin([-1, 0, 1, 2]))
    codes = codes.fillna(text).fillna(-1)
    return codes.astype(int)


def _chr_label(value):
    """'1', 1, 1.0 -> '1'; non-numeric names unchanged."""
    s = str(value).strip()
    try:
        f = float(s)
        return str(int(f)) if f == int(f) else s
    except ValueError:
        return s


def _chr_sort_key(label):
    try:
        return (0, float(label), "")
    except ValueError:
        return (1, 0.0, label)


def _safe_stem(path):
    return re.sub(r"[^\w\-]+", "_", Path(path).stem).strip("_") or "linkage"


def export_qtl_input(input_file, linkage_xlsx, out_dir):
    """Write <name>_for_QTL.csv / .xlsx into out_dir. Returns a summary dict."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    linkage = pd.read_excel(linkage_xlsx, sheet_name="Linkage Map", dtype={"Marker": str})
    linkage = linkage.dropna(subset=["Marker", "Position"]).copy()
    linkage["Chromosome"] = linkage["Chromosome"].map(_chr_label)

    raw = _read_genotypes(input_file)
    marker_col, chrom_col = "Marker", "Chromosome"
    if marker_col not in raw.columns:
        raise ValueError("The input file has no 'Marker' column.")
    geno_cols = [c for c in raw.columns if c not in (marker_col, chrom_col)]

    codes = pd.DataFrame({c: _to_codes(raw[c]) for c in geno_cols})
    # Columns with no real genotype call (status/notes columns) aren't individuals
    individuals = [c for c in geno_cols if (codes[c] != -1).any()]
    codes = codes[individuals]
    codes.insert(0, "__key", raw[marker_col].astype(str).str.strip())

    n_dup = int(codes["__key"].duplicated().sum())
    codes = codes.drop_duplicates("__key", keep="first")

    linkage["__key"] = linkage["Marker"].astype(str).str.strip()
    linkage["__order"] = np.arange(len(linkage))
    linkage["__chr_key"] = linkage["Chromosome"].map(_chr_sort_key)

    merged = linkage.merge(codes, on="__key", how="inner", sort=False)
    if merged.empty:
        raise ValueError("No marker IDs in the linkage map match the genotype file.")
    # Natural chromosome order (1, 2, ... 10, 11), then map position
    merged = merged.sort_values(["__chr_key", "Position", "__order"], kind="stable")

    markers = merged["Marker"].astype(str).tolist()
    chroms = merged["Chromosome"].tolist()
    positions = merged["Position"].astype(float).tolist()
    matrix = merged[individuals].to_numpy(dtype=int)        # markers x individuals

    # ── audit: what didn't make it, and why ───────────────────────────────
    codes_keys = set(codes["__key"])
    unmatched_linkage = linkage.loc[~linkage["__key"].isin(set(merged["__key"])),
                                    ["Marker", "Chromosome", "Position"]]
    in_linkage = set(linkage["__key"])
    unmatched_input = pd.DataFrame({
        "Marker": raw[marker_col].astype(str).str.strip(),
        "Chromosome": raw[chrom_col].map(_chr_label) if chrom_col in raw.columns else "",
    })
    unmatched_input = unmatched_input[~unmatched_input["Marker"].isin(in_linkage)].drop_duplicates("Marker")
    reasons = {}
    try:
        log = pd.read_excel(linkage_xlsx, sheet_name="Removal Log", dtype=str)
        reasons = {str(m).strip(): (r, d) for m, r, d in
                   zip(log["Marker"], log["Reason"], log.get("Detail", [""] * len(log)))}
    except Exception:
        pass
    unmatched_input["Reason"] = unmatched_input["Marker"].map(lambda m: reasons.get(m, ("Not in linkage map", ""))[0])
    unmatched_input["Detail"] = unmatched_input["Marker"].map(lambda m: reasons.get(m, ("", ""))[1])
    reason_counts = unmatched_input["Reason"].value_counts().to_dict()

    # ── the QTL file: transposed, three header rows ───────────────────────
    grid = [["Marker", *markers], ["Chromosome", *chroms], ["Position", *positions]]
    for j, ind in enumerate(individuals):
        grid.append([ind, *matrix[:, j].tolist()])
    qtl = pd.DataFrame(grid)

    stem = _safe_stem(input_file)
    csv_path = out_dir / f"{stem}_for_QTL.csv"
    xlsx_path = out_dir / f"{stem}_for_QTL.xlsx"
    qtl.to_csv(csv_path, header=False, index=False)

    merged_view = merged[["Marker", "Chromosome", "Position", *individuals]].reset_index(drop=True)
    notes = []
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as w:
        if qtl.shape[1] <= XLSX_MAX_COLS and qtl.size <= XLSX_MAX_CELLS:
            qtl.to_excel(w, sheet_name="QTL_Input", header=False, index=False)
        else:
            notes.append("QTL_Input sheet skipped (too large for Excel) — use the CSV.")
        if merged_view.shape[1] <= XLSX_MAX_COLS and merged_view.size <= XLSX_MAX_CELLS:
            merged_view.to_excel(w, sheet_name="Linkage_Merged", index=False)
        else:
            notes.append("Linkage_Merged sheet skipped (too large for Excel).")
        unmatched_linkage.to_excel(w, sheet_name="Unmatched_Linkage_IDs", index=False)
        unmatched_input.to_excel(w, sheet_name="Unmatched_Input_IDs", index=False)
        pd.DataFrame([
            ("Input genotype file", Path(input_file).name),
            ("Linkage markers", len(linkage)),
            ("Input markers", int(raw[marker_col].nunique())),
            ("Matched markers (in the QTL file)", len(markers)),
            ("Linkage markers not found in input", len(unmatched_linkage)),
            ("Input markers not in the linkage map", len(unmatched_input)),
            ("Individuals", len(individuals)),
            ("Duplicate input marker IDs dropped", n_dup),
            ("Genotype codes", "0 = hom. P1, 1 = heterozygous, 2 = hom. P2, -1 = missing"),
            *[(f"  not in map: {k}", v) for k, v in reason_counts.items()],
            *[("Note", n) for n in notes],
        ], columns=["Item", "Value"]).to_excel(w, sheet_name="Summary", index=False)

    n_prev_m, n_prev_i = min(12, len(markers)), min(8, len(individuals))
    return {
        "qtl_input_csv": str(csv_path).replace("\\", "/"),
        "qtl_input_xlsx": str(xlsx_path).replace("\\", "/"),
        "qtl_export": {
            "markers": len(markers),
            "individuals": len(individuals),
            "chromosomes": len(set(chroms)),
            "linkage_markers": len(linkage),
            "unmatched_linkage": len(unmatched_linkage),
            "unmatched_input": len(unmatched_input),
            "unmatched_input_reasons": {str(k): int(v) for k, v in reason_counts.items()},
            "duplicates_dropped": n_dup,
            "missing_pct": round(100 * float((matrix == -1).mean()), 2),
            "notes": notes,
            "preview": {
                "markers": markers[:n_prev_m],
                "chromosome": chroms[:n_prev_m],
                "position": [round(p, 3) for p in positions[:n_prev_m]],
                "rows": [{"id": individuals[j], "codes": matrix[:n_prev_m, j].tolist()} for j in range(n_prev_i)],
            },
        },
    }
