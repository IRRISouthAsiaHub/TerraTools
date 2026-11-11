"""
ocs_usefulness.py — optimal cross selection on the USEFULNESS of every cross (marker mode).

Inputs
  allele dosage matrix   genotypes x markers, 0 / 1 / 2 (inbred parents: 0 and 2; a heterozygous call is set to missing)
  marker effects         Marker, Effect (additive)
  genetic map            Marker, Chromosome, cM  - or, when there is no map, a VCF: its base-pair positions are turned
                         into cM with a constant recombination rate (cM per Mb) given by the user
  relationship matrix    always calculated from the dosage matrix (VanRaden 2008); nothing to upload

What is calculated
  1. the genomic value of every candidate  = dosage x marker effects
  2. the relationship (kinship) matrix      = VanRaden G from the dosage matrix, unless one is uploaded
  3. for every cross of the best candidates, SimpleMating::getUsefA() (ocs_uc.R) gives the mean, the progeny
     variance (from the effects and the recombination map, Lehermeier et al. 2017) and the usefulness
         UC = mean + i * sd           i = selection intensity of the selected proportion
  4. the plan of crosses that maximises the mean usefulness while keeping the coancestry of the parents and the
     inbreeding of the progeny down (ocs_optimizer.py).
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from tools.common.r_locator import find_rscript, get_r_env
from tools.common.subprocess_utils import no_window_kwargs

SCRIPT = Path(__file__).resolve().parent / "ocs_uc.R"
DEFAULT_POOL = 60


class UsefulnessError(Exception):
    pass


# ------------------------------------------------------------------ reading the files

def _clean(s):
    return str(s).strip()


def read_effects(path):
    df = pd.read_csv(path, sep=None, engine="python")
    if df.shape[1] < 2:
        raise UsefulnessError("The marker effects file needs two columns: Marker and Effect.")
    names = df.iloc[:, 0].map(_clean)
    lower = {str(c).strip().lower(): c for c in df.columns[1:]}
    col = next((lower[k] for k in ("effect", "effects", "additive", "add", "add_effect", "a", "beta") if k in lower), None)
    note = ""
    if col is None:
        num = [c for c in df.columns[1:] if pd.to_numeric(df[c], errors="coerce").notna().mean() > 0.9]
        if not num:
            raise UsefulnessError("No numeric effect column was found in the marker effects file (expected: Marker, Effect).")
        col = num[0]
        if len(num) > 1:
            note = f"The marker effects file has several numeric columns; '{col}' was used."
    eff = pd.Series(pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float), index=names)
    eff = eff[~eff.index.duplicated()].dropna()
    if eff.empty:
        raise UsefulnessError("The marker effects file has no usable effect values.")
    return eff, note


def read_map(path):
    df = pd.read_csv(path, sep=None, engine="python")
    if df.shape[1] < 3:
        raise UsefulnessError("The genetic map needs three columns: Marker, Chromosome and cM position.")
    low = {str(c).strip().lower(): c for c in df.columns}

    def pick(keys, fallback):
        for k in keys:
            if k in low:
                return low[k]
        return df.columns[fallback]

    mk = pick(("marker", "snp", "mkr", "id", "name", "marker_name"), 0)
    ch = pick(("chromosome", "chr", "chrom", "lg", "linkage_group", "linkage group"), 1)
    ps = pick(("cm", "position_cm", "pos_cm", "position", "pos", "map", "genetic_position"), 2)
    out = pd.DataFrame({"mkr": df[mk].map(_clean), "chr": df[ch].astype(str).str.strip(),
                        "pos": pd.to_numeric(df[ps], errors="coerce")})
    out = out.dropna(subset=["pos"]).drop_duplicates("mkr")
    if out.empty:
        raise UsefulnessError("The genetic map has no marker with a numeric cM position.")
    return out


def is_vcf(path, kind=None):
    return str(kind or "").lower() == "vcf" or str(path).lower().endswith((".vcf", ".vcf.gz"))


def read_map_vcf(path, marker_names, cm_per_mb):
    """Genetic map from the CHROM / POS / ID columns of a VCF: cM = bp / 1e6 x cM per Mb (one constant rate).
    A variant is matched by its ID, or by CHROM_POS / CHROM:POS when the ID is empty."""
    from tools.gene_extractor.backend.gff_genes import read_marker_positions

    try:
        rate = float(cm_per_mb)
    except (TypeError, ValueError):
        rate = 0.0
    if not rate > 0:
        raise UsefulnessError("Set the recombination rate (cM per Mb) used to turn the VCF positions into a genetic map.")
    names = {str(m).strip().lower(): str(m) for m in marker_names}
    try:
        pos = read_marker_positions(path, set(names))
    except Exception as e:
        raise UsefulnessError(f"The VCF could not be read: {e}")
    rows = [(names[k], str(chrom), bp / 1e6 * rate) for k, (chrom, bp) in pos.items() if k in names]
    if not rows:
        raise UsefulnessError("None of the marker names of the effects file is in the VCF (its ID column, or CHROM_POS). "
                              "The marker names must be the same as in the dosage matrix and the marker effects.")
    return pd.DataFrame(rows, columns=["mkr", "chr", "pos"]).drop_duplicates("mkr")


def read_dosage(path, marker_names):
    """Genotypes x markers dataframe of floats (rows = genotypes). A matrix with markers in the rows is turned round."""
    df = pd.read_csv(path, index_col=0, low_memory=False)
    df.index = df.index.map(_clean)
    df.columns = [_clean(c) for c in df.columns]
    names = set(marker_names)
    if len(names & set(df.index)) > len(names & set(df.columns)):
        df = df.T
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.loc[~df.index.duplicated()]
    df = df.loc[:, ~pd.Index(df.columns).duplicated()]
    return df


def _chrom_codes(chroms):
    """Chromosome labels -> 1, 2, 3 ... in natural order (chr02 and 2 are the same chromosome)."""
    import re

    def key(s):
        t = re.sub(r"^(chr|chromosome|lg|group)[_\-\s]*", "", s, flags=re.I)
        m = re.search(r"\d+", t)
        return (0, int(m.group()), t) if m else (1, 0, t)

    uniq = {}
    for s in chroms:
        uniq.setdefault(re.sub(r"^(chr|chromosome|lg|group)[_\-\s]*", "", s, flags=re.I).lstrip("0") or "0", s)
    order = sorted(uniq.values(), key=key)
    rank = {}
    for i, s in enumerate(order, 1):
        rank[re.sub(r"^(chr|chromosome|lg|group)[_\-\s]*", "", s, flags=re.I).lstrip("0") or "0"] = i
    return [rank[re.sub(r"^(chr|chromosome|lg|group)[_\-\s]*", "", s, flags=re.I).lstrip("0") or "0"] for s in chroms]


# ------------------------------------------------------------------ preparing the data

def prepare(geno_file, effects_file, map_file, map_kind=None, cm_per_mb=None):
    """Aligned dosage (0/2/NA), effects, map and relationship matrix, plus notes.
    map_file is a genetic map table (cM) or, with map_kind 'vcf' or a .vcf / .vcf.gz name, a VCF (positions in bp)."""
    notes = []
    eff, n = read_effects(effects_file)
    if n:
        notes.append(n)
    from_vcf = is_vcf(map_file, map_kind)
    if from_vcf:
        gmap = read_map_vcf(map_file, list(eff.index), cm_per_mb)
        notes.append(f"No genetic map was given: the cM positions were calculated from the base-pair positions of the VCF at "
                     f"{float(cm_per_mb):g} cM per Mb (one rate for the whole genome). Recombination is not even along a chromosome "
                     "(it is rare around the centromere), so the progeny variance is approximate; a real genetic map is more accurate.")
    else:
        gmap = read_map(map_file)
    M = read_dosage(geno_file, list(eff.index))

    common = [m for m in gmap["mkr"] if m in eff.index and m in M.columns]
    if len(common) < 10:
        raise UsefulnessError(
            f"Only {len(common)} marker(s) are in the dosage matrix, the marker effects and the genetic map at the same time. "
            "The marker names must be identical in the three files.")
    dropped = {"dosage": len(M.columns) - len(common), "effects": len(eff) - len(common), "map": len(gmap) - len(common)}
    if any(dropped.values()):
        notes.append("Markers that are not in all three files were left out "
                     f"(dosage {dropped['dosage']}, effects {dropped['effects']}, map {dropped['map']}); {len(common)} are used.")
    gmap = gmap.set_index("mkr").loc[common].reset_index()
    gmap["chr_code"] = _chrom_codes(gmap["chr"].tolist())
    gmap = gmap.sort_values(["chr_code", "pos"], kind="stable").reset_index(drop=True)
    order = gmap["mkr"].tolist()
    M = M[order]
    e = eff.loc[order].to_numpy(dtype=float)

    vals = M.to_numpy(dtype=float)
    lo, hi = np.nanmin(vals), np.nanmax(vals)
    if lo < 0 and hi <= 1:                       # -1 / 0 / 1 coding
        vals = vals + 1
        notes.append("The dosage matrix was coded -1 / 0 / 1; it was read as 0 / 1 / 2.")
        lo, hi = lo + 1, hi + 1
    if lo < 0 or hi > 2:
        raise UsefulnessError("The allele dosage matrix must hold 0, 1 and 2 (the number of copies of the allele).")
    vals = np.where(np.isfinite(vals), vals, np.nan)
    M = pd.DataFrame(vals, index=M.index, columns=M.columns)

    allna = M.columns[M.isna().all()]
    if len(allna):
        keep = [c for c in M.columns if c not in set(allna)]
        e = e[[order.index(c) for c in keep]]
        gmap = gmap[gmap["mkr"].isin(keep)].reset_index(drop=True)
        M = M[keep]
        notes.append(f"{len(allna)} marker(s) without any genotype call were left out.")

    n_het = int(np.nansum(M.to_numpy() == 1))
    miss = float(M.isna().mean().mean())
    if miss > 0.2:
        notes.append(f"{miss:.0%} of the dosage matrix is missing; check the file.")

    # relationship matrix: VanRaden (2008) on the full 0/1/2 data (missing calls take the marker mean)
    K = vanraden(M)
    source = "calculated from the dosage matrix (VanRaden)"
    ids = list(M.index)
    if len(ids) < 2:
        raise UsefulnessError("The dosage matrix needs at least 2 genotypes.")

    # heterozygous calls: SimpleMating takes inbred parents (0 / 2); a call of 1 is treated as missing
    M01 = M.copy()
    M02 = M.where(M != 1, np.nan)
    if n_het:
        notes.append(f"{n_het:,} heterozygous call(s) (dosage 1) were set to missing: the parents are taken as inbred lines "
                     "(0 / 2), which is what the usefulness model needs.")
    return {"M": M01, "M02": M02, "eff": e, "map": gmap, "K": K, "notes": notes, "kmat_source": source,
            "n_het": n_het, "missing": miss, "map_source": "VCF (cM from bp)" if from_vcf else "genetic map"}


def vanraden(M):
    """VanRaden (2008) genomic relationship matrix of a genotypes x markers 0/1/2 table."""
    X = M.to_numpy(dtype=float)
    p = np.nanmean(X, axis=0) / 2.0
    X = np.where(np.isnan(X), (2 * p)[None, :], X)
    Z = X - (2 * p)[None, :]
    denom = 2.0 * float(np.sum(p * (1 - p)))
    if denom <= 0:
        raise UsefulnessError("All markers are monomorphic: a relationship matrix cannot be computed.")
    G = Z @ Z.T / denom
    return pd.DataFrame(G, index=M.index, columns=M.index)


def genomic_values(M02, eff):
    """Dosage x effects with missing calls taking the marker mean (the genomic value of every candidate)."""
    X = M02.to_numpy(dtype=float)
    mu = np.nanmean(X, axis=0)
    X = np.where(np.isnan(X), mu[None, :], X)
    return X @ eff


# ------------------------------------------------------------------ usefulness of the crosses (R)

def usefulness_table(prep, pool, cfg, workdir):
    """SimpleMating::getUsefA for every pair of the `pool` genotypes."""
    rscript = find_rscript()
    if rscript is None:
        raise UsefulnessError("R was not found. The usefulness of the crosses is calculated with the R package SimpleMating.")
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    tag = tempfile.mkdtemp(prefix="ocs_uc_", dir=str(workdir))
    t = Path(tag)
    try:
        M = prep["M02"].loc[pool]
        M.to_csv(t / "geno.csv", index_label="Genotype")
        pd.DataFrame({"Marker": prep["map"]["mkr"], "Effect": prep["eff"]}).to_csv(t / "effects.csv", index=False)
        prep["map"][["chr_code", "pos", "mkr"]].rename(columns={"chr_code": "chr"}).to_csv(t / "map.csv", index=False)
        prep["K"].loc[pool, pool].to_csv(t / "kmat.csv", index_label="Genotype")
        ii, jj = np.triu_indices(len(pool), 1)
        pd.DataFrame({"Parent1": [pool[i] for i in ii], "Parent2": [pool[j] for j in jj]}).to_csv(t / "pairs.csv", index=False)
        with open(t / "cfg.json", "w", encoding="utf-8") as fh:
            json.dump({"dir": str(t).replace("\\", "/"), "prop_sel": float(cfg["prop_sel"]), "type": cfg["cross_type"],
                       "generation": int(cfg["generation"])}, fh)
        proc = subprocess.run([rscript, str(SCRIPT), str(t / "cfg.json")], capture_output=True, text=True, encoding="utf-8",
                              errors="replace", env=get_r_env(rscript), cwd=str(SCRIPT.parent), **no_window_kwargs())
        if proc.returncode == 3:
            raise UsefulnessError("The R package SimpleMating is not installed, and it is needed to calculate the usefulness of the crosses.")
        out = t / "usefulness.csv"
        if proc.returncode != 0 or not out.exists():
            tail = (proc.stderr or proc.stdout or "SimpleMating failed.").strip().splitlines()
            raise UsefulnessError("SimpleMating could not calculate the usefulness: " + (tail[-1] if tail else "unknown error"))
        df = pd.read_csv(out, dtype={"Parent1": str, "Parent2": str})
    finally:
        import shutil
        shutil.rmtree(tag, ignore_errors=True)
        try:
            workdir.rmdir()                      # only when empty
        except OSError:
            pass
    return df


# ------------------------------------------------------------------ public entry points

def _settings(config):
    cross_type = str(config.get("cross_type") or "DH").upper()
    if cross_type not in ("DH", "RIL"):
        raise UsefulnessError("The population type must be DH or RIL.")
    prop = float(config.get("prop_sel") if config.get("prop_sel") not in (None, "") else 0.05)
    if not 0 < prop < 1:
        raise UsefulnessError("The selected proportion must be between 0 and 1 (for example 0.05 for the best 5 %).")
    gen = int(config.get("generation") or 1)
    pool = int(config.get("pool_size") or DEFAULT_POOL)
    if pool < 2:
        raise UsefulnessError("The candidate pool needs at least 2 genotypes.")
    return {"cross_type": cross_type, "prop_sel": prop, "generation": gen, "pool_size": pool}


def preview(geno_file, effects_file, map_file, map_kind=None, cm_per_mb=None):
    """Counts shown on the screen before a run."""
    p = prepare(geno_file, effects_file, map_file, map_kind, cm_per_mb)
    return {"genotypes": int(p["M"].shape[0]), "markers": int(p["M"].shape[1]),
            "chromosomes": int(p["map"]["chr_code"].nunique()), "heterozygous_calls": p["n_het"],
            "missing": round(p["missing"], 4), "kinship": p["kmat_source"], "notes": p["notes"]}


def run_usefulness(config, results_dir):
    from tools.cross_selection.backend.ocs_optimizer import OCSError, plan_and_write

    for k, label in (("geno_file", "the allele dosage matrix"), ("effects_file", "the marker effects"),
                     ("map_file", "the genetic map (or the VCF)")):
        if not config.get(k):
            raise OCSError(f"Please select {label}.")
    try:
        st = _settings(config)
        prep = prepare(config["geno_file"], config["effects_file"], config["map_file"],
                       config.get("map_kind"), config.get("cm_per_mb"))
    except UsefulnessError as e:
        raise OCSError(str(e))

    nc = int(config.get("n_crosses") or 20)
    m = int(config.get("max_crosses_per_parent") or 4)
    angle = float(config.get("penalty_degree", 30))
    if nc < 1 or m < 1:
        raise OCSError("The number of crosses and the maximum crosses per parent must be at least 1.")
    if not 0 <= angle <= 90:
        raise OCSError("The penalty degree must be between 0 and 90.")

    ids = list(prep["M"].index)
    gv_all = genomic_values(prep["M02"], prep["eff"])
    pool_n = min(st["pool_size"], len(ids))
    pool_idx = np.argsort(-gv_all)[:pool_n]
    pool = [ids[i] for i in pool_idx]
    n = len(pool)
    if n * m < 2 * nc:
        raise OCSError(f"{n} candidates cannot make {nc} crosses with at most {m} crosses per parent "
                       f"(needs at least {-(-2 * nc // m)} candidates; raise the candidate pool or the limit per parent).")
    if nc > n * (n - 1) // 2:
        raise OCSError(f"There are only {n * (n - 1) // 2} different crosses among {n} candidates.")

    try:
        tab = usefulness_table(prep, pool, st, Path(results_dir) / "_work")
    except UsefulnessError as e:
        raise OCSError(str(e))

    pos = {g: i for i, g in enumerate(pool)}
    U = np.full((n, n), np.nan)
    Mean = np.zeros((n, n))
    SD = np.zeros((n, n))
    for r in tab.itertuples(index=False):
        i, j = pos[str(r.Parent1)], pos[str(r.Parent2)]
        U[i, j] = U[j, i] = r.Usefulness
        Mean[i, j] = Mean[j, i] = r.Mean
        SD[i, j] = SD[j, i] = r.sd
    if np.isnan(U[np.triu_indices(n, 1)]).any():
        raise OCSError("SimpleMating did not return the usefulness of every cross.")
    np.fill_diagonal(U, np.nan_to_num(np.nanmin(U)))

    K = prep["K"].loc[pool, pool].to_numpy(dtype=float)
    C = np.clip(0.5 * (K + K.T) / 2.0, 0.0, None)         # kinship (coancestry) scale
    a = gv_all[pool_idx]

    out = plan_and_write(pool, a, C, config, results_dir, U=U, pairinfo={"Mean": Mean, "sd": SD},
                         gain_label="mean usefulness of the crosses", crit_label="Genomic")

    results_dir = Path(results_dir)
    prep["K"].round(5).to_csv(
        results_dir / "ocs_relationship_matrix.csv", index_label="Genotype")
    tab = tab.rename(columns={"Mean": "Cross_Mean", "Variance": "Progeny_Variance", "sd": "Progeny_SD"})
    tab.sort_values("Usefulness", ascending=False).round(6).to_csv(results_dir / "ocs_all_crosses.csv", index=False)

    s = out["summary"]
    s.update({
        "mode": "markers", "n_markers": int(prep["M"].shape[1]), "n_genotypes": int(len(ids)),
        "pool_size": n, "cross_type": st["cross_type"], "generation": st["generation"], "prop_sel": st["prop_sel"],
        "kinship_source": prep["kmat_source"], "map_source": prep["map_source"],
        "engine": "SimpleMating getUsefA for the usefulness of the crosses; built-in simulated annealing for the plan",
    })
    notes = list(prep["notes"])
    if len(ids) > n:
        notes.append(f"The {n} genotypes with the highest genomic value (of {len(ids)}) were used as the candidates; "
                     f"all {n * (n - 1) // 2:,} crosses among them were evaluated.")
    if notes:
        s["note"] = " ".join(notes)
    out["usefulness_csv"] = str((results_dir / "ocs_all_crosses.csv").resolve()).replace("\\", "/")
    out["kinship_csv"] = str((results_dir / "ocs_relationship_matrix.csv").resolve()).replace("\\", "/")
    return out
