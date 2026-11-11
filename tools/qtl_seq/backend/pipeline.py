"""
QTL-seq pipeline: inputs (BAM or VCF) -> informative SNPs -> filters -> SNP-index / Delta -> simulation limits ->
sliding windows -> QTL intervals -> files, figures and report.

`run_analysis` is the one entry point the GUI uses; `analyze_snps` is the part that starts from a SNP table (it is
also what the validation mode and the tests call).
"""

from __future__ import annotations

import datetime as _dt
import json
import platform
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

from tools.common.seqtools import Fasta
from tools.qtl_seq import config as cfgmod
from tools.qtl_seq.backend import annotation, bcftools, plots, stats, vcfio
from tools.qtl_seq.backend.bamio import BamPipelineError, build_bam_snps, validate_bam_inputs
from tools.qtl_seq.backend.vcfio import VcfError

ASSUMPTIONS = [
    "SNP-index = reads of the index allele / (reads of both alleles). Delta = SNP-index(high bulk) - SNP-index(low bulk).",
    "Null model: each bulk is a random sample of plants of the chosen population type (F2: 2n alleles, RIL: n homozygous plants, "
    "backcross: n plants with the recurrent parent chosen in the settings); reads are sampled binomially from the bulk's allele share.",
    "The two bulks are independent. Residual heterozygosity of RILs, sequencing error and mapping error are not part of the simulation.",
    "Window limits: 'snp_average' averages the simulated per-SNP limits (each from its own pair of read depths) over the SNPs of the "
    "window; 'window' simulates the window as a whole with one allele share per bulk.",
    "The limits are not corrected for the large number of overlapping windows that are compared with them.",
    "Without both parents the markers are not validated as parental polymorphisms.",
]


class AnalysisError(Exception):
    pass


class Cancelled(Exception):
    pass


class Log:
    def __init__(self, path, progress=None):
        self.path = path
        self.fh = open(path, "a", encoding="utf-8")
        self.progress = progress
        self.pct = 0

    def __call__(self, msg, pct=None):
        if pct is not None:
            self.pct = pct
        self.fh.write(f"{_dt.datetime.now().strftime('%H:%M:%S')}  {msg}\n")
        self.fh.flush()
        if self.progress:
            self.progress(self.pct, msg)

    def close(self):
        self.fh.close()


# ============================================================
# SETTINGS
# ============================================================

def check_settings(p):
    err = []
    if p["population"] not in cfgmod.POPULATIONS:
        err.append("Population type must be F2, RIL or backcross (BC).")
    for k, label in (("n_high", "high bulk"), ("n_low", "low bulk")):
        try:
            if int(p[k]) < 1:
                raise ValueError
        except Exception:
            err.append(f"Enter the number of plants in the {label} (a whole number of at least 1).")
    if p["window_size"] < 1 or p["window_step"] < 1:
        err.append("Window size and step must be positive.")
    elif p["window_step"] > p["window_size"]:
        err.append("The window step is larger than the window size, so SNPs between windows would be skipped.")
    if p["min_depth"] < 1:
        err.append("The minimum read depth must be at least 1.")
    if p["replicates"] < 100:
        err.append("Use at least 100 simulation replicates (10,000 recommended).")
    if p["depth_mode"] == "fixed" and int(p["fixed_depth"]) < 1:
        err.append("The fixed sequencing depth must be at least 1.")
    if p["qtl_min_level"] not in cfgmod.LEVELS:
        err.append("The QTL confidence level must be 90, 95 or 99.")
    return err


def coerce(cfg):
    p = cfgmod.merged(cfg)
    for k in ("n_high", "n_low", "min_depth", "max_depth", "min_parent_depth", "min_parent_gq", "min_mapq", "min_baseq",
              "min_alt_reads", "window_size", "window_step", "min_snps_per_window", "replicates", "seed", "fixed_depth",
              "qtl_min_level", "merge_gap"):
        try:
            p[k] = int(float(p[k]))
        except (TypeError, ValueError):
            pass
    for k in ("min_snp_index", "min_qual", "base_error", "min_alt_fraction"):
        try:
            p[k] = float(p[k])
        except (TypeError, ValueError):
            pass
    p["population"] = {"BACKCROSS": "BC"}.get(str(p["population"]).upper(), str(p["population"]).upper())
    return p


# ============================================================
# ANALYSIS OF A SNP TABLE
# ============================================================

def analyze_snps(snps, chrom_lengths, p, log=lambda *a, **k: None, cancelled=None):
    """Filters, indices, simulation, windows and QTLs for a SNP table (columns of stats.SNP_COLUMNS).

    Returns a dict: kept (SNP table with indices), ledger, limits (LimitTable), windows, qtls, sig_windows."""
    kept, ledger = stats.apply_filters(snps, p["min_depth"], p["max_depth"], p["min_snp_index"], p["apply_snp_index_filter"])
    log(f"{ledger['retained']:,} of {ledger['input']:,} SNPs pass the depth and SNP-index filters "
        f"(low depth {ledger['removed_low_depth']:,}, high depth {ledger['removed_high_depth']:,}, "
        f"SNP-index {ledger['removed_snp_index']:,})", 58)
    if kept.empty:
        raise AnalysisError("No SNP is left after the depth and SNP-index filters. Lower the minimum depth or SNP-index threshold, "
                            "or check the bulks and parents.")
    # chromosome lengths: fall back to the largest SNP position
    maxpos = kept.groupby("Chromosome")["Position"].max()
    lengths = {c: int(max(chrom_lengths.get(c, 0), maxpos[c])) for c in maxpos.index}
    order = [c for c in chrom_lengths if c in lengths] + [c for c in lengths if c not in chrom_lengths]
    lengths = {c: lengths[c] for c in order}

    log("Simulating the null distribution of the Delta SNP-index", 60)
    sim_args = (p["population"], p["n_high"], p["n_low"], p["replicates"], p["seed"], p.get("bc_recurrent", "non_index"))
    if p["depth_mode"] == "fixed":
        limits = stats.fixed_depth_limits(p["fixed_depth"], *sim_args, cancelled=cancelled)
        dh = dl = np.full(len(kept), int(p["fixed_depth"]))
    else:
        limits = stats.simulate_limits(kept["Depth_High"].to_numpy(), kept["Depth_Low"].to_numpy(), *sim_args,
                                       progress=lambda f: log("Simulating the null distribution of the Delta SNP-index", 60 + 15 * f),
                                       cancelled=cancelled)
        dh, dl = kept["Depth_High"].to_numpy(), kept["Depth_Low"].to_numpy()
    lim = limits.limits_for(dh, dl)
    for c, arr in lim.items():
        kept[c] = arr
    log(f"Limits simulated for {limits.meta['pairs']:,} pair(s) of read depths ({p['replicates']:,} replicates each, seed {p['seed']})", 76)

    log("Sliding-window averages", 78)
    kept = kept.sort_values(["Chromosome", "Position"], kind="stable").reset_index(drop=True)
    windows = stats.build_windows(kept, lengths, p, limits if p["ci_method"] == "snp_average" else None, cancelled=cancelled)
    log(f"{len(windows):,} windows with at least {p['min_snps_per_window']} SNPs", 82)
    positions = {c: g["Position"].to_numpy() for c, g in kept.groupby("Chromosome")}
    qtls, sigw = stats.call_qtls(windows, positions, p["window_step"], p["qtl_min_level"], p["merge_gap"], p["interval_basis"]) \
        if len(windows) else (pd.DataFrame(columns=stats.QTL_COLUMNS + ["QTL_ID"]), pd.DataFrame())
    counts = stats.level_counts(qtls)
    log(f"{len(qtls)} QTL interval(s): {counts['99']} at 99%, {counts['95']} at 95%, {counts['90']} at 90%", 84)
    return {"kept": kept, "ledger": ledger, "limits": limits, "windows": windows, "qtls": qtls, "sig_windows": sigw,
            "chrom_lengths": lengths, "counts": counts}


# ============================================================
# OUTPUT FILES
# ============================================================

QTL_FILE_COLUMNS = stats.QTL_COLUMNS + ["Direction", "QTL_ID", "Window_Count", "Window_Extent_Start", "Window_Extent_End",
                                        "Peak_Window_Start", "Peak_Window_End", "P_Threshold", "Interval_Basis"]
WINDOW_FILE_COLUMNS = ["Chromosome", "Window_Start", "Window_End", "SNP_Count", "Mean_SNP_index_High", "Mean_SNP_index_Low",
                       "Mean_Delta_SNP_index", "CI90_Upper", "CI90_Lower", "CI95_Upper", "CI95_Lower", "CI99_Upper",
                       "CI99_Lower", "Significant_90", "Significant_95", "Significant_99", "Mean_Depth_High", "Mean_Depth_Low"]


def _csv(df, path, **kw):
    df.to_csv(path, index=False, float_format="%.6g", **kw)


def write_outputs(out, res, p, log):
    kept, windows, qtls, limits = res["kept"], res["windows"], res["qtls"], res["limits"]
    files = {}

    def put(key, name):
        files[key] = str(out / name)
        return out / name

    q = qtls.reindex(columns=QTL_FILE_COLUMNS)
    _csv(q, put("qtls", "significant_QTLs.csv"))
    _csv(windows.reindex(columns=WINDOW_FILE_COLUMNS), put("windows", "all_sliding_windows.csv"))
    _csv(res["sig_windows"].drop(columns=[c for c in ("Sign",) if c in res["sig_windows"]]), put("sig_windows", "significant_windows.csv"))
    cols = ["Chromosome", "Position", "Ref", "Alt", "Ref_High", "Alt_High", "Depth_High", "SNP_index_High",
            "Ref_Low", "Alt_Low", "Depth_Low", "SNP_index_Low"]
    _csv(kept[cols], put("snp_index", "SNP_index_per_position.csv"))
    dcols = ["Chromosome", "Position", "Delta_SNP_index", "SNP_index_High", "SNP_index_Low", "Depth_High", "Depth_Low"] + \
        [f"CI{lv}_{s}" for lv in stats.LEVELS for s in ("Lower", "Upper")]
    d = kept[dcols].copy()
    d["Outside_CI95"] = (d["Delta_SNP_index"] > d["CI95_Upper"]) | (d["Delta_SNP_index"] < d["CI95_Lower"])
    _csv(d, put("delta", "delta_SNP_index.csv"))
    _csv(limits.frame, put("simulation", "simulation_confidence_intervals.csv"))
    vcfio.write_filtered_vcf(put("vcf", "filtered_informative_SNPs.vcf.gz"), kept, p["high_name"].replace(" ", "_"),
                             p["low_name"].replace(" ", "_"))
    log("Tables written", 87)
    chk = bcftools.index_and_check(files["vcf"], len(kept))
    res["vcf_check"] = chk
    log(("" if chk["ok"] else "WARNING: " if chk["available"] else "") + chk["message"], 88)
    if chk["indexed"]:
        files["vcf_index"] = files["vcf"] + ".tbi"
    return files


def make_figures(out, res, p, log):
    files = {}
    kept, win, qtls, L = res["kept"], res["windows"], res["qtls"], res["chrom_lengths"]
    import matplotlib.pyplot as plt

    def done(fig, key, stem):
        files[key] = plots.save_figure(fig, out / stem)
        plt.close(fig)

    done(plots.fig_genomewide(kept, win, qtls, L, p["high_name"], p["low_name"]), "genomewide", "QTLseq_genomewide_plot")
    log("Genome-wide plot", 90)
    done(plots.fig_bulk_index(kept, L, "High", p["high_name"]), "high_plot", "high_bulk_SNP_index")
    done(plots.fig_bulk_index(kept, L, "Low", p["low_name"]), "low_plot", "low_bulk_SNP_index")
    done(plots.fig_simulation(res["limits"], p["ci_method"]), "simulation_plot", "simulation_distribution")
    log("Bulk and simulation plots", 92)
    cdir = out / "chromosome_wise_QTL_plots"
    cdir.mkdir(exist_ok=True)
    chroms = [c for c in L if c in set(kept["Chromosome"])][:60]
    files["chrom_plots"] = {}
    for i, c in enumerate(chroms):
        fig = plots.fig_chromosome(c, kept, win, qtls, L, p["high_name"], p["low_name"])
        files["chrom_plots"][c] = plots.save_figure(fig, cdir / f"QTLseq_{plots.safe_name(c)}")
        plt.close(fig)
        log(f"Chromosome plot {c}", 92 + 5 * (i + 1) / max(len(chroms), 1))
    return files


# ============================================================
# UI DATA
# ============================================================

def _spread(n, cap):
    if n <= cap:
        return np.arange(n)
    return np.unique(np.linspace(0, n - 1, cap).astype(int))


def plot_payload(res, max_snps=45000, max_windows=1400):
    """Compact data for the interactive plots of the page: thinned SNPs and windows per chromosome."""
    kept, win = res["kept"], res["windows"]
    total = max(len(kept), 1)
    chroms = []
    for c, length in res["chrom_lengths"].items():
        s = kept[kept["Chromosome"] == c]
        w = win[win["Chromosome"] == c]
        if s.empty:
            continue
        si = _spread(len(s), max(1500, int(max_snps * len(s) / total)))
        wi = _spread(len(w), max_windows)
        s, w = s.iloc[si], w.iloc[wi]
        r3 = lambda a: [round(float(x), 3) for x in a]
        chroms.append({
            "name": c, "length": int(length), "n_snps": int((kept["Chromosome"] == c).sum()),
            "snp": {"pos": s["Position"].astype(int).tolist(), "high": r3(s["SNP_index_High"]), "low": r3(s["SNP_index_Low"]),
                    "delta": r3(s["Delta_SNP_index"])},
            "win": {"mid": w["Window_Mid"].astype(int).tolist(), "delta": r3(w["Mean_Delta_SNP_index"]),
                    "high": r3(w["Mean_SNP_index_High"]), "low": r3(w["Mean_SNP_index_Low"]),
                    **{f"ci{lv}{s_}": r3(w[f"CI{lv}_{side}"]) for lv in stats.LEVELS for s_, side in (("u", "Upper"), ("l", "Lower"))},
                    "n": w["SNP_Count"].astype(int).tolist()},
        })
    lim = res["limits"].frame
    eq = lim[lim["Depth_High"] == lim["Depth_Low"]].sort_values("Depth_High")
    use = eq if len(eq) >= 2 else lim.assign(_m=np.minimum(lim["Depth_High"], lim["Depth_Low"])).sort_values("_m")
    xs = use["Depth_High"] if len(eq) >= 2 else use["_m"]
    sim = {"depth": xs.astype(int).tolist(), **{f"ci{lv}{s_}": r3(use[f"CI{lv}_{side}"]) for lv in stats.LEVELS
                                                for s_, side in (("u", "Upper"), ("l", "Lower"))}}
    return {"chromosomes": chroms, "limits_by_depth": sim, "meta": res["limits"].meta,
            "null": stats.null_histogram(res["limits"])}


def chromosome_payload(res, chrom, max_snps=20000):
    kept, win = res["kept"], res["windows"]
    s = kept[kept["Chromosome"] == chrom]
    w = win[win["Chromosome"] == chrom]
    if s.empty:
        return {"status": "error", "message": f"No data for {chrom}."}
    s = s.iloc[_spread(len(s), max_snps)]
    r3 = lambda a: [round(float(x), 3) for x in a]
    q = res["qtls"][res["qtls"]["Chromosome"] == chrom] if len(res["qtls"]) else res["qtls"]
    return {"status": "success", "name": chrom, "length": int(res["chrom_lengths"][chrom]),
            "snp": {"pos": s["Position"].astype(int).tolist(), "high": r3(s["SNP_index_High"]), "low": r3(s["SNP_index_Low"]),
                    "delta": r3(s["Delta_SNP_index"])},
            "win": {"mid": w["Window_Mid"].astype(int).tolist(), "delta": r3(w["Mean_Delta_SNP_index"]),
                    "high": r3(w["Mean_SNP_index_High"]), "low": r3(w["Mean_SNP_index_Low"]),
                    **{f"ci{lv}{s_}": r3(w[f"CI{lv}_{side}"]) for lv in stats.LEVELS for s_, side in (("u", "Upper"), ("l", "Lower"))},
                    "n": w["SNP_Count"].astype(int).tolist()},
            "qtls": json.loads(q.to_json(orient="records")) if len(q) else []}


def _records(df):
    return json.loads(df.to_json(orient="records", double_precision=4)) if len(df) else []


# ============================================================
# ENTRY POINT
# ============================================================

def new_output_dir(root, name="QTLseq_results"):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    out = root / name
    n = 1
    while out.exists():
        n += 1
        out = root / f"{name}_{n}"
    out.mkdir()
    return out


def run_analysis(cfg, results_root, progress=None, cancelled=None):
    """Runs everything. Returns (result dict for the GUI, in-memory results or None)."""
    p = coerce(cfg)
    t0 = time.time()
    out_root = p.get("output_dir") or results_root
    problems = check_settings(p)
    if problems:
        return {"status": "error", "message": "\n".join(problems)}, None
    try:
        out = new_output_dir(out_root)
    except OSError as e:
        return {"status": "error", "message": f"The output folder cannot be created: {e}"}, None
    log = Log(out / "analysis.log", progress)
    log(f"Output folder: {out}", 1)
    fasta = None
    try:
        if p["mode"] == "vcf":
            snps, ledger0, lengths, warnings, info = vcfio.build_vcf_snps(p, log=log, cancelled=cancelled)
            if p.get("ref_fasta"):
                fasta = Fasta(p["ref_fasta"])
                v = vcfio.validate_vcf_inputs(p, fasta)
                warnings += [w for w in v["warnings"] if w not in warnings]
                lengths = {**{c: l for c, l in fasta.contigs().items() if c in set(snps["Chromosome"])}, **lengths}
        else:
            if not p.get("ref_fasta"):
                raise AnalysisError("Select the reference genome FASTA.")
            fasta = Fasta(p["ref_fasta"])
            snps, ledger0, lengths, warnings, info = build_bam_snps(p, fasta, log=log, cancelled=cancelled)
        log(f"{len(snps):,} informative SNPs ({ledger0['mode'].replace('_', ' ')})", 56)
        for w in warnings:
            log("WARNING: " + w)
        res = analyze_snps(snps, lengths, p, log=log, cancelled=cancelled)
        led = {**ledger0, **res["ledger"]}
        res["ledger"] = led
        files = write_outputs(out, res, p, log)
        vc = res["vcf_check"]
        if vc["available"] and not vc["ok"]:
            warnings.append(vc["message"])
        files.update(make_figures(out, res, p, log))

        genes = None
        if p.get("annotation"):
            try:
                g = annotation.load_genes(p["annotation"])
                if g.empty:
                    warnings.append("The annotation file has no gene features (a GFF3/GTF needs lines of type 'gene'), so no "
                                    "candidate genes could be listed.")
                genes = annotation.candidate_genes(res["qtls"], g)
                _csv(genes, out / "candidate_genes.csv")
                files["genes"] = str(out / "candidate_genes.csv")
                log(f"{len(genes):,} gene(s) overlap the QTL intervals ({len(g):,} genes read)", 97)
            except Exception as e:
                warnings.append(f"The annotation file could not be used: {e}")
                log("WARNING: " + warnings[-1])
        marker_text = {"two_parents": "SNPs homozygous for different alleles in both parents",
                       "one_parent": "SNPs from one parent only (not fully validated)",
                       "no_parent": "no parents (bulk data only)"}[led["mode"]]
        index_text = ("the VCF ALT allele (no parent)" if led["mode"] == "no_parent" else
                      f"{p['parent2_name'] if p['index_parent'] == 'parent2' else p['parent1_name']}")
        ctx = {"params": p, "ledger": led, "warnings": warnings, "snps": res["kept"], "windows": res["windows"],
               "qtls": res["qtls"], "genes": genes, "chrom_lengths": res["chrom_lengths"], "limits": res["limits"],
               "n_windows": len(res["windows"]), "marker_mode_text": marker_text, "index_text": index_text}
        plots.write_report(out / "QTLseq_analysis_report.pdf", ctx)
        files["report"] = str(out / "QTLseq_analysis_report.pdf")
        log("Report written", 99)

        record = {k: v for k, v in p.items()}
        record.update({"finished": _dt.datetime.now().isoformat(timespec="seconds"), "seconds": round(time.time() - t0, 1),
                       "python": platform.python_version(), "platform": platform.platform(), "marker_mode": led["mode"],
                       "ledger": led, "warnings": warnings, "input_info": json.loads(json.dumps(info, default=str)),
                       "simulation": res["limits"].meta, "assumptions": ASSUMPTIONS, "vcf_check": res.get("vcf_check")})
        (out / "analysis_parameters.json").write_text(json.dumps(record, indent=2, default=str), encoding="utf-8")
        files["params"], files["log"] = str(out / "analysis_parameters.json"), str(out / "analysis.log")
        log(f"Done in {time.time() - t0:.0f} s", 100)
        result = {
            "status": "success", "out_dir": str(out), "files": files, "warnings": warnings, "ledger": led,
            "summary": {"informative_snps": int(len(res["kept"])), "qtls": int(len(res["qtls"])),
                        "qtl_99": res["counts"]["99"], "qtl_95": res["counts"]["95"], "qtl_90": res["counts"]["90"],
                        "windows": int(len(res["windows"])), "chromosomes": len(res["chrom_lengths"]),
                        "mode": led["mode"], "marker_text": marker_text, "index_text": index_text,
                        "population": p["population"], "n_high": p["n_high"], "n_low": p["n_low"],
                        "high_name": p["high_name"], "low_name": p["low_name"], "ci_method": p["ci_method"],
                        "replicates": p["replicates"], "seed": p["seed"], "seconds": round(time.time() - t0, 1),
                        "vcf_check": vc["message"], "engine": info.get("engine", "")},
            "qtls": _records(res["qtls"]), "genes": _records(genes) if genes is not None else [],
            "plot": plot_payload(res), "assumptions": ASSUMPTIONS,
            "params": {k: p[k] for k in ("window_size", "window_step", "min_snps_per_window", "min_depth", "max_depth",
                                         "min_snp_index", "qtl_min_level", "merge_gap", "interval_basis", "depth_mode",
                                         "fixed_depth")},
        }
        return result, res
    except InterruptedError:
        log("Cancelled by the user")
        return {"status": "cancelled", "message": "Cancelled."}, None
    except (AnalysisError, VcfError, BamPipelineError) as e:
        log("ERROR: " + str(e))
        return {"status": "error", "message": str(e)}, None
    except Exception as e:                                           # a bug or an unexpected file: keep the details in the log
        log("ERROR: " + traceback.format_exc())
        return {"status": "error", "message": f"The analysis failed: {e}. Details are in {out / 'analysis.log'}."}, None
    finally:
        log.close()
        if fasta is not None:
            fasta.close()


def check_inputs(cfg):
    """The 'Check inputs' button: validation without running the analysis."""
    p = coerce(cfg)
    errors = check_settings(p)
    warnings, info = [], {}
    fasta = None
    try:
        if p.get("ref_fasta") and Path(p["ref_fasta"]).is_file():
            try:
                fasta = Fasta(p["ref_fasta"])
            except Exception as e:
                errors.append(f"The reference FASTA cannot be read: {e}")
        if p["mode"] == "vcf":
            v = vcfio.validate_vcf_inputs(p, fasta)
        else:
            v = validate_bam_inputs(p, fasta)
        errors += v["errors"]
        warnings += v["warnings"]
        info = v["files"]
        if p["mode"] == "vcf" and p.get("parent_vcf"):
            try:
                h = vcfio.read_header(p["parent_vcf"])
                s1, s2 = vcfio.parent_samples(h, p)
                info["parents"] = {"parent1_sample": s1, "parent2_sample": s2, "samples": h["samples"]}
            except VcfError as e:
                errors.append(str(e))
        have = [bool(p.get(k)) for k in (("parent_vcf",) if p["mode"] == "vcf" else ("parent1_bam", "parent2_bam"))]
        if not any(have):
            warnings.append("No parental data: parental polymorphism cannot be validated and the analysis is a bulk-only "
                            "allele-frequency comparison.")
        elif p["mode"] == "bam" and sum(have) == 1:
            warnings.append("Only one parent was supplied: parental polymorphism cannot be completely validated.")
    except Exception as e:
        errors.append(str(e))
    finally:
        if fasta is not None:
            fasta.close()
    return {"status": "success", "ok": not errors, "errors": errors, "warnings": warnings, "files": info}
