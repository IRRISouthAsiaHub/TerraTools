"""
GWAS runner for the TerraTools app: reads the genotypes, runs the analysis (analysis.py), draws the figures (plots.py) and writes the
tables (reporting.py). The statistics live in genotype_io, qc, population, mixed and inference.
"""

import os
import tempfile

import numpy as np
import pandas as pd

from tools.gwas.backend import analysis as an, genotype_io, inference as inf, plots, population, reporting
from tools.gwas.backend.genotype_io import InputError


# ─────────────────────── DATA PARSING (kept for the Et-GWAS tool) ───────────────────────────
def parse_data(file_obj=None, text_input=None, encoding="auto", extra_missing=()):
    """(DataFrame with ID, Phenotype and one dosage column per SNP, SNP names, chromosomes, positions).

    Genotypes are converted to 0/1/2 dosage with NaN for missing. The encoding is never guessed from a matrix of 0 and 1: pass
    encoding="dosage" or "centred" for such files."""
    path = file_obj.name if file_obj is not None else None
    data = genotype_io.read_gwas_file(path=path, text=text_input, encoding=encoding, extra_missing=extra_missing)
    head = pd.DataFrame({"ID": data.ids, "Phenotype": data.y})
    geno = pd.DataFrame(data.G, columns=list(data.snp_ids))
    return pd.concat([head, geno], axis=1), list(data.snp_ids), np.asarray(data.chrom), np.asarray(data.pos)


class _SimpleFileObj:
    """Tiny shim so parse_data() also works with a plain path string (it expects an object with a `.name`)."""
    def __init__(self, path):
        self.name = path


def _jsonable(o):
    """Plain Python objects only (the result is sent to the page as JSON)."""
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return [_jsonable(v) for v in o.tolist()]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, (float, np.floating)):
        o = float(o)
        return o if np.isfinite(o) else None
    return o


def _write_csv(path, rows, columns):
    pd.DataFrame(rows, columns=columns).to_csv(path, index=False, encoding="utf-8-sig")
    return path


class GWASEngine:
    """Runs one GWAS from a GWASConfig and returns the result dictionary shown by the app."""

    def __init__(self, config):
        self.config = config

    def run(self):
        cfg = self.config
        out_dir = tempfile.mkdtemp(prefix="terratools_gwas_")
        log = lambda m: print(m, flush=True)                                                  # noqa: E731
        try:
            from tools.gwas.backend.gwas_inputs import describe, prepare_input
            input_path, input_info = prepare_input(cfg, out_dir, log=log)
            data = an.read_genotypes(cfg, input_path, generated=bool(input_info))
            labels = None
            alleles = None
            if input_info:
                alleles = input_info.get("alleles")
                if input_info.get("subpopulation") is not None:
                    lab = dict(zip(input_info["sample_ids"], input_info["subpopulation"]))
                    labels = np.array([lab.get(i, "unknown") for i in data.ids])
            A = an.analyse(cfg, data, labels=labels, alleles=alleles, log=log)
            if input_info and input_info.get("source") == "VCF":      # what the reader dropped before the QC step, for the ledger
                dropped = int(input_info["variants_read"] - input_info["markers"])
                A.qc.ledger.insert(1, {"stage": "VCF reading", "samples_removed": 0, "snps_removed": dropped,
                                       "samples_left": int(data.G.shape[0]), "snps_left": int(data.G.shape[1]),
                                       "rule": f"{input_info['variants_read']:,} records read; biallelic SNPs with GT, SNP missing <= "
                                               f"{input_info['prefilter_missing']:.0%} and MAF >= {input_info['prefilter_maf']:g} kept"})
        except InputError as e:
            return {"status": "error", "message": str(e)}
        except Exception as e:
            import traceback
            traceback.print_exc()
            return {"status": "error", "message": str(e)}

        try:
            return self._finish(A, out_dir, input_path, input_info, describe)
        except Exception as e:
            import traceback
            traceback.print_exc()
            return {"status": "error", "message": f"The analysis finished but the report could not be written: {e}"}

    # ------------------------------------------------------------------------------------
    def _finish(self, A, out_dir, input_path, input_info, describe):
        import matplotlib.pyplot as plt
        cfg = A.cfg
        res = A.results
        summary = an.summary_rows(A)
        mta, dropped_pve = an.mta_rows(A)
        diag = reporting.diagnostics_rows(A)

        # ── figures ──
        items, circ_p, circ_thr = [], {}, {}
        for name, r in res.items():
            sig = r.info["sig"]
            ex = np.zeros(len(r.p), dtype=bool)
            thr = dict(sig["thresholds"])
            if "permutation_p" in thr and not thr["permutation_p"]:
                thr.pop("permutation_p")
            items.append((name, np.asarray(r.p), thr, sig["genome_wide"], ex & ~sig["genome_wide"]))
            circ_p[name] = np.asarray(r.p)
            circ_thr[name] = sig["thresholds"]["bonferroni_p"]
        figs = {}
        if len(items) == 1:
            name, p, thr, gw, ex = items[0]
            figs["01_Manhattan"] = plots.plot_manhattan(p, A.chrom, A.pos, name, thr, gw, ex, show_bonferroni=cfg.show_bonferroni,
                                                        show_suggestive=cfg.show_suggestive)
        else:
            figs["01_Manhattan"] = plots.plot_model_comparison_manhattan(items, A.chrom, A.pos)
        lam = {n: r.info["diag"]["lambda"] for n, r in res.items()}
        figs["02_QQ"] = plots.plot_qq({n: r.p for n, r in res.items()}, lam, "QQ plot")
        try:
            hits = [r for r in mta if r["Category"] in ("Genome-wide significant", "QTN")]
            figs["03_Circos"] = plots.plot_circos_mta(circ_p, A.chrom, A.pos, A.snp, list(res), hits, thresholds=circ_thr)
        except Exception as ex:
            import traceback
            traceback.print_exc()
            print(f"Circos plot failed: {ex}")
        grp = A.structure.get("groups", {})

        png_paths = {}
        for name, fig in figs.items():
            path = os.path.join(out_dir, f"{name}.png")
            fig.savefig(path, format="png", dpi=plots.RENDER_DPI, bbox_inches="tight", facecolor=fig.get_facecolor())
            png_paths[name] = path
        zip_path = plots.save_all_plots_zip(figs, out_dir)
        plt.close("all")

        # ── tables ──
        mta_csv = _write_csv(os.path.join(out_dir, "gwas_mtas.csv"), mta, an.MTA_COLUMNS)
        qc_csv = _write_csv(os.path.join(out_dir, "gwas_qc_stages.csv"), A.qc.ledger, reporting.QC_COLUMNS)
        pd.DataFrame(A.qc.snp_table).to_csv(os.path.join(out_dir, "gwas_qc_snps.csv"), index=False, encoding="utf-8-sig")
        pd.DataFrame(A.qc.sample_table).to_csv(os.path.join(out_dir, "gwas_qc_samples.csv"), index=False, encoding="utf-8-sig")
        excel_path = reporting.export_excel(A, summary, mta, os.path.join(out_dir, "gwas_results.xlsx"))

        # ── message ──
        any_hit = any(r.info["diag"]["n_genome_wide"] for r in res.values())
        message = inf.NO_HITS if not any_hit else ""
        lines = [f"Samples: {len(A.y):,} analysed of {A.data.G.shape[0]:,}  |  SNPs: {len(A.snp):,} analysed of {A.data.G.shape[1]:,}",
                 A.data.encoding["message"]]
        lines += [f"QC {s['stage']}: {s['snps_removed']:,} SNPs and {s['samples_removed']:,} samples removed" for s in A.qc.ledger[1:]
                  if s["snps_removed"] or s["samples_removed"]]
        lines.append(f"Kinship: {A.kinship_info['method']} ({A.kinship_info['n_snps']:,} SNPs)" + ("; leave-one-chromosome-out" if cfg.loco else ""))
        for row in summary:
            lines.append(f"{row['Model']}: {row['Result']}; λ = {row['Genomic_Lambda']}; Bonferroni p <= {row['Bonferroni_threshold_p']:.2e}")
        for name, msg in A.errors.items():
            lines.append(f"{name} did not run: {msg}")
        if message:
            lines.append(message)
        text = "\n".join(lines)
        if input_info:
            text = describe(input_info) + "\n\n" + text

        result = {
            "status": "success", "output_dir": out_dir, "plots": png_paths, "excel": excel_path, "mta_csv": mta_csv, "qc_csv": qc_csv,
            "gwas_input": input_path if input_info else None, "input_info": input_info, "zip": zip_path, "summary": text,
            "message": message, "model_table": summary, "mta_table": mta, "qc_table": A.qc.ledger, "diagnostics": diag,
            "errors": A.errors,
            "notes": [k + ": " + v for k, v in reporting.info_rows(A)[:6]] + list(A.notes),
            "mta_dropped_by_pve": dropped_pve,
            "structure": {k: v for k, v in A.structure.items() if k != "groups"} | (
                {"groups": {k: v for k, v in grp.items() if k != "assignment"}} if grp else {}),
        }
        return _jsonable(result)
