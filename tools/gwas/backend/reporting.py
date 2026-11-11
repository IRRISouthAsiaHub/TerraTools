"""
Result tables, the Notes sheet that defines every reported quantity, and the Excel workbook.
"""

from __future__ import annotations

import numpy as np
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

from tools.gwas.backend import analysis as an, inference as inf

EXCEL_MAX_ROWS = 1_048_575          # worksheet limit minus the header row

SUMMARY_COLUMNS = ["Model", "SNPs_Tested", "Genomic_Lambda", "Min_P", "Bonferroni_threshold_p", "Bonferroni_hits", "BH_hits",
                   "Permutation_threshold_p", "Permutation_hits", "Category_rule", "Genome_wide", "Suggestive", "Invalid_Statistics",
                   "Result"]
DIAG_COLUMNS = ["Model", "SNPs_Tested", "Invalid_Statistics", "Genomic_Lambda", "Min_P", "Genome_wide", "Suggestive", "Bonferroni_hits",
                "BH_hits", "Permutation_hits", "Convergence", "REML_check", "Heritability_null", "Variance_ratio_at_edge", "P_histogram_20bins"]
QC_COLUMNS = ["stage", "rule", "samples_removed", "snps_removed", "samples_left", "snps_left"]
PC_COLUMNS = ["PC", "Eigenvalue", "Variance_explained_%", "Cumulative_%"]
PC_COMPARE_COLUMNS = ["Model", "PCs", "Covariates", "Lambda", "Min_P", "Bonferroni_hits", "SNPs_Tested"]


def _f(v, nd=None):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(v):
        return None
    return round(v, nd) if nd is not None else v


def diagnostics_rows(A):
    rows = []
    for name, res in A.results.items():
        d = res.info["diag"]
        rows.append({"Model": name, "SNPs_Tested": d["snps_tested"], "Invalid_Statistics": d["invalid_statistics"],
                     "Genomic_Lambda": _f(d["lambda"], 3), "Min_P": _f(d["min_p"]), "Genome_wide": d["n_genome_wide"],
                     "Suggestive": d["n_suggestive"], "Bonferroni_hits": d["n_bonferroni"], "BH_hits": d["n_bh"],
                     "Permutation_hits": d["n_permutation"], "Convergence": d.get("convergence", "not applicable"),
                     "REML_check": d.get("reml_check", "not applicable"), "Heritability_null": _f(d.get("h2"), 3),
                     "Variance_ratio_at_edge": d.get("markers_with_variance_ratio_at_edge"),
                     "P_histogram_20bins": " ".join(str(x) for x in d["p_histogram"])})
    return rows


def pc_rows(A, n=20):
    ev = np.asarray(A.pc["eigenvalues"], float)[:n]
    ve = np.asarray(A.pc["variance_explained"], float)[:n] * 100
    return [{"PC": f"PC{i + 1}", "Eigenvalue": _f(ev[i], 5), "Variance_explained_%": _f(ve[i], 3), "Cumulative_%": _f(ve[:i + 1].sum(), 3)}
            for i in range(len(ev))]


def info_rows(A):
    cfg, qc = A.cfg, A.qc
    G = A.data
    n, m = len(A.y), len(A.snp)
    rows = []
    add = lambda k, v: rows.append((k, v))                                                    # noqa: E731
    add("Samples / SNPs", f"{G.G.shape[0]:,} samples and {G.G.shape[1]:,} SNPs read; {n:,} samples and {m:,} SNPs analysed after quality control.")
    add("Genotype encoding", G.encoding["message"] + (" The matrix was converted to dosage." if G.encoding.get("converted") else ""))
    add("Quality control", "Stages in order: duplicate positions, SNP call rate (<= %.0f%% missing), sample call rate (<= %.0f%% missing), "
                           "monomorphic SNPs, MAF >= %g. The QC sheet gives the numbers removed at every stage; the QC SNP and sample "
                           "tables list the status of every marker and sample. The raw genotype matrix is never overwritten."
        % (100 * cfg.max_snp_missing, 100 * cfg.max_sample_missing, cfg.maf))
    add("Missing genotypes", {"mean": "Each remaining missing call was replaced by the mean dosage of its marker (2 x allele frequency from the "
                                      "observed calls). This is mean imputation, not haplotype-based imputation: it ignores linkage "
                                      "disequilibrium and shrinks the variance of markers with many missing calls. %s calls were replaced."
                                      % f"{qc.imputed_calls:,}",
                              "remove": "Markers with any missing call were removed; no value was imputed.",
                              "none": "No imputation was applied; the genotypes were required to be complete (externally imputed)."}[cfg.imputation])
    k = A.kinship_info
    add("Relationship matrix", f"{k['method']} computed from {k['n_snps']:,} markers"
                               + (f" (LD-pruned from {k['n_snps_available']:,}: r² > {cfg.ld_prune_r2:g} within a sliding window)"
                                  if k["kinship_method"] == "ld_pruned" else "")
                               + f". Mean diagonal {k['mean_diagonal']:.3f}; smallest eigenvalue {k['eigen']['min_eigenvalue']:.2e} "
                                 f"({'positive semi-definite' if k['eigen']['psd'] else 'NOT positive semi-definite'}); rank {k['eigen']['rank']}. "
                                 "The matrix is built from the observed calls, before imputation, so the imputation choice does not "
                                 "enter it. Nothing is added to the diagonal and no correction is applied."
                               + (" Mixed models used leave-one-chromosome-out matrices: the chromosome of the tested marker is left out."
                                  if cfg.loco else ""))
    pc = A.pc
    add("Principal components", f"Leading eigenvectors of the relationship matrix; used only where a model asks for PCs. PC1 and PC2 "
                                f"explain {100 * A.structure['pc1_pc2_variance']:.1f}% of the relationship variance.")
    if "groups" in A.structure:
        g = A.structure["groups"]
        add("Subpopulation structure", f"Groups from {g['source']}: " + ", ".join(f"{a} ({b})" for a, b in g["sizes"].items())
            + ". Share of PC1 / PC2 / PC3 variance lying between groups: "
            + " / ".join(f"{100 * x:.0f}%" for x in g["pc_variance_between_groups"]) + ".")
    for name, res in A.results.items():
        add(f"Model: {name}", res.label)
        for note in res.notes:
            add(f"  {name} note", note)
        if res.permutation:
            pm = res.permutation
            t = pm["thresholds"][f"{float(cfg.alpha):g}"]
            add(f"  {name} simulated threshold",
                f"{pm['null'].replace('_', ' ')} under the model's own fitted null (covariates and covariance kept, no marker effect): "
                f"{pm['n_successful']} of {pm['n_requested']} replicates valid; threshold p = {t['p']:.3e} "
                f"(interval {t['p_stricter']:.2e} to {t['p_looser']:.2e}). " + " ".join(pm["limitations"]))
    for name, msg in A.errors.items():
        add(f"Model failed: {name}", msg)
    add("Significance", f"Bonferroni: p <= {cfg.alpha:g}/m, m = SNPs tested. BH: Benjamini-Hochberg adjusted p <= {cfg.alpha:g}. Simulated "
                        "threshold: see the model notes (single-locus models only). Only markers that pass the "
                        f"genome-wide rule ({cfg.significance_rule}) are listed as MTAs; the other criteria are reported beside them. The "
                        "suggestive threshold (p <= 1/m, one false positive expected per scan) is a reference line and a count in the "
                        "Summary, not a list of MTAs. When nothing passes the answer is: \"" + inf.NO_HITS + "\".")
    add("Multi-locus models", "FarmCPU and BLINK apply the chosen rule to the p-values of their final iteration. mrMLM and FASTmrMLM "
                              "report QTNs by their own criterion (LOD in the joint model); the p-value shown is the first-step scan "
                              "p-value, and the LOD is not converted into a p-value.")
    add("Effect / SE", "Change in the trait per copy of the counted allele (the allele whose copies are counted by the 0/1/2 dosage), "
                       "and its standard error from the model that produced the p-value.")
    add("LOD", "Single-locus models: the 1-df chi-square statistic implied by the p-value divided by 2 ln 10 (a LOD-scale "
               "re-expression of the same test). Multi-locus models: the model's own LOD.")
    add("PVE_%", "Share of the phenotypic variance, always as a percentage of the total variance. OLS: (RSS0 - RSS1)/SS_total, covariates "
                 "in both models. Mixed models: b² v / (b² v + s2_g mean(diag K) + s2_e) from the null-model variance components. "
                 "Multi-locus and iterative models: conditional on the other selected markers. Empty = unavailable. These are not "
                 "claimed to equal the values another program reports.")
    add("Genomic_Lambda", "Median of the chi-square(1) statistics implied by the p-values divided by 0.455. It is reported as a number; "
                          "no p-value was rescaled or corrected by it. For mrMLM / FASTmrMLM the first-step p-values come from a "
                          "random-effect scan and are not on a comparable scale.")
    return rows


# ============================================================================================ workbook
# Layout of the original TerraTools GWAS report: Summary, Threshold_Explanation, one MTA sheet per model, Significant_MTAs,
# QTN_Effects, SNP_QC. Only genome-wide significant markers (QTNs for the multi-locus models) are listed as MTAs.
# Diagnostics, QC_Stages and Notes follow at the end.

_HEAD_FILL = PatternFill("solid", start_color="4A90D9", end_color="4A90D9")
_SIG_FILL = PatternFill("solid", start_color="FFF3CD", end_color="FFF3CD")
_ALT_FILL = PatternFill("solid", start_color="F8F9FA", end_color="F8F9FA")
_WHITE = PatternFill("solid", start_color="FFFFFF", end_color="FFFFFF")
_HEAD_FONT = Font(bold=True, color="FFFFFF", name="Calibri", size=11)
_TITLE_FONT = Font(bold=True, color="2C5F8A", name="Calibri", size=13)
_EDGE = Side(style="thin", color="CCCCCC")
_BORDER = Border(left=_EDGE, right=_EDGE, top=_EDGE, bottom=_EDGE)
_CENTER = Alignment(horizontal="center")


def _clean(v):
    if v is None:
        return None
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, (float, np.floating)):
        v = float(v)
        return v if np.isfinite(v) else None
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    return v if isinstance(v, (int, str)) else str(v)


def _title(ws, text, last_col):
    ws["A1"] = text
    ws["A1"].font = _TITLE_FONT
    ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=last_col)


def _header(ws, row, names, widths):
    for c, (name, w) in enumerate(zip(names, widths), 1):
        cell = ws.cell(row=row, column=c, value=name)
        cell.font, cell.fill, cell.alignment, cell.border = _HEAD_FONT, _HEAD_FILL, _CENTER, _BORDER
        ws.column_dimensions[get_column_letter(c)].width = w


def _put(ws, row, values, fill, sci_col=None):
    for c, v in enumerate(values, 1):
        cell = ws.cell(row=row, column=c, value=_clean(v))
        cell.fill, cell.border, cell.alignment = fill, _BORDER, _CENTER
        if c == sci_col and isinstance(cell.value, float):
            cell.number_format = "0.00E+00"


def _bonf_flag(res, row, alpha):
    """YES / NO from the Bonferroni-adjusted p; multi-locus models have no comparable first-step p-value, so n/a."""
    if res.kind == "multilocus":
        return "n/a"
    pb = row["P_Bonferroni"]
    return "YES" if (pb is not None and pb <= alpha) else "NO"


def _t_stat(effect, se):
    try:
        return round(float(effect) / float(se), 4) if float(se) > 0 else None
    except (TypeError, ValueError):
        return None


def _plain_table(ws, rows, columns, widths=None, sci=()):
    """A table starting at A1 (used for the sheets added after the original ones)."""
    _header(ws, 1, columns, [(widths or {}).get(n, max(12, len(n) + 4)) for n in columns])
    for r, row in enumerate(rows, 2):
        for c, name in enumerate(columns, 1):
            cell = ws.cell(row=r, column=c, value=_clean(row.get(name)))
            cell.border = _BORDER
            if name in sci and isinstance(cell.value, float):
                cell.number_format = "0.00E+00"
    ws.freeze_panes = "A2"
    if rows:
        ws.auto_filter.ref = f"A1:{get_column_letter(len(columns))}{len(rows) + 1}"


def _threshold_rows(A):
    cfg = A.cfg
    rows = [
        ("Bonferroni Threshold", "Formula", f"{cfg.alpha:g} / m  where m = number of SNPs tested"),
        ("Bonferroni Threshold", "Rationale", f"Controls the genome-wide Type I error rate at {100 * cfg.alpha:g}%."),
        ("Bonferroni Threshold", "Assumption", "All SNP tests are independent (conservative when SNPs are in linkage disequilibrium)."),
        ("Benjamini-Hochberg", "Formula", f"Adjusted p <= {cfg.alpha:g} (false discovery rate)."),
        ("Suggestive Threshold", "Formula", "1 / m  where m = number of SNPs tested"),
        ("Suggestive Threshold", "Rationale", "Expects about one false positive per genome scan. Suggestive markers are not genome-wide significant."),
        ("Simulated Threshold", "Formula", "Single-locus models only: the genome-wide minimum p-value of data simulated under the model's own fitted null."),
        ("MTA list", "Rule", f"Only genome-wide significant markers ({cfg.significance_rule}) are listed as MTAs; QTNs for the multi-locus models."),
    ]
    for name, res in A.results.items():
        rows.append((f"{name} Note", "Correction", res.label))
    return rows


def export_excel(A, summary, mta, path):
    cfg = A.cfg
    wb = Workbook()
    n_models = len(A.results)

    # ---- Summary
    ws = wb.active
    ws.title = "Summary"
    ws["A1"] = f"GWAS Complete Analysis Summary - {n_models} model{'s' if n_models != 1 else ''}"
    ws["A1"].font = Font(bold=True, color="2C5F8A", name="Calibri", size=16)
    ws.merge_cells("A1:H1")
    _header(ws, 3, ["Model", "N_SNPs_Tested", "Bonferroni_Threshold", "Suggestive_Threshold", "Sig_SNPs_Bonferroni",
                    "Sig_SNPs_Suggestive", "Genomic_Lambda", "Min_Pvalue"], [24] * 8)
    for r, row in enumerate(summary, 4):
        res = A.results[row["Model"]]
        thr = res.info["sig"]["thresholds"]
        with np.errstate(invalid="ignore"):
            n_sug = int(np.sum(np.asarray(res.p, float) <= thr["suggestive_p"]))
        _put(ws, r, [row["Model"], row["SNPs_Tested"], f"{thr['bonferroni_p']:.2e}", f"{thr['suggestive_p']:.2e}", row["Bonferroni_hits"],
                     n_sug, row["Genomic_Lambda"], row["Min_P"]], _ALT_FILL if r % 2 == 0 else _WHITE)
        ws.cell(row=r, column=9, value=row["Result"])
    ws.cell(row=3, column=9, value="Result").font = Font(bold=True)
    ws.column_dimensions["I"].width = 70

    # ---- Threshold_Explanation
    ws = wb.create_sheet("Threshold_Explanation")
    _title(ws, "Threshold Calculation Explanation", 3)
    _header(ws, 3, ["Threshold", "Aspect", "Details"], [22, 16, 110])
    for r, (t, a, d) in enumerate(_threshold_rows(A), 4):
        for c, v in enumerate((t, a, d), 1):
            cell = ws.cell(row=r, column=c, value=v)
            cell.border = _BORDER
            cell.alignment = Alignment(wrap_text=True, vertical="top")

    by_model = {name: [r for r in mta if r["Model"] == name] for name in A.results}
    n_hits = {name: len(rows) for name, rows in by_model.items()}

    # ---- one MTA sheet per model
    headers = ["Rank", "SNP", "Chr", "Position", "MAF", "Beta", "SE", "T_statistic", "P_value", "neg_log10P", "LOD", "FDR_BH",
               "Bonferroni_sig"]
    widths = [8, 25, 8, 12, 10, 12, 12, 12, 14, 12, 10, 14, 16]
    for name, rows in by_model.items():
        ws = wb.create_sheet(f"{name}_MTA(LOD)"[:31])
        _title(ws, f"GWAS Results - {name} Model - {n_hits[name]} genome-wide significant MTA(s)", len(headers))
        _header(ws, 3, headers, widths)
        if not rows:
            ws.cell(row=4, column=1, value=inf.NO_HITS + ".")
        for rank, r in enumerate(rows, 1):
            _put(ws, rank + 3, [rank, r["SNP"], r["Chr"], r["Position"], r["MAF"], r["Effect"], r["SE"], _t_stat(r["Effect"], r["SE"]),
                                r["P_value"], r["neg_log10P"], r["LOD"], r["P_BH"],
                                _bonf_flag(A.results[name], r, cfg.alpha)], _SIG_FILL, sci_col=9)
        ws.freeze_panes = "A4"

    # ---- Significant_MTAs (all models)
    headers = ["Model", "Rank", "SNP", "Chr", "Position", "MAF", "Beta", "SE", "T_stat", "P_value", "neg_log10P", "LOD", "FDR_BH_p",
               "Bonferroni_sig", "PVE_%"]
    ws = wb.create_sheet("Significant_MTAs")
    _title(ws, "Significant MTAs - All Models (genome-wide significant markers only)", len(headers))
    _header(ws, 3, headers, [16] * len(headers))
    row_i = 4
    for name, rows in by_model.items():
        if not rows:
            ws.cell(row=row_i, column=1, value=f"{name}: no MTAs passed threshold")
            row_i += 1
            continue
        for rank, r in enumerate(rows, 1):
            _put(ws, row_i, [name, rank, r["SNP"], r["Chr"], r["Position"], r["MAF"], r["Effect"], r["SE"], _t_stat(r["Effect"], r["SE"]),
                             r["P_value"], r["neg_log10P"], r["LOD"], r["P_BH"], _bonf_flag(A.results[name], r, cfg.alpha),
                             r["PVE_%"]], _SIG_FILL, sci_col=10)
            row_i += 1
    ws.freeze_panes = "A4"

    # ---- QTN_Effects
    headers = ["Model", "Rank", "SNP", "Chr", "Position", "MAF", "QTN_Effect_Beta", "SE", "95%_CI_Lower", "95%_CI_Upper", "P_value", "LOD",
               "PVE_%"]
    ws = wb.create_sheet("QTN_Effects")
    _title(ws, "Quantitative Trait Nucleotides (QTNs) - genome-wide significant markers per model", len(headers))
    _header(ws, 3, headers, [17] * len(headers))
    row_i = 4
    for name, rows in by_model.items():
        hits = rows
        if not hits:
            ws.cell(row=row_i, column=1, value=f"{name}: no QTNs passed threshold")
            row_i += 1
            continue
        for rank, r in enumerate(hits, 1):
            b, se = r["Effect"], r["SE"]
            lo = None if b is None or se is None else round(b - 1.96 * se, 6)
            hi = None if b is None or se is None else round(b + 1.96 * se, 6)
            _put(ws, row_i, [name, rank, r["SNP"], r["Chr"], r["Position"], r["MAF"], b, se, lo, hi, r["P_value"], r["LOD"], r["PVE_%"]],
                 _ALT_FILL if row_i % 2 == 0 else _WHITE, sci_col=11)
            row_i += 1
    ws.freeze_panes = "A4"

    # ---- SNP_QC (the markers that were analysed)
    ws = wb.create_sheet("SNP_QC")
    _title(ws, "SNP Quality Control Statistics (markers analysed after quality control)", 6)
    _header(ws, 3, ["SNP", "Chr", "Position", "MAF", "Call_Rate", "Category"], [20] * 6)
    maf, cr = np.asarray(A.maf, float), np.asarray(A.qc.call_rate, float)
    for i in range(min(len(A.snp), EXCEL_MAX_ROWS - 3)):
        _put(ws, i + 4, [str(A.snp[i]), int(A.chrom[i]), int(A.pos[i]), round(float(maf[i]), 4) if np.isfinite(maf[i]) else None,
                         round(float(cr[i]), 4) if np.isfinite(cr[i]) else None,
                         "Rare (<0.05)" if np.isfinite(maf[i]) and maf[i] < 0.05 else "Common"], _ALT_FILL if i % 2 == 0 else _WHITE)
    ws.freeze_panes = "A4"

    # ---- added after the original sheets
    _plain_table(wb.create_sheet("Diagnostics"), diagnostics_rows(A), DIAG_COLUMNS, {"Convergence": 34, "P_histogram_20bins": 60})
    _plain_table(wb.create_sheet("QC_Stages"), A.qc.ledger, QC_COLUMNS, {"stage": 24, "rule": 62})
    ws = wb.create_sheet("Notes")
    _header(ws, 1, ["Item", "Details"], [30, 110])
    for r, (k, v) in enumerate(info_rows(A) + [("Analysis notes", n) for n in A.notes], 2):
        ws.cell(row=r, column=1, value=_clean(k)).alignment = Alignment(vertical="top")
        ws.cell(row=r, column=2, value=_clean(v)).alignment = Alignment(wrap_text=True, vertical="top")

    wb.save(path)
    load_workbook(path).close()          # raises if the file just written is unreadable
    return path
