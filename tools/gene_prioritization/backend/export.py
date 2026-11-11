"""
Result files: ranked CSV, an Excel workbook (ranked genes + evidence details + publications +
pathways + run information) and a self-contained HTML report.
"""

from __future__ import annotations

import html
from pathlib import Path

import pandas as pd

from tools.gene_prioritization import config

RANK_COLUMNS = ["Rank", "Gene_ID", "Symbol", "Description", "Final_Score", "Evidence_Level", "Categories_Used",
                "E_association", "E_function", "E_pathway", "E_ortholog", "E_literature", "Chromosome", "Start", "End", "Strand"]
DISCLAIMER = ("The score ranks genes by how much evidence in public databases relates them to the trait. It is not a "
              "probability that a gene is causal. A category without data is left out of the score, not counted as negative "
              "evidence. PubMed co-mentions show that a gene and the trait are named together, not that a relationship was demonstrated.")


def _evidence_rows(res):
    rows = []
    for gid, d in res["detail"].items():
        for k, s in d["scores"].items():
            if not s["items"] and s["score"] is None:        # nothing found is not listed; a category that could not be queried is
                rows.append({"Gene_ID": gid, "Symbol": d["symbol"], "Category": s["label"], "Weight": s["weight"],
                             "Category_Score": s["score"], "Evidence": s["note"] or "no evidence", "Evidence_Type": "",
                             "Strength": None, "Source": "", "URL": "", "Retrieved": d["retrieved"]})
            for it in s["items"]:
                rows.append({"Gene_ID": gid, "Symbol": d["symbol"], "Category": s["label"], "Weight": s["weight"],
                             "Category_Score": s["score"], "Evidence": it["text"], "Evidence_Type": it.get("kind", ""),
                             "Strength": it.get("strength"), "Source": it.get("source", ""), "URL": it.get("url", ""),
                             "Retrieved": d["retrieved"]})
    return rows


def _pub_rows(res):
    rows = []
    for gid, d in res["detail"].items():
        for p in d["publications"]:
            rows.append({"Gene_ID": gid, "Symbol": d["symbol"], "PMID": p["pmid"], "Year": p["year"], "Title": p["title"],
                         "Journal": p["journal"], "Authors": p["authors"], "DOI": p["doi"], "Relation": p["relation"], "URL": p["url"]})
    return rows


def _pathway_rows(res):
    return [{"Gene_ID": gid, "Symbol": d["symbol"], "Pathway_ID": p["id"], "Pathway": p["name"], "URL": p["url"]}
            for gid, d in res["detail"].items() for p in d["pathways"]]


def _go_rows(res):
    return [{"Gene_ID": gid, "Symbol": d["symbol"], "GO_ID": g["id"], "Term": g.get("name"), "Aspect": g.get("namespace"),
             "Evidence_Code": g.get("evidence_code"), "Source": g.get("source")} for gid, d in res["detail"].items() for g in d["go"]]


def _run_info(res):
    info = [("Species", res["species"]), ("Trait", res["trait"]), ("Input", "gene list"), ("Run date", res["run_date"]),
            ("Ranked by", res["rank_by"]),
            ("Trait words (direct)", ", ".join(res["trait_terms"]["direct"])), ("Trait words (related)", ", ".join(res["trait_terms"]["related"])),
            ("Score", "S(g) = [sum(w_i x E_gi) / sum(w_i)] x 100 over the evidence categories available for the gene")]
    info += [(f"Weight: {v['label']}", v["weight"]) for v in res["weights"].values()]
    info += [(f"Service: {s['source']}", f"{s['state']} ({s['calls']} requests, {s['cached']} from cache)") for s in res["api_status"]]
    info += [("Warning", w) for w in res["warnings"]]
    info += [("Unresolved input", f"{u['input']}: {u['reason']}") for u in res["unresolved"]]
    info.append(("Note", DISCLAIMER))
    return [{"Item": a, "Value": b} for a, b in info]


def write_csv(res, path):
    df = pd.DataFrame(res["genes"])
    cols = [c for c in RANK_COLUMNS if c in df.columns]
    df[cols].to_csv(path, index=False, encoding="utf-8-sig")


def write_xlsx(res, path):
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        df = pd.DataFrame(res["genes"])
        df[[c for c in RANK_COLUMNS if c in df.columns]].to_excel(xw, sheet_name="Ranked genes", index=False)
        for name, rows in (("Evidence details", _evidence_rows(res)),
                           ("Publications", _pub_rows(res)), ("Pathways", _pathway_rows(res)), ("GO terms", _go_rows(res)),
                           ("Run information", _run_info(res))):
            if rows or name == "Run information":
                pd.DataFrame(rows).to_excel(xw, sheet_name=name, index=False)
        for ws in xw.book.worksheets:
            ws.freeze_panes = "A2"
            for col in ws.columns:
                width = min(60, max(10, max((len(str(c.value)) for c in col[:60] if c.value is not None), default=10) + 2))
                ws.column_dimensions[col[0].column_letter].width = width


# ======================================================================= HTML report
CSS = """
body{font-family:Segoe UI,Arial,sans-serif;color:#1b2b21;max-width:1100px;margin:28px auto;padding:0 18px;line-height:1.45}
h1{color:#1f5d3a;margin-bottom:2px} h2{color:#1f5d3a;border-bottom:2px solid #b7dcc6;padding-bottom:4px;margin-top:30px}
h3{margin:18px 0 4px} .sub{color:#5c6f63;margin-top:0} table{border-collapse:collapse;width:100%;font-size:13px;margin:8px 0}
th{background:#e4f3ea;text-align:left} th,td{border:1px solid #cfe3d6;padding:5px 8px;vertical-align:top}
.note{background:#f4f9f6;border-left:4px solid #52b788;padding:8px 12px;margin:12px 0;font-size:13px}
.bar{background:#e4f3ea;border-radius:3px;height:12px;width:160px;display:inline-block;vertical-align:middle}
.bar i{display:block;background:#2d8a5a;height:12px;border-radius:3px} .muted{color:#6b7c72} a{color:#1f6f9f}
.gene{page-break-inside:avoid;margin-bottom:20px} .tag{display:inline-block;background:#e4f3ea;border-radius:10px;padding:1px 8px;font-size:12px}
"""


def _e(x):
    return html.escape("" if x is None else str(x))


def _fmt(x):
    return "—" if x is None else (f"{x:.1f}" if isinstance(x, float) else str(x))


def _bar(v):
    return f'<span class="bar"><i style="width:{max(0, min(100, v or 0))}%"></i></span> {_fmt(v)}'


def write_html(res, path):
    s = res["summary"]
    out = [f"<!doctype html><html><head><meta charset='utf-8'><title>Gene Prioritization — {_e(res['trait'])}</title><style>{CSS}</style></head><body>",
           f"<h1>Gene Prioritization</h1><p class='sub'>{_e(res['species'])} · trait: <b>{_e(res['trait'])}</b> · {_e(res['run_date'])} · TerraTools</p>",
           f"<div class='note'>{_e(DISCLAIMER)}</div>", "<h2>Summary</h2><table>"]
    for label, key in (("Genes in the input list", "input_genes"), ("Candidate genes", "candidate_genes"),
                       ("Successfully annotated", "annotated"), ("With a known trait association", "with_association"),
                       ("With literature support (co-mention)", "with_literature"),
                       ("Identifiers not matched", "unresolved")):
        out.append(f"<tr><th>{label}</th><td>{s[key]}</td></tr>")
    out.append(f"</table><p>Genes are ranked by the <b>{_e(res['rank_by'].lower())}</b>.</p>")
    out.append("<h2>Ranked genes</h2><table><tr><th>Rank</th><th>Gene</th><th>Symbol</th><th>Description</th><th>Final score</th><th>Evidence level</th></tr>")
    for r in res["genes"][:200]:
        out.append(f"<tr><td>{r['Rank']}</td><td>{_e(r['Gene_ID'])}</td><td>{_e(r['Symbol'])}</td><td>{_e(r['Description'])}</td><td>{_bar(r['Final_Score'])}</td>"
                   f"<td>{_e(r['Evidence_Level'])}</td></tr>")
    out.append("</table>")
    if len(res["genes"]) > 200:
        out.append(f"<p class='muted'>The first 200 of {len(res['genes'])} genes are shown; the CSV and Excel files have all of them.</p>")
    out.append("<h2>Evidence for the top genes</h2>")
    for r in res["genes"][:25]:
        d = res["detail"][r["Gene_ID"]]
        out.append(f"<div class='gene'><h3>{r['Rank']}. {_e(d['symbol'] or d['id'])} <span class='tag'>{_e(d['evidence_level'])}</span></h3>"
                   f"<div class='muted'>{_e(d['id'])} · {_e(d['location'])} · <a href='{_e(d['links']['Ensembl Plants'])}'>Ensembl Plants</a> · "
                   f"<a href='{_e(d['links']['Gramene'])}'>Gramene</a></div><p>{_e(d['description']) or '<i>no description</i>'}</p>"
                   "<table><tr><th>Evidence category</th><th>Weight</th><th>Score (0–1)</th><th>Evidence</th></tr>")
        for k, sc in d["scores"].items():
            if not sc["items"] and sc["score"] is not None:      # nothing found: not shown
                continue
            ev = "<br>".join(f"{_e(it['text'])} <span class='muted'>— {_e(it['source'])}</span>"
                             + (f" <a href='{_e(it['url'])}'>link</a>" if it.get("url") else "") for it in sc["items"][:5])
            shown = "unavailable" if sc["score"] is None else "%.2f" % sc["score"]
            ev = ev or "<span class='muted'>" + _e(sc["note"]) + "</span>"
            out.append(f"<tr><td>{_e(sc['label'])}</td><td>{int(sc['weight'] * 100)}%</td><td>{shown}</td><td>{ev}</td></tr>")
        out.append("</table></div>")
    out.append("<h2>Methods</h2><ul>"
               "<li><b>Score.</b> S(g) = [Σ w<sub>i</sub> · E(g,i) / Σ w<sub>i</sub>] × 100 over the evidence categories available for the gene. Weights: "
               + ", ".join(f"{_e(v['label'])} {int(v['weight'] * 100)}%" for v in res["weights"].values()) + ".</li>"
               f"<li><b>Trait words.</b> direct: {_e(', '.join(res['trait_terms']['direct']))}; related: {_e(', '.join(res['trait_terms']['related']))}.</li>"
               "<li><b>Sources.</b> Ensembl Plants REST, Gramene, Plant Reactome, NCBI E-utilities (PubMed, Gene) and, for rice, Oryzabase (curated genes keyed by RAP-DB ID). Each evidence item lists its source and link; "
               "retrieval dates are in the Excel workbook.</li></ul>")
    if res["warnings"]:
        out.append("<h2>Notes</h2><ul>" + "".join(f"<li>{_e(w)}</li>" for w in res["warnings"]) + "</ul>")
    if res["unresolved"]:
        out.append("<h2>Identifiers not matched</h2><ul>" + "".join(f"<li>{_e(u['input'])}: {_e(u['reason'])}</li>" for u in res["unresolved"]) + "</ul>")
    out.append("<h2>Services</h2><table><tr><th>Service</th><th>State</th><th>Requests</th><th>From cache</th></tr>"
               + "".join(f"<tr><td>{_e(a['source'])}</td><td>{_e(a['state'])}</td><td>{a['calls']}</td><td>{a['cached']}</td></tr>" for a in res["api_status"])
               + "</table></body></html>")
    Path(path).write_text("".join(out), encoding="utf-8")


def write_all(res, directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    files = {"csv": directory / "gene_prioritization_ranked_genes.csv",
             "xlsx": directory / "gene_prioritization_results.xlsx",
             "html": directory / "gene_prioritization_report.html"}
    write_csv(res, files["csv"])
    write_xlsx(res, files["xlsx"])
    write_html(res, files["html"])
    return {k: str(v) for k, v in files.items()}
