"""
Builds the TerraTools user manual:  docs/manual/TerraTools_User_Manual.html -> .pdf

    python docs/manual/build_manual.py
Needs the screenshots (python docs/make_screenshots.py) and Microsoft Edge + Node for the PDF.

    python docs/manual/build_manual.py --help-only
only refreshes the in-app Help page (frontend/public/help + frontend/src/lib/manualChapters.json), no PDF.
"""
import json
import shutil
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import content_a as A   # noqa: E402
import content_b as B   # noqa: E402
import content_c as C   # noqa: E402
import mdkit            # noqa: E402

HTML = HERE / "TerraTools_User_Manual.html"
PDF = HERE / "TerraTools_User_Manual.pdf"

# the in-app Help page shows this copy of the manual (frontend/public is copied into frontend/dist by the build)
FRONTEND = HERE.parent.parent / "frontend"
HELP_DIR = FRONTEND / "public" / "help"
HELP_CHAPTERS = FRONTEND / "src" / "lib" / "manualChapters.json"

# the manual is laid out for A4 paper; on screen it becomes one continuous page
SCREEN_CSS = """
@media screen {
  html { font-size: 15px; }
  body { max-width: 900px; margin: 0 auto; padding: 22px 28px 90px; background: #fff; }
  .cover { height: auto; margin: 0 0 18px; padding: 34px 30px; border-radius: 14px; break-after: auto; }
  .cover .leaf { font-size: 46px; margin-bottom: 10px; }
  .cover h1 { font-size: 38px; }
  .toc { display: none; }
  .part { height: auto; padding: 46px 0 6px; break-before: auto; }
  .part h1 { font-size: 30px; }
  .chapter { break-before: auto; border-top: 1px solid #d8eadf; margin-top: 34px; padding-top: 30px; }
  figure.shot img { max-height: none; }
  h1 { font-size: 28px; }
  h2 { font-size: 20px; }
  h3 { font-size: 16px; }
  code { font-size: 13px; }
  pre.code { font-size: 12.5px; }
  table { font-size: 13.5px; }
  table.ex { font-size: 12px; }
  .excap, .exnote, .samples, ol.legend li, figure.shot figcaption { font-size: 13px; }
  .box { font-size: 14px; }
}
"""

CSS = """
@page { size: A4; margin: 17mm 16mm 20mm 16mm; }
* { box-sizing: border-box; }
html { font-family: "Segoe UI", "Helvetica Neue", Arial, sans-serif; color: #1c2b24; font-size: 10pt; line-height: 1.5; }
body { margin: 0; }
h1, h2, h3, .cover h1 { font-family: Georgia, "Times New Roman", serif; color: #1b4332; line-height: 1.2; }
h1 { font-size: 24pt; margin: 0 0 4mm; padding-bottom: 3mm; border-bottom: 3px solid #52b788; }
h2 { font-size: 15pt; margin: 8mm 0 3mm; }
h3 { font-size: 12pt; margin: 6mm 0 2mm; color: #2d6a4f; }
h1, h2, h3 { break-after: avoid; }
p { margin: 0 0 3mm; }
.chapter { break-before: page; }
.part { break-before: page; height: 230mm; display: flex; flex-direction: column; justify-content: center; }
.part h1 { border: 0; font-size: 34pt; }
.part p { font-size: 13pt; color: #40594c; }
.lead { font-size: 11pt; color: #29443a; border-left: 4px solid #52b788; padding: 1mm 0 1mm 4mm; background: #f2f9f5; }
.chnum { font: 600 9pt "Segoe UI", sans-serif; letter-spacing: .12em; text-transform: uppercase; color: #52b788; margin: 0 0 1mm; }
code { font-family: Consolas, "Courier New", monospace; font-size: 8.8pt; background: #eef5f1; padding: 0 1.2mm; border-radius: 1mm; color: #1b4332; }
pre.code { font-family: Consolas, monospace; font-size: 8.4pt; background: #14352a; color: #d8f3dc; padding: 3mm 4mm; border-radius: 2mm; white-space: pre-wrap; break-inside: avoid; }
ul, ol { margin: 0 0 3mm 5mm; padding-left: 4mm; }
li { margin-bottom: 1.2mm; }
ol.steps { counter-reset: s; list-style: none; padding-left: 0; margin-left: 0; }
ol.steps > li { counter-increment: s; position: relative; padding-left: 9mm; margin-bottom: 2mm; break-inside: avoid; }
ol.steps > li::before { content: counter(s); position: absolute; left: 0; top: 0.2mm; width: 6mm; height: 6mm; border-radius: 50%; background: #2d6a4f; color: #fff; font: 700 8.5pt "Segoe UI"; text-align: center; line-height: 6mm; }
figure.shot { margin: 4mm 0 5mm; }
figure.shot img { display: block; width: auto; max-width: 100%; max-height: 205mm; margin: 0 auto; border: 1px solid #cfe3d8; border-radius: 2mm; box-shadow: 0 1px 4px rgba(0,0,0,.12); break-inside: avoid; break-after: avoid; }
figure.shot figcaption { font-size: 8.8pt; color: #3d5a4d; margin: 1.5mm 0 1mm; break-before: avoid; }
ol.legend { counter-reset: l; list-style: none; margin: 1mm 0 0; padding: 0; columns: 1; }
ol.legend li { counter-increment: l; position: relative; padding-left: 7mm; font-size: 8.8pt; margin-bottom: 1mm; break-inside: avoid; }
ol.legend li::before { content: counter(l); position: absolute; left: 0; top: .2mm; width: 5mm; height: 5mm; line-height: 5mm; text-align: center; border-radius: 50%; background: #e65a14; color: #fff; font: 700 8pt "Segoe UI"; }
table { border-collapse: collapse; width: 100%; margin: 2mm 0 4mm; font-size: 9pt; }
table.tbl th { background: #2d6a4f; color: #fff; text-align: left; padding: 1.6mm 2.4mm; font-weight: 600; }
table.tbl td { padding: 1.6mm 2.4mm; vertical-align: top; border-bottom: 1px solid #dbe9e1; }
table.tbl tbody tr:nth-child(even) td { background: #f6faf7; }
table.tbl tr { break-inside: avoid; }
.example { margin: 3mm 0 4mm; break-inside: avoid; }
.excap { font-size: 8.6pt; color: #3d5a4d; margin-bottom: 1mm; }
table.ex { font-family: Consolas, monospace; font-size: 8pt; width: auto; max-width: 100%; margin: 0; border: 1px solid #b9d3c4; }
table.ex th { background: #d8f3dc; color: #1b4332; padding: 1mm 2.2mm; text-align: left; border-right: 1px solid #b9d3c4; white-space: nowrap; }
table.ex td { padding: .8mm 2.2mm; border-top: 1px solid #e3efe8; border-right: 1px solid #e3efe8; white-space: nowrap; background: #fff; }
.exnote { font-size: 8pt; color: #5a7566; margin-top: 1mm; }
.box { border-radius: 2mm; padding: 2.4mm 3.5mm; margin: 3mm 0 4mm; font-size: 9.4pt; break-inside: avoid; border-left: 4px solid; }
.box.tip { background: #eaf7ef; border-color: #40a56d; }
.box.warn { background: #fff3e6; border-color: #e8832a; }
.box.note { background: #eef3fb; border-color: #5a86c8; }
.box.key { background: #fdf6d8; border-color: #d4a400; }
.samples { font-size: 8.6pt; color: #4a6558; }
.missing { background: #fde8e8; color: #9b1c1c; padding: 3mm; border-radius: 2mm; }
/* cover */
.cover { height: 297mm; display: flex; flex-direction: column; justify-content: space-between; color: #fff; padding: 26mm 20mm 22mm; margin: -17mm -16mm -37mm; break-after: page; background: linear-gradient(160deg,#0b2b1f 0%,#1b4332 55%,#2d6a4f 100%); }
.cover h1 { color: #fff; font-size: 40pt; border: 0; margin: 0; }
.cover .sub { font-size: 15pt; color: #b7e4c7; margin-top: 4mm; max-width: 130mm; line-height: 1.35; }
.cover .brand { font: 600 11pt "Segoe UI"; letter-spacing: .2em; text-transform: uppercase; color: #95d5b2; }
.cover .meta { font-size: 10.5pt; color: #d8f3dc; }
.cover .leaf { font-size: 70pt; line-height: 1; margin-bottom: 8mm; }
/* contents */
.toc { break-before: page; }
.toc h1 { border: 0; }
.toc ol { list-style: none; margin: 0; padding: 0; }
.toc li { margin: 0; }
.toc a { color: inherit; text-decoration: none; display: flex; align-items: baseline; gap: 2mm; padding: 1mm 0; border-bottom: 1px dotted #b9d3c4; }
.toc a .t { flex: 1; }
.toc a .n { width: 9mm; color: #52b788; font-weight: 700; }
.toc a .p { color: #2d6a4f; font-variant-numeric: tabular-nums; }
.toc .grp { font: 700 9pt "Segoe UI"; letter-spacing: .1em; text-transform: uppercase; color: #52b788; margin: 5mm 0 1mm; }
"""

CHAPTERS = []   # (id, number, title) for the contents page


def chapter(cid, num, title, body, grp=None):
    CHAPTERS.append((cid, num, title, grp))
    return (f'<section class="chapter" id="{cid}"><div class="chnum">{grp + " · " if grp else ""}{("Chapter " + str(num)) if num else ""}</div>'
            f'<h1>{title}</h1>{body}</section>')


def part(pid, title, text):
    return f'<section class="part" id="{pid}"><div class="chnum">TerraTools User Manual</div><h1>{title}</h1><p>{text}</p></section>'


def build(pages=None, screen=False):
    mdkit._fig["n"] = 0
    CHAPTERS.clear()
    out = []
    out.append(part("p1", "Part 1 — Getting started", "Installing TerraTools, finding your way around, and preparing your data."))
    out.append(chapter("c-welcome", 1, "Welcome to TerraTools", A.front_welcome(), "Getting started"))
    out.append(chapter("c-install", 2, "Installing and starting", A.front_install(), "Getting started"))
    out.append(chapter("c-tour", 3, "Finding your way around", A.front_tour(), "Getting started"))
    out.append(chapter("c-data", 4, "Preparing your data", A.front_data(), "Getting started"))
    out.append(part("p2", "Part 2 — The tools", "One chapter per tool: the files it needs, what to click, what the settings mean and how to read the results."))
    tools = [("genofile", C.ch_genofile), ("phenotypic", B.ch_phenotypic), ("mendelian", A.ch_mendelian), ("linkage", A.ch_linkage), ("genotype", C.ch_genotype), ("qtl", A.ch_qtl),
             ("qtlseq", A.ch_qtlseq), ("diversity", A.ch_diversity), ("gwas", B.ch_gwas), ("pybsa", A.ch_pybsa), ("gene", A.ch_gene), ("prioritize", A.ch_prioritize),
             ("clustering", B.ch_clustering), ("additive", B.ch_additive), ("breedsight", B.ch_breedsight), ("cross", B.ch_cross),
             ("heterosis", B.ch_heterosis), ("favallele", C.ch_favallele)]
    names = {k: (n, g) for k, n, g, _d in A.TOOLS}
    for i, (k, fn) in enumerate(tools):
        n, g = names[k]
        out.append(chapter(f"c-{k}", 5 + i, n, fn(), g))
    out.append(part("p3", "Part 3 — Reference", "Typical workflows, troubleshooting, a glossary and a one-page summary of every input format."))
    k = 5 + len(tools)
    out.append(chapter("c-work", k, "Typical workflows", B.appendix_workflows(), "Reference"))
    out.append(chapter("c-trouble", k + 1, "Troubleshooting", B.appendix_troubleshooting(), "Reference"))
    out.append(chapter("c-glossary", k + 2, "Glossary", B.appendix_glossary(), "Reference"))
    out.append(chapter("c-formats", k + 3, "Input formats at a glance", B.appendix_formats(), "Reference"))

    toc = ['<section class="toc"><h1>Contents</h1><ol>']
    last = None
    for cid, num, title, grp in CHAPTERS:
        if grp != last:
            toc.append(f'<li class="grp">{grp}</li>')
            last = grp
        pg = (pages or {}).get(cid, "")
        toc.append(f'<li><a href="#{cid}"><span class="n">{num}</span><span class="t">{title}</span><span class="p">{pg}</span></a></li>')
    toc.append("</ol></section>")

    cover = """<section class="cover"><div><div class="brand">IRRI · SAH</div></div>
<div><div class="leaf">🌱</div><h1>TerraTools</h1><div class="sub">User manual — plant genetics and breeding analyses in one desktop application</div></div>
<div class="meta">Nineteen tools · input formats · step-by-step guides<br>Version for TerraTools 1.0 · October 2026</div></section>"""
    return (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><title>TerraTools User Manual</title><style>{CSS}{SCREEN_CSS if screen else ""}</style></head><body>'
            + cover + "".join(toc) + "".join(out) + "</body></html>")


def render():
    r = subprocess.run(["node", str(HERE / "render_pdf.mjs"), str(HTML), str(PDF)], capture_output=True, text=True)
    print(r.stdout.strip(), r.stderr.strip()[:400])
    return r.returncode == 0


def page_numbers():
    """Find the PDF page on which each chapter starts (by its heading text)."""
    from pypdf import PdfReader
    rd = PdfReader(str(PDF))
    texts = [(p.extract_text() or "") for p in rd.pages]
    pages = {}
    for cid, num, title, grp in CHAPTERS:
        marker = re.sub(r"\s+", " ", title).lower()
        for i, t in enumerate(texts):
            if i < 2:
                continue            # cover and contents
            flat = re.sub(r"\s+", " ", t).lower()
            if ("chapter %s" % num) in flat and marker in flat:
                pages[cid] = i + 1
                break
    return pages, len(texts)


def write_help():
    """The manual as one continuous page for the in-app Help window, with its screenshots and the chapter list."""
    html = build(screen=True)
    HELP_DIR.mkdir(parents=True, exist_ok=True)
    (HELP_DIR / "index.html").write_text(html, encoding="utf-8")
    shots = HELP_DIR / "shots"
    if shots.exists():
        shutil.rmtree(shots)
    shutil.copytree(HERE / "shots", shots, ignore=shutil.ignore_patterns("*.json"))
    chapters = [{"id": cid, "n": num, "title": title, "group": grp} for cid, num, title, grp in CHAPTERS]
    HELP_CHAPTERS.write_text(json.dumps(chapters, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"help page: {len(chapters)} chapters -> {HELP_DIR}")


if __name__ == "__main__":
    if "--help-only" in sys.argv:
        write_help()
        sys.exit(0)
    HTML.write_text(build(), encoding="utf-8")
    if not render():
        sys.exit("PDF rendering failed")
    try:
        pages, total = page_numbers()
        HTML.write_text(build(pages), encoding="utf-8")
        render()
        print(f"pages: {total}; chapters located: {len(pages)}/{len(CHAPTERS)}")
    except ImportError:
        print("pypdf not installed: contents page has no page numbers (pip install pypdf)")
    write_help()
