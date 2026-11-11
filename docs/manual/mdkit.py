"""Small HTML helpers for the manual: figures with numbered legends, example-file previews, boxes."""
from __future__ import annotations

import csv
import gzip
import html
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SAMPLES = ROOT / "docs" / "sample_data"
SHOTS = Path(__file__).resolve().parent / "shots"
MANIFEST = json.loads((SHOTS / "manifest.json").read_text(encoding="utf-8")) if (SHOTS / "manifest.json").exists() else {}

esc = html.escape
_fig = {"n": 0}


def figure(name, caption=None, width=None):
    """Screenshot with the numbered legend taken from the capture manifest."""
    path = SHOTS / f"{name}.jpg"
    if not path.exists():
        return f'<div class="missing">[screenshot "{esc(name)}" was not captured]</div>'
    m = MANIFEST.get(name, {})
    _fig["n"] += 1
    cap = caption or m.get("title", name)
    legend = "".join(f"<li>{esc(t)}</li>" for t in m.get("notes", []) if t)
    style = f' style="max-width:{width}"' if width else ""
    return (f'<figure class="shot"><img src="shots/{name}.jpg"{style} alt="{esc(cap)}">'
            f'<figcaption><b>Figure {_fig["n"]}.</b> {esc(cap)}</figcaption>'
            f'{"<ol class=legend>" + legend + "</ol>" if legend else ""}</figure>')


def _read_rows(path: Path, n_rows: int):
    suf = "".join(path.suffixes).lower()
    rows = []
    if suf.endswith(".xlsx"):
        import openpyxl
        ws = openpyxl.load_workbook(path, read_only=True).worksheets[0]
        for r in ws.iter_rows(values_only=True):
            rows.append(["" if v is None else (f"{v:g}" if isinstance(v, float) else str(v)) for v in r])
            if len(rows) >= n_rows:
                break
        return rows
    opener = gzip.open if suf.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig", errors="replace", newline="") as f:
        if ".vcf" in suf:
            for line in f:
                if line.startswith("##"):
                    continue
                rows.append(line.rstrip("\n").split("\t"))
                if len(rows) >= n_rows:
                    break
            return rows
        dialect = "excel-tab" if suf.endswith((".tsv", ".txt")) else "excel"
        for r in csv.reader(f, dialect=dialect):
            rows.append(r)
            if len(rows) >= n_rows:
                break
    return rows


def example(rel, rows=5, cols=7, caption=None, skip_cols=None, header_rows=1, tail_cols=0, note=None):
    """Preview of an example file (first rows / columns) as a table; the file lives in docs/sample_data."""
    path = SAMPLES / rel
    data = _read_rows(path, rows + header_rows)
    if skip_cols:
        data = [[c for i, c in enumerate(r) if i not in skip_cols] for r in data]
    width = max(len(r) for r in data)
    show = cols
    out = []
    for r in data:
        r = r + [""] * (width - len(r))
        cut = r[:show]
        if width > show:
            cut = cut + ["…"] + (r[-tail_cols:] if tail_cols else [])
        out.append(cut)
    head = "".join(f"<th>{esc(c)}</th>" for c in out[0])
    body = "".join("<tr>" + "".join(f"<td>{esc(c)}</td>" for c in r) + "</tr>" for r in out[1:])
    cap = caption or f"Example: <code>{esc(path.name)}</code>"
    nt = f'<div class="exnote">{esc(note)}</div>' if note else ""
    return f'<div class="example"><div class="excap">{cap}</div><table class="ex"><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>{nt}</div>'


def table(header, rows, cls="", widths=None):
    th = "".join(f"<th>{h}</th>" for h in header)
    tb = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f'<table class="tbl {cls}"><thead><tr>{th}</tr></thead><tbody>{tb}</tbody></table>'


def ul(items):
    return "<ul>" + "".join(f"<li>{i}</li>" for i in items) + "</ul>"


def ol(items):
    return "<ol class=steps>" + "".join(f"<li>{i}</li>" for i in items) + "</ol>"


def box(kind, text, title=None):
    t = {"tip": "Tip", "warn": "Watch out", "note": "Note", "key": "Key point"}[kind]
    return f'<div class="box {kind}"><b>{title or t}.</b> {text}</div>'


def sample_files(rel_dir):
    files = sorted((SAMPLES / rel_dir).iterdir())
    return "<p class=samples>Example files: " + ", ".join(f"<code>sample_data/{rel_dir}/{f.name}</code>" for f in files) + "</p>"
