"""
Evidence scoring.

Each of the five evidence categories gets a score E between 0 and 1, or NONE when the evidence cannot be assessed
(no data for that gene). The final score is

    S(g) = [ sum_i( w_i * E(g,i) ) / sum_i( w_i ) ] * 100        over the categories that are available for g

so a missing category is left out of the denominator and is never counted as negative evidence. A category that was
assessed and found unrelated to the trait scores 0 (that is information). The score ranks genes by how much
database evidence relates them to the trait; it is not a probability that a gene is causal.

No evidence item is counted twice: GO terms are merged by GO id, pathways by Reactome id, papers by PMID, orthologs by
gene id, and one item is scored in one category only (trait-ontology annotations in "association", GO text in
"function", Reactome pathways in "pathway", ortholog texts in "ortholog", papers in "literature").

"""

from __future__ import annotations

import math

from tools.gene_prioritization import config
from tools.gene_prioritization.backend import traits

EXPERIMENTAL = {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HDA", "HMP", "HGI", "HEP"}
CURATED = {"TAS", "IC", "NAS", "ISS", "ISO", "ISA", "ISM", "IGC", "IBA", "IBD", "IKR", "IRD", "RCA"}
KIND_STRENGTH = {"direct": 0.6, "related": 0.4}


def _code_weight(code):
    code = (code or "").upper()
    if code in EXPERIMENTAL:
        return 1.0
    if code in CURATED:
        return 0.8
    return 0.6                                     # IEA, unknown


def _noisy_or(strengths):
    p = 1.0
    for s in strengths:
        p *= 1.0 - max(0.0, min(1.0, s))
    return 1.0 - p


def _result(score, items, note=""):
    return {"score": None if score is None else round(float(score), 4), "available": score is not None,
            "items": items, "note": note}


def _unavailable(note):
    return _result(None, [], note)


# ======================================================================= association
def association(rec, tt):
    """Known gene-trait association: the gene is annotated with a trait term, or has a phenotype annotation, matching the trait.
    A trait term of a QTL that contains the gene counts less: it is positional."""
    items = []
    for t in rec.get("to", []):
        for term, kind in traits.match(t["name"], tt)[:1]:
            items.append({"text": f"{t['name']} ({t['id']})", "source": "Gramene (Trait Ontology annotation of the gene)",
                          "url": rec.get("gramene_url", ""), "strength": 1.0 if kind == "direct" else 0.6,
                          "kind": "gene-level annotation", "matched": term, "code": t.get("evidence_code", "")})
    for p in rec.get("phenotypes", []):
        for term, kind in traits.match(p["description"], tt)[:1]:
            items.append({"text": p["description"], "source": f"Ensembl Plants phenotype ({p.get('source') or 'n/a'})",
                          "url": rec.get("ensembl_url", ""), "strength": 1.0 if kind == "direct" else 0.6,
                          "kind": "gene-level annotation", "matched": term, "code": ""})
    oz = rec.get("oryzabase_entry")
    if oz:                                           # Oryzabase: curated trait ontology terms and trait classes of the gene
        seen = set()
        for t in oz.get("to", []):
            for term, kind in traits.match(t["name"], tt)[:1]:
                seen.add(t["id"])
                items.append({"text": f"{t['name']} ({t['id']})", "source": "Oryzabase (curated Trait Ontology annotation of the gene)",
                              "url": oz["url"], "strength": 1.0 if kind == "direct" else 0.6, "kind": "gene-level annotation",
                              "matched": term, "code": "curated"})
        for c in oz.get("classes", []):
            for term, kind in traits.match(c, tt)[:1]:
                items.append({"text": f"trait class: {c}", "source": "Oryzabase (trait class of the gene)", "url": oz["url"],
                              "strength": 0.8 if kind == "direct" else 0.5, "kind": "gene-level annotation", "matched": term, "code": "curated"})
    for t in rec.get("qtl_to", []):
        for term, kind in traits.match(t["name"], tt)[:1]:
            items.append({"text": f"inside a QTL for {t['name']} ({t['id']})", "source": "Gramene (QTL trait term)",
                          "url": rec.get("gramene_url", ""), "strength": 0.5 if kind == "direct" else 0.3,
                          "kind": "QTL-based (positional)", "matched": term, "code": ""})
    if not items:
        return _unavailable("No annotation links this gene to the trait in the databases queried.")
    return _result(max(i["strength"] for i in items), items)


# ======================================================================= function
def function(rec, tt):
    """Functional relevance: description, GO terms and protein domains whose text names the trait or its vocabulary."""
    oz = rec.get("oryzabase_entry") or {}
    has_text = bool(rec.get("description") or rec.get("go") or rec.get("domains") or oz.get("name"))
    if not has_text:
        return _unavailable("The gene has no description, GO terms or domains in the databases queried.")
    items = []

    def add(text, where, base, url="", key=None):
        for term, kind in traits.match(text, tt)[:1]:
            items.append({"text": text, "source": where, "url": url, "strength": KIND_STRENGTH[kind] * base,
                          "kind": kind, "matched": term, "key": key or text})

    add(rec.get("description", ""), "Gene description", 1.0, rec.get("gramene_url", ""))
    if oz:
        add(f"{oz.get('symbol', '')} {oz.get('name', '')}".strip(), "Oryzabase gene name", 0.8, oz["url"], key="oryzabase:name")
    seen = set()
    for g in rec.get("go", []):
        if g["id"] in seen:
            continue
        seen.add(g["id"])
        w = _code_weight(g.get("evidence_code"))
        url = f"https://amigo.geneontology.org/amigo/term/{g['id']}"
        label = f"{g['id']} {g.get('name') or ''} [{g.get('evidence_code') or 'IEA'}]"
        add(label, f"GO annotation ({g.get('source', 'Gramene')})", w, url, key=g["id"])
        if g.get("definition"):
            add(f"{g['id']} {g.get('name') or ''}: {g['definition']}", "GO term definition", 0.5 * w, url, key=g["id"])
    for d in rec.get("domains", []):
        add(f"{d['id']} {d.get('name', '')} {d.get('description', '')}", "InterPro domain (via Gramene)", 0.5,
            f"https://www.ebi.ac.uk/interpro/entry/InterPro/{d['id']}/", key=d["id"])
    # one item per GO term / domain: the strongest of its name and its definition
    best = {}
    for it in items:
        if it["key"] not in best or it["strength"] > best[it["key"]]["strength"]:
            best[it["key"]] = it
    items = sorted(best.values(), key=lambda i: -i["strength"])
    return _result(_noisy_or(i["strength"] for i in items), items)


# ======================================================================= pathway
def pathway(rec, tt, trait_pathways):
    """Pathway relevance: the gene's Plant Reactome pathways, judged by their name and by whether Plant Reactome
    itself returns them for the trait words."""
    pws = rec.get("pathways", [])
    if not pws:
        return _unavailable("The gene is not in a Plant Reactome pathway.")
    items = []
    for p in pws:
        hits = traits.match(p["name"], tt)
        in_search = p["id"] in trait_pathways
        if hits:
            term, kind = hits[0]
            s, why = (0.9 if kind == "direct" else 0.6), f"pathway name mentions '{term}'"
        elif in_search:
            s, why = 0.7, f"Plant Reactome returns this pathway for '{trait_pathways[p['id']]['query']}'"
        else:
            s, why = 0.0, "pathway not related to the trait"
        items.append({"text": p["name"], "source": "Plant Reactome", "url": p.get("url", ""), "strength": s, "kind": why,
                      "id": p["id"]})
    rel = [i for i in items if i["strength"] > 0]
    score = max((i["strength"] for i in items), default=0.0)
    if len(rel) > 1:
        score = min(1.0, score + 0.05 * (len(rel) - 1))
    return _result(score, sorted(items, key=lambda i: -i["strength"]))


# ======================================================================= orthologs
TYPE_FACTOR = {"ortholog_one2one": 1.0, "ortholog_one2many": 0.8, "ortholog_many2many": 0.7}


def ortholog(rec, tt):
    """Ortholog evidence: orthologs (Ensembl Plants, Gramene) that are annotated with the trait."""
    ors = [o for o in rec.get("orthologs", []) if o.get("name") or o.get("description") or o.get("phenotypes")]
    if not ors:
        return _unavailable("No annotated orthologs were found.")
    own = (rec.get("description") or "").strip().lower()
    items = []
    for o in ors:
        f = TYPE_FACTOR.get(o.get("type"), 0.7)
        for ph in o.get("phenotypes", []):
            for term, kind in traits.match(ph["description"], tt)[:1]:
                items.append({"text": f"{o['id']} ({o.get('species', '')}): phenotype '{ph['description']}'", "source": "Ensembl Plants phenotype of an ortholog",
                              "url": o.get("url", ""), "strength": (1.0 if kind == "direct" else 0.6) * f, "kind": "ortholog phenotype", "matched": term})
        text = f"{o.get('name', '')} {o.get('description', '')}".strip()
        if text and text.lower() != own:
            for term, kind in traits.match(text, tt)[:1]:
                items.append({"text": f"{o['id']} ({o.get('species', '')}): {text}", "source": o.get("source", "Ensembl Plants / Gramene"),
                              "url": o.get("url", ""), "strength": KIND_STRENGTH[kind] * f, "kind": "ortholog annotation", "matched": term})
    return _result(_noisy_or(i["strength"] for i in items), sorted(items, key=lambda i: -i["strength"]),
                   "" if items else "Orthologs exist but none is annotated with the trait.")


# ======================================================================= literature
def literature(pubs):
    """Literature support from papers whose title/abstract name the gene and the trait. Keyword-only hits are listed, not scored."""
    co = [p for p in pubs if p["relation"] in ("title co-mention", "co-mention")]
    if not co:
        return _unavailable("No paper in the PubMed results names this gene together with the trait.")
    n_title = sum(p["relation"] == "title co-mention" for p in co)
    score = 1.0 - math.exp(-(len(co) + n_title) / 4.0)
    items = [{"text": f"{p['title']} ({p['year']}) PMID {p['pmid']}", "source": "PubMed", "url": p["url"], "strength": 1.0 if p["relation"] == "title co-mention" else 0.6,
              "kind": p["relation"]} for p in co]
    return _result(score, items, "Co-mentions in a title or abstract: they do not by themselves prove a gene-trait relationship.")


# ======================================================================= combine
def combine(cats):
    """cats: {key: result}. -> (final score or None, {key: contribution in score points}, used keys)."""
    used = {k: config.EVIDENCE_WEIGHTS[k][1] for k, r in cats.items() if r["available"]}
    den = sum(used.values())
    if not den:
        return None, {k: 0.0 for k in cats}, []
    num = sum(w * cats[k]["score"] for k, w in used.items())
    contrib = {k: (used[k] * cats[k]["score"] / den * 100.0 if k in used else 0.0) for k in cats}
    return round(num / den * 100.0, 2), {k: round(v, 2) for k, v in contrib.items()}, sorted(used)


def evidence_level(cats):
    assoc = cats["association"]
    if any(i["kind"] == "gene-level annotation" and i["strength"] >= 0.6 for i in assoc["items"]):
        return "Direct association"
    n = sum(1 for r in cats.values() if r["available"] and (r["score"] or 0) > 0)
    return {0: "No trait evidence", 1: "Single line", 2: "Two lines"}.get(n, "Multiple lines")
