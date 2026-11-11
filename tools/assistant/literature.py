"""Live literature search for the assistant.

Queries Europe PMC (covers PubMed/MEDLINE, PMC, preprints and agricola-type
records) and falls back to Crossref. Every record returned carries the
identifiers exactly as the service reported them (DOI / PMID), so the model
can only cite papers that really exist. No API keys are needed.
"""

import html
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CROSSREF = "https://api.crossref.org/works"
UA = "TerraTools-Assistant/1.0 (mailto:terratools@irri.org)"
TIMEOUT = 20


def _get_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:      # these public APIs return a transient 502/503 now and then
            if e.code < 500 or attempt == 2:
                raise
            time.sleep(1.5 * (attempt + 1))


def _clean(text, limit=None):
    text = html.unescape(text or "")            # Europe PMC sends titles like "&lt;i&gt;Oryza&lt;/i&gt;"
    text = re.sub(r"<[^>]+>", "", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text[:limit] + "…" if limit and len(text) > limit else text


def _epmc(query, n, year_from, open_access_only, sort="relevance"):
    q = f"({query})"
    if year_from:
        q += f" AND (PUB_YEAR:[{int(year_from)} TO 3000])"
    if open_access_only:
        q += " AND OPEN_ACCESS:y"
    params = {
        "query": q, "format": "json", "resultType": "core", "pageSize": max(1, min(int(n), 10)),
    }                                           # (relevance is the default order; "sort=RELEVANCE" gives a 503)
    if sort == "cited":                         # foundational method papers are the most cited ones
        params["sort"] = "CITED desc"
    data = _get_json(EPMC + "?" + urllib.parse.urlencode(params))
    out = []
    for r in data.get("resultList", {}).get("result", []):
        pmid, doi = r.get("pmid"), r.get("doi")
        out.append({
            "title": _clean(r.get("title")),
            "authors": _clean(r.get("authorString")),
            "journal": _clean(r.get("journalTitle") or (r.get("journalInfo") or {}).get("journal", {}).get("title")),
            "year": r.get("pubYear"),
            "doi": doi,
            "pmid": pmid,
            "url": (f"https://doi.org/{doi}" if doi else
                    f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else
                    f"https://europepmc.org/article/{r.get('source')}/{r.get('id')}"),
            "cited_by": r.get("citedByCount"),
            "abstract": _clean(r.get("abstractText"), 900),
            "source": "Europe PMC",
        })
    return out


def _crossref(query, n, year_from):
    params = {"query.bibliographic": query, "rows": max(1, min(int(n), 10)),
              "select": "DOI,title,author,container-title,issued,abstract,is-referenced-by-count"}
    if year_from:
        params["filter"] = f"from-pub-date:{int(year_from)}"
    data = _get_json(CROSSREF + "?" + urllib.parse.urlencode(params))
    out = []
    for r in data.get("message", {}).get("items", []):
        authors = ", ".join(
            f"{a.get('family', '')} {''.join(w[0] for w in a.get('given', '').split())}".strip()
            for a in (r.get("author") or [])[:8])
        year = ((r.get("issued") or {}).get("date-parts") or [[None]])[0][0]
        out.append({
            "title": _clean((r.get("title") or [""])[0]),
            "authors": authors,
            "journal": _clean((r.get("container-title") or [""])[0]),
            "year": str(year) if year else None,
            "doi": r.get("DOI"),
            "pmid": None,
            "url": f"https://doi.org/{r.get('DOI')}",
            "cited_by": r.get("is-referenced-by-count"),
            "abstract": _clean(r.get("abstract"), 900),
            "source": "Crossref",
        })
    return out


def search(query, max_results=6, year_from=None, open_access_only=False, sort="relevance"):
    """Return {'results': [...], 'note': str}. Never raises."""
    query = (query or "").strip()
    if not query:
        return {"results": [], "note": "Empty query."}
    notes = []
    results = []
    try:
        results = _epmc(query, max_results, year_from, open_access_only, sort)
    except Exception as e:                       # offline, rate-limited, bad query syntax
        notes.append(f"Europe PMC unavailable ({type(e).__name__}).")
    if len(results) < 3 and not open_access_only:
        try:
            seen = {r["doi"].lower() for r in results if r.get("doi")}
            for r in _crossref(query, max_results, year_from):
                if r["doi"] and r["doi"].lower() not in seen:
                    results.append(r)
        except Exception as e:
            notes.append(f"Crossref unavailable ({type(e).__name__}).")
    if not results:
        notes.append("No papers found. Try fewer or different keywords; do not cite from memory.")
    return {"results": results[: max(1, int(max_results))], "note": " ".join(notes)}


def lookup_doi(doi):
    """Verify a DOI a user (or the model) mentions and return its bibliographic record."""
    doi = re.sub(r"^(https?://(dx\.)?doi\.org/|doi:)", "", (doi or "").strip(), flags=re.I)
    if not doi:
        return {"found": False, "note": "No DOI given."}
    try:
        r = _get_json(f"{CROSSREF}/{urllib.parse.quote(doi)}")["message"]
    except Exception as e:
        return {"found": False, "note": f"Could not verify {doi} ({type(e).__name__}); treat it as unverified."}
    year = ((r.get("issued") or {}).get("date-parts") or [[None]])[0][0]
    return {
        "found": True, "doi": r.get("DOI"),
        "title": _clean((r.get("title") or [""])[0]),
        "journal": _clean((r.get("container-title") or [""])[0]),
        "year": year, "url": f"https://doi.org/{r.get('DOI')}",
        "authors": ", ".join(f"{a.get('family', '')} {a.get('given', '')[:1]}" for a in (r.get("author") or [])[:8]),
    }


if __name__ == "__main__":
    import sys
    print(json.dumps(search(" ".join(sys.argv[1:]) or "QTL-seq SNP-index", 3), indent=1)[:3000])
