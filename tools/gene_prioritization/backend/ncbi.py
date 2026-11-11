"""
NCBI E-utilities (https://www.ncbi.nlm.nih.gov/books/NBK25501/): ESearch, ESummary and EFetch on PubMed and Gene.

A key in the environment variable NCBI_API_KEY raises the allowed rate from 3 to 10 requests per second.

A PubMed hit is NOT a verified gene-trait relationship. Each paper is classified from its own title and abstract:
  "co-mention"       the gene (symbol or ID) and a trait word both appear in the title or abstract
  "title co-mention" ... and both appear in the title
  "keyword only"     PubMed matched the query but the text does not name the gene
Only co-mentions count towards the literature score, and they are reported as co-mentions.
"""

from __future__ import annotations

import os
import re
import xml.etree.ElementTree as ET

from tools.gene_prioritization import config
from tools.gene_prioritization.backend.http import ApiClient, today

SOURCE = "NCBI (PubMed / Gene)"


class Ncbi:
    def __init__(self, ctx):
        key = os.environ.get(config.NCBI_API_KEY_ENV, "").strip()
        extra = {"tool": "TerraTools", "retmode": "json"}
        if key:
            extra["api_key"] = key
        if os.environ.get(config.NCBI_EMAIL_ENV):
            extra["email"] = os.environ[config.NCBI_EMAIL_ENV]
        self.has_key = bool(key)
        self.api = ApiClient(SOURCE, config.NCBI_URL, ctx, min_interval=0.11 if key else 0.36, extra_params=extra)

    @staticmethod
    def pubmed_url(pmid):
        return f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

    @staticmethod
    def gene_url(gene_id):
        return f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}"

    def ping(self):
        r = self.api.request("/einfo.fcgi", ttl=0, retries=2, timeout=12)
        return bool(r)

    # ------------------------------------------------------------ ESearch / ESummary
    def search(self, db, term, retmax=20, sort=None):
        params = {"db": db, "term": term, "retmax": retmax}
        if sort:
            params["sort"] = sort
        r = self.api.request("/esearch.fcgi", params=params)
        res = (r or {}).get("esearchresult", {})
        return [str(i) for i in res.get("idlist", [])], int(res.get("count", 0) or 0)

    def summary(self, db, ids):
        out = {}
        ids = [str(i) for i in dict.fromkeys(ids)]
        for i in range(0, len(ids), 150):
            chunk = ids[i:i + 150]
            r = self.api.request("/esummary.fcgi", params={"db": db, "id": ",".join(chunk)}) or {}
            res = r.get("result", {})
            for uid in res.get("uids", []):
                out[str(uid)] = res.get(str(uid), {})
        return out

    # ------------------------------------------------------------ Gene database (identifier mapping)
    def gene_ids_for(self, term, taxon):
        ids, _ = self.search("gene", f"{term}[All Fields] AND txid{taxon}[Organism:exp]", retmax=5)
        return ids

    # ------------------------------------------------------------ PubMed records
    def articles(self, pmids):
        """{pmid: {pmid, title, journal, year, authors, doi, abstract, url}} through ESummary and EFetch."""
        pmids = [str(p) for p in dict.fromkeys(pmids)]
        summ = self.summary("pubmed", pmids)
        abstracts = {}
        for i in range(0, len(pmids), 80):
            chunk = pmids[i:i + 80]
            xml = self.api.request("/efetch.fcgi", params={"db": "pubmed", "id": ",".join(chunk), "retmode": "xml",
                                                          "rettype": "abstract"}, text=True)
            if xml:
                abstracts.update(self._parse_abstracts(xml))
        out = {}
        for p in pmids:
            s = summ.get(p)
            if not s:
                continue
            doi = next((a.get("value") for a in s.get("articleids", []) if a.get("idtype") == "doi"), "")
            out[p] = {"pmid": p, "title": (s.get("title") or "").rstrip("."), "journal": s.get("fulljournalname") or s.get("source") or "",
                      "year": (s.get("pubdate") or "")[:4], "authors": ", ".join(a.get("name", "") for a in (s.get("authors") or [])[:4])
                      + (" et al." if len(s.get("authors") or []) > 4 else ""),
                      "doi": doi, "abstract": abstracts.get(p, ""), "url": self.pubmed_url(p)}
        return out

    @staticmethod
    def _parse_abstracts(xml):
        out = {}
        try:
            root = ET.fromstring(xml)
        except ET.ParseError:
            return out
        for art in root.iter("PubmedArticle"):
            pmid = art.findtext(".//MedlineCitation/PMID")
            text = " ".join("".join(t.itertext()) for t in art.findall(".//Abstract/AbstractText"))
            if pmid:
                out[pmid] = text
        return out


def classify(article, gene_names, trait_terms):
    """'title co-mention' | 'co-mention' | 'keyword only' for one article (see the module notes)."""
    def has(text, names):
        return any(re.search(r"(?<![A-Za-z0-9])" + re.escape(n) + r"(?![A-Za-z0-9])", text, re.I) for n in names if n)

    def has_trait(text):
        return any(re.search(r"(?<![a-z])" + re.escape(t.lower()), text.lower()) for t in trait_terms)

    title, abstract = article.get("title", ""), article.get("abstract", "")
    if has(title, gene_names) and has_trait(title + " " + abstract):
        return "title co-mention"
    if has(title + " " + abstract, gene_names) and has_trait(title + " " + abstract):
        return "co-mention"
    return "keyword only"
