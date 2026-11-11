"""
Plant Reactome ContentService (https://plantreactome.gramene.org/ContentService). Endpoints used, checked live:

  GET /search/query?query=...&species=...&types=Pathway&cluster=true   pathways whose name or summary matches a word
  GET /search/query?query=<gene identifier>&species=...                  the gene product's entry (stable id R-OSA-…)
  GET /data/pathways/low/entity/{stId}?species={taxon}                   the lowest-level pathways that contain it

The gene product is found by a UniProt accession or by a locus name (MSU or RAP ID, symbol), whichever the record has.
The site https://plantreactome.org does not answer API calls; the ContentService host above does.
"""

from __future__ import annotations

import re

from tools.gene_prioritization import config
from tools.gene_prioritization.backend.http import ApiClient, today

SOURCE = "Plant Reactome"
TAG = re.compile(r"<[^>]+>")


def clean(text):
    return TAG.sub("", text or "").strip()


class PlantReactome:
    def __init__(self, ctx):
        self.api = ApiClient(SOURCE, config.REACTOME_URL, ctx, min_interval=0.1)

    @staticmethod
    def pathway_url(st_id):
        return f"https://plantreactome.gramene.org/content/detail/{st_id}"

    def prov(self, st_id, kind="pathway"):
        return {"source": SOURCE, "url": self.pathway_url(st_id), "retrieved": today(), "type": kind}

    def ping(self):
        r = self.api.request("/data/species/main", ttl=0, retries=2, timeout=12)
        return bool(r)

    # ------------------------------------------------------------ pathways about a trait
    def trait_pathways(self, species_name, words):
        """{stId: {name, summary, query}} of pathways whose name or summary mentions one of the words."""
        found = {}
        for w in words:
            r = self.api.request("/search/query", params={"query": w, "species": species_name, "types": "Pathway",
                                                          "cluster": "true"})
            for grp in (r or {}).get("results", []):
                for e in grp.get("entries", []):
                    if e.get("exactType") != "Pathway" or not e.get("stId"):
                        continue
                    found.setdefault(e["stId"], {"name": clean(e.get("name")), "summary": clean(e.get("summation")),
                                                 "query": w})
        return found

    # ------------------------------------------------------------ pathways of a gene
    def gene_entity(self, species_name, identifiers):
        """stId of the gene product's entry, searching each identifier in turn. None when Plant Reactome has none."""
        for ident in identifiers:
            if not ident:
                continue
            r = self.api.request("/search/query", params={"query": ident, "species": species_name, "cluster": "true"})
            for grp in (r or {}).get("results", []):
                for e in grp.get("entries", []):
                    if e.get("exactType") not in ("ReferenceGeneProduct", "ReferenceIsoform") or not e.get("stId"):
                        continue
                    ref = clean(e.get("referenceIdentifier")).lower()
                    names = f"{e.get('referenceName') or ''} {clean(e.get('name'))}".lower()
                    if ident.lower() == ref or ident.lower() in re.split(r"[\s,;]+", names):
                        return e["stId"]
        return None

    def pathways_of(self, st_id, taxon):
        r = self.api.request(f"/data/pathways/low/entity/{st_id}", params={"species": taxon}) or []
        return [{"id": p.get("stId"), "name": p.get("displayName") or "", "url": self.pathway_url(p.get("stId"))}
                for p in r if p.get("stId")]
