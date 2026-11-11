"""
Gramene REST API (https://data.gramene.org, documented at https://data.gramene.org/swagger). Endpoints used:

  GET /genes?idList=a,b,c   gene documents: description, location, cross-references, GO terms (with evidence codes),
                            protein domains, QTL trait-ontology terms (genes inside QTL), homology (orthologs, closest
                            representative)
  GET /search?q=...&fq=...  full-text gene search (Solr): resolves symbols, UniProt accessions and NCBI Gene IDs to genes

Two kinds of trait terms (Trait Ontology, TO) sit on a gene record:
  annotations.TO      the gene itself is annotated with the trait term (a curated gene-trait association)
  annotations.QTL_TO  the gene lies inside a QTL annotated with the term (positional, not a gene-level result)
The tool keeps them apart. Plant Reactome pathways of the gene are also listed in annotations.pathways.
"""

from __future__ import annotations

from tools.gene_prioritization import config
from tools.gene_prioritization.backend.http import ApiClient, today

SOURCE = "Gramene"
GENE_FIELDS = "_id,name,description,biotype,system_name,location,xrefs,annotations,homology,MSU_id"


class Gramene:
    def __init__(self, ctx):
        self.api = ApiClient(SOURCE, config.GRAMENE_URL, ctx, min_interval=0.1)

    @staticmethod
    def gene_url(gene_id):
        return f"https://www.gramene.org/genes?q={gene_id}"

    def prov(self, gene_id, kind):
        return {"source": SOURCE, "url": self.gene_url(gene_id), "retrieved": today(), "type": kind}

    def ping(self):
        r = self.api.request("/genes", params={"idList": "Os01g0100100", "fl": "_id"}, ttl=0, retries=2, timeout=12)
        return bool(r)

    # ------------------------------------------------------------ genes
    def genes(self, gene_ids):
        """{id: gene document} for every id Gramene knows."""
        out = {}
        ids = list(dict.fromkeys(gene_ids))
        for i in range(0, len(ids), 40):
            chunk = ids[i:i + 40]
            r = self.api.request("/genes", params={"idList": ",".join(chunk), "fl": GENE_FIELDS}) or []
            for d in r:
                out[d["_id"]] = d
        return out

    def search(self, term, system_name, rows=6):
        """Gene ids whose record matches `term` (a symbol, UniProt accession, NCBI Gene ID …) in one species."""
        r = self.api.request("/search", params={"q": f'"{term}"', "fq": f"system_name:{system_name}", "fl": "id,name",
                                                "rows": rows, "wt": "json"})
        if not r:
            return []
        return [d["id"] for d in r.get("response", {}).get("docs", [])]

    # ------------------------------------------------------------ parsing helpers
    @staticmethod
    def xref(doc, db):
        for x in doc.get("xrefs", []):
            if x.get("db") == db:
                return list(x.get("ids", []))
        return []

    @staticmethod
    def identifiers(doc):
        """Every identifier the record carries, lower-cased: used to confirm that a search hit really is the input."""
        ids = {str(doc.get("_id", "")).lower(), str(doc.get("name", "")).lower(), str(doc.get("MSU_id", "")).lower()}
        for x in doc.get("xrefs", []):
            ids.update(str(i).lower() for i in x.get("ids", []))
        ids.discard("")
        return ids

    @staticmethod
    def parse(doc):
        """Evidence of one Gramene gene document in the tool's own form."""
        ann = doc.get("annotations", {}) or {}
        go = [{"id": e.get("id"), "name": e.get("name"), "namespace": e.get("namespace") or "",
               "definition": e.get("def") or "", "evidence_code": e.get("evidence_code") or "", "source": SOURCE}
              for e in (ann.get("GO", {}) or {}).get("entries", []) if e.get("id")]
        domains = [{"id": e.get("id"), "name": e.get("name") or "", "description": e.get("description") or ""}
                   for e in (ann.get("domains", {}) or {}).get("entries", []) if e.get("id")]
        def terms(key):
            return [{"id": e.get("id"), "name": e.get("name") or "", "definition": e.get("def") or "",
                     "evidence_code": e.get("evidence_code") or ""}
                    for e in (ann.get(key, {}) or {}).get("entries", []) if e.get("id")]

        pathways = [{"id": e.get("id"), "name": e.get("name") or "", "url": e.get("link") or ""}
                    for e in (ann.get("pathways", {}) or {}).get("entries", []) if e.get("id")]
        hom = doc.get("homology", {}) or {}
        rep = ((hom.get("gene_tree") or {}).get("representative") or {})
        closest = rep.get("closest") or rep.get("model") or {}
        loc = doc.get("location") or {}
        return {
            "name": doc.get("name") or "", "description": (doc.get("description") or "").split(" [Source")[0],
            "biotype": doc.get("biotype") or "", "system_name": doc.get("system_name") or "",
            "chromosome": loc.get("region"), "start": loc.get("start"), "end": loc.get("end"), "strand": loc.get("strand"),
            "msu_id": doc.get("MSU_id") or "",
            "uniprot": [*Gramene.xref(doc, "Uniprot/SWISSPROT"), *Gramene.xref(doc, "Uniprot/SPTREMBL")],
            "entrez": Gramene.xref(doc, "EntrezGene"), "oryzabase": Gramene.xref(doc, "Oryzabase"),
            "go": go, "domains": domains, "to": terms("TO"), "qtl_to": terms("QTL_TO"), "pathways": pathways,
            "closest_ortholog": ({"id": closest.get("id"), "name": closest.get("name") or "",
                                  "description": closest.get("description") or "",
                                  "identity": closest.get("percent_identity")} if closest.get("id") else None),
        }
