"""
Ensembl Plants REST API (https://rest.ensembl.org/documentation/). Endpoints used, each checked against the live service:

  GET  /info/ping                                  service up?
  GET  /info/assembly/{species}                    chromosome names and lengths, to map the input's chromosome names
  POST /lookup/id                                  gene coordinates, biotype, description (batch)
  GET  /overlap/region/{species}/{region}          genes overlapping a region (GWAS window, QTL interval)
  GET  /homology/id/{species}/{id}                 orthologues in other species
  GET  /xrefs/id/{id}?external_db=GO               GO cross-references
  GET  /phenotype/gene/{species}/{id}              gene-to-phenotype annotations (empty for most plant genes)
"""

from __future__ import annotations

from tools.gene_prioritization import config
from tools.gene_prioritization.backend.http import ApiClient, today

SOURCE = "Ensembl Plants"
SKIP_BIOTYPES = {"tRNA", "rRNA", "snoRNA", "snRNA", "miRNA", "Mt_tRNA", "Mt_rRNA", "misc_RNA", "scaRNA"}


class Ensembl:
    def __init__(self, ctx):
        self.api = ApiClient(SOURCE, config.ENSEMBL_URL, ctx, min_interval=0.08,
                             extra_params={"content-type": "application/json"})

    # ------------------------------------------------------------ links
    @staticmethod
    def gene_url(species, gene_id):
        return f"https://plants.ensembl.org/{species}/Gene/Summary?g={gene_id}"

    def prov(self, url, kind):
        return {"source": SOURCE, "url": url, "retrieved": today(), "type": kind}

    # ------------------------------------------------------------ service
    def ping(self):
        r = self.api.request("/info/ping", ttl=0, retries=3, timeout=30)      # the service is sometimes slow to answer; a false "down" aborts the run
        return bool(r and r.get("ping") == 1)

    def chromosomes(self, species):
        """{region name: length}. None when Ensembl Plants does not know the species."""
        r = self.api.request(f"/info/assembly/{species}", ttl=30)
        if not r:
            return None
        return {x["name"]: x["length"] for x in r.get("top_level_region", []) if x.get("coord_system") == "chromosome"}

    # ------------------------------------------------------------ genes
    def lookup(self, gene_ids):
        """{gene id: info or None}. Info keeps the coordinates, biotype and description."""
        out = {}
        ids = list(dict.fromkeys(gene_ids))
        for i in range(0, len(ids), 200):
            chunk = ids[i:i + 200]
            r = self.api.request("/lookup/id", body={"ids": chunk}) or {}
            for gid in chunk:
                info = r.get(gid)
                out[gid] = None if not info else {
                    "id": info.get("id", gid), "species": info.get("species"), "chromosome": info.get("seq_region_name"),
                    "start": info.get("start"), "end": info.get("end"), "strand": info.get("strand"),
                    "biotype": info.get("biotype"), "symbol": info.get("display_name") or "",
                    "description": (info.get("description") or "").split(" [Source")[0],
                    "assembly": info.get("assembly_name"), "source": info.get("source")}
        return out

    def genes_in_region(self, species, chromosome, start, end):
        """Genes that overlap chromosome:start-end (Ensembl allows regions up to 5 Mb per call)."""
        genes, step = [], 4_900_000
        s = int(max(1, start))
        end = int(end)
        while s <= end:
            e = min(end, s + step - 1)
            # a region overlap is a heavy query for Ensembl: it can need far longer than the default timeout when the service is busy
            r = self.api.request(f"/overlap/region/{species}/{chromosome}:{s}-{e}", params={"feature": "gene"}, timeout=120) or []
            for g in r:
                if g.get("biotype") in SKIP_BIOTYPES:
                    continue
                genes.append({"id": g["id"], "chromosome": g.get("seq_region_name"), "start": g.get("start"),
                              "end": g.get("end"), "strand": g.get("strand"), "biotype": g.get("biotype"),
                              "symbol": g.get("external_name") or "", "description": g.get("description") or ""})
            s = e + 1
        return list({g["id"]: g for g in genes}.values())

    # ------------------------------------------------------------ evidence
    def orthologs(self, species, gene_id, targets):
        """[{id, species, type, taxonomy_level}] in the target species."""
        params = {"type": "orthologues", "format": "condensed"}
        url = f"/homology/id/{species}/{gene_id}"
        # several target_species values: requests repeats the key for a list
        params["target_species"] = list(targets)
        r = self.api.request(url, params=params)
        if not r or not r.get("data"):
            return []
        out = []
        for h in r["data"][0].get("homologies", []):
            out.append({"id": h["id"], "species": h.get("species"), "type": h.get("type"),
                        "taxonomy_level": h.get("taxonomy_level")})
        return out

    def go_terms(self, gene_id):
        r = self.api.request(f"/xrefs/id/{gene_id}", params={"external_db": "GO", "all_levels": 0}) or []
        return [{"id": x.get("primary_id"), "name": x.get("description") or x.get("display_id"),
                 "namespace": "", "evidence_code": "", "source": SOURCE}
                for x in r if x.get("primary_id")]

    def phenotypes(self, species, gene_id):
        """Phenotype annotations of a gene: [{description, source, accession}]."""
        r = self.api.request(f"/phenotype/gene/{species}/{gene_id}") or []
        out = []
        for x in r:
            for p in x.get("phenotypes_associated", [x]):
                d = p.get("description")
                if d:
                    out.append({"description": d, "source": p.get("source") or "", "accession": p.get("ontology_accessions") or []})
        return out
