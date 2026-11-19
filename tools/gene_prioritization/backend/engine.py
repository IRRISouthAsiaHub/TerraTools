"""
The Gene Prioritization pipeline.

  1  inputs         the gene ID list (for example the Gene_ID_list.txt of Genes Extraction)
  2  services       ping Ensembl Plants, Gramene, Plant Reactome, NCBI (and Oryzabase for rice)
  3  genes          resolve the user's identifiers
  4  annotation     coordinates, description, GO, domains, trait annotations, pathways, orthologs (cached)
  5  scoring        five evidence categories -> S(g); literature is searched for the best-ranked genes
  6  results        ranked table, per-gene evidence with sources, chart data

Nothing is invented: when a service fails, the evidence from it is marked unavailable for the genes concerned.
"""

from __future__ import annotations

import re
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

from tools.gene_prioritization import config
from tools.gene_prioritization.backend import inputs, scoring, traits
from tools.gene_prioritization.backend.cache import Cache
from tools.gene_prioritization.backend.ensembl import Ensembl
from tools.gene_prioritization.backend.gramene import Gramene
from tools.gene_prioritization.backend.http import Cancelled, NetworkError, ApiError, RunContext, today
from tools.gene_prioritization.backend.ncbi import Ncbi, classify
from tools.gene_prioritization.backend.oryzabase import Oryzabase, rapdb_url
from tools.gene_prioritization.backend.plant_reactome import PlantReactome

CATS = ["association", "function", "pathway", "ortholog", "literature"]
MSU_ID = re.compile(r"^LOC_Os\d{2}g\d{5}", re.I)


class PrioritizationError(Exception):
    pass


class Engine:
    def __init__(self, cfg, progress=None, ctx=None):
        self.cfg = cfg
        self.sp = config.SPECIES.get(cfg.get("species"))
        self._progress = progress or (lambda pct, msg, extra=None: None)
        self.cache = Cache(config.CACHE_DB)
        self.ctx = ctx or RunContext(self.cache, refresh=bool(cfg.get("refresh_cache")))
        self.ctx.cache = self.cache
        self.ctx.refresh = bool(cfg.get("refresh_cache"))
        self.ens, self.gra = Ensembl(self.ctx), Gramene(self.ctx)
        self.rea, self.ncbi = PlantReactome(self.ctx), Ncbi(self.ctx)
        self.ory = Oryzabase(self.ctx) if (self.sp or {}).get("oryzabase") else None
        self.up = {}
        self.warnings = []
        self._lock = threading.Lock()

    # ------------------------------------------------------------ small helpers
    def progress(self, pct, msg, extra=None):
        self.ctx.check()
        self._progress(int(pct), msg, extra)

    def warn(self, msg):
        with self._lock:
            if msg not in self.warnings:
                self.warnings.append(msg)

    def api_status(self):
        return [c.api.status() for c in (self.ens, self.gra, self.rea, self.ncbi, *([self.ory] if self.ory else []))]

    # ============================================================ 1. validation
    def validate(self):
        cfg = self.cfg
        if not self.sp:
            raise PrioritizationError("Choose a species.")
        if not self.sp["enabled"]:
            raise PrioritizationError(self.sp.get("reason") or f"{self.sp['label']} is not supported yet.")
        trait = (cfg.get("trait") or "").strip()
        if len(trait) < 3:
            raise PrioritizationError("Enter the trait of interest (for example Drought tolerance).")
        gene_text, gene_file = cfg.get("gene_text") or "", cfg.get("gene_file") or ""
        if not (gene_text.strip() or gene_file):
            raise PrioritizationError("Provide the genes to rank: a gene ID list file or pasted gene IDs.")
        import os
        if gene_file and not os.path.exists(gene_file):
            raise PrioritizationError(f"The gene ID list file was not found: {gene_file}")
        if gene_file and not gene_file.lower().endswith((".txt", ".csv", ".tsv", ".xlsx", ".xls", ".tab")):
            raise PrioritizationError("The gene ID list file must be TXT, CSV, TSV or XLSX.")
        return trait

    # ============================================================ 2. services
    def check_services(self):
        """Ping every service. Raises when no service can be reached."""
        up = {}
        for name, conn in (("Ensembl Plants", self.ens), ("Gramene", self.gra), ("Plant Reactome", self.rea), ("NCBI", self.ncbi),
                           *((("Oryzabase", self.ory),) if self.ory else ())):
            try:
                up[name] = bool(conn.ping())
            except Cancelled:
                raise
            except (NetworkError, ApiError):
                up[name] = False
        self.up = up
        if not any(up.values()):
            raise PrioritizationError(
                "No internet connection: none of Ensembl Plants, Gramene, Plant Reactome or NCBI could be reached. "
                "Check the connection and run again. No annotations or scores were produced.")
        for name, ok in up.items():
            if not ok:
                self.warn(f"{name} could not be reached; evidence from it is unavailable in this run.")
        return up

    # ============================================================ 3. genes
    def resolve(self, tokens):
        """-> ({input: gene id}, [{input, reason}], {gene id: Gramene document})"""
        sp = self.sp
        found, unresolved, docs = {}, [], {}

        def accept(tok, gid, doc):
            if doc.get("system_name") and doc["system_name"] != sp["gramene"]:
                unresolved.append({"input": tok, "reason": f"{gid} belongs to {doc['system_name']}, not {sp['scientific']}."})
                return
            found[tok] = gid
            docs[gid] = doc

        if self.up.get("Gramene"):
            direct, n = {}, len(tokens)
            for i in range(0, n, 40):                     # the same 40-gene requests as Gramene.genes(), with a progress report per request
                self.progress(10 + 8 * i / n, f"Looking up genes {i + 1}-{min(i + 40, n)} of {n} in Gramene…", {"api": self.api_status()})
                direct.update(self.gra.genes(tokens[i:i + 40]))
            todo = []
            for tok in tokens:
                if tok in direct:
                    accept(tok, tok, direct[tok])
                else:
                    todo.append(tok)
            for k, tok in enumerate(todo):
                self.ctx.check()
                self.progress(18 + 19 * k / len(todo), f"Searching for identifier {k + 1} of {len(todo)} that Gramene does not know as a gene ID: {tok}",
                              {"api": self.api_status()})
                gid = None
                base = re.sub(r"\.\d+$", "", tok)                            # Os01g0100100.1 -> Os01g0100100
                if base != tok:
                    d = self.gra.genes([base])
                    if base in d:
                        accept(tok, base, d[base])
                        continue
                cands = self.gra.search(tok, sp["gramene"])
                if cands:
                    d = self.gra.genes(cands)
                    exact = [i for i in cands if i in d and tok.lower() in Gramene.identifiers(d[i])]
                    if len(exact) == 1:
                        gid = exact[0]
                        accept(tok, gid, d[gid])
                        continue
                    if len(exact) > 1:
                        unresolved.append({"input": tok, "reason": f"ambiguous: matches {len(exact)} genes ({', '.join(exact[:4])}…). Use the gene ID."})
                        continue
                if self.up.get("NCBI") and (tok.isdigit() or tok.upper().startswith("LOC")):
                    gid = self._via_ncbi(tok, docs)
                    if gid:
                        accept(tok, gid, docs[gid])
                        continue
                if MSU_ID.match(tok):
                    rap = self._rap_from_msu(tok)
                    if rap:
                        d = self.gra.genes([rap])
                        if rap in d:
                            accept(tok, rap, d[rap])
                            continue
                    unresolved.append({"input": tok, "reason": "MSU (LOC_Os…) locus ID not translated: only the genes listed in Oryzabase can be mapped to "
                                                               "RAP-DB IDs. Use the RAP-DB ID (Os01g0100100) or the gene symbol."})
                else:
                    unresolved.append({"input": tok, "reason": f"not found in {sp['scientific']} (Gramene / Ensembl Plants / NCBI Gene)."})
        elif self.up.get("Ensembl Plants"):
            info = self.ens.lookup(tokens)
            for tok in tokens:
                if info.get(tok) and info[tok].get("species") == sp["ensembl"]:
                    found[tok] = tok
                else:
                    unresolved.append({"input": tok, "reason": "not an Ensembl Plants gene ID of this species (Gramene is unreachable, so symbols cannot be searched)."})
        return found, unresolved, docs

    def _rap_from_msu(self, msu):
        if not (self.ory and self.up.get("Oryzabase")):
            return None
        try:
            return self.ory.rap_from_msu(msu)
        except NetworkError:
            self.up["Oryzabase"] = False
            self.warn("Oryzabase stopped answering; its evidence is unavailable in this run.")
            return None

    def _via_ncbi(self, tok, docs):
        """NCBI Gene ID / LOC number -> the species' own gene ID, through the other names NCBI lists for the gene."""
        num = re.sub(r"^LOC", "", tok, flags=re.I)
        summ = self.ncbi.summary("gene", [num]).get(num)
        if not summ:
            return None
        names = re.split(r"[|,]\s*", f"{summ.get('otherdesignations', '')}|{summ.get('otheraliases', '')}|{summ.get('name', '')}")
        names = [n.strip() for n in names if re.fullmatch(r"[A-Za-z0-9_.\-]{5,}", n.strip())]
        d = self.gra.genes(names)
        for n in names:
            if n in d and d[n].get("system_name") == self.sp["gramene"]:
                docs[n] = d[n]
                return n
        return None

    # ============================================================ 4. annotation
    def annotate(self, gene_ids, gdocs):
        """Annotation record per gene (cached). Returns {gene id: record}."""
        sp = self.sp
        records, todo = {}, []
        for gid in gene_ids:
            rec = None if self.ctx.refresh else self.cache.get_annotation(sp["ensembl"], gid, config.CACHE_DAYS)
            if rec:
                records[gid] = rec
            else:
                todo.append(gid)
        self.n_cached = len(records)
        self._load_oryzabase()
        for rec in records.values():                     # records cached before Oryzabase was added get its entry now
            if self.ory and self.up.get("Oryzabase") and "oryzabase_entry" not in rec:
                self._add_oryzabase(rec)
                self.cache.put_annotation(sp["ensembl"], rec["id"], rec)
        if not todo:
            return records
        self.progress(40, f"Retrieving annotations for {len(todo)} gene(s)…")
        try:
            ens = self.ens.lookup(todo) if self.up.get("Ensembl Plants") else {}
        except NetworkError:
            ens = {}
        need_docs = [g for g in todo if g not in gdocs]
        if need_docs and self.up.get("Gramene"):
            try:
                gdocs.update(self.gra.genes(need_docs))
            except NetworkError:
                pass

        done, total = [0], len(todo)

        def one(gid):
            self.ctx.check()
            rec = self._build_record(gid, ens.get(gid), gdocs.get(gid))
            with self._lock:
                done[0] += 1
                self.progress(40 + 30 * done[0] / total, f"Annotated {done[0]} of {total} genes", None)
            return gid, rec

        with ThreadPoolExecutor(max_workers=config.WORKERS) as ex:
            futs = [ex.submit(one, g) for g in todo]
            try:
                for f in as_completed(futs):
                    gid, rec = f.result()
                    records[gid] = rec
            except Cancelled:
                for f in futs:
                    f.cancel()
                raise
        self._fill_orthologs([records[g] for g in todo])
        for gid in todo:
            rec = records[gid]
            if not rec["incomplete"]:
                self.cache.put_annotation(sp["ensembl"], gid, rec)
        return records

    def _load_oryzabase(self):
        """Read the Oryzabase gene list once (downloaded the first time, then cached). A failure only removes that evidence."""
        if not (self.ory and self.up.get("Oryzabase")):
            return
        self.progress(39, "Reading the Oryzabase gene list…")
        try:
            self.ory.load()
        except NetworkError:
            self.up["Oryzabase"] = False
            self.warn("The Oryzabase gene list could not be retrieved; evidence from Oryzabase is unavailable in this run.")

    def _add_oryzabase(self, rec):
        try:
            e = self.ory.entry(rec["id"], rec.get("msu_id", ""))
        except NetworkError:
            rec["incomplete"] = True
            rec["unavailable"]["Oryzabase"] = "stopped answering"
            return
        rec["oryzabase_entry"] = e
        if e:
            rec["provenance"].append(self.ory.prov(e))
            if not rec.get("symbol") and e["symbol"]:
                rec["symbol"] = e["symbol"]

    def _build_record(self, gid, einfo, gdoc):
        sp = self.sp
        rec = {"id": gid, "species": sp["scientific"], "symbol": "", "description": "", "chromosome": None, "start": None,
               "end": None, "strand": None, "biotype": "", "go": [], "domains": [], "to": [], "qtl_to": [], "pathways": [],
               "orthologs": [], "phenotypes": [], "uniprot": [], "entrez": [], "msu_id": "", "provenance": [],
               "unavailable": {}, "incomplete": False, "retrieved": today(),
               "ensembl_url": Ensembl.gene_url(sp["ensembl"], gid), "gramene_url": Gramene.gene_url(gid)}
        if einfo:
            rec.update({"chromosome": einfo["chromosome"], "start": einfo["start"], "end": einfo["end"], "strand": einfo["strand"],
                        "biotype": einfo["biotype"] or "", "symbol": einfo["symbol"] or "", "description": einfo["description"] or ""})
            rec["assembly"] = einfo.get("assembly")
            rec["provenance"].append(self.ens.prov(rec["ensembl_url"], "gene coordinates and identifiers"))
        elif self.up.get("Ensembl Plants"):
            rec["unavailable"]["Ensembl Plants"] = "gene not found"
        if gdoc:
            p = Gramene.parse(gdoc)
            for k in ("go", "domains", "to", "qtl_to", "uniprot", "entrez", "msu_id"):
                rec[k] = p[k]
            if not rec["description"]:
                rec["description"] = p["description"]
            if not rec["symbol"] and p["name"] and p["name"] != gid:
                rec["symbol"] = p["name"]
            if rec["chromosome"] is None:
                rec.update({"chromosome": p["chromosome"], "start": p["start"], "end": p["end"], "strand": p["strand"], "biotype": p["biotype"]})
            rec["gramene_pathways"] = p["pathways"]
            rec["closest_ortholog"] = p["closest_ortholog"]
            rec["oryzabase"] = p["oryzabase"]
            rec["provenance"].append(self.gra.prov(gid, "annotation: description, GO, domains, trait terms, orthologs"))
        elif self.up.get("Gramene"):
            rec["unavailable"]["Gramene"] = "gene not found"
        # ---- Ensembl: GO cross-references, phenotypes, orthologs
        if self.up.get("Ensembl Plants"):
            try:
                have = {g["id"] for g in rec["go"]}
                for t in self.ens.go_terms(gid):
                    if t["id"] not in have:
                        rec["go"].append(t)
                rec["phenotypes"] = self.ens.phenotypes(sp["ensembl"], gid)
                rec["orthologs"] = [{**o, "source": "Ensembl Plants", "url": ""} for o in
                                    self.ens.orthologs(sp["ensembl"], gid, sp["orthologs_in"])][:40]
                rec["provenance"].append(self.ens.prov(f"{config.ENSEMBL_URL}/homology/id/{sp['ensembl']}/{gid}", "orthologs"))
            except NetworkError:
                rec["incomplete"] = True
                rec["unavailable"]["Ensembl Plants"] = "stopped answering"
        # ---- Oryzabase (rice): curated gene symbol/name, trait classes, Trait Ontology
        if self.ory and self.up.get("Oryzabase"):
            self._add_oryzabase(rec)
        # ---- Plant Reactome
        if self.up.get("Plant Reactome") and sp["reactome"]:
            try:
                idents = [*rec["uniprot"], rec["msu_id"], gid, rec["symbol"]]
                st = self.rea.gene_entity(sp["reactome"], idents)
                pws = self.rea.pathways_of(st, sp["reactome_taxon"]) if st else []
                rec["reactome_entity"] = st
                have = {p["id"] for p in pws}
                for p in rec.get("gramene_pathways", []):                      # the same pathways, as listed on the Gramene record
                    if p["id"] not in have:
                        pws.append({"id": p["id"], "name": p["name"], "url": self.rea.pathway_url(p["id"])})
                rec["pathways"] = pws
                if pws:
                    rec["provenance"].append(self.rea.prov(pws[0]["id"], "pathway membership"))
            except NetworkError:
                rec["incomplete"] = True
                rec["unavailable"]["Plant Reactome"] = "stopped answering"
        else:
            rec["pathways"] = [{"id": p["id"], "name": p["name"], "url": self.rea.pathway_url(p["id"])} for p in rec.get("gramene_pathways", [])]
        return rec

    def _fill_orthologs(self, recs):
        """Names and descriptions of the orthologs, one batch lookup; plus the closest ortholog Gramene lists."""
        ids = sorted({o["id"] for r in recs for o in r["orthologs"]})
        info = {}
        if ids and self.up.get("Ensembl Plants"):
            try:
                info = self.ens.lookup(ids)
            except NetworkError:
                self.warn("Ortholog names could not be retrieved from Ensembl Plants.")
        arabidopsis = "arabidopsis_thaliana"
        for r in recs:
            ors, seen = [], set()
            rank = {"ortholog_one2one": 0, "ortholog_one2many": 1}
            for o in sorted(r["orthologs"], key=lambda o: (o["species"] != arabidopsis, rank.get(o["type"], 2))):
                i = info.get(o["id"]) or {}
                o = {**o, "name": i.get("symbol", ""), "description": i.get("description", ""),
                     "url": Ensembl.gene_url(o["species"], o["id"])}
                if o["id"] not in seen:
                    seen.add(o["id"])
                    ors.append(o)
            c = r.get("closest_ortholog")
            if c and c["id"] not in seen:
                ors.append({"id": c["id"], "species": "closest representative in the Gramene gene tree", "type": "closest", "name": c["name"],
                            "description": c["description"], "source": "Gramene", "url": Gramene.gene_url(c["id"])})
            elif c:
                for o in ors:
                    if o["id"] == c["id"] and not o.get("description"):
                        o["description"] = c["description"]
                    if o["id"] == c["id"] and not o.get("name"):
                        o["name"] = c["name"]
            r["orthologs"] = ors[:14]
            # phenotype annotations of the two best Arabidopsis orthologs (Ensembl Plants)
            if self.up.get("Ensembl Plants"):
                for o in [x for x in r["orthologs"] if x.get("species") == arabidopsis][:2]:
                    try:
                        o["phenotypes"] = self.ens.phenotypes(arabidopsis, o["id"])
                    except NetworkError:
                        break

    # ============================================================ 5. literature
    def literature_for(self, genes, recs, tt):
        """{gene id: [paper dicts]} for the genes given."""
        sp = self.sp
        words = [f'"{w}"' for w in traits.search_words(tt, 5)]
        trait_q = "(" + " OR ".join(words) + ")"
        species_q = f'("{sp["scientific"]}"[Organism] OR {sp["common"]})'
        per_gene, all_pmids = {}, []
        for n, gid in enumerate(genes, 1):
            self.progress(80 + 8 * n / max(len(genes), 1), f"Searching PubMed for gene {n} of {len(genes)}")
            r = recs[gid]
            oz = r.get("oryzabase_entry") or {}
            names = [x for x in dict.fromkeys([gid, r.get("symbol", ""), *r.get("oryzabase", []), oz.get("symbol", "")]) if x and len(x) >= 3]
            if not names:
                continue
            gene_q = "(" + " OR ".join(f'"{x}"' for x in names) + ")"
            try:
                pmids, _ = self.ncbi.search("pubmed", f"{gene_q} AND {species_q} AND {trait_q}", retmax=config.PUBS_PER_GENE, sort="relevance")
            except NetworkError:
                self.warn("NCBI stopped answering; literature evidence is missing for some genes.")
                break
            per_gene[gid] = (names, pmids)
            all_pmids += pmids
        arts = {}
        if all_pmids:
            self.progress(89, "Reading PubMed titles and abstracts")
            try:
                arts = self.ncbi.articles(all_pmids)
            except NetworkError:
                self.warn("NCBI stopped answering; paper details are missing.")
        terms = tt["direct"] + tt["related"]
        out = {}
        for gid, (names, pmids) in per_gene.items():
            pubs = []
            for p in pmids:
                a = arts.get(p)
                if a:
                    pubs.append({**a, "relation": classify(a, names, terms)})
            out[gid] = pubs
        return out

    # ============================================================ main
    def run(self):
        cfg, sp = self.cfg, self.sp
        trait = self.validate()

        self.progress(2, "Reading the gene list…")
        try:
            list_ids, notes = inputs.read_gene_ids(cfg.get("gene_file") or "", cfg.get("gene_text") or "")
        except inputs.InputError as e:
            raise PrioritizationError(str(e))
        except Exception as e:
            raise PrioritizationError(f"The gene list could not be read: {e}")
        for n in notes:
            self.warn(n)
        if not list_ids:
            raise PrioritizationError("The gene ID list is empty.")

        self.progress(6, "Checking the databases…")
        self.check_services()
        self.progress(10, "Resolving gene identifiers…", {"api": self.api_status()})

        found, unresolved, gdocs = self.resolve(list_ids)
        if unresolved:
            self.warn(f"{len(unresolved)} input identifier(s) could not be matched to a {sp['scientific']} gene (see the list below the summary).")
        candidates = list(dict.fromkeys(found[tok] for tok in list_ids if tok in found))
        if not candidates:
            raise PrioritizationError("No gene could be found: none of the input identifiers matched a gene.")

        gene_ids = candidates
        self.progress(38, f"{len(gene_ids)} candidate gene(s)")
        recs = self.annotate(gene_ids, gdocs)

        # ---- trait vocabulary and Plant Reactome pathways for the trait
        tt = traits.expand(trait)
        if tt["note"]:
            self.warn(tt["note"])
        trait_pw = {}
        if self.up.get("Plant Reactome") and sp["reactome"]:
            self.progress(72, "Searching Plant Reactome for trait pathways…")
            try:
                trait_pw = self.rea.trait_pathways(sp["reactome"], traits.search_words(tt, 4))
            except NetworkError:
                self.warn("Plant Reactome stopped answering; pathway relevance uses pathway names only.")

        # ---- scoring: four categories first, then literature for the best-ranked genes
        self.progress(74, "Scoring the evidence…")
        cats = {}
        for gid in gene_ids:
            r = recs[gid]
            cats[gid] = {"association": scoring.association(r, tt), "function": scoring.function(r, tt),
                         "pathway": scoring.pathway(r, tt, trait_pw), "ortholog": scoring.ortholog(r, tt),
                         "literature": scoring.literature([])}
            cats[gid]["literature"]["note"] = "PubMed was not searched for this gene (only the best-ranked genes are searched)."

        def prelim(g):
            s = scoring.combine(cats[g])[0]
            return s if s is not None else -1

        top = sorted(gene_ids, key=prelim, reverse=True)[:config.LITERATURE_TOP_N]
        pubs = {}
        if self.up.get("NCBI"):
            pubs = self.literature_for(top, recs, tt)
            for gid in top:
                if gid in pubs:
                    cats[gid]["literature"] = scoring.literature(pubs[gid])
                else:
                    cats[gid]["literature"]["note"] = "The PubMed search failed or the gene has no usable name."
        else:
            for gid in gene_ids:
                cats[gid]["literature"]["note"] = "NCBI could not be reached."

        self.progress(92, "Ranking the genes…")
        result = self._assemble(trait, tt, gene_ids, recs, cats, pubs, unresolved, list_ids)
        self.progress(97, "Writing the result files…")
        return result

    # ============================================================ results
    def _gene_links(self, r):
        links = {"Ensembl Plants": r["ensembl_url"], "Gramene": r["gramene_url"]}
        if self.sp.get("oryzabase"):
            links["RAP-DB"] = rapdb_url(r["id"])
            e = r.get("oryzabase_entry")
            if e:
                links["Oryzabase"] = e["url"]
        return links

    @staticmethod
    def _oryzabase_detail(r):
        e = r.get("oryzabase_entry")
        if not e:
            return None
        return {"symbol": e["symbol"], "name": e["name"], "synonyms": e["synonyms"][:8], "explanation": e["explanation"],
                "classes": e["classes"], "to": e["to"][:40], "url": e["url"]}

    def _assemble(self, trait, tt, gene_ids, recs, cats, pubs, unresolved, list_ids):
        sp = self.sp
        rows, detail = [], {}
        for gid in gene_ids:
            r, c = recs[gid], cats[gid]
            final, contrib, used = scoring.combine(c)
            row = {"Gene_ID": gid, "Symbol": r["symbol"], "Description": r["description"],
                   "Final_Score": final, "Evidence_Level": scoring.evidence_level(c),
                   "Chromosome": r["chromosome"], "Start": r["start"], "End": r["end"], "Strand": r["strand"],
                   "Categories_Used": len(used)}
            for k in CATS:
                row[f"E_{k}"] = c[k]["score"]
            rows.append(row)
            detail[gid] = {"id": gid, "symbol": r["symbol"], "description": r["description"], "species": r["species"],
                           "location": f"{r['chromosome']}:{r['start']}-{r['end']} ({'+' if r['strand'] == 1 else '-' if r['strand'] == -1 else '?'})"
                                       if r["chromosome"] else "", "assembly": r.get("assembly", ""), "biotype": r["biotype"],
                           "retrieved": r["retrieved"], "links": self._gene_links(r),
                           "oryzabase_entry": self._oryzabase_detail(r),
                           "uniprot": r["uniprot"], "ncbi_gene": r["entrez"],
                           "go": r["go"][:80], "domains": r["domains"], "to": r["to"], "qtl_to": r["qtl_to"], "pathways": r["pathways"],
                           "orthologs": [{k: v for k, v in o.items() if k != "phenotypes"} for o in r["orthologs"]],
                           "phenotypes": r["phenotypes"], "publications": [{k: v for k, v in p.items() if k != "abstract"} for p in pubs.get(gid, [])],
                           "scores": {k: {"score": c[k]["score"], "available": c[k]["available"], "note": c[k]["note"], "items": c[k]["items"][:10],
                                          "weight": config.EVIDENCE_WEIGHTS[k][1], "label": config.EVIDENCE_WEIGHTS[k][0],
                                          "contribution": contrib[k]} for k in CATS},
                           "final": final, "evidence_level": row["Evidence_Level"], "used": used,
                           "provenance": r["provenance"], "unavailable": r["unavailable"]}
        rows.sort(key=lambda x: (-(x["Final_Score"] or 0), x["Gene_ID"]))
        for i, row in enumerate(rows, 1):
            row["Rank"] = i
            detail[row["Gene_ID"]]["rank"] = i

        summary = {
            "input_genes": len(set(list_ids)), "candidate_genes": len(gene_ids),
            "annotated": sum(1 for g in gene_ids if recs[g]["description"] or recs[g]["go"] or recs[g]["domains"] or recs[g]["pathways"] or recs[g]["chromosome"]),
            "with_association": sum(1 for g in gene_ids if cats[g]["association"]["available"]),
            "with_literature": sum(1 for g in gene_ids if cats[g]["literature"]["available"]),
            "with_pathway": sum(1 for g in gene_ids if cats[g]["pathway"]["available"]),
            "unresolved": len(unresolved), "from_cache": getattr(self, "n_cached", 0),
        }
        top = [r for r in rows if detail[r["Gene_ID"]]["pathways"]][:25]      # the best-ranked genes that belong to a pathway
        charts = {
            "bars": [{"gene": r["Symbol"] or r["Gene_ID"], "id": r["Gene_ID"], "final": r["Final_Score"]} for r in rows[:20]],
            "contributions": [{"gene": r["Symbol"] or r["Gene_ID"], "id": r["Gene_ID"],
                               **{k: detail[r["Gene_ID"]]["scores"][k]["contribution"] for k in CATS}}
                              for r in sorted((x for x in rows if x["Final_Score"]), key=lambda x: -x["Final_Score"])[:20]],
            "pathways": _pathway_matrix([detail[r["Gene_ID"]] for r in top]),
        }
        return {
            "status": "success", "species": sp["label"], "species_key": self.cfg["species"], "trait": trait,
            "trait_terms": {"direct": tt["direct"], "related": tt["related"], "note": tt["note"]},
            "rank_by": "Final score",
            "summary": summary, "genes": rows, "detail": detail, "charts": charts, "unresolved": unresolved,
            "warnings": self.warnings, "api_status": self.api_status(),
            "weights": {k: {"label": v[0], "weight": v[1]} for k, v in config.EVIDENCE_WEIGHTS.items()},
            "run_date": today(), "ncbi_key": bool(getattr(self.ncbi, "has_key", False)),
        }


def _pathway_matrix(details):
    """Gene x pathway membership of the top genes (only pathways shared by or specific to them), for the chart."""
    names = {}
    for d in details:
        for p in d["pathways"]:
            names.setdefault(p["id"], p["name"])
    count = {pid: sum(any(p["id"] == pid for p in d["pathways"]) for d in details) for pid in names}
    pids = sorted(names, key=lambda p: (-count[p], names[p]))[:14]
    return {"pathways": [{"id": p, "name": names[p]} for p in pids],
            "genes": [{"id": d["id"], "gene": d["symbol"] or d["id"], "in": [p for p in pids if any(x["id"] == p for x in d["pathways"])]}
                      for d in details if any(x["id"] in pids for x in d["pathways"])]}
