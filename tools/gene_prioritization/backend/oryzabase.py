"""
Oryzabase (NBRP-Rice, National Institute of Genetics, Japan: https://shigen.nig.ac.jp/rice/oryzabase) — curated rice genes
keyed by their RAP-DB identifier.

Oryzabase is the curated gene catalogue behind the RAP-DB gene names: for every characterised gene it lists the CGSNL gene
symbol and name, the trait classes, Trait Ontology and Gene Ontology terms, and the RAP-DB and MSU identifiers. It is used here
as one more evidence source (known gene-trait association, gene name) and to translate MSU (LOC_Os...) identifiers.

RAP-DB's own web server (rapdb.dna.affrc.go.jp) is not queried: its HTTPS certificate does not verify.

Nothing is stored in the application. The gene list (about 11 MB of text, downloaded once) is reduced to the genes that have a RAP-DB
identifier and kept in the user's result cache for CACHE_DAYS days, like every other response.

  GET https://shigen.nig.ac.jp/rice/oryzabase/gene/download?classtag=GENE_EN_LIST     tab-separated gene list
  https://shigen.nig.ac.jp/rice/oryzabase/gene/detail/<Trait Gene Id>                  page of one gene
"""

from __future__ import annotations

import csv
import io
import re
import threading

import requests

from tools.gene_prioritization import config
from tools.gene_prioritization.backend.http import ApiClient, ApiError, NetworkError, today

SOURCE = "Oryzabase"
BASE = "https://shigen.nig.ac.jp"
LIST_URL = BASE + "/rice/oryzabase/gene/download?classtag=GENE_EN_LIST"
CACHE_KEY = "oryzabase:gene_list:v1"
RAP_ID = re.compile(r"Os\d{2}g\d{7}", re.I)
TERM = re.compile(r"((?:TO|GO):\d+) - (.*?)(?=, (?:TO|GO):\d+ - |$)")


def gene_url(trait_gene_id):
    return f"{BASE}/rice/oryzabase/gene/detail/{trait_gene_id}"


def rapdb_url(rap_id):
    return f"https://rapdb.dna.affrc.go.jp/viewer/gene_detail/irgsp1?name={rap_id}"


def msu_key(msu_id):
    """LOC_Os01g01010.1 / loc_os01g01010 -> LOC_Os01g01010"""
    return re.sub(r"\.\d+$", "", (msu_id or "").strip()).upper().replace("LOC_OS", "LOC_Os")


def _clean(v):
    v = (v or "").strip()
    return "" if v in ("_", "-") else v


def _terms(text):
    return [{"id": i, "name": n.strip()} for i, n in TERM.findall(text or "")]


def parse_list(raw):
    """Tab-separated Oryzabase gene list -> {"genes": {RAP ID: entry}, "msu": {MSU ID: RAP ID}} for the genes that have a RAP ID."""
    rows = csv.reader(io.StringIO(raw), delimiter="\t")
    head = next(rows)
    col = {name: i for i, name in enumerate(head)}
    need = ("Trait Gene Id", "CGSNL Gene Symbol", "CGSNL Gene Name", "Explanation", "Trait Class", "RAP ID", "MSU ID", "Trait Ontology")
    missing = [n for n in need if n not in col]
    if missing:
        raise ApiError(f"Oryzabase: the gene list has no column {missing[0]!r} (the file format changed).")
    genes, msu = {}, {}
    for r in rows:
        if len(r) < len(head):
            continue
        raps = list(dict.fromkeys(RAP_ID.findall(r[col["RAP ID"]])))
        if not raps:
            continue
        symbol, name = _clean(r[col["CGSNL Gene Symbol"]]), _clean(r[col["CGSNL Gene Name"]])
        classes = [c.strip() for c in r[col["Trait Class"]].split(",") if c.strip()]
        to = _terms(r[col["Trait Ontology"]])
        go = _terms(r[col["Gene Ontology"]]) if "Gene Ontology" in col else []
        syn = [s.strip(" []*") for s in re.split(r"[,;]", r[col["Gene symbol synonym(s)"]]) if s.strip(" []*")] if "Gene symbol synonym(s)" in col else []
        for rap in raps:
            e = genes.setdefault(rap, {"id": [], "symbol": "", "name": "", "synonyms": [], "explanation": "", "classes": [], "to": [], "go": [], "msu": []})
            e["id"].append(r[col["Trait Gene Id"]])
            e["symbol"] = e["symbol"] or symbol
            e["name"] = e["name"] or name
            e["explanation"] = (e["explanation"] or r[col["Explanation"]].strip())[:400]
            for key, items in (("synonyms", syn), ("classes", classes)):
                e[key] += [x for x in items if x not in e[key]]
            for key, items in (("to", to), ("go", go)):
                have = {t["id"] for t in e[key]}
                e[key] += [t for t in items if t["id"] not in have]
        for m in re.findall(r"LOC_Os\d{2}g\d{5}", r[col["MSU ID"]]):
            msu.setdefault(msu_key(m), raps[0])
            for rap in raps:
                if m not in genes[rap]["msu"]:
                    genes[rap]["msu"].append(m)
    return {"genes": genes, "msu": msu}


class Oryzabase:
    def __init__(self, ctx):
        self.ctx = ctx
        self.api = ApiClient(SOURCE, BASE, ctx, min_interval=0.2)
        self._data = None
        self._lock = threading.Lock()

    def ping(self):
        try:
            self.api.calls += 1
            r = self.api.session.head(BASE + "/rice/oryzabase/", timeout=12, allow_redirects=True)
            ok = r.status_code < 400
        except requests.RequestException:
            self.api.failures += 1
            return False
        if ok:
            self.api.ok_once = True
        return ok

    # ------------------------------------------------------------ the gene list
    def load(self):
        """Download (or read from the cache) the reduced gene list. Raises NetworkError when it cannot be had."""
        with self._lock:
            if self._data is not None:
                return self._data
            self.ctx.check()
            cache = self.ctx.cache
            if not self.ctx.refresh:
                hit, payload, _ = cache.get(CACHE_KEY, config.CACHE_DAYS)
                if hit and payload:
                    self.api.cached += 1
                    self._data = payload
                    return payload
            try:
                self.api.calls += 1
                r = self.api.session.get(LIST_URL, timeout=180, headers={"Accept": "text/plain,*/*"})
            except requests.RequestException as e:
                self.api.failures += 1
                self.api.last_error = type(e).__name__
                raise NetworkError(f"{SOURCE} could not be reached ({type(e).__name__}).")
            if r.status_code != 200 or len(r.content) < 100000:
                self.api.failures += 1
                self.api.last_error = f"HTTP {r.status_code}"
                raise ApiError(f"{SOURCE}: the gene list was not delivered (HTTP {r.status_code}).")
            try:
                raw = r.content.decode("utf-8")
            except UnicodeDecodeError:
                raw = r.content.decode("cp932", errors="replace")
            data = parse_list(raw)
            data["retrieved"] = today()
            self.api.ok_once = True
            cache.put(CACHE_KEY, SOURCE, LIST_URL, 200, data)
            self._data = data
            return data

    # ------------------------------------------------------------ lookups
    def entry(self, rap_id, msu_id=""):
        """Oryzabase record of a gene by RAP-DB ID (or, failing that, MSU ID), or None when Oryzabase does not list the gene."""
        d = self.load()
        key = re.sub(r"\.\d+$", "", rap_id or "")
        e = d["genes"].get(key)
        if e is None and msu_id:
            rap = d["msu"].get(msu_key(msu_id))
            e = d["genes"].get(rap) if rap else None
        if e is None:
            return None
        return {**e, "url": gene_url(e["id"][0]), "retrieved": d.get("retrieved", "")}

    def rap_from_msu(self, msu_id):
        """RAP-DB ID of an MSU (LOC_Os...) locus, for the genes Oryzabase lists. None otherwise."""
        return self.load()["msu"].get(msu_key(msu_id))

    @staticmethod
    def prov(entry):
        return {"source": SOURCE, "url": entry["url"], "retrieved": entry.get("retrieved") or today(),
                "type": "curated gene name, trait classes and Trait Ontology terms (keyed by RAP-DB ID)"}
