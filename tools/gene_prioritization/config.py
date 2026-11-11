"""
Gene Prioritization — configuration.

SPECIES is the one place that says which crops are supported. A species is ENABLED only when its identifiers and the
API calls the tool needs were checked against the live services (see "verified" below). To add a crop, add an entry
here: nothing else in the code names a species.

Verified on 2026-10-30 (live requests):
  Ensembl Plants REST  https://rest.ensembl.org      oryza_sativa, triticum_aestivum, zea_mays — lookup, overlap, homology, xrefs, phenotype
  Gramene REST         https://data.gramene.org/v69  same three — genes (annotation, GO, domains, QTL trait ontology, orthologs)
  Plant Reactome       https://plantreactome.gramene.org/ContentService  — search, pathways of a gene product
  NCBI E-utilities     https://eutils.ncbi.nlm.nih.gov/entrez/eutils     — gene, pubmed
  Oryzabase (rice)     https://shigen.nig.ac.jp/rice/oryzabase           — curated gene list (RAP-DB ID, MSU ID, trait classes, Trait Ontology)
  Cicer arietinum is in neither Ensembl Plants nor Gramene, so it is listed but not enabled.
"""

from pathlib import Path

from tools.common.paths import results_dir

TOOL_DIR = Path(__file__).resolve().parent
RESULT_DIR = results_dir("gene_prioritization", TOOL_DIR / "ui" / "results")
CACHE_DB = RESULT_DIR / "api_cache.sqlite"

# ------------------------------------------------------------------ species
SPECIES = {
    "oryza_sativa": {
        "label": "Oryza sativa (Rice)",
        "common": "rice",
        "scientific": "Oryza sativa",
        "ensembl": "oryza_sativa",          # Ensembl Plants species name
        "gramene": "oryza_sativa",          # Gramene system_name
        "reactome": "Oryza sativa",         # Plant Reactome species name
        "reactome_taxon": "4530",
        "ncbi_taxon": "4530",
        "assembly": "IRGSP-1.0 (RAP-DB gene models)",
        "id_help": "RAP-DB / Ensembl IDs such as Os01g0100100, gene symbols, UniProt accessions or NCBI Gene IDs. "
                   "MSU (LOC_Os…) IDs are translated only for the genes Oryzabase lists.",
        "id_example": "Os02g0766700",
        "orthologs_in": ["arabidopsis_thaliana", "zea_mays"],
        "oryzabase": True,                  # Oryzabase (curated rice genes keyed by RAP-DB ID) is an extra evidence source
        "enabled": True,
    },
    "triticum_aestivum": {
        "label": "Triticum aestivum (Wheat)",
        "common": "wheat",
        "scientific": "Triticum aestivum",
        "ensembl": "triticum_aestivum",
        "gramene": "triticum_aestivum",
        "reactome": "Triticum aestivum",
        "reactome_taxon": "4565",
        "ncbi_taxon": "4565",
        "assembly": "IWGSC RefSeq v1.0 (Ensembl Plants 'triticum_aestivum')",
        "id_help": "Ensembl Plants IDs such as TraesCS4B02G300000, gene symbols, UniProt accessions or NCBI Gene IDs.",
        "id_example": "TraesCS4B02G300000",
        "orthologs_in": ["arabidopsis_thaliana", "oryza_sativa", "zea_mays"],
        "enabled": True,
    },
    "zea_mays": {
        "label": "Zea mays (Maize)",
        "common": "maize",
        "scientific": "Zea mays",
        "ensembl": "zea_mays",
        "gramene": "zea_mays",
        "reactome": "Zea mays",
        "reactome_taxon": "4577",
        "ncbi_taxon": "4577",
        "assembly": "Zm-B73-REFERENCE-NAM-5.0",
        "id_help": "B73 NAM-5.0 IDs such as Zm00001eb000010, gene symbols, UniProt accessions or NCBI Gene IDs.",
        "id_example": "Zm00001eb000010",
        "orthologs_in": ["arabidopsis_thaliana", "oryza_sativa"],
        "enabled": True,
    },
    "cicer_arietinum": {
        "label": "Cicer arietinum (Chickpea)",
        "common": "chickpea",
        "scientific": "Cicer arietinum",
        "ensembl": None,
        "gramene": None,
        "reactome": None,
        "reactome_taxon": None,
        "ncbi_taxon": "3827",
        "assembly": "",
        "id_help": "",
        "id_example": "",
        "orthologs_in": [],
        "enabled": False,
        "reason": "Chickpea is not in Ensembl Plants or Gramene, so gene coordinates, regions, orthologs and "
                  "annotations cannot be retrieved. It will be enabled when a verified source exists.",
    },
}

# ------------------------------------------------------------------ evidence
# category key -> (label, weight)
EVIDENCE_WEIGHTS = {
    "association": ("Known gene-trait association", 0.30),
    "function": ("Functional relevance to the trait", 0.25),
    "pathway": ("Pathway relevance", 0.20),
    "ortholog": ("Ortholog evidence", 0.15),
    "literature": ("Literature support", 0.10),
}

# ------------------------------------------------------------------ run limits
LITERATURE_TOP_N = 80            # genes that get a PubMed search (ranked after the other evidence)
PUBS_PER_GENE = 15
WORKERS = 8                       # parallel annotation requests; each service is still held to its own request rate

# ------------------------------------------------------------------ http
USER_AGENT = "TerraTools-GenePrioritization/1.0 (IRRI-SAH plant breeding toolkit; +https://www.irri.org)"
TIMEOUT = 30                     # seconds
RETRIES = 4
BACKOFF = 1.5                    # seconds, doubled after each failed attempt
BREAKER_FAILURES = 8             # consecutive failed attempts after which a service is skipped for the rest of the run
CACHE_DAYS = 30
NEGATIVE_CACHE_DAYS = 7

ENSEMBL_URL = "https://rest.ensembl.org"
GRAMENE_URL = "https://data.gramene.org/v69"
REACTOME_URL = "https://plantreactome.gramene.org/ContentService"
NCBI_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# NCBI asks for at most 3 requests/s without an API key and 10/s with one (environment variable NCBI_API_KEY).
NCBI_API_KEY_ENV = "NCBI_API_KEY"
NCBI_EMAIL_ENV = "NCBI_EMAIL"


def species_options():
    """What the page needs for its drop-down."""
    return [
        {"value": k, "label": v["label"], "enabled": bool(v["enabled"]), "reason": v.get("reason", ""),
         "id_help": v["id_help"], "id_example": v["id_example"], "assembly": v["assembly"]}
        for k, v in SPECIES.items()
    ]
