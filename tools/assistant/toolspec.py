"""Tool schemas the model may call. The proxy sends these to the API; the
desktop client executes them locally (tools/assistant/agent.py). Keep this
file free of heavy imports."""

TOOL_IDS = [
    "mendelian", "linkage", "qtl", "gene_extractor", "gene_prioritization", "pybsaseq", "qtl_seq", "diversity",
    "gwas", "et_gwas", "clustering", "additive_dominance", "genomic_prediction", "cross_selection",
    "heterosis", "phenotypic",
]

TOOLS = [
    {
        "name": "search_manual",
        "description": "Search the TerraTools user manual: input file formats, settings, how outputs are read, workflows, troubleshooting and glossary. Returns the best-matching passages with their chapter and section.",
        "input_schema": {
            "type": "object",
            "properties": {
                "query": {"type": "string", "description": "Keywords or a short question, e.g. 'QTL genotype file header rows' or 'PVE meaning GWAS'."},
                "chapter": {"type": "string", "description": "Optional chapter filter, e.g. 'QTL Mapping', 'GWAS', 'Troubleshooting', 'Glossary'."},
                "max_results": {"type": "integer", "minimum": 1, "maximum": 8},
            },
            "required": ["query"],
            "additionalProperties": False,
        },
    },
    {
        "name": "list_runs",
        "description": "List the analyses the user has finished in this session (tool, status, time) and the output files each produced, plus any files the user attached. Call this first when the user asks about their results.",
        "input_schema": {"type": "object", "properties": {}, "additionalProperties": False},
    },
    {
        "name": "get_run_result",
        "description": "Get the result summary of the user's latest run of one tool (the same numbers the tool showed on screen: counts, thresholds, tables, file names). Long lists are shortened; use read_run_file for full tables.",
        "input_schema": {
            "type": "object",
            "properties": {"tool": {"type": "string", "enum": TOOL_IDS}},
            "required": ["tool"],
            "additionalProperties": False,
        },
    },
    {
        "name": "read_run_file",
        "description": "Read one output file from the user's runs: CSV/TSV/Excel tables (columns, numeric summary, a window of rows), text/JSON, or a PNG/JPG plot, which you will be able to see. Use a file name returned by list_runs.",
        "input_schema": {
            "type": "object",
            "properties": {
                "name": {"type": "string", "description": "File name from list_runs, e.g. 'qtl_significant_CIM.csv'."},
                "tool": {"type": "string", "enum": TOOL_IDS, "description": "Optional: restrict to this tool's run."},
                "offset": {"type": "integer", "minimum": 0, "description": "First row to show (tables)."},
                "limit": {"type": "integer", "minimum": 1, "maximum": 60, "description": "Rows to show (tables), default 20."},
                "sheet": {"type": "string", "description": "Sheet name for Excel files."},
            },
            "required": ["name"],
            "additionalProperties": False,
        },
    },
    {
        "name": "search_literature",
        "description": "Live search of the scientific literature (Europe PMC, with Crossref as fallback). Returns real records with title, authors, journal, year, DOI/PMID, citation count and abstract. This is the ONLY allowed source for citations.",
        "input_schema": {
            "type": "object",
            "properties": {
                "query": {"type": "string", "description": "Keywords, e.g. 'QTL-seq SNP-index bulked segregant rice' or 'composite interval mapping Zeng'. Europe PMC syntax such as AUTH:\"Zeng\" or TITLE:\"QTL-seq\" is accepted."},
                "max_results": {"type": "integer", "minimum": 1, "maximum": 10},
                "year_from": {"type": "integer", "description": "Only papers published in or after this year."},
                "sort": {"type": "string", "enum": ["relevance", "cited"], "description": "'cited' surfaces foundational method papers."},
                "open_access_only": {"type": "boolean"},
            },
            "required": ["query"],
            "additionalProperties": False,
        },
    },
    {
        "name": "verify_doi",
        "description": "Check that a DOI exists and get its title, authors, journal and year. Use it before repeating any reference the user gives you.",
        "input_schema": {
            "type": "object",
            "properties": {"doi": {"type": "string"}},
            "required": ["doi"],
            "additionalProperties": False,
        },
    },
]

TOOL_NAMES = {t["name"] for t in TOOLS}
