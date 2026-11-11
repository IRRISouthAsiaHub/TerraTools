"""
Entry point of the Gene Prioritization tool (the shape every TerraTools runner has: a config dict in, a result dict out).

config: {
    "species": "oryza_sativa" | "triticum_aestivum" | "zea_mays",
    "trait": "Drought tolerance",
    "gene_text": "pasted gene IDs",              optional   (the gene_file or the gene_text is required)
    "gene_file": "path.(txt|csv|tsv|xlsx)",       optional   for example the Gene_ID_list.txt of Genes Extraction
    "refresh_cache": false                        optional   ignore cached API responses and annotations
}
Returns {"status": "success", ...} (see engine.Engine._assemble), {"status": "cancelled"} or {"status": "error", "message": ...}.
"""

from __future__ import annotations

import traceback

from tools.gene_prioritization import config
from tools.gene_prioritization.backend import export
from tools.gene_prioritization.backend.engine import Engine, PrioritizationError
from tools.gene_prioritization.backend.http import Cancelled, RunContext


def run_prioritization(cfg, results_dir=None, progress=None, ctx=None):
    engine = None
    try:
        engine = Engine(cfg, progress=progress, ctx=ctx)
        result = engine.run()
        result["files"] = export.write_all(result, results_dir or config.RESULT_DIR)
        return result
    except PrioritizationError as e:
        return {"status": "error", "message": str(e), "api_status": engine.api_status() if engine else []}
    except Cancelled:
        return {"status": "cancelled", "message": "The run was cancelled."}
    except Exception as e:                                              # never take the app down
        traceback.print_exc()
        return {"status": "error", "message": f"Gene Prioritization failed: {e}"}
    finally:
        if engine is not None:
            engine.cache.close()


def check_services():
    """Ping the four services (the 'API connection status' of the page)."""
    import tempfile
    from tools.gene_prioritization.backend.cache import Cache
    from tools.gene_prioritization.backend.ensembl import Ensembl
    from tools.gene_prioritization.backend.gramene import Gramene
    from tools.gene_prioritization.backend.http import ApiError, NetworkError
    from tools.gene_prioritization.backend.ncbi import Ncbi
    from tools.gene_prioritization.backend.plant_reactome import PlantReactome

    cache = Cache(tempfile.gettempdir() + "/tt_gp_ping.sqlite")
    ctx = RunContext(cache)
    out = []
    for name, cls in (("Ensembl Plants", Ensembl), ("Gramene", Gramene), ("Plant Reactome", PlantReactome), ("NCBI", Ncbi)):
        try:
            ok = bool(cls(ctx).ping())
            out.append({"source": name, "state": "ok" if ok else "unavailable"})
        except (NetworkError, ApiError) as e:
            out.append({"source": name, "state": "unavailable", "error": str(e)})
    cache.close()
    return out
