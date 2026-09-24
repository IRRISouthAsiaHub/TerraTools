"""Every tool's Python<->JS bridge, hosted inside the single main window.

The imports below are deliberately written out as plain `from ... import`
statements (not importlib strings): PyInstaller discovers which libraries
to bundle by scanning import statements, so a string-based import would
leave the tools' dependencies (pandas, statsmodels, ...) out of the exe.
"""

import traceback


def _mendelian():
    from tools.mendelian.ui.mendelian_window import MendelianBridge
    return MendelianBridge


def _linkage():
    from tools.linkage.ui.linkage_window import LinkageBridge
    return LinkageBridge


def _qtl():
    from tools.qtl.ui.qtl_window import QTLBridge
    return QTLBridge


def _gene_extractor():
    from tools.gene_extractor.ui.gene_extractor_window import geneExtractorBridge
    return geneExtractorBridge


def _pybsaseq():
    from tools.pybsaseq.ui.pybsaseq_window import PyBSASeqBridge
    return PyBSASeqBridge


def _qtl_seq():
    from tools.qtl_seq.ui.qtl_seq_window import QTLSeqBridge
    return QTLSeqBridge


def _diversity():
    from tools.diversity.ui.diversity_window import DiversityBridge
    return DiversityBridge


def _gwas():
    from tools.gwas.ui.gwas_window import GWASBridge
    return GWASBridge


def _clustering():
    from tools.clustering.ui.clustering_window import ClusteringBridge
    return ClusteringBridge


def _additive_dominance():
    from tools.additive_dominance.ui.additive_dominance_window import AdditiveDominanceBridge
    return AdditiveDominanceBridge


def _genomic_prediction():
    from tools.genomic_prediction.ui.genomic_prediction_window import GenomicPredictionBridge
    return GenomicPredictionBridge


def _cross_selection():
    from tools.cross_selection.ui.cross_window import CrossBridge
    return CrossBridge


def _heterosis():
    from tools.heterosis.ui.heterosis_window import HeterosisBridge
    return HeterosisBridge


def _phenotypic():
    from tools.phenotypic.ui.phenotypic_window import PhenotypicBridge
    return PhenotypicBridge


# (loader, QWebChannel object name). The object name is what the React page
# asks for (frontend/src/pages/*Tool.jsx).
TOOL_BRIDGES = [
    (_mendelian, "mendelianBridge"),
    (_linkage, "linkageBridge"),
    (_qtl, "qtlBridge"),
    (_gene_extractor, "geneExtractorBridge"),
    (_pybsaseq, "pybsaseqBridge"),
    (_qtl_seq, "qtlSeqBridge"),
    (_diversity, "diversityBridge"),
    (_gwas, "gwasBridge"),
    (_clustering, "clusteringBridge"),
    (_additive_dominance, "additiveDominanceBridge"),
    (_genomic_prediction, "genomicPredictionBridge"),
    (_cross_selection, "crossBridge"),
    (_heterosis, "heterosisBridge"),
    (_phenotypic, "phenotypicBridge"),
]


def register_tool_bridges(channel, window):
    """Create every tool bridge with `window` as its parent and register it
    on `channel`. The bridges use the window for file dialogs, and the
    linkage/QTL bridges push some results via window.browser.

    A tool that fails to import is logged and skipped, so one broken tool
    never takes down the whole app. Returns {object name: bridge}."""
    bridges = {}
    for load, name in TOOL_BRIDGES:
        try:
            bridge = load()(window)
            channel.registerObject(name, bridge)
            bridges[name] = bridge
        except Exception:
            print(f"[TerraTools] could not load the {name} tool:", flush=True)
            traceback.print_exc()
    return bridges
