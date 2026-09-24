from PySide6.QtWidgets import QMessageBox


class ToolManager:

    def __init__(self, parent=None):

        self.parent = parent

        # Keep references to all open tool windows.
        self.open_windows = {}

    # ========================================================
    # TOOL ROUTER
    # ========================================================

    def open_tool(self, tool_name):

        tool_name = str(
            tool_name
        ).strip().lower()

        print(
            f"[ToolManager] Opening tool: "
            f"{tool_name!r}"
        )

        if tool_name == "gwas":

            self.open_gwas()

        elif tool_name == "qtl":

            self.open_qtl()

        elif tool_name == "mendelian":

            self.open_mendelian()

        elif tool_name == "diversity":

            self.open_diversity()
        
        elif tool_name == "qtl_seq":
            self.open_qtl_seq()
        
        elif tool_name == "gene_extractor":

            self.open_gene_extractor()

        elif tool_name == "linkage":

            self.open_linkage()
            
        elif tool_name == "phenotypic":
 
            self.open_phenotypic()
            
        elif tool_name == "pybsaseq":
 
            self.open_pybsaseq()

        elif tool_name == "clustering":
 
            self.open_clustering()
            
        elif tool_name == "genomic_prediction":
 
            self.open_genomic_prediction()

        elif tool_name == "additive_dominance":

            self.open_additive_dominance()

        elif tool_name == "cross_selection":
            
            self.open_cross_selection()

        elif tool_name == "heterosis":

            self.open_heterosis()
 

        else:

            self.show_message(
                "TerraTools",
                f"Unknown tool: {tool_name}"
            )

    # ========================================================
    # GWAS
    # ========================================================

    def open_gwas(self):

        print(
            "[ToolManager] Loading GWAS..."
        )

        try:

            from tools.gwas.ui.gwas_window import GWASWindow

            if "gwas" in self.open_windows:

                window = self.open_windows["gwas"]

                window.show()
                window.raise_()
                window.activateWindow()

                return

            window = GWASWindow(
                self.parent
            )

            self.open_windows[
                "gwas"
            ] = window

            window.show()

            window.raise_()
            window.activateWindow()

        except Exception:

            import traceback

            traceback.print_exc()

            raise

    # ========================================================
    # QTL
    # ========================================================

    def open_qtl(self):

        print(
            "[ToolManager] Loading QTL..."
        )

        try:

            from tools.qtl.ui.qtl_window import QTLWindow

            print(
                "[ToolManager] QTLWindow imported"
            )

            if "qtl" in self.open_windows:

                print(
                    "[ToolManager] "
                    "Reusing existing QTL window"
                )

                window = self.open_windows[
                    "qtl"
                ]

                window.show()
                window.raise_()
                window.activateWindow()

                return

            print(
                "[ToolManager] "
                "Creating QTLWindow"
            )

            window = QTLWindow(
                self.parent
            )

            # IMPORTANT:
            # Keep a Python reference so Qt does not
            # destroy the window when this function returns.
            self.open_windows[
                "qtl"
            ] = window

            window.show()

            window.raise_()
            window.activateWindow()

            print(
                "[ToolManager] QTL window shown"
            )

        except Exception:

            import traceback

            print(
                "\n========== QTL ERROR =========="
            )

            traceback.print_exc()

            print(
                "===============================\n"
            )

            raise

    # ========================================================
    # MENDELIAN
    # ========================================================

    def open_mendelian(self):

        print(
            "[ToolManager] Loading Mendelian..."
        )

        try:

            from tools.mendelian.ui.mendelian_window import MendelianWindow

            if "mendelian" in self.open_windows:

                window = self.open_windows["mendelian"]

                window.show()
                window.raise_()
                window.activateWindow()

                return

            window = MendelianWindow(
                self.parent
            )

            self.open_windows[
                "mendelian"
            ] = window

            window.show()

            window.raise_()
            window.activateWindow()

        except Exception:

            import traceback

            traceback.print_exc()

            raise

    # ========================================================
    # DIVERSITY
    # ========================================================

    def open_diversity(self):

        print(
            "[ToolManager] Loading Diversity..."
        )

        try:

            from tools.diversity.ui.diversity_window import DiversityWindow

            if "diversity" in self.open_windows:

                window = self.open_windows["diversity"]

                window.show()
                window.raise_()
                window.activateWindow()

                return

            window = DiversityWindow(
                self.parent
            )

            self.open_windows[
                "diversity"
            ] = window

            window.show()

            window.raise_()
            window.activateWindow()

        except Exception:

            import traceback

            traceback.print_exc()

            raise
       # ========================================================
    # QTL-SEQ (BSA)
    # ========================================================
 
    def open_qtl_seq(self):
 
        print(
            "[ToolManager] Loading QTL-seq..."
        )
 
        try:
 
            from tools.qtl_seq.ui.qtl_seq_window import (
                QTLSeqWindow
            )
 
            print(
                "[ToolManager] QTLSeqWindow imported"
            )
 
            if "qtl_seq" in self.open_windows:
 
                print(
                    "[ToolManager] "
                    "Reusing existing QTL-seq window"
                )
 
                window = self.open_windows[
                    "qtl_seq"
                ]
 
                window.show()
                window.raise_()
                window.activateWindow()
 
                return
 
            print(
                "[ToolManager] "
                "Creating QTLSeqWindow"
            )
 
            window = QTLSeqWindow(
                self.parent
            )
 
            self.open_windows[
                "qtl_seq"
            ] = window
 
            window.show()
 
            window.raise_()
            window.activateWindow()
 
            print(
                "[ToolManager] QTL-seq window shown"
            )
 
        except Exception:
 
            import traceback
 
            print(
                "\n========== QTL-SEQ ERROR =========="
            )
 
            traceback.print_exc()
 
            print(
                "====================================\n"
            )
 
            raise
            
        # ========================================================
# gene EXTRACTOR
# ========================================================

    def open_gene_extractor(self):

        print(
            "[ToolManager] Loading Gene Extractor..."
        )

        try:

            from tools.gene_extractor.ui.gene_extractor_window import (
                geneExtractorWindow
            )

            print(
                "[ToolManager] geneExtractorWindow imported"
            )

            if "gene_extractor" in self.open_windows:

                print(
                    "[ToolManager] "
                    "Reusing existing gene Extractor window"
                )

                window = self.open_windows[
                    "gene_extractor"
                ]

                window.show()
                window.raise_()
                window.activateWindow()

                return

            print(
                "[ToolManager] "
                "Creating geneExtractorWindow"
            )

            window = geneExtractorWindow(
                self.parent
            )

            self.open_windows[
                "gene_extractor"
            ] = window

            window.show()

            window.raise_()
            window.activateWindow()

            print(
                "[ToolManager] gene Extractor window shown"
            )

        except Exception:

            import traceback

            print(
                "\n========== gene EXTRACTOR ERROR =========="
            )

            traceback.print_exc()

            print(
                "===========================================\n"
            )

            raise

    # ========================================================
# LINKAGE
# ========================================================

    def open_linkage(self):

        print(
            "[ToolManager] Loading Linkage..."
        )

        try:

            from tools.linkage.ui.linkage_window import LinkageWindow

            print(
                "[ToolManager] LinkageWindow imported"
            )

            if "linkage" in self.open_windows:

                print(
                    "[ToolManager] "
                    "Reusing existing Linkage window"
                )

                window = self.open_windows[
                    "linkage"
                ]

                window.show()
                window.raise_()
                window.activateWindow()

                return

            print(
                "[ToolManager] "
                "Creating LinkageWindow"
            )

            window = LinkageWindow(
                self.parent
            )

            # Keep Python reference so Qt does not
            # destroy the window after this function returns.
            self.open_windows[
                "linkage"
            ] = window

            window.show()

            window.raise_()
            window.activateWindow()

            print(
                "[ToolManager] Linkage window shown"
            )

        except Exception:

            import traceback

            print(
                "\n========== LINKAGE ERROR =========="
            )

            traceback.print_exc()

            print(
                "===================================\n"
            )

            raise
        # ========================================================
    # PHENOTYPIC DATA ANALYSIS
    # ========================================================
 
    def open_phenotypic(self):
 
        print(
            "[ToolManager] Loading Phenotypic Data Analysis..."
        )
 
        try:
 
            from tools.phenotypic.ui.phenotypic_window import PhenotypicWindow
 
            if "phenotypic" in self.open_windows:
 
                window = self.open_windows["phenotypic"]
 
                window.show()
                window.raise_()
                window.activateWindow()
 
                return
 
            window = PhenotypicWindow(
                self.parent
            )
 
            self.open_windows[
                "phenotypic"
            ] = window
 
            window.show()
 
            window.raise_()
            window.activateWindow()
 
        except Exception:
 
            import traceback
 
            traceback.print_exc()
 
            raise
            
        # ========================================================
    # PYBSASEQ
    # ========================================================
 
    def open_pybsaseq(self):
 
        print(
            "[ToolManager] Loading PyBSA-Seq..."
        )
 
        try:
 
            from tools.pybsaseq.ui.pybsaseq_window import PyBSASeqWindow
 
            if "pybsaseq" in self.open_windows:
 
                window = self.open_windows["pybsaseq"]
 
                window.show()
                window.raise_()
                window.activateWindow()
 
                return
 
            window = PyBSASeqWindow(
                self.parent
            )
 
            self.open_windows[
                "pybsaseq"
            ] = window
 
            window.show()
 
            window.raise_()
            window.activateWindow()
 
        except Exception:
 
            import traceback
 
            traceback.print_exc()
 
            raise

    # ========================================================
    # CLUSTERING / GENOMIC SELECTION
    # ========================================================
 
    def open_clustering(self):
 
        print(
            "[ToolManager] Loading Clustering..."
        )
 
        try:
 
            from tools.clustering.ui.clustering_window import ClusteringWindow
 
            if "clustering" in self.open_windows:
 
                window = self.open_windows["clustering"]
 
                window.show()
                window.raise_()
                window.activateWindow()
 
                return
 
            window = ClusteringWindow(
                self.parent
            )
 
            self.open_windows[
                "clustering"
            ] = window
 
            window.show()
 
            window.raise_()
            window.activateWindow()
 
        except Exception:
 
            import traceback
 
            traceback.print_exc()
 
            raise
 
        # ========================================================
    # GENOMIC PREDICTION (BREEDSIGHT)
    # ========================================================
 
    def open_genomic_prediction(self):
 
        print(
            "[ToolManager] Loading Genomic Prediction..."
        )
 
        try:
 
            from tools.genomic_prediction.ui.genomic_prediction_window import GenomicPredictionWindow
 
            if "genomic_prediction" in self.open_windows:
 
                window = self.open_windows["genomic_prediction"]
 
                window.show()
                window.raise_()
                window.activateWindow()
 
                return
 
            window = GenomicPredictionWindow(
                self.parent
            )
 
            self.open_windows[
                "genomic_prediction"
            ] = window
 
            window.show()
 
            window.raise_()
            window.activateWindow()
 
        except Exception:
 
            import traceback
 
            traceback.print_exc()
 
            raise

       # ========================================================
    # ADDITIVE / DOMINANCE
    # ========================================================
 
    def open_additive_dominance(self):
 
        print(
            "[ToolManager] Loading Additive/Dominance..."
        )
 
        try:
 
            from tools.additive_dominance.ui.additive_dominance_window import (
                AdditiveDominanceWindow
            )
 
            print(
                "[ToolManager] AdditiveDominanceWindow imported"
            )
 
            if "additive_dominance" in self.open_windows:
 
                print(
                    "[ToolManager] "
                    "Reusing existing Additive/Dominance window"
                )
 
                window = self.open_windows[
                    "additive_dominance"
                ]
 
                window.show()
                window.raise_()
                window.activateWindow()
 
                return
 
            print(
                "[ToolManager] "
                "Creating AdditiveDominanceWindow"
            )
 
            window = AdditiveDominanceWindow(
                self.parent
            )
 
            self.open_windows[
                "additive_dominance"
            ] = window
 
            window.show()
 
            window.raise_()
            window.activateWindow()
 
            print(
                "[ToolManager] Additive/Dominance window shown"
            )
 
        except Exception:
 
            import traceback
 
            print(
                "\n========== ADDITIVE/DOMINANCE ERROR =========="
            )
 
            traceback.print_exc()
 
            print(
                "================================================\n"
            )
 
            raise
    # ========================================================
# cross_selection
# ========================================================

    def open_cross_selection(self):

        print(
            "[ToolManager] Loading cross_selection..."
        )

        try:

            from tools.cross_selection.ui.cross_window import CrossWindow

            print(
                "[ToolManager] CrossWindow imported"
            )

            if "cross_selection" in self.open_windows:

                print(
                    "[ToolManager] "
                    "Reusing existing cross_selection window"
                )

                window = self.open_windows["cross_selection"]

                window.show()
                window.raise_()
                window.activateWindow()

                return

            print(
                "[ToolManager] "
                "Creating CrossWindow"
            )

            window = CrossWindow(
                self.parent
            )

            # Keep a Python reference so Qt does not
            # destroy the window after this method returns.
            self.open_windows["cross_selection"] = window

            window.show()
            window.raise_()
            window.activateWindow()

            print(
                "[ToolManager] cross_selection window shown"
            )

        except Exception:

            import traceback

            print(
                "\n========== cross_selection ERROR =========="
            )

            traceback.print_exc()

            print(
                "===============================\n"
            )

            raise

    # ========================================================
    # HETEROSIS
    # ========================================================

    def open_heterosis(self):

        print(
            "[ToolManager] Loading Heterosis..."
        )

        try:

            from tools.heterosis.ui.heterosis_window import HeterosisWindow

            print(
                "[ToolManager] HeterosisWindow imported"
            )

            if "heterosis" in self.open_windows:

                print(
                    "[ToolManager] "
                    "Reusing existing Heterosis window"
                )

                window = self.open_windows[
                    "heterosis"
                ]

                window.show()
                window.raise_()
                window.activateWindow()

                return

            print(
                "[ToolManager] "
                "Creating HeterosisWindow"
            )

            window = HeterosisWindow(
                self.parent
            )

            # Keep Python reference so Qt does not
            # destroy the window after this function returns.
            self.open_windows[
                "heterosis"
            ] = window

            window.show()

            window.raise_()
            window.activateWindow()

            print(
                "[ToolManager] Heterosis window shown"
            )

        except Exception:

            import traceback

            print(
                "\n========== HETEROSIS ERROR =========="
            )

            traceback.print_exc()

            print(
                "=======================================\n"
            )

            raise
            # ========================================================
            # MESSAGE
            # ========================================================

    def show_message(
        self,
        title,
        message
    ):

        QMessageBox.information(
            self.parent,
            title,
            message
        )