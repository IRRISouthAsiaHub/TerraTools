import os
import sys
import json
import shutil
from pathlib import Path

from PySide6.QtCore import QObject, Qt, Signal, Slot, QThread, QUrl
from PySide6.QtWidgets import QFileDialog
from PySide6.QtWebChannel import QWebChannel
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QMainWindow

from tools.gene_extractor.backend import gene_extractor_runner
from ui.frontend import load_frontend
from tools.common.paths import results_dir


# ============================================================
# PATHS
# ============================================================

gene_EXTRACTOR_DIR = Path(__file__).resolve().parent

RESULT_DIR = results_dir("gene_extractor", gene_EXTRACTOR_DIR / "results")


# ============================================================
# BRIDGE
# ============================================================

class geneExtractorBridge(QObject):

    resultReady = Signal(str)
    statusUpdate = Signal(str)
    errorOccurred = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)

        self.last_result = None

    # --------------------------------------------------------
    # FILE BROWSING
    # --------------------------------------------------------

    @Slot(result=str)
    def browse_extracted_markers(self):

        path, _ = QFileDialog.getOpenFileName(
            None,
            "Select Extracted Markers Excel File",
            "",
            "Excel Files (*.xlsx *.xls);;All Files (*)"
        )

        return path or ""

    @Slot(result=str)
    def browse_gene_markers(self):

        path, _ = QFileDialog.getOpenFileName(
            None,
            "Select gene Markers Excel File",
            "",
            "Excel Files (*.xlsx *.xls);;All Files (*)"
        )

        return path or ""

    @Slot(result=str)
    def browse_rice_genes(self):

        path, _ = QFileDialog.getOpenFileName(
            None,
            "Select Rice Gene List Excel File",
            "",
            "Excel Files (*.xlsx *.xls);;All Files (*)"
        )

        return path or ""

    # --------------------------------------------------------
    # RUN EXTRACTION
    # --------------------------------------------------------

    @Slot(str)
    def run_extraction(self, config_json):

        print("[gene EXTRACTOR BRIDGE] run_extraction() called", flush=True)

        try:

            config = json.loads(config_json)

            print("\n========== EXTRACTION CONFIGURATION ==========")
            print(config)
            print("================================================\n")

            # Same validation the Shiny observeEvent(input$extract, ...) did
            if (not config.get("extracted_markers_path")
                    or not config.get("gene_markers_path")
                    or not config.get("rice_gene_path")):

                self._analysis_error(
                    "Please upload all three required files."
                )

                return

            self.statusUpdate.emit(
                "Extracting markers and genes..."
            )

            thread = QThread()

            worker = geneExtractorWorker(
                config,
                RESULT_DIR
            )

            worker.moveToThread(thread)

            thread.started.connect(
                worker.run, Qt.QueuedConnection
            )

            worker.finished.connect(
                self._analysis_finished, Qt.QueuedConnection
            )

            worker.error.connect(
                self._analysis_error, Qt.QueuedConnection
            )

            worker.finished.connect(
                thread.quit, Qt.QueuedConnection
            )

            worker.finished.connect(
                worker.deleteLater, Qt.QueuedConnection
            )

            thread.finished.connect(
                thread.deleteLater, Qt.QueuedConnection
            )

            self._thread = thread
            self._worker = worker

            print("[gene EXTRACTOR BRIDGE] about to call self._thread.start()", flush=True)
            thread.start()
            print("[gene EXTRACTOR BRIDGE] self._thread.start() returned", flush=True)

        except Exception as e:

            self._analysis_error(
                str(e)
            )

    @Slot(dict)
    def _analysis_finished(self, result):

        print("[gene EXTRACTOR BRIDGE] _analysis_finished() called, status:", result.get("status"), flush=True)

        self.last_result = result

        if result.get("status") == "error":
            self.statusUpdate.emit(
                "Extraction failed."
            )
        else:
            self.statusUpdate.emit(
                "Data extraction completed successfully! Navigate to the tabs above to view and download the results."
            )

        self.resultReady.emit(
            json.dumps(result)
        )

    @Slot(str)
    def _analysis_error(self, message):

        print("[gene EXTRACTOR BRIDGE] _analysis_error() called:", message, flush=True)

        result = {
            "status": "error",
            "message": message
        }

        self.last_result = result

        self.statusUpdate.emit(
            "Extraction failed."
        )

        self.resultReady.emit(
            json.dumps(result)
        )

    # --------------------------------------------------------
    # DOWNLOADS
    # --------------------------------------------------------

    @Slot(str)
    def download_result(self, result_type):

        print("[gene EXTRACTOR BRIDGE] download_result() called for:", result_type, flush=True)

        if not self.last_result or self.last_result.get("status") != "success":

            self.statusUpdate.emit(
                "Run extraction first."
            )

            return

        mapping = {

            "markers":
                self.last_result.get(
                    "markers_file"
                ),

            "genes":
                self.last_result.get(
                    "genes_file"
                )
        }

        default_names = {
            "markers": "Markers_with_physical_position.xlsx",
            "genes": "Genes_in_gene_regions.xlsx"
        }

        source = mapping.get(
            result_type
        )

        if not source or not os.path.exists(source):

            self.statusUpdate.emit(
                "Requested result is not available."
            )

            return

        destination, _ = QFileDialog.getSaveFileName(
            None,
            "Save Result",
            default_names.get(result_type, os.path.basename(source)),
            "Excel Files (*.xlsx);;All Files (*)"
        )

        if not destination:
            return

        try:

            shutil.copy2(
                source,
                destination
            )

            print("[gene EXTRACTOR BRIDGE] saved", result_type, "to", destination, flush=True)

            self.statusUpdate.emit(
                "Result saved."
            )

        except Exception as e:

            import traceback
            traceback.print_exc()

            self.statusUpdate.emit(
                f"Download failed: {e}"
            )


# ============================================================
# EXTRACTION WORKER
# ============================================================

class geneExtractorWorker(QObject):

    finished = Signal(dict)
    error = Signal(str)

    def __init__(
        self,
        config,
        output_dir
    ):

        super().__init__()

        self.config = config
        self.output_dir = output_dir

    @Slot()
    def run(self):

        print("[gene EXTRACTOR WORKER] run() called on worker thread", flush=True)

        try:

            print("[gene EXTRACTOR WORKER] calling gene_extractor_runner.run_extraction()...", flush=True)

            result = gene_extractor_runner.run_extraction(
                self.config,
                results_dir=self.output_dir
            )

            print("[gene EXTRACTOR WORKER] gene_extractor_runner.run_extraction() returned, status:", result.get("status"), flush=True)

            self.finished.emit(
                result
            )
            print("[gene EXTRACTOR WORKER] finished signal emitted", flush=True)

        except Exception as e:

            print("[gene EXTRACTOR WORKER] error:", e, flush=True)

            self.error.emit(
                str(e)
            )


# ============================================================
# WINDOW
# ============================================================

class geneExtractorWindow(QMainWindow):

    def __init__(self, parent=None):

        super().__init__(
            parent
        )

        self.setWindowTitle(
            "TerraTools — gene Marker and Gene Extractor"
        )

        self.showMaximized()

        self.browser = QWebEngineView()

        self.setCentralWidget(
            self.browser
        )

        # ----------------------------------------------------
        # WebChannel
        # ----------------------------------------------------

        self.channel = QWebChannel(
            self.browser.page()
        )

        self.bridge = geneExtractorBridge(
            self
        )

        self.channel.registerObject(
            "geneExtractorBridge",
            self.bridge
        )

        self.browser.page().setWebChannel(
            self.channel
        )

        # ----------------------------------------------------
        # Load HTML
        # ----------------------------------------------------

        load_frontend(self.browser, "gene_extractor")


# ============================================================
# STANDALONE TEST
# ============================================================

if __name__ == "__main__":

    from PySide6.QtWidgets import QApplication

    app = QApplication(
        sys.argv
    )

    window = geneExtractorWindow()

    window.show()

    sys.exit(
        app.exec()
    )
