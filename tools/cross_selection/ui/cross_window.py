import os
import sys
import json
from pathlib import Path

from PySide6.QtCore import QObject, Qt, Signal, Slot, QThread, QUrl
from tools.common.single_run import single_run
from PySide6.QtWidgets import QFileDialog
from PySide6.QtWebChannel import QWebChannel
from PySide6.QtWebEngineCore import QWebEngineSettings
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QMainWindow

from tools.cross_selection.backend import cross_runner
from ui.frontend import load_frontend
from tools.common.paths import results_dir


# ============================================================
# PATHS
# ============================================================

CROSS_DIR = Path(__file__).resolve().parent

RESULT_DIR = results_dir("cross_selection", CROSS_DIR / "results")


# ============================================================
# BRIDGE
# ============================================================

class CrossBridge(QObject):

    resultReady = Signal(str)
    statusUpdate = Signal(str)
    errorOccurred = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)

        self.last_result = None
        self.last_config = None

    # --------------------------------------------------------
    # FILE BROWSING
    # --------------------------------------------------------

    @Slot(result=str)
    def browse_geno(self):
        return self._browse("Select Allele Dosage Matrix", "CSV Files (*.csv *.txt *.tsv);;All Files (*)")

    @Slot(result=str)
    def browse_effects(self):
        return self._browse("Select Marker Effects", "CSV Files (*.csv *.txt *.tsv);;All Files (*)")

    @Slot(result=str)
    def browse_map(self):
        return self._browse("Select Genetic Map", "CSV Files (*.csv *.txt *.tsv);;All Files (*)")

    @staticmethod
    def _browse(title, flt):
        path, _ = QFileDialog.getOpenFileName(None, title, "", flt)
        return path or ""

    @Slot(result=str)
    def browse_vcf(self):
        return self._browse("Select VCF", "VCF Files (*.vcf *.vcf.gz);;All Files (*)")

    @Slot(str, str, str, str, str, result=str)
    def preview_markers(self, geno_path, effects_path, map_path, map_kind, cm_per_mb):
        try:
            return json.dumps(cross_runner.preview_markers(geno_path or None, effects_path or None, map_path or None,
                                                           map_kind or None, cm_per_mb or None))
        except Exception as e:
            return json.dumps({"status": "error", "message": str(e)})

    # --------------------------------------------------------
    # RUN CROSS SELECTION
    # --------------------------------------------------------

    @Slot(str)
    @single_run
    def run_cross_selection(self, config_json):

        print("[CROSS BRIDGE] run_cross_selection() called", flush=True)

        try:

            config = json.loads(config_json)

            self.last_config = config

            print("\n========== CROSS SELECTION CONFIGURATION ==========")
            print(config)
            print("=====================================================\n")

            self.statusUpdate.emit(
                "Preparing cross selection analysis..."
            )

            thread = QThread()

            worker = CrossWorker(
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

            print("[CROSS BRIDGE] about to call self._thread.start()", flush=True)
            thread.start()
            print("[CROSS BRIDGE] self._thread.start() returned", flush=True)

        except Exception as e:

            self._analysis_error(
                str(e)
            )

    @Slot(dict)
    def _analysis_finished(self, result):

        print("[CROSS BRIDGE] _analysis_finished() called, status:", result.get("status"), flush=True)

        self.last_result = result

        self.statusUpdate.emit(
            "Cross selection analysis completed successfully."
        )

        self.resultReady.emit(
            json.dumps(result)
        )

    @Slot(str)
    def _analysis_error(self, message):

        print("[CROSS BRIDGE] _analysis_error() called:", message, flush=True)

        result = {
            "status": "error",
            "message": message
        }

        self.last_result = result

        self.statusUpdate.emit(
            "Cross selection analysis failed."
        )

        self.resultReady.emit(
            json.dumps(result)
        )

    # --------------------------------------------------------
    # DOWNLOADS
    # --------------------------------------------------------

    @Slot(str)
    def download_result(self, result_type):

        print("[CROSS BRIDGE] download_result() called for:", result_type, flush=True)

        if not self.last_result or self.last_result.get("status") != "success":

            self.statusUpdate.emit(
                "Run an analysis first."
            )

            return

        ocs = self.last_result.get("ocs", {})

        mapping = {
            "ocs_crosses_csv": ocs.get("crosses_csv"),
            "ocs_contributions_csv": ocs.get("contributions_csv"),
            "ocs_frontier_csv": ocs.get("frontier_csv"),
            "ocs_frontier_plot": ocs.get("frontier_plot"),
            "ocs_usefulness_csv": ocs.get("usefulness_csv"),
            "ocs_kinship_csv": ocs.get("kinship_csv")
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
            os.path.basename(source),
            "All Files (*)"
        )

        if not destination:
            return

        try:

            import shutil

            shutil.copy2(
                source,
                destination
            )

            print("[CROSS BRIDGE] saved", result_type, "to", destination, flush=True)

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
# WORKER
# ============================================================

class CrossWorker(QObject):

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

        print("[CROSS WORKER] run() called on worker thread", flush=True)

        try:

            print("[CROSS WORKER] calling cross_runner.run_cross_selection()...", flush=True)

            result = cross_runner.run_cross_selection(
                self.config,
                results_dir=self.output_dir
            )

            print("[CROSS WORKER] run_cross_selection() returned, status:", result.get("status"), flush=True)

            self.finished.emit(
                result
            )

        except Exception as e:

            print("[CROSS WORKER] error:", e, flush=True)

            self.error.emit(
                str(e)
            )


# ============================================================
# WINDOW
# ============================================================

class CrossWindow(QMainWindow):

    def __init__(self, parent=None):

        super().__init__(
            parent
        )

        self.setWindowTitle(
            "TerraTools — Cross Selection"
        )

        self.showMaximized()

        self.browser = QWebEngineView()

        self.setCentralWidget(
            self.browser
        )

        # Local file:// results (plots, CSVs) live outside the html's own
        # directory, so local content needs permission to load them.
        # This was missing here (present in every other tool window),
        # which blocked plot images from ever rendering at all.
        settings = self.browser.settings()
        settings.setAttribute(QWebEngineSettings.LocalContentCanAccessFileUrls, True)
        settings.setAttribute(QWebEngineSettings.LocalContentCanAccessRemoteUrls, False)

        # ----------------------------------------------------
        # WebChannel
        # ----------------------------------------------------

        self.channel = QWebChannel(
            self.browser.page()
        )

        self.bridge = CrossBridge(
            self
        )

        self.channel.registerObject(
            "crossBridge",
            self.bridge
        )

        self.browser.page().setWebChannel(
            self.channel
        )

        # ----------------------------------------------------
        # Load Cross Selection HTML
        # ----------------------------------------------------

        load_frontend(self.browser, "cross_selection")


# ============================================================
# STANDALONE TEST
# ============================================================

if __name__ == "__main__":

    from PySide6.QtWidgets import QApplication

    app = QApplication(
        sys.argv
    )

    window = CrossWindow()

    window.show()

    sys.exit(
        app.exec()
    )
