import os
import sys
import json
from pathlib import Path

from PySide6.QtCore import QObject, Qt, Signal, Slot, QThread, QUrl
from PySide6.QtWidgets import QFileDialog
from PySide6.QtWebChannel import QWebChannel
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QMainWindow

from tools.qtl.backend import qtl_runner


# ============================================================
# PATHS
# ============================================================

QTL_DIR = Path(__file__).resolve().parent
QTL_HTML = QTL_DIR / "qtl.html"

RESULT_DIR = QTL_DIR / "results"
RESULT_DIR.mkdir(exist_ok=True)


# ============================================================
# BRIDGE
# ============================================================

class QTLBridge(QObject):

    resultReady = Signal(str)
    statusUpdate = Signal(str)
    errorOccurred = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)

        self.last_result = None
        self.last_config = None
        self.last_compare_result = None
        self.last_common_result = None

    # --------------------------------------------------------
    # FILE BROWSING
    # --------------------------------------------------------

    @Slot(result=str)
    def browse_genotype(self):

        path, _ = QFileDialog.getOpenFileName(
            None,
            "Select Genotype CSV",
            "",
            "CSV Files (*.csv);;All Files (*)"
        )

        return path or ""

    @Slot(result=str)
    def browse_phenotype(self):

        path, _ = QFileDialog.getOpenFileName(
            None,
            "Select Phenotype CSV",
            "",
            "CSV Files (*.csv);;All Files (*)"
        )

        return path or ""

    # --------------------------------------------------------
    # PHENOTYPE COLUMNS
    # --------------------------------------------------------

    @Slot(str, result=str)
    def get_phenotype_columns(self, path):

        try:

            import pandas as pd

            df = pd.read_csv(path)

            if df.shape[1] <= 1:
                return json.dumps([])

            columns = [
                str(c)
                for c in df.columns[1:]
            ]

            return json.dumps(columns)

        except Exception as e:

            self.errorOccurred.emit(
                f"Could not read phenotype file: {e}"
            )

            return json.dumps([])

    # --------------------------------------------------------
    # RUN QTL
    # --------------------------------------------------------

    @Slot(str)
    def run_qtl(self, config_json):

        print("[QTL BRIDGE] run_qtl() called", flush=True)

        try:

            config = json.loads(config_json)

            self.last_config = config

            print("\n========== QTL CONFIGURATION ==========")
            print(config)
            print("========================================\n")

            self.statusUpdate.emit(
                "Preparing QTL analysis..."
            )

            thread = QThread()

            worker = QTLWorker(
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

            print("[QTL BRIDGE] about to call self._thread.start()", flush=True)
            thread.start()
            print("[QTL BRIDGE] self._thread.start() returned", flush=True)

        except Exception as e:

            self._analysis_error(
                str(e)
            )

    @Slot(dict)
    def _analysis_finished(self, result):

        print("[QTL BRIDGE] _analysis_finished() called, status:", result.get("status"), flush=True)

        self.last_result = result

        self.statusUpdate.emit(
            "QTL analysis completed successfully."
        )

        print("[QTL BRIDGE] about to emit resultReady...", flush=True)
        self.resultReady.emit(
            json.dumps(result)
        )
        print("[QTL BRIDGE] resultReady emitted", flush=True)

    @Slot(str)
    def _analysis_error(self, message):

        print("[QTL BRIDGE] _analysis_error() called:", message, flush=True)

        result = {
            "status": "error",
            "message": message
        }

        self.last_result = result

        self.statusUpdate.emit(
            "QTL analysis failed."
        )

        self.resultReady.emit(
            json.dumps(result)
        )

    # --------------------------------------------------------
    # MODEL COMPARISON
    # --------------------------------------------------------

    @Slot(str)
    def run_compare(self, config_json):

        print("[QTL BRIDGE] run_compare() called", flush=True)

        try:

            config = json.loads(
                config_json
            )

            models = config.get(
                "models",
                []
            )

            if len(models) < 2:

                self._send_compare_error(
                    "Select at least two models."
                )

                return

            self.statusUpdate.emit(
                "Comparing QTL models..."
            )

            thread = QThread()

            worker = QTLCompareWorker(
                config,
                RESULT_DIR
            )

            worker.moveToThread(thread)

            thread.started.connect(
                worker.run, Qt.QueuedConnection
            )

            worker.finished.connect(
                self._compare_finished, Qt.QueuedConnection
            )

            worker.error.connect(
                self._compare_error, Qt.QueuedConnection
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

            self._compare_thread = thread
            self._compare_worker = worker

            print("[QTL BRIDGE] about to call self._compare_thread.start()", flush=True)
            thread.start()
            print("[QTL BRIDGE] self._compare_thread.start() returned", flush=True)

        except Exception as e:

            self._compare_error(
                str(e)
            )

    @Slot(dict)
    def _compare_finished(self, result):

        print("[QTL BRIDGE] _compare_finished() called, status:", result.get("status"), flush=True)

        self.last_compare_result = result

        self.statusUpdate.emit(
            "Model comparison completed."
        )

        # Tell the HTML to render comparison results.
        payload = json.dumps(result)

        print("[QTL BRIDGE] pushing showCompareResult(...) via runJavaScript", flush=True)
        self.parent().browser.page().runJavaScript(
            f"showCompareResult({payload});"
        )

    @Slot(str)
    def _compare_error(self, message):

        print("[QTL BRIDGE] _compare_error() called:", message, flush=True)

        self.statusUpdate.emit(
            "Model comparison failed."
        )

        self.parent().browser.page().runJavaScript(
            "showCompareResult(" +
            json.dumps({
                "status": "error",
                "message": message
            }) +
            ");"
        )

    def _send_compare_error(self, message):

        self._compare_error(message)

    # --------------------------------------------------------
    # COMMON QTL
    # --------------------------------------------------------

    @Slot(str)
    def find_common(self, config_json):

        print("[QTL BRIDGE] find_common() called", flush=True)

        try:

            config = json.loads(
                config_json
            )

            if not self.last_compare_result:

                print("[QTL BRIDGE] find_common() called without a prior comparison result", flush=True)

                self.parent().browser.page().runJavaScript(
                    "showCommonResult(" +
                    json.dumps({
                        "status": "error",
                        "message":
                            "Run Compare Models first."
                    }) +
                    ");"
                )

                return

            # Common-QTL calculation is handled by the
            # shared backend, from the comparison result.
            result = qtl_runner.find_common_qtl(
                config,
                self.last_compare_result,
                results_dir=RESULT_DIR
            )

            print("[QTL BRIDGE] find_common_qtl() returned, status:", result.get("status"), flush=True)

            self.last_common_result = result

            self.parent().browser.page().runJavaScript(
                "showCommonResult(" +
                json.dumps(result) +
                ");"
            )

        except Exception as e:

            print("[QTL BRIDGE] find_common() error:", e, flush=True)

            self.parent().browser.page().runJavaScript(
                "showCommonResult(" +
                json.dumps({
                    "status": "error",
                    "message": str(e)
                }) +
                ");"
            )

    # --------------------------------------------------------
    # PLOT SETTINGS
    # --------------------------------------------------------

    @Slot(str)
    def save_plot_settings(self, config_json):

        try:

            config = json.loads(
                config_json
            )

            self.last_plot_settings = config

            self.statusUpdate.emit(
                "Plot settings saved."
            )

        except Exception as e:

            self.statusUpdate.emit(
                f"Could not save settings: {e}"
            )

    # --------------------------------------------------------
    # DOWNLOADS
    # --------------------------------------------------------

    @Slot(str)
    def download_result(self, result_type):

        print("[QTL BRIDGE] download_result() called for:", result_type, flush=True)

        if not self.last_result:

            self.statusUpdate.emit(
                "Run QTL analysis first."
            )

            return

        mapping = {

            "manhattan":
                self.last_result.get(
                    "manhattan"
                ),

            "ideogram":
                self.last_result.get(
                    "ideogram"
                ),

            "qtl_csv":
                self.last_result.get(
                    "qtl_csv"
                ),

            "full_csv":
                self.last_result.get(
                    "full_csv"
                )
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

            print("[QTL BRIDGE] saved", result_type, "to", destination, flush=True)

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
# QTL WORKER
# ============================================================

class QTLWorker(QObject):

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

        print("[QTL WORKER] run() called on worker thread", flush=True)

        try:

            print("[QTL WORKER] calling qtl_runner.run_qtl()...", flush=True)

            result = qtl_runner.run_qtl(
                self.config,
                results_dir=self.output_dir
            )

            print("[QTL WORKER] qtl_runner.run_qtl() returned, status:", result.get("status"), flush=True)

            self.finished.emit(
                result
            )
            print("[QTL WORKER] finished signal emitted", flush=True)

        except Exception as e:

            print("[QTL WORKER] error:", e, flush=True)

            self.error.emit(
                str(e)
            )


# ============================================================
# COMPARISON WORKER
# ============================================================

class QTLCompareWorker(QObject):

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

        print("[QTL COMPARE WORKER] run() called on worker thread", flush=True)

        try:

            print("[QTL COMPARE WORKER] calling qtl_runner.run_comparison()...", flush=True)

            result = qtl_runner.run_comparison(
                self.config,
                results_dir=self.output_dir
            )

            print("[QTL COMPARE WORKER] qtl_runner.run_comparison() returned, status:", result.get("status"), flush=True)

            self.finished.emit(
                result
            )
            print("[QTL COMPARE WORKER] finished signal emitted", flush=True)

        except Exception as e:

            print("[QTL COMPARE WORKER] error:", e, flush=True)

            self.error.emit(
                str(e)
            )


# ============================================================
# WINDOW
# ============================================================

class QTLWindow(QMainWindow):

    def __init__(self, parent=None):

        super().__init__(
            parent
        )

        self.setWindowTitle(
            "TerraTools — QTL Mapping"
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

        self.bridge = QTLBridge(
            self
        )

        self.channel.registerObject(
            "qtlBridge",
            self.bridge
        )

        self.browser.page().setWebChannel(
            self.channel
        )

        # ----------------------------------------------------
        # Load QTL HTML
        # ----------------------------------------------------

        self.browser.load(
            QUrl.fromLocalFile(
                str(QTL_HTML)
            )
        )


# ============================================================
# STANDALONE TEST
# ============================================================

if __name__ == "__main__":

    from PySide6.QtWidgets import QApplication

    app = QApplication(
        sys.argv
    )

    window = QTLWindow()

    window.show()

    sys.exit(
        app.exec()
    )
