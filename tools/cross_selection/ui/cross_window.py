import os
import sys
import json
from pathlib import Path

from PySide6.QtCore import QObject, Qt, Signal, Slot, QThread, QUrl
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
    def browse_criterion(self):

        path, _ = QFileDialog.getOpenFileName(
            None,
            "Select Criterion CSV",
            "",
            "CSV Files (*.csv);;All Files (*)"
        )

        return path or ""

    @Slot(result=str)
    def browse_kmat(self):

        path, _ = QFileDialog.getOpenFileName(
            None,
            "Select Kinship Matrix CSV",
            "",
            "CSV Files (*.csv);;All Files (*)"
        )

        return path or ""

    # --------------------------------------------------------
    # DATA PREVIEW (pandas, no R needed for this)
    # --------------------------------------------------------

    @Slot(str, str, result=str)
    def preview_data(self, criterion_path, kmat_path):

        try:

            result = cross_runner.preview_data(
                criterion_path or None,
                kmat_path or None
            )

            return json.dumps(result)

        except Exception as e:

            return json.dumps({
                "status": "error",
                "message": str(e)
            })

    # --------------------------------------------------------
    # RUN CROSS SELECTION
    # --------------------------------------------------------

    @Slot(str)
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
    # TABLE FILTERING (pandas, re-uses the last full CSV)
    # --------------------------------------------------------

    @Slot(str, result=str)
    def filter_crosses(self, filter_json):

        try:

            filters = json.loads(filter_json)

            if not self.last_result or self.last_result.get("status") != "success":

                return json.dumps({
                    "status": "error",
                    "message": "Run an analysis first."
                })

            crosses_csv = self.last_result.get("crosses_csv")

            result = cross_runner.filter_crosses(
                crosses_csv,
                relatedness=filters.get("relatedness", "All"),
                performance=filters.get("performance", "All"),
                kinship_max=filters.get("kinship_max", 1.0),
                kinship_min=filters.get("kinship_min", 0.0)
            )

            return json.dumps(result)

        except Exception as e:

            return json.dumps({
                "status": "error",
                "message": str(e)
            })

    # --------------------------------------------------------
    # TABLE EXPORT (writes straight to a user-chosen file — no
    # need to round-trip through R again, the full crosses CSV
    # written by the last analysis already has everything needed)
    # --------------------------------------------------------

    @Slot(str)
    def export_filtered_crosses(self, filter_json):

        print("[CROSS BRIDGE] export_filtered_crosses() called", flush=True)

        try:

            filters = json.loads(filter_json)

            if not self.last_result or self.last_result.get("status") != "success":
                self.statusUpdate.emit("Run an analysis first.")
                return

            crosses_csv = self.last_result.get("crosses_csv")

            result = cross_runner.filter_crosses(
                crosses_csv,
                relatedness=filters.get("relatedness", "All"),
                performance=filters.get("performance", "All"),
                kinship_max=filters.get("kinship_max", 1.0),
                kinship_min=filters.get("kinship_min", 0.0)
            )

            if result.get("status") != "success":
                self.statusUpdate.emit(result.get("message", "Filter failed."))
                return

            self._save_rows_as_csv(
                result.get("rows", []),
                "filtered_crosses.csv"
            )

        except Exception as e:

            import traceback
            traceback.print_exc()

            self.statusUpdate.emit(f"Export failed: {e}")

    @Slot(str)
    def export_top_crosses(self, band):

        print("[CROSS BRIDGE] export_top_crosses() called for:", band, flush=True)

        try:

            if not self.last_result or self.last_result.get("status") != "success":
                self.statusUpdate.emit("Run an analysis first.")
                return

            crosses_csv = self.last_result.get("crosses_csv")

            # Same kinship bands already used to populate the on-screen
            # Top Crosses tables — kept in sync with those, not
            # re-derived, so an export always matches what's shown.
            if band == "low":
                kinship_min, kinship_max = 0.0, 0.125
                suggested_name = "top_crosses_low_kinship.csv"
            elif band == "medium":
                kinship_min, kinship_max = 0.125, 0.25
                suggested_name = "top_crosses_medium_kinship.csv"
            else:
                self.statusUpdate.emit(f"Unknown band: {band}")
                return

            result = cross_runner.filter_crosses(
                crosses_csv,
                relatedness="All",
                performance="All",
                kinship_max=kinship_max,
                kinship_min=kinship_min
            )

            if result.get("status") != "success":
                self.statusUpdate.emit(result.get("message", "Filter failed."))
                return

            self._save_rows_as_csv(
                result.get("rows", []),
                suggested_name
            )

        except Exception as e:

            import traceback
            traceback.print_exc()

            self.statusUpdate.emit(f"Export failed: {e}")

    def _save_rows_as_csv(self, rows, suggested_name):

        destination, _ = QFileDialog.getSaveFileName(
            None,
            "Save Table",
            suggested_name,
            "CSV Files (*.csv);;All Files (*)"
        )

        if not destination:
            return

        try:

            import pandas as pd

            pd.DataFrame(rows).to_csv(destination, index=False)

            print("[CROSS BRIDGE] saved table to", destination, flush=True)

            self.statusUpdate.emit(f"Saved to {destination}")

        except Exception as e:

            import traceback
            traceback.print_exc()

            self.statusUpdate.emit(f"Save failed: {e}")

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

        plots = self.last_result.get("plots", {})

        mapping = {

            "hist_pov": plots.get("hist_pov"),
            "hist_kinship": plots.get("hist_kinship"),
            "bar_relatedness": plots.get("bar_relatedness"),
            "scatter_vk": plots.get("scatter_vk"),
            "parent_bar": plots.get("parent_bar"),

            "crosses_csv": self.last_result.get("crosses_csv"),
            "parent_csv": self.last_result.get("parent_csv")
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
