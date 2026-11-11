import os
import sys
import json
import shutil

from PySide6.QtCore import QObject, Qt, Signal, Slot, QThread, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import QFileDialog, QMainWindow

from tools.gene_prioritization import config
from tools.gene_prioritization.backend import gene_prioritization_runner as runner
from tools.gene_prioritization.backend.cache import Cache
from tools.gene_prioritization.backend.http import RunContext


RESULT_DIR = config.RESULT_DIR


# ============================================================
# WORKERS (the network never runs on the UI thread)
# ============================================================

class RunWorker(QObject):

    finished = Signal(dict)
    progress = Signal(str)

    def __init__(self, cfg, ctx_holder):
        super().__init__()
        self.cfg = cfg
        self.ctx_holder = ctx_holder

    @Slot()
    def run(self):
        try:
            # the engine attaches its own cache connection; the context carries the cancel flag
            ctx = RunContext(None, refresh=bool(self.cfg.get("refresh_cache")))
            self.ctx_holder["ctx"] = ctx
            result = runner.run_prioritization(
                self.cfg,
                results_dir=RESULT_DIR,
                progress=lambda pct, msg, extra=None: self.progress.emit(
                    json.dumps({"pct": pct, "message": msg, **({"api": extra["api"]} if extra and "api" in extra else {})})
                ),
                ctx=ctx,
            )
        except Exception as e:
            result = {"status": "error", "message": f"Gene Prioritization failed: {e}"}
        self.finished.emit(result)


class PingWorker(QObject):

    finished = Signal(str)

    @Slot()
    def run(self):
        try:
            self.finished.emit(json.dumps(runner.check_services()))
        except Exception as e:
            self.finished.emit(json.dumps([{"source": "Services", "state": "unavailable", "error": str(e)}]))


# ============================================================
# BRIDGE
# ============================================================

class GenePrioritizationBridge(QObject):

    resultReady = Signal(str)
    statusUpdate = Signal(str)
    progressUpdate = Signal(str)
    apiStatus = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._window = parent
        self.last_result = None
        self._threads = []
        self._ctx_holder = {}

    # --------------------------------------------------------
    # FILES AND LOOKUPS
    # --------------------------------------------------------

    @Slot(result=str)
    def browse_gene_list(self):
        path, _ = QFileDialog.getOpenFileName(
            self._window, "Select Gene ID List", "",
            "Gene lists (*.txt *.csv *.tsv *.xlsx *.xls);;All Files (*)")
        return path or ""

    @Slot(str)
    def open_link(self, url):
        """Source links of the evidence open in the user's browser, never inside the app window."""
        if isinstance(url, str) and url.startswith(("https://", "http://")):
            QDesktopServices.openUrl(QUrl(url))

    @Slot(result=str)
    def get_species(self):
        return json.dumps(config.species_options())

    @Slot(result=str)
    def cache_info(self):
        try:
            cache = Cache(config.CACHE_DB)
            info = cache.stats()
            cache.close()
            return json.dumps({"status": "success", **info, "path": str(config.CACHE_DB)})
        except Exception as e:
            return json.dumps({"status": "error", "message": str(e)})

    @Slot(result=str)
    def clear_cache(self):
        try:
            cache = Cache(config.CACHE_DB)
            cache.clear()
            cache.close()
            self.statusUpdate.emit("Cached API responses were deleted.")
            return json.dumps({"status": "success"})
        except Exception as e:
            return json.dumps({"status": "error", "message": str(e)})

    # --------------------------------------------------------
    # SERVICE CHECK
    # --------------------------------------------------------

    @Slot()
    def check_apis(self):
        self.statusUpdate.emit("Checking the connection to Ensembl Plants, Gramene, Plant Reactome and NCBI…")
        worker = PingWorker()
        self._start(worker, worker.finished, self._on_ping)

    @Slot(str)
    def _on_ping(self, payload):
        self.apiStatus.emit(payload)
        self.statusUpdate.emit("Connection check finished.")

    # --------------------------------------------------------
    # RUN / CANCEL
    # --------------------------------------------------------

    @Slot(str)
    def run_prioritization(self, config_json):
        try:
            cfg = json.loads(config_json)
        except Exception as e:
            self._failed(f"Could not read the settings: {e}")
            return

        self._ctx_holder.clear()
        self.statusUpdate.emit("Starting…")
        worker = RunWorker(cfg, self._ctx_holder)
        worker.progress.connect(self._on_progress, Qt.QueuedConnection)
        self._start(worker, worker.finished, self._on_finished)

    @Slot()
    def cancel(self):
        ctx = self._ctx_holder.get("ctx")
        if ctx is not None:
            ctx.cancel.set()
            self.statusUpdate.emit("Cancelling… the requests in progress are finishing.")

    def _start(self, worker, done_signal, done_slot):
        thread = QThread()
        worker.moveToThread(thread)
        thread.started.connect(worker.run, Qt.QueuedConnection)
        done_signal.connect(done_slot, Qt.QueuedConnection)
        done_signal.connect(thread.quit, Qt.QueuedConnection)
        thread.finished.connect(thread.deleteLater, Qt.QueuedConnection)
        thread.finished.connect(lambda t=thread, w=worker: self._release(t, w), Qt.QueuedConnection)
        self._threads.append((thread, worker))
        thread.start()

    def _release(self, thread, worker):
        self._threads = [(t, w) for t, w in self._threads if w is not worker]

    @Slot(str)
    def _on_progress(self, payload):
        self.progressUpdate.emit(payload)
        try:
            self.statusUpdate.emit(json.loads(payload).get("message", ""))
        except Exception:
            pass

    @Slot(dict)
    def _on_finished(self, result):
        self.last_result = result if result.get("status") == "success" else self.last_result
        status = result.get("status")
        if status == "success":
            self.statusUpdate.emit("Gene prioritization finished. Open the tabs above to see the ranking and the evidence.")
        elif status == "cancelled":
            self.statusUpdate.emit("Cancelled.")
        else:
            self.statusUpdate.emit("Gene prioritization failed.")
        self.resultReady.emit(json.dumps(result))

    def _failed(self, message):
        self.statusUpdate.emit(message)
        self.resultReady.emit(json.dumps({"status": "error", "message": message}))

    # --------------------------------------------------------
    # DOWNLOADS
    # --------------------------------------------------------

    @Slot(str)
    def download_result(self, result_type):
        if not self.last_result or self.last_result.get("status") != "success":
            self.statusUpdate.emit("Run the prioritization first.")
            return

        source = (self.last_result.get("files") or {}).get(result_type)
        if not source or not os.path.exists(source):
            self.statusUpdate.emit("Requested result is not available.")
            return

        filters = {"csv": "CSV Files (*.csv);;All Files (*)",
                   "xlsx": "Excel Files (*.xlsx);;All Files (*)",
                   "html": "HTML Files (*.html);;All Files (*)"}
        destination, _ = QFileDialog.getSaveFileName(
            self._window, "Save Result", os.path.basename(source), filters.get(result_type, "All Files (*)"))
        if not destination:
            return

        try:
            shutil.copy2(source, destination)
            self.statusUpdate.emit("Result saved.")
        except Exception as e:
            self.statusUpdate.emit(f"Download failed: {e}")


# ============================================================
# STANDALONE WINDOW (not used by the main app, kept for testing)
# ============================================================

class GenePrioritizationWindow(QMainWindow):

    def __init__(self, parent=None):
        super().__init__(parent)

        from PySide6.QtWebChannel import QWebChannel
        from PySide6.QtWebEngineWidgets import QWebEngineView
        from ui.frontend import load_frontend

        self.setWindowTitle("TerraTools — Gene Prioritization")
        self.showMaximized()

        self.browser = QWebEngineView()
        self.setCentralWidget(self.browser)

        self.channel = QWebChannel(self.browser.page())
        self.bridge = GenePrioritizationBridge(self)
        self.channel.registerObject("genePrioritizationBridge", self.bridge)
        self.browser.page().setWebChannel(self.channel)

        load_frontend(self.browser, "gene_prioritization")


if __name__ == "__main__":
    from PySide6.QtWidgets import QApplication

    app = QApplication(sys.argv)
    window = GenePrioritizationWindow()
    window.show()
    sys.exit(app.exec())
