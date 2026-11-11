import json
import shutil
from pathlib import Path

from PySide6.QtCore import QObject, Qt, QThread, Signal, Slot
from PySide6.QtWebChannel import QWebChannel
from PySide6.QtWebEngineCore import QWebEngineSettings
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QFileDialog, QMainWindow

from tools.common.save_dialog import save_dialog
from tools.phenotypic.config import PhenotypicConfig
from tools.phenotypic.backend.phenotypic_runner import PhenotypicEngine
from ui.frontend import load_frontend


class _Worker(QObject):
    """Runs one PhenotypicEngine call (it shells out to R) on a background thread so the window stays responsive."""

    # The result travels as a JSON string: a dict crossing a Qt signal becomes a QVariantMap, which re-sorts every key
    # alphabetically and would scramble the column order of the tables.
    finished = Signal(int, str, str)      # job id, kind, result as JSON

    def __init__(self, job_id: int, kind: str, fn):
        super().__init__()
        self.job_id = job_id
        self.kind = kind
        self.fn = fn

    @Slot()
    def run(self):
        try:
            result = self.fn()
        except Exception as e:            # never leave the page waiting
            import traceback
            traceback.print_exc()
            result = {"status": "error", "message": str(e)}
        self.finished.emit(self.job_id, self.kind, json.dumps(result))


class PhenotypicBridge(QObject):
    """Python <-> JavaScript bridge for the Phenotypic Analysis frontend.

    The data file is read once (inspect_data); every analysis then has its own run_module call. The results of the
    analyses run so far stay in one folder (outdir) until the data settings change.
    """

    resultReady = Signal(str)       # result of run_module
    statusUpdate = Signal(str)
    dataLoaded = Signal(str)        # result of inspect_data
    exportReady = Signal(str)       # result of export_files

    def __init__(self, parent_window=None):
        super().__init__()
        self._parent_window = parent_window
        self._jobs = {}             # job id -> (thread, worker); kept alive until the thread has finished
        self._next_job = 0
        self._outdir = ""           # folder of the current results

    # ── dialogs ───────────────────────────────────────────────────────────────────────────
    @Slot(result=str)
    def browse_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self._parent_window, "Select Phenotypic Data File", "",
            "Data Files (*.xlsx *.xls *.csv);;All Files (*)")
        return file_path or ""

    @Slot(str, str)
    def save_file(self, source_path: str, suggested_name: str):
        """Native Save As dialog + copy, since <a download> is unreliable
        for local file:// URLs inside QWebEngineView."""
        try:
            dest_path, _ = save_dialog(self._parent_window, suggested_name)
            if not dest_path:
                return
            shutil.copy2(source_path, dest_path)
            self.statusUpdate.emit(f"Saved to {dest_path}")
        except Exception as e:
            import traceback
            traceback.print_exc()
            self.statusUpdate.emit(f"Save failed: {e}")

    # ── background jobs ───────────────────────────────────────────────────────────────────
    def _start(self, kind: str, fn):
        # Finished threads are dropped here, never inside the completion slot: that runs while the thread is still
        # winding down, and destroying a running QThread aborts the program.
        self._jobs = {k: v for k, v in self._jobs.items() if not v[0].isFinished()}
        self._next_job += 1
        thread = QThread()
        worker = _Worker(self._next_job, kind, fn)
        worker.moveToThread(thread)
        self._jobs[self._next_job] = (thread, worker)

        # Explicit QueuedConnection on every cross-thread connection (auto-detection failed silently in the GWAS tool).
        thread.started.connect(worker.run, Qt.QueuedConnection)
        worker.finished.connect(self._on_finished, Qt.QueuedConnection)
        worker.finished.connect(thread.quit, Qt.QueuedConnection)
        thread.start()

    def _config(self, config_json: str):
        try:
            d = json.loads(config_json)
        except Exception as e:
            return None, None, f"Invalid configuration: {e}"
        if not d.get("input_file"):
            return None, None, "Please select an input file."
        cfg = PhenotypicConfig.from_dict(d)
        if d.get("design") == "auto":            # inspect only: pick the design from the columns
            cfg.design = "auto"
        return cfg, d.get("req"), None

    @Slot(str)
    def inspect_data(self, config_json: str):
        """Reads and validates the data file; the page gets the answer through dataLoaded."""
        cfg, req, err = self._config(config_json)
        if err:
            return self.dataLoaded.emit(json.dumps({"status": "error", "message": err}))
        self.statusUpdate.emit("Reading the data file...")
        # the request id comes back with the answer, so the page can ignore an answer to an older request
        self._start("inspect", lambda: {**PhenotypicEngine(cfg).inspect(), "req": req})

    @Slot(str)
    def run_module(self, config_json: str):
        """Runs one analysis (config.module) with its settings; the page gets the answer through resultReady."""
        cfg, _, err = self._config(config_json)
        if err:
            return self.resultReady.emit(json.dumps({"status": "error", "message": err}))
        if not cfg.module:
            return self.resultReady.emit(json.dumps({"status": "error", "message": "No analysis chosen."}))
        if not cfg.outdir:                       # new data settings: the earlier results are discarded
            self._discard_results()
        self.statusUpdate.emit(f"Running {cfg.module} ...")
        self._start("run", lambda: PhenotypicEngine(cfg).run())

    @Slot(str)
    def export_files(self, request_json: str):
        """Builds a download (workbook, figures zip, everything zip); the page gets the file path through exportReady."""
        try:
            req = json.loads(request_json)
        except Exception as e:
            return self.exportReady.emit(json.dumps({"status": "error", "message": f"Invalid request: {e}"}))
        outdir, kind = req.get("outdir", ""), req.get("kind", "all")
        self.statusUpdate.emit("Preparing the file...")
        self._start("export", lambda: PhenotypicEngine().export(outdir, kind))

    @Slot(int, str, str)
    def _on_finished(self, job_id: int, kind: str, payload: str):
        result = json.loads(payload)
        if kind == "inspect":
            if result.get("status") == "success":
                self.statusUpdate.emit(
                    f"Loaded {result.get('n_plots', 0)} plots, {len(result.get('traits', []))} traits, "
                    f"{result.get('n_environments', 0)} environments.")
            else:
                self.statusUpdate.emit(f"Error: {result.get('message', 'Unknown error')}")
            self.dataLoaded.emit(payload)
        elif kind == "export":
            if result.get("status") != "success":
                self.statusUpdate.emit(f"Error: {result.get('message', 'Unknown error')}")
            self.exportReady.emit(payload)
        else:
            if result.get("status") == "success":
                self._outdir = result.get("outdir", self._outdir)
                self.statusUpdate.emit("Analysis completed.")
            else:
                self.statusUpdate.emit(f"Error: {result.get('message', 'Unknown error')}")
            self.resultReady.emit(payload)

    # ── temporary results ─────────────────────────────────────────────────────────────────
    def _discard_results(self):
        if self._outdir:
            shutil.rmtree(self._outdir, ignore_errors=True)
            self._outdir = ""

    def cleanup(self):
        self._discard_results()


class PhenotypicWindow(QMainWindow):

    def __init__(self, parent=None):
        super().__init__(parent)

        self.setWindowTitle("TerraTools — Phenotypic Data Analysis")
        self.resize(1500, 950)

        self.browser = QWebEngineView()
        self.setCentralWidget(self.browser)

        settings = self.browser.settings()
        settings.setAttribute(QWebEngineSettings.LocalContentCanAccessFileUrls, True)
        settings.setAttribute(QWebEngineSettings.LocalContentCanAccessRemoteUrls, False)

        self.bridge = PhenotypicBridge(parent_window=self)
        self.channel = QWebChannel()
        self.channel.registerObject("phenotypicBridge", self.bridge)
        self.browser.page().setWebChannel(self.channel)

        load_frontend(self.browser, "phenotypic")

    def closeEvent(self, event):
        self.bridge.cleanup()
        super().closeEvent(event)
