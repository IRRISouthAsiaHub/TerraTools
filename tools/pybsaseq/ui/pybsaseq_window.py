import json
import shutil
from pathlib import Path

from PySide6.QtCore import QObject, Qt, QThread, QUrl, Signal, Slot
from PySide6.QtWebChannel import QWebChannel
from PySide6.QtWebEngineCore import QWebEngineSettings
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QFileDialog, QMainWindow

from tools.common.save_dialog import save_dialog
from tools.pybsaseq.config import PyBSASeqConfig
from tools.pybsaseq.backend.pybsaseq_runner import PyBSASeqEngine
from ui.frontend import load_frontend


class _PyBSASeqWorker(QObject):
    """Runs PyBSASeqEngine.run() (shells out to the bundled script) on a
    background thread. This can be slow — thousands of simulation
    replications plus genome-wide sliding-window scans — so the UI must
    stay responsive throughout."""

    finished = Signal(dict)

    def __init__(self, config: PyBSASeqConfig):
        super().__init__()
        self.config = config

    @Slot()
    def run(self):
        engine = PyBSASeqEngine(self.config)
        result = engine.run()
        self.finished.emit(result)


class PyBSASeqBridge(QObject):
    """Python <-> JavaScript bridge for the PyBSA-Seq frontend."""

    resultReady = Signal(str)
    statusUpdate = Signal(str)

    def __init__(self, parent_window=None):
        super().__init__()
        self._parent_window = parent_window
        self._thread = None
        self._worker = None

    @Slot(result=str)
    def browse_parents_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self._parent_window, "Select Parents File", "",
            "TSV/CSV Files (*.tsv *.csv *.txt);;All Files (*)"
        )
        return file_path or ""

    @Slot(result=str)
    def browse_bulks_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self._parent_window, "Select Bulks File", "",
            "TSV/CSV Files (*.tsv *.csv *.txt);;All Files (*)"
        )
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

    @Slot(str)
    def run_analysis(self, config_json: str):
        try:
            config_dict = json.loads(config_json)
        except Exception as e:
            self.resultReady.emit(json.dumps({"status": "error", "message": f"Invalid configuration: {e}"}))
            return

        if not config_dict.get("parents_file") or not config_dict.get("bulks_file"):
            self.resultReady.emit(json.dumps({
                "status": "error",
                "message": "Please select both a parents file and a bulks file."
            }))
            return

        config = PyBSASeqConfig.from_dict(config_dict)
        print("\n===== PYBSASEQ CONFIGURATION =====")
        print(config)
        print("===================================\n")

        self.statusUpdate.emit("Running BSA-Seq analysis (this can take a while for large replication counts)...")

        self._thread = QThread()
        self._worker = _PyBSASeqWorker(config)
        self._worker.moveToThread(self._thread)

        # Explicit QueuedConnection on every cross-thread connection —
        # relying on Qt's auto-detection for a plain (non-@Slot) receiver
        # method caused a silent delivery failure in the GWAS tool.
        self._thread.started.connect(self._worker.run, Qt.QueuedConnection)
        self._worker.finished.connect(self._on_finished, Qt.QueuedConnection)
        self._worker.finished.connect(self._thread.quit, Qt.QueuedConnection)
        self._worker.finished.connect(self._worker.deleteLater, Qt.QueuedConnection)
        self._thread.finished.connect(self._thread.deleteLater, Qt.QueuedConnection)

        self._thread.start()

    @Slot(dict)
    def _on_finished(self, result: dict):
        if result.get("status") == "success":
            n = len(result.get("output_files", []))
            self.statusUpdate.emit(f"Analysis completed successfully — {n} output file(s) generated.")
        else:
            self.statusUpdate.emit(f"Error: {result.get('message', 'Unknown error')}")
        self.resultReady.emit(json.dumps(result))


class PyBSASeqWindow(QMainWindow):

    def __init__(self, parent=None):
        super().__init__(parent)

        self.setWindowTitle("TerraTools — PyBSA-Seq (BSA-Seq / G-statistic)")
        self.resize(1450, 950)

        self.browser = QWebEngineView()
        self.setCentralWidget(self.browser)

        settings = self.browser.settings()
        settings.setAttribute(QWebEngineSettings.LocalContentCanAccessFileUrls, True)
        settings.setAttribute(QWebEngineSettings.LocalContentCanAccessRemoteUrls, False)

        self.bridge = PyBSASeqBridge(parent_window=self)
        self.channel = QWebChannel()
        self.channel.registerObject("pybsaseqBridge", self.bridge)
        self.browser.page().setWebChannel(self.channel)

        load_frontend(self.browser, "pybsaseq")
