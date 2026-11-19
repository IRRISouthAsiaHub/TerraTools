import json
import shutil
from pathlib import Path

from PySide6.QtCore import QObject, Qt, QThread, QUrl, Signal, Slot
from tools.common.single_run import single_run
from PySide6.QtWebChannel import QWebChannel
from PySide6.QtWebEngineCore import QWebEngineSettings
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QFileDialog, QMainWindow

from tools.common.save_dialog import save_dialog
from tools.diversity.config import DiversityConfig
from ui.frontend import load_frontend


class _DiversityWorker(QObject):
    """Runs DiversityEngine.run() on a background thread so the UI never
    freezes (the PCA/clustering/NJ-tree computation can take a while on
    larger sample sets)."""

    finished = Signal(dict)

    def __init__(self, config: DiversityConfig):
        super().__init__()
        self.config = config

    @Slot()
    def run(self):
        # Imported here, not at module level: the engine pulls in heavy
        # plotting/stats libraries and every tool bridge is created at
        # app startup.
        from tools.diversity.backend.diversity_runner import DiversityEngine
        engine = DiversityEngine(self.config)
        result = engine.run()
        self.finished.emit(result)


class DiversityBridge(QObject):
    """Python <-> JavaScript bridge for the Diversity HTML frontend."""

    resultReady = Signal(str)
    statusUpdate = Signal(str)

    def __init__(self, parent_window=None):
        super().__init__()
        self._parent_window = parent_window
        self._thread = None
        self._worker = None

    @Slot(result=str)
    def browse_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self._parent_window,
            "Select Input Matrix CSV",
            "",
            "CSV Files (*.csv);;All Files (*)"
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
    @single_run
    def run_diversity(self, config_json: str):
        try:
            config_dict = json.loads(config_json)
        except Exception as e:
            self.resultReady.emit(json.dumps({
                "status": "error",
                "message": f"Invalid configuration: {e}"
            }))
            return

        if not config_dict.get("input_file"):
            self.resultReady.emit(json.dumps({
                "status": "error",
                "message": "Please select an input file."
            }))
            return

        config = DiversityConfig.from_dict(config_dict)
        print("\n===== DIVERSITY CONFIGURATION =====")
        print(config)
        print("====================================\n")

        self.statusUpdate.emit("Running diversity analysis...")

        self._thread = QThread()
        self._worker = _DiversityWorker(config)
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
            self.statusUpdate.emit("Diversity analysis completed successfully.")
        else:
            self.statusUpdate.emit(f"Error: {result.get('message', 'Unknown error')}")
        self.resultReady.emit(json.dumps(result))


class DiversityWindow(QMainWindow):

    def __init__(self, parent=None):
        super().__init__(parent)

        self.setWindowTitle("TerraTools — Diversity Analysis")
        self.resize(1400, 900)

        self.browser = QWebEngineView()
        self.setCentralWidget(self.browser)

        settings = self.browser.settings()
        settings.setAttribute(QWebEngineSettings.LocalContentCanAccessFileUrls, True)
        settings.setAttribute(QWebEngineSettings.LocalContentCanAccessRemoteUrls, False)

        self.bridge = DiversityBridge(parent_window=self)
        self.channel = QWebChannel()
        self.channel.registerObject("diversityBridge", self.bridge)
        self.browser.page().setWebChannel(self.channel)

        load_frontend(self.browser, "diversity")
