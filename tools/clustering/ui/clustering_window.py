import json
import shutil
from pathlib import Path

from PySide6.QtCore import QObject, Qt, QThread, QUrl, Signal, Slot
from PySide6.QtWebChannel import QWebChannel
from PySide6.QtWebEngineCore import QWebEngineSettings
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QFileDialog, QMainWindow

from tools.common.save_dialog import save_dialog
from tools.clustering.config import ClusteringConfig
from tools.clustering.backend.clustering_runner import ClusteringEngine
from ui.frontend import load_frontend


class _ClusteringWorker(QObject):
    """Runs ClusteringEngine.run() on a background thread so the UI never
    freezes during PCA/K-Means.

    The result travels as a JSON string: a Python dict passed through a Qt
    Signal becomes a QVariantMap, which re-sorts every key alphabetically and
    would scramble the column order of the membership table."""

    finished = Signal(str)

    def __init__(self, config: ClusteringConfig):
        super().__init__()
        self.config = config

    @Slot()
    def run(self):
        engine = ClusteringEngine(self.config)
        result = engine.run()
        self.finished.emit(json.dumps(result))


class ClusteringBridge(QObject):
    """Python <-> JavaScript bridge for the Clustering HTML frontend."""

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
            "Select Genotype Data File",
            "",
            "Data Files (*.csv *.tsv *.txt);;All Files (*)"
        )
        return file_path or ""

    @Slot(result=str)
    def browse_pheno_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self._parent_window,
            "Select Phenotype File (sample IDs + trait)",
            "",
            "Phenotype Files (*.csv *.tsv *.txt *.xlsx);;All Files (*)"
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
    def run_clustering(self, config_json: str):
        try:
            config_dict = json.loads(config_json)
        except Exception as e:
            self.resultReady.emit(json.dumps({"status": "error", "message": f"Invalid configuration: {e}"}))
            return

        if not config_dict.get("input_file"):
            self.resultReady.emit(json.dumps({
                "status": "error",
                "message": "Please select an input file."
            }))
            return

        config = ClusteringConfig.from_dict(config_dict)
        print("\n===== CLUSTERING CONFIGURATION =====")
        print(config)
        print("=====================================\n")

        self.statusUpdate.emit("Running PCA + K-Means clustering and train/test split...")

        self._thread = QThread()
        self._worker = _ClusteringWorker(config)
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

    @Slot(str)
    def _on_finished(self, result_json: str):
        try:
            result = json.loads(result_json)
        except Exception:
            result = {"status": "error", "message": "The clustering result could not be read."}
            result_json = json.dumps(result)
        if result.get("status") == "success":
            self.statusUpdate.emit(
                f"Done — {result.get('n_clusters_used')} clusters, "
                f"{result.get('n_train')} train / {result.get('n_test')} test samples."
            )
        else:
            self.statusUpdate.emit(f"Error: {result.get('message', 'Unknown error')}")
        self.resultReady.emit(result_json)


class ClusteringWindow(QMainWindow):

    def __init__(self, parent=None):
        super().__init__(parent)

        self.setWindowTitle("TerraTools — Genotypic Clustering & Genomic Selection")
        self.resize(1400, 900)

        self.browser = QWebEngineView()
        self.setCentralWidget(self.browser)

        settings = self.browser.settings()
        settings.setAttribute(QWebEngineSettings.LocalContentCanAccessFileUrls, True)
        settings.setAttribute(QWebEngineSettings.LocalContentCanAccessRemoteUrls, False)

        self.bridge = ClusteringBridge(parent_window=self)
        self.channel = QWebChannel()
        self.channel.registerObject("clusteringBridge", self.bridge)
        self.browser.page().setWebChannel(self.channel)

        load_frontend(self.browser, "clustering")
