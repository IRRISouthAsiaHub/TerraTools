import json
import shutil
from pathlib import Path

from PySide6.QtCore import QObject, Qt, QThread, QUrl, Signal, Slot
from tools.common.single_run import single_run
from PySide6.QtWebChannel import QWebChannel
from PySide6.QtWebEngineCore import QWebEngineSettings
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QFileDialog, QMainWindow

from tools.common import json_safe
from tools.common.save_dialog import save_dialog
from tools.genomic_prediction.config import GenomicPredictionConfig
from tools.genomic_prediction.backend.genomic_prediction_runner import GenomicPredictionEngine
from ui.frontend import load_frontend


class _GenomicPredictionWorker(QObject):
    """Runs GenomicPredictionEngine.run() on a background thread. This is
    the slowest pipeline in TerraTools by a wide margin — K-fold CV, each
    fold training a DNN + 3 other models, then a final full retrain — so
    it streams log lines live via logLine rather than only reporting at
    the very end."""

    finished = Signal(str)   # JSON text: a dict through a Qt signal loses its key order (columns get sorted)
    logLine = Signal(str)

    def __init__(self, config: GenomicPredictionConfig):
        super().__init__()
        self.config = config

    @Slot()
    def run(self):
        engine = GenomicPredictionEngine(self.config, on_log=self._on_log)
        result = engine.run()
        self.finished.emit(json_safe.dumps(result))

    def _on_log(self, msg: str):
        # Called synchronously from inside engine.run() on this worker
        # thread; Signal.emit() is thread-safe and will queue delivery to
        # the bridge's thread automatically.
        self.logLine.emit(msg)


class GenomicPredictionBridge(QObject):
    """Python <-> JavaScript bridge for the Genomic Prediction frontend."""

    resultReady = Signal(str)
    statusUpdate = Signal(str)
    logUpdate = Signal(str)

    def __init__(self, parent_window=None):
        super().__init__()
        self._parent_window = parent_window
        self._thread = None
        self._worker = None

    def _browse(self, title: str) -> str:
        file_path, _ = QFileDialog.getOpenFileName(
            self._parent_window, title, "", "CSV Files (*.csv);;All Files (*)"
        )
        return file_path or ""

    @Slot(result=str)
    def browse_train_pheno(self):
        return self._browse("Select Training Phenotypic Data")

    @Slot(result=str)
    def browse_train_geno(self):
        return self._browse("Select Training Genotypic (Additive) Data")

    @Slot(result=str)
    def browse_test_pheno(self):
        return self._browse("Select Testing Phenotypic Data")

    @Slot(result=str)
    def browse_test_geno(self):
        return self._browse("Select Testing Genotypic (Additive) Data")

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
    def run_prediction(self, config_json: str):
        try:
            config_dict = json.loads(config_json)
        except Exception as e:
            self.resultReady.emit(json.dumps({"status": "error", "message": f"Invalid configuration: {e}"}))
            return

        required = ["train_pheno_file", "train_geno_file", "test_pheno_file", "test_geno_file"]
        if not all(config_dict.get(k) for k in required):
            self.resultReady.emit(json.dumps({
                "status": "error",
                "message": "Please select all four input files (training/testing x phenotypic/genotypic)."
            }))
            return

        config = GenomicPredictionConfig.from_dict(config_dict)
        print("\n===== GENOMIC PREDICTION CONFIGURATION =====")
        print(config)
        print("==============================================\n")

        self.statusUpdate.emit(
            "Running BreedSight genomic prediction pipeline "
            "(this trains a neural network per fold and can take a long time)..."
        )

        self._thread = QThread()
        self._worker = _GenomicPredictionWorker(config)
        self._worker.moveToThread(self._thread)

        # Explicit QueuedConnection on every cross-thread connection —
        # relying on Qt's auto-detection for a plain (non-@Slot) receiver
        # method caused a silent delivery failure in the GWAS tool.
        self._thread.started.connect(self._worker.run, Qt.QueuedConnection)
        self._worker.logLine.connect(self._on_log_line, Qt.QueuedConnection)
        self._worker.finished.connect(self._on_finished, Qt.QueuedConnection)
        self._worker.finished.connect(self._thread.quit, Qt.QueuedConnection)
        self._worker.finished.connect(self._worker.deleteLater, Qt.QueuedConnection)
        self._thread.finished.connect(self._thread.deleteLater, Qt.QueuedConnection)

        self._thread.start()

    @Slot(str)
    def _on_log_line(self, msg: str):
        self.logUpdate.emit(msg)

    @Slot(str)
    def _on_finished(self, payload: str):
        result = json.loads(payload)
        if result.get("status") == "success":
            self.statusUpdate.emit("Analysis complete.")
        else:
            self.statusUpdate.emit(f"Error: {result.get('message', 'Unknown error')}")
        self.resultReady.emit(payload)


class GenomicPredictionWindow(QMainWindow):

    def __init__(self, parent=None):
        super().__init__(parent)

        self.setWindowTitle("TerraTools — BreedSight Genomic Prediction")
        self.resize(1500, 950)

        self.browser = QWebEngineView()
        self.setCentralWidget(self.browser)

        settings = self.browser.settings()
        settings.setAttribute(QWebEngineSettings.LocalContentCanAccessFileUrls, True)
        settings.setAttribute(QWebEngineSettings.LocalContentCanAccessRemoteUrls, False)

        self.bridge = GenomicPredictionBridge(parent_window=self)
        self.channel = QWebChannel()
        self.channel.registerObject("genomicPredictionBridge", self.bridge)
        self.browser.page().setWebChannel(self.channel)

        load_frontend(self.browser, "genomic_prediction")
