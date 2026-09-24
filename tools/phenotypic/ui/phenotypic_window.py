import json
import shutil
from pathlib import Path

from PySide6.QtCore import QObject, Qt, QThread, QUrl, Signal, Slot
from PySide6.QtWebChannel import QWebChannel
from PySide6.QtWebEngineCore import QWebEngineSettings
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QFileDialog, QMainWindow

from tools.phenotypic.config import PhenotypicConfig
from tools.phenotypic.backend.phenotypic_runner import PhenotypicEngine
from ui.frontend import load_frontend


class _PhenotypicWorker(QObject):
    """Runs PhenotypicEngine.run() (shells out to R) on a background thread.
    This pipeline fits many mixed models in one pass, so it can take a
    while — the UI must stay responsive throughout."""

    finished = Signal(dict)

    def __init__(self, config: PhenotypicConfig):
        super().__init__()
        self.config = config

    @Slot()
    def run(self):
        engine = PhenotypicEngine(self.config)
        result = engine.run()
        self.finished.emit(result)


class PhenotypicBridge(QObject):
    """Python <-> JavaScript bridge for the Phenotypic Analysis frontend."""

    resultReady = Signal(str)
    statusUpdate = Signal(str)
    dataLoaded = Signal(str)

    def __init__(self, parent_window=None):
        super().__init__()
        self._parent_window = parent_window
        self._thread = None
        self._worker = None

    @Slot(result=str)
    def browse_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self._parent_window,
            "Select Phenotypic Data File",
            "",
            "Data Files (*.xlsx *.xls *.csv);;All Files (*)"
        )
        if not file_path:
            return ""
        preview = self._preview_file(file_path)
        self.dataLoaded.emit(json.dumps(preview))
        return file_path

    def _preview_file(self, file_path: str) -> dict:
        """Reads just enough of the file to populate the trait/environment
        dropdowns before the user clicks Run — mirrors app.R's file-upload
        observer, which updates selectInputs immediately on upload."""
        try:
            import pandas as pd
            if file_path.lower().endswith(".csv"):
                df = pd.read_csv(file_path, nrows=5000)
            else:
                df = pd.read_excel(file_path, nrows=5000)

            envs = sorted(map(str, df["Environment"].dropna().unique())) if "Environment" in df.columns else []
            trait_cols = list(df.columns[9:]) if len(df.columns) > 9 else []
            numeric_traits = [c for c in trait_cols if pd.api.types.is_numeric_dtype(df[c])]

            return {
                "status": "ok",
                "n_vars": int(df.shape[1]),
                "n_obs": int(df.shape[0]),
                "n_geno": int(df["Genotype"].nunique()) if "Genotype" in df.columns else None,
                "n_env": len(envs),
                "environments": envs,
                "traits": numeric_traits,
            }
        except Exception as e:
            return {"status": "error", "message": str(e)}

    @Slot(str, str)
    def save_file(self, source_path: str, suggested_name: str):
        """Native Save As dialog + copy, since <a download> is unreliable
        for local file:// URLs inside QWebEngineView."""
        try:
            dest_path, _ = QFileDialog.getSaveFileName(
                self._parent_window,
                "Save File",
                suggested_name,
                "All Files (*)"
            )
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
        if not config_dict.get("trait"):
            self.resultReady.emit(json.dumps({
                "status": "error",
                "message": "Please select a trait to analyze."
            }))
            return

        config = PhenotypicConfig.from_dict(config_dict)
        print("\n===== PHENOTYPIC CONFIGURATION =====")
        print(config)
        print("=====================================\n")

        self.statusUpdate.emit("Running phenotypic analysis (this covers 10+ models, may take a while)...")

        self._thread = QThread()
        self._worker = _PhenotypicWorker(config)
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
            self.statusUpdate.emit("Phenotypic analysis completed successfully.")
        else:
            self.statusUpdate.emit(f"Error: {result.get('message', 'Unknown error')}")
        self.resultReady.emit(json.dumps(result))


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
