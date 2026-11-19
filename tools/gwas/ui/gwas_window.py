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
from tools.gwas.config import GWASConfig
from ui.frontend import load_frontend


class _GWASWorker(QObject):
    """Runs GWASEngine.run() on a background thread so the UI never freezes."""

    finished = Signal(dict)

    def __init__(self, config: GWASConfig):
        super().__init__()
        self.config = config

    @Slot()
    def run(self):
        print("[GWAS WORKER] run() called on worker thread", flush=True)
        # Imported here, not at module level: the engine is heavy and every
        # tool bridge is created at app startup.
        from tools.gwas.backend.gwas_engine import GWASEngine
        engine = GWASEngine(self.config)
        print("[GWAS WORKER] GWASEngine created, calling engine.run()...", flush=True)
        result = engine.run()
        print("[GWAS WORKER] engine.run() returned, status:", result.get("status"), flush=True)
        self.finished.emit(result)
        print("[GWAS WORKER] finished signal emitted", flush=True)


class GWASBridge(QObject):
    """Python <-> JavaScript bridge for the GWAS HTML frontend."""

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
            "Select GWAS Input File",
            "",
            "Data Files (*.csv *.tsv *.txt);;All Files (*)"
        )
        return file_path or ""

    @Slot(result=str)
    def browse_phenotype(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self._parent_window, "Select phenotype file", "",
            "Phenotype files (*.csv *.tsv *.txt *.xlsx *.xls);;All Files (*)"
        )
        return file_path or ""

    @Slot(result=str)
    def browse_vcf(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self._parent_window, "Select VCF file", "",
            "VCF (*.vcf *.vcf.gz);;All Files (*)"
        )
        return file_path or ""

    @Slot(str, result=str)
    def get_phenotype_columns(self, path: str):
        """Trait columns (numeric columns after the ID column) of a phenotype file, as JSON."""
        try:
            from tools.gwas.backend.gwas_inputs import phenotype_columns
            return json.dumps(phenotype_columns(path))
        except Exception as e:
            self.statusUpdate.emit(f"Could not read the phenotype file: {e}")
            return json.dumps([])

    @Slot(result=str)
    def rice3k_info(self):
        """Size of the bundled rice 3K panel, example IDs, and a ready-to-fill phenotype template
        (all accession IDs + an empty Trait column), as JSON."""
        try:
            import tempfile
            from tools.gwas.backend.gwas_inputs import rice3k_info, rice3k_template
            info = rice3k_info()
            info.pop("ids", None)
            template = Path(tempfile.mkdtemp(prefix="terratools_3k_")) / "rice3k_phenotype_template.csv"
            info["template"] = rice3k_template(template)
            return json.dumps(info)
        except Exception as e:
            self.statusUpdate.emit(f"Rice 3K panel unavailable: {e}")
            return json.dumps({})

    @Slot(str, str)
    def save_file(self, source_path: str, suggested_name: str):
        """Opens a native Save As dialog and copies source_path to the chosen location.
        Used for the Excel report / plots ZIP download buttons, since <a download>
        does not reliably trigger a save dialog for local file:// URLs inside
        QWebEngineView."""
        try:
            dest_path, _ = save_dialog(self._parent_window, suggested_name)
            if not dest_path:
                return  # user cancelled
            shutil.copy2(source_path, dest_path)
            self.statusUpdate.emit(f"Saved to {dest_path}")
        except Exception as e:
            import traceback
            traceback.print_exc()
            self.statusUpdate.emit(f"Save failed: {e}")

    @Slot(str)
    @single_run
    def run_gwas(self, config_json: str):
        try:
            config_dict = json.loads(config_json)
        except Exception as e:
            self.resultReady.emit(json.dumps({
                "status": "error",
                "message": f"Invalid configuration: {e}"
            }))
            return

        mode = config_dict.get("input_mode", "file")
        if mode == "file" and not config_dict.get("input_file"):
            self.resultReady.emit(json.dumps({
                "status": "error",
                "message": "Please select an input file."
            }))
            return
        if mode in ("3k", "vcf") and not (config_dict.get("phenotype_file") and config_dict.get("trait")):
            self.resultReady.emit(json.dumps({
                "status": "error",
                "message": "Please select the phenotype file and the trait to analyse."
            }))
            return
        if mode == "vcf" and not config_dict.get("vcf_file"):
            self.resultReady.emit(json.dumps({
                "status": "error",
                "message": "Please select the VCF file."
            }))
            return

        if not config_dict.get("models"):
            self.resultReady.emit(json.dumps({
                "status": "error",
                "message": "Please select at least one GWAS model."
            }))
            return

        config = GWASConfig.from_dict(config_dict)
        print("\n========== GWAS CONFIGURATION ==========")
        print(config)
        print("========================================\n")

        self.statusUpdate.emit("Running GWAS...")

        self._thread = QThread()
        self._worker = _GWASWorker(config)
        self._worker.moveToThread(self._thread)

        self._thread.started.connect(self._worker.run, Qt.QueuedConnection)
        self._worker.finished.connect(self._on_finished, Qt.QueuedConnection)
        self._worker.finished.connect(self._thread.quit, Qt.QueuedConnection)
        self._worker.finished.connect(self._worker.deleteLater, Qt.QueuedConnection)
        self._thread.finished.connect(self._thread.deleteLater, Qt.QueuedConnection)

        print("[GWAS BRIDGE] about to call self._thread.start()", flush=True)
        self._thread.start()
        print("[GWAS BRIDGE] self._thread.start() returned", flush=True)

    @Slot(dict)
    def _on_finished(self, result: dict):
        print("[GWAS BRIDGE] _on_finished() called, status:", result.get("status"), flush=True)
        if result.get("status") == "success":
            self.statusUpdate.emit("GWAS completed successfully.")
        else:
            self.statusUpdate.emit(f"Error: {result.get('message', 'Unknown error')}")
        print("[GWAS BRIDGE] about to emit resultReady...", flush=True)
        self.resultReady.emit(json.dumps(result))
        print("[GWAS BRIDGE] resultReady emitted", flush=True)


class GWASWindow(QMainWindow):

    def __init__(self, parent=None):
        super().__init__(parent)

        self.setWindowTitle("TerraTools — GWAS Analysis")
        self.resize(1400, 900)

        self.browser = QWebEngineView()
        self.setCentralWidget(self.browser)

        # Local file:// results (plots, excel, zip) live outside the html's own
        # directory, so local content needs permission to load them.
        settings = self.browser.settings()
        settings.setAttribute(QWebEngineSettings.LocalContentCanAccessFileUrls, True)
        settings.setAttribute(QWebEngineSettings.LocalContentCanAccessRemoteUrls, False)

        self.bridge = GWASBridge(parent_window=self)
        self.channel = QWebChannel()
        self.channel.registerObject("gwasBridge", self.bridge)
        self.browser.page().setWebChannel(self.channel)

        load_frontend(self.browser, "gwas")
