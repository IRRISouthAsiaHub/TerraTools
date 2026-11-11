import json
import shutil
from pathlib import Path

from PySide6.QtCore import QObject, Qt, QThread, Signal, Slot
from PySide6.QtWidgets import QFileDialog

from tools.common import json_safe
from tools.common.save_dialog import save_dialog
from tools.et_gwas.config import EtGwasConfig


class _EtGwasWorker(QObject):
    """Runs the analysis on a background thread so the UI never freezes."""

    finished = Signal(dict)

    def __init__(self, config: EtGwasConfig):
        super().__init__()
        self.config = config

    @Slot()
    def run(self):
        # Imported here, not at module level: every tool bridge is created at app start-up.
        from tools.et_gwas.backend.et_gwas_engine import run_et_gwas
        try:
            result = run_et_gwas(self.config, log=lambda m: print(m, flush=True))
        except Exception as e:                 # never leave the page spinning
            import traceback
            traceback.print_exc()
            result = {"status": "error", "message": f"Et-GWAS failed: {e}"}
        self.finished.emit(result)


class EtGwasBridge(QObject):
    """Python <-> JavaScript bridge for the ET-GWAS page."""

    resultReady = Signal(str)
    statusUpdate = Signal(str)

    def __init__(self, parent_window=None):
        super().__init__()
        self._parent_window = parent_window
        self._thread = None
        self._worker = None
        self._running = False

    # ---- file dialogs
    @Slot(result=str)
    def browse_file(self):
        p, _ = QFileDialog.getOpenFileName(self._parent_window, "Select the ready-made genotype + phenotype file", "",
                                           "Data Files (*.csv *.tsv *.txt);;All Files (*)")
        return p or ""

    @Slot(result=str)
    def browse_phenotype(self):
        p, _ = QFileDialog.getOpenFileName(self._parent_window, "Select phenotype file", "",
                                           "Phenotype files (*.csv *.tsv *.txt *.xlsx *.xls);;All Files (*)")
        return p or ""

    @Slot(result=str)
    def browse_vcf(self):
        p, _ = QFileDialog.getOpenFileName(self._parent_window, "Select VCF file", "", "VCF (*.vcf *.vcf.gz);;All Files (*)")
        return p or ""

    # ---- helpers for the page
    @Slot(str, result=str)
    def get_phenotype_columns(self, path: str):
        try:
            from tools.gwas.backend.gwas_inputs import phenotype_columns
            return json.dumps(phenotype_columns(path))
        except Exception as e:
            self.statusUpdate.emit(f"Could not read the phenotype file: {e}")
            return json.dumps([])

    @Slot(result=str)
    def rice3k_info(self):
        """Size of the bundled rice 3K panel, example IDs and an ID template, as JSON."""
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
        try:
            dest, _ = save_dialog(self._parent_window, suggested_name)
            if not dest:
                return
            shutil.copy2(source_path, dest)
            self.statusUpdate.emit(f"Saved to {dest}")
        except Exception as e:
            import traceback
            traceback.print_exc()
            self.statusUpdate.emit(f"Save failed: {e}")

    # ---- run
    @Slot(str)
    def run_et_gwas(self, config_json: str):
        try:
            d = json.loads(config_json)
        except Exception as e:
            self.resultReady.emit(json.dumps({"status": "error", "message": f"Invalid configuration: {e}"}))
            return
        mode = d.get("input_mode", "3k")
        problem = None
        if mode == "file" and not d.get("input_file"):
            problem = "Please select the genotype + phenotype file."
        elif mode in ("3k", "vcf") and not (d.get("phenotype_file") and d.get("trait")):
            problem = "Please select the phenotype file and the trait to analyse."
        elif mode == "vcf" and not d.get("vcf_file"):
            problem = "Please select the VCF file."
        if problem:
            self.resultReady.emit(json.dumps({"status": "error", "message": problem}))
            return
        if self._running:                      # (not QThread.isRunning(): the thread object is deleted after each run)
            self.statusUpdate.emit("An Et-GWAS run is already in progress.")
            return

        config = EtGwasConfig.from_dict(d)
        self._running = True
        self.statusUpdate.emit("Running Et-GWAS...")
        self._thread = QThread()
        self._worker = _EtGwasWorker(config)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run, Qt.QueuedConnection)
        self._worker.finished.connect(self._on_finished, Qt.QueuedConnection)
        self._worker.finished.connect(self._thread.quit, Qt.QueuedConnection)
        self._worker.finished.connect(self._worker.deleteLater, Qt.QueuedConnection)
        self._thread.finished.connect(self._thread.deleteLater, Qt.QueuedConnection)
        self._thread.start()

    @Slot(dict)
    def _on_finished(self, result: dict):
        self._running = False                  # (the thread/worker are replaced by the next run, never touched again)
        if result.get("status") == "success":
            self.statusUpdate.emit("Et-GWAS completed successfully.")
        else:
            self.statusUpdate.emit(f"Error: {result.get('message', 'Unknown error')}")
        self.resultReady.emit(json_safe.dumps(result))
