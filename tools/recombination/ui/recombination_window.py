import json
import shutil
import tempfile

from PySide6.QtCore import QObject, Qt, QThread, Signal, Slot
from PySide6.QtWidgets import QFileDialog

from tools.common.save_dialog import save_dialog
from tools.common.single_run import single_run


class _RunWorker(QObject):
    """Runs recombination.run() off the GUI thread (smoothing and the circos figure take a few seconds)."""

    finished = Signal(dict)

    def __init__(self, config):
        super().__init__()
        self.config = config

    @Slot()
    def run(self):
        try:
            from tools.recombination.backend.recombination import run
            result = run(self.config)
        except Exception as e:
            import traceback
            traceback.print_exc()
            result = {"status": "error", "message": str(e)}
        self.finished.emit(result)


class RecombinationBridge(QObject):
    """Python <-> JavaScript bridge of the Recombination Analysis page."""

    resultReady = Signal(str)
    statusUpdate = Signal(str)

    def __init__(self, parent_window=None):
        super().__init__()
        self._parent_window = parent_window
        self._thread = None
        self._worker = None

    @Slot(result=str)
    def browse_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self._parent_window, "Select the genotype file of the mapping population", "",
            "Genotype tables (*.csv *.tsv *.txt *.xlsx *.xls);;All Files (*)")
        return path or ""

    @Slot(str, str, result=str)
    def inspect_file(self, path, sheet):
        """Samples, marker count, a guess of the recurrent parent and a preview (the file is read once)."""
        try:
            from tools.recombination.backend.recombination import RecombError, inspect_file
            return json.dumps(inspect_file(path, sheet or None), ensure_ascii=False, default=str)
        except Exception as e:
            from tools.recombination.backend.recombination import RecombError
            if not isinstance(e, RecombError):
                import traceback
                traceback.print_exc()
            return json.dumps({"status": "error", "message": str(e) if isinstance(e, RecombError) else f"Could not read the file: {e}"})

    @Slot(str)
    @single_run
    def run_recombination(self, config_json):
        try:
            config = json.loads(config_json)
        except Exception as e:
            self.resultReady.emit(json.dumps({"status": "error", "message": f"Invalid configuration: {e}"}))
            return
        config["output_dir"] = tempfile.mkdtemp(prefix="terratools_recomb_")
        self.statusUpdate.emit("Detecting crossovers…")
        self._thread = QThread()
        self._worker = _RunWorker(config)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run, Qt.QueuedConnection)
        self._worker.finished.connect(self._finished, Qt.QueuedConnection)
        self._worker.finished.connect(self._thread.quit, Qt.QueuedConnection)
        self._worker.finished.connect(self._worker.deleteLater, Qt.QueuedConnection)
        self._thread.finished.connect(self._thread.deleteLater, Qt.QueuedConnection)
        self._thread.start()

    @Slot(dict)
    def _finished(self, result):
        self.resultReady.emit(json.dumps(result, ensure_ascii=False, default=str))

    @Slot(str, str)
    def save_file(self, source_path, suggested_name):
        """Native Save As + copy (downloads of the tables and figures)."""
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
