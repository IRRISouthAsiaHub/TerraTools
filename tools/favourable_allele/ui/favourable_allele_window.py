import json
import shutil

from PySide6.QtCore import QObject, Qt, QThread, Signal, Slot
from tools.common.single_run import single_run
from PySide6.QtWidgets import QFileDialog

from tools.common.save_dialog import save_dialog
from tools.favourable_allele.backend import balloon_runner as runner
from tools.favourable_allele.backend import matcher


class _Worker(QObject):
    """Runs the R script on a background thread. The result travels as a JSON string (a dict through a Signal loses its order)."""

    finished = Signal(str)

    def __init__(self, payload: dict):
        super().__init__()
        self.payload = payload

    @Slot()
    def run(self):
        self.finished.emit(json.dumps(runner.draw(self.payload)))


class FavourableAlleleBridge(QObject):
    """Python <-> JavaScript bridge for the Heatmap tab of the Favourable Allele Analysis page."""

    resultReady = Signal(str)
    statusUpdate = Signal(str)

    def __init__(self, parent_window=None):
        super().__init__()
        self._parent_window = parent_window
        self._thread = None
        self._worker = None

    @Slot(result=str)
    def browse_file(self):
        p, _ = QFileDialog.getOpenFileName(self._parent_window, "Select your 4K haplopanel genotype file", "",
                                           "Genotype files (*.csv *.txt *.tsv *.vcf *.vcf.gz);;All Files (*)")
        return p or ""

    @Slot(result=str)
    def get_info(self):
        return json.dumps(matcher.info())

    @Slot(str, result=str)
    def load_genotypes(self, path: str):
        """Match the user's file against the built-in trait-gene markers; returns the matched data as JSON."""
        self.statusUpdate.emit("Matching the marker IDs…")
        result = matcher.load(path)
        if result.get("status") == "success":
            m = result["match"]
            self.statusUpdate.emit(f"Matched {m['matched_markers']} of {m['panel_markers']} trait-gene markers in {m['lines']} lines.")
        else:
            self.statusUpdate.emit("Error: " + result.get("message", "")[:200])
        return json.dumps(result)

    @Slot(str, str)
    def save_file(self, source_path: str, suggested_name: str):
        try:
            dest, _ = save_dialog(self._parent_window, suggested_name)
            if not dest:
                return
            shutil.copy2(source_path, dest)
            self.statusUpdate.emit(f"Saved to {dest}")
        except Exception as e:
            self.statusUpdate.emit(f"Save failed: {e}")

    @Slot(str)
    @single_run
    def run_heatmap(self, payload_json: str):
        try:
            payload = json.loads(payload_json)
        except Exception as e:
            self.resultReady.emit(json.dumps({"status": "error", "message": f"Invalid request: {e}"}))
            return
        self.statusUpdate.emit("Drawing the balloon heatmap in R…")
        self._thread = QThread()
        self._worker = _Worker(payload)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run, Qt.QueuedConnection)
        self._worker.finished.connect(self._on_finished, Qt.QueuedConnection)
        self._worker.finished.connect(self._thread.quit, Qt.QueuedConnection)
        self._worker.finished.connect(self._worker.deleteLater, Qt.QueuedConnection)
        self._thread.finished.connect(self._thread.deleteLater, Qt.QueuedConnection)
        self._thread.start()

    @Slot(str)
    def _on_finished(self, payload: str):
        try:
            result = json.loads(payload)
        except Exception:
            result = {"status": "error", "message": "The result could not be read."}
            payload = json.dumps(result)
        if result.get("status") == "success":
            self.statusUpdate.emit(f"Done — {result.get('n_genes')} genes × {result.get('n_lines')} lines.")
        else:
            self.statusUpdate.emit(f"Error: {result.get('message', 'Unknown error')[:200]}")
        self.resultReady.emit(payload)
