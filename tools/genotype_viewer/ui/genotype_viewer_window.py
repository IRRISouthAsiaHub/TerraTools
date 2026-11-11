import json
import shutil
from pathlib import Path

from PySide6.QtCore import QObject, Qt, QThread, Signal, Slot
from PySide6.QtWidgets import QFileDialog

from tools.common.save_dialog import save_dialog
from tools.genotype_viewer.config import GenotypeViewerConfig
from tools.genotype_viewer.backend import genotype_viewer_runner as runner


class _ViewerWorker(QObject):
    """Draws the images on a background thread (a big file takes a while).

    The result travels as a JSON string: a dict passed through a Qt Signal becomes a QVariantMap and loses its order."""

    finished = Signal(str)

    def __init__(self, config: GenotypeViewerConfig):
        super().__init__()
        self.config = config

    @Slot()
    def run(self):
        self.finished.emit(json.dumps(runner.GenotypeViewerEngine(self.config).run()))


class _SpeciesWorker(QObject):
    finished = Signal(str)

    def __init__(self, fasta: str, name: str):
        super().__init__()
        self.fasta, self.name = fasta, name

    @Slot()
    def run(self):
        self.finished.emit(json.dumps(runner.add_species(self.fasta, self.name)))


class GenotypeViewerBridge(QObject):
    """Python <-> JavaScript bridge for the Genotype Viewer page."""

    resultReady = Signal(str)
    statusUpdate = Signal(str)
    speciesAdded = Signal(str)

    def __init__(self, parent_window=None):
        super().__init__()
        self._parent_window = parent_window
        self._thread = None
        self._worker = None
        self._sthread = None
        self._sworker = None

    # ── dialogs ──
    @Slot(result=str)
    def browse_file(self):
        p, _ = QFileDialog.getOpenFileName(self._parent_window, "Select the genotype file (CSV or VCF)",
                                           "", "Genotype files (*.csv *.txt *.tsv *.vcf *.vcf.gz);;All Files (*)")
        return p or ""

    @Slot(result=str)
    def browse_fasta(self):
        p, _ = QFileDialog.getOpenFileName(self._parent_window, "Select the genome FASTA file", "",
                                           "FASTA (*.fa *.fasta *.fna *.fa.gz *.fasta.gz *.fna.gz);;All Files (*)")
        return p or ""

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

    # ── options and file checks ──
    @Slot(result=str)
    def get_options(self):
        return json.dumps({"species": runner.list_species(), "palettes": runner.list_palettes()})

    @Slot(str, result=str)
    def get_samples(self, path: str):
        try:
            return json.dumps({"status": "success", "samples": runner.sample_names(path)})
        except Exception as e:
            return json.dumps({"status": "error", "message": str(e)})

    # ── draw ──
    @Slot(str)
    def run_viewer(self, config_json: str):
        try:
            config = GenotypeViewerConfig.from_dict(json.loads(config_json))
        except Exception as e:
            self.resultReady.emit(json.dumps({"status": "error", "message": f"Invalid configuration: {e}"}))
            return
        self.statusUpdate.emit("Drawing the graphical genotypes…")
        self._thread = QThread()
        self._worker = _ViewerWorker(config)
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
            self.statusUpdate.emit(f"Done — {result.get('n_images')} image(s) drawn.")
        else:
            self.statusUpdate.emit(f"Error: {result.get('message', 'Unknown error')}")
        self.resultReady.emit(payload)

    # ── a new species from its genome FASTA (GenoSee's add_chromosome_length_to_database.py) ──
    @Slot(str, str)
    def add_species(self, fasta: str, name: str):
        self.statusUpdate.emit("Reading the chromosome lengths from the FASTA file…")
        self._sthread = QThread()
        self._sworker = _SpeciesWorker(fasta, name)
        self._sworker.moveToThread(self._sthread)
        self._sthread.started.connect(self._sworker.run, Qt.QueuedConnection)
        self._sworker.finished.connect(self._on_species, Qt.QueuedConnection)
        self._sworker.finished.connect(self._sthread.quit, Qt.QueuedConnection)
        self._sworker.finished.connect(self._sworker.deleteLater, Qt.QueuedConnection)
        self._sthread.finished.connect(self._sthread.deleteLater, Qt.QueuedConnection)
        self._sthread.start()

    @Slot(str)
    def _on_species(self, payload: str):
        try:
            r = json.loads(payload)
        except Exception:
            r = {"status": "error", "message": "The species could not be added."}
            payload = json.dumps(r)
        self.statusUpdate.emit(f"Added {r['added']}." if r.get("status") == "success" else f"Error: {r.get('message', '')[:200]}")
        self.speciesAdded.emit(payload)
