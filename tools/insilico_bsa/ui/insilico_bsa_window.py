import json
import shutil
import tempfile

from PySide6.QtCore import QObject, Qt, QThread, Signal, Slot
from PySide6.QtWidgets import QFileDialog

from tools.common.save_dialog import save_dialog
from tools.common.single_run import single_run
from tools.insilico_bsa.backend.json_util import json_safe


class _RunWorker(QObject):
    """Runs insilico_bsa.run() off the GUI thread (the permutations take a few seconds to a minute)."""

    finished = Signal(dict)
    progress = Signal(str)

    def __init__(self, config):
        super().__init__()
        self.config = config

    @Slot()
    def run(self):
        try:
            from tools.insilico_bsa.backend.insilico_bsa import run
            result = run(self.config, log=self.progress.emit)
        except Exception as e:
            import traceback
            traceback.print_exc()
            result = {"status": "error", "message": str(e)}
        self.finished.emit(result)


class InSilicoBsaBridge(QObject):
    """Python <-> JavaScript bridge of the In silico BSA page."""

    resultReady = Signal(str)
    statusUpdate = Signal(str)

    def __init__(self, parent_window=None):
        super().__init__()
        self._parent_window = parent_window
        self._thread = None
        self._worker = None

    @Slot(result=str)
    def browse_geno(self):
        path, _ = QFileDialog.getOpenFileName(
            self._parent_window, "Select the genotype file of the F2 / mapping population", "",
            "Genotype tables (*.csv *.tsv *.txt *.xlsx *.xls);;All Files (*)")
        return path or ""

    @Slot(result=str)
    def browse_pheno(self):
        path, _ = QFileDialog.getOpenFileName(
            self._parent_window, "Select the phenotype file (sample name, then the trait columns)", "",
            "Phenotype tables (*.csv *.tsv *.txt *.xlsx *.xls);;All Files (*)")
        return path or ""

    @Slot(result=str)
    def browse_gff(self):
        path, _ = QFileDialog.getOpenFileName(
            self._parent_window, "Select the genome annotation (GFF3 / GTF)", "",
            "Annotation (*.gff *.gff3 *.gtf *.gff.gz *.gff3.gz *.gtf.gz);;All Files (*)")
        return path or ""

    @Slot(str, result=str)
    def annotate_regions(self, config_json):
        """Separate step after the analysis: the genes of the annotation that overlap each significant region."""
        try:
            config = json.loads(config_json)
            config["output_dir"] = tempfile.mkdtemp(prefix="terratools_insilicobsa_genes_")
            from tools.insilico_bsa.backend.annotation import run
            return json.dumps(json_safe(run(config)), ensure_ascii=False, default=str, allow_nan=False)
        except Exception as e:
            import traceback
            traceback.print_exc()
            return json.dumps({"status": "error", "message": f"The annotation failed: {e}"})

    @Slot(str, str, result=str)
    def inspect_inputs(self, geno_path, pheno_path):
        """Samples, markers, population type, the traits found and a preview (read once, on choosing the files)."""
        from tools.insilico_bsa.backend.insilico_bsa import BsaError, inspect_inputs
        try:
            return json.dumps(json_safe(inspect_inputs(geno_path, pheno_path or None)), ensure_ascii=False, default=str, allow_nan=False)
        except Exception as e:
            if not isinstance(e, BsaError):
                import traceback
                traceback.print_exc()
            return json.dumps({"status": "error", "message": str(e) if isinstance(e, BsaError) else f"Could not read the file: {e}"})

    @Slot(str)
    @single_run
    def run_insilico_bsa(self, config_json):
        try:
            config = json.loads(config_json)
        except Exception as e:
            self.resultReady.emit(json.dumps({"status": "error", "message": f"Invalid configuration: {e}"}))
            return
        config["output_dir"] = tempfile.mkdtemp(prefix="terratools_insilicobsa_")
        self.statusUpdate.emit("Pooling the lines and computing the statistics…")
        self._thread = QThread()
        self._worker = _RunWorker(config)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run, Qt.QueuedConnection)
        self._worker.progress.connect(self.statusUpdate, Qt.QueuedConnection)
        self._worker.finished.connect(self._finished, Qt.QueuedConnection)
        self._worker.finished.connect(self._thread.quit, Qt.QueuedConnection)
        self._worker.finished.connect(self._worker.deleteLater, Qt.QueuedConnection)
        self._thread.finished.connect(self._thread.deleteLater, Qt.QueuedConnection)
        self._thread.start()

    @Slot(dict)
    def _finished(self, result):
        self.resultReady.emit(json.dumps(json_safe(result), ensure_ascii=False, default=str, allow_nan=False))

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
