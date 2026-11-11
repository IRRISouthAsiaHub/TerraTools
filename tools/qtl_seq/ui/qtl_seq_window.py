import json
import os
import shutil
import sys
import threading

from PySide6.QtCore import QObject, Qt, QThread, QUrl, Signal, Slot
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import QFileDialog, QMainWindow

from tools.qtl_seq import config
from tools.qtl_seq.backend import deps, pipeline, validation

RESULT_DIR = config.RESULT_DIR


# ============================================================
# WORKERS (nothing heavy runs on the UI thread)
# ============================================================

class RunWorker(QObject):

    finished = Signal(dict)
    progress = Signal(str)

    def __init__(self, cfg, holder):
        super().__init__()
        self.cfg = cfg
        self.holder = holder

    @Slot()
    def run(self):
        try:
            result, res = pipeline.run_analysis(
                self.cfg, RESULT_DIR,
                progress=lambda pct, msg: self.progress.emit(json.dumps({"pct": round(pct, 1), "message": msg})),
                cancelled=self.holder["cancel"].is_set)
            self.holder["res"] = res
        except Exception as e:
            result = {"status": "error", "message": f"QTL-seq failed: {e}"}
        self.finished.emit(result)


class CheckWorker(QObject):

    finished = Signal(str)

    def __init__(self, cfg):
        super().__init__()
        self.cfg = cfg

    @Slot()
    def run(self):
        try:
            self.finished.emit(json.dumps(pipeline.check_inputs(self.cfg)))
        except Exception as e:
            self.finished.emit(json.dumps({"status": "error", "ok": False, "errors": [str(e)], "warnings": [], "files": {}}))


class ValidationWorker(QObject):

    finished = Signal(str)

    def __init__(self, key, seed, holder):
        super().__init__()
        self.key, self.seed, self.holder = key, seed, holder

    @Slot()
    def run(self):
        try:
            self.finished.emit(json.dumps(validation.run_validation(self.key, self.seed)))
        except Exception as e:
            self.finished.emit(json.dumps({"status": "error", "message": f"Validation failed: {e}"}))


# ============================================================
# BRIDGE
# ============================================================

class QTLSeqBridge(QObject):

    resultReady = Signal(str)
    statusUpdate = Signal(str)
    progressUpdate = Signal(str)
    inputsChecked = Signal(str)
    validationReady = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._window = parent
        self._threads = []
        self._holder = {"cancel": threading.Event(), "res": None}
        self.last_result = None

    # ---------------------------------------------------------- files
    def _open(self, title, filt):
        path, _ = QFileDialog.getOpenFileName(self._window, title, "", filt)
        return path or ""

    @Slot(result=str)
    def browse_fasta(self):
        return self._open("Select the reference genome", "FASTA (*.fa *.fasta *.fna *.fa.gz *.fasta.gz);;All files (*)")

    @Slot(result=str)
    def browse_bam(self):
        return self._open("Select a BAM file", "BAM (*.bam);;All files (*)")

    @Slot(result=str)
    def browse_vcf(self):
        return self._open("Select a VCF file", "VCF (*.vcf *.vcf.gz *.bcf);;All files (*)")

    @Slot(result=str)
    def browse_annotation(self):
        return self._open("Select a gene annotation", "GFF3 / GTF (*.gff *.gff3 *.gtf *.gff.gz *.gff3.gz *.gtf.gz);;All files (*)")

    @Slot(result=str)
    def browse_folder(self):
        return QFileDialog.getExistingDirectory(self._window, "Select the output folder") or ""

    @Slot(str)
    def open_folder(self, path):
        if path and os.path.isdir(path):
            QDesktopServices.openUrl(QUrl.fromLocalFile(path))

    # ---------------------------------------------------------- information
    @Slot(result=str)
    def get_defaults(self):
        return json.dumps({"defaults": config.DEFAULTS, "output_default": str(RESULT_DIR)})

    @Slot(result=str)
    def get_dependencies(self):
        return json.dumps(deps.check())

    @Slot(result=str)
    def get_validation_experiments(self):
        return json.dumps(validation.reference_summary())

    # ---------------------------------------------------------- threads
    def _start(self, worker, done_signal, done_slot):
        thread = QThread()
        worker.moveToThread(thread)
        thread.started.connect(worker.run, Qt.QueuedConnection)
        done_signal.connect(done_slot, Qt.QueuedConnection)
        done_signal.connect(thread.quit, Qt.QueuedConnection)
        thread.finished.connect(thread.deleteLater, Qt.QueuedConnection)
        thread.finished.connect(lambda w=worker: self._release(w), Qt.QueuedConnection)
        self._threads.append((thread, worker))
        thread.start()

    def _release(self, worker):
        self._threads = [(t, w) for t, w in self._threads if w is not worker]

    # ---------------------------------------------------------- checks
    @Slot(str)
    def check_inputs(self, config_json):
        try:
            cfg = json.loads(config_json)
        except Exception as e:
            self.inputsChecked.emit(json.dumps({"status": "error", "ok": False, "errors": [f"Could not read the settings: {e}"],
                                                 "warnings": [], "files": {}}))
            return
        self.statusUpdate.emit("Checking the input files…")
        worker = CheckWorker(cfg)
        self._start(worker, worker.finished, self._on_checked)

    @Slot(str)
    def _on_checked(self, payload):
        self.inputsChecked.emit(payload)
        self.statusUpdate.emit("Input check finished.")

    # ---------------------------------------------------------- run
    @Slot(str)
    def run_qtlseq(self, config_json):
        try:
            cfg = json.loads(config_json)
        except Exception as e:
            self._fail(f"Could not read the settings: {e}")
            return
        self._holder["cancel"] = threading.Event()
        self._holder["res"] = None
        self.statusUpdate.emit("Starting…")
        worker = RunWorker(cfg, self._holder)
        worker.progress.connect(self._on_progress, Qt.QueuedConnection)
        self._start(worker, worker.finished, self._on_finished)

    @Slot()
    def cancel(self):
        self._holder["cancel"].set()
        self.statusUpdate.emit("Cancelling… the current step is finishing.")

    @Slot(str)
    def _on_progress(self, payload):
        self.progressUpdate.emit(payload)
        try:
            self.statusUpdate.emit(json.loads(payload).get("message", ""))
        except Exception:
            pass

    @Slot(dict)
    def _on_finished(self, result):
        status = result.get("status")
        if status == "success":
            self.last_result = result
            self.statusUpdate.emit("QTL-seq analysis finished.")
        elif status == "cancelled":
            self.statusUpdate.emit("Cancelled.")
        else:
            self.statusUpdate.emit("QTL-seq analysis failed.")
        self.resultReady.emit(json.dumps(result))

    def _fail(self, message):
        self.statusUpdate.emit(message)
        self.resultReady.emit(json.dumps({"status": "error", "message": message}))

    # ---------------------------------------------------------- one chromosome
    @Slot(str, result=str)
    def get_chromosome(self, chrom):
        res = self._holder.get("res")
        if res is None:
            return json.dumps({"status": "error", "message": "Run the analysis first."})
        return json.dumps(pipeline.chromosome_payload(res, chrom))

    # ---------------------------------------------------------- validation
    @Slot(str, int)
    def run_validation(self, key, seed):
        self.statusUpdate.emit("Running the validation on simulated data…")
        worker = ValidationWorker(key, int(seed), self._holder)
        self._start(worker, worker.finished, self._on_validated)

    @Slot(str)
    def _on_validated(self, payload):
        self.validationReady.emit(payload)
        self.statusUpdate.emit("Validation finished.")

    # ---------------------------------------------------------- downloads
    @Slot(str, str)
    def save_file(self, source_path, suggested_name):
        """Native Save As + copy (an <a download> link is unreliable for file:// pages)."""
        try:
            if not source_path or not os.path.exists(source_path):
                self.statusUpdate.emit("That file is not available.")
                return
            dest, _ = QFileDialog.getSaveFileName(self._window, "Save", suggested_name or os.path.basename(source_path))
            if not dest:
                return
            shutil.copy2(source_path, dest)
            self.statusUpdate.emit(f"Saved to {dest}")
        except Exception as e:
            self.statusUpdate.emit(f"Save failed: {e}")

    @Slot(str, str)
    def download_result(self, key, fmt):
        """key: a name from result['files']; fmt 'png' / 'pdf' for figures (both exist for every figure)."""
        if not self.last_result:
            self.statusUpdate.emit("Run the analysis first.")
            return
        files = self.last_result.get("files", {})
        src = files.get(key)
        if isinstance(src, dict):
            return
        if src and fmt in ("pdf", "png") and src.lower().endswith(".png"):
            src = os.path.splitext(src)[0] + "." + fmt
        self.save_file(src or "", os.path.basename(src or ""))

    @Slot(str, str)
    def download_chromosome_plot(self, chrom, fmt):
        if not self.last_result:
            return
        src = (self.last_result.get("files", {}).get("chrom_plots") or {}).get(chrom)
        if src:
            src = os.path.splitext(src)[0] + "." + (fmt if fmt in ("png", "pdf") else "png")
            self.save_file(src, os.path.basename(src))

    @Slot(str)
    def download_all(self, _unused=""):
        """Copy the whole result folder to a place of the user's choice."""
        if not self.last_result:
            self.statusUpdate.emit("Run the analysis first.")
            return
        src = self.last_result["out_dir"]
        dest = QFileDialog.getExistingDirectory(self._window, "Copy the result folder to…")
        if not dest:
            return
        try:
            target = os.path.join(dest, os.path.basename(src))
            shutil.copytree(src, target, dirs_exist_ok=True)
            self.statusUpdate.emit(f"Copied to {target}")
        except Exception as e:
            self.statusUpdate.emit(f"Copy failed: {e}")


# ============================================================
# STANDALONE WINDOW (ToolManager.open_qtl_seq, kept for testing)
# ============================================================

class QTLSeqWindow(QMainWindow):

    def __init__(self, parent=None):
        super().__init__(parent)

        from PySide6.QtWebChannel import QWebChannel
        from PySide6.QtWebEngineWidgets import QWebEngineView
        from ui.frontend import load_frontend

        self.setWindowTitle("TerraTools — QTL-seq")
        self.resize(1400, 900)
        self.browser = QWebEngineView()
        self.setCentralWidget(self.browser)
        self.channel = QWebChannel(self.browser.page())
        self.bridge = QTLSeqBridge(self)
        self.channel.registerObject("qtlSeqBridge", self.bridge)
        self.browser.page().setWebChannel(self.channel)
        load_frontend(self.browser, "qtl_seq")


if __name__ == "__main__":
    from PySide6.QtWidgets import QApplication

    app = QApplication(sys.argv)
    w = QTLSeqWindow()
    w.show()
    sys.exit(app.exec())
