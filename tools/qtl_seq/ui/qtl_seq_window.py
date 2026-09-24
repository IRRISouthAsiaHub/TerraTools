import json
from pathlib import Path

from PySide6.QtCore import QObject, Qt, QThread, QUrl, Signal, Slot
import PySide6.QtWebChannel
from PySide6.QtWebEngineCore import QWebEngineSettings
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QFileDialog, QMainWindow

from tools.qtl_seq.backend import qtlseq_runner


# ============================================================
# PATHS
# ============================================================

QTLSEQ_DIR = Path(__file__).resolve().parent
QTLSEQ_HTML = QTLSEQ_DIR / "qtl_seq.html"

RESULT_DIR = QTLSEQ_DIR / "results"
RESULT_DIR.mkdir(exist_ok=True)


# ============================================================
# WORKER
# ============================================================

class QTLSeqWorker(QObject):

    finished = Signal(dict)
    error = Signal(str)

    def __init__(self, config, output_dir):

        super().__init__()

        self.config = config
        self.output_dir = output_dir

    @Slot()
    def run(self):

        print(
            "[QTL-SEQ WORKER] run() called on worker thread",
            flush=True
        )

        try:

            result = qtlseq_runner.run_qtlseq(
                self.config,
                results_dir=self.output_dir
            )

            print(
                "[QTL-SEQ WORKER] run_qtlseq() returned, status:",
                result.get("status"),
                flush=True
            )

            self.finished.emit(result)

            print(
                "[QTL-SEQ WORKER] finished signal emitted",
                flush=True
            )

        except Exception as e:

            print(
                "[QTL-SEQ WORKER] error:", e, flush=True
            )

            self.error.emit(str(e))


# ============================================================
# BRIDGE
# ============================================================

class QTLSeqBridge(QObject):

    resultReady = Signal(str)
    statusUpdate = Signal(str)

    def __init__(self, parent_window=None):

        super().__init__()

        self._parent_window = parent_window
        self._thread = None
        self._worker = None
        self.last_result = None

    # --------------------------------------------------------
    # FILE BROWSING
    # --------------------------------------------------------

    @Slot(result=str)
    def browse_ref_file(self):

        path, _ = QFileDialog.getOpenFileName(
            self._parent_window,
            "Select Reference FASTA",
            "",
            "FASTA Files (*.fa *.fasta *.fna);;All Files (*)"
        )

        return path or ""

    @Slot(result=str)
    def browse_reads_file(self):

        path, _ = QFileDialog.getOpenFileName(
            self._parent_window,
            "Select FASTQ / BAM File",
            "",
            "Reads (*.fastq *.fq *.fastq.gz *.fq.gz *.bam);;"
            "All Files (*)"
        )

        return path or ""

    @Slot(str, str)
    def save_file(self, source_path: str, suggested_name: str):
        """Native Save As dialog + copy, since <a download> is
        unreliable for local file:// URLs inside QWebEngineView."""

        try:

            dest_path, _ = QFileDialog.getSaveFileName(
                self._parent_window,
                "Save File",
                suggested_name,
                "All Files (*)"
            )

            if not dest_path:
                return

            import shutil
            shutil.copy2(source_path, dest_path)

            print(
                "[QTL-SEQ BRIDGE] saved",
                source_path, "to", dest_path,
                flush=True
            )

            self.statusUpdate.emit(f"Saved to {dest_path}")

        except Exception as e:

            import traceback
            traceback.print_exc()

            self.statusUpdate.emit(f"Save failed: {e}")

    # --------------------------------------------------------
    # RUN
    # --------------------------------------------------------

    @Slot(str)
    def run_qtlseq(self, config_json: str):

        print(
            "[QTL-SEQ BRIDGE] run_qtlseq() called",
            flush=True
        )

        try:
            config = json.loads(config_json)
        except Exception as e:

            self.resultReady.emit(json.dumps({
                "status": "error",
                "message": f"Invalid configuration: {e}"
            }))

            return

        for key in ("ref_file", "parent_file_r1", "bulk1_file_r1", "bulk2_file_r1"):

            if not config.get(key):

                self.resultReady.emit(json.dumps({
                    "status": "error",
                    "message":
                        "Please select all required input files "
                        "(reference, and at least R1 reads for "
                        "parent, bulk1, and bulk2)."
                }))

                return

        print("\n===== QTL-SEQ CONFIGURATION =====")
        print(config)
        print("==================================\n")

        self.statusUpdate.emit("Running QTL-seq analysis...")

        self._thread = QThread()
        self._worker = QTLSeqWorker(config, RESULT_DIR)
        self._worker.moveToThread(self._thread)

        self._thread.started.connect(
            self._worker.run, Qt.QueuedConnection
        )

        self._worker.finished.connect(
            self._on_finished, Qt.QueuedConnection
        )

        self._worker.error.connect(
            self._on_error, Qt.QueuedConnection
        )

        self._worker.finished.connect(
            self._thread.quit, Qt.QueuedConnection
        )

        self._worker.finished.connect(
            self._worker.deleteLater, Qt.QueuedConnection
        )

        self._thread.finished.connect(
            self._thread.deleteLater, Qt.QueuedConnection
        )

        print(
            "[QTL-SEQ BRIDGE] about to call self._thread.start()",
            flush=True
        )

        self._thread.start()

        print(
            "[QTL-SEQ BRIDGE] self._thread.start() returned",
            flush=True
        )

    @Slot(dict)
    def _on_finished(self, result: dict):

        print(
            "[QTL-SEQ BRIDGE] _on_finished() called, status:",
            result.get("status"),
            flush=True
        )

        self.last_result = result

        if result.get("status") == "success":
            self.statusUpdate.emit(
                "QTL-seq analysis completed successfully."
            )
        else:
            self.statusUpdate.emit(
                f"Error: {result.get('message', 'Unknown error')}"
            )

        print(
            "[QTL-SEQ BRIDGE] about to emit resultReady...",
            flush=True
        )

        self.resultReady.emit(json.dumps(result))

        print(
            "[QTL-SEQ BRIDGE] resultReady emitted",
            flush=True
        )

    @Slot(str)
    def _on_error(self, message: str):

        print(
            "[QTL-SEQ BRIDGE] _on_error() called:",
            message,
            flush=True
        )

        result = {"status": "error", "message": message}

        self.last_result = result

        self.statusUpdate.emit("QTL-seq analysis failed.")

        self.resultReady.emit(json.dumps(result))

    # --------------------------------------------------------
    # DOWNLOAD
    # --------------------------------------------------------

    @Slot(str)
    def download_result(self, result_type: str):

        print(
            "[QTL-SEQ BRIDGE] download_result() called for:",
            result_type,
            flush=True
        )

        if not self.last_result:

            self.statusUpdate.emit(
                "Run the analysis first."
            )

            return

        key_map = {
            "snp": "snp_csv",
            "window": "window_csv",
            "plot": "plot",
        }

        key = key_map.get(result_type)

        source = self.last_result.get(key) if key else None

        if not source or not Path(source).exists():

            self.statusUpdate.emit(
                "Result file not found."
            )

            return

        suggested = Path(source).name

        try:

            dest_path, _ = QFileDialog.getSaveFileName(
                self._parent_window,
                "Save File",
                suggested,
                "All Files (*)"
            )

            if not dest_path:
                return

            import shutil
            shutil.copy2(source, dest_path)

            print(
                "[QTL-SEQ BRIDGE] saved", result_type,
                "to", dest_path,
                flush=True
            )

            self.statusUpdate.emit("Result saved.")

        except Exception as e:

            import traceback
            traceback.print_exc()

            self.statusUpdate.emit(f"Save failed: {e}")


# ============================================================
# WINDOW
# ============================================================

class QTLSeqWindow(QMainWindow):

    def __init__(self, parent=None):

        super().__init__(parent)

        self.setWindowTitle("TerraTools — QTL-seq (BSA)")
        self.resize(1400, 900)

        self.browser = QWebEngineView()
        self.setCentralWidget(self.browser)

        # Local file:// results (plots, CSVs) live outside the
        # html's own directory, so local content needs
        # permission to load them.
        settings = self.browser.settings()
        settings.setAttribute(
            QWebEngineSettings.LocalContentCanAccessFileUrls, True
        )
        settings.setAttribute(
            QWebEngineSettings.LocalContentCanAccessRemoteUrls, False
        )

        self.bridge = QTLSeqBridge(parent_window=self)
        self.channel = PySide6.QtWebChannel.QWebChannel()
        self.channel.registerObject("qtlSeqBridge", self.bridge)
        self.browser.page().setWebChannel(self.channel)

        print(
            "[QTL-SEQ WINDOW] loading:", QTLSEQ_HTML, flush=True
        )

        self.browser.load(
            QUrl.fromLocalFile(str(QTLSEQ_HTML))
        )

        self.browser.loadFinished.connect(self._page_loaded)

    def _page_loaded(self, success):

        if success:

            print(
                "[QTL-SEQ WINDOW] HTML loaded successfully",
                flush=True
            )

        else:

            print(
                "[QTL-SEQ WINDOW] ERROR: HTML failed to load",
                flush=True
            )