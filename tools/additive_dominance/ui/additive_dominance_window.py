import json
import threading
from pathlib import Path

from PySide6.QtCore import QObject, Qt, QThread, QUrl, Signal, Slot
from PySide6.QtWebChannel import QWebChannel
from PySide6.QtWebEngineCore import QWebEngineSettings
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QFileDialog, QMainWindow

from tools.common.r_locator import find_rscript, get_r_env, rscript_not_found_message
from tools.common.subprocess_utils import no_window_kwargs
from ui.frontend import load_frontend
from tools.common.paths import results_dir


# ============================================================
# PATHS
# ============================================================

AD_DIR = Path(__file__).resolve().parent

BACKEND_DIR = AD_DIR.parent / "backend"
CLI_FILE = BACKEND_DIR / "additive_dominance_cli.R"

RESULT_DIR = results_dir("additive_dominance", AD_DIR / "results")


# ============================================================
# WORKER
# ============================================================

class AdditiveDominanceWorker(QObject):

    finished = Signal(str)      # JSON text, see the note in the bridge
    error = Signal(str)

    def __init__(self, config, results_dir):

        super().__init__()

        self.config = config
        self.results_dir = Path(results_dir)
        self.process = None
        self._process_lock = threading.Lock()

    @Slot()
    def run(self):

        print(
            "[ADDITIVE-DOMINANCE WORKER] run() called on "
            "worker thread",
            flush=True
        )

        import subprocess

        config_file = None

        try:

            self.results_dir.mkdir(parents=True, exist_ok=True)

            config_file = (
                self.results_dir
                / "additive_dominance_config.json"
            )

            config_file.write_text(
                json.dumps(self.config, indent=2),
                encoding="utf-8"
            )

            rscript = find_rscript()
            if rscript is None:
                raise RuntimeError(rscript_not_found_message())

            print(
                "[ADDITIVE-DOMINANCE] Rscript:", rscript,
                flush=True
            )

            command = [
                str(rscript),
                str(CLI_FILE),
                str(config_file),
                str(self.results_dir),
            ]

            print(
                "[ADDITIVE-DOMINANCE] Command:", command,
                flush=True
            )

            proc = subprocess.Popen(
                command,
                cwd=str(BACKEND_DIR),
                env=get_r_env(rscript),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                stdin=subprocess.DEVNULL,
                text=True,
                encoding="utf-8",
                errors="replace",
                bufsize=1,
                **no_window_kwargs(),
            )

            with self._process_lock:
                self.process = proc

            output_lines = []

            if proc.stdout is not None:

                for line in proc.stdout:

                    line = line.rstrip()

                    if not line:
                        continue

                    output_lines.append(line)

                    print(
                        "[ADDITIVE-DOMINANCE R]", line,
                        flush=True
                    )

            return_code = proc.wait()

            with self._process_lock:

                if self.process is proc:
                    self.process = None

            if return_code != 0:

                output = "\n".join(output_lines)

                if output.strip():

                    raise RuntimeError(
                        "Matrix computation failed "
                        f"(exit code {return_code}).\n\n"
                        + output[-10000:]
                    )

                else:

                    raise RuntimeError(
                        "Matrix computation failed "
                        f"(exit code {return_code}) and "
                        "produced no output at all. Try "
                        "running this command manually in a "
                        "terminal:\n\n"
                        + " ".join(f'"{c}"' for c in command)
                    )

            result_json = (
                self.results_dir
                / "additive_dominance_results.json"
            )

            result = None

            if result_json.exists():

                result = json.loads(
                    result_json.read_text(encoding="utf-8")
                )

            else:

                output = "\n".join(output_lines)
                json_start = output.rfind('{"status"')

                if json_start != -1:

                    try:
                        result = json.loads(output[json_start:])
                    except Exception:
                        pass

            if result is None:

                raise FileNotFoundError(
                    "R completed successfully, but no "
                    "result JSON was returned."
                )

            print(
                "[ADDITIVE-DOMINANCE WORKER] parsed result, "
                "status:", result.get("status"),
                flush=True
            )

            self.finished.emit(json.dumps(result))

        except Exception as e:

            import traceback
            traceback.print_exc()

            self.error.emit(str(e))

        finally:

            with self._process_lock:
                self.process = None

            if config_file is not None:

                try:
                    config_file.unlink()
                except OSError:
                    pass

    @Slot()
    def stop(self):

        with self._process_lock:
            proc = self.process

        if proc is not None:

            try:

                if proc.poll() is None:

                    proc.terminate()

                    try:
                        proc.wait(timeout=5)
                    except Exception:
                        proc.kill()

            except Exception:
                pass

            with self._process_lock:

                if self.process is proc:
                    self.process = None


# ============================================================
# BRIDGE
# ============================================================

class AdditiveDominanceBridge(QObject):

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
    def browse_snp_file(self):

        path, _ = QFileDialog.getOpenFileName(
            self._parent_window,
            "Select SNP File",
            "",
            "SNP Files (*.csv *.txt *.tsv);;All Files (*)"
        )

        return path or ""

    @Slot(result=str)
    def browse_train_file(self):

        path, _ = QFileDialog.getOpenFileName(
            self._parent_window,
            "Select the Training file (from Clustering)",
            "",
            "Genotype Files (*.csv *.txt *.tsv);;All Files (*)"
        )

        return path or ""

    @Slot(result=str)
    def browse_test_file(self):

        path, _ = QFileDialog.getOpenFileName(
            self._parent_window,
            "Select the Testing file (from Clustering)",
            "",
            "Genotype Files (*.csv *.txt *.tsv);;All Files (*)"
        )

        return path or ""

    # --------------------------------------------------------
    # RUN
    # --------------------------------------------------------

    @Slot(str)
    def run_matrices(self, config_json: str):

        print(
            "[ADDITIVE-DOMINANCE BRIDGE] run_matrices() called",
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

        if config.get("train_file") or config.get("test_file"):

            if not (config.get("train_file") and config.get("test_file")):

                self.resultReady.emit(json.dumps({
                    "status": "error",
                    "message": "Please select both the training and the testing file."
                }))

                return

        elif not config.get("snp_file"):

            self.resultReady.emit(json.dumps({
                "status": "error",
                "message": "Please select a SNP file."
            }))

            return

        if not config.get("matrix_types"):

            self.resultReady.emit(json.dumps({
                "status": "error",
                "message":
                    "Please select at least one matrix type."
            }))

            return

        print("\n===== ADDITIVE/DOMINANCE CONFIGURATION =====")
        print(config)
        print("==============================================\n")

        self.statusUpdate.emit("Building relationship matrices...")

        self._thread = QThread()
        self._worker = AdditiveDominanceWorker(config, RESULT_DIR)
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
            "[ADDITIVE-DOMINANCE BRIDGE] about to call "
            "self._thread.start()",
            flush=True
        )

        self._thread.start()

    @Slot(str)
    def _on_finished(self, result_json: str):

        result = json.loads(result_json)

        print(
            "[ADDITIVE-DOMINANCE BRIDGE] _on_finished() "
            "called, status:", result.get("status"),
            flush=True
        )

        self.last_result = result

        if result.get("status") == "success":
            self.statusUpdate.emit(
                "Matrices computed successfully."
            )
        else:
            self.statusUpdate.emit(
                f"Error: {result.get('message', 'Unknown error')}"
            )

        self.resultReady.emit(result_json)

    @Slot(str)
    def _on_error(self, message: str):

        print(
            "[ADDITIVE-DOMINANCE BRIDGE] _on_error() called:",
            message, flush=True
        )

        result = {"status": "error", "message": message}

        self.last_result = result

        self.statusUpdate.emit("Matrix computation failed.")

        self.resultReady.emit(json.dumps(result))

    # --------------------------------------------------------
    # STOP (called from window.closeEvent)
    # --------------------------------------------------------

    @Slot()
    def stop_matrices(self):

        if self._worker is not None:

            print(
                "[ADDITIVE-DOMINANCE BRIDGE] stopping "
                "running worker...",
                flush=True
            )

            try:
                self._worker.stop()
            except Exception:
                import traceback
                traceback.print_exc()

        thread = self._thread

        if thread is None:
            return

        try:
            if thread.isRunning():

                print(
                    "[ADDITIVE-DOMINANCE BRIDGE] waiting for "
                    "worker thread to stop...",
                    flush=True
                )

                thread.quit()

                if not thread.wait(5000):

                    print(
                        "[ADDITIVE-DOMINANCE BRIDGE] worker "
                        "thread did not stop in time, "
                        "terminating it.",
                        flush=True
                    )

                    thread.terminate()
                    thread.wait(2000)

        except RuntimeError as e:
            # The underlying C++ QThread has already been deleted.
            print(
                "[ADDITIVE-DOMINANCE BRIDGE] thread already "
                f"deleted: {e}",
                flush=True
            )

        finally:
            self._thread = None
            self._worker = None
    # --------------------------------------------------------
    # DOWNLOAD
    # --------------------------------------------------------

    @Slot(str)
    def download_result(self, result_type: str):

        print(
            "[ADDITIVE-DOMINANCE BRIDGE] download_result() "
            "called for:", result_type,
            flush=True
        )

        if not self.last_result:

            self.statusUpdate.emit(
                "Run the matrix computation first."
            )

            return

        key_map = {
            "additive": "additive_csv",
            "dominance": "dominance_csv",
            "train_additive": "train_additive_csv",
            "test_additive": "test_additive_csv",
            "train_dominance": "train_dominance_csv",
            "test_dominance": "test_dominance_csv",
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
                "CSV Files (*.csv);;All Files (*)"
            )

            if not dest_path:
                return

            import shutil
            shutil.copy2(source, dest_path)

            print(
                "[ADDITIVE-DOMINANCE BRIDGE] saved",
                result_type, "to", dest_path,
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

class AdditiveDominanceWindow(QMainWindow):

    def __init__(self, parent=None):

        super().__init__(parent)

        self.setWindowTitle(
            "TerraTools — Genomic Relationship Matrix Builder"
        )
        self.resize(1400, 900)

        self.browser = QWebEngineView()
        self.setCentralWidget(self.browser)

        settings = self.browser.settings()
        settings.setAttribute(
            QWebEngineSettings.LocalContentCanAccessFileUrls, True
        )
        settings.setAttribute(
            QWebEngineSettings.LocalContentCanAccessRemoteUrls, False
        )

        self.bridge = AdditiveDominanceBridge(parent_window=self)
        self.channel = QWebChannel()
        self.channel.registerObject(
            "additiveDominanceBridge", self.bridge
        )
        self.browser.page().setWebChannel(self.channel)

        print(
            "[ADDITIVE-DOMINANCE WINDOW] loading React frontend",
            flush=True
        )

        load_frontend(self.browser, "additive_dominance")

        self.browser.loadFinished.connect(self._page_loaded)

    def _page_loaded(self, success):

        if success:

            print(
                "[ADDITIVE-DOMINANCE WINDOW] HTML loaded "
                "successfully",
                flush=True
            )

        else:

            print(
                "[ADDITIVE-DOMINANCE WINDOW] ERROR: HTML "
                "failed to load",
                flush=True
            )

    def closeEvent(self, event):

        print(
            "[ADDITIVE-DOMINANCE WINDOW] closeEvent() called",
            flush=True
        )

        try:
            self.bridge.stop_matrices()
        except Exception:
            import traceback
            traceback.print_exc()

        try:
            super().closeEvent(event)
        except Exception:
            import traceback
            traceback.print_exc()

        event.accept()

        print(
            "[ADDITIVE-DOMINANCE WINDOW] closeEvent() finished",
            flush=True
        )
