import os
import sys
import json
import tempfile
import subprocess
from pathlib import Path

from PySide6.QtCore import QObject, Qt, Signal, Slot, QThread, QUrl
from tools.common.single_run import single_run
from PySide6.QtWidgets import QFileDialog, QMainWindow
from PySide6.QtWebChannel import QWebChannel
from PySide6.QtWebEngineWidgets import QWebEngineView

# IMPORTANT: heterosis_runner (which imports TensorFlow) is
# deliberately NOT imported here. TensorFlow's bundled LLVM
# conflicts with QtWebEngine's bundled LLVM if both are loaded in
# the same process ("CommandLine Error: ... registered more than
# once!" -> abort). The heavy computation runs in a separate
# subprocess via heterosis_process.py instead. Only the lightweight,
# TensorFlow-free preview helper is imported directly.
from tools.heterosis.backend import heterosis_preview
from tools.common.subprocess_utils import no_window_kwargs
from ui.frontend import load_frontend
from tools.common import json_safe
from tools.common.paths import results_dir


# ============================================================
# PATHS
# ============================================================

HETEROSIS_DIR = Path(__file__).resolve().parent

BACKEND_DIR = HETEROSIS_DIR.parent / "backend"
PROCESS_SCRIPT = BACKEND_DIR / "heterosis_process.py"

# Project root = parent of the top-level "tools" package, so the
# subprocess can resolve `tools.heterosis.backend.heterosis_runner`.
PROJECT_ROOT = HETEROSIS_DIR.parents[2]

RESULT_DIR = results_dir("heterosis", HETEROSIS_DIR / "results")


# ============================================================
# BRIDGE
# ============================================================

class HeterosisBridge(QObject):

    resultReady = Signal(str)
    statusUpdate = Signal(str)
    progressUpdate = Signal(float, str)
    errorOccurred = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)

        self.last_result = None
        self.last_config = None

        self.file_paths = {
            "train_pheno_file": None,
            "train_geno_file": None,
            "parent_geno_file": None,
            "test_geno_file": None,
        }

    # --------------------------------------------------------
    # FILE BROWSING
    # --------------------------------------------------------

    @Slot(str, result=str)
    def browse_file(self, key):

        path, _ = QFileDialog.getOpenFileName(
            None,
            "Select CSV File",
            "",
            "CSV Files (*.csv);;All Files (*)"
        )

        if not path:
            return ""

        if key in self.file_paths:
            self.file_paths[key] = path

        return path

    @Slot(str, result=str)
    def get_file_preview(self, path):

        print("[HETEROSIS BRIDGE] get_file_preview() called for:", path, flush=True)

        try:

            preview = heterosis_preview.load_preview(path)

            return json.dumps(preview)

        except Exception as e:

            return json.dumps({
                "status": "error",
                "message": str(e)
            })

    # --------------------------------------------------------
    # RUN HETEROSIS PREDICTION
    # --------------------------------------------------------

    @Slot(str)
    @single_run
    def run_heterosis(self, config_json):

        print("[HETEROSIS BRIDGE] run_heterosis() called", flush=True)

        try:

            config = json.loads(config_json)

            config["train_pheno_file"] = self.file_paths.get("train_pheno_file")
            config["train_geno_file"] = self.file_paths.get("train_geno_file")
            config["parent_geno_file"] = self.file_paths.get("parent_geno_file")
            config["test_geno_file"] = self.file_paths.get("test_geno_file")

            self.last_config = config

            print("\n========== HETEROSIS CONFIGURATION ==========")
            print({k: v for k, v in config.items() if not k.endswith("_file")})
            print("===============================================\n")

            missing = [
                k for k in ("train_pheno_file", "train_geno_file", "parent_geno_file")
                if not config.get(k)
            ]

            if missing:
                self._run_error(
                    "Missing required file(s): " + ", ".join(missing)
                )
                return

            self.statusUpdate.emit(
                "Starting cross-validation..."
            )

            thread = QThread()

            worker = HeterosisWorker(
                config,
                RESULT_DIR
            )

            worker.moveToThread(thread)

            thread.started.connect(
                worker.run, Qt.QueuedConnection
            )

            worker.finished.connect(
                self._run_finished, Qt.QueuedConnection
            )

            worker.error.connect(
                self._run_error, Qt.QueuedConnection
            )

            worker.progress.connect(
                self._run_progress, Qt.QueuedConnection
            )

            worker.finished.connect(
                thread.quit, Qt.QueuedConnection
            )

            worker.finished.connect(
                worker.deleteLater, Qt.QueuedConnection
            )

            thread.finished.connect(
                thread.deleteLater, Qt.QueuedConnection
            )

            self._thread = thread
            self._worker = worker

            print("[HETEROSIS BRIDGE] about to call self._thread.start()", flush=True)
            thread.start()
            print("[HETEROSIS BRIDGE] self._thread.start() returned", flush=True)

        except Exception as e:

            self._run_error(
                str(e)
            )

    @Slot(float, str)
    def _run_progress(self, fraction, message):

        self.progressUpdate.emit(
            fraction if fraction is not None else -1.0,
            message
        )

    @Slot(dict)
    def _run_finished(self, result):

        print("[HETEROSIS BRIDGE] _run_finished() called, status:", result.get("status"), flush=True)

        self.last_result = result

        self.statusUpdate.emit(
            "Heterosis prediction completed successfully."
            if result.get("status") == "success"
            else "Heterosis prediction failed."
        )

        self.resultReady.emit(
            json_safe.dumps(result)
        )

    @Slot(str)
    def _run_error(self, message):

        print("[HETEROSIS BRIDGE] _run_error() called:", message, flush=True)

        result = {
            "status": "error",
            "message": message
        }

        self.last_result = result

        self.statusUpdate.emit(
            "Heterosis prediction failed."
        )

        self.resultReady.emit(
            json_safe.dumps(result)
        )

    # --------------------------------------------------------
    # DOWNLOADS
    # --------------------------------------------------------

    @Slot(str)
    def download_result(self, result_type):

        print("[HETEROSIS BRIDGE] download_result() called for:", result_type, flush=True)

        if not self.last_result or self.last_result.get("status") != "success":

            self.statusUpdate.emit(
                "Run heterosis prediction first."
            )

            return

        mapping = {
            "f1_hybrids": self.last_result.get("f1_csv"),
            "top_parents": self.last_result.get("parents_csv"),
            "all_hybrids": self.last_result.get("all_csv"),
            "test_predictions": self.last_result.get("test_csv"),
            "plot": self.last_result.get("plot_path"),
        }

        source = mapping.get(
            result_type
        )

        if not source or not os.path.exists(source):

            self.statusUpdate.emit(
                "Requested result is not available."
            )

            return

        destination, _ = QFileDialog.getSaveFileName(
            None,
            "Save Result",
            os.path.basename(source),
            "All Files (*)"
        )

        if not destination:
            return

        try:

            import shutil

            shutil.copy2(
                source,
                destination
            )

            print("[HETEROSIS BRIDGE] saved", result_type, "to", destination, flush=True)

            self.statusUpdate.emit(
                "Result saved."
            )

        except Exception as e:

            import traceback
            traceback.print_exc()

            self.statusUpdate.emit(
                f"Download failed: {e}"
            )


# ============================================================
# WORKER  (launches an isolated subprocess for the TF-heavy work)
# ============================================================

class HeterosisWorker(QObject):

    finished = Signal(dict)
    error = Signal(str)
    progress = Signal(float, str)

    def __init__(
        self,
        config,
        output_dir
    ):

        super().__init__()

        self.config = config
        self.output_dir = output_dir

    @Slot()
    def run(self):

        print("[HETEROSIS WORKER] run() called on worker thread", flush=True)

        config_path = None

        try:

            os.makedirs(self.output_dir, exist_ok=True)

            # Write the config to a temp file the subprocess reads back.
            with tempfile.NamedTemporaryFile(
                mode="w",
                suffix=".json",
                dir=str(self.output_dir),
                delete=False
            ) as f:
                json.dump(self.config, f)
                config_path = f.name

            if getattr(sys, "frozen", False):
                # Frozen build: sys.executable IS TerraTools.exe itself.
                # Passing it a script path (as in the dev-mode branch
                # below) does NOT run that script — it just launches
                # another instance of the whole app, which can exit
                # quickly without ever running heterosis_process.py's
                # logic. This was the actual cause of "subprocess exited
                # without a result" in the packaged build. Instead,
                # re-invoke the exe with a --worker flag, which main.py
                # intercepts BEFORE creating a QApplication and
                # dispatches straight to heterosis_process.main() in
                # this same (separate OS) process.
                cmd = [
                    sys.executable, "--worker", "heterosis",
                    config_path, str(self.output_dir),
                ]
                cwd = None
                env = None  # the frozen exe already has every bundled module importable
            else:
                # Dev mode: sys.executable is a real python.exe, safe to
                # point directly at the standalone script.
                env = os.environ.copy()
                existing = env.get("PYTHONPATH", "")
                env["PYTHONPATH"] = str(PROJECT_ROOT) + (os.pathsep + existing if existing else "")

                cmd = [
                    sys.executable, "-u",
                    str(PROCESS_SCRIPT),
                    config_path,
                    str(self.output_dir)
                ]
                cwd = str(PROJECT_ROOT)

            print("[HETEROSIS WORKER] launching subprocess:", " ".join(cmd), flush=True)

            # ── TEMPORARY DEBUG — remove once the frozen-mode issue is found ──
            try:
                with open(r"C:\Users\baboo\Desktop\terratools_heterosis_cmd_debug.txt", "w") as _dbg:
                    _dbg.write(f"frozen={getattr(sys, 'frozen', False)}\n")
                    _dbg.write(f"sys.executable={sys.executable!r}\n")
                    _dbg.write(f"cmd={cmd!r}\n")
            except Exception:
                pass
            # ────────────────────────────────────────────────────────────────

            process = subprocess.Popen(
                cmd,
                cwd=cwd,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,  # merged into stdout — reading only
                                            # stdout while a separate stderr pipe
                                            # fills up (TensorFlow is very chatty
                                            # on stderr at import time) deadlocks
                                            # both the child (blocked writing) and
                                            # this loop (blocked reading)
                stdin=subprocess.DEVNULL,  # a frozen console=False app has no
                                            # valid stdin to inherit — without this
                                            # the child can hang waiting on input
                                            # it will never receive
                text=True,
                bufsize=1,
                **no_window_kwargs(),
            )

            result = None

            for line in iter(process.stdout.readline, ""):

                line = line.rstrip("\n")

                if not line:
                    continue

                if line.startswith("PROGRESS "):

                    rest = line[len("PROGRESS "):]
                    parts = rest.split(" ", 1)

                    try:
                        fraction = float(parts[0])
                    except ValueError:
                        fraction = -1.0

                    message = parts[1] if len(parts) > 1 else ""

                    self.progress.emit(fraction, message)

                elif line.startswith("RESULT "):

                    payload = line[len("RESULT "):]

                    try:
                        result = json.loads(payload)
                    except Exception as e:
                        result = {
                            "status": "error",
                            "message": f"Could not parse subprocess result: {e}"
                        }

                else:

                    print("[HETEROSIS SUBPROCESS]", line, flush=True)

            process.stdout.close()
            process.wait()

            if result is None:

                message = (
                    "Heterosis subprocess exited without a result "
                    f"(exit code {process.returncode})."
                )

                self.error.emit(message)
                return

            print("[HETEROSIS WORKER] subprocess finished, status:", result.get("status"), flush=True)

            self.finished.emit(result)

        except Exception as e:

            print("[HETEROSIS WORKER] error:", e, flush=True)

            import traceback
            traceback.print_exc()

            self.error.emit(
                str(e)
            )

        finally:

            if config_path and os.path.exists(config_path):
                try:
                    os.remove(config_path)
                except OSError:
                    pass


# ============================================================
# WINDOW
# ============================================================

class HeterosisWindow(QMainWindow):

    def __init__(self, parent=None):

        super().__init__(
            parent
        )

        self.setWindowTitle(
            "TerraTools — Heterosis Prediction"
        )

        self.showMaximized()

        self.browser = QWebEngineView()

        self.setCentralWidget(
            self.browser
        )

        # ----------------------------------------------------
        # WebChannel
        # ----------------------------------------------------

        self.channel = QWebChannel(
            self.browser.page()
        )

        self.bridge = HeterosisBridge(
            self
        )

        self.channel.registerObject(
            "heterosisBridge",
            self.bridge
        )

        self.browser.page().setWebChannel(
            self.channel
        )

        # ----------------------------------------------------
        # Load Heterosis HTML
        # ----------------------------------------------------

        load_frontend(self.browser, "heterosis")


# ============================================================
# STANDALONE TEST
# ============================================================

if __name__ == "__main__":

    from PySide6.QtWidgets import QApplication

    app = QApplication(
        sys.argv
    )

    window = HeterosisWindow()

    window.show()

    sys.exit(
        app.exec()
    )