import os
import sys
import json
import tempfile
from pathlib import Path

from PySide6.QtCore import QObject, Qt, Signal, Slot, QThread, QUrl
from tools.common.single_run import single_run
from PySide6.QtWidgets import QFileDialog
from PySide6.QtWebChannel import QWebChannel
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QMainWindow

import json
import os
import subprocess
import shutil

from tools.common.save_dialog import save_dialog
from tools.common.r_locator import find_rscript, get_r_env, rscript_not_found_message
from tools.common.subprocess_utils import no_window_kwargs


# ============================================================
# PATHS
# ============================================================
from pathlib import Path
from ui.frontend import load_frontend

# tools/linkage/ui
UI_DIR = Path(__file__).resolve().parent

# tools/linkage
LINKAGE_DIR = UI_DIR.parent

# tools/linkage/ui/linkage.html

# tools/linkage/backend
BACKEND_DIR = LINKAGE_DIR / "backend"

# tools/linkage/backend/linkage_runner.R
LINKAGE_RUNNER = BACKEND_DIR / "linkage_runner.R"

# Results go into a short, OS-managed temp folder rather than nested
# inside the app's own installed location. A path like
# C:\Program Files\TerraTools\_internal\tools\linkage\results\... can
# combine with a long Windows username or a redirected/networked profile
# path (common on institutional, AD-joined machines) to exceed Windows'
# 260-character path limit, producing WinError 206. Every other tool in
# this app already uses tempfile.mkdtemp() for exactly this reason —
# Linkage was the one exception.
RESULT_DIR = Path(tempfile.mkdtemp(prefix="terratools_linkage_"))


def _is_chr0(value):
    try:
        return float(value) == 0
    except (TypeError, ValueError):
        return False


def run_r_backend(config):
    if not LINKAGE_RUNNER.exists():
        raise FileNotFoundError(
            "Linkage R backend not found:\n"
            f"{LINKAGE_RUNNER}"
        )

    rscript = find_rscript()
    if rscript is None:
        raise RuntimeError(rscript_not_found_message())

    # The configuration goes to R as a file, never on the command line:
    # Windows caps a command line at 32,767 characters, and configs that
    # carry real data exceed it (WinError 206 "filename or extension is
    # too long" / R's "command line too long").
    fd, config_path = tempfile.mkstemp(prefix="linkage_cfg_", suffix=".json", dir=str(RESULT_DIR))
    with os.fdopen(fd, "w", encoding="utf-8") as f:
        json.dump(config, f, ensure_ascii=False)

    try:
        process = subprocess.run(
            [rscript, str(LINKAGE_RUNNER), config_path],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            cwd=str(BACKEND_DIR),
            env=get_r_env(rscript),
            **no_window_kwargs(),
        )
    finally:
        try:
            os.remove(config_path)
        except OSError:
            pass

    stdout = process.stdout.strip()
    stderr = process.stderr.strip()

    if stderr:
        print("[LINKAGE R] stderr:\n" + stderr, flush=True)

    if process.returncode != 0:
        raise RuntimeError(
            f"R linkage backend failed (exit {process.returncode}).\n\n"
            f"{stderr or stdout or 'No R output.'}"
        )

    if not stdout:
        raise RuntimeError("R linkage backend returned no output.")

    try:
        result = json.loads(stdout)
    except json.JSONDecodeError:
        start = stdout.find("{")
        end = stdout.rfind("}")
        if start < 0 or end <= start:
            raise RuntimeError(
                "Could not parse R backend JSON output:\n" + stdout
            )
        result = json.loads(stdout[start:end + 1])

    # R reports its own failures as {"ok": false, "error": "..."}; those
    # used to fall through and be treated as a successful result.
    if isinstance(result, dict) and result.get("ok") is False:
        raise RuntimeError(result.get("error") or "R linkage backend reported an error.")
    return result


# ============================================================
# PLOT RE-RENDER WORKER ("Update" buttons)
# ============================================================

class _RenderWorker(QObject):

    finished = Signal(str, dict)

    def __init__(self, operation, payload):
        super().__init__()
        self.operation = operation
        self.payload = payload

    @Slot()
    def run(self):
        try:
            result = run_r_backend(self.payload)
            result.setdefault("status", "success")
        except Exception as e:
            import traceback
            traceback.print_exc()
            result = {"status": "error", "message": str(e)}
        self.finished.emit(self.operation, result)


# ============================================================
# BRIDGE
# ============================================================

class LinkageBridge(QObject):

    resultReady = Signal(str)
    statusUpdate = Signal(str)
    errorOccurred = Signal(str)
    cytoDataReady = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)

        self.last_result = None
        self.last_config = None
        self.last_plot_settings = None
        self.uploaded_path = None
        self._render_jobs = {}

    @Slot(str, str)
    def save_file(self, source_path, suggested_name):
        """Native Save As + copy."""
        try:
            dest, _ = save_dialog(None, suggested_name)
            if not dest:
                return
            shutil.copy2(source_path, dest)
            self.statusUpdate.emit(f"Saved to {dest}")
        except Exception as e:
            import traceback
            traceback.print_exc()
            self.statusUpdate.emit(f"Save failed: {e}")

    # --------------------------------------------------------
    # FILE BROWSING
    # --------------------------------------------------------

    @Slot(result=str)
    def browse_genotype(self):

        path, _ = QFileDialog.getOpenFileName(
            None,
            "Select Genotype File",
            "",
            "Genotype Files (*.xlsx *.xls *.csv *.tsv);;All Files (*)"
        )

        return path or ""

    # --------------------------------------------------------
    # UPLOAD / PREVIEW
    # --------------------------------------------------------

    @Slot(str, result=str)
    def get_upload_preview(self, path):

        print("[LINKAGE BRIDGE] get_upload_preview():", path, flush=True)

        try:
            import pandas as pd

            if not path:
                raise ValueError("No genotype file selected.")

            path = os.path.abspath(path)

            if not os.path.exists(path):
                raise FileNotFoundError(path)

            from tools.common.genotype_layout import is_qtl_layout, to_marker_rows
            if is_qtl_layout(path):
                # the QTL Mapping layout (one row per individual): use it as one row per marker
                path = to_marker_rows(path, tempfile.mkdtemp(prefix="terratools_linkage_upload_"))

            if path.lower().endswith((".csv", ".tsv", ".txt")):
                # This runs on the GUI thread (the page waits for the answer), so a text file is not parsed whole:
                # the first rows give the preview and the column list, and only Marker + Chromosome are read in full.
                sep = "	" if path.lower().endswith((".tsv", ".txt")) else ","
                head_df = pd.read_csv(path, sep=sep, nrows=10)
                if "Marker" not in head_df.columns:
                    raise ValueError("Column 'Marker' not found.")
                if "Chromosome" not in head_df.columns:
                    raise ValueError("Column 'Chromosome' not found.")
                df = pd.read_csv(path, sep=sep, usecols=["Marker", "Chromosome"])
                all_columns = list(head_df.columns)
            else:
                df = pd.read_excel(path, sheet_name=0)
                if "Marker" not in df.columns:
                    raise ValueError("Column 'Marker' not found.")
                if "Chromosome" not in df.columns:
                    raise ValueError("Column 'Chromosome' not found.")
                head_df = df.head(10)
                all_columns = list(df.columns)

            self.uploaded_path = path

            chromosomes = [
                str(x) for x in df["Chromosome"].dropna().unique()
            ]

            # Per-chromosome marker counts — the frontend's chrSummaryTable
            # expects this and never received it, since this key didn't
            # exist before.
            chr_counts = (
                df["Chromosome"].dropna().astype(str).value_counts().sort_index()
            )
            chr_summary = [
                {"Chromosome": f"{chrom} (unplaced, excluded)" if _is_chr0(chrom) else chrom,
                 "N_Markers": int(count)}
                for chrom, count in chr_counts.items()
            ]

            n_individuals = max(int(len(all_columns) - 2), 0)

            # Chromosome 0 = unplaced markers; the analysis leaves them out.
            chr0 = pd.to_numeric(df["Chromosome"], errors="coerce") == 0
            n_chr0 = int(chr0.sum())
            real_chromosomes = [c for c in chromosomes if not _is_chr0(c)]
            log = (
                f"Loaded {len(df)} markers across {len(real_chromosomes)} "
                f"chromosomes ({n_individuals} individuals)."
            )
            if n_chr0:
                log += (f"\n{n_chr0} marker(s) on chromosome 0 (unplaced) will be "
                        "left out of the map and listed in the QC report.")

            return json.dumps({
                "status": "success",
                "filename": os.path.basename(path),
                "path": path,
                "markers": int(len(df)),
                "individuals": n_individuals,
                "chromosomes": len(real_chromosomes),
                "chromosome_values": real_chromosomes,
                "unplaced_markers": n_chr0,
                "log": log,
                # Kept both naming conventions so any other frontend code
                # relying on the original "preview"/"columns" keys (if any)
                # keeps working, while matching what linkage.html actually
                # reads (preview_rows/preview_columns/chr_summary).
                "preview_rows": head_df.fillna("").to_dict(orient="records"),
                "preview_columns": [str(x) for x in all_columns],
                "chr_summary": chr_summary,
                "preview": head_df.fillna("").to_dict(orient="records"),
                "columns": [str(x) for x in all_columns],
            }, ensure_ascii=False)

        except Exception as e:

            import traceback
            traceback.print_exc()

            self.errorOccurred.emit(
                f"Could not read genotype file: {e}"
            )

            return json.dumps({
                "status": "error",
                "message": str(e)
            })

    # --------------------------------------------------------
    # RUN LINKAGE ANALYSIS
    # --------------------------------------------------------

    @Slot(str)
    @single_run
    def run_linkage(self, config_json):

        print("[LINKAGE BRIDGE] run_linkage() called", flush=True)

        try:

            config = json.loads(config_json)
            config["genotype_file"] = self.uploaded_path

            self.last_config = config

            print("\n========== LINKAGE CONFIGURATION ==========")
            print(config)
            print("=============================================\n")

            if not self.uploaded_path:
                self._analysis_error(
                    "Upload a genotype file first."
                )
                return

            self.statusUpdate.emit(
                "Building linkage map..."
            )

            thread = QThread()

            worker = LinkageWorker(
                config,
                RESULT_DIR
            )

            worker.moveToThread(thread)

            thread.started.connect(
                worker.run, Qt.QueuedConnection
            )

            worker.finished.connect(
                self._analysis_finished, Qt.QueuedConnection
            )

            worker.error.connect(
                self._analysis_error, Qt.QueuedConnection
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

            print("[LINKAGE BRIDGE] about to call self._thread.start()", flush=True)
            thread.start()
            print("[LINKAGE BRIDGE] self._thread.start() returned", flush=True)

        except Exception as e:

            self._analysis_error(
                str(e)
            )

    @Slot(dict)
    def _analysis_finished(self, result):

        print("[LINKAGE BRIDGE] _analysis_finished() called, status:", result.get("status"), flush=True)

        self.last_result = result

        self.statusUpdate.emit(
            "Linkage map analysis completed successfully."
            if result.get("status") == "success"
            else "Linkage map analysis failed."
        )

        print("[LINKAGE BRIDGE] about to emit resultReady...", flush=True)
        self.resultReady.emit(
            json.dumps(result)
        )
        print("[LINKAGE BRIDGE] resultReady emitted", flush=True)

    @Slot(str)
    def _analysis_error(self, message):

        print("[LINKAGE BRIDGE] _analysis_error() called:", message, flush=True)

        result = {
            "status": "error",
            "message": message
        }

        self.last_result = result

        self.statusUpdate.emit(
            "Linkage map analysis failed."
        )

        self.resultReady.emit(
            json.dumps(result)
        )

    # --------------------------------------------------------
    # MAP / DENSITY / HEATMAP PLOT RE-RENDER
    # --------------------------------------------------------
    # These re-render a figure from the *already computed*
    # linkage map using new display settings (theme, palette,
    # chromosome selection, sizing, etc.) without re-running
    # the analysis itself.

    # Each "Update" button redraws one plot. R reloads the saved analysis
    # (result_rds) instead of receiving the whole result again, and the
    # work runs on a background thread so the window never freezes.

    _RENDER_TARGETS = {
        "render_map_plot": ("showMapPlot", "map_plot"),
        "render_density_plot": ("showDensityPlot", "density_plot"),
        "render_heatmap_plot": ("showHeatmapPlot", "heatmap_plot"),
    }

    @Slot(str)
    def render_map_plot(self, config_json):
        self._start_render("render_map_plot", config_json)

    @Slot(str)
    def render_density_plot(self, config_json):
        self._start_render("render_density_plot", config_json)

    @Slot(str)
    def render_heatmap_plot(self, config_json):
        self._start_render("render_heatmap_plot", config_json)

    def _start_render(self, operation, config_json):
        print(f"[LINKAGE BRIDGE] {operation}() called", flush=True)
        js_fn, _key = self._RENDER_TARGETS[operation]

        if not self.last_result or self.last_result.get("status") != "success":
            self._push_plot(js_fn, {"status": "error", "message": "Run the analysis first."})
            return
        if not self.last_result.get("result_rds"):
            self._push_plot(js_fn, {"status": "error", "message": "Please run the analysis again (results are from an older version)."})
            return

        try:
            payload = dict(json.loads(config_json))
        except Exception as e:
            self._push_plot(js_fn, {"status": "error", "message": f"Invalid plot settings: {e}"})
            return
        payload.update(
            operation=operation,
            result_rds=self.last_result["result_rds"],
            output_dir=str(RESULT_DIR.resolve()),
        )

        thread = QThread()
        worker = _RenderWorker(operation, payload)
        worker.moveToThread(thread)
        thread.started.connect(worker.run, Qt.QueuedConnection)
        worker.finished.connect(self._render_finished, Qt.QueuedConnection)
        worker.finished.connect(thread.quit, Qt.QueuedConnection)
        worker.finished.connect(worker.deleteLater, Qt.QueuedConnection)
        thread.finished.connect(thread.deleteLater, Qt.QueuedConnection)
        # one slot per plot type, so updating two plots at once is fine
        self._render_jobs[operation] = (thread, worker)
        thread.start()

    @Slot(str, dict)
    def _render_finished(self, operation, result):
        js_fn, key = self._RENDER_TARGETS[operation]
        if result.get("status") == "success" and result.get(key) and self.last_result:
            self.last_result[key] = result[key]   # Download buttons save the latest image
        self._push_plot(js_fn, result)

    def _push_plot(self, js_fn, result):
        self.parent().browser.page().runJavaScript(
            f"window.{js_fn} && window.{js_fn}({json.dumps(result)});"
        )

    # --------------------------------------------------------
    # CYTOGENETIC IDEOGRAM DATA
    # --------------------------------------------------------
    # Only the underlying data needs to come from Python;
    # zoom / color scheme / sizing are handled client-side in
    # the canvas renderer without another round trip.

    @Slot(str)
    def get_cyto_data(self, config_json):

        print("[LINKAGE BRIDGE] get_cyto_data() called", flush=True)

        if not self.last_result or self.last_result.get("status") != "success":
            self.statusUpdate.emit("Run analysis first.")
            return

        try:

            config = json.loads(config_json)

            cyto = self.last_result.get("cyto_data", {})

            # The R side returns cyto_data as a flat object keyed
            # directly by chromosome name ({"1": {...}, "2": {...}}).
            # The frontend's renderCytogeneticMap() checks for
            # cytoData.chromosomes specifically and returns immediately
            # if that key is missing — which it always was, since
            # nothing wrapped it. This is the direct cause of the
            # ideogram never rendering. Wrapping it here, at the single
            # point already converting this data for the frontend, is
            # the minimal fix.
            self.cytoDataReady.emit(
                json.dumps({"chromosomes": cyto})
            )

        except Exception as e:

            self.cytoDataReady.emit(
                json.dumps({
                    "status": "error",
                    "message": str(e)
                })
            )

    # --------------------------------------------------------
    # PLOT SETTINGS
    # --------------------------------------------------------

    @Slot(str)
    def save_plot_settings(self, config_json):

        try:

            config = json.loads(
                config_json
            )

            self.last_plot_settings = config

            self.statusUpdate.emit(
                "Plot settings saved."
            )

        except Exception as e:

            self.statusUpdate.emit(
                f"Could not save settings: {e}"
            )

    # --------------------------------------------------------
    # DOWNLOADS
    # --------------------------------------------------------

    @Slot(str)
    def download_result(self, result_type):

        print("[LINKAGE BRIDGE] download_result() called for:", result_type, flush=True)

        if not self.last_result or self.last_result.get("status") != "success":

            self.statusUpdate.emit(
                "Run linkage analysis first."
            )

            return

        mapping = {

            "full_results":
                self.last_result.get("result_xlsx")
                or self.last_result.get("full_results_xlsx"),

            "all_markers":
                self.last_result.get("all_markers_xlsx"),

            "retained_markers":
                self.last_result.get("retained_markers_xlsx"),

            "removal_log":
                self.last_result.get("removal_log_xlsx"),

            "recombination":
                self.last_result.get("recombination_xlsx"),

            "qtl_input_csv":
                self.last_result.get("qtl_input_csv"),

            "qtl_input_xlsx":
                self.last_result.get("qtl_input_xlsx"),

            "map_plot":
                self.last_result.get("map_plot"),

            "density_plot":
                self.last_result.get("density_plot"),

            "heatmap_plot":
                self.last_result.get("heatmap_plot"),
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

            print("[LINKAGE BRIDGE] saved", result_type, "to", destination, flush=True)

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
# LINKAGE WORKER
# ============================================================

class LinkageWorker(QObject):

    finished = Signal(dict)
    error = Signal(str)

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

        print("[LINKAGE WORKER] run() called on worker thread", flush=True)

        try:

            self.config["input_path"] = self.config.get(
                "input_path",
                self.config.get("genotype_file")
            )
            self.config["genotype_file"] = self.config["input_path"]
            self.config["output_dir"] = str(self.output_dir.resolve())

            print(
                "[LINKAGE WORKER] calling R linkage backend...",
                flush=True
            )

            result = run_r_backend(self.config)

            if result.get("ok") is True:
                result["status"] = "success"

            if "status" not in result:
                result["status"] = "success"

            # File for QTL Mapping: the linkage map matched to the input
            # genotypes and transposed. A failure here must not hide an
            # otherwise good linkage result.
            if result.get("status") == "success" and result.get("result_xlsx"):
                try:
                    from tools.linkage.backend.linkage_qtl_export import export_qtl_input
                    result.update(export_qtl_input(
                        self.config["input_path"], result["result_xlsx"], self.output_dir
                    ))
                except Exception as e:
                    import traceback
                    traceback.print_exc()
                    result["qtl_export_error"] = str(e)

            print(
                "[LINKAGE WORKER] R backend returned:",
                result.get("status"),
                flush=True
            )

            self.finished.emit(
                result
            )
            print("[LINKAGE WORKER] finished signal emitted", flush=True)

        except Exception as e:

            print("[LINKAGE WORKER] error:", e, flush=True)

            self.error.emit(
                str(e)
            )


# ============================================================
# WINDOW
# ============================================================

class LinkageWindow(QMainWindow):

    def __init__(self, parent=None):

        super().__init__(
            parent
        )

        self.setWindowTitle(
            "TerraTools — Linkage Map Creator"
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

        self.bridge = LinkageBridge(
            self
        )

        self.channel.registerObject(
            "linkageBridge",
            self.bridge
        )

        self.browser.page().setWebChannel(
            self.channel
        )

        # ----------------------------------------------------
        # Load Linkage HTML
        # ----------------------------------------------------

        load_frontend(self.browser, "linkage")


# ============================================================
# STANDALONE TEST
# ============================================================

if __name__ == "__main__":

    from PySide6.QtWidgets import QApplication

    app = QApplication(
        sys.argv
    )

    window = LinkageWindow()

    window.show()

    sys.exit(
        app.exec()
    )