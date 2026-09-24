"""
pybsaseq_runner.py — Python side of the PyBSA-Seq (BSA-Seq / G-statistic)
tool. Runs the bundled pybsaseq_script.py, exactly the way the original
Gradio app did (script runs with cwd set to a working directory
containing the input files, and writes outputs into a Results/ subfolder
there).

CRITICAL — how the script actually gets run differs between dev mode and
a frozen (PyInstaller) build, for the exact same reason documented in
genomic_prediction_runner.py: in a frozen build, sys.executable IS
TerraTools.exe itself, not a real python.exe. Passing it a script path
does NOT run that script — it just launches another instance of the
whole app, which never runs pybsaseq_script.py's logic and can exit
without producing any output (this was the actual cause of "process
appears to start, but no results are generated"). In a frozen build we
instead re-invoke the exe with `--worker pybsaseq --cwd <dir> ...`,
which main.py intercepts BEFORE creating a QApplication and executes the
bundled script directly via runpy in that separate process.

Same deadlock-safe / timeout-safe subprocess pattern used for the
R-backed tools: capture_output=True, an explicit timeout, stdin closed
and no console window on Windows, and collecting results by walking the
output directory afterward rather than parsing stdout.
"""

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from tools.pybsaseq.config import PyBSASeqConfig
from tools.common.subprocess_utils import no_window_kwargs

TIMEOUT_SECONDS = 3600  # matches the original Gradio app's 1-hour timeout
VALID_EXTENSIONS = (".csv", ".pdf", ".png", ".eps", ".svg", ".txt", ".html")


class PyBSASeqEngine:
    def __init__(self, config: PyBSASeqConfig):
        self.config = config
        self._script_path = Path(__file__).resolve().parent / "pybsaseq_script.py"

    def _build_flags(self, parents_name: str, bulks_name: str) -> list:
        c = self.config
        flags = ["-i", f"{parents_name},{bulks_name}"]
        if c.bulksizes:
            flags += ["-b", c.bulksizes]
        if c.popstruct:
            flags += ["-p", c.popstruct]
        if c.pvalues:
            flags += ["-v", c.pvalues]
        if c.replication:
            flags += ["-r", str(c.replication)]
        if c.slidingwindow:
            flags += ["-s", c.slidingwindow]
        if c.smoothing:
            flags += ["-m", c.smoothing]
        if c.chromosome_order:
            flags += ["--chromosome_order", c.chromosome_order]
        if c.read_length:
            flags += ["-d", str(c.read_length)]
        if c.region:
            flags += ["-e", c.region]
        if c.misc_params:
            flags += ["-c", c.misc_params]
        if c.adjust_gap:
            flags += ["-a", "true"]
        if c.parent_ref:
            flags += ["--parent", c.parent_ref]
        return flags

    def run(self) -> dict:
        if not self.config.parents_file or not self.config.bulks_file:
            return {"status": "error", "message": "Both a parents file and a bulks file are required."}

        work_dir = Path(tempfile.mkdtemp(prefix="terratools_pybsaseq_work_"))
        outdir = Path(tempfile.mkdtemp(prefix="terratools_pybsaseq_out_"))

        try:
            parents_name = Path(self.config.parents_file).name
            bulks_name = Path(self.config.bulks_file).name
            shutil.copyfile(self.config.parents_file, work_dir / parents_name)
            shutil.copyfile(self.config.bulks_file, work_dir / bulks_name)
        except Exception as e:
            return {"status": "error", "message": f"Failed to stage input files: {e}"}

        flags = self._build_flags(parents_name, bulks_name)

        if getattr(sys, "frozen", False):
            # Frozen build: re-invoke this same exe with a worker flag —
            # main.py intercepts this before touching PySide6/QtWebEngine
            # and runs the bundled script directly via runpy, chdir'ing
            # to work_dir itself before doing so.
            cmd = [sys.executable, "--worker", "pybsaseq", "--cwd", str(work_dir)] + flags
            popen_cwd = None
        else:
            # Dev mode: sys.executable is a real python.exe, safe to
            # point directly at the script file. The script needs to be
            # physically present in work_dir since it's run with that as
            # cwd (matching the original Gradio app's approach).
            try:
                shutil.copyfile(self._script_path, work_dir / self._script_path.name)
            except Exception as e:
                return {"status": "error", "message": f"Failed to stage the analysis script: {e}"}
            cmd = [sys.executable, self._script_path.name] + flags
            popen_cwd = str(work_dir)

        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=popen_cwd,
                timeout=TIMEOUT_SECONDS,
                stdin=subprocess.DEVNULL,
                **no_window_kwargs(),
            )
        except subprocess.TimeoutExpired:
            return {
                "status": "error",
                "message": f"Analysis did not finish within {TIMEOUT_SECONDS}s (timed out). "
                           "Try reducing the number of replications.",
            }
        except Exception as e:
            return {"status": "error", "message": f"Failed to launch the analysis script: {e}"}

        results_dir = work_dir / "Results"
        output_files = []
        html_report = None
        if results_dir.exists():
            for src in sorted(results_dir.rglob("*")):
                if src.is_file() and src.suffix.lower() in VALID_EXTENSIONS:
                    rel = src.relative_to(results_dir)
                    dest_name = "__".join(rel.parts)  # flatten subfolders into a readable unique name
                    dest = outdir / dest_name
                    shutil.copyfile(src, dest)
                    output_files.append({
                        "name": dest_name,
                        "path": str(dest),
                        "ext": src.suffix.lower().lstrip("."),
                    })
                    if src.suffix.lower() == ".html" and html_report is None:
                        html_report = str(dest)

        # Clean up the working directory now that outputs are collected —
        # only the copies in outdir need to persist.
        shutil.rmtree(work_dir, ignore_errors=True)

        if proc.returncode != 0:
            return {
                "status": "error",
                "message": (
                    f"Script exited with code {proc.returncode}.\n\n"
                    f"stderr:\n{proc.stderr[-3000:]}\n\n"
                    f"stdout:\n{proc.stdout[-1500:]}"
                ),
                "output_files": output_files,  # partial results, if any
            }

        return {
            "status": "success",
            "stdout": proc.stdout[-4000:],
            "output_files": output_files,
            "html_report": html_report,
        }
