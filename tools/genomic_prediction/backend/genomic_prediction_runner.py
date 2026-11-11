"""
genomic_prediction_runner.py — Python side of the BreedSight genomic
prediction tool.

IMPORTANT: this runner launches the heavy computation as a SEPARATE
PROCESS, rather than calling the engine in-process like GWAS/Diversity/
Clustering do. TensorFlow cannot be imported into the same process as Qt
WebEngine (see genomic_prediction_cli.py's docstring for why) — doing so
causes a hard, uncatchable LLVM_ERROR abort() that kills the entire
application, not just this tool.

CRITICAL — how the subprocess is launched differs between dev mode and a
frozen (PyInstaller) build:
  - Running from source, sys.executable is a real python.exe, so it can
    be pointed directly at genomic_prediction_cli.py as a script to run.
  - In a frozen build, sys.executable IS TerraTools.exe itself — passing
    it a script path does NOT run that script, it just launches another
    instance of the whole app, which never runs the pipeline and can
    exit quickly without producing a result (this was the actual cause
    of "subprocess exited (code 0) without producing a result").
    Instead, the frozen exe must be re-invoked with a special
    `--worker genomic_prediction ...` flag, which main.py intercepts
    BEFORE creating a QApplication and dispatches straight to
    genomic_prediction_cli.main() in that separate process.

Uses subprocess.Popen with line-by-line stdout reading (not
subprocess.run(capture_output=True)) so log lines can stream to the UI
live as the pipeline runs, and so a multi-minute run's output never has
to sit fully buffered in a pipe.
"""

import os
import subprocess
import sys
import json
import tempfile
from pathlib import Path

from tools.genomic_prediction.config import GenomicPredictionConfig
from tools.common.subprocess_utils import no_window_kwargs

TIMEOUT_SECONDS = 4 * 3600  # 4 hours — this pipeline's runtime varies enormously
                            # with dataset size / epochs / folds; a hard ceiling
                            # exists only so a genuinely broken run fails loudly
                            # instead of hanging forever.


class GenomicPredictionEngine:
    def __init__(self, config: GenomicPredictionConfig, on_log=None):
        self.config = config
        self.on_log = on_log
        self._cli_path = Path(__file__).resolve().parent / "genomic_prediction_cli.py"

    def _build_flags(self, outdir: Path) -> list:
        c = self.config
        return [
            "--train_pheno", c.train_pheno_file,
            "--train_geno", c.train_geno_file,
            "--test_pheno", c.test_pheno_file,
            "--test_geno", c.test_geno_file,
            "--outer_n_splits", str(c.outer_n_splits),
            "--epochs", str(c.epochs),
            "--batch_size", str(c.batch_size),
            "--learning_rate", str(c.learning_rate),
            "--l2_reg", str(c.l2_reg),
            "--dropout_rate", str(c.dropout_rate),
            "--rf_n_estimators", str(c.rf_n_estimators),
            "--rf_max_depth", str(c.rf_max_depth),
            "--alpha", str(c.alpha),
            "--heritability", str(c.heritability),
            "--use_pca", "true" if c.use_pca else "false",
            "--n_components", str(c.n_components),
            "--feature_selection", "true" if c.feature_selection else "false",
            "--rfe_n_features", str(c.rfe_n_features),
            "--use_raw_genotypes", "true" if c.use_raw_genotypes else "false",
            "--baselines", ",".join(c.baselines),
            "--outdir", str(outdir),
        ]

    def run(self) -> dict:
        outdir = Path(tempfile.mkdtemp(prefix="terratools_genomicpred_"))
        flags = self._build_flags(outdir)

        if getattr(sys, "frozen", False):
            # Frozen build: re-invoke this same exe with a worker flag —
            # main.py intercepts this before touching PySide6/QtWebEngine.
            cmd = [sys.executable, "--worker", "genomic_prediction"] + flags
            cwd = None
            env = None
        else:
            # Dev mode: sys.executable is a real python.exe, safe to point
            # directly at the CLI script.
            cmd = [sys.executable, str(self._cli_path)] + flags
            cwd = str(self._cli_path.parent)  # (the CLI imports the `tools` package by absolute name)
            # the CLI imports the `tools` package, which lives two levels above its folder
            project_root = str(Path(__file__).resolve().parents[3])
            env = dict(os.environ)
            env["PYTHONPATH"] = project_root + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")

        try:
            proc = subprocess.Popen(
                cmd,
                cwd=cwd,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,  # merge stderr into stdout for a single live stream
                stdin=subprocess.DEVNULL,  # a frozen console=False app has no
                                            # valid stdin to inherit — without this
                                            # the child can hang waiting on input
                                            # it will never receive
                text=True,
                encoding="utf-8",
                errors="replace",
                bufsize=1,  # line-buffered
                **no_window_kwargs(),
            )
        except Exception as e:
            return {"status": "error", "message": f"Failed to launch the prediction subprocess: {e}"}

        try:
            for line in proc.stdout:
                line = line.rstrip("\n")
                if line and self.on_log:
                    try:
                        self.on_log(line)
                    except Exception:
                        pass
            proc.wait(timeout=TIMEOUT_SECONDS)
        except subprocess.TimeoutExpired:
            proc.kill()
            return {
                "status": "error",
                "message": f"Analysis did not finish within {TIMEOUT_SECONDS}s (timed out).",
            }

        summary_path = outdir / "summary.json"
        if not summary_path.exists():
            return {
                "status": "error",
                "message": f"Subprocess exited (code {proc.returncode}) without producing a result. "
                           "Check the log above for the actual error.",
            }

        with open(summary_path, "r", encoding="utf-8") as f:
            return json.load(f)
