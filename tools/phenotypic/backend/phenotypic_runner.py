"""
phenotypic_runner.py — Python side of the Phenotypic Data Analysis tool's R backend.

Three calls, each one Rscript process (phenotypic_cli.R) that leaves its answer in <outdir>/result_<action>.json:
  inspect  reads and validates the data file (columns, environments, usable traits, design check)
  run      runs ONE analysis module; the results of earlier runs in the same folder are kept
  export   writes the Excel workbooks of everything run so far

Same deadlock-safe / timeout-safe subprocess pattern as the other tools: capture_output=True, an explicit timeout,
and a file (not stdout) for the result.
"""

import csv
import json
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

from tools.phenotypic.config import PhenotypicConfig
from tools.common.r_locator import find_rscript, rscript_not_found_message
from tools.common.subprocess_utils import no_window_kwargs

R_TIMEOUT_SECONDS = 900  # the mixed models of a large trial can take a while
OUT_PREFIX = "terratools_phenotypic_"


def _csv_to_records(path: str) -> list:
    # written with a byte-order mark so Excel shows H² and G×E correctly
    with open(path, "r", newline="", encoding="utf-8-sig") as f:
        return [row for row in csv.DictReader(f)]


def _p(path) -> str:
    return str(path).replace(chr(92), "/")


class PhenotypicEngine:
    def __init__(self, config: PhenotypicConfig = None):
        self.config = config
        self._cli_path = Path(__file__).resolve().parent / "phenotypic_cli.R"

    # ── one Rscript call ──────────────────────────────────────────────────────────────────
    def _call(self, action: str, outdir: Path, params: dict = None) -> dict:
        rscript = find_rscript()
        if rscript is None:
            return {"status": "error", "message": rscript_not_found_message()}
        outdir.mkdir(parents=True, exist_ok=True)
        result_path = outdir / f"result_{action}.json"
        result_path.unlink(missing_ok=True)

        cmd = [rscript, str(self._cli_path), "--action", action, "--outdir", _p(outdir)]
        if params is not None:
            params_path = outdir / f"params_{action}.json"
            params_path.write_text(json.dumps(params), encoding="utf-8")
            cmd += ["--params", _p(params_path)]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace",
                                  timeout=R_TIMEOUT_SECONDS, stdin=subprocess.DEVNULL, **no_window_kwargs())
        except subprocess.TimeoutExpired:
            return {"status": "error", "message": f"R process did not finish within {R_TIMEOUT_SECONDS}s (timed out)."}
        except Exception as e:
            return {"status": "error", "message": f"Failed to launch Rscript: {e}"}

        if not result_path.exists():
            return {"status": "error", "message": ("R process exited without producing a result "
                                                    f"(exit code {proc.returncode}). stderr:\n{proc.stderr[-2000:]}")}
        with open(result_path, "r", encoding="utf-8") as f:
            return json.load(f)

    # ── inspect ───────────────────────────────────────────────────────────────────────────
    def inspect(self) -> dict:
        tmp = Path(tempfile.mkdtemp(prefix="terratools_phenotypic_inspect_"))
        try:
            return self._call("inspect", tmp, self.config.params())
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    # ── run one analysis ──────────────────────────────────────────────────────────────────
    def run(self) -> dict:
        cfg = self.config
        given = Path(cfg.outdir) if cfg.outdir else None
        if given is not None and given.name.startswith(OUT_PREFIX) and given.is_dir():
            outdir = given
        else:
            outdir = Path(tempfile.mkdtemp(prefix=OUT_PREFIX))
        res = self._call("run", outdir, cfg.params())
        res["outdir"] = _p(outdir)
        if res.get("status") != "success":
            return res

        tables = {}
        for name, path in (res.get("tables") or {}).items():
            try:
                tables[name] = _csv_to_records(path)
            except Exception as e:
                tables[name] = []
                res.setdefault("msgs", []).append(
                    {"module": cfg.module, "type": "warning", "text": f"{name}: could not be read ({e})"})
        res["tables_data"] = tables

        wide = res.get("anova_wide")
        if wide:
            try:
                with open(wide, "r", encoding="utf-8") as f:
                    res["anova_wide_data"] = json.load(f)
            except Exception:
                res["anova_wide_data"] = []
        return res

    # ── downloads ─────────────────────────────────────────────────────────────────────────
    def export(self, outdir: str, kind: str) -> dict:
        out = Path(outdir) if outdir else None
        if out is None or not out.is_dir() or not out.name.startswith(OUT_PREFIX):
            return {"status": "error", "message": "Run an analysis first."}
        if kind == "figures":
            plots = out / "plots"
            files = sorted(p for p in plots.rglob("*") if p.suffix.lower() in (".png", ".pdf")) if plots.is_dir() else []
            if not files:
                return {"status": "error", "message": "There are no figures yet."}
            zp = out / "figures_300dpi.zip"
            with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED) as z:
                for p in files:
                    z.write(p, p.relative_to(plots).as_posix())
            return {"status": "success", "kind": kind, "path": _p(zp), "name": "figures_300dpi.zip"}

        wb = self._call("export", out)
        if wb.get("status") != "success":
            return wb
        if kind == "all":
            return {"status": "success", "kind": kind, "path": wb["excel_all"], "name": "Phenotypic_Analysis.xlsx"}
        if kind == "anova":
            if not wb.get("excel_anova"):
                return {"status": "error", "message": "Run the ANOVA first."}
            return {"status": "success", "kind": kind, "path": wb["excel_anova"], "name": "ANOVA_tables.xlsx"}
        # everything: tables, workbooks and figures in one zip (the internal state file stays out)
        zp = out / "phenotypic_results.zip"
        zp.unlink(missing_ok=True)
        with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED) as z:
            for p in sorted(out.rglob("*")):
                if (p.is_file() and p.name != "state.rds" and not p.name.startswith(("params_", "result_"))
                        and p.suffix.lower() != ".zip"):
                    z.write(p, p.relative_to(out).as_posix())
        return {"status": "success", "kind": kind, "path": _p(zp), "name": "phenotypic_results.zip"}
