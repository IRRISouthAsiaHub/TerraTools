from __future__ import annotations

import json
import subprocess
from pathlib import Path

from tools.common.r_locator import find_rscript, get_r_env, rscript_not_found_message
from tools.common.subprocess_utils import no_window_kwargs
from tools.qtl.backend.qtl_compare_plots import plot_common, plot_comparison, to_native
from tools.qtl.backend.qtl_regions import cluster_regions, load_regions

# One R run per model; ML models on big panels can take a while.
R_TIMEOUT_S = 3600


def _r_error_message(stderr, stdout):
    """R prints package start-up chatter before any error, so the actual
    reason used to be cut off in the app's status line. Lead with it."""
    text = (stderr or "").strip() or (stdout or "").strip()
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.startswith(("Error", "Fehler")):
            return "\n".join(l for l in lines[i:] if l.strip() != "Execution halted").strip()
    return "\n".join(lines[-15:]) or "R process failed."


BACKEND_DIR = Path(__file__).resolve().parent

CLI_FILE = BACKEND_DIR / "qtl_cli.R"


def _run_r(
    config: dict,
    results_dir: Path,
) -> dict:

    results_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    config_file = (
        results_dir /
        "_qtl_config.json"
    )

    config_file.write_text(
        json.dumps(
            config,
            indent=2
        ),
        encoding="utf-8"
    )

    try:

        rscript = find_rscript()
        if rscript is None:
            raise RuntimeError(rscript_not_found_message())

        try:
            process = subprocess.run(
                [
                    rscript,
                    str(CLI_FILE),
                    str(config_file),
                    str(results_dir),
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                encoding="utf-8",
                errors="replace",
                cwd=str(BACKEND_DIR),
                env=get_r_env(rscript),
                timeout=R_TIMEOUT_S,
                **no_window_kwargs(),
            )
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(
                f"R process timed out after {R_TIMEOUT_S}s. This usually means Rscript "
                "is stuck (e.g. trying to install a missing package with no "
                "internet access) rather than actually computing. Check that "
                "all required R packages are pre-installed.\n\n"
                f"Partial stdout:\n{(exc.stdout or '').strip()}\n\n"
                f"Partial stderr:\n{(exc.stderr or '').strip()}"
            ) from exc

        if process.returncode != 0:

            raise RuntimeError(
                _r_error_message(process.stderr, process.stdout)
            )

        output = (
            process.stdout.strip()
        )

        # R packages/functions may print
        # messages before the JSON.
        json_start = output.rfind(
            '{"status"'
        )

        if json_start == -1:

            raise RuntimeError(
                "R completed, but no JSON result was returned.\n\n"
                f"R output:\n{output}"
            )

        result = json.loads(
            output[json_start:]
        )

        return result

    finally:

        if config_file.exists():

            config_file.unlink()


def run_qtl(
    config: dict,
    backend_dir: Path | None = None,
    results_dir: Path | None = None,
) -> dict:

    if results_dir is None:

        results_dir = BACKEND_DIR / "results"

    result = _run_r(
        config,
        Path(results_dir)
    )

    return result


def _write_table(rows, path):
    """Save a comparison table as CSV (opens cleanly in Excel); returns its path."""
    import pandas as pd

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    if not df.empty:
        first = [c for c in ("Region_ID", "Chromosome") if c in df.columns]
        df = df[first + [c for c in df.columns if c not in first]]
    df.to_csv(path, index=False, encoding="utf-8-sig")
    return str(path.resolve()).replace("\\", "/")


def run_comparison(
    config: dict,
    backend_dir: Path | None = None,
    results_dir: Path | None = None,
) -> dict:

    if results_dir is None:

        results_dir = BACKEND_DIR / "results"

    results_dir = Path(
        results_dir
    )

    models = config.get(
        "models",
        []
    )

    if len(models) < 2:

        raise ValueError(
            "Select at least two QTL models."
        )

    model_results = {}
    failed = {}

    # --------------------------------------------------------
    # Run each selected model through the same R engine.
    # One model failing must not throw away the others.
    # --------------------------------------------------------

    for model in models:

        model_config = dict(
            config
        )

        model_config[
            "mapping_method"
        ] = model

        try:
            result = run_qtl(
                model_config,
                backend_dir,
                results_dir
            )
        except Exception as e:
            failed[model] = str(e)
            continue

        if result.get(
            "status"
        ) != "success":

            failed[model] = result.get("message", "did not complete")
            continue

        model_results[
            model
        ] = result

    if not model_results:

        raise RuntimeError(
            "None of the selected QTL models completed.\n\n"
            + "\n\n".join(f"{m}: {msg}" for m, msg in failed.items())
        )

    # --------------------------------------------------------
    # One row per QTL region found by at least one model; regions
    # from different models that overlap are the same QTL.
    # --------------------------------------------------------

    regions = load_regions(model_results)
    table = cluster_regions(regions)
    table.sort(key=lambda x: (-x["N_Models_Detected"], -x["Max_LOD"]))

    ideogram, bar = plot_comparison(model_results, results_dir, config.get("plot"), config.get("phenotype_column"))
    csv_path = _write_table(table, results_dir / "compare_regions.csv")

    return to_native({

        "csv":
            csv_path,

        "status":
            "success",

        "models":
            list(
                model_results.keys()
            ),

        "failed_models":
            failed,

        "model_results":
            model_results,

        "thresholds": {
            m: r.get("summary", {}).get("threshold_lod") for m, r in model_results.items()
        },

        "table":
            table,

        "ideogram":
            ideogram,

        "bar":
            bar,
    })


def find_common_qtl(
    config: dict,
    comparison_result: dict,
    backend_dir: Path | None = None,
    results_dir: Path | None = None,
) -> dict:

    min_models = int(
        config.get(
            "min_models",
            2
        )
    )

    model_results = comparison_result.get("model_results", {})
    out_dir = results_dir or BACKEND_DIR / "results"

    regions = load_regions(model_results)
    common = [c for c in cluster_regions(regions) if c["N_Models_Detected"] >= min_models]
    common.sort(key=lambda x: (-x["N_Models_Detected"], -x["Max_LOD"]))

    ideogram, bar = plot_common(common, model_results, min_models, out_dir, config.get("plot"))
    csv_path = _write_table(common, Path(out_dir) / "common_regions.csv")

    return to_native({
        "status": "success",
        "csv": csv_path,
        "table": common,
        "ideogram": ideogram,
        "bar": bar,
    })
