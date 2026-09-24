from __future__ import annotations

import json
import subprocess
from pathlib import Path

from tools.common.r_locator import find_rscript, rscript_not_found_message
from tools.common.subprocess_utils import no_window_kwargs


BACKEND_DIR = Path(__file__).resolve().parent

CLI_FILE = BACKEND_DIR / "cross_cli.R"


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
        "_cross_config.json"
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
                cwd=str(BACKEND_DIR),
                timeout=300,  # seconds; tune to your largest expected parent set
                **no_window_kwargs(),
            )
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(
                "R process timed out after 300s. This usually means Rscript "
                "is stuck (e.g. trying to install a missing package with no "
                "internet access) rather than actually computing. Check that "
                "all required R packages are pre-installed.\n\n"
                f"Partial stdout:\n{(exc.stdout or '').strip()}\n\n"
                f"Partial stderr:\n{(exc.stderr or '').strip()}"
            ) from exc

        if process.returncode != 0:

            raise RuntimeError(
                process.stderr.strip()
                or process.stdout.strip()
                or "R process failed."
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


def run_cross_selection(
    config: dict,
    backend_dir: Path | None = None,
    results_dir: Path | None = None,
) -> dict:

    if results_dir is None:

        results_dir = BACKEND_DIR / "results"

    return _run_r(
        config,
        Path(results_dir)
    )


# ================================================================
# FAST PYTHON-SIDE PREVIEW
# Mirrors the Shiny app's "Data Validation" panel, but runs
# straight off pandas so the UI gets feedback before the (much
# heavier) R thinning + diallel pipeline is ever invoked.
# ================================================================

def preview_data(
    criterion_file: str,
    kmat_file: str | None = None,
) -> dict:

    import pandas as pd

    result: dict = {
        "status": "success"
    }

    if criterion_file:

        crit = pd.read_csv(criterion_file)

        if "Genotype" not in crit.columns or "Criterion" not in crit.columns:

            return {
                "status": "error",
                "message": (
                    "Criterion file must have columns: Genotype, Criterion. "
                    f"Found: {', '.join(crit.columns.astype(str))}"
                )
            }

        crit["Genotype"] = crit["Genotype"].astype(str).str.strip()

        result["criterion"] = {
            "rows": int(len(crit)),
            "unique_genotypes": int(crit["Genotype"].nunique()),
            "preview": crit.head(8).to_dict(orient="records")
        }

    if kmat_file:

        kmat_raw = pd.read_csv(kmat_file)

        genotypes = (
            kmat_raw.iloc[:, 0]
            .astype(str)
            .str.strip()
        )

        values = kmat_raw.iloc[:, 1:].apply(
            pd.to_numeric,
            errors="coerce"
        )

        result["kmat"] = {
            "rows": int(kmat_raw.shape[0]),
            "cols": int(kmat_raw.shape[1] - 1),
            "min": None if values.empty else float(values.min().min()),
            "max": None if values.empty else float(values.max().max())
        }

        if criterion_file:

            overlap = int(
                crit["Genotype"].isin(genotypes).sum()
            )

            result["overlap"] = overlap

    return result


# ================================================================
# FAST PYTHON-SIDE TABLE FILTERING
# The R engine already wrote the full cross table to disk once;
# re-filtering it for the "All Crosses" tab doesn't need another
# R invocation.
# ================================================================

def filter_crosses(
    crosses_csv: str,
    relatedness: str = "All",
    performance: str = "All",
    kinship_max: float = 1.0,
    kinship_min: float = 0.0,
) -> dict:

    import pandas as pd

    path = Path(crosses_csv)

    if not path.exists():

        return {
            "status": "error",
            "message": "Run an analysis first."
        }

    df = pd.read_csv(path)

    if relatedness and relatedness != "All":
        df = df[df["Relatedness_Category"] == relatedness]

    if performance and performance != "All":
        df = df[df["Performance_Category"] == performance]

    kinship = pd.to_numeric(df["Kinship_Coefficient"], errors="coerce")

    df = df[
        (kinship <= float(kinship_max)) &
        (kinship >= float(kinship_min))
    ]

    return {
        "status": "success",
        "count": int(len(df)),
        "rows": df.to_dict(orient="records")
    }
