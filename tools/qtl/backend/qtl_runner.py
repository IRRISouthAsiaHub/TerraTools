from __future__ import annotations

import json
import subprocess
from pathlib import Path

from tools.common.r_locator import find_rscript, rscript_not_found_message
from tools.common.subprocess_utils import no_window_kwargs


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
                cwd=str(BACKEND_DIR),
                timeout=300,  # seconds; tune to your largest expected dataset
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


def run_qtl(
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

    # --------------------------------------------------------
    # Run each selected model through the same R engine.
    # --------------------------------------------------------

    for model in models:

        model_config = dict(
            config
        )

        model_config[
            "mapping_method"
        ] = model

        result = run_qtl(
            model_config,
            backend_dir,
            results_dir
        )

        if result.get(
            "status"
        ) != "success":

            continue

        model_results[
            model
        ] = result

    if not model_results:

        raise RuntimeError(
            "None of the selected QTL models completed."
        )

    # --------------------------------------------------------
    # Build comparison table from each model's full scan.
    # --------------------------------------------------------

    comparison_rows = []

    for model, result in model_results.items():

        full_csv = result.get(
            "full_csv"
        )

        if not full_csv:

            continue

        full_path = Path(
            full_csv
        )

        if not full_path.exists():

            continue

        import pandas as pd

        df = pd.read_csv(
            full_path
        )

        if df.empty:

            continue

        df = df.copy()

        df["Model"] = model

        threshold = float(
            config.get(
                "lod_threshold",
                2.5
            )
        )

        df = df[
            pd.to_numeric(
                df["LOD"],
                errors="coerce"
            ) >= threshold
        ]

        comparison_rows.append(
            df
        )

    if comparison_rows:

        import pandas as pd

        combined = pd.concat(
            comparison_rows,
            ignore_index=True
        )

        table = []

        for marker, group in combined.groupby(
            "Marker",
            sort=False
        ):

            lod = pd.to_numeric(
                group["LOD"],
                errors="coerce"
            )

            pve = pd.to_numeric(
                group.get(
                    "PVE_Marker_Percent"
                ),
                errors="coerce"
            )

            row = group.iloc[
                lod.argmax()
            ]

            table.append({

                "Marker":
                    marker,

                "Chromosome":
                    row.get(
                        "Chromosome"
                    ),

                "Position":
                    row.get(
                        "Position"
                    ),

                "Models_List":
                    ", ".join(
                        sorted(
                            group[
                                "Model"
                            ].unique()
                        )
                    ),

                "Mean_LOD":
                    round(
                        lod.mean(),
                        3
                    ),

                "Max_LOD":
                    round(
                        lod.max(),
                        3
                    ),

                "Mean_PVE":
                    round(
                        pve.mean(),
                        3
                    ),
            })

        table.sort(
            key=lambda x:
                -float(
                    x["Max_LOD"]
                    or 0
                )
        )

    else:

        table = []

    return {

        "status":
            "success",

        "models":
            list(
                model_results.keys()
            ),

        "model_results":
            model_results,

        "table":
            table
    }


def find_common_qtl(
    config: dict,
    comparison_result: dict,
    backend_dir: Path | None = None,
    results_dir: Path | None = None,
) -> dict:

    import pandas as pd

    threshold = float(
        config.get(
            "lod_threshold",
            2.5
        )
    )

    min_models = int(
        config.get(
            "min_models",
            2
        )
    )

    rows = []

    for model, result in comparison_result.get(
        "model_results",
        {}
    ).items():

        csv_file = result.get(
            "full_csv"
        )

        if not csv_file:
            continue

        path = Path(
            csv_file
        )

        if not path.exists():
            continue

        df = pd.read_csv(
            path
        )

        if "LOD" not in df.columns:
            continue

        df["LOD"] = pd.to_numeric(
            df["LOD"],
            errors="coerce"
        )

        df = df[
            df["LOD"] >= threshold
        ].copy()

        if df.empty:
            continue

        df["Model"] = model

        rows.append(
            df
        )

    if not rows:

        return {
            "status": "success",
            "table": []
        }

    all_qtl = pd.concat(
        rows,
        ignore_index=True
    )

    common = []

    for marker, group in all_qtl.groupby(
        "Marker",
        sort=False
    ):

        models = sorted(
            group["Model"].unique()
        )

        if len(models) < min_models:
            continue

        lod = pd.to_numeric(
            group["LOD"],
            errors="coerce"
        )

        row = group.iloc[
            lod.argmax()
        ]

        common.append({

            "Marker":
                marker,

            "Chromosome":
                row.get(
                    "Chromosome"
                ),

            "Position":
                row.get(
                    "Position"
                ),

            "N_Models_Detected":
                len(models),

            "Models_List":
                ", ".join(models),

            "Mean_LOD":
                round(
                    lod.mean(),
                    3
                ),

            "Max_LOD":
                round(
                    lod.max(),
                    3
                ),

            "Mean_PVE":
                round(
                    pd.to_numeric(
                        group.get(
                            "PVE_Marker_Percent"
                        ),
                        errors="coerce"
                    ).mean(),
                    3
                )
        })

    common.sort(
        key=lambda x: (
            -x["N_Models_Detected"],
            -x["Max_LOD"]
        )
    )

    return {

        "status":
            "success",

        "table":
            common
    }