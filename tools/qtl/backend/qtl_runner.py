from __future__ import annotations

import json
import subprocess
from pathlib import Path

from tools.common.r_locator import find_rscript, get_r_env, rscript_not_found_message
from tools.common.subprocess_utils import no_window_kwargs
from tools.qtl.backend.qtl_compare_plots import plot_common, plot_comparison, to_native
from tools.qtl.backend.qtl_regions import cluster_regions, load_regions

# One R run per model; ML models on big panels can take a while.
R_TIMEOUT_S = 1800


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


def _read_linkage_map(path):
    """Marker, Chromosome, Position (cM) of a Linkage Mapping output (Linkage Map / Retained Markers /
    All Markers sheet of linkage_results.xlsx, retained_markers.xlsx, all_markers.xlsx, or a CSV);
    markers that were removed from the map are dropped."""
    import pandas as pd

    if str(path).lower().endswith((".csv", ".tsv", ".txt")):
        frames = [pd.read_csv(path, sep=None, engine="python")]
    else:
        xl = pd.ExcelFile(path)
        names = sorted(xl.sheet_names, key=lambda n: n not in ("Linkage Map", "Retained Markers"))
        frames = [xl.parse(n) for n in names]
    for df in frames:
        cols = {str(c).strip().lower(): c for c in df.columns}
        if {"marker", "chromosome", "position"} <= set(cols):
            df = df.rename(columns={cols["marker"]: "Marker", cols["chromosome"]: "Chromosome",
                                    cols["position"]: "Position"})
            if "status" in cols:
                df = df[df[cols["status"]].astype(str).str.lower() == "retained"]
            df["Position"] = pd.to_numeric(df["Position"], errors="coerce")
            return df.dropna(subset=["Position"])[["Marker", "Chromosome", "Position"]]
    raise ValueError("The linkage map needs the columns Marker, Chromosome and Position.")


def _physical_sheet(path):
    """{marker: (chromosome, bp)} from the 'Physical Map' sheet that Linkage Mapping writes (bp read from the
    marker names, the rest interpolated from cM); {} for any other file."""
    import pandas as pd

    if not path or str(path).lower().endswith((".csv", ".tsv", ".txt")):
        return {}
    try:
        df = pd.read_excel(path, sheet_name="Physical Map")
        out = {}
        src = df["Physical_Source"] if "Physical_Source" in df.columns else [""] * len(df)
        for m, c, b, how in zip(df["Marker"], df["Physical_Chromosome"], df["Physical_bp"], src):
            # "cM x population factor" has no known chromosome (the group number is only a placeholder)
            if pd.notna(c) and pd.notna(b) and str(how) != "cM x population factor":
                out[str(m).strip()] = (int(c), int(b))
        return out
    except Exception:
        return {}


def _marker_map(genotype_file, linkage_map_file=None):
    """({linkage group: [(cM, chromosome, bp, marker), ...]}, note): the markers whose NAME carries a
    physical position (chr02_16919033 ...). The chromosome comes from the name, because a linkage group
    number is not always the chromosome number. The cM positions come from the Linkage Mapping output when
    one is given, and otherwise from the first three rows of the genotype file. A linkage map that does not
    agree with the genotype file's positions is not used: the QTL intervals are in the genotype file's cM."""
    import csv

    from tools.qtl.backend.qtl_regions import physical_position

    with open(genotype_file, newline="", encoding="utf-8-sig", errors="replace") as fh:
        rows = []
        for row in csv.reader(fh):
            if any(c.strip() for c in row):
                rows.append([c.strip() for c in row])
            if len(rows) == 3:
                break
    if len(rows) < 3:
        return {}, ""

    def collect(items):
        out = {}
        for name, group, cm in items:
            try:
                cm = float(cm)
                key = str(int(float(group)))
            except (TypeError, ValueError):
                continue
            pos = phys.get(str(name).strip()) or physical_position(name)
            if pos is not None:
                out.setdefault(key, []).append((cm, pos[0], pos[1], str(name)))
        return {k: sorted(v) for k, v in out.items()}

    phys = _physical_sheet(linkage_map_file)
    geno = list(zip(rows[0][1:], rows[1][1:], rows[2][1:]))
    note = ""
    if linkage_map_file:
        try:
            lm = _read_linkage_map(linkage_map_file)
            have = {n: c for n, _ch, c in geno}
            shared = [(have[m], float(p)) for m, p in zip(lm["Marker"], lm["Position"]) if m in have]
            same = sum(1 for g, p in shared if _close(g, p))
            if shared and same >= 0.8 * len(shared):
                return collect(zip(lm["Marker"], lm["Chromosome"], lm["Position"])), "Marker positions from the linkage map."
            note = ("The linkage map does not match the genotype file (" +
                    ("no shared markers" if not shared else f"{len(shared) - same} of {len(shared)} shared markers "
                     "have a different cM position") + "), so the genotype file's own marker positions were used.")
        except Exception as e:
            note = f"The linkage map could not be used ({e}); the genotype file's own marker positions were used."
    return collect(geno), note


def _close(a, b):
    try:
        return abs(float(a) - float(b)) <= 0.05
    except (TypeError, ValueError):
        return False


def _cm_to_bp(points, cm_start, cm_end):
    """Physical span (bp) of a cM interval: the ends are placed by linear interpolation between the
    neighbouring markers, and every marker inside the interval is kept so that a map which runs
    against the genome still gives start <= end."""
    import numpy as np

    cms = np.array([p[0] for p in points], dtype=float)
    bps = np.array([p[1] for p in points], dtype=float)
    lo, hi = sorted((float(cm_start), float(cm_end)))
    ends = [float(np.interp(lo, cms, bps)), float(np.interp(hi, cms, bps))]
    inside = bps[(cms >= lo) & (cms <= hi)]
    values = ends + [float(v) for v in inside]
    return int(round(min(values))), int(round(max(values)))


def _qtl_interval(q, pts):
    """(chromosome, start bp, end bp, remark) of one QTL, or (chromosome or None, None, None, reason).

    The chromosome is read from the marker names, not taken from the linkage group number: first from the
    QTL's own flanking/peak markers, otherwise the one most of the group's markers belong to. The cM
    interval is converted with the markers of that chromosome on this group; when the group has too few of
    them around the QTL, the flanking markers' own positions give the span."""
    from collections import Counter

    from tools.qtl.backend.qtl_regions import physical_position

    group = str(int(float(q.get("Chromosome"))))
    lo, hi = sorted((float(q.get("CI_Start_cM")), float(q.get("CI_End_cM"))))
    by_name = {m[3]: (m[1], m[2]) for m in pts}
    flank = [(k, physical_position(q.get(k)) or by_name.get(str(q.get(k)).strip())) for k in ("Left_Marker", "Peak_Marker", "Right_Marker")]
    flank = [(k, p) for k, p in flank if p]
    votes = Counter(p[0] for _k, p in flank) or Counter(m[1] for m in pts)
    if not votes:
        return None, None, None, "none of its markers has a position in its name"
    chrom = votes.most_common(1)[0][0]
    remark = []
    if chrom != int(group):
        remark.append(f"linkage group {group} is chromosome {chrom} by its marker names")

    on_chr = [m for m in pts if m[1] == chrom]
    cms = [m[0] for m in on_chr]
    if len(on_chr) >= 2 and max(cms) >= lo and min(cms) <= hi:
        start, end = _cm_to_bp([(m[0], m[2]) for m in on_chr], lo, hi)
        if end > start:
            f = [p[1] for _k, p in flank if p[0] == chrom]
            if len(f) == 3 and not (f[0] <= f[1] <= f[2] or f[0] >= f[1] >= f[2]):
                remark.append("the map order of its markers differs from the physical order, so the interval may be wide")
            return chrom, start, end, "; ".join(remark)

    f = sorted({p[1] for _k, p in flank if p[0] == chrom})
    if len(f) >= 2:
        remark.append("taken from the positions of its flanking markers")
        return chrom, f[0], f[-1], "; ".join(remark)
    names = [str(q.get(k)) for k in ("Left_Marker", "Peak_Marker", "Right_Marker")]
    return chrom, None, None, (f"chromosome {chrom} has fewer than two markers with a physical position around it "
                               f"(flanking markers: {', '.join(names)})")


def write_interval_file(config, result, results_dir):
    """QTL_ID, Chromosome, QTL_Start, QTL_End (bp) for every QTL of this run, from its confidence interval.

    The mapping works in cM, a genome annotation in bp, so the CI is converted through the markers'
    physical positions (see _qtl_interval). Chromosome is the physical chromosome named by the markers.
    QTLs that cannot be converted are listed in 'interval_skipped', with the reason in 'interval_remarks'.
    Adds interval_csv / interval_skipped / interval_remarks to the result."""
    import csv

    table = result.get("qtl_table") or []
    method = (result.get("summary") or {}).get("method") or config.get("mapping_method") or "qtl"
    path = Path(results_dir) / f"qtl_intervals_{method}.csv"

    try:
        markers, note = _marker_map(config.get("genotype_file"), config.get("linkage_map_file"))
    except Exception:
        markers, note = {}, ""

    rows, skipped, remarks = [], [], []
    for q in table:
        qid = q.get("QTL_ID")
        try:
            chrom, start, end, why = _qtl_interval(q, markers.get(str(int(float(q.get("Chromosome")))), []))
        except (TypeError, ValueError):
            skipped.append(qid)
            remarks.append(f"{qid} is not in the file: no usable chromosome or confidence interval")
            continue
        if start is None:
            skipped.append(qid)
            remarks.append(f"{qid} is not in the file: {why}")
            continue
        rows.append((qid, chrom, start, end))
        if why:
            remarks.append(f"{qid} (chromosome {chrom}): {why}")

    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["QTL_ID", "Chromosome", "QTL_Start", "QTL_End"])
        w.writerows(rows)

    result["interval_csv"] = str(path.resolve()).replace("\\", "/")
    result["interval_skipped"] = skipped
    result["interval_remarks"] = remarks
    result["interval_note"] = note
    return result


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

    if result.get("status") == "success":
        write_interval_file(config, result, Path(results_dir))

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

    ideogram, bar = plot_comparison(model_results, results_dir)
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

    ideogram, bar = plot_common(common, model_results, min_models, out_dir)
    csv_path = _write_table(common, Path(out_dir) / "common_regions.csv")

    return to_native({
        "status": "success",
        "csv": csv_path,
        "table": common,
        "ideogram": ideogram,
        "bar": bar,
    })
