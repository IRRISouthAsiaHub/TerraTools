from __future__ import annotations

import subprocess
from pathlib import Path


BACKEND_DIR = Path(__file__).resolve().parent
RESULTS_DIR_DEFAULT = BACKEND_DIR / "results"


# ============================================================
# TOOL DISCOVERY
# ============================================================

def _check_tool(name: str) -> None:

    try:

        subprocess.run(
            [name, "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=15,
        )

    except FileNotFoundError:

        raise RuntimeError(
            f"'{name}' was not found on PATH. QTL-seq requires "
            "bwa and samtools to be installed/bundled."
        )

    except Exception:

        # Some tools (e.g. bwa with no args) exit non-zero on
        # --version but still prove they exist.
        pass


def _run(cmd, **kwargs):

    print(
        "[QTL-SEQ] running:",
        " ".join(str(c) for c in cmd),
        flush=True,
    )

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        **kwargs,
    )

    if result.returncode != 0:

        raise RuntimeError(
            f"Command failed ({' '.join(str(c) for c in cmd)}):\n"
            + (result.stderr.strip() or result.stdout.strip())
        )

    return result


# ============================================================
# ALIGNMENT (FASTQ -> sorted, indexed BAM)
# ============================================================

def _align_sample(
    label: str,
    fastq: Path,
    ref: Path,
    work_dir: Path,
    threads: int,
    fastq_r2: Path | None = None,
) -> Path:

    sorted_bam = work_dir / f"{label}.sorted.bam"

    if fastq_r2 is not None:
        print(
            f"[QTL-SEQ] aligning {label} (paired-end): "
            f"{fastq} + {fastq_r2} -> {sorted_bam}",
            flush=True,
        )
    else:
        print(
            f"[QTL-SEQ] aligning {label} (single-end): "
            f"{fastq} -> {sorted_bam}",
            flush=True,
        )

    bwa_cmd = ["bwa", "mem", "-t", str(threads), str(ref), str(fastq)]
    if fastq_r2 is not None:
        bwa_cmd.append(str(fastq_r2))

    bwa = subprocess.Popen(
        bwa_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    sort = subprocess.Popen(
        [
            "samtools", "sort",
            "-@", str(threads),
            "-o", str(sorted_bam),
        ],
        stdin=bwa.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if bwa.stdout is not None:
        bwa.stdout.close()

    _, sort_err = sort.communicate()
    _, bwa_err = bwa.communicate()

    if sort.returncode != 0:

        raise RuntimeError(
            f"Alignment/sort failed for {label}:\n"
            f"{bwa_err.strip()}\n{sort_err.strip()}"
        )

    _run(["samtools", "index", str(sorted_bam)])

    return sorted_bam


def _prepare_bam(
    label: str,
    path: Path,
    ref: Path,
    work_dir: Path,
    threads: int,
    path_r2: Path | None = None,
) -> Path:
    """Accept either a raw FASTQ (aligned here, single- or paired-end)
    or an already-sorted BAM. path_r2, if given, is the mate/reverse
    FASTQ for paired-end alignment — ignored if path is a .bam."""

    if path.suffix.lower() in (".bam",):

        if path_r2 is not None:
            print(
                f"[QTL-SEQ] {label}: R2 file provided alongside a .bam "
                "input — ignoring R2, the BAM is already aligned.",
                flush=True,
            )

        if not (
            path.with_suffix(path.suffix + ".bai").exists()
            or Path(str(path) + ".bai").exists()
        ):

            _run(["samtools", "index", str(path)])

        return path

    return _align_sample(
        label, path, ref, work_dir, threads, fastq_r2=path_r2
    )


# ============================================================
# MPILEUP PARSING
# ============================================================

def _parse_mpileup_line(fields, min_depth):

    if len(fields) < 12:
        return None

    chrom, pos, ref_base = (
        fields[0],
        int(fields[1]),
        fields[2].upper(),
    )

    samples = []

    for i in range(3):

        start_idx = 3 + i * 3

        if start_idx + 2 >= len(fields):
            return None

        depth = int(fields[start_idx])
        reads = fields[start_idx + 1]

        if depth < min_depth:
            return None

        ref_count = reads.count(".") + reads.count(",")
        alt_count = depth - ref_count

        ref_freq = ref_count / depth if depth else 0.0
        alt_freq = alt_count / depth if depth else 0.0

        samples.append({
            "depth": depth,
            "ref_freq": ref_freq,
            "alt_freq": alt_freq,
        })

    return {
        "chrom": chrom,
        "pos": pos,
        "ref_base": ref_base,
        "parent": samples[0],
        "bulk1": samples[1],
        "bulk2": samples[2],
    }


def _sliding_window(snp_rows, window_size, step_size):

    import math

    by_chrom: dict[str, list[dict]] = {}

    for row in snp_rows:
        by_chrom.setdefault(row["CHROM"], []).append(row)

    windows = []

    for chrom, rows in by_chrom.items():

        rows = sorted(rows, key=lambda r: r["POS"])

        min_pos = rows[0]["POS"]
        max_pos = rows[-1]["POS"]

        start = min_pos
        made_any = False

        while start <= max_pos:

            end = start + window_size

            in_window = [
                r for r in rows
                if start <= r["POS"] < end
            ]

            if in_window:

                made_any = True

                deltas = [r["deltaSNP"] for r in in_window]
                n = len(deltas)
                mean_delta = sum(deltas) / n

                if n > 1:

                    variance = sum(
                        (d - mean_delta) ** 2 for d in deltas
                    ) / (n - 1)

                    std = math.sqrt(variance)
                    se = std / math.sqrt(n)

                else:

                    std = 0.0
                    se = 0.0

                windows.append({
                    "CHROM": chrom,
                    "START": start,
                    "END": end,
                    "POS": start + window_size // 2,
                    "N_SNPS": n,
                    "deltaSNP": round(mean_delta, 5),
                    "CI_95": round(1.96 * se, 5),
                    "CI_99": round(2.58 * se, 5),
                })

            start += step_size

        if not made_any:

            windows.append({
                "CHROM": chrom,
                "START": min_pos,
                "END": min_pos + window_size,
                "POS": min_pos + window_size // 2,
                "N_SNPS": 0,
                "deltaSNP": 0,
                "CI_95": 0,
                "CI_99": 0,
            })

    return windows


# ============================================================
# PLOT
# ============================================================

def _save_plot(snp_rows, window_rows, results_dir: Path) -> Path | None:

    if not window_rows:
        return None

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    chroms = sorted(
        {r["CHROM"] for r in window_rows}
    )

    n = len(chroms)
    ncols = min(4, n) or 1
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(4.2 * ncols, 3.2 * nrows),
        dpi=140,
        squeeze=False,
    )

    for i, chrom in enumerate(chroms):

        ax = axes[i // ncols][i % ncols]

        snp_chrom = [
            r for r in snp_rows if r["CHROM"] == chrom
        ]

        win_chrom = [
            r for r in window_rows if r["CHROM"] == chrom
        ]

        if snp_chrom:

            ax.scatter(
                [r["POS"] / 1e6 for r in snp_chrom],
                [r["deltaSNP"] for r in snp_chrom],
                s=2, color="#95d5b2", alpha=0.5,
                label="Raw",
            )

        if win_chrom:

            win_sorted = sorted(
                win_chrom, key=lambda r: r["POS"]
            )

            ax.plot(
                [r["POS"] / 1e6 for r in win_sorted],
                [r["deltaSNP"] for r in win_sorted],
                color="#d62828", linewidth=2,
                label="Smoothed",
            )

            ax.plot(
                [r["POS"] / 1e6 for r in win_sorted],
                [r["CI_95"] for r in win_sorted],
                color="#f4a261", linewidth=1,
                linestyle="--", label="95% CI",
            )

            ax.plot(
                [r["POS"] / 1e6 for r in win_sorted],
                [-r["CI_95"] for r in win_sorted],
                color="#f4a261", linewidth=1,
                linestyle="--",
            )

        ax.axhline(0, color="black", linewidth=1)
        ax.set_ylim(-1.1, 1.1)
        ax.set_title(f"Chr {chrom}", fontsize=10)
        ax.set_xlabel("Position (Mb)", fontsize=8)
        ax.set_ylabel("\u0394(SNP-index)", fontsize=8)
        ax.tick_params(labelsize=7)

        if i == 0:
            ax.legend(fontsize=7, loc="upper right")

    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")

    fig.suptitle(
        "QTL-seq: \u0394(SNP-index) across chromosomes",
        fontsize=12,
    )
    fig.tight_layout()

    plot_file = results_dir / "qtlseq_manhattan.png"
    fig.savefig(plot_file)
    plt.close(fig)

    return plot_file


# ============================================================
# MAIN ENTRY POINT
# ============================================================

def run_qtlseq(
    config: dict,
    results_dir: Path | None = None,
) -> dict:

    ref = Path(config.get("ref_file", ""))

    # R1 is required for each sample; R2 is optional (single-end reads
    # or an already-aligned .bam remain fully supported when omitted).
    def _optional_path(key: str) -> Path | None:
        val = config.get(key, "")
        return Path(val) if val else None

    parent_path = Path(config.get("parent_file_r1", ""))
    parent_path_r2 = _optional_path("parent_file_r2")

    bulk1_path = Path(config.get("bulk1_file_r1", ""))
    bulk1_path_r2 = _optional_path("bulk1_file_r2")

    bulk2_path = Path(config.get("bulk2_file_r1", ""))
    bulk2_path_r2 = _optional_path("bulk2_file_r2")

    for label, p in (
        ("reference", ref),
        ("parent", parent_path),
        ("bulk1", bulk1_path),
        ("bulk2", bulk2_path),
    ):

        if not p.exists():

            return {
                "status": "error",
                "message": f"{label} file not found: {p}",
            }

    for label, p in (
        ("parent R2", parent_path_r2),
        ("bulk1 R2", bulk1_path_r2),
        ("bulk2 R2", bulk2_path_r2),
    ):

        if p is not None and not p.exists():

            return {
                "status": "error",
                "message": f"{label} file not found: {p}",
            }

    min_depth = int(config.get("min_depth", 10))
    min_total_depth = int(config.get("min_total_depth", 50))
    ref_allele_freq = float(config.get("ref_allele_freq", 0.1))
    window_size = int(config.get("window_size", 1_000_000))
    step_size = int(config.get("step_size", 100_000))
    threads = int(config.get("threads", 2))

    if results_dir is None:
        results_dir = RESULTS_DIR_DEFAULT

    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    try:

        _check_tool("bwa")
        _check_tool("samtools")

        # ------------------------------------------------
        # Index reference (skip if already indexed).
        # ------------------------------------------------

        if not Path(str(ref) + ".bwt").exists():

            print("[QTL-SEQ] indexing reference...", flush=True)
            _run(["bwa", "index", str(ref)])

        if not Path(str(ref) + ".fai").exists():
            _run(["samtools", "faidx", str(ref)])

        # ------------------------------------------------
        # Align (or reuse) BAMs.
        # ------------------------------------------------

        parent_bam = _prepare_bam(
            "parent", parent_path, ref, results_dir, threads,
            path_r2=parent_path_r2,
        )

        bulk1_bam = _prepare_bam(
            "bulk1", bulk1_path, ref, results_dir, threads,
            path_r2=bulk1_path_r2,
        )

        bulk2_bam = _prepare_bam(
            "bulk2", bulk2_path, ref, results_dir, threads,
            path_r2=bulk2_path_r2,
        )

        # ------------------------------------------------
        # mpileup -> parse -> filter -> delta-SNP-index.
        # ------------------------------------------------

        print("[QTL-SEQ] running mpileup...", flush=True)

        mpileup = subprocess.Popen(
            [
                "samtools", "mpileup",
                "-f", str(ref),
                "-q", "20", "-Q", "20",
                str(parent_bam), str(bulk1_bam), str(bulk2_bam),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        snp_rows = []

        if mpileup.stdout is not None:

            for line_num, line in enumerate(mpileup.stdout):

                if line_num % 50_000 == 0 and line_num:

                    print(
                        f"[QTL-SEQ] processed {line_num} "
                        "mpileup lines...",
                        flush=True,
                    )

                parsed = _parse_mpileup_line(
                    line.rstrip("\n").split("\t"),
                    min_depth,
                )

                if parsed is None:
                    continue

                bulk1_idx = parsed["bulk1"]["alt_freq"]
                bulk2_idx = parsed["bulk2"]["alt_freq"]
                total_depth = (
                    parsed["bulk1"]["depth"]
                    + parsed["bulk2"]["depth"]
                )

                if total_depth < min_total_depth:
                    continue

                if parsed["parent"]["ref_freq"] < ref_allele_freq:
                    continue

                snp_rows.append({
                    "CHROM": parsed["chrom"],
                    "POS": parsed["pos"],
                    "REF": parsed["ref_base"],
                    "PARENT_DEPTH": parsed["parent"]["depth"],
                    "BULK1_DEPTH": parsed["bulk1"]["depth"],
                    "BULK1_SNP_INDEX": round(bulk1_idx, 5),
                    "BULK2_DEPTH": parsed["bulk2"]["depth"],
                    "BULK2_SNP_INDEX": round(bulk2_idx, 5),
                    "deltaSNP": round(bulk1_idx - bulk2_idx, 5),
                    "TOTAL_DEPTH": total_depth,
                })

        mpileup_err = ""

        if mpileup.stderr is not None:
            mpileup_err = mpileup.stderr.read()

        return_code = mpileup.wait()

        if return_code != 0:

            return {
                "status": "error",
                "message":
                    "samtools mpileup failed:\n"
                    + mpileup_err.strip(),
            }

        if not snp_rows:

            return {
                "status": "error",
                "message":
                    "No SNPs passed the quality/depth filters. "
                    "Try lowering the depth or allele-frequency "
                    "thresholds.",
            }

        print(
            f"[QTL-SEQ] {len(snp_rows)} SNPs passed filters",
            flush=True,
        )

        # ------------------------------------------------
        # Sliding window analysis.
        # ------------------------------------------------

        window_rows = _sliding_window(
            snp_rows, window_size, step_size
        )

        print(
            f"[QTL-SEQ] {len(window_rows)} windows",
            flush=True,
        )

        # ------------------------------------------------
        # Save CSVs.
        # ------------------------------------------------

        import csv

        snp_csv = results_dir / "qtlseq_snp_data.csv"

        with open(snp_csv, "w", newline="") as f:

            writer = csv.DictWriter(
                f, fieldnames=list(snp_rows[0].keys())
            )
            writer.writeheader()
            writer.writerows(snp_rows)

        window_csv = results_dir / "qtlseq_window_data.csv"

        with open(window_csv, "w", newline="") as f:

            writer = csv.DictWriter(
                f, fieldnames=list(window_rows[0].keys())
            )
            writer.writeheader()
            writer.writerows(window_rows)

        # ------------------------------------------------
        # Plot.
        # ------------------------------------------------

        plot_file = _save_plot(
            snp_rows, window_rows, results_dir
        )

        # ------------------------------------------------
        # Summary + capped table for the UI (top by |delta|).
        # ------------------------------------------------

        chroms = sorted({r["CHROM"] for r in snp_rows})

        top_snps = sorted(
            snp_rows,
            key=lambda r: abs(r["deltaSNP"]),
            reverse=True,
        )[:1000]

        candidate_windows = [
            r for r in window_rows
            if abs(r["deltaSNP"]) > r.get("CI_95", 0)
            and abs(r["deltaSNP"]) > 0.3
        ]

        return {

            "status": "success",

            "summary": {
                "n_snps": len(snp_rows),
                "n_windows": len(window_rows),
                "n_chromosomes": len(chroms),
                "n_candidate_windows": len(candidate_windows),
                "window_size": window_size,
                "step_size": step_size,
            },

            "chromosomes": chroms,

            "snp_table": top_snps,

            "window_table": window_rows,

            "candidate_windows": candidate_windows,

            "plot": str(plot_file) if plot_file else None,

            "snp_csv": str(snp_csv),

            "window_csv": str(window_csv),
        }

    except Exception as e:

        import traceback
        traceback.print_exc()

        return {
            "status": "error",
            "message": str(e),
        }