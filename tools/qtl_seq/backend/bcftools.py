"""bcftools on the VCF that the analysis writes: index it (tabix index, .tbi), count its records and read the SNP statistics.

Indexing is what makes the file usable in a genome browser (IGV) or with `bcftools view -r chr:start-end`; the record count
and the transition / transversion ratio are an independent check of the file. Everything here is optional: without bcftools
the VCF is written as usual."""

from __future__ import annotations

import subprocess

from tools.qtl_seq.backend import deps


def _run(*args, timeout=600):
    return subprocess.run(deps.command("bcftools", *args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=timeout,
                          **deps.run_kwargs())


def parse_stats(text):
    """{'records', 'snps', 'ts', 'tv', 'ts_tv'} from the output of `bcftools stats`."""
    out = {}
    for line in text.splitlines():
        f = line.split("\t")
        if f[0] == "SN" and len(f) >= 4:
            key = {"number of records:": "records", "number of SNPs:": "snps"}.get(f[2])
            if key:
                out[key] = int(f[3])
        elif f[0] == "TSTV" and len(f) >= 5:
            out["ts"], out["tv"] = int(f[2]), int(f[3])
            out["ts_tv"] = None if f[4] in ("nan", "inf") else float(f[4])
    return out


def index_and_check(path, expected_records=None):
    """Indexes `path` and checks it. Never raises: returns a dict with 'available', 'ok', 'message' and the numbers found."""
    info = {"available": False, "ok": False, "indexed": False, "message": ""}
    if not deps.where("bcftools"):
        info["message"] = "bcftools was not found, so the filtered VCF was not indexed or checked."
        return info
    info["available"] = True
    info["version"] = deps.version("bcftools")
    try:
        r = _run("index", "-f", "-t", path)
        if r.returncode != 0:
            info["message"] = "bcftools could not index the filtered VCF: " + (r.stderr or "").strip()[-300:]
            return info
        info["indexed"] = True
        r = _run("index", "-n", path)
        n = int(r.stdout.strip()) if r.returncode == 0 and r.stdout.strip().isdigit() else None
        r = _run("stats", path)
        stats = parse_stats(r.stdout) if r.returncode == 0 else {}
        info.update(stats)
        info["records_indexed"] = n
        if n is None or (expected_records is not None and n != expected_records):
            info["message"] = f"bcftools counts {n} records in the filtered VCF but {expected_records} SNPs were retained."
            return info
        info["ok"] = True
        tstv = f", Ts/Tv {info['ts_tv']:.2f}" if info.get("ts_tv") is not None else ""
        info["message"] = f"{info['version'] or 'bcftools'}: the filtered VCF is indexed and holds {n:,} records{tstv}."
    except (OSError, subprocess.SubprocessError, ValueError) as e:
        info["message"] = f"bcftools could not check the filtered VCF: {e}"
    return info


def region_count(path, region):
    """Number of records of an indexed VCF inside `region` ('chr:start-end'), or None."""
    try:
        r = _run("view", "-H", "-r", region, path)
        return len(r.stdout.splitlines()) if r.returncode == 0 else None
    except (OSError, subprocess.SubprocessError):
        return None

