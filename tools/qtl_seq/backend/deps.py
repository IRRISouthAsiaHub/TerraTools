"""Which external programs are available (samtools, bcftools, pysam) and how to call them.

On Windows a program that is not on PATH is looked for inside the default WSL distribution (through a login shell, so a
tool installed under ~/.local/bin is found too); Windows paths are translated to /mnt/<drive>/... for it. pysam has no
Windows build: there it is used from the Python of the WSL distribution."""

from __future__ import annotations

import re
import shutil
import subprocess
import sys

_WIN = sys.platform.startswith("win")
NO_WINDOW = getattr(subprocess, "CREATE_NO_WINDOW", 0)     # keeps a console from flashing up in the windowed app
_WSL: dict = {}


def run_kwargs():
    return {"creationflags": NO_WINDOW} if _WIN else {}


def _wsl_shell(script, *args):
    """wsl.exe command running `script` in a login shell; $0, $@ are the arguments after it."""
    return ["wsl.exe", "-e", "bash", "-lc", script, *args]


def _wsl_has(name):
    if not _WIN:
        return False
    if name not in _WSL:
        try:
            r = subprocess.run(_wsl_shell('command -v "$0"', name), stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60,
                               **run_kwargs())
            _WSL[name] = r.returncode == 0
        except Exception:
            _WSL[name] = False
    return _WSL[name]


def where(name):
    """'path', 'wsl' or ''."""
    if shutil.which(name):
        return "path"
    return "wsl" if _wsl_has(name) else ""


def wsl_path(arg):
    s = str(arg)
    m = re.match(r"^([A-Za-z]):[\\/](.*)$", s)
    return f"/mnt/{m.group(1).lower()}/{m.group(2).replace(chr(92), '/')}" if m else s


def command(name, *args):
    """Command line list for `name` (from PATH, otherwise through WSL)."""
    if shutil.which(name) is None and _wsl_has(name):
        return _wsl_shell('exec "$0" "$@"', name, *[wsl_path(a) for a in args])
    return [name, *[str(a) for a in args]]


def version(name):
    """First line of `name --version`, or ''."""
    try:
        r = subprocess.run(command(name, "--version"), stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60, text=True,
                           **run_kwargs())
        return (r.stdout or "").splitlines()[0].strip() if r.returncode == 0 and r.stdout else ""
    except Exception:
        return ""


# ---------------------------------------------------------------- pysam
def have_pysam():
    """pysam importable by this Python (Linux / macOS)."""
    try:
        import pysam  # noqa: F401
        return True
    except Exception:
        return False


def pysam_where():
    """'python' (importable here), 'wsl' (importable by the Python of WSL) or ''."""
    if have_pysam():
        return "python"
    if not _WIN:
        return ""
    if "pysam" not in _WSL:
        try:
            r = subprocess.run(_wsl_shell('exec python3 -c "import pysam"'), stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=90,
                               **run_kwargs())
            _WSL["pysam"] = r.returncode == 0
        except Exception:
            _WSL["pysam"] = False
    return "wsl" if _WSL["pysam"] else ""


def pysam_version():
    w = pysam_where()
    if w == "python":
        import pysam
        return pysam.__version__
    if w == "wsl":
        try:
            r = subprocess.run(_wsl_shell('exec python3 -c "import pysam;print(pysam.__version__)"'), stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, timeout=90, text=True, **run_kwargs())
            return r.stdout.strip()
        except Exception:
            return ""
    return ""


def pysam_script_command(script, *args):
    """Command line running a stand-alone script with the Python that has pysam (WSL's python3 on Windows)."""
    if pysam_where() == "wsl":
        return _wsl_shell('exec python3 "$0" "$@"', wsl_path(script), *[wsl_path(a) for a in args])
    return [sys.executable, str(script), *[str(a) for a in args]]


def check():
    """Status of every optional dependency, with a plain message for each."""
    out = []
    uses = {"samtools": "counting reads of BAM files (samtools mpileup, the fastest route) and indexing BAM files",
            "bcftools": "indexing and checking the filtered VCF that is written (bcftools index / stats)"}
    for name in ("samtools", "bcftools"):
        w = where(name)
        ver = version(name) if w else ""
        out.append({"name": name, "found": bool(w), "where": {"path": "on PATH", "wsl": "inside WSL"}.get(w, ""),
                    "version": ver, "needed_for": uses[name],
                    "message": (f"{ver or name} found ({'inside WSL' if w == 'wsl' else 'on PATH'}): {uses[name]}." if w else
                                f"{name} was not found. " + (
                                    "The built-in BAM reader (or pysam) is used instead; samtools is faster for whole-genome BAM files. "
                                    "Install it with: conda install -c bioconda samtools, or apt install samtools." if name == "samtools" else
                                    "The filtered VCF is still written, but it is not indexed or checked. "
                                    "Install it with: conda install -c bioconda bcftools, or apt install bcftools."))})
    pw = pysam_where()
    out.append({"name": "pysam", "found": bool(pw), "where": {"python": "python", "wsl": "inside WSL"}.get(pw, ""),
                "version": pysam_version() if pw else "",
                "needed_for": "counting reads of BAM files (alternative to samtools mpileup)",
                "message": (f"pysam {pysam_version()} found ({'inside WSL' if pw == 'wsl' else 'python'}): an alternative to "
                            "samtools for counting reads of BAM files." if pw else
                            "pysam is not installed. It is optional: the built-in reader or samtools count the reads instead. "
                            "Install it with: pip install pysam (Linux / macOS), or inside WSL on Windows.")})
    return out
