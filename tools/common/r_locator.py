"""
r_locator.py — Shared helper for finding the Rscript executable that
TerraTools' R-backed tools (Mendelian, QTL, Phenotypic, OCS) should use.

Lookup order:
  1. TERRATOOLS_RSCRIPT env var, if set (explicit override)
  2. A bundled R runtime shipped inside the frozen .exe build
     (<exe_dir>/r_runtime/bin/Rscript.exe) — this is what end users hit
  3. The project-local conda env created by setup_r_env.sh, for
     developers running from source (<project_root>/r_env/bin/Rscript)
  4. Rscript on the system PATH (fallback)

Returns None if nothing is found, so callers can show a clear,
actionable error instead of hanging or crashing.
"""

import os
import shutil
import sys
from pathlib import Path


def _is_frozen() -> bool:
    """True when running inside a PyInstaller-built .exe."""
    return getattr(sys, "frozen", False)


def _app_root() -> Path:
    """Directory containing the app's own resources.
    - Frozen build: the folder the .exe lives in (or _MEIPASS for --onefile
      temp extraction; PyInstaller sets both, sys.executable is preferred
      here since bundled *data* like r_runtime should be placed next to
      the exe by the Inno Setup installer, not inside the onefile temp
      extraction folder).
    - Running from source: tools/common -> tools -> project root.
    """
    if _is_frozen():
        return Path(sys.executable).resolve().parent
    return Path(__file__).resolve().parents[2]


def find_rscript() -> str | None:
    override = os.environ.get("TERRATOOLS_RSCRIPT")
    if override and Path(override).exists():
        return override

    root = _app_root()
    is_windows = sys.platform.startswith("win")

    candidates = []
    if _is_frozen():
        # Bundled runtime shipped by the installer, e.g.:
        #   C:\Program Files\TerraTools\r_runtime\bin\Rscript.exe
        if is_windows:
            candidates.append(root / "r_runtime" / "bin" / "Rscript.exe")
        else:
            candidates.append(root / "r_runtime" / "bin" / "Rscript")
    else:
        # Developer setup via setup_r_env.sh
        if is_windows:
            candidates.append(root / "r_env" / "Scripts" / "Rscript.exe")
        else:
            candidates.append(root / "r_env" / "bin" / "Rscript")

    for candidate in candidates:
        if candidate.exists():
            return str(candidate)

    system_r = shutil.which("Rscript")
    if system_r:
        return system_r

    return None


def rscript_not_found_message() -> str:
    if _is_frozen():
        return (
            "Rscript was not found.\n\n"
            "This looks like a corrupted or incomplete installation — the "
            "bundled R runtime should be included automatically. Try "
            "reinstalling TerraTools, or contact support if the problem "
            "persists."
        )
    return (
        "Rscript was not found.\n\n"
        "TerraTools looks for R in this order:\n"
        "  1. A project-local R environment at ./r_env (recommended)\n"
        "  2. R installed system-wide on your PATH\n\n"
        "To set up the project-local environment (isolated, no admin/sudo "
        "needed), run this once from the project root:\n"
        "  bash setup_r_env.sh\n\n"
        "This requires conda or mamba to be installed first — see "
        "setup_r_env.sh for a one-time Miniconda install link if you don't "
        "have it yet."
    )


def get_r_env(rscript) -> dict:
    """
    Build a subprocess environment for the R runtime, IF R_HOME can be
    confidently resolved — and leave the environment untouched otherwise.

    This matters specifically for a portable/copied R install (like the
    r_runtime folder the installer ships): a normal Windows R install
    relies on the registry or its install path to auto-detect R_HOME,
    but a bare copied folder run from an arbitrary location may not
    auto-detect correctly without R_HOME set explicitly.

    Only sets R_HOME when real evidence is found (a library/ or
    etc/Renviron marker) within a bounded number of parent levels of the
    Rscript executable — never guessed for a normal system/conda R
    install, to avoid breaking a perfectly good existing R setup.
    """
    import os

    env = os.environ.copy()
    rscript = Path(rscript)

    r_home = None
    candidates = list(rscript.parents)[:4]
    candidates.append(rscript.parent.parent / "lib" / "R")

    for parent in candidates:
        if (parent / "library").exists() or (parent / "etc" / "Renviron").exists():
            r_home = parent
            break

    if r_home is None:
        return env

    env["R_HOME"] = str(r_home)
    env["PATH"] = str(rscript.parent) + os.pathsep + env.get("PATH", "")
    return env
