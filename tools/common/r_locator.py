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
import threading
from pathlib import Path


# ── macOS packaged app: bundled R runtime ────────────────────────────────
# A conda R environment has its install path baked into its files, so it
# only works at the exact path it was built at. build_macos.sh therefore
# builds it at this fixed, world-writable location, and the .app ships it
# as Resources/r_runtime.tar.gz. On first use it is unpacked back to the
# same path — no files are rewritten, so macOS code signatures stay valid.
# MUST match R_PREFIX in build_macos/build_macos.sh.
MACOS_R_PREFIX = Path("/Users/Shared/TerraTools/r_runtime")
_MACOS_R_LOCK = threading.Lock()


def prepare_macos_r() -> Path | None:
    """Unpack the bundled R runtime if it isn't unpacked (or is from an
    older app version). Returns the Rscript path, or None when this isn't a
    packaged macOS build. Safe to call from several threads; the first run
    takes a minute or so, later calls return immediately."""
    if sys.platform != "darwin" or not _is_frozen():
        return None
    bundle = Path(getattr(sys, "_MEIPASS", Path(sys.executable).parent))
    archive = bundle / "r_runtime.tar.gz"
    stamp_src = bundle / "r_runtime.stamp"
    rscript = MACOS_R_PREFIX / "bin" / "Rscript"
    if not archive.is_file():
        return rscript if rscript.exists() else None

    wanted = stamp_src.read_text().strip() if stamp_src.is_file() else str(archive.stat().st_size)
    stamp = MACOS_R_PREFIX / ".terratools-stamp"

    with _MACOS_R_LOCK:
        if rscript.exists() and stamp.is_file() and stamp.read_text().strip() == wanted:
            return rscript

        import tarfile
        print(f"[TerraTools] Unpacking bundled R runtime to {MACOS_R_PREFIX} ...", flush=True)
        staging = MACOS_R_PREFIX.with_name(".r_runtime.unpacking")
        shutil.rmtree(staging, ignore_errors=True)
        staging.mkdir(parents=True)
        with tarfile.open(archive, "r:gz") as tar:
            if sys.version_info >= (3, 12):
                tar.extractall(staging, filter="tar")
            else:
                tar.extractall(staging)
        # Nobody but the unpacking user may modify these files: the runtime is
        # executed by every account that uses the app, so it must not be
        # replaceable by another local user.
        for root, dirs, files in os.walk(staging):
            for name in dirs + files:
                p = os.path.join(root, name)
                if not os.path.islink(p):
                    os.chmod(p, os.stat(p).st_mode & ~0o022)
        os.chmod(staging, os.stat(staging).st_mode & ~0o022)
        # Contents only reference MACOS_R_PREFIX, so moving the finished
        # tree into place is safe (and a half-unpacked tree never is used).
        shutil.rmtree(MACOS_R_PREFIX, ignore_errors=True)
        try:
            staging.rename(MACOS_R_PREFIX)
        except OSError as e:
            shutil.rmtree(staging, ignore_errors=True)
            raise RuntimeError(
                f"Could not install the R runtime at {MACOS_R_PREFIX}: {e}. Another user account on this Mac "
                "probably installed an older copy; ask them to update TerraTools, or remove that folder as an administrator."
            ) from e
        stamp.write_text(wanted)
        # Shared location like /tmp: every account can create its entry, but
        # (sticky bit) nobody can delete or replace another account's files.
        try:
            os.chmod(MACOS_R_PREFIX.parent, 0o1777)
        except OSError:
            pass
        print("[TerraTools] R runtime ready.", flush=True)
        return rscript


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
    if _is_frozen() and sys.platform == "darwin":
        try:
            mac_r = prepare_macos_r()
        except Exception as e:
            print(f"[TerraTools] Could not unpack bundled R runtime: {e}", flush=True)
            mac_r = None
        if mac_r is not None:
            candidates.append(mac_r)
    elif _is_frozen():
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
            # portable R copied into the project (the folder the installer is built from)
            for folder in ("r_runtime_build", "r_runtime"):
                candidates.append(root / folder / "bin" / "x64" / "Rscript.exe")
                candidates.append(root / folder / "bin" / "Rscript.exe")
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
