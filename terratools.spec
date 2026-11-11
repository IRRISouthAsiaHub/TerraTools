# terratools.spec — PyInstaller build spec for TerraTools.
#
# Build with:
#   pyinstaller terratools.spec
#
# Output goes to dist/TerraTools/ (a folder build, NOT --onefile).
# A folder build is strongly recommended over --onefile here: --onefile
# extracts everything to a temp dir on every launch, which is painfully
# slow with TensorFlow + Qt WebEngine's Chromium data files, and makes it
# awkward to place the bundled r_runtime folder next to the exe. The
# Inno Setup script below just installs this whole dist/TerraTools/
# folder as-is.
#
# IMPORTANT — run this from your project root (where main.py lives).
# Adjust the paths below if your layout differs.

# -*- mode: python ; coding: utf-8 -*-

import sys
from pathlib import Path

block_cipher = None

PROJECT_ROOT = Path(SPECPATH)  # SPECPATH is injected by PyInstaller

# ── React frontend: rebuild frontend/dist if it's missing or stale ────────
# The exe only ships the built frontend/dist/index.html (Node is NOT needed
# on end-user machines), but it must be built on the packaging machine.
# This rebuilds automatically whenever any file under frontend/src (or the
# frontend config) is newer than the last build.
import shutil
import subprocess

FRONTEND = PROJECT_ROOT / "frontend"
FRONTEND_INDEX = FRONTEND / "dist" / "index.html"
_sources = [p for p in (FRONTEND / "src").rglob("*") if p.is_file()]
_sources += [FRONTEND / n for n in ("index.html", "package.json", "vite.config.js")]
_stale = not FRONTEND_INDEX.exists() or max(p.stat().st_mtime for p in _sources) > FRONTEND_INDEX.stat().st_mtime

if _stale:
    npm = shutil.which("npm")
    if not npm:
        raise SystemExit(
            "frontend/dist is missing or out of date and npm was not found.\n"
            "Install Node.js (https://nodejs.org), then run:\n"
            "    cd frontend && npm install && npm run build"
        )
    print("[terratools.spec] Building React frontend (npm run build)...")
    if not (FRONTEND / "node_modules").exists():
        subprocess.run([npm, "ci"], cwd=FRONTEND, check=True)
    subprocess.run([npm, "run", "build"], cwd=FRONTEND, check=True)
else:
    print("[terratools.spec] React frontend is up to date.")

# ── Data files: the React frontend and tool assets ────────────────────────
def collect_tree(src, dest):
    """Every file under `src` as a PyInstaller data entry, minus things that
    must not ship: results/ folders (output from running the tools on the
    developer's own data — genotype matrices, hybrid results, even BAM files —
    which the tools recreate when they run), Python caches and R history."""
    skip_dirs = {"results", "__pycache__"}
    entries = []
    for f in sorted(Path(src).rglob("*")):
        rel = f.relative_to(src)
        if not f.is_file() or skip_dirs & set(rel.parts[:-1]):
            continue
        if f.suffix in (".pyc", ".pyo") or f.name == ".Rhistory":
            continue
        entries.append((str(f), str(Path(dest) / rel.parent)))
    return entries


datas = [
    # React frontend — a single self-contained index.html used by every
    # window. Run `npm run build` in frontend/ before packaging.
    (str(PROJECT_ROOT / "frontend" / "dist"), "frontend/dist"),
    *collect_tree(PROJECT_ROOT / "tools", "tools"),   # backends, R scripts, configs (and OS_model.dat if bundled)
    *collect_tree(PROJECT_ROOT / "ui", "ui"),
    (str(PROJECT_ROOT / "terratools.ico"), "."),   # main window icon
]

# ── Hidden imports PyInstaller's static analysis tends to miss ────────────
hidden_imports = [
    # TensorFlow pulls in a lot dynamically
    "tensorflow",
    "tensorflow.keras",
    "tensorflow.python.keras.engine.base_layer_v1",
    # sklearn/scipy dynamic submodules
    "sklearn.utils._typedefs",
    "sklearn.neighbors._partition_nodes",
    "scipy.special.cython_special",
    # imbalanced-learn (gene prioritization)
    "imblearn.over_sampling",
    # Genotype Viewer (GenoSee): Biopython is imported only when a species is added from a FASTA file
    "Bio.SeqIO",
    "tqdm",
    "vcf",
    # PySide6 WebEngine
    "PySide6.QtWebEngineWidgets",
    "PySide6.QtWebEngineCore",
    "PySide6.QtWebChannel",
    # numpy submodules referenced only via pickle's dynamic find_class()
    # lookup (e.g. joblib.load() on a pickled DataFrame/array) — these
    # are invisible to PyInstaller's static import scanning since
    # nothing in the source code does `import numpy.core.multiarray`
    # directly. Missing these produces:
    #   ModuleNotFoundError: No module named 'numpy.core.multiarray'
    # Covers both numpy <2.0 (numpy.core.*) and >=2.0 (numpy._core.*,
    # with numpy.core kept only as a pickle-compatibility shim) naming,
    # since the pickling environment's numpy version may not match the
    # build environment's.
    "numpy.core.multiarray",
    "numpy.core._multiarray_umath",
    "numpy.core.multiarray_umath",
    "numpy._core.multiarray",
    "numpy._core._multiarray_umath",
]

# ── Collect entire packages whose static analysis is unreliable ───────────
# (PyInstaller's collect_all pulls in all submodules + data files for these)
from PyInstaller.utils.hooks import collect_all

collect_names = ["numpy", "tensorflow", "sklearn", "scipy", "matplotlib", "plotly", "PySide6"]
binaries = []
for name in collect_names:
    d, b, h = collect_all(name)
    datas += d
    binaries += b
    hidden_imports += h

a = Analysis(
    ["main.py"],  # your app's entrypoint — adjust if it's named differently
    pathex=[str(PROJECT_ROOT)],
    binaries=binaries,
    datas=datas,
    hiddenimports=hidden_imports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="TerraTools",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,          # UPX compression + TensorFlow's DLLs are a known
                         # source of false-positive antivirus flags — leave off
    console=False,       # no terminal window; set True temporarily if you
                          # need to debug a launch crash
    icon='terratools.ico',            # point at a .ico file here if you have one
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name="TerraTools",
)
