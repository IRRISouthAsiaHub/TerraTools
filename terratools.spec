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

# ── Data files: everything the HTML frontends and tool assets need ────────
datas = [
    (str(PROJECT_ROOT / "dashboard"), "dashboard"),
    # React frontend (run `npm run build` in frontend/ before packaging).
    (str(PROJECT_ROOT / "frontend" / "dist"), "frontend/dist"),
    (str(PROJECT_ROOT / "tools"), "tools"),   # includes every .html, config, and (if you bundled it) OS_model.dat
    (str(PROJECT_ROOT / "ui"), "ui"),
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
