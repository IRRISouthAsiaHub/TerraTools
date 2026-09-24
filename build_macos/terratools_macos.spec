# terratools_macos.spec — PyInstaller spec for the macOS app (TerraTools.app).
#
# Don't run this directly — use build_macos/build_macos.sh, which first
# builds the React frontend and the packed R runtime this spec bundles.
# (The Windows build uses terratools.spec in the project root.)

# -*- mode: python ; coding: utf-8 -*-

from pathlib import Path

from PyInstaller.utils.hooks import collect_all

block_cipher = None

PROJECT_ROOT = Path(SPECPATH).parent          # SPECPATH = build_macos/
WORK = PROJECT_ROOT / "build_macos" / "work"  # created by build_macos.sh

for required in (
    PROJECT_ROOT / "frontend" / "dist" / "index.html",
    WORK / "r_runtime.tar.gz",
    WORK / "r_runtime.stamp",
):
    if not required.exists():
        raise SystemExit(f"Missing {required} — run build_macos/build_macos.sh instead of this spec directly.")

datas = [
    (str(PROJECT_ROOT / "frontend" / "dist"), "frontend/dist"),
    (str(PROJECT_ROOT / "tools"), "tools"),
    (str(PROJECT_ROOT / "ui"), "ui"),
    # Packed R runtime, unpacked on first use by tools/common/r_locator.py
    (str(WORK / "r_runtime.tar.gz"), "."),
    (str(WORK / "r_runtime.stamp"), "."),
]

hidden_imports = [
    "tensorflow",
    "tensorflow.keras",
    "sklearn.utils._typedefs",
    "sklearn.neighbors._partition_nodes",
    "scipy.special.cython_special",
    "imblearn.over_sampling",
    "PySide6.QtWebEngineWidgets",
    "PySide6.QtWebEngineCore",
    "PySide6.QtWebChannel",
    "numpy.core.multiarray",
    "numpy.core._multiarray_umath",
    "numpy._core.multiarray",
    "numpy._core._multiarray_umath",
    "tarfile",   # R runtime unpacking
]

# PySide6 is intentionally NOT collected wholesale here: on macOS that drags
# in every Qt framework as loose data and breaks code signing. PyInstaller's
# own PySide6 hooks bundle exactly the Qt modules the app imports (including
# QtWebEngine and its helper app).
binaries = []
for name in ["numpy", "tensorflow", "sklearn", "scipy", "matplotlib", "plotly"]:
    d, b, h = collect_all(name)
    datas += d
    binaries += b
    hidden_imports += h

a = Analysis(
    [str(PROJECT_ROOT / "main.py")],
    pathex=[str(PROJECT_ROOT)],
    binaries=binaries,
    datas=datas,
    hiddenimports=hidden_imports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
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
    upx=False,
    console=False,
    argv_emulation=False,
    target_arch=None,        # build machine's architecture (arm64 on Apple Silicon)
    codesign_identity=None,  # ad-hoc signed; see build_macos.sh for Developer ID signing
    entitlements_file=None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    name="TerraTools",
)

app = BUNDLE(
    coll,
    name="TerraTools.app",
    icon=str(PROJECT_ROOT / "terratools.ico"),   # converted to .icns automatically (needs Pillow)
    bundle_identifier="org.irri.terratools",
    info_plist={
        "CFBundleName": "TerraTools",
        "CFBundleDisplayName": "TerraTools",
        "CFBundleShortVersionString": "1.0.0",
        "CFBundleVersion": "1.0.0",
        "NSHighResolutionCapable": True,
        "LSMinimumSystemVersion": "12.0",
        "NSRequiresAquaSystemAppearance": False,
    },
)
