#!/usr/bin/env bash
# build_macos.sh — builds TerraTools.app and TerraTools.dmg on macOS.
#
# The macOS counterpart of `pyinstaller terratools.spec` + Inno Setup.
# Must run ON a Mac (PyInstaller cannot cross-compile from Windows); the
# result runs on Macs with the same CPU type as the build Mac (Apple
# Silicon build -> Apple Silicon Macs).
#
#   bash build_macos/build_macos.sh
#
# Requirements on the build Mac:
#   - conda (Miniforge recommended: https://github.com/conda-forge/miniforge)
#   - Node.js / npm (https://nodejs.org)
#   - codesign and hdiutil ship with macOS
#
# Output:
#   dist/TerraTools.app
#   installer_output/TerraTools.dmg   (drag-to-Applications disk image)

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

WORK="$ROOT/build_macos/work"
PY_ENV="$WORK/py_env"
PYTHON_VERSION="3.12"
# The R runtime is built at its final install path (see
# tools/common/r_locator.py). MUST match MACOS_R_PREFIX there.
R_PREFIX="/Users/Shared/TerraTools/r_runtime"

# Every R package the tools' R scripts load.
R_PACKAGES='c("agricolae","broom","car","caret","cluster","corrplot","dplyr","emmeans","factoextra","GGally","ggplot2","ggrepel","gridExtra","htmlwidgets","jsonlite","lme4","moments","openxlsx","pheatmap","plotly","purrr","RColorBrewer","readr","readxl","reshape2","rrBLUP","scales","sommer","stringr","tidyr","tidyverse","viridisLite","writexl")'

step() { printf '\n\033[1;32m==> %s\033[0m\n' "$*"; }
die()  { printf '\033[1;31mError: %s\033[0m\n' "$*" >&2; exit 1; }

[[ "$(uname)" == "Darwin" ]] || die "This script must run on macOS."
command -v conda >/dev/null || die "conda not found. Install Miniforge: https://github.com/conda-forge/miniforge"
command -v npm   >/dev/null || die "npm not found. Install Node.js: https://nodejs.org"
mkdir -p "$WORK"

# ---------------------------------------------------------------------------
step "1/6  React frontend"
(cd frontend && npm ci && npm run build)

# ---------------------------------------------------------------------------
step "2/6  Python build environment"
if [[ ! -x "$PY_ENV/bin/python" ]]; then
  conda create -y -p "$PY_ENV" -c conda-forge "python=$PYTHON_VERSION"
fi
PY="$PY_ENV/bin/python"
"$PY" -m pip install --upgrade pip
"$PY" -m pip install -r requirements.txt pyinstaller pillow   # pillow: .ico -> .icns

# ---------------------------------------------------------------------------
step "3/6  R runtime at $R_PREFIX"
mkdir -p "$(dirname "$R_PREFIX")"
if [[ -d "$R_PREFIX/conda-meta" ]]; then
  conda env update -p "$R_PREFIX" -f environment.yml --prune
else
  rm -rf "$R_PREFIX"
  conda env create -p "$R_PREFIX" -f environment.yml
fi
"$R_PREFIX/bin/Rscript" -e "
  pk <- $R_PACKAGES
  missing <- pk[!vapply(pk, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) { message('Missing R packages: ', paste(missing, collapse = ', ')); quit(status = 1) }
  cat('R', as.character(getRversion()), '- all', length(pk), 'packages load\n')
"

# ---------------------------------------------------------------------------
step "4/6  Pack R runtime"
rm -f "$WORK/r_runtime.tar.gz"
COPYFILE_DISABLE=1 tar --no-mac-metadata -czf "$WORK/r_runtime.tar.gz" \
  --exclude ./.terratools-stamp -C "$R_PREFIX" .
shasum -a 256 "$WORK/r_runtime.tar.gz" | cut -c1-16 > "$WORK/r_runtime.stamp"
# Mark this Mac's copy as current so running the app here doesn't re-unpack it.
cp "$WORK/r_runtime.stamp" "$R_PREFIX/.terratools-stamp"
du -sh "$WORK/r_runtime.tar.gz"

# ---------------------------------------------------------------------------
step "5/6  PyInstaller -> dist/TerraTools.app"
rm -rf "$ROOT/dist/TerraTools" "$ROOT/dist/TerraTools.app"
"$PY" -m PyInstaller build_macos/terratools_macos.spec --noconfirm \
  --distpath "$ROOT/dist" --workpath "$WORK/pyinstaller"
APP="$ROOT/dist/TerraTools.app"
[[ -d "$APP" ]] || die "PyInstaller did not produce $APP"
codesign --verify --deep --strict "$APP" && echo "Signature OK (ad-hoc)" \
  || echo "WARNING: codesign verification failed — the app may be reported as damaged."

# ---------------------------------------------------------------------------
step "6/6  Disk image -> installer_output/TerraTools.dmg"
STAGE="$WORK/dmg"
rm -rf "$STAGE"
mkdir -p "$STAGE" "$ROOT/installer_output"
ditto "$APP" "$STAGE/TerraTools.app"
ln -s /Applications "$STAGE/Applications"
rm -f "$ROOT/installer_output/TerraTools.dmg"
hdiutil create -volname "TerraTools" -srcfolder "$STAGE" -ov -format UDZO \
  "$ROOT/installer_output/TerraTools.dmg"
rm -rf "$STAGE"

step "Done"
echo "App: $APP"
echo "DMG: $ROOT/installer_output/TerraTools.dmg"
