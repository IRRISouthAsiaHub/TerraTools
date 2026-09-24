"""
paths.py — Where tools write their result files.

Windows installer / running from source: results stay next to the tool
(tools/<tool>/ui/results), exactly as before — the Windows installer grants
users write access to the install folder for this.

Packaged macOS app: the .app bundle is read-only and code-signed, so
writing inside it fails (or breaks the signature). Results go to
~/Library/Application Support/TerraTools/results/<tool> instead.
"""

import sys
from pathlib import Path


def results_dir(tool: str, default: Path) -> Path:
    if sys.platform == "darwin" and getattr(sys, "frozen", False):
        path = Path.home() / "Library" / "Application Support" / "TerraTools" / "results" / tool
    else:
        path = Path(default)
    path.mkdir(parents=True, exist_ok=True)
    return path
