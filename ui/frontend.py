import sys
from pathlib import Path

from PySide6.QtCore import QUrl


def _dist_index() -> Path:
    # PyInstaller folder builds unpack data files under sys._MEIPASS.
    base = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent.parent))
    return base / "frontend" / "dist" / "index.html"


_MISSING_BUILD_HTML = """
<html><body style="font-family:sans-serif;background:#f5faf7;color:#1a3a2a;padding:48px">
<h2>TerraTools frontend has not been built</h2>
<p>Expected <code>{path}</code>.</p>
<p>Build it once with:</p>
<pre style="background:#1a3a2a;color:#d8f3dc;padding:14px;border-radius:8px">cd frontend
npm install
npm run build</pre>
</body></html>
"""


def load_frontend(browser, route: str = "") -> None:
    """Load the React frontend (frontend/dist/index.html) into a
    QWebEngineView.

    route "" is the launcher dashboard; any other value selects a tool page
    and matches the ids used by ToolManager (see frontend/src/App.jsx)."""
    index = _dist_index()
    if not index.is_file():
        print(f"[TerraTools] React build not found: {index}", flush=True)
        browser.setHtml(_MISSING_BUILD_HTML.format(path=index))
        return
    url = QUrl.fromLocalFile(str(index))
    if route:
        url.setFragment(f"/{route}")
    browser.load(url)
