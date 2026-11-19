import os
import sys
import time
from pathlib import Path

from PySide6.QtCore import QThread, QTimer, Qt, QUrl, Signal
from PySide6.QtGui import QColor, QIcon
from PySide6.QtWidgets import QFileDialog, QMainWindow
from PySide6.QtWebEngineCore import QWebEnginePage, QWebEngineSettings
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWebChannel import QWebChannel

from ui.bridge import TerraToolsBridge
from ui.file_bridge import FileBridge
from ui.frontend import load_frontend
from ui.tool_registry import TOOL_BRIDGES, register_tool_bridges
from tools.common import diagnostics as diag


# Shown the moment the window opens, while the tool modules load (pandas, scipy, sklearn ... are large; on a freshly
# installed machine with antivirus scanning every DLL this takes tens of seconds). Plain CSS, no script.
_SPLASH = """<!doctype html><meta charset="utf-8"><body style="margin:0;height:100vh;display:grid;place-items:center;
background:#0f2419;color:#d8f3dc;font:500 15px 'Segoe UI',sans-serif"><div style="text-align:center">
<div style="font:700 34px Georgia,serif;color:#fff;margin-bottom:6px">TerraTools</div>
<div style="opacity:.75">Loading the analysis tools&hellip;</div>
<div style="margin:22px auto 0;width:180px;height:4px;border-radius:4px;background:rgba(255,255,255,.12);overflow:hidden">
<div style="width:40%;height:100%;background:#52b788;animation:m 1.2s ease-in-out infinite alternate"></div></div></div>
<style>@keyframes m{from{transform:translateX(-30%)}to{transform:translateX(200%)}}</style>"""


class _Preload(QThread):
    """Imports every tool module off the GUI thread, so the window stays responsive while pandas / scipy / sklearn load.
    Creating the bridge objects (QObjects) stays on the GUI thread afterwards and is fast."""

    done = Signal()

    def run(self):
        for load, name in TOOL_BRIDGES:
            try:
                with diag.timed("import_tool", tool=name):
                    load()
            except Exception:
                pass          # reported again, with a traceback, by register_tool_bridges
        try:
            from tools.assistant.ui import assistant_window  # noqa: F401
        except Exception:
            pass
        self.done.emit()


class AppPage(QWebEnginePage):
    """Keeps the top-level page on the TerraTools frontend.

    Without this, the view behaves like a browser tab: dropping a file onto
    the window or following a stray link replaces the whole app with that
    file/page. Embedded frames (interactive plotly/HTML reports) are still
    allowed to load anything."""

    def acceptNavigationRequest(self, url, nav_type, is_main_frame):
        if not is_main_frame:
            return True
        current = self.url()
        if current.isEmpty() or current.scheme() == "data":
            return True
        # Same document (hash routing) is fine; anything else is refused.
        return url.scheme() == current.scheme() and url.path() == current.path()


class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()

        self.setWindowTitle("TerraTools")
        self.setMinimumSize(1100, 720)
        self.resize(1440, 900)
        icon = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent.parent)) / "terratools.ico"
        if icon.is_file():
            self.setWindowIcon(QIcon(str(icon)))

        # -------------------------------------------------
        # Browser — configured to behave like an app surface
        # -------------------------------------------------

        self.browser = QWebEngineView()
        self.browser.setPage(AppPage(self.browser))
        self.browser.setContextMenuPolicy(Qt.NoContextMenu)   # no Back/Reload/View source menu
        self.browser.page().setBackgroundColor(QColor("#0f2419"))  # no white flash on start
        self.setCentralWidget(self.browser)

        settings = self.browser.settings()
        # Result plots/tables live under tools/*/results, outside the frontend folder.
        settings.setAttribute(QWebEngineSettings.LocalContentCanAccessFileUrls, True)
        settings.setAttribute(QWebEngineSettings.LocalContentCanAccessRemoteUrls, False)
        settings.setAttribute(QWebEngineSettings.ScrollAnimatorEnabled, True)

        # "Save as PNG" style downloads triggered from the page.
        self.browser.page().profile().downloadRequested.connect(self._on_download)

        # A crashed or killed Chromium renderer leaves a blank white view: reload the app instead.
        self.browser.page().renderProcessTerminated.connect(self._on_renderer_gone)
        self._renderer_restarts = 0

        self.channel = QWebChannel()
        self.bridge = TerraToolsBridge(self)
        self.channel.registerObject("terraTools", self.bridge)
        self.file_bridge = FileBridge(self)           # image downloads, shared by every tool
        self.channel.registerObject("terraFiles", self.file_bridge)
        self.tool_bridges = {}
        self.browser.page().setWebChannel(self.channel)

        # The window appears at once with a small splash; tool modules load on a worker thread, then the bridges are
        # registered (all on one channel, so tools open inside this window) and the React app is loaded.
        self.browser.setHtml(_SPLASH)
        self._t_open = time.perf_counter()
        self._preload = _Preload(self)
        self._preload.done.connect(self._finish_startup, Qt.QueuedConnection)
        QTimer.singleShot(0, self._preload.start)

    def _finish_startup(self):
        with diag.timed("register_bridges"):
            self.tool_bridges = register_tool_bridges(self.channel, self)
        for name, bridge in self.tool_bridges.items():
            diag.observe_bridge(name, bridge)
        self.browser.loadFinished.connect(self._after_load)
        load_frontend(self.browser)
        diag.event("ui_loading", since_window_s=round(time.perf_counter() - self._t_open, 2))

    def _after_load(self, ok):
        """Tell the page which graphics mode is active (software rendering starts with the decorative animations off)."""
        if not ok or self.browser.url().scheme() != "file":
            return
        software = "--disable-gpu" in os.environ.get("QTWEBENGINE_CHROMIUM_FLAGS", "")
        if software:
            self.browser.page().runJavaScript("document.documentElement.dataset.gfx='software';document.documentElement.dataset.motion='reduced'")
        diag.event("ui_loaded", since_window_s=round(time.perf_counter() - self._t_open, 2), graphics="software" if software else "hardware")

    def _on_renderer_gone(self, status, code):
        diag.event("renderer_terminated", status=str(status), code=code, restarts=self._renderer_restarts)
        if self._renderer_restarts < 3:
            self._renderer_restarts += 1
            QTimer.singleShot(300, lambda: load_frontend(self.browser))

    def closeEvent(self, event):
        """Ask every running tool to stop, give them a moment, then quit. The job object (tools.common.win_job)
        removes any process that is still running."""
        for name, bridge in self.tool_bridges.items():
            for slot in ("cancel", "cancel_run", "stop", "stop_matrices"):
                fn = getattr(bridge, slot, None)
                if callable(fn):
                    try:
                        fn()
                    except Exception:
                        pass
                    break
        super().closeEvent(event)

    def _on_download(self, download):
        path, _ = QFileDialog.getSaveFileName(self, "Save File", download.downloadFileName())
        if not path:
            download.cancel()
            return
        p = Path(path)
        download.setDownloadDirectory(str(p.parent))
        download.setDownloadFileName(p.name)
        download.accept()
