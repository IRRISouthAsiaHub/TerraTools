import sys
from pathlib import Path

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QColor, QIcon
from PySide6.QtWidgets import QFileDialog, QMainWindow
from PySide6.QtWebEngineCore import QWebEnginePage, QWebEngineSettings
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWebChannel import QWebChannel

from ui.bridge import TerraToolsBridge
from ui.frontend import load_frontend
from ui.tool_registry import register_tool_bridges


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

        # -------------------------------------------------
        # Python ↔ JavaScript bridges: the launcher plus every tool,
        # all on one channel, so tools open inside this window.
        # -------------------------------------------------

        self.channel = QWebChannel()
        self.bridge = TerraToolsBridge(self)
        self.channel.registerObject("terraTools", self.bridge)
        self.tool_bridges = register_tool_bridges(self.channel, self)
        self.browser.page().setWebChannel(self.channel)

        load_frontend(self.browser)

    def _on_download(self, download):
        path, _ = QFileDialog.getSaveFileName(self, "Save File", download.downloadFileName())
        if not path:
            download.cancel()
            return
        p = Path(path)
        download.setDownloadDirectory(str(p.parent))
        download.setDownloadFileName(p.name)
        download.accept()
