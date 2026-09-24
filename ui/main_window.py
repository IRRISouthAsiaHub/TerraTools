from pathlib import Path

from PySide6.QtCore import QUrl
from PySide6.QtWidgets import QMainWindow
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWebChannel import QWebChannel

from ui.bridge import TerraToolsBridge
from ui.frontend import frontend_url


class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()

        self.setWindowTitle("TerraTools")
        self.resize(1400, 900)

        # -------------------------------------------------
        # Browser
        # -------------------------------------------------

        self.browser = QWebEngineView()
        self.setCentralWidget(self.browser)

        # -------------------------------------------------
        # Python ↔ JavaScript bridge
        # -------------------------------------------------

        self.bridge = TerraToolsBridge(self)

        self.channel = QWebChannel()
        self.channel.registerObject("terraTools", self.bridge)

        self.browser.page().setWebChannel(self.channel)

        # -------------------------------------------------
        # Load dashboard
        # -------------------------------------------------

        project_dir = Path(__file__).resolve().parent.parent
        html_file = project_dir / "dashboard" / "index.html"

        # React frontend if built, otherwise the original static dashboard.
        self.browser.load(
            frontend_url() or QUrl.fromLocalFile(str(html_file))
        )