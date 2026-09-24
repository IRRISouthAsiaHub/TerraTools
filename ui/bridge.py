import json

from PySide6.QtCore import QObject, Slot


class TerraToolsBridge(QObject):
    """Launcher-level bridge. Tools now open inside the main window, so this
    only handles in-app navigation requested from Python-side callers."""

    def __init__(self, parent_window=None):
        super().__init__(parent_window)
        self.parent_window = parent_window

    @Slot(str)
    def open_tool(self, tool_name):
        # Kept for compatibility with older callers: navigate in place.
        tool_name = str(tool_name).strip().lower()
        print(f"[TerraTools] open_tool({tool_name!r}) -> in-app navigation", flush=True)
        self.parent_window.browser.page().runJavaScript(
            f"location.replace({json.dumps('#/' + tool_name)});"
        )
