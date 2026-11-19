import json

from PySide6.QtCore import QObject, Slot

from tools.common import diagnostics as diag


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

    @Slot(str)
    def diag_log(self, payload):
        """Frontend reports (long tasks, slow frames, result sizes, render times) go to the same structured log as the
        Python side. Cheap: the page batches them and sends at most one call every few seconds."""
        try:
            for rec in json.loads(payload):
                diag.event("fe_" + str(rec.pop("ev", "event")), **rec)
        except Exception:
            pass

    @Slot(str)
    def set_active_tool(self, tool):
        diag.set_tool(tool or None)
