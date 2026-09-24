from PySide6.QtCore import QObject, Slot

from ui.tool_manager import ToolManager


class TerraToolsBridge(QObject):

    def __init__(self, parent_window=None):

        super().__init__(parent_window)

        self.parent_window = parent_window

        self.tool_manager = ToolManager(
            parent_window
        )

    @Slot(str)
    def open_tool(self, tool_name):

        print(
            f"[TerraTools] Requested tool: "
            f"{tool_name!r}"
        )

        try:

            # Normalize the value coming from HTML.
            tool_name = str(
                tool_name
            ).strip().lower()

            print(
                f"[TerraTools] Normalized tool: "
                f"{tool_name!r}"
            )

            self.tool_manager.open_tool(
                tool_name
            )

        except Exception as e:

            import traceback

            print(
                "\n========== TOOL ERROR =========="
            )

            traceback.print_exc()

            print(
                "================================\n"
            )