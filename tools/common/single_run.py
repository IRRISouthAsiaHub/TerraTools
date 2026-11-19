"""single_run: one analysis per tool at a time.

Every tool bridge keeps its worker in self._thread / self._worker and replaces them on each run. A second run request
while the first is still computing (a double click, a shortcut, a page that calls twice) dropped the last Python
reference to the running QThread, and Qt aborts the whole process with "QThread: Destroyed while thread is still
running" - the window just vanished. This decorator refuses the second request instead and tells the user.

    @Slot(str)
    @single_run
    def run_something(self, config_json): ...
"""

import functools

try:
    from shiboken6 import isValid
except Exception:                       # pragma: no cover
    def isValid(_obj):
        return True

from tools.common import diagnostics as diag


def thread_busy(thread) -> bool:
    """True while `thread` exists and runs. A QThread that finished and was deleteLater()-ed counts as not busy."""
    try:
        return thread is not None and isValid(thread) and thread.isRunning()
    except RuntimeError:
        return False


def single_run(fn):
    @functools.wraps(fn)
    def wrapper(self, *args, **kwargs):
        if thread_busy(getattr(self, "_thread", None)):
            diag.event("run_rejected", tool=type(self).__name__, slot=fn.__name__)
            try:
                self.statusUpdate.emit("A run of this tool is already in progress. Wait for it to finish.")
            except Exception:
                pass
            return None
        return fn(self, *args, **kwargs)
    return wrapper
