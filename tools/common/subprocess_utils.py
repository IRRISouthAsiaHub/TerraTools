"""
subprocess_utils.py — Shared helper for spawning subprocesses without a
visible console window on Windows.

On Windows, subprocess.run()/Popen() launching a console application
(Rscript.exe, python.exe, etc.) from a windowed (console=False) app opens
a brand-new console window for that child process. That window is the
process's actual controlling window — closing it kills the process, not
just hides output. subprocess.CREATE_NO_WINDOW suppresses that window
entirely while leaving stdout/stderr piping completely unaffected.

This is a no-op on Linux/macOS, where no such per-process console window
exists in the first place.

Usage:
    from tools.common.subprocess_utils import no_window_kwargs
    subprocess.run(cmd, capture_output=True, text=True, **no_window_kwargs())
    subprocess.Popen(cmd, stdout=subprocess.PIPE, **no_window_kwargs())
"""

import sys
import subprocess


def no_window_kwargs() -> dict:
    if sys.platform.startswith("win"):
        return {"creationflags": subprocess.CREATE_NO_WINDOW}
    return {}
