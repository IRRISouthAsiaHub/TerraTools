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
    """creationflags for every analysis subprocess: no console window, and BELOW_NORMAL priority.

    The second flag matters on small machines: R / TensorFlow / scikit-learn workers use every core (n_jobs=-1, BLAS
    threads), and at normal priority they starve the app's own UI thread and the Chromium processes, so the whole
    window stops repainting while an analysis runs. Below-normal priority costs the analysis nothing when the machine
    is otherwise idle and keeps the interface responsive when it is not."""
    if sys.platform.startswith("win"):
        return {"creationflags": subprocess.CREATE_NO_WINDOW | subprocess.BELOW_NORMAL_PRIORITY_CLASS}
    return {}
