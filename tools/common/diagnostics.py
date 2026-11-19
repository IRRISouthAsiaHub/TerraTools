"""Lightweight, always-safe diagnostics for TerraTools (development and packaged Windows builds).

Writes JSON lines to  %LOCALAPPDATA%\\TerraTools\\logs\\diagnostics.log  (rotated at 1 MB, 3 files kept).
Cost is a few microseconds per event and one 250 ms timer; nothing is written unless something is worth noting:

  * startup phases (imports, bridge registration, first page load) with durations
  * GUI-thread stalls: a timer on the GUI thread measures its own lateness; when the GUI thread has not run
    for more than STALL_MS a watchdog thread records the Python stack of every thread, which names the blocking call
  * every tool result that crosses the bridge: tool, status, payload size, time since the run was requested
  * frontend reports (long tasks, slow frames, slow renders, JSON parse time) sent through TerraToolsBridge.diag_log
  * the graphics mode chosen at start (see graphics_mode) and why

Set TERRATOOLS_DIAG=verbose to also log every status update and a CPU / memory sample every 5 s.
"""

from __future__ import annotations

import json
import logging
import logging.handlers
import os
import sys
import threading
import time
import traceback
from contextlib import contextmanager
from pathlib import Path

STALL_MS = 400          # a GUI thread that does not run for this long is a stall worth logging
HEARTBEAT_MS = 100      # GUI-thread timer period
WATCHDOG_S = 2.0        # the GUI thread silent this long -> dump all thread stacks (once per stall)

_log = logging.getLogger("terratools.diag")
_t0 = time.perf_counter()
_ready = False
_verbose = os.environ.get("TERRATOOLS_DIAG", "").lower() == "verbose"
_context = {"tool": None}       # the tool page the user is on / last started, for stall reports


def log_dir() -> Path:
    base = os.environ.get("LOCALAPPDATA") or os.path.expanduser("~")
    p = Path(base) / "TerraTools" / "logs"
    p.mkdir(parents=True, exist_ok=True)
    return p


def _setup():
    global _ready
    if _ready:
        return
    _ready = True
    try:
        h = logging.handlers.RotatingFileHandler(log_dir() / "diagnostics.log", maxBytes=1_000_000, backupCount=3, encoding="utf-8")
        h.setFormatter(logging.Formatter("%(message)s"))
        _log.addHandler(h)
        _log.setLevel(logging.INFO)
        _log.propagate = False
    except Exception:
        _log.addHandler(logging.NullHandler())


def event(name: str, **fields):
    """One structured log line: {"t": seconds since start, "ev": name, ...}."""
    _setup()
    try:
        rec = {"t": round(time.perf_counter() - _t0, 3), "ev": name}
        rec.update(fields)
        _log.info(json.dumps(rec, default=str, ensure_ascii=False))
    except Exception:
        pass


@contextmanager
def timed(name: str, **fields):
    s = time.perf_counter()
    try:
        yield
    finally:
        event(name, ms=round((time.perf_counter() - s) * 1000, 1), **fields)


def set_tool(tool: str | None):
    _context["tool"] = tool


def all_stacks() -> dict:
    names = {t.ident: t.name for t in threading.enumerate()}
    out = {}
    for tid, frame in sys._current_frames().items():
        out[names.get(tid, str(tid))] = "".join(traceback.format_stack(frame)[-8:])
    return out


class GuiWatchdog:
    """Detects a frozen GUI thread. A QTimer on the GUI thread stamps a heartbeat; a plain Python thread checks the stamp.
    Whatever the GUI thread is stuck in is visible in the stack dump (e.g. a pandas read inside a bridge slot)."""

    def __init__(self):
        from PySide6.QtCore import QTimer
        self._last = time.perf_counter()
        self._expected = self._last
        self._dumped = False
        self._stop = False
        self._main = threading.main_thread().ident
        self._timer = QTimer()
        self._timer.setInterval(HEARTBEAT_MS)
        self._timer.timeout.connect(self._beat)
        self._timer.start()
        self._thread = threading.Thread(target=self._watch, name="tt-watchdog", daemon=True)
        self._thread.start()
        self._proc = None
        self._tick = 0

    def _beat(self):
        now = time.perf_counter()
        late = (now - self._last) * 1000 - HEARTBEAT_MS
        if late > STALL_MS:
            event("gui_stall", ms=round(late), tool=_context["tool"])
        self._last = now
        self._dumped = False
        self._tick += 1
        if _verbose and self._tick % 50 == 0:
            self._sample()

    def _sample(self):
        try:
            import psutil
            if self._proc is None:
                self._proc = psutil.Process()
            kids = self._proc.children(recursive=True)
            event("sample", cpu=self._proc.cpu_percent(None), rss_mb=round(self._proc.memory_info().rss / 1e6),
                  children=[(c.name(), round(c.cpu_percent(None))) for c in kids][:8])
        except Exception:
            pass

    def _watch(self):
        while not self._stop:
            time.sleep(0.5)
            silent = time.perf_counter() - self._last
            if silent > WATCHDOG_S and not self._dumped:
                self._dumped = True
                event("gui_frozen", silent_s=round(silent, 1), tool=_context["tool"], stacks=all_stacks())


_watchdog = None


def start_gui_watchdog():
    """Call once, after the QApplication exists."""
    global _watchdog
    if _watchdog is None:
        try:
            _watchdog = GuiWatchdog()
        except Exception:
            traceback.print_exc()


def observe_bridge(name: str, bridge):
    """Log payload size / duration of every result a tool bridge emits. Uses only the bridge's existing signals."""
    marks = {"start": None}
    try:
        if hasattr(bridge, "statusUpdate"):
            def _status(msg, n=name):
                if marks["start"] is None:
                    marks["start"] = time.perf_counter()
                    _context["tool"] = n
                if _verbose:
                    event("status", tool=n, msg=str(msg)[:120])
            bridge.statusUpdate.connect(_status)
        if hasattr(bridge, "resultReady"):
            def _result(payload, n=name):
                s = marks["start"]
                size = len(payload) if isinstance(payload, str) else -1
                status = None
                if isinstance(payload, str) and size < 200_000:
                    try:
                        status = json.loads(payload).get("status")
                    except Exception:
                        pass
                event("result", tool=n, bytes=size, status=status,
                      since_first_status_s=None if s is None else round(time.perf_counter() - s, 2),
                      big=size > 5_000_000)
                marks["start"] = None
            bridge.resultReady.connect(_result)
    except Exception:
        pass


# ───────────────────────── graphics mode ─────────────────────────

def _state_file() -> Path:
    return log_dir() / "unclean_exits.txt"


def unclean_exits() -> int:
    try:
        return int(_state_file().read_text().strip() or 0)
    except Exception:
        return 0


def mark_start():
    """A run that reaches mark_clean_exit() resets the counter; one that does not leaves it one higher."""
    try:
        _state_file().write_text(str(unclean_exits() + 1))
    except Exception:
        pass


def mark_clean_exit():
    try:
        _state_file().write_text("0")
    except Exception:
        pass


def _is_remote_session() -> bool:
    try:
        import ctypes
        return bool(ctypes.windll.user32.GetSystemMetrics(0x1000))     # SM_REMOTESESSION
    except Exception:
        return False


def _gpu_names() -> list[str]:
    """Display adapters of the machine (EnumDisplayDevices: no subprocess, a few milliseconds)."""
    try:
        import ctypes
        from ctypes import wintypes

        class DISPLAY_DEVICE(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("DeviceName", wintypes.WCHAR * 32), ("DeviceString", wintypes.WCHAR * 128),
                        ("StateFlags", wintypes.DWORD), ("DeviceID", wintypes.WCHAR * 128), ("DeviceKey", wintypes.WCHAR * 128)]

        names, i = [], 0
        while i < 16:
            d = DISPLAY_DEVICE()
            d.cb = ctypes.sizeof(d)
            if not ctypes.windll.user32.EnumDisplayDevicesW(None, i, ctypes.byref(d), 0):
                break
            if d.DeviceString and d.DeviceString not in names:
                names.append(d.DeviceString)
            i += 1
        return names
    except Exception:
        return []


def graphics_mode(argv: list[str]) -> tuple[str, str]:
    """('hardware' | 'software', reason). Decided before QApplication exists.

    software (Chromium's WARP / SwiftShader rasteriser, no GPU process involved) is chosen when
      * the user asks for it:  --safe-graphics   or   TERRATOOLS_GPU=off
      * the session is a remote desktop (RDP exposes no real GPU; hardware compositing flickers there)
      * the only display adapter is Microsoft Basic Render Driver / a virtual adapter (VM, missing driver)
    TERRATOOLS_GPU=on forces hardware. Two consecutive runs that did not end normally (the window frozen or the
    process killed) also switch the next start to software rendering, so a machine whose graphics driver hangs
    the app recovers by itself."""
    env = os.environ.get("TERRATOOLS_GPU", "").lower()
    if "--safe-graphics" in argv or env in ("off", "0", "software"):
        return "software", "requested (--safe-graphics / TERRATOOLS_GPU=off)"
    if env in ("on", "1", "hardware"):
        return "hardware", "forced (TERRATOOLS_GPU=on)"
    if unclean_exits() >= 2:
        return "software", "the last %d runs did not end normally (frozen or crashed); delete %s to retry hardware" % (unclean_exits(), _state_file())
    if sys.platform.startswith("win"):
        if _is_remote_session():
            return "software", "remote desktop session"
        names = [n.lower() for n in _gpu_names()]
        real = [n for n in names if not any(k in n for k in ("basic render", "basic display", "vmware", "virtualbox", "hyper-v", "remote display", "parsec", "citrix", "microsoft remote"))]
        if names and not real:
            return "software", "no hardware display adapter (%s)" % ", ".join(names)
    return "hardware", "default"
