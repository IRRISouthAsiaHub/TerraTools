import sys
import os

if len(sys.argv) > 2 and sys.argv[1] == "--worker":
    worker_name = sys.argv[2]
    # Rewrite sys.argv so the target CLI's own argparse sees exactly what
    # it expects (its own argv[0] placeholder + its real flags), exactly
    # as if it had been launched standalone.
    sys.argv = [f"{worker_name}_worker"] + sys.argv[3:]

    if worker_name == "genomic_prediction":
        from tools.genomic_prediction.backend.genomic_prediction_cli import main as _worker_main
    elif worker_name == "heterosis":
        from tools.heterosis.backend.heterosis_process import main as _worker_main
    elif worker_name == "pybsaseq":
        # Located via Python's own import system (works whether frozen
        # or not) rather than a hand-built path, then executed in-process
        # via runpy exactly as if it had been run standalone as
        # `python pybsaseq_script.py <args>` from the given directory.
        import runpy
        import importlib
        from pathlib import Path
 
        if len(sys.argv) < 3 or sys.argv[1] != "--cwd":
            print("pybsaseq worker: expected --cwd <path> as first argument", file=sys.stderr)
            sys.exit(1)
 
        work_dir = sys.argv[2]
        script_args = sys.argv[3:]
 
        backend_pkg = importlib.import_module("tools.pybsaseq.backend")
        script_path = str(Path(backend_pkg.__file__).resolve().parent / "pybsaseq_script.py")
 
        os.chdir(work_dir)
        sys.argv = ["pybsaseq_script.py"] + script_args
        runpy.run_path(script_path, run_name="__main__")
        _worker_main = None

    else:
        print(f"Unknown worker: {worker_name}", file=sys.stderr)
        sys.exit(1)

    # The imported entry point has to be CALLED: importing it only defines it, and the process
    # then exited with code 0 without ever running the analysis (BreedSight and Heterosis in a packaged build).
    # The BreedSight CLI exits by itself (0 on success, 1 on error); the Heterosis one returns.
    if _worker_main is not None:
        _worker_main()
    sys.exit(0)

# ── Normal GUI startup ───────────────────────────────────────────────────
from tools.common import diagnostics as diag

diag.event("start", frozen=bool(getattr(sys, "frozen", False)), argv=sys.argv[1:], python=sys.version.split()[0])
diag.mark_start()          # reset by a normal exit; two runs in a row that do not end normally switch to software graphics

# Graphics: hardware compositing by default; software rendering on remote desktops, machines without a real display
# adapter, after repeated abnormal exits, or when asked (--safe-graphics / TERRATOOLS_GPU=off).
GRAPHICS_MODE, GRAPHICS_WHY = diag.graphics_mode(sys.argv)
if GRAPHICS_MODE == "software":
    _flags = os.environ.get("QTWEBENGINE_CHROMIUM_FLAGS", "")
    if "--disable-gpu" not in _flags:
        os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = (_flags + " --disable-gpu").strip()
diag.event("graphics", mode=GRAPHICS_MODE, why=GRAPHICS_WHY, flags=os.environ.get("QTWEBENGINE_CHROMIUM_FLAGS", ""))

# The Chromium remote-debugging port is a development aid: an open local port in an installed app is both a
# security exposure and a source of conflicts (another Chromium app, a second TerraTools). Packaged builds leave
# it off unless TERRATOOLS_DEBUG is set.
if not getattr(sys, "frozen", False) or os.environ.get("TERRATOOLS_DEBUG"):
    os.environ["QTWEBENGINE_REMOTE_DEBUGGING"] = "9222"
if sys.platform.startswith("win"):
    # ANGLE (Direct3D) backend selection — Windows-only; macOS/Linux Qt
    # builds use native OpenGL and ignore this variable.
    os.environ["QT_ANGLE_PLATFORM"] = "d3d11"
    # Every process TerraTools starts (Rscript, workers) dies with it, however it exits.
    from tools.common.win_job import kill_children_on_exit
    kill_children_on_exit()
from PySide6.QtWidgets import QApplication

from ui.main_window import MainWindow


def main():
    if sys.platform.startswith("win"):
        # Own taskbar identity, so Windows groups/icons TerraTools as its
        # own app rather than as "python.exe".
        import ctypes
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("IRRI.TerraTools")

    app = QApplication(sys.argv)
    app.setApplicationName("TerraTools")
    app.setApplicationDisplayName("TerraTools")
    app.setOrganizationName("IRRI")
    diag.start_gui_watchdog()

    window = MainWindow()
    window.showMaximized()

    if sys.platform == "darwin" and getattr(sys, "frozen", False):
        # Packaged Mac app: unpack the bundled R runtime in the background
        # (first launch only) so it's ready before an R tool is used.
        import threading
        from tools.common.r_locator import prepare_macos_r
        threading.Thread(target=prepare_macos_r, daemon=True).start()

    code = app.exec()
    diag.event("exit", code=code)
    diag.mark_clean_exit()
    sys.exit(code)


if __name__ == "__main__":
    main()
