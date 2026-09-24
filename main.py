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
 
    else:
        print(f"Unknown worker: {worker_name}", file=sys.stderr)
        sys.exit(1)
 
    sys.exit(0)

# ── Normal GUI startup continues below, unchanged ──────────────────────
# (your existing os.environ["QTWEBENGINE_REMOTE_DEBUGGING"] = "9222" line,
#  the PySide6 imports, QApplication creation, install_polish(), etc.)
os.environ["QTWEBENGINE_REMOTE_DEBUGGING"] = "9222"
if sys.platform.startswith("win"):
    # ANGLE (Direct3D) backend selection — Windows-only; macOS/Linux Qt
    # builds use native OpenGL and ignore this variable.
    os.environ["QT_ANGLE_PLATFORM"] = "d3d11"
from PySide6.QtWidgets import QApplication

from ui.main_window import MainWindow


def main():
    app = QApplication(sys.argv)

    window = MainWindow()
    window.show()

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
