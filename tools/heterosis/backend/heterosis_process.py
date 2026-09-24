# -*- coding: utf-8 -*-
"""
heterosis_process.py

Standalone subprocess entry point for the heavy heterosis
computation (TensorFlow/Keras training).

WHY THIS FILE EXISTS
---------------------
TensorFlow bundles its own copy of LLVM/XLA/MLIR. PySide6's
QtWebEngine embeds Chromium, which also bundles LLVM-based
components. Each registers its own global command-line flags at
import time. If TensorFlow and QtWebEngine are ever imported into
the *same* process, LLVM's global flag registry sees a flag such as
`disable-auto-upgrade-debug-info` registered twice and the whole
process aborts:

    CommandLine Error: Option 'disable-auto-upgrade-debug-info'
    registered more than once!
    LLVM ERROR: inconsistency in registered CommandLine options
    Aborted (core dumped)

The fix is process isolation, not an import-order trick: this
script runs as a *separate* `python` process (launched via
subprocess by HeterosisWorker in heterosis_window.py), so
TensorFlow only ever loads in a process that never touches
QtWebEngine.

USAGE
-----
    python heterosis_process.py <config_json_path> <results_dir>

Writes newline-delimited status lines to stdout so the parent
process can stream progress without blocking on the full run:

    PROGRESS <fraction_or_-1> <message>
    RESULT <json>

Exactly one RESULT line is written, always last.
"""

import sys
import json
import traceback


def main():

    if len(sys.argv) < 3:

        print("RESULT " + json.dumps({
            "status": "error",
            "message": "Usage: heterosis_process.py <config.json> <results_dir>"
        }), flush=True)

        return

    config_path = sys.argv[1]
    results_dir = sys.argv[2]

    try:

        with open(config_path, "r") as f:
            config = json.load(f)

    except Exception as e:

        print("RESULT " + json.dumps({
            "status": "error",
            "message": f"Could not read config file: {e}"
        }), flush=True)

        return

    # Imported here, inside the subprocess, so TensorFlow never
    # loads in the parent (GUI) process.
    try:

        from tools.heterosis.backend import heterosis_runner

    except Exception as e:

        traceback.print_exc()

        print("RESULT " + json.dumps({
            "status": "error",
            "message": f"Could not import heterosis_runner: {e}"
        }), flush=True)

        return

    def progress_callback(fraction, message):

        frac = fraction if fraction is not None else -1

        print(f"PROGRESS {frac} {message or ''}", flush=True)

    try:

        result = heterosis_runner.run_heterosis(
            config,
            results_dir=results_dir,
            progress_callback=progress_callback
        )

    except Exception as e:

        traceback.print_exc()

        result = {
            "status": "error",
            "message": str(e)
        }

    print("RESULT " + json.dumps(result), flush=True)


if __name__ == "__main__":
    main()