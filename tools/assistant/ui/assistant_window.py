import json
import threading
from pathlib import Path

from PySide6.QtCore import QObject, Signal, Slot
from PySide6.QtWidgets import QFileDialog

from tools.assistant import config
from tools.assistant.agent import Agent
from tools.assistant.proxy_client import AssistantError
from tools.assistant.run_files import RunStore

# QWebChannel object name of each tool bridge -> the tool id used in the assistant's tools
OBJECT_TO_TOOL = {
    "mendelianBridge": "mendelian", "linkageBridge": "linkage", "qtlBridge": "qtl",
    "geneExtractorBridge": "gene_extractor", "genePrioritizationBridge": "gene_prioritization", "pybsaseqBridge": "pybsaseq", "qtlSeqBridge": "qtl_seq",
    "diversityBridge": "diversity", "gwasBridge": "gwas", "etGwasBridge": "et_gwas", "clusteringBridge": "clustering",
    "additiveDominanceBridge": "additive_dominance", "genomicPredictionBridge": "genomic_prediction",
    "crossBridge": "cross_selection", "heterosisBridge": "heterosis", "phenotypicBridge": "phenotypic",
}


class AssistantBridge(QObject):
    """Chat backend for the Assistant page.

    JS -> Python: ask(text), cancel(), new_chat(), attach_file(), status(), files_read()
    Python -> JS: chunk(text), toolUsed(json), notice(text), answerDone(json), answerFailed(text)
    """

    chunk = Signal(str)
    toolUsed = Signal(str)
    notice = Signal(str)
    answerDone = Signal(str)
    answerFailed = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.store = RunStore()
        self.agent = Agent(self.store)
        self._busy = False
        self._stop = threading.Event()
        self._attached = []            # names attached since the last question

    # ---- fed by the tool bridges (see tool_registry.register_tool_bridges)
    def note_result(self, tool, payload):
        try:
            self.store.note_result(tool, payload)
        except Exception:               # never let bookkeeping disturb a tool run
            pass

    # ---- chat
    @Slot(str)
    def ask(self, text):
        text = (text or "").strip()
        if not text or self._busy:
            return
        if self._attached:
            names = ", ".join(self._attached)
            text += f"\n\n(The user attached: {names}. Use list_runs / read_run_file to look at them.)"
            self._attached = []
        self._busy = True
        self._stop.clear()
        threading.Thread(target=self._work, args=(text,), daemon=True).start()

    def _work(self, text):
        def emit(kind, payload):
            if kind == "delta":
                self.chunk.emit(payload)
            elif kind == "tool":
                self.toolUsed.emit(json.dumps(payload))
            else:
                self.notice.emit(str(payload))
        try:
            answer = self.agent.ask(text, emit, self._stop.is_set)
            self.answerDone.emit(json.dumps({"text": answer}))
        except AssistantError as e:
            self.answerFailed.emit(str(e))
        except Exception as e:          # keep the page usable whatever happens
            print("[ASSISTANT] unexpected error:", repr(e), flush=True)
            self.answerFailed.emit("Something went wrong in the assistant. Please try again.")
        finally:
            self._busy = False

    @Slot()
    def cancel(self):
        self._stop.set()

    @Slot()
    def new_chat(self):
        if self._busy:
            self._stop.set()
        self.agent.reset()
        self._attached = []

    @Slot(result=str)
    def attach_file(self):
        path, _ = QFileDialog.getOpenFileName(
            None, "Attach a file for the assistant to read", "",
            "Data and plots (*.csv *.tsv *.txt *.xlsx *.json *.png *.jpg *.jpeg *.log);;All files (*)")
        if not path:
            return ""
        name = self.store.attach(path)
        if name:
            self._attached.append(name)
        return name or ""

    @Slot(result=str)
    def status(self):
        cfg = config.load()
        runs = self.store.list_runs()
        return json.dumps({
            "configured": bool(cfg["proxy_url"]),
            "runs": [r["tool"] for r in runs.get("runs", [])],
        })

    @Slot(result=str)
    def files_read(self):
        """Names of the files/results the assistant has opened in this chat (for transparency)."""
        return json.dumps([name for _t, name in self.store.log[-200:]])
