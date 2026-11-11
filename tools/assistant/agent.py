"""The assistant's conversation loop, running inside the desktop app.

The model lives behind the proxy; the tools it calls (manual search, run
files, literature search) run here, on the user's machine. The conversation
history is append-only: assistant turns are stored exactly as the service sent
them (including thinking blocks), and a failed turn is removed from the end
as a whole.
"""

import json

from . import literature, proxy_client
from .knowledge import Manual
from .proxy_client import AssistantError
from .run_files import RunStore
from .toolspec import TOOL_NAMES

MAX_ROUNDS = 12                 # tool calls chained for one question
MAX_RESULT_CHARS = 14000
CONTEXT_WARN_TOKENS = 600_000   # the model's window is far larger; this just keeps chats snappy


def _text_block(s):
    return [{"type": "text", "text": s[:MAX_RESULT_CHARS]}]


class Agent:
    def __init__(self, store: RunStore):
        self.store = store
        self.history = []
        self.manual = Manual()
        self.last_input_tokens = 0

    def reset(self):
        self.history = []
        self.last_input_tokens = 0

    # ------------------------------------------------------------ tools
    def _run_tool(self, name, args):
        """-> (content blocks, is_error, short label for the UI)"""
        if name not in TOOL_NAMES:
            return _text_block(f"Unknown tool {name}"), True, name
        try:
            if name == "search_manual":
                hits = self.manual.search(args["query"], k=min(int(args.get("max_results") or 5), 8), chapter=args.get("chapter"))
                if not hits:
                    return _text_block("No passage in the manual matches. Try other keywords."), False, f"manual: {args['query']}"
                out = "\n\n".join(f"[{h['chapter']}{' > ' + h['section'] if h['section'] else ''}]\n{h['text']}" for h in hits)
                return _text_block(out), False, f"manual: {args['query']}"
            if name == "search_literature":
                res = literature.search(
                    args["query"], int(args.get("max_results") or 6), args.get("year_from"),
                    bool(args.get("open_access_only")), args.get("sort") or "relevance")
                return _text_block(json.dumps(res, ensure_ascii=False)), False, f"papers: {args['query']}"
            if name == "verify_doi":
                return _text_block(json.dumps(literature.lookup_doi(args["doi"]), ensure_ascii=False)), False, f"doi: {args['doi']}"
            if name == "list_runs":
                return _text_block(json.dumps(self.store.list_runs(), ensure_ascii=False)), False, "your runs"
            if name == "get_run_result":
                return _text_block(json.dumps(self.store.get_run_result(args["tool"]), ensure_ascii=False, default=str)), False, f"results: {args['tool']}"
            if name == "read_run_file":
                blocks, err = self.store.read_file(args["name"], args.get("tool"), args.get("offset", 0), args.get("limit", 20), args.get("sheet"))
                return blocks, err, f"file: {args['name']}"
        except KeyError as e:
            return _text_block(f"Missing argument {e}"), True, name
        except Exception as e:                      # a broken tool must not break the chat
            return _text_block(f"{name} failed: {type(e).__name__}: {e}"), True, name
        return _text_block("Unhandled tool"), True, name

    # ------------------------------------------------------------ one question
    def ask(self, text, emit, cancelled=lambda: False):
        """emit(kind, payload): kinds 'delta', 'tool', 'info'. Returns the final answer text.
        Raises AssistantError (message is user-safe)."""
        if self.last_input_tokens > CONTEXT_WARN_TOKENS:
            raise AssistantError("This chat has become very long. Please start a new chat (the New chat button) to keep answers fast and accurate.")
        mark = len(self.history)
        self.history.append({"role": "user", "content": text})
        answer = []
        try:
            for _ in range(MAX_ROUNDS):
                if cancelled():
                    raise AssistantError("Stopped.")
                final = proxy_client.stream_chat(self.history, lambda t: (answer.append(t), emit("delta", t)), cancelled)
                msg = final["message"]
                self.history.append(msg)
                self.last_input_tokens = (final.get("usage") or {}).get("input") or self.last_input_tokens
                calls = [b for b in msg["content"] if b.get("type") == "tool_use"]
                stop = final.get("stop_reason")
                if stop == "refusal":
                    emit("info", "The assistant declined to answer this one.")
                    break
                if not calls or stop not in ("tool_use", None):
                    if stop == "max_tokens":
                        emit("info", "The answer was cut short; ask me to continue.")
                    break
                results = []
                for call in calls:
                    blocks, is_error, label = self._run_tool(call["name"], call.get("input") or {})
                    emit("tool", {"name": call["name"], "label": label, "error": is_error})
                    results.append({"type": "tool_result", "tool_use_id": call["id"], "content": blocks, "is_error": is_error})
                self.history.append({"role": "user", "content": results})
                if answer:                  # visual gap between the text before and after tool calls
                    answer.append("\n\n")
                    emit("delta", "\n\n")
            else:
                emit("info", "I stopped after several lookups; ask again to continue.")
        except Exception:
            del self.history[mark:]       # drop the whole failed turn; earlier history is untouched
            raise
        return "".join(answer).strip()
