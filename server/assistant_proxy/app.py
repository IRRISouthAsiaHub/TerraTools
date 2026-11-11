"""TerraTools Assistant proxy.

Holds the Anthropic API key so the desktop app never has to. The app sends the
conversation (messages only); this service adds the system prompt and the tool
schemas from tools/assistant/{prompts,toolspec}.py, pins the model, limits the
size and rate of requests, and streams the reply back as server-sent events:

    event: delta   data: {"text": "..."}          streamed answer text
    event: final   data: {"message": {...}, "usage": {...}}   the full assistant message
    event: error   data: {"message": "..."}

Run:   uvicorn server.assistant_proxy.app:app --host 0.0.0.0 --port 8787
Needs: ANTHROPIC_API_KEY (or another credential the SDK resolves) in the environment.
See README.md for deployment and the other settings (TT_*).
"""

import asyncio
import collections
import json
import logging
import os
import sys
import time
from pathlib import Path

from fastapi import FastAPI, Header, HTTPException, Request
from fastapi.responses import JSONResponse, StreamingResponse

# The prompt and tool schemas live with the client code so both stay in step.
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from tools.assistant.prompts import SYSTEM_PROMPT          # noqa: E402
from tools.assistant.toolspec import TOOLS, TOOL_NAMES      # noqa: E402

log = logging.getLogger("tt-assistant")
logging.basicConfig(level=os.environ.get("TT_LOG", "INFO"))

MODEL = os.environ.get("TT_MODEL", "claude-opus-5-5")
EFFORT = os.environ.get("TT_EFFORT", "medium")
MAX_TOKENS = int(os.environ.get("TT_MAX_TOKENS", "8000"))
MAX_BODY = int(os.environ.get("TT_MAX_BODY_BYTES", str(8 * 1024 * 1024)))   # images travel in tool results
MAX_MESSAGES = int(os.environ.get("TT_MAX_MESSAGES", "160"))
PER_INSTALL_HOUR = int(os.environ.get("TT_PER_INSTALL_PER_HOUR", "150"))    # one question = several calls
PER_INSTALL_DAY = int(os.environ.get("TT_PER_INSTALL_PER_DAY", "600"))
PER_IP_HOUR = int(os.environ.get("TT_PER_IP_PER_HOUR", "400"))
DAILY_TOKENS = int(os.environ.get("TT_DAILY_TOKEN_BUDGET", "20000000"))     # all users together
TOKENS = {t.strip() for t in os.environ.get("TT_APP_TOKENS", "").split(",") if t.strip()}
MOCK = os.environ.get("TT_MOCK") == "1"                                       # no upstream call (plumbing tests)

USER_BLOCKS = {"text", "image", "tool_result"}
ASSISTANT_BLOCKS = {"text", "thinking", "redacted_thinking", "tool_use", "fallback"}

app = FastAPI(title="TerraTools Assistant proxy", docs_url=None, redoc_url=None)


# ------------------------------------------------------------------ limits
class Window:
    """Sliding-window counter (in memory: resets when the service restarts)."""

    def __init__(self, span):
        self.span = span
        self.hits = collections.defaultdict(collections.deque)

    def take(self, key, limit):
        now = time.time()
        q = self.hits[key]
        while q and q[0] < now - self.span:
            q.popleft()
        if len(q) >= limit:
            return False
        q.append(now)
        return True


_install_h, _install_d, _ip_h = Window(3600), Window(86400), Window(3600)
_spent = {"day": time.strftime("%Y-%m-%d"), "tokens": 0}


def _budget_ok():
    today = time.strftime("%Y-%m-%d")
    if _spent["day"] != today:
        _spent.update(day=today, tokens=0)
    return _spent["tokens"] < DAILY_TOKENS


def _client_ip(request: Request):
    fwd = request.headers.get("x-forwarded-for")       # behind a reverse proxy / platform router
    return (fwd.split(",")[0].strip() if fwd else (request.client.host if request.client else "?"))


# ------------------------------------------------------------------ validation
def _validate(body):
    msgs = body.get("messages")
    if not isinstance(msgs, list) or not msgs:
        raise HTTPException(400, "messages must be a non-empty list")
    if len(msgs) > MAX_MESSAGES:
        raise HTTPException(413, "This conversation is too long - please start a new chat.")
    for i, m in enumerate(msgs):
        role, content = m.get("role"), m.get("content")
        if role not in ("user", "assistant"):
            raise HTTPException(400, f"message {i}: role must be user or assistant")
        if isinstance(content, str):
            continue
        if not isinstance(content, list):
            raise HTTPException(400, f"message {i}: bad content")
        allowed = USER_BLOCKS if role == "user" else ASSISTANT_BLOCKS
        for b in content:
            t = b.get("type") if isinstance(b, dict) else None
            if t not in allowed:
                raise HTTPException(400, f"message {i}: block type {t!r} not allowed for {role}")
            if t == "tool_use" and b.get("name") not in TOOL_NAMES:
                raise HTTPException(400, f"message {i}: unknown tool {b.get('name')!r}")
    if msgs[-1]["role"] != "user":
        raise HTTPException(400, "the last message must be from the user")
    return msgs


def _sse(event, data):
    return f"event: {event}\ndata: {json.dumps(data, ensure_ascii=False)}\n\n"


# ------------------------------------------------------------------ routes
@app.get("/health")
async def health():
    return {"ok": True, "model": MODEL, "mock": MOCK, "token_budget_left": max(0, DAILY_TOKENS - _spent["tokens"])}


@app.post("/v1/chat")
async def chat(request: Request, x_tt_token: str = Header(default=""), x_tt_install: str = Header(default="")):
    if TOKENS and x_tt_token not in TOKENS:
        raise HTTPException(401, "Unauthorized")
    install = x_tt_install.strip()[:64] or "anonymous"
    ip = _client_ip(request)

    raw = await request.body()
    if len(raw) > MAX_BODY:
        raise HTTPException(413, "Request too large")
    try:
        body = json.loads(raw)
    except ValueError:
        raise HTTPException(400, "Invalid JSON")
    messages = _validate(body)

    if not _budget_ok():
        raise HTTPException(429, "The assistant has reached its daily usage limit. Please try again tomorrow.")
    if not (_ip_h.take(ip, PER_IP_HOUR) and _install_h.take(install, PER_INSTALL_HOUR) and _install_d.take(install, PER_INSTALL_DAY)):
        raise HTTPException(429, "You are sending requests too quickly. Please wait a little and try again.")

    stream = _mock_stream(messages) if MOCK else _claude_stream(messages)
    return StreamingResponse(stream, media_type="text/event-stream", headers={"Cache-Control": "no-store", "X-Accel-Buffering": "no"})


# ------------------------------------------------------------------ upstream
_client = None


def _anthropic():
    global _client
    if _client is None:
        import anthropic                       # imported lazily so mock mode needs no SDK
        _client = anthropic.AsyncAnthropic()   # reads ANTHROPIC_API_KEY / profile from the environment
    return _client


async def _claude_stream(messages):
    import anthropic
    try:
        async with _anthropic().beta.messages.stream(
            model=MODEL,
            max_tokens=MAX_TOKENS,
            system=[{"type": "text", "text": SYSTEM_PROMPT, "cache_control": {"type": "ephemeral"}}],
            tools=TOOLS,
            tool_choice={"type": "auto"},          # forced tool use is rejected by this model family
            thinking={"type": "adaptive"},
            output_config={"effort": EFFORT},
            messages=messages,
            betas=["server-side-fallback-2026-07-01"],
            fallbacks="default",                  # re-run on another model if a safety classifier declines
        ) as stream:
            async for text in stream.text_stream:
                yield _sse("delta", {"text": text})
            final = await stream.get_final_message()
        usage = final.usage
        used = (usage.input_tokens or 0) + (usage.output_tokens or 0) + (getattr(usage, "cache_creation_input_tokens", 0) or 0)
        _spent["tokens"] += used
        yield _sse("final", {
            "message": {"role": "assistant", "content": [b.model_dump(mode="json", exclude_none=True) for b in final.content]},
            "stop_reason": final.stop_reason,
            "usage": {"input": usage.input_tokens, "output": usage.output_tokens,
                      "cache_read": getattr(usage, "cache_read_input_tokens", 0)},
        })
    except anthropic.RateLimitError:
        yield _sse("error", {"message": "The assistant is busy right now. Please try again in a minute."})
    except anthropic.AuthenticationError:
        log.error("Anthropic credentials rejected - check ANTHROPIC_API_KEY on the proxy")
        yield _sse("error", {"message": "The assistant is temporarily unavailable (service configuration)."})
    except anthropic.BadRequestError as e:
        log.warning("bad request upstream: %s", e)
        yield _sse("error", {"message": "The assistant could not process this conversation. Try starting a new chat."})
    except anthropic.APIConnectionError:
        yield _sse("error", {"message": "The assistant could not reach its model service. Please try again."})
    except anthropic.APIStatusError as e:
        log.error("upstream %s: %s", e.status_code, e)
        yield _sse("error", {"message": "The assistant is temporarily unavailable. Please try again later."})
    except Exception:
        log.exception("unexpected error")
        yield _sse("error", {"message": "Something went wrong in the assistant service."})


async def _mock_stream(messages):
    """Canned replies: exercises streaming, the tool loop and the error path without an API key."""
    last = messages[-1]["content"]
    text = last if isinstance(last, str) else " ".join(b.get("text", "") for b in last if b.get("type") == "text")
    has_result = isinstance(last, list) and any(b.get("type") == "tool_result" for b in last)
    if "BOOM" in text:
        yield _sse("error", {"message": "mock error"})
        return
    if "TOOLTEST" in text and not has_result:
        content = [{"type": "text", "text": "Let me check the manual. "},
                   {"type": "tool_use", "id": "toolu_mock1", "name": "search_manual", "input": {"query": "QTL genotype file"}}]
        for c in content[:1]:
            yield _sse("delta", {"text": c["text"]})
        yield _sse("final", {"message": {"role": "assistant", "content": content}, "stop_reason": "tool_use", "usage": {"input": 10, "output": 5}})
        return
    reply = "Mock answer after tool result." if has_result else f"Mock answer to: {text[:80]}"
    for w in reply.split(" "):
        yield _sse("delta", {"text": w + " "})
        await asyncio.sleep(0.01)
    yield _sse("final", {"message": {"role": "assistant", "content": [{"type": "text", "text": reply}]}, "stop_reason": "end_turn", "usage": {"input": 10, "output": 5}})


@app.exception_handler(HTTPException)
async def _http_error(_request, exc: HTTPException):
    return JSONResponse({"error": exc.detail}, status_code=exc.status_code)
