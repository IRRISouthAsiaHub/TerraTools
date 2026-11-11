"""Talks to the assistant proxy (HTTP + server-sent events), stdlib only."""

import json
import urllib.error
import urllib.request

from . import config

TIMEOUT = 300           # a long answer with thinking can take minutes; this is per socket read


class AssistantError(Exception):
    """Message is safe to show to the user."""


def _parse_sse(lines):
    event, data = None, []
    for raw in lines:
        line = raw.decode("utf-8", errors="replace").rstrip("\r\n")
        if not line:
            if event and data:
                yield event, json.loads("\n".join(data))
            event, data = None, []
        elif line.startswith("event:"):
            event = line[6:].strip()
        elif line.startswith("data:"):
            data.append(line[5:].lstrip())


def stream_chat(messages, on_delta, cancelled=lambda: False):
    """Send the conversation, call on_delta(text) for streamed text, and return
    {"message": assistant message, "stop_reason", "usage"}."""
    cfg = config.load()
    if not cfg["proxy_url"]:
        raise AssistantError(
            "The assistant is not set up in this copy of TerraTools (no service address). "
            "Ask whoever installed TerraTools to configure it - see server/assistant_proxy/README.md.")
    req = urllib.request.Request(
        cfg["proxy_url"] + "/v1/chat",
        data=json.dumps({"messages": messages}).encode("utf-8"),
        headers={
            "Content-Type": "application/json", "Accept": "text/event-stream",
            "X-TT-Token": cfg["app_token"], "X-TT-Install": config.install_id(),
        },
        method="POST",
    )
    try:
        resp = urllib.request.urlopen(req, timeout=TIMEOUT)
    except urllib.error.HTTPError as e:
        try:
            detail = json.loads(e.read().decode("utf-8", "replace")).get("error")
        except Exception:
            detail = None
        raise AssistantError(detail if isinstance(detail, str) else f"The assistant service returned an error ({e.code}).")
    except (urllib.error.URLError, TimeoutError, OSError):
        raise AssistantError("Could not reach the assistant service. Check your internet connection and try again.")

    final = None
    try:
        with resp:
            for event, data in _parse_sse(resp):
                if cancelled():
                    raise AssistantError("Stopped.")
                if event == "delta":
                    on_delta(data.get("text", ""))
                elif event == "final":
                    final = data
                elif event == "error":
                    raise AssistantError(data.get("message", "The assistant failed."))
    except (TimeoutError, OSError):
        raise AssistantError("The connection to the assistant was interrupted. Please try again.")
    if final is None:
        raise AssistantError("The assistant did not finish its answer. Please try again.")
    return final
