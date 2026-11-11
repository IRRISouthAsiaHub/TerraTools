"""Where the assistant finds its service.

The desktop app never holds an Anthropic key: it talks to the TerraTools
Assistant proxy (server/assistant_proxy), which does. The proxy address ships
with the build in `assistant_config.json`; the environment or a per-user file
can override it (useful for testing against a local proxy).

    {"proxy_url": "https://assistant.example.org", "app_token": "..."}

`app_token` is not a secret in the cryptographic sense (it is inside every
installed copy); it only keeps casual traffic away. The proxy's rate limits
and daily budget are what protect the quota.
"""

import json
import os
import sys
import uuid
from pathlib import Path

HERE = Path(__file__).resolve().parent
USER_DIR = Path.home() / ".terratools"


def _bundled_dir():
    if getattr(sys, "frozen", False):
        return Path(getattr(sys, "_MEIPASS", HERE)) / "tools" / "assistant"
    return HERE


def load():
    """-> {"proxy_url": str, "app_token": str}; proxy_url is "" when not configured."""
    cfg = {"proxy_url": "", "app_token": ""}
    for path in (_bundled_dir() / "assistant_config.json", USER_DIR / "assistant.json"):
        try:
            cfg.update({k: v for k, v in json.loads(path.read_text(encoding="utf-8")).items() if k in cfg and v})
        except (OSError, ValueError):
            pass
    cfg["proxy_url"] = (os.environ.get("TERRATOOLS_ASSISTANT_URL") or cfg["proxy_url"]).rstrip("/")
    cfg["app_token"] = os.environ.get("TERRATOOLS_ASSISTANT_TOKEN") or cfg["app_token"]
    return cfg


def install_id():
    """Random per-installation id (no personal data); lets the proxy rate-limit fairly."""
    path = USER_DIR / "install_id"
    try:
        return path.read_text().strip() or _new_id(path)
    except OSError:
        return _new_id(path)


def _new_id(path):
    value = uuid.uuid4().hex
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(value)
    except OSError:
        pass
    return value
