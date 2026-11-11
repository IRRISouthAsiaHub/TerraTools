"""
One HTTP client per database: User-Agent, timeout, retries with exponential backoff, a minimum interval between
requests (rate limit), the SQLite cache and a cancel flag. A failed service never produces data: the caller gets an
exception (NetworkError / ApiError) and records the evidence as unavailable.
"""

from __future__ import annotations

import hashlib
import json
import random
import threading
import time
from datetime import datetime, timezone

import requests

from tools.gene_prioritization import config


class Cancelled(Exception):
    pass


class NetworkError(Exception):
    """The service could not be reached (no internet, DNS, timeout after all retries)."""


class ApiError(NetworkError):
    """The service answered with an error that retrying does not fix. It is handled like an unreachable service:
    the evidence that needed the answer is reported as unavailable and the run goes on."""


def today():
    return datetime.now(timezone.utc).strftime("%Y-%m-%d")


class RunContext:
    """Shared by every client of one run: the cancel flag, the cache, and whether cached entries are ignored."""

    def __init__(self, cache, refresh=False):
        self.cache = cache
        self.refresh = refresh
        self.cancel = threading.Event()

    def check(self):
        if self.cancel.is_set():
            raise Cancelled()


class ApiClient:
    def __init__(self, source, base_url, ctx, min_interval=0.1, headers=None, extra_params=None):
        self.source = source
        self.base = base_url.rstrip("/")
        self.ctx = ctx
        self.min_interval = min_interval
        self.extra = dict(extra_params or {})
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": config.USER_AGENT, "Accept": "application/json"})
        if headers:
            self.session.headers.update(headers)
        self._gate = threading.Lock()
        self._last = 0.0
        self.calls = 0
        self.cached = 0
        self.failures = 0
        self.last_error = ""
        self.ok_once = False
        self.consecutive = 0          # failed attempts in a row; a service that keeps failing is skipped instead of waited for
        self.tripped = False

    # ------------------------------------------------------------ rate limit
    def _wait(self):
        with self._gate:
            delay = self._last + self.min_interval - time.monotonic()
            if delay > 0:
                time.sleep(delay)
            self._last = time.monotonic()

    # ------------------------------------------------------------ request
    def request(self, path, params=None, body=None, ttl=None, text=False, retries=None, timeout=None):
        """JSON (or text) of a GET (POST when `body` is given). Returns None when the service says 'not found'.
        Responses are cached with their retrieval time; `ttl` is in days."""
        self.ctx.check()
        url = path if path.startswith("http") else f"{self.base}/{path.lstrip('/')}"
        params = {**self.extra, **(params or {})}
        ident = json.dumps([url, sorted((k, str(v)) for k, v in params.items()), body], sort_keys=True)
        key = hashlib.sha1(ident.encode("utf-8")).hexdigest()
        max_age = config.CACHE_DAYS if ttl is None else ttl
        if not self.ctx.refresh and max_age > 0:
            hit, payload, _ = self.ctx.cache.get(key, max_age, config.NEGATIVE_CACHE_DAYS)
            if hit:
                self.cached += 1
                return payload

        if self.tripped:
            raise NetworkError(f"{self.source} stopped answering during this run ({self.last_error or 'unreachable'}); skipped.")
        tries = retries or config.RETRIES
        wait_for = timeout or config.TIMEOUT
        delay, last = config.BACKOFF, None
        for attempt in range(tries):
            self.ctx.check()
            self._wait()
            try:
                self.calls += 1
                if body is None:
                    r = self.session.get(url, params=params, timeout=wait_for)
                else:
                    r = self.session.post(url, params=params, json=body, timeout=wait_for)
            except (requests.ConnectionError, requests.Timeout) as e:
                last = f"{type(e).__name__}"
                self.failures += 1
                self.consecutive += 1
            else:
                if r.status_code == 200:
                    try:
                        # strict=False: NCBI puts raw line breaks inside JSON strings (the translated query)
                        payload = r.text if text else json.loads(r.text, strict=False)
                    except ValueError:
                        raise ApiError(f"{self.source}: the answer was not {'text' if text else 'JSON'}")
                    self.ok_once = True
                    self.consecutive = 0
                    self.ctx.cache.put(key, self.source, r.url, 200, payload)
                    return payload
                if r.status_code in (400, 404):                 # nothing for this identifier
                    self.ok_once = True
                    self.consecutive = 0
                    self.ctx.cache.put(key, self.source, r.url, r.status_code, None)
                    return None
                if r.status_code in (429, 500, 502, 503, 504):
                    last = f"HTTP {r.status_code}"
                    self.failures += 1
                    self.consecutive += 1
                    retry_after = r.headers.get("Retry-After")
                    if retry_after and retry_after.isdigit():
                        delay = max(delay, min(float(retry_after), 30))
                else:
                    self.last_error = f"HTTP {r.status_code}"
                    raise ApiError(f"{self.source}: HTTP {r.status_code} for {r.url}")
            if self.consecutive >= config.BREAKER_FAILURES:
                self.tripped = True
                self.last_error = last or "unreachable"
                raise NetworkError(f"{self.source} stopped answering ({self.last_error}); it is skipped for the rest of this run.")
            if attempt < tries - 1:
                time.sleep(delay + random.uniform(0, 0.4))
                delay *= 2
        self.last_error = last or "unreachable"
        raise NetworkError(f"{self.source} could not be reached ({self.last_error}).")

    def status(self):
        if self.failures and not self.ok_once:
            state = "unavailable"
        elif self.failures:
            state = "partial"
        elif self.calls or self.cached:
            state = "ok"
        else:
            state = "not used"
        return {"source": self.source, "state": state, "calls": self.calls, "cached": self.cached,
                "failures": self.failures, "error": self.last_error}
