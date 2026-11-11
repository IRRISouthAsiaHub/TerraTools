"""
SQLite cache of API responses and of assembled gene annotations.

Both tables keep the retrieval time, so every piece of evidence can be traced and old entries can be refreshed.
A single connection is shared by the worker threads and guarded by a lock.
"""

from __future__ import annotations

import json
import sqlite3
import threading
import time


class Cache:
    def __init__(self, path):
        self.path = str(path)
        self._lock = threading.Lock()
        self._db = sqlite3.connect(self.path, check_same_thread=False)
        self._db.execute("PRAGMA journal_mode=WAL")
        self._db.execute(
            "CREATE TABLE IF NOT EXISTS responses (key TEXT PRIMARY KEY, source TEXT, url TEXT, "
            "fetched REAL, status INTEGER, body TEXT)")
        self._db.execute(
            "CREATE TABLE IF NOT EXISTS annotations (species TEXT, gene TEXT, fetched REAL, body TEXT, "
            "PRIMARY KEY (species, gene))")
        self._db.commit()

    # ---------------------------------------------------------------- responses
    def get(self, key, max_age_days, miss_age_days=7):
        """(found, payload, fetched_epoch). A stored miss (404) comes back as (True, None, t) and expires sooner."""
        with self._lock:
            row = self._db.execute("SELECT status, body, fetched FROM responses WHERE key=?", (key,)).fetchone()
        if not row:
            return False, None, 0.0
        status, body, fetched = row
        limit = max_age_days if status == 200 else miss_age_days
        if time.time() - fetched > limit * 86400:
            return False, None, fetched
        return True, (None if status != 200 or body is None else json.loads(body)), fetched

    def put(self, key, source, url, status, payload):
        body = None if payload is None else json.dumps(payload)
        with self._lock:
            self._db.execute("REPLACE INTO responses VALUES (?,?,?,?,?,?)", (key, source, url, time.time(), status, body))
            self._db.commit()

    # ---------------------------------------------------------------- annotations
    def get_annotation(self, species, gene, max_age_days):
        with self._lock:
            row = self._db.execute("SELECT body, fetched FROM annotations WHERE species=? AND gene=?",
                                   (species, gene)).fetchone()
        if not row or time.time() - row[1] > max_age_days * 86400:
            return None
        return json.loads(row[0])

    def put_annotation(self, species, gene, record):
        with self._lock:
            self._db.execute("REPLACE INTO annotations VALUES (?,?,?,?)", (species, gene, time.time(), json.dumps(record)))
            self._db.commit()

    # ---------------------------------------------------------------- housekeeping
    def clear(self):
        with self._lock:
            self._db.execute("DELETE FROM responses")
            self._db.execute("DELETE FROM annotations")
            self._db.commit()

    def stats(self):
        with self._lock:
            n_r = self._db.execute("SELECT COUNT(*) FROM responses").fetchone()[0]
            n_a = self._db.execute("SELECT COUNT(*) FROM annotations").fetchone()[0]
        return {"responses": n_r, "annotations": n_a}

    def close(self):
        with self._lock:
            try:
                self._db.close()
            except sqlite3.Error:
                pass
