"""Search over the TerraTools user manual (the assistant's grounding).

The manual HTML is turned into small heading-delimited chunks once
(`python -m tools.assistant.knowledge`, which writes data/manual_chunks.json);
at run time the chunks are ranked with a small BM25 so no extra library or
network access is needed.
"""

import json
import math
import re
import sys
from html.parser import HTMLParser
from pathlib import Path

HERE = Path(__file__).resolve().parent
CHUNKS_FILE = HERE / "data" / "manual_chunks.json"
MANUAL_HTML = HERE.parents[1] / "docs" / "manual" / "TerraTools_User_Manual.html"

_SKIP = {"script", "style", "figure", "svg"}   # img is a void tag: it never closes
_BLOCK = {"p", "li", "tr", "h1", "h2", "h3", "h4", "div", "pre", "table", "ul", "ol", "br"}
MAX_CHARS = 1800


class _Flatten(HTMLParser):
    """Walk the manual and cut it at every h2/h3 heading."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.chunks = []
        self.h2 = ""
        self.h3 = ""
        self._buf = []
        self._skip = 0
        self._in_head = None
        self._head = []
        self._cell = False

    def _flush(self):
        text = re.sub(r"[ \t]+", " ", "".join(self._buf))
        text = re.sub(r"\n\s*\n+", "\n", text).strip()
        self._buf = []
        if len(text) > 40:
            self.chunks.append({"chapter": self.h2, "section": self.h3, "text": text})

    def handle_starttag(self, tag, attrs):
        if tag in _SKIP:
            self._skip += 1
            return
        if tag in ("h1", "h2", "h3"):
            self._flush()
            self._in_head, self._head = tag, []
        elif tag in _BLOCK:
            self._buf.append("\n")
        elif tag in ("td", "th"):
            self._buf.append(" | ")

    def handle_endtag(self, tag):
        if tag in _SKIP:
            self._skip = max(0, self._skip - 1)
            return
        if tag == self._in_head:
            title = re.sub(r"\s+", " ", "".join(self._head)).strip()
            if tag in ("h1", "h2"):
                self.h2, self.h3 = title, ""
            else:
                self.h3 = title
            self._in_head = None
        elif tag in ("li",):
            self._buf.append("\n")

    def handle_data(self, data):
        if self._skip:
            return
        if self._in_head:
            self._head.append(data)
        else:
            self._buf.append(data)

    def close(self):
        super().close()
        self._flush()


def _split(chunk):
    """Keep chunks small enough that several fit in one tool result."""
    text = chunk["text"]
    if len(text) <= MAX_CHARS:
        return [chunk]
    out, cur = [], ""
    for line in text.split("\n"):
        if cur and len(cur) + len(line) > MAX_CHARS:
            out.append({**chunk, "text": cur.strip()})
            cur = ""
        cur += line + "\n"
    if cur.strip():
        out.append({**chunk, "text": cur.strip()})
    return out


def build(html_path=MANUAL_HTML, out_path=CHUNKS_FILE):
    parser = _Flatten()
    parser.feed(Path(html_path).read_text(encoding="utf-8"))
    parser.close()
    chunks = [c for raw in parser.chunks for c in _split(raw)]
    for i, c in enumerate(chunks):
        c["id"] = i
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(chunks, ensure_ascii=False), encoding="utf-8")
    return chunks


# ---------------------------------------------------------------- BM25 ----

_WORD = re.compile(r"[a-z0-9][a-z0-9_\-\.]*[a-z0-9]|[a-z0-9]")
_STOP = set("the a an of and or to in on for is are be it this that with as by at from how do does what which can i my you your".split())


def _tokens(text):
    return [t for t in _WORD.findall(text.lower()) if t not in _STOP]


class Manual:
    def __init__(self, chunks=None):
        if chunks is None:
            chunks = self._load()
        self.chunks = chunks
        self._docs = []
        df = {}
        for c in chunks:
            # the heading words are worth more than body text: repeat them
            toks = _tokens(f"{c['chapter']} {c['section']} " * 3 + c["text"])
            tf = {}
            for t in toks:
                tf[t] = tf.get(t, 0) + 1
            self._docs.append((tf, len(toks)))
            for t in tf:
                df[t] = df.get(t, 0) + 1
        n = max(len(chunks), 1)
        self._idf = {t: math.log(1 + (n - d + 0.5) / (d + 0.5)) for t, d in df.items()}
        self._avg = (sum(l for _, l in self._docs) / n) or 1.0

    @staticmethod
    def _load():
        path = CHUNKS_FILE
        if getattr(sys, "frozen", False):
            path = Path(getattr(sys, "_MEIPASS", HERE)) / "tools" / "assistant" / "data" / "manual_chunks.json"
        if not path.is_file():
            return []
        return json.loads(path.read_text(encoding="utf-8"))

    def search(self, query, k=5, chapter=None):
        q = _tokens(query)
        if not q or not self.chunks:
            return []
        scored = []
        for c, (tf, ln) in zip(self.chunks, self._docs):
            if chapter and chapter.lower() not in c["chapter"].lower():
                continue
            s = 0.0
            for t in q:
                f = tf.get(t)
                if f:
                    s += self._idf.get(t, 0) * f * 2.2 / (f + 1.2 * (0.25 + 0.75 * ln / self._avg))
            if s > 0:
                scored.append((s, c))
        scored.sort(key=lambda x: -x[0])
        return [dict(c, score=round(s, 2)) for s, c in scored[:k]]

    def chapters(self):
        seen = []
        for c in self.chunks:
            if c["chapter"] and c["chapter"] not in seen:
                seen.append(c["chapter"])
        return seen


if __name__ == "__main__":
    ch = build()
    print(f"{len(ch)} chunks -> {CHUNKS_FILE}")
