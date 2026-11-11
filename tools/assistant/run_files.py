"""Read-only access to the files a user's TerraTools runs produced.

Every tool bridge announces its finished run through `resultReady`; the
assistant listens to those and remembers (a) the result JSON and (b) which
files it points to. The model can only touch files that were either named in
a result, live in the same tool output folder, or were attached by the user
through a file dialog. Nothing else on disk is reachable - the path check in
`_resolve` is the only door, and it is deliberately strict.
"""

import base64
import io
import json
import os
import tempfile
import time
from pathlib import Path

MAX_TEXT = 14000          # characters returned to the model per call
MAX_ROWS = 60
IMG_EDGE = 1400
TEXT_EXT = {".txt", ".csv", ".tsv", ".json", ".log", ".md", ".r", ".py", ".vcf"}
IMG_EXT = {".png", ".jpg", ".jpeg"}
XLS_EXT = {".xlsx", ".xlsm"}


def _output_roots():
    """Folders in which every file is a tool output (siblings are readable)."""
    roots = [Path(tempfile.gettempdir())]
    here = Path(__file__).resolve().parents[1]
    roots += list(here.glob("*/ui/results")) + list(here.glob("*/results")) + list(here.glob("*/backend/results"))
    roots.append(Path.home() / "Library" / "Application Support" / "TerraTools")
    return [r.resolve() for r in roots if r.exists()]


def _inside(p: Path, root: Path) -> bool:
    try:
        p.relative_to(root)
        return True
    except ValueError:
        return False


def _shrink(obj, depth=0):
    """Compact a result JSON: long lists / long strings are cut with a count."""
    if isinstance(obj, dict):
        return {k: _shrink(v, depth + 1) for k, v in list(obj.items())[:80]}
    if isinstance(obj, list):
        head = [_shrink(v, depth + 1) for v in obj[:12]]
        return head + [f"... {len(obj) - 12} more items"] if len(obj) > 12 else head
    if isinstance(obj, str) and len(obj) > 600:
        if obj.startswith("data:"):
            return f"<embedded image {len(obj)} chars>"
        return obj[:600] + f"... ({len(obj)} chars)"
    return obj


class RunStore:
    def __init__(self):
        self.runs = {}                 # tool -> {"time", "status", "result", "files": {abs path}}
        self.attached = {}             # abs path -> display name
        self._roots = None
        self.log = []                  # [(time, path)] of files given to the model, shown to the user

    # ------------------------------------------------------------ intake
    def note_result(self, tool, payload):
        try:
            result = json.loads(payload) if isinstance(payload, str) else payload
        except Exception:
            return
        if not isinstance(result, dict):
            return
        files = set()
        self._collect_paths(result, files)
        self.runs[tool] = {
            "time": time.time(), "status": result.get("status", "?"), "result": result,
            "files": files,
        }

    def _collect_paths(self, obj, out, depth=0):
        if depth > 6:
            return
        if isinstance(obj, dict):
            for v in obj.values():
                self._collect_paths(v, out, depth + 1)
        elif isinstance(obj, list):
            for v in obj[:200]:
                self._collect_paths(v, out, depth + 1)
        elif isinstance(obj, str) and 3 < len(obj) < 600 and not obj.startswith(("data:", "http")):
            try:
                p = Path(obj)
                if p.is_absolute() and p.is_file():
                    out.add(str(p.resolve()))
            except (OSError, ValueError):
                pass

    def attach(self, path):
        p = Path(path).resolve()
        if p.is_file():
            self.attached[str(p)] = p.name
            return p.name
        return None

    # ------------------------------------------------------------ policy
    def _roots_cached(self):
        if self._roots is None:
            self._roots = _output_roots()
        return self._roots

    def _in_output_dir(self, fp: Path) -> bool:
        """True when fp sits in a folder *below* a tool-output root (so its siblings
        are all outputs). A file lying directly in %TEMP% does not qualify: its
        siblings are other programs' files."""
        return any(_inside(fp.parent, r) and fp.parent != r for r in self._roots_cached())

    def _allowed(self, p: Path) -> bool:
        s = str(p)
        if s in self.attached:
            return True
        for run in list(self.runs.values()):
            if s in run["files"]:
                return True
        # siblings of a recorded file are fine only when its folder is a tool output folder
        for run in list(self.runs.values()):
            for f in run["files"]:
                fp = Path(f)
                if p.parent == fp.parent and self._in_output_dir(fp):
                    return True
        return False

    def _resolve(self, tool, name):
        """name is a file name from list_runs, or an absolute path that is allowed."""
        name = (name or "").strip()
        if not name:
            raise ValueError("No file given.")
        cands = []
        if os.path.isabs(name):
            cands.append(Path(name))
        run = self.runs.get(tool or "")
        pools = [run] if run else list(self.runs.values())
        for r in pools:
            for f in r["files"]:
                if Path(f).name.lower() == name.lower():
                    cands.append(Path(f))
                    d = Path(f).parent
                    cands.append(d / name)
        for a in self.attached:
            if Path(a).name.lower() == name.lower():
                cands.append(Path(a))
        for c in cands:
            try:
                c = c.resolve()
            except OSError:
                continue
            if c.is_file() and self._allowed(c):
                return c
        raise ValueError(f"'{name}' is not one of the files from the user's runs. Call list_runs to see the available files.")

    # ------------------------------------------------------------ tools
    def list_runs(self):
        if not self.runs and not self.attached:
            return {"runs": [], "note": "The user has not finished any analysis in this session yet, and has attached no files."}
        out = []
        for tool, r in sorted(list(self.runs.items()), key=lambda kv: -kv[1]["time"]):
            names = {}
            for f in r["files"]:
                p = Path(f)
                names[p.name] = p
                if self._in_output_dir(p):
                    for sib in p.parent.iterdir():
                        if sib.is_file() and not sib.name.startswith("."):
                            names.setdefault(sib.name, sib)
            files = []
            for n, p in sorted(names.items()):
                try:
                    files.append({"name": n, "kb": round(p.stat().st_size / 1024, 1)})
                except OSError:
                    pass
            out.append({
                "tool": tool, "status": r["status"],
                "finished": time.strftime("%H:%M:%S", time.localtime(r["time"])),
                "files": files[:60],
            })
        att = [{"name": n} for n in self.attached.values()]
        return {"runs": out, "attached_by_user": att}

    def get_run_result(self, tool):
        r = self.runs.get(tool)
        if not r:
            return {"error": f"No finished run for '{tool}' in this session. Known: {sorted(self.runs)}"}
        text = json.dumps(_shrink(r["result"]), ensure_ascii=False, default=str)
        self.log.append((time.time(), f"{tool} (run summary)"))
        return {"tool": tool, "status": r["status"], "result": json.loads(text) if len(text) < MAX_TEXT else text[:MAX_TEXT] + "…"}

    def read_file(self, name, tool=None, offset=0, limit=20, sheet=None):
        """Return a list of content blocks (text, or text + image)."""
        try:
            p = self._resolve(tool, name)
        except ValueError as e:
            return [{"type": "text", "text": f"Error: {e}"}], True
        self.log.append((time.time(), p.name))
        ext = p.suffix.lower()
        try:
            if ext in IMG_EXT:
                return self._image(p), False
            if ext in XLS_EXT:
                return [{"type": "text", "text": self._excel(p, sheet, offset, limit)}], False
            if ext in (".csv", ".tsv") or (ext == ".txt" and _looks_tabular(p)):
                return [{"type": "text", "text": self._table(p, offset, limit)}], False
            if ext in TEXT_EXT:
                return [{"type": "text", "text": _read_text(p)}], False
            return [{"type": "text", "text": f"{p.name}: {p.stat().st_size} bytes, type '{ext}' - cannot be shown as text."}], False
        except Exception as e:
            return [{"type": "text", "text": f"Error reading {p.name}: {type(e).__name__}: {e}"}], True

    # ------------------------------------------------------------ readers
    @staticmethod
    def _image(p):
        from PIL import Image
        with Image.open(p) as im:
            im = im.convert("RGB")
            w, h = im.size
            im.thumbnail((IMG_EDGE, IMG_EDGE))
            buf = io.BytesIO()
            im.save(buf, "PNG", optimize=True)
            media = "image/png"
            if buf.tell() > 1_300_000:
                buf = io.BytesIO()
                im.save(buf, "JPEG", quality=82)
                media = "image/jpeg"
        data = base64.standard_b64encode(buf.getvalue()).decode()
        return [
            {"type": "text", "text": f"Image {p.name} ({w}x{h} px, shown at {im.size[0]}x{im.size[1]})."},
            {"type": "image", "source": {"type": "base64", "media_type": media, "data": data}},
        ]

    @staticmethod
    def _table(p, offset, limit):
        import pandas as pd
        limit = max(1, min(int(limit or 20), MAX_ROWS))
        offset = max(0, int(offset or 0))
        sep = "\t" if p.suffix.lower() == ".tsv" else None
        df = pd.read_csv(p, sep=sep, engine="python", dtype=object, on_bad_lines="skip", encoding="utf-8-sig", encoding_errors="replace")
        return _describe_frame(p.name, df, offset, limit)

    @staticmethod
    def _excel(p, sheet, offset, limit):
        import pandas as pd
        limit = max(1, min(int(limit or 20), MAX_ROWS))
        offset = max(0, int(offset or 0))
        xl = pd.ExcelFile(p)
        names = xl.sheet_names
        if not sheet:
            if len(names) > 1:
                parts = []
                for n in names[:12]:
                    df = xl.parse(n, dtype=object)
                    parts.append(f"sheet '{n}': {len(df)} rows x {len(df.columns)} cols; columns: {', '.join(map(str, df.columns[:14]))}")
                return (f"{p.name} has {len(names)} sheets. Call read_file again with sheet=<name> to see rows.\n" + "\n".join(parts))[:MAX_TEXT]
            sheet = names[0]
        if sheet not in names:
            return f"No sheet '{sheet}'. Sheets: {names}"
        return _describe_frame(f"{p.name} [{sheet}]", xl.parse(sheet, dtype=object), offset, limit)


def _looks_tabular(p):
    try:
        head = p.read_text(errors="replace")[:2000]
    except OSError:
        return False
    first = head.splitlines()[0] if head else ""
    return first.count(",") >= 2 or first.count("\t") >= 2


def _read_text(p):
    raw = p.read_bytes()[: MAX_TEXT * 4]
    text = raw.decode("utf-8", errors="replace")
    size = p.stat().st_size
    if len(text) > MAX_TEXT:
        text = text[:MAX_TEXT] + f"\n... (truncated; file is {size} bytes)"
    return text


def _describe_frame(label, df, offset, limit):
    import pandas as pd
    n_rows, n_cols = df.shape
    lines = [f"{label}: {n_rows} rows x {n_cols} columns", "columns: " + ", ".join(map(str, df.columns[:60]))]
    num = df.apply(pd.to_numeric, errors="coerce")
    stats = []
    for c in df.columns[:40]:
        col = num[c]
        if col.notna().sum() >= max(3, 0.8 * df[c].notna().sum()) and df[c].notna().sum() > 0:
            stats.append(f"{c}: min {col.min():.4g}, median {col.median():.4g}, max {col.max():.4g}, missing {int(df[c].isna().sum())}")
    if stats:
        lines.append("numeric summary:\n  " + "\n  ".join(stats[:25]))
    part = df.iloc[offset: offset + limit]
    lines.append(f"rows {offset + 1}-{offset + len(part)} (use offset/limit for more, max {MAX_ROWS} per call):")
    lines.append(part.to_csv(index=False))
    return "\n".join(lines)[:MAX_TEXT]

