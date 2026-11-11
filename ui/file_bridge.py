"""
file_bridge.py — Saving result images to a place the user picks.

One bridge for the whole app ("terraFiles" on the shared web channel), so
every tool's images can be downloaded the same way without each tool
needing its own save code:

    save_image(source, suggested_name)  -> Save As dialog, copies one image
    save_images(items_json)             -> Folder dialog, copies many images

`source` is a local file path, a file:// URL, or a data: URL.
"""

import base64
import json
import re
import shutil
from pathlib import Path
from urllib.parse import unquote, urlparse

from PySide6.QtCore import QObject, QSettings, Slot
from PySide6.QtWidgets import QFileDialog

_BAD_NAME_CHARS = re.compile(r'[<>:"/\\|?*\x00-\x1f]+')


def _clean_name(name, default="image.png"):
    name = _BAD_NAME_CHARS.sub("_", str(name or "")).strip(" ._")
    if not name:
        name = default
    if not Path(name).suffix:
        name += ".png"
    return name


def _local_path(source):
    """file:// URL (with or without a cache-busting ?query) -> filesystem path."""
    source = str(source or "")
    if source.lower().startswith("file:"):
        path = unquote(urlparse(source).path)
        if re.match(r"^/[A-Za-z]:", path):      # file:///C:/x -> C:/x
            path = path[1:]
        return path
    return source.split("?t=")[0] if re.search(r"\?t=\d+$", source) else source


def _write_to(source, dest):
    """Copy a file, or decode a data: URL, to dest."""
    source = str(source or "")
    dest = Path(dest)
    if source.startswith("data:"):
        header, _, payload = source.partition(",")
        raw = base64.b64decode(payload) if ";base64" in header else unquote(payload).encode("utf-8")
        dest.write_bytes(raw)
        return
    src = Path(_local_path(source))
    if not src.is_file():
        raise FileNotFoundError(f"Image file not found: {src}")
    shutil.copy2(src, dest)


def _unique(folder, name):
    candidate = Path(folder) / name
    stem, suffix, n = candidate.stem, candidate.suffix, 2
    while candidate.exists():
        candidate = Path(folder) / f"{stem} ({n}){suffix}"
        n += 1
    return candidate


class FileBridge(QObject):

    def __init__(self, parent_window=None):
        super().__init__(parent_window)
        self._window = parent_window
        self._settings = QSettings("IRRI", "TerraTools")

    def _last_dir(self):
        d = self._settings.value("last_save_dir", "", type=str)
        return d if d and Path(d).is_dir() else str(Path.home() / "Pictures") if (Path.home() / "Pictures").is_dir() else str(Path.home())

    def _remember(self, folder):
        self._settings.setValue("last_save_dir", str(folder))

    @Slot(str, str, result=str)
    def save_image(self, source, suggested_name):
        """Save one image. Returns the saved path, "" if cancelled, or
        "ERROR: <reason>"."""
        try:
            name = _clean_name(suggested_name)
            dest, _ = QFileDialog.getSaveFileName(
                self._window, "Save image", str(Path(self._last_dir()) / name),
                "PNG image (*.png);;All files (*)",
            )
            if not dest:
                return ""
            _write_to(source, dest)
            self._remember(Path(dest).parent)
            print(f"[FILES] saved {source!r} -> {dest}", flush=True)
            return dest
        except Exception as e:
            import traceback
            traceback.print_exc()
            return f"ERROR: {e}"

    @Slot(str, result=str)
    def save_images(self, items_json):
        """Save many images into one folder. items_json = [{"src","name"}, ...].
        Returns JSON: {"status": "success"|"cancelled"|"error", "saved", "folder", "failed"}."""
        try:
            items = json.loads(items_json)
        except Exception as e:
            return json.dumps({"status": "error", "message": f"Bad request: {e}"})
        if not items:
            return json.dumps({"status": "error", "message": "No images to save."})

        folder = QFileDialog.getExistingDirectory(self._window, "Choose a folder for the images", self._last_dir())
        if not folder:
            return json.dumps({"status": "cancelled"})

        saved, failed, used = 0, [], set()
        for item in items:
            try:
                name = _clean_name(item.get("name"))
                if name.lower() in used:           # two images with the same label
                    name = _unique(folder, name).name
                dest = _unique(folder, name)
                _write_to(item.get("src"), dest)
                used.add(dest.name.lower())
                saved += 1
            except Exception as e:
                failed.append({"name": item.get("name"), "error": str(e)})
        self._remember(folder)
        print(f"[FILES] saved {saved}/{len(items)} images to {folder}", flush=True)
        return json.dumps({"status": "success" if saved else "error", "saved": saved,
                           "folder": folder, "failed": failed,
                           "message": None if saved else "None of the images could be saved."})
