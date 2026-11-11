"""
save_dialog.py — one Save As dialog for every tool's result downloads.

The dialog filter comes from the suggested file name (an .xlsx offers
"Excel workbook (*.xlsx)") and the extension is added back if the user types a
name without one. Without this a report saved as "results" has no extension, and
Excel / Windows then call the file invalid.
"""

from pathlib import Path

from PySide6.QtWidgets import QFileDialog

_KINDS = {
    ".xlsx": "Excel workbook", ".csv": "CSV file", ".zip": "ZIP archive", ".png": "PNG image",
    ".html": "HTML file", ".txt": "Text file", ".pdf": "PDF document", ".json": "JSON file",
    ".tsv": "TSV file", ".docx": "Word document",
}


def save_dialog(parent, suggested_name, title="Save File"):
    """Returns (path, selected_filter); path is "" when cancelled."""
    suffix = Path(str(suggested_name)).suffix.lower()
    kind = _KINDS.get(suffix, f"{suffix[1:].upper()} file" if suffix else "")
    flt = f"{kind} (*{suffix});;All Files (*)" if suffix else "All Files (*)"
    path, chosen = QFileDialog.getSaveFileName(parent, title, str(suggested_name), flt)
    if path and suffix and not Path(path).suffix:
        path += suffix
    return path, chosen
