from pathlib import Path
from typing import Optional

from PySide6.QtCore import QUrl

_DIST = Path(__file__).resolve().parent.parent / "frontend" / "dist" / "index.html"


def frontend_url(route: str = "") -> Optional[QUrl]:
    """URL of the built React frontend (frontend/dist/index.html) for a route,
    or None if it hasn't been built (`npm run build` in frontend/), in which
    case callers fall back to the original standalone HTML page.

    route "" is the launcher dashboard; "gwas" / "clustering" / ... select a
    tool page (see frontend/src/App.jsx)."""
    if not _DIST.is_file():
        return None
    url = QUrl.fromLocalFile(str(_DIST))
    if route:
        url.setFragment(f"/{route}")
    return url
