"""JSON helper of the In silico BSA bridge."""

import numpy as np


def json_safe(obj):
    """Result -> something json.dumps writes as VALID JSON: NaN / inf (empty cells, windows that were not tested) become null, numpy numbers
    become Python numbers. Python's own dumps would write the bare word NaN, which JSON.parse in the page rejects."""
    if isinstance(obj, dict):
        return {str(k): json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [json_safe(v) for v in obj]
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        return float(obj) if np.isfinite(obj) else None
    return obj
