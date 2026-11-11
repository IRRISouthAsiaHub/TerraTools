"""
json_safe.py — json.dumps that the web page can always parse.

Python writes NaN / Infinity for non-finite floats, which is not valid JSON: JSON.parse() in the page
then throws and the results of a perfectly good analysis never appear (for example a correlation of
constant predictions is NaN). Here every non-finite number becomes null, and numpy numbers become
plain Python ones.
"""
import json
import math


def clean(obj):
    if isinstance(obj, dict):
        return {str(k): clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean(v) for v in obj]
    if hasattr(obj, "item") and not isinstance(obj, (str, bytes)) and getattr(obj, "shape", None) == ():
        obj = obj.item()                    # numpy scalar
    if hasattr(obj, "tolist") and not isinstance(obj, (str, bytes)):
        return clean(obj.tolist())          # numpy array
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    return obj


def dumps(obj, **kw):
    return json.dumps(clean(obj), allow_nan=False, default=str, **kw)
