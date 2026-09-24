# -*- coding: utf-8 -*-
"""
heterosis_preview.py

Lightweight CSV preview used by the GUI process directly (file
upload -> quick "N rows x N cols" preview). Deliberately has NO
TensorFlow / Keras import — only pandas/numpy — so it is safe to
call from the main Qt process without pulling TensorFlow's bundled
LLVM into the same process as QtWebEngine (see heterosis_process.py
for why that matters).
"""

import os

import numpy as np
import pandas as pd


def load_preview(path, n_rows=10):

    try:

        df = pd.read_csv(path)

        preview = df.head(n_rows).replace({np.nan: None})

        return {
            "status": "success",
            "filename": os.path.basename(path),
            "n_rows": int(df.shape[0]),
            "n_cols": int(df.shape[1]),
            "preview_columns": list(df.columns),
            "preview_rows": preview.to_dict(orient="records")
        }

    except Exception as e:

        return {
            "status": "error",
            "message": str(e)
        }