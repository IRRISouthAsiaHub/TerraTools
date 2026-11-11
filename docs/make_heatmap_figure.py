"""Draws the balloon-heatmap figure of the manual with the real backend (the R script) from the demo 4K file.
Run after make_screenshots.py:  python docs/make_heatmap_figure.py"""
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from PIL import Image  # noqa: E402

from tools.favourable_allele.backend import balloon_runner as runner  # noqa: E402
from tools.favourable_allele.backend import matcher  # noqa: E402

DEMO = ROOT / "Demo_Input_Files" / "18_Favourable_Allele"
GENES = ["DROT1", "DRG9", "Xa26|Xa3", "Pi25|Pid3", "HIS1", "Wx|qFC6", "TOND1", "OsNRAMP5", "GS5",
         "NAL1|qFLW4|LVPA4", "sd1|GA20ox2|OsGA20ox2", "Gn1a|OsCKX2"]
d = matcher.load(str(DEMO / "input" / "rice_4K_haplopanel_demo.csv"))
if d["status"] != "success":
    sys.exit(d["message"])
order = {g: i for i, g in enumerate(GENES)}
idx = sorted((i for i, s in enumerate(d["snps"]) if (s["name"] or s["gene"]) in order), key=lambda i: (order[d["snps"][i]["name"] or d["snps"][i]["gene"]], i))
vi = list(range(len(d["varieties"])))
res = runner.draw({
    "snps": [d["snps"][i] for i in idx],
    "varieties": [{"name": d["varieties"][j]["name"], "group": "Variety"} for j in vi],
    "calls": ["".join(d["calls"][j][i] for i in idx) for j in vi],
    "group_label": "Lines",
    "dpi": 300,
})
if res["status"] != "success":
    sys.exit(res["message"])
out = ROOT / "docs" / "manual" / "shots" / "favallele_heatmap.jpg"
im = Image.open(res["image"]).convert("RGB")
im.thumbnail((2000, 2000))
im.save(out, quality=90)
exp = DEMO / "expected_output"
exp.mkdir(exist_ok=True)
shutil.copy2(res["image"], exp / "gene_balloon.png")
for f in res["files"]:
    shutil.copy2(f["path"], exp / f["name"])
print("wrote", out, im.size)
