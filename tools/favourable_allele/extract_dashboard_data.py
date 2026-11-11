"""
Extract the trait-gene annotation of the Favourable Allele dashboard into tools/favourable_allele/backend/data/trait_genes_4k.json.

The source is the standalone dashboard page (Dashboard_inputfile_embedded.html), which carries the annotation (ANNOT_DB: gene, trait, marker, chromosome,
position, REF, ALT, favourable allele) and the list of genes without a scorable marker (LOST_GENES). Only the annotation is kept: the genotypes of the
dashboard's own panels are not part of the app, the user uploads a 4K haplopanel file that is matched against this list (backend/matcher.py).

    python tools/favourable_allele/extract_dashboard_data.py path/to/Dashboard_inputfile_embedded.html
"""
import json
import sys
from pathlib import Path

src = Path(sys.argv[1])
out = Path(__file__).resolve().parent / "backend" / "data" / "trait_genes_4k.json"

lines = src.read_text(encoding="utf-8").split("\n")


def grab(name):
    for ln in lines:
        if ln.startswith(f"const {name} = "):
            return json.loads(ln[len(f"const {name} = "):].rstrip().rstrip(";"))
    raise SystemExit(f"{name} not found")


annot, gt, lost = grab("ANNOT_DB"), grab("GT_DATA"), grab("LOST_GENES")
typed = set(gt["markers"])      # markers that are on the 4K haplopanel
snps = [a for a in annot if a["m"] in typed]
data = {
    "name": "Rice 4K haplopanel - trait-gene markers",
    "snps": [{"gene": a["g"], "name": a["n"], "trait": a["t"], "sub": a["s"], "marker": a["m"], "chrom": a["c"],
              "pos": a["p"], "ref": a["r"], "alt": a["a"], "fav": a["f"], "unfav": a["u"]} for a in snps],
    "lostGenes": lost,
}
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(data, separators=(",", ":")), encoding="utf-8")
print(f"{len(snps)} SNPs in {len({a['g'] for a in snps})} genes -> {out}")
