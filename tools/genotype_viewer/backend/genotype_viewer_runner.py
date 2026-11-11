"""
genotype_viewer_runner.py — runs GenoSee (backend/genosee, https://github.com/hashimotoshumpei/GenoSee) for the
Genotype Viewer screen.

GenoSee.py's main() is a command line; this module does what main() does (load the species and colour databases,
read the CSV or VCF (the first three columns are taken by position), call create_normal_plot / create_comparison_plot / create_zoomed_plot
with the same checks) and returns the written images. The drawing itself is GenoSee's own code.
"""

import contextlib
import gzip
import io
import json
import re
import shutil
import tempfile
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")   # images are written to files; no window, and safe on a worker thread

import pandas as pd

from tools.genotype_viewer.backend.genosee import plotting as genosee_plotting
from tools.genotype_viewer.backend.genosee.plotting import create_comparison_plot, create_normal_plot, create_zoomed_plot
from tools.genotype_viewer.backend.genosee.utilities import get_value_from_database
from tools.genotype_viewer.config import GenotypeViewerConfig

GENOSEE_DIR = Path(__file__).resolve().parent / "genosee"
COLOR_JSON = GENOSEE_DIR / "color_set.json"
# palettes added by TerraTools, written in the format of GenoSee's color_set.json
EXTRA_COLOR_JSON = Path(__file__).resolve().parent / "terratools_color_set.json"
SPECIES_JSON = GENOSEE_DIR / "chromosome_length_database.json"
# species added from a genome FASTA are kept here, so the installed (read-only) database is never changed
USER_SPECIES_JSON = Path.home() / ".terratools" / "genotype_viewer" / "chromosome_length_database.json"


def species_db_path() -> Path:
    return USER_SPECIES_JSON if USER_SPECIES_JSON.is_file() else SPECIES_JSON


def list_species() -> list:
    with open(species_db_path(), "r", encoding="utf-8") as f:
        return list(json.load(f))


def list_palettes() -> dict:
    with open(COLOR_JSON, "r", encoding="utf-8") as f:
        palettes = json.load(f)
    with open(EXTRA_COLOR_JSON, "r", encoding="utf-8") as f:
        palettes.update(json.load(f))
    return palettes


def add_species(fasta_file: str, name: str) -> dict:
    """GenoSee's add_chromosome_length_to_database.py: read the chromosome lengths of a genome FASTA into the database."""
    try:
        if not fasta_file or not Path(fasta_file).is_file():
            return {"status": "error", "message": "Please choose the genome FASTA file."}
        if not name.strip():
            return {"status": "error", "message": "Please enter the species name."}
        from tools.genotype_viewer.backend.genosee.add_chromosome_length_to_database import add_chromosome_lengths_to_json
        if not USER_SPECIES_JSON.is_file():
            USER_SPECIES_JSON.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(SPECIES_JSON, USER_SPECIES_JSON)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            add_chromosome_lengths_to_json(fasta_file, str(USER_SPECIES_JSON), name.strip().replace(" ", "_"))
        return {"status": "success", "species": list_species(), "added": name.strip().replace(" ", "_")}
    except Exception as e:
        return {"status": "error", "message": f"{e}\n\n{traceback.format_exc()}"}


def is_vcf(path: str) -> bool:
    n = str(path).lower()
    return n.endswith(".vcf") or n.endswith(".vcf.gz")


def _delimiter(path: str) -> str:
    with open(path, "r", encoding="utf-8-sig", errors="replace") as f:
        head = f.readline()
    return max([",", "\t", ";"], key=head.count)


def _vcf_samples(vcf_file: str) -> list:
    opener = gzip.open if str(vcf_file).lower().endswith(".gz") else open
    with opener(vcf_file, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#CHROM"):
                return line.rstrip("\r\n").split("\t")[9:]
            if not line.startswith("#"):
                break
    raise ValueError("This VCF file has no #CHROM header line.")


def vcf_to_csv(vcf_file: str, marker_names: bool = True) -> str:
    """GenoSee's make_input_file_from_VCF.py (PyVCF): one row per site, chr, marker_name, pos, then each sample's GT.
    marker_names=True names the markers CHROM_POS (GenoSee's --add_marker_names); otherwise the VCF ID column is used."""
    from tools.genotype_viewer.backend.genosee.make_input_file_from_VCF import convert_vcf_to_csv
    work = Path(tempfile.mkdtemp(prefix="terratools_genotype_viewer_vcf_"))
    src = Path(vcf_file)
    if src.name.lower().endswith(".gz"):          # GenoSee's script opens a plain-text VCF: unpack a .vcf.gz first
        plain = work / src.name[:-3]
        with gzip.open(src, "rb") as zin, open(plain, "wb") as zout:
            shutil.copyfileobj(zin, zout)
        src = plain
    out = work / (src.stem + ".csv")
    convert_vcf_to_csv(str(src), str(out), marker_names)
    return str(out)


def _norm_chr(name) -> str:
    return re.sub(r"^chr", "", str(name).strip(), flags=re.I).lstrip("0").lower() or "0"


def read_genotypes(path: str, vcf_marker_names: bool = True) -> pd.DataFrame:
    """The input file as GenoSee wants it. The columns are taken by position, whatever they are called: the first is the
    chromosome, then the marker name, then the position, then one column per sample. A VCF is converted first."""
    if is_vcf(path):
        path = vcf_to_csv(path, vcf_marker_names)
    data = pd.read_csv(path, sep=_delimiter(path), dtype=str, encoding="utf-8-sig")
    if data.shape[1] < 4:
        raise ValueError("The file needs at least four columns: chromosome, marker name, position, then one column per sample.")
    data.columns = ["chr", "marker_name", "pos"] + [str(c).strip() for c in data.columns[3:]]
    data["chr"] = data["chr"].astype(str).str.strip()
    pos = pd.to_numeric(data["pos"], errors="coerce")
    bad = int(pos.isna().sum())
    if bad:
        raise ValueError(f"The third column must be the position (a number). {bad} row(s) have none, "
                         f"for example '{data.loc[pos.isna(), 'pos'].iloc[0]}'.")
    data["pos"] = pos.astype(int)
    return data


def sample_names(path: str) -> list:
    """The sample columns (everything after the first three) of a genotype file."""
    if is_vcf(path):
        return _vcf_samples(path)
    cols = list(pd.read_csv(path, sep=_delimiter(path), nrows=0, encoding="utf-8-sig").columns)
    if len(cols) < 4:
        raise ValueError("The file needs at least four columns: chromosome, marker name, position, then one column per sample.")
    return [str(c).strip() for c in cols[3:]]


class GenotypeViewerEngine:
    def __init__(self, config: GenotypeViewerConfig):
        self.config = config

    def run(self) -> dict:
        log = io.StringIO()
        try:
            return self._run(log)
        except SystemExit:
            # GenoSee's plotting code ends the program (sys.exit) on a genotype it cannot draw; it has printed why
            return {"status": "error", "message": log.getvalue().strip()[-1500:] or "GenoSee could not draw this file."}
        except Exception as e:
            return {"status": "error", "message": f"{e}\n\nDebug info:\n{traceback.format_exc()}"}

    def _run(self, log) -> dict:
        c = self.config
        if not c.input_file or not Path(c.input_file).is_file():
            return {"status": "error", "message": "Please select the input CSV file."}

        chrs_dict = get_value_from_database(str(species_db_path()), c.species) if c.species else None
        if not chrs_dict:
            return {"status": "error", "message": f"Species '{c.species}' is not in the chromosome length database. "
                                                  "Choose another species or add it from its genome FASTA."}
        color_dict = list_palettes().get(c.color_set)
        if not color_dict:
            return {"status": "error", "message": f"Colour set '{c.color_set}' does not exist."}

        data = read_genotypes(c.input_file, c.vcf_marker_names)

        # the chromosome names of the file may be written differently (1, chr1, Chr01 ...): use the species' own names
        by_norm = {}
        for name in chrs_dict:
            by_norm.setdefault(_norm_chr(name), []).append(name)
        by_norm = {k: v[0] for k, v in by_norm.items() if len(v) == 1}
        data["chr"] = data["chr"].map(lambda v: v if v in chrs_dict else by_norm.get(_norm_chr(v), v))
        on_species = data["chr"].isin(chrs_dict)
        if not on_species.any():
            return {"status": "error", "message": f"None of the chromosome names in the first column ({', '.join(map(str, data['chr'].unique()[:6]))}) "
                                                  f"match the species {c.species} ({', '.join(chrs_dict)}). Choose the right species or add it from its genome FASTA."}
        n_off = int((~on_species).sum())
        data = data[on_species]
        data.sort_values(["chr", "pos"], inplace=True)

        if c.samples:
            missing = [s for s in c.samples if s not in data.columns[3:]]
            if missing:
                return {"status": "error", "message": f"Sample(s) not in the file: {', '.join(missing[:5])}"}
            data = data[list(data.columns[:3]) + [s for s in data.columns[3:] if s in c.samples]]

        if c.drawing_mode == "zoomed":
            if not c.chr or not c.end:
                return {"status": "error", "message": "The zoomed view needs a chromosome, a start and an end position."}
            if c.chr not in chrs_dict:
                return {"status": "error", "message": f"Chromosome name {c.chr} doesn't match the species "
                                                      f"(it has: {', '.join(chrs_dict)})."}
            if c.start > c.end:
                return {"status": "error", "message": "Start position is larger than end position."}
            if c.end > chrs_dict[c.chr]:
                return {"status": "error", "message": f"End position is larger than the chromosome length ({chrs_dict[c.chr]})."}
        elif c.drawing_mode not in ("normal", "compare"):
            return {"status": "error", "message": "Invalid drawing mode."}

        # the legend texts: GenoSee's own words unless the user named the parents / alleles
        for key in genosee_plotting.LEGEND_LABELS:
            genosee_plotting.LEGEND_LABELS[key] = c.legend_labels.get(key) or key
        genosee_plotting.GUIDE_LINES = c.guide_lines
        genosee_plotting.GRAPHIC = c.graphic
        outdir = Path(tempfile.mkdtemp(prefix="terratools_genotype_viewer_"))
        out = str(outdir)
        with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            if c.drawing_mode == "normal":
                for column in range(3, len(data.columns)):
                    create_normal_plot(data, chrs_dict, color_dict, c.coloring_mode, c.fill, c.display_marker_names,
                                       data.columns[column], column, out, c.dpi, c.pdf)
            elif c.drawing_mode == "compare":
                create_comparison_plot(data, chrs_dict, color_dict, c.coloring_mode, c.fill, c.display_marker_names,
                                       out, c.dpi, c.pdf)
            else:
                create_zoomed_plot(data, chrs_dict, color_dict, c.coloring_mode, c.fill, c.display_marker_names,
                                   c.chr, c.start, c.end, out, c.dpi, c.pdf)

        images = []
        for png in sorted(outdir.glob("*.png")):
            pdf = png.with_suffix(".pdf")
            images.append({"name": png.stem, "png": str(png), "pdf": str(pdf) if pdf.is_file() else None})
        if not images:
            return {"status": "error", "message": "No image was written.\n" + log.getvalue().strip()[-800:]}
        return {
            "status": "success",
            "mode": c.drawing_mode,
            "images": images,
            "n_images": len(images),
            "n_samples": len(data.columns) - 3,
            "n_markers": int(len(data)),
            "n_chromosomes": int(data["chr"].nunique()),
            "n_off_species": n_off,
            "outdir": out,
        }
