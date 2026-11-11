from dataclasses import dataclass, field


@dataclass
class GenotypeViewerConfig:
    """The options of GenoSee's command line (GenoSee.py --help), plus an optional list of samples to draw."""
    input_file: str
    species: str = ""
    drawing_mode: str = "normal"          # normal | compare | zoomed
    coloring_mode: str = "2-color"        # 2-color | 3-color
    fill: str = "on"                      # on | off
    display_marker_names: str = "on"      # on | off
    color_set: str = "normal"             # a key of color_set.json
    chr: str = ""                         # zoomed mode: chromosome, start and end of the region
    start: int = 0
    end: int = 0
    dpi: int = 200
    pdf: bool = False
    samples: list = field(default_factory=list)   # empty = every sample column of the file
    legend_labels: dict = field(default_factory=dict)   # REF / ALT / HET / MISSING -> text of the legend (empty = GenoSee's words)
    graphic: bool = True                  # graphic look: rounded, shaded chromosomes without the black lines inside
    guide_lines: bool = True              # compare / zoomed: the dashed line from each marker name (needs marker names)
    vcf_marker_names: bool = True         # VCF input: name the markers CHROM_POS (else the VCF ID column)

    @staticmethod
    def from_dict(d: dict) -> "GenotypeViewerConfig":
        return GenotypeViewerConfig(
            input_file=d.get("input_file", ""),
            species=d.get("species", "") or "",
            drawing_mode=d.get("drawing_mode", "normal"),
            coloring_mode=d.get("coloring_mode", "2-color"),
            fill="on" if str(d.get("fill", "on")).lower() in ("on", "true", "1") else "off",
            display_marker_names="on" if str(d.get("display_marker_names", "on")).lower() in ("on", "true", "1") else "off",
            color_set=d.get("color_set", "normal"),
            chr=str(d.get("chr", "") or ""),
            start=int(d.get("start") or 0),
            end=int(d.get("end") or 0),
            dpi=int(d.get("dpi", 200)),
            pdf=bool(d.get("pdf", False)),
            samples=[str(s) for s in (d.get("samples") or [])],
            legend_labels={k: str(v).strip() for k, v in (d.get("legend_labels") or {}).items() if str(v).strip()},
            graphic=bool(d.get("graphic", True)),
            guide_lines=bool(d.get("guide_lines", True)),
            vcf_marker_names=bool(d.get("vcf_marker_names", True)),
        )
