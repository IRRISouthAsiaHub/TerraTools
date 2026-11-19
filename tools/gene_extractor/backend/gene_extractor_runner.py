from pathlib import Path


def _error(message):
    return {"status": "error", "message": message}


# ============================================================
# PUBLIC ENTRYPOINT
#
# Mirrors the shape of gene_runner.run_gene(config, results_dir) from
# the gene mapping tool, so the Qt bridge/worker pattern is identical.
# ============================================================

def run_extraction(config, results_dir):
    """
    config: {
        "marker_path": str,   # the QTL table of QTL Mapping (qtl_significant_*.csv, positions in cM)
        "gff_path": str,      # GFF3 gene annotation
        "positions_path": str,  # REQUIRED: the marker position file with the base-pair position of the markers of the QTL table (a table
                              # Marker, Chromosome, bp; also the 'Physical Map' sheet of a Linkage Mapping workbook) or a VCF (CHROM, POS, ID)
        "interval": str,      # optional, QTL table in cM: "ci" (confidence interval, default) or "flank" (left to right marker)
        "fasta_path": str,    # optional reference FASTA: adds each gene's sequence
        "match": str          # optional: "inside" (gene completely inside the region, default) or "overlap"
    }

    results_dir: Path or str, directory Genes_in_gene_regions.xlsx, Gene_ID_list.txt (and Gene_sequences.fasta)
    are written to.

    Returns a dict:
        {
          "status": "success", "input_kind": "interval",
          "genes": [ {Chromosome, Start, Stop, Strand, Gene_id, Description, Region_ID, Region_Start, Region_End,
                       Name, Feature, Source, Length_bp,
                       <other GFF attributes>}, ... ],
          "columns": [...], "region_summary": [ {Region_ID, Chromosome, Region_Start, Region_End, Genes} ],
          "region_columns": [...], "gene_ids": int,
          "genes_file": "<path>", "gene_list_file": "<path>", "fasta_file": "<path>" (with a FASTA), "warnings": [...]
        }
    or:
        { "status": "error", "message": "..." }
    """
    from tools.gene_extractor.backend.gff_genes import GeneExtractionError, run_gff_extraction

    for key, label in (("marker_path", "QTL table"), ("gff_path", "Annotation (GFF3)")):
        if not config.get(key):
            return _error("Please upload the QTL table, the marker position file (or a VCF) and the annotation (GFF3) file.")
        if not Path(config[key]).exists():
            return _error(f"{label} file not found: {config[key]}")
    if not config.get("positions_path"):
        return _error("Please select the marker position file (Marker, Chromosome, bp) or a VCF file: it gives the base-pair position of the QTL markers.")
    if not Path(config["positions_path"]).exists():
        return _error(f"Marker position file / VCF not found: {config['positions_path']}")
    if config.get("fasta_path") and not Path(config["fasta_path"]).exists():
        return _error(f"Reference FASTA file not found: {config['fasta_path']}")

    try:
        return run_gff_extraction(config, results_dir)
    except GeneExtractionError as e:
        return _error(str(e))
    except Exception as e:
        return _error(f"Gene extraction failed: {e}")
