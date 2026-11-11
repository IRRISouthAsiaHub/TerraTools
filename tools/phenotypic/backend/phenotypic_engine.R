# phenotypic_engine.R — loads the analysis engine (the files in engine/, in order).
# Usage: source("phenotypic_engine.R")   (phenotypic_cli.R does this)
#
# The engine is the Advanced Phenotypic Analysis app: run_phenotypic() runs the analyses a caller asks for
# (descriptives, environments, BLUEs, ANOVA, heritability, correlation, PCA, clustering, stability) and writes tables
# (CSV), figures (HTML, 300-dpi PNG, vector PDF) and an Excel workbook. See engine/01_setup.R for the formulas.
.engine_dir <- local({
  d <- if (exists(".tt_backend_dir", inherits = TRUE)) get(".tt_backend_dir", inherits = TRUE) else "."
  for (i in rev(seq_len(sys.nframe()))) {          # the folder of the file being sourced (innermost source() call)
    of <- tryCatch(get("ofile", envir = sys.frame(i), inherits = FALSE), error = function(e) NULL)
    if (is.character(of)) { d <- dirname(of); break }
  }
  normalizePath(d)
})
for (.f in sort(list.files(file.path(.engine_dir, "engine"), pattern = "\\.R$", full.names = TRUE))) source(.f, encoding = "UTF-8")
rm(.f)
