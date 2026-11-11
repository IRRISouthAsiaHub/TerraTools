# Runs the mrMLM R package (mrMLM and/or FASTmrMLM) on files written by mrmlm_package.py.
#
# Rscript mrmlm_run.R <genotype.csv> <phenotype.csv> <popstr.csv|none> <outdir> <methods, comma separated> <search radius> <critical LOD>
#
# The package does all the statistics; this script only passes its arguments and reports whether it finished.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) stop("expected 7 arguments")
gen <- args[1]; phe <- args[2]; ps <- args[3]; out <- args[4]
methods <- strsplit(args[5], ",")[[1]]
radius <- as.numeric(args[6]); lod <- as.numeric(args[7])

if (!requireNamespace("mrMLM", quietly = TRUE)) {
  writeLines("mrMLM package is not installed", file.path(out, "error.txt"))
  quit(status = 3)
}
suppressPackageStartupMessages(library(mrMLM))
dir.create(out, showWarnings = FALSE, recursive = TRUE)

res <- try(
  if (identical(ps, "none")) {
    mrMLM(fileGen = gen, filePhe = phe, Genformat = "Num", method = methods, trait = 1, SearchRadius = radius,
          CriLOD = lod, DrawPlot = FALSE, dir = out, PC = FALSE, RAM = 4)
  } else {
    mrMLM(fileGen = gen, filePhe = phe, filePS = ps, PopStrType = "PCA", Genformat = "Num", method = methods, trait = 1,
          SearchRadius = radius, CriLOD = lod, DrawPlot = FALSE, dir = out, PC = FALSE, RAM = 4)
  },
  silent = TRUE)
if (inherits(res, "try-error")) {
  writeLines(as.character(res), file.path(out, "error.txt"))
  quit(status = 2)
}
writeLines(as.character(packageVersion("mrMLM")), file.path(out, "version.txt"))
