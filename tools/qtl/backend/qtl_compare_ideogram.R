#!/usr/bin/env Rscript

# ============================================================
# TerraTools QTL: the "Compare Models" ideogram.
#   Rscript qtl_compare_ideogram.R <config.json> <out.png>
#
# Draws the SAME chromosome map as the single-model run (make_chromosome_ideogram_reference
# of qtl_engine.R: same chromosome bands, same flanking-marker bars, peak dots and LOD labels),
# with one lane per model and one colour per model.
#
# config: list(models = list(list(name, qtl_csv, full_csv, color), ...), plot = list(<plot settings of the page>))
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage:\nRscript qtl_compare_ideogram.R <config.json> <out.png>")

suppressPackageStartupMessages(library(jsonlite))
source(normalizePath(file.path(getwd(), "qtl_engine.R"), mustWork = TRUE), local = TRUE)

cfg  <- fromJSON(normalizePath(args[1], mustWork = TRUE), simplifyVector = FALSE)
out  <- args[2]
pc   <- cfg$plot
num  <- function(x, d) { v <- suppressWarnings(as.numeric(x)); if (length(v) != 1 || is.na(v)) d else v }
txt  <- function(x, d) if (is.null(x) || !nzchar(as.character(x))) d else as.character(x)

qtl_list <- list(); colours <- character(0); markers <- NULL
for (m in cfg$models) {
  nm <- toupper(as.character(m$name))
  colours[nm] <- m$color
  q <- tryCatch(read.csv(m$qtl_csv, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  need <- c("Chromosome", "Left_Position", "Right_Position", "Peak_Position", "LOD")
  qtl_list[[nm]] <- if (!is.null(q) && nrow(q) > 0 && all(need %in% names(q))) q else q[0, , drop = FALSE]
  sc <- tryCatch(read.csv(m$full_csv, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (!is.null(sc) && all(c("Chromosome", "Position") %in% names(sc))) {
    markers <- rbind(markers, sc[, c("Chromosome", "Position")])
  }
}
if (is.null(markers) || nrow(markers) == 0) stop("No marker positions were found in the full-scan files of the models.")

# one lane per model: widen the spacing between chromosomes so that the lanes of one chromosome never touch the next
lane_step <- 0.36
n_mod     <- length(qtl_list)
x_gap     <- max(1.8, 0.55 / 2 + 0.12 + lane_step * (n_mod - 1) + 0.75)
w         <- num(pc$ideogram_width, 1400) * x_gap / 1.8
h         <- num(pc$ideogram_height, 900)

png(out, width = w, height = h, res = 150, bg = "white")
make_chromosome_ideogram_reference(
  qtl_list, markers, 0,
  model_colors   = colours,
  title_str      = paste("Model comparison -", cfg$trait),
  band_dark_col  = txt(pc$chr_dark, "#2E6B35"),
  band_light_col = txt(pc$chr_light, "#6AAB6E"),
  border_col     = txt(pc$chr_border, "#1A4A1F"),
  qtl_size_min   = min(num(pc$qtl_marker_size_min, 0.7), 0.7),     # up to six lanes side by side: smaller dots so that they do not cover each other
  qtl_size_max   = min(num(pc$qtl_marker_size_max, 2), 1.0),
  x_gap          = x_gap,
  lane_step      = lane_step
)
invisible(dev.off())
cat(normalizePath(out, winslash = "/"), "\n")
