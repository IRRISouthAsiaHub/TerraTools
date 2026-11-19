#!/usr/bin/env Rscript

# ============================================================
# TerraTools QTL: the "Common QTL" ideogram.
#   Rscript qtl_common_ideogram.R <config.json> <out.png>
#
# The same chromosome map as the single-model run and the Compare Models figure
# (make_chromosome_ideogram_reference of qtl_engine.R: banded chromosomes, flanking-marker
# bars, peak dots), with one lane per "found by N models" level and one colour per level.
#
# config: list(regions_csv = common_regions.csv, full_csvs = c(<full-scan csv of each model>),
#              plot = list(<plot settings of the page>), min_models = n, n_models = number of models run)
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage:\nRscript qtl_common_ideogram.R <config.json> <out.png>")

suppressPackageStartupMessages(library(jsonlite))
source(normalizePath(file.path(getwd(), "qtl_engine.R"), mustWork = TRUE), local = TRUE)

cfg <- fromJSON(normalizePath(args[1], mustWork = TRUE), simplifyVector = FALSE)
out <- args[2]
pc  <- cfg$plot
num <- function(x, d) { v <- suppressWarnings(as.numeric(x)); if (length(v) != 1 || is.na(v)) d else v }
txt <- function(x, d) if (is.null(x) || !nzchar(as.character(x))) d else as.character(x)

markers <- NULL
for (f in unlist(cfg$full_csvs)) {
  sc <- tryCatch(read.csv(f, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (!is.null(sc) && all(c("Chromosome", "Position") %in% names(sc))) markers <- rbind(markers, sc[, c("Chromosome", "Position")])
}
if (is.null(markers) || nrow(markers) == 0) stop("No marker positions were found in the full-scan files of the models.")

common <- read.csv(cfg$regions_csv, stringsAsFactors = FALSE, check.names = FALSE)
need   <- c("Chromosome", "Left_Position", "Right_Position", "Peak_Position", "N_Models_Detected")
min_m  <- as.integer(num(cfg$min_models, 2))
max_m  <- max(as.integer(num(cfg$n_models, 2)), min_m)
levels <- seq(max(2, min_m), max_m)

# one lane per detection level; the dot is sized and labelled by the number of models
palette <- c("#F2A541", "#E4572E", "#B5179E", "#4361EE", "#2B2D42")
lane_names <- paste0("Found by ", levels, " models")
colours <- setNames(palette[seq_along(levels)], lane_names)
qtl_list <- list()
for (k in seq_along(levels)) {
  q <- if (nrow(common) && all(need %in% names(common))) common[common$N_Models_Detected == levels[k], , drop = FALSE] else common[0, , drop = FALSE]
  if (nrow(q)) q$LOD <- q$N_Models_Detected * 1.5          # only used for the size of the dot
  if (nrow(q) || length(levels) == 1) qtl_list[[lane_names[k]]] <- q          # a level without a region gets no lane
}
if (!length(qtl_list)) qtl_list[[lane_names[1]]] <- common[0, , drop = FALSE]

lane_step <- 0.36
n_lane    <- length(qtl_list)
x_gap     <- max(1.8, 0.55 / 2 + 0.12 + lane_step * (n_lane - 1) + 0.95)
w         <- num(pc$ideogram_width, 1400) * x_gap / 1.8
h         <- num(pc$ideogram_height, 900)

png(out, width = w, height = h, res = 150, bg = "white")
make_chromosome_ideogram_reference(
  qtl_list, markers, 0,
  model_colors   = colours,
  title_str      = paste0("Common QTL regions (found by at least ", min_m, " models)"),
  band_dark_col  = txt(pc$chr_dark, "#2E6B35"),
  band_light_col = txt(pc$chr_light, "#6AAB6E"),
  border_col     = txt(pc$chr_border, "#1A4A1F"),
  qtl_size_min   = min(num(pc$qtl_marker_size_min, 0.7), 0.7),
  qtl_size_max   = min(num(pc$qtl_marker_size_max, 2), 1.0),
  x_gap          = x_gap,
  lane_step      = lane_step,
  label_fun      = function(q) paste0("x", q$N_Models_Detected),
  legend_title   = "Detected by",
  legend_note    = "bar = flanking-marker interval\ndot = QTL peak\nxN = found by N models"
)
invisible(dev.off())
cat(normalizePath(out, winslash = "/"), "\n")
