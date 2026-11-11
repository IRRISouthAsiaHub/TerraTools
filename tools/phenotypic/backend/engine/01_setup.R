# =============================================================================
# Advanced Phenotypic Analysis: statistical engine
# =============================================================================
# Loaded by phenotypic_engine.R (it sources the files of this folder in order). run_phenotypic() in 07_pipeline.R runs the
# analyses a caller asks for and writes tables (CSV), figures (HTML, PNG, PDF) and an Excel workbook.
#
#   01_setup.R        packages, constants, colours, plot themes, figure writers
#   02_data_design.R  reading and validating data; design models (formulas, augmented RCBD)
#   03_analysis.R     the statistical analyses (each returns NULL when the data cannot support it)
#   04_plots.R        ggplot builders
#   05_tables_excel.R readable tables and the Excel workbook
#   06_modules.R      one function per analysis tab
#   07_pipeline.R     run_phenotypic()
#
# ANALYSES AND THEIR FORMULAS
# ---------------------------------------------------------------------------
#  Design models (design_formula)
#    Per-env ANOVA   Augmented: Federer fixed-block analysis (Value ~ Block + Genotype),
#                    otherwise  RCBD: Value ~ Rep + Genotype
#                               Alpha: Value ~ Rep + RepBlock + Genotype
#    Per-env BLUE    Value ~ Genotype + (1|Rep) [+ (1|RepBlock)]         (emmeans)
#    Across-env BLUE Value ~ Environment + Genotype [+ Environment:Genotype] + (1|EnvRep ...)
#    Pooled ANOVA    Value ~ Environment + blocks nested in Environment + Genotype [+ Genotype:Environment]
#  ANOVA             every trait x season (Federer table for augmented RCBD, aov otherwise) and pooled over seasons;
#                    CV% = sqrt(MS_error) / mean * 100
#  Descriptives      CV% = SD / mean * 100;  Kurtosis = EXCESS kurtosis (normal = 0)
#  Heritability      Checks are fixed effects, Genotype and Genotype:Environment random
#                    H2_plot  = Vg / (Vg + Vge + Ve)
#                    H2_entry = Vg / (Vg + Vge/ne + Ve/np)    (entry-mean basis; ne, np = harmonic means per entry)
#  Correlation, PCA, clustering   on genotype means (mean of environment means)
#  Stability (genotype x environment means x_ij, only genotypes present in every environment)
#                    GE_ij       = x_ij - x_i. - x_.j + x_..          (interaction effect)
#                    Ecovalence  = Wi = sum_j GE_ij^2                  (Wricke 1962)
#                    Slope_b     = sum_j (x_ij - x_i.) * I_j / sum_j I_j^2,  I_j = x_.j - x_..  (Finlay-Wilkinson)
#                    AMMI        = SVD of the GE matrix; scores = U*sqrt(d), V*sqrt(d)
#                    ASV         = sqrt( (d1^2/d2^2 * IPCA1)^2 + IPCA2^2 )
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(ggplot2)   # data handling and plots
  library(lme4); library(emmeans)                                     # mixed models and adjusted means
  library(plotly); library(htmlwidgets)                               # interactive figures
  library(openxlsx)                                                   # Excel workbooks
})
# Used as pkg::fun(): readxl, readr, stringr, corrplot, broom, moments, ragg, pbkrtest, viridisLite

`%||%` <- function(a, b) if (is.null(a)) b else a

# =============================================================================
# 1. CONSTANTS AND SMALL HELPERS
# =============================================================================

# Columns are found by NAME (case and surrounding spaces do not matter). Design columns:
# Genotype, Environment, Rep, Block. Every other numeric column is a trait, except layout
# columns such as Plot, Row, Column or Year. Keep LAYOUT_NAMES in step with tools/phenotypic/columns.py.
DESIGN_COLS  <- c("Genotype", "Environment", "Rep", "Block")
LAYOUT_NAMES <- c("plot", "plotno", "plotid", "row", "column", "col", "range", "year", "season", "serial",
                  "sno", "entry", "entryno", "position", "type", "linetype", "check", "remarks", "notes")

# augmented : checks repeated in every block, new entries once      (Rep + Block)
# rcbd      : complete blocks, every genotype once per replicate    (Rep)
# alpha     : alpha-lattice, incomplete blocks nested in replicates (Rep + Block within Rep)
DESIGN_LABELS <- c(augmented = "Augmented RCBD", rcbd = "RCBD", alpha = "Alpha-lattice")

canonical_names <- function(nm) {
  nm <- trimws(nm)
  for (cn in DESIGN_COLS) nm[tolower(nm) == tolower(cn)] <- cn
  nm
}

is_layout_name <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x)) %in% LAYOUT_NAMES

# Numeric columns that can be traits; `chosen` (if given) restricts them to the columns the user ticked.
find_traits <- function(data, chosen = character(0)) {
  num <- names(data)[vapply(data, is.numeric, logical(1))]
  num <- setdiff(num, DESIGN_COLS)
  if (length(chosen) > 0) intersect(chosen, num) else num[!is_layout_name(num)]
}

# Traits that have at least one real value (0 counts as missing): the only ones the analyses and selectors offer
usable_traits <- function(data) {
  tr <- find_traits(data)
  tr[vapply(tr, function(x) any(!is.na(data[[x]]) & data[[x]] != 0), logical(1))]
}

normalize_design <- function(x) {
  x <- tolower(gsub("[^a-z]", "", if (is.null(x) || length(x) == 0) "" else x[[1]]))
  if (x %in% c("rcbd", "randomizedcompleteblock")) "rcbd"
  else if (x %in% c("alpha", "alphalattice", "lattice")) "alpha"
  else "augmented"
}

design_required_cols <- function(design) {
  switch(design, rcbd = c("Genotype", "Environment", "Rep"), c("Genotype", "Environment", "Rep", "Block"))
}

harmonic_mean <- function(x) length(x) / sum(1 / x)

# Sets 0 to NA in every trait column; returns the data and one note per trait listing the seasons affected
zeros_to_missing <- function(d, traits) {
  notes <- character(0)
  for (tr in traits) {
    z <- !is.na(d[[tr]]) & d[[tr]] == 0
    if (!any(z)) next
    by_env <- tapply(z, as.character(d$Environment), sum)
    by_env <- by_env[by_env > 0]
    notes <- c(notes, sprintf("%s: %d zero value(s) treated as missing (%s).", tr, sum(z),
                              paste0(names(by_env), ": ", by_env, collapse = ", ")))
    d[[tr]][z] <- NA
  }
  list(data = d, notes = notes)
}

# ── Colours and plot styling ─────────────────────────────────────────────────
get_color_palette <- function(scheme, n) {
  n <- max(n, 1)
  if (scheme == "blue") {
    return(colorRampPalette(c("#E6F0FA", "#1E3A8A"))(n))
  } else if (scheme %in% c("viridis", "plasma", "magma") && requireNamespace("viridisLite", quietly = TRUE)) {
    fn <- switch(scheme, viridis = viridisLite::viridis, plasma = viridisLite::plasma, magma = viridisLite::magma)
    return(fn(n))
  }
  colorRampPalette(c("#1E3A8A", "#6B8E23", "#FF8C00", "#DC143C", "#4169E1"))(n)
}

# x-axis label rotation for environment-labelled axes
label_geom <- function(n_env) {
  angle <- if (n_env > 6) 90 else if (n_env > 3) 45 else 0
  list(angle = angle,
       hjust = if (angle > 0) 1 else 0.5,
       vjust = if (angle == 90) 0.5 else 1,
       bmargin = if (angle >= 45) 60 else 20)
}

# One theme for every ggplot (replaces the copy-pasted theme() blocks)
theme_pheno <- function(lg = NULL, legend = "right") {
  th <- theme_minimal() +
    theme(axis.text  = element_text(size = 11, color = "#000"),
          axis.title = element_text(size = 12, color = "#000"),
          plot.title = element_text(size = 13, color = "#000", face = "bold"),
          legend.position = legend)
  if (!is.null(lg))
    th <- th + theme(axis.text.x = element_text(angle = lg$angle, hjust = lg$hjust, vjust = lg$vjust),
                     plot.margin = margin(t = 10, r = 10, b = lg$bmargin, l = 10))
  th
}

text_size_for <- function(n) if (n > 40) 6 else if (n > 20) 8 else 10

# selfcontained = FALSE avoids a pandoc dependency; the plotly JavaScript (~4 MB) is written once, to a "lib" folder
# next to the figures, and shared by every interactive figure instead of being copied for each one
save_plotly <- function(p, path) htmlwidgets::saveWidget(plotly::as_widget(p), path, selfcontained = FALSE, libdir = "lib")

png_device <- function(file, width, height) {
  if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png(file, width = width, height = height, units = "in", res = 300, background = "white")
  else png(file, width = width, height = height, units = "in", res = 300, bg = "white")
}

# 300-dpi PNG and vector PDF of a ggplot, at manuscript size (inches). ggsave() measures the text on the device it draws
# on, so each file gets its own build: drawing one pre-built plot on two devices shifts the axis labels.
save_figure <- function(p_gg, png_fp, pdf_fp, size_in) {
  tryCatch(ggsave(png_fp, plot = p_gg, width = size_in[1], height = size_in[2], dpi = 300, units = "in", bg = "white",
                  device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"), error = function(e) NULL)
  tryCatch(ggsave(pdf_fp, plot = p_gg, width = size_in[1], height = size_in[2], units = "in", device = "pdf"), error = function(e) NULL)
}

# Figure sizes in inches for manuscripts: a 7 in wide figure fits a journal page; every figure is written at 300 dpi
# (PNG) and as vector PDF.
FIG_STD <- c(7, 5); FIG_WIDE <- c(8, 5.5); FIG_SQUARE <- c(6.5, 6); FIG_SMALL <- c(6.5, 4.5)
fig_tall <- function(n, width = 7) c(width, min(12, max(5, 1.8 + 0.17 * n)))

