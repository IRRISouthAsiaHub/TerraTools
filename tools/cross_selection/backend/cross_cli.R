#!/usr/bin/env Rscript

# ============================================================
# TerraTools Cross Selection CLI
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage:\n",
    "Rscript cross_cli.R <config.json> <results_dir>"
  )
}

config_file <- normalizePath(
  args[1],
  mustWork = TRUE
)

results_dir <- normalizePath(
  args[2],
  mustWork = FALSE
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# PACKAGES
# ============================================================

suppressPackageStartupMessages({
  library(jsonlite)
})

# ============================================================
# LOAD CROSS SELECTION ENGINE
# ============================================================

engine_file <- file.path(
  getwd(),
  "cross_engine.R"
)

engine_file <- normalizePath(
  engine_file,
  mustWork = TRUE
)

source(
  engine_file,
  local = TRUE
)

# ============================================================
# CONFIG
# ============================================================

config <- fromJSON(
  config_file,
  simplifyVector = FALSE
)

criterion_file <- config$criterion_file
kmat_file <- config$kmat_file

threshold <- as.numeric(
  config$threshold
)

max_per_cluster <- as.numeric(
  config$max_per_cluster
)

plot_cfg <- config$plot

# ============================================================
# VALIDATE INPUTS
# ============================================================

if (!file.exists(criterion_file)) {
  stop(
    paste(
      "Criterion file not found:",
      criterion_file
    )
  )
}

if (!file.exists(kmat_file)) {
  stop(
    paste(
      "Kinship matrix file not found:",
      kmat_file
    )
  )
}

# ============================================================
# LOAD CRITERION
# ============================================================

crit <- load_data(
  criterion_file,
  c("Genotype", "Criterion")
)

crit$Genotype <- cn(crit$Genotype)

n_criterion_rows <- nrow(crit)

# ============================================================
# LOAD KINSHIP MATRIX
# ============================================================

kmat_raw <- load_data(
  kmat_file
)

Kmat <- process_kmat(
  kmat_raw
)

n_kmat_genotypes <- nrow(Kmat)

# ============================================================
# OVERLAP CHECK
# ============================================================

crit_filtered <- crit[
  crit$Genotype %in% rownames(Kmat),
]

n_overlap <- nrow(crit_filtered)

if (n_overlap == 0) {

  stop(
    "No matching genotypes between the criterion file and the kinship matrix. ",
    "Check that Genotype names match exactly (case, whitespace, IDs)."
  )
}

# ============================================================
# RELATEDNESS THINNING
# ============================================================

keep <- relate_thinning(
  Kmat,
  crit_filtered,
  threshold = threshold,
  max_per_cluster = max_per_cluster
)

all_parents <- intersect(
  keep,
  rownames(Kmat)
)

if (length(all_parents) < 2) {

  stop(
    "Fewer than 2 parents remain after relatedness thinning. ",
    "Try a higher kinship threshold or a larger max-per-cluster value."
  )
}

# ============================================================
# HALF-DIALLEL CROSS GENERATION
# ============================================================

crosses <- generate_half_diallel(
  all_parents,
  crit_filtered,
  Kmat
)

parent_usage <- build_parent_usage(
  crosses,
  crit_filtered
)

# ============================================================
# OUTPUT PATHS
# ============================================================

crosses_file <- file.path(
  results_dir,
  "crosses_full.csv"
)

parent_file <- file.path(
  results_dir,
  "parent_usage.csv"
)

hist_pov_file <- file.path(
  results_dir,
  "hist_pov.png"
)

hist_kinship_file <- file.path(
  results_dir,
  "hist_kinship.png"
)

bar_relatedness_file <- file.path(
  results_dir,
  "bar_relatedness.png"
)

scatter_vk_file <- file.path(
  results_dir,
  "scatter_vk.png"
)

parent_bar_file <- file.path(
  results_dir,
  "parent_bar.png"
)

# ============================================================
# PLOTS
# ============================================================

hist_width  <- as.numeric(get_color(plot_cfg$hist_width, 900))
hist_height <- as.numeric(get_color(plot_cfg$hist_height, 550))

parent_bar_width  <- as.numeric(get_color(plot_cfg$parent_bar_width, 900))
parent_bar_height <- as.numeric(get_color(plot_cfg$parent_bar_height, 900))

bins <- as.numeric(get_color(plot_cfg$bins, 40))

pov_color        <- get_color(plot_cfg$pov_color, "#3A7D44")
kinship_color    <- get_color(plot_cfg$kinship_color, "#C8760A")
mean_line_color  <- get_color(plot_cfg$mean_line_color, "#1A4A1F")
mod_line_color   <- get_color(plot_cfg$mod_line_color, "#E8834A")
close_line_color <- get_color(plot_cfg$close_line_color, "#990000")
parent_bar_color <- get_color(plot_cfg$parent_bar_color, "#3A7D44")

ggsave(
  hist_pov_file,
  make_pov_histogram(crosses, bar_color = pov_color, mean_line_color = mean_line_color, bins = bins),
  width = hist_width / 100, height = hist_height / 100, dpi = 100, bg = "white"
)

ggsave(
  hist_kinship_file,
  make_kinship_histogram(crosses, bar_color = kinship_color, mod_line_color = mod_line_color,
                          close_line_color = close_line_color, bins = bins),
  width = hist_width / 100, height = hist_height / 100, dpi = 100, bg = "white"
)

ggsave(
  bar_relatedness_file,
  make_relatedness_bar(crosses),
  width = hist_width / 100, height = hist_height / 100, dpi = 100, bg = "white"
)

ggsave(
  scatter_vk_file,
  make_value_kinship_scatter(crosses),
  width = hist_width / 100, height = hist_height / 100, dpi = 100, bg = "white"
)

ggsave(
  parent_bar_file,
  make_parent_bar(parent_usage, bar_color = parent_bar_color),
  width = parent_bar_width / 100,
  height = max(3, length(all_parents) * 0.22) ,
  dpi = 100, bg = "white", limitsize = FALSE
)

# ============================================================
# CSV
# ============================================================

write.csv(
  crosses,
  crosses_file,
  row.names = FALSE
)

write.csv(
  parent_usage,
  parent_file,
  row.names = FALSE
)

# ============================================================
# RESULT JSON
# ============================================================

result <- list(

  status = "success",

  summary = list(
    n_criterion_rows = n_criterion_rows,
    n_kmat_genotypes = n_kmat_genotypes,
    n_overlap = n_overlap,
    n_parents = length(all_parents),
    n_crosses = nrow(crosses),
    mean_pov = round(mean(crosses$Predicted_Offspring_Value), 4),
    max_pov = round(max(crosses$Predicted_Offspring_Value), 4),
    mean_kinship = round(mean(crosses$Kinship_Coefficient), 6),
    n_low_kinship = sum(crosses$Kinship_Coefficient < 0.125)
  ),

  plots = list(
    hist_pov = normalizePath(hist_pov_file, winslash = "/"),
    hist_kinship = normalizePath(hist_kinship_file, winslash = "/"),
    bar_relatedness = normalizePath(bar_relatedness_file, winslash = "/"),
    scatter_vk = normalizePath(scatter_vk_file, winslash = "/"),
    parent_bar = normalizePath(parent_bar_file, winslash = "/")
  ),

  crosses_csv = normalizePath(crosses_file, winslash = "/"),

  parent_csv = normalizePath(parent_file, winslash = "/"),

  top_crosses = jsonlite::fromJSON(
    jsonlite::toJSON(
      head(crosses, 100),
      dataframe = "rows",
      na = "null",
      auto_unbox = TRUE
    )
  ),

  parent_table = jsonlite::fromJSON(
    jsonlite::toJSON(
      parent_usage,
      dataframe = "rows",
      na = "null",
      auto_unbox = TRUE
    )
  )
)

cat(
  jsonlite::toJSON(
    result,
    auto_unbox = TRUE,
    na = "null"
  )
)
