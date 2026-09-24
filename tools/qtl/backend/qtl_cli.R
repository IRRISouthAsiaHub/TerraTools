#!/usr/bin/env Rscript

# ============================================================
# TerraTools QTL CLI
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage:\n",
    "Rscript qtl_cli.R <config.json> <results_dir>"
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
# LOAD QTL ENGINE
# ============================================================

engine_file <- file.path(
  getwd(),
  "qtl_engine.R"
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

geno_file <- config$genotype_file
pheno_file <- config$phenotype_file
trait <- config$phenotype_column

method <- config$mapping_method

lod_threshold <- as.numeric(
  config$lod_threshold
)

validation_fraction <- as.numeric(
  config$validation_fraction
)

n_cofactors <- as.numeric(
  config$n_cofactors
)

window_size <- as.numeric(
  config$window_size
)

lod_drop <- as.numeric(
  config$lod_drop
)

# ============================================================
# VALIDATE INPUTS
# ============================================================

if (!file.exists(geno_file)) {
  stop(
    paste(
      "Genotype file not found:",
      geno_file
    )
  )
}

if (!file.exists(pheno_file)) {
  stop(
    paste(
      "Phenotype file not found:",
      pheno_file
    )
  )
}

# ============================================================
# LOAD GENOTYPE
# ============================================================

geno_res <- load_genotype_file(
  geno_file
)

geno <- geno_res$geno

marker_info <- geno_res$marker_info

all_markers <- geno_res$marker_info

# ============================================================
# LOAD PHENOTYPE
# ============================================================

pheno <- read.csv(
  pheno_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rownames(pheno) <- pheno[, 1]

pheno <- pheno[
  ,
  -1,
  drop = FALSE
]

for (col in colnames(pheno)) {

  pheno[[col]] <- suppressWarnings(
    as.numeric(
      as.character(
        pheno[[col]]
      )
    )
  )
}

if (!(trait %in% colnames(pheno))) {

  stop(
    paste(
      "Trait not found:",
      trait,
      "\nAvailable traits:",
      paste(
        colnames(pheno),
        collapse = ", "
      )
    )
  )
}

# ============================================================
# MATCH SAMPLES
# ============================================================

common_samples <- intersect(
  rownames(geno),
  rownames(pheno)
)

if (length(common_samples) < 8) {

  stop(
    paste(
      "Only",
      length(common_samples),
      "matching samples found."
    )
  )
}

X <- geno[
  common_samples,
  ,
  drop = FALSE
]

y <- as.numeric(
  pheno[
    common_samples,
    trait
  ]
)

# ============================================================
# TRAINING SET
# ============================================================

set.seed(42)

if (validation_fraction > 0) {

  train_idx <- caret::createDataPartition(
    y,
    p = 1 - validation_fraction,
    list = FALSE
  )

  train_idx <- as.vector(
    train_idx
  )

} else {

  train_idx <- seq_along(y)
}

X_train <- X[
  train_idx,
  ,
  drop = FALSE
]

y_train <- y[
  train_idx
]

# ============================================================
# QTL SCAN
# ============================================================

cat(
  "\nRunning:",
  method,
  "\n"
)

scan_all <- qtl_scan(
  X_train,
  y_train,
  marker_info,
  method = method,
  n_cofactors = n_cofactors,
  window_size = window_size
)

# ============================================================
# CHECK RESULT
# ============================================================

if ("Message" %in% colnames(scan_all)) {

  stop(
    scan_all$Message[1]
  )
}

# ============================================================
# SIGNIFICANT QTL
# ============================================================

sig <- scan_all[
  !is.na(scan_all$LOD) &
  scan_all$LOD >= lod_threshold,
  ,
  drop = FALSE
]

if (nrow(sig) > 0) {

  sig <- add_lod_support_interval(
    sig,
    scan_all,
    lod_drop
  )

  sig <- sig[
    order(
      -sig$LOD
    ),
    ,
    drop = FALSE
  ]
}

# ============================================================
# OUTPUT PATHS
# ============================================================

manhattan_file <- file.path(
  results_dir,
  paste0(
    "manhattan_",
    method,
    ".png"
  )
)

ideogram_file <- file.path(
  results_dir,
  paste0(
    "ideogram_",
    method,
    ".png"
  )
)

qtl_file <- file.path(
  results_dir,
  paste0(
    "qtl_significant_",
    method,
    ".csv"
  )
)

full_file <- file.path(
  results_dir,
  paste0(
    "qtl_full_scan_",
    method,
    ".csv"
  )
)

# ============================================================
# PLOT SETTINGS
# ============================================================

plot_cfg <- config$plot

# ============================================================
# MANHATTAN PNG
# ============================================================

png(
  manhattan_file,
  width = as.numeric(
    plot_cfg$manhattan_width
  ),
  height = as.numeric(
    plot_cfg$manhattan_height
  ),
  res = 150,
  bg = "white"
)

pd <- scan_all[
  !is.na(scan_all$LOD) &
  !is.na(scan_all$Position) &
  !is.na(scan_all$Chromosome),
  ,
  drop = FALSE
]

if (nrow(pd) > 0) {

  pd$Chromosome <- as.character(
    pd$Chromosome
  )

  pd$sig <- pd$LOD >= lod_threshold

  chroms <- sort(
    unique(
      pd$Chromosome
    )
  )

  n_chroms <- length(
    chroms
  )

  n_cols <- ceiling(
    sqrt(n_chroms)
  )

  n_rows <- ceiling(
    n_chroms / n_cols
  )

  par(
    mfrow = c(
      n_rows,
      n_cols
    ),
    mar = c(
      4,
      4,
      3,
      1
    ),
    bg = "white"
  )

  for (chr in chroms) {

    chr_data <- pd[
      pd$Chromosome == chr,
      ,
      drop = FALSE
    ]

    chr_data <- chr_data[
      order(
        chr_data$Position
      ),
      ,
      drop = FALSE
    ]

    y_max <- max(
      chr_data$LOD,
      na.rm = TRUE
    )

    if (!is.finite(y_max) ||
        y_max <= lod_threshold) {

      y_max <-
        lod_threshold * 1.2

    } else {

      y_max <-
        y_max * 1.1
    }

    plot(
      chr_data$Position,
      chr_data$LOD,
      type = "n",
      xlab = paste(
        "Chr",
        chr,
        "Position (cM)"
      ),
      ylab = "LOD Score",
      main = paste(
        "Chr",
        chr
      ),
      ylim = c(
        0,
        y_max
      ),
      frame.plot = FALSE
    )

    nonsig <- chr_data[
      !chr_data$sig,
      ,
      drop = FALSE
    ]

    sig_data <- chr_data[
      chr_data$sig,
      ,
      drop = FALSE
    ]

    if (nrow(nonsig) > 0) {

      points(
        nonsig$Position,
        nonsig$LOD,
        pch = 16,
        cex =
          as.numeric(
            plot_cfg$nonsig_point_size
          ) / 5,
        col =
          plot_cfg$nonsig_color
      )
    }

    if (nrow(sig_data) > 0) {

      points(
        sig_data$Position,
        sig_data$LOD,
        pch = 18,
        cex =
          as.numeric(
            plot_cfg$sig_point_size
          ) / 5,
        col =
          plot_cfg$sig_color
      )
    }

    lty <- switch(
      plot_cfg$threshold_style,
      solid = 1,
      dash = 2,
      dot = 3,
      dashdot = 4,
      2
    )

    abline(
      h = lod_threshold,
      lty = lty,
      col = plot_cfg$sig_color,
      lwd =
        as.numeric(
          plot_cfg$threshold_width
        )
    )
  }
}

dev.off()

# ============================================================
# IDEOGRAM
# ============================================================

png(
  ideogram_file,
  width =
    as.numeric(
      plot_cfg$ideogram_width
    ),
  height =
    as.numeric(
      plot_cfg$ideogram_height
    ),
  res = 150,
  bg = "white"
)

qtl_list <- list()

qtl_list[[method]] <- scan_all

make_chromosome_ideogram_reference(
  qtl_list,
  all_markers,
  lod_threshold,

  title_str =
    paste(
      method,
      "— Trait:",
      trait
    ),

  band_dark_col =
    plot_cfg$chr_dark,

  band_light_col =
    plot_cfg$chr_light,

  border_col =
    plot_cfg$chr_border,

  qtl_size_min =
    as.numeric(
      plot_cfg$qtl_marker_size_min
    ),

  qtl_size_max =
    as.numeric(
      plot_cfg$qtl_marker_size_max
    )
)

dev.off()

# ============================================================
# CSV
# ============================================================

write.csv(
  sig,
  qtl_file,
  row.names = FALSE
)

write.csv(
  scan_all,
  full_file,
  row.names = FALSE
)

# ============================================================
# RESULT JSON
# ============================================================

result <- list(

  status = "success",

  summary = list(
    samples = nrow(X),
    markers = ncol(X),
    chromosomes =
      length(
        unique(
          marker_info$Chromosome
        )
      ),
    significant_qtl =
      nrow(sig),
    method = method,
    trait = trait
  ),

  manhattan =
    normalizePath(
      manhattan_file,
      winslash = "/"
    ),

  ideogram =
    normalizePath(
      ideogram_file,
      winslash = "/"
    ),

  qtl_csv =
    normalizePath(
      qtl_file,
      winslash = "/"
    ),

  full_csv =
    normalizePath(
      full_file,
      winslash = "/"
    ),

  qtl_table =
    if (nrow(sig) > 0) {

      jsonlite::fromJSON(
        jsonlite::toJSON(
          sig,
          dataframe = "rows",
          na = "null",
          auto_unbox = TRUE
        )
      )

    } else {

      list()
    }
)

cat(
  jsonlite::toJSON(
    result,
    auto_unbox = TRUE,
    na = "null"
  )
)
