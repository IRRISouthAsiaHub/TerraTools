#!/usr/bin/env Rscript

# ============================================================
# TerraTools QTL CLI
#   Rscript qtl_cli.R <config.json> <results_dir>
#
# Reports QTL as regions (left flanking marker, peak, right flanking
# marker, LOD-drop support interval).  See qtl_engine.R for the methods.
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage:\nRscript qtl_cli.R <config.json> <results_dir>")
}

config_file <- normalizePath(args[1], mustWork = TRUE)
results_dir <- normalizePath(args[2], mustWork = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(jsonlite)
})

engine_file <- normalizePath(file.path(getwd(), "qtl_engine.R"), mustWork = TRUE)
source(engine_file, local = TRUE)

# ============================================================
# CONFIG
# ============================================================

config <- fromJSON(config_file, simplifyVector = FALSE)

geno_file  <- config$genotype_file
pheno_file <- config$phenotype_file
trait      <- config$phenotype_column
method     <- config$mapping_method

get_config_numeric <- function(config, key, default) {
  value <- config[[key]]
  if (is.null(value) || length(value) == 0 || (length(value) == 1 && is.na(value))) {
    return(default)
  }
  value <- as.numeric(unlist(value))
  if (length(value) != 1 || is.na(value) || !is.finite(value)) {
    stop(paste("Invalid configuration parameter:", key))
  }
  value
}

lod_threshold       <- get_config_numeric(config, "lod_threshold", 3)
validation_fraction <- get_config_numeric(config, "validation_fraction", 0)
n_cofactors         <- get_config_numeric(config, "n_cofactors", 5)
window_size         <- get_config_numeric(config, "window_size", 10)     # cM
lod_drop            <- get_config_numeric(config, "lod_drop", 1.5)
peak_drop           <- get_config_numeric(config, "peak_drop", 2)
p_threshold         <- get_config_numeric(config, "p_threshold", 0.05)
n_perm              <- get_config_numeric(config, "n_perm", 1000)
perm_alpha          <- get_config_numeric(config, "perm_alpha", 0.05)
n_perm_ml           <- get_config_numeric(config, "n_perm_ml", 100)          # permutations of the ML methods (RF / SVR / XGBoost)
ml_budget_s         <- get_config_numeric(config, "ml_time_budget_s", 1200)  # the ML null is cut short rather than run longer than this
bayes_prob          <- get_config_numeric(config, "bayes_prob", 0.95)
impute_missing      <- get_config_numeric(config, "impute_missing", 1) != 0
threshold_mode      <- tolower(if (is.null(config$threshold_mode)) "lod" else as.character(config$threshold_mode))

if (validation_fraction < 0 || validation_fraction >= 1) stop("validation_fraction must be >= 0 and < 1")
if (lod_threshold < 0) stop("lod_threshold must be non-negative")
if (n_cofactors < 0 || n_cofactors != floor(n_cofactors)) stop("n_cofactors must be a non-negative integer")
if (window_size <= 0) stop("window_size (cM) must be greater than 0")
if (lod_drop < 0) stop("lod_drop must be non-negative")
if (peak_drop < 0) stop("peak_drop must be non-negative")
if (p_threshold <= 0 || p_threshold > 1) stop("p_threshold must be in (0, 1]")
if (!(threshold_mode %in% c("lod", "permutation"))) stop("threshold_mode must be 'lod' or 'permutation'")
ml_method <- method %in% ML_METHODS
if (ml_method) {
  # the importance of a marker has no scale of its own: the ML methods are always tested by permutation
  if (n_perm_ml < ML_MIN_PERM || n_perm_ml != floor(n_perm_ml)) stop(paste("n_perm_ml must be a whole number >=", ML_MIN_PERM))
  if (perm_alpha <= 0 || perm_alpha >= 1) stop("perm_alpha must be between 0 and 1")
}
if (threshold_mode == "permutation") {
  if (n_perm < 100 || n_perm != floor(n_perm)) stop("n_perm must be a whole number >= 100")
  if (perm_alpha <= 0 || perm_alpha >= 1) stop("perm_alpha must be between 0 and 1")
}

if (ml_method) threshold_mode <- "permutation"          # (after the checks of the linear permutation settings)
if (is.null(geno_file) || !file.exists(geno_file)) stop(paste("Genotype file not found:", geno_file))
if (is.null(pheno_file) || !file.exists(pheno_file)) stop(paste("Phenotype file not found:", pheno_file))
if (is.null(trait) || !nzchar(trait)) stop("Phenotype column is not specified in config")
if (is.null(method) || !nzchar(method)) stop("Mapping method is not specified in config")

cat("Genotype file:", geno_file, "\n")
cat("Phenotype file:", pheno_file, "\n")
cat("Trait:", trait, "\n")
cat("Mapping method:", method, "\n")
cat("Threshold mode:", threshold_mode, "\n")
if (!ml_method) cat("LOD threshold:", lod_threshold, "\n")
cat("Peak p-value cut-off:", p_threshold, "\n")
cat("Support-interval LOD drop:", lod_drop, " Peak drop:", peak_drop, "\n")
cat("Cofactors:", n_cofactors, " CIM window (cM):", window_size, "\n")
cat("Validation fraction:", validation_fraction, "\n")

# ============================================================
# LOAD GENOTYPE / PHENOTYPE
# ============================================================

geno_res    <- load_genotype_file(geno_file)
geno        <- geno_res$geno
marker_info <- geno_res$marker_info
all_markers <- geno_res$marker_info
cat(sprintf("Missing genotypes: %.1f%%\n", 100 * geno_res$missing_fraction))

pheno <- read.csv(pheno_file, stringsAsFactors = FALSE, check.names = FALSE)
rownames(pheno) <- pheno[, 1]
pheno <- pheno[, -1, drop = FALSE]
for (col in colnames(pheno)) {
  pheno[[col]] <- suppressWarnings(as.numeric(as.character(pheno[[col]])))
}

if (!(trait %in% colnames(pheno))) {
  stop(paste("Trait not found:", trait, "\nAvailable traits:", paste(colnames(pheno), collapse = ", ")))
}

common_samples <- intersect(rownames(geno), rownames(pheno))
if (length(common_samples) < 8) {
  stop(paste("Only", length(common_samples), "matching samples found."))
}

X <- geno[common_samples, , drop = FALSE]
y <- as.numeric(pheno[common_samples, trait])

# ============================================================
# TRAINING SET (validation_fraction = 0 uses every individual)
# ============================================================

set.seed(42)

if (validation_fraction > 0) {
  train_idx <- as.vector(caret::createDataPartition(y, p = 1 - validation_fraction, list = FALSE))
} else {
  train_idx <- seq_along(y)
}

X_train <- X[train_idx, , drop = FALSE]
y_train <- y[train_idx]

# fill missing genotypes from the flanking markers (keeps the LOD curve smooth)
n_imputed <- 0L
if (impute_missing && anyNA(X_train)) {
  imp <- impute_flanking(X_train, marker_info)
  X_train <- imp$X
  n_imputed <- imp$n_imputed
  cat("Missing genotypes imputed from flanking markers:", n_imputed, "\n")
}

# ============================================================
# QTL SCAN
# ============================================================

cat("\nRunning:", method, "\n")

scan_all <- qtl_scan(X_train, y_train, marker_info, method = method,
                     n_cofactors = n_cofactors, window_size = window_size,
                     ml_perm = n_perm_ml, ml_alpha = perm_alpha, ml_budget_s = ml_budget_s)
cofactors_used <- attr(scan_all, "cofactors")

if ("Message" %in% colnames(scan_all)) {
  stop(scan_all$Message[1])
}

# ============================================================
# THRESHOLD
#   lod         -> the LOD cut-off typed by the user
#   permutation -> genome-wide cut-off from shuffling the phenotype
#                  (Churchill & Doerge 1994).  SIM: single-marker scan.  CIM / MQM:
#                  the cofactor selection and the composite scan are repeated on every
#                  permutation (Doerge & Churchill 1996).
# ============================================================

linear_method <- method %in% LINEAR_METHODS
perm_null <- NULL
threshold_note <- ""
thr <- lod_threshold
ml_info <- NULL

if (ml_method) {
  # ML: score = importance / genome-wide permutation threshold of the importance, so the threshold is 1
  thr_imp   <- attr(scan_all, "ml_threshold")
  perm_null <- attr(scan_all, "ml_null") / thr_imp
  thr <- 1
  ml_info <- list(importance_threshold = thr_imp, inflation = attr(scan_all, "ml_inflation"), permutations = attr(scan_all, "ml_nperm"),
                  permutations_requested = attr(scan_all, "ml_nperm_requested"))
  threshold_note <- sprintf("Genome-wide permutation threshold of the %s importance (%d permutations, alpha = %g); score = importance / threshold, a QTL needs a score of at least 1",
                            method, attr(scan_all, "ml_nperm"), perm_alpha)
  if (attr(scan_all, "ml_nperm") < attr(scan_all, "ml_nperm_requested"))
    threshold_note <- paste0(threshold_note, "; permutations reduced from ", attr(scan_all, "ml_nperm_requested"), " to stay within the time budget")
  cat(sprintf("%s importance threshold (alpha=%g): %.4g  (score 1)\n", method, perm_alpha, thr_imp))
} else if (threshold_mode == "permutation") {
  if (linear_method) {
    cat("Permutation test:", n_perm, "permutations ...\n")
    if (method == "SIM") {
      perm_null <- perm_max_lod(X_train, y_train, n_perm = n_perm)
    } else {
      cat("Each permutation repeats the cofactor selection and the", method, "scan; this takes longer than for SIM.
")
      perm_null <- perm_max_lod_cim(X_train, y_train, marker_info, method = method, n_perm = n_perm,
                                    window_cM = window_size, n_cofactors = n_cofactors)
    }
    thr <- as.numeric(quantile(perm_null, 1 - perm_alpha, type = 7))
    threshold_note <- sprintf("Genome-wide permutation threshold (%d permutations, alpha = %g)", n_perm, perm_alpha)
    if (method != "SIM") threshold_note <- paste0(threshold_note, "; cofactors re-selected and the ", method, " scan repeated in every permutation")
    cat(sprintf("Permutation LOD threshold (alpha=%g): %.3f\n", perm_alpha, thr))
  } else {
    threshold_note <- "Permutation threshold is only available for SIM, CIM and MQM; the typed LOD threshold was used."
    cat(threshold_note, "\n")
  }
} else {
  threshold_note <- "Manual LOD threshold"
}

# ============================================================
# QTL REGIONS
# ============================================================

if (ml_method) {
  # score units: a peak must stand out, two peaks are separate if the score falls by at least half a threshold between them
  sig <- call_qtl_regions(scan_all, thr, lod_drop = 0, peak_drop = 0.5, p_max = 1, bayes = FALSE,
                          perm_null = perm_null, interval_floor = 1, min_sep = window_size)
} else {
  sig <- call_qtl_regions(scan_all, thr, lod_drop = lod_drop, peak_drop = peak_drop,
                          p_max = p_threshold, bayes_prob = bayes_prob,
                          bayes = linear_method, perm_null = perm_null)
}
scan_all <- tag_scan_regions(scan_all, sig)

cat("QTL regions found:", nrow(sig), "\n")
if (nrow(sig) > 0) {
  print(sig[, c("QTL_ID", "Left_Marker", "Peak_Marker", "Right_Marker", "Interval_Size_cM", "LOD", "P_Value")], row.names = FALSE)
}

# ============================================================
# HELD-OUT VALIDATION
#   ML methods: the fitted model predicts the lines that were held out.
#   Linear methods: a regression on the peak markers of the QTL found (fitted on the training lines).
#   R2 = 1 - SSE / SS about the training mean (it can be negative: worse than predicting the mean).
# ============================================================

validation <- NULL
oob_r2 <- NULL
if (ml_method && !is.na(attr(scan_all, "ml_model")$oob_r2)) oob_r2 <- attr(scan_all, "ml_model")$oob_r2
test_idx <- setdiff(seq_along(y), train_idx)
if (validation_fraction > 0 && length(test_idx) >= 5) {
  X_te <- X[test_idx, , drop = FALSE]; y_te <- y[test_idx]
  okt <- !is.na(y_te); X_te <- X_te[okt, , drop = FALSE]; y_te <- y_te[okt]
  if (impute_missing && anyNA(X_te)) X_te <- impute_flanking(X_te, marker_info)$X
  pred <- NULL; what <- ""
  if (ml_method) {
    cols <- attr(scan_all, "ml_cols")
    pred <- attr(scan_all, "ml_model")$predict(impute_mode(X_te[, cols, drop = FALSE]))
    what <- paste(method, "model")
  } else if (nrow(sig) > 0) {
    pk <- unique(sig$Peak_Marker)
    Xp <- impute_mean(X_train[, pk, drop = FALSE]); fitp <- lm.fit(cbind(1, Xp), y_train)
    cf <- fitp$coefficients; cf[is.na(cf)] <- 0
    pred <- as.numeric(cbind(1, impute_mean(X_te[, pk, drop = FALSE])) %*% cf)
    what <- paste("regression on the", length(pk), "QTL peak marker(s)")
  }
  if (!is.null(pred) && length(pred) == length(y_te)) {
    sse <- sum((y_te - pred)^2); sst <- sum((y_te - mean(y_train))^2)
    validation <- list(model = what, n_train = length(y_train), n_test = length(y_te),
                       r2 = if (sst > 0) 1 - sse / sst else NA_real_, rmse = sqrt(sse / length(y_te)),
                       correlation = suppressWarnings(cor(y_te, pred)))
    cat(sprintf("Held-out validation (%s): n = %d, R2 = %.3f, RMSE = %.3f, r = %.3f\n", what, length(y_te),
                validation$r2, validation$rmse, validation$correlation))
  } else if (!ml_method) {
    cat("Held-out validation: no QTL found, so there is no model to validate.\n")
  }
}

# ============================================================
# OUTPUT PATHS
# ============================================================

manhattan_file <- file.path(results_dir, paste0("manhattan_", method, ".png"))
ideogram_file  <- file.path(results_dir, paste0("ideogram_", method, ".png"))
qtl_file       <- file.path(results_dir, paste0("qtl_significant_", method, ".csv"))
full_file      <- file.path(results_dir, paste0("qtl_full_scan_", method, ".csv"))

plot_cfg <- config$plot

# ============================================================
# MANHATTAN PNG
# ============================================================

png(manhattan_file, width = as.numeric(plot_cfg$manhattan_width),
    height = as.numeric(plot_cfg$manhattan_height), res = 150, bg = "white")
draw_manhattan(scan_all, sig, thr, plot_cfg,
               title = if (ml_method) paste0(method, " - ", trait, "  |  importance score (1 = permutation threshold)  |  ", nrow(sig), " QTL")
                       else paste0(method, " - ", trait, "  |  threshold LOD ", signif(thr, 3), "  |  ", nrow(sig), " QTL"),
               ylab = if (ml_method) "Importance / threshold" else "LOD")
dev.off()

# ============================================================
# IDEOGRAM
# ============================================================

png(ideogram_file, width = as.numeric(plot_cfg$ideogram_width),
    height = as.numeric(plot_cfg$ideogram_height), res = 150, bg = "white")

qtl_list <- list()
qtl_list[[method]] <- sig

make_chromosome_ideogram_reference(
  qtl_list, all_markers, thr,
  title_str = paste(method, "- Trait:", trait),
  band_dark_col = plot_cfg$chr_dark, band_light_col = plot_cfg$chr_light,
  border_col = plot_cfg$chr_border,
  qtl_size_min = as.numeric(plot_cfg$qtl_marker_size_min),
  qtl_size_max = as.numeric(plot_cfg$qtl_marker_size_max)
)

dev.off()

# ============================================================
# CSV
# ============================================================

write.csv(sig, qtl_file, row.names = FALSE)
write.csv(scan_all, full_file, row.names = FALSE)

# ============================================================
# RESULT JSON
# ============================================================

result <- list(
  status = "success",
  summary = list(
    samples = nrow(X_train),
    markers = ncol(X),
    chromosomes = length(unique(marker_info$Chromosome)),
    significant_qtl = nrow(sig),
    method = method,
    trait = trait,
    threshold_mode = threshold_mode,
    threshold_lod = thr,
    permutations = if (!is.null(ml_info)) ml_info$permutations else if (threshold_mode == "permutation" && linear_method) n_perm else NULL,
    threshold_note = threshold_note,
    lod_drop = lod_drop,
    peak_drop = max(peak_drop, lod_drop),
    p_threshold = p_threshold,
    cofactors = if (is.null(cofactors_used)) list() else as.list(cofactors_used),
    score_type = if (ml_method) "importance / permutation threshold (not a LOD)" else "LOD",
    ml = ml_info,
    oob_r2 = oob_r2,
    validation = validation,
    missing_percent = round(100 * geno_res$missing_fraction, 2),
    imputed_genotypes = n_imputed
  ),
  manhattan = normalizePath(manhattan_file, winslash = "/"),
  ideogram  = normalizePath(ideogram_file, winslash = "/"),
  qtl_csv   = normalizePath(qtl_file, winslash = "/"),
  full_csv  = normalizePath(full_file, winslash = "/"),
  qtl_table = if (nrow(sig) > 0) {
    jsonlite::fromJSON(jsonlite::toJSON(sig, dataframe = "rows", na = "null", auto_unbox = TRUE))
  } else {
    list()
  }
)

cat(jsonlite::toJSON(result, auto_unbox = TRUE, na = "null"))
