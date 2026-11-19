# ============================================================
# TerraTools QTL engine
#
# QTL are reported as *regions*: one row per QTL peak, with the
# flanking markers on the left and right and a LOD-drop support
# interval, following the standard QTL-mapping literature:
#
#   Lander & Botstein 1989   interval mapping, 1-LOD / 1.5-LOD support interval
#   Haley & Knott 1992       regression mapping, LOD = (n/2) log10(RSS0/RSS1)
#   Churchill & Doerge 1994  permutation test for the genome-wide LOD threshold
#   Zeng 1994, Jansen & Stam 1994
#                            composite interval / multiple-QTL mapping (cofactors,
#                            cofactors inside the scan window are dropped)
#   Manichaikul et al. 2006  Bayesian credible interval
#   Broman & Sen 2009 (R/qtl lodint, bayesint) and Broman et al. 2019
#   (R/qtl2 find_peaks: threshold + peakdrop)
#                            peak calling, intervals expanded to flanking markers
# ============================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

required_packages <- c(
  "randomForest", "caret", "ggplot2",
  "dplyr", "glmnet", "scales", "e1071", "xgboost",
  "writexl", "jsonlite"
)

missing_packages <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them once ahead of time (e.g. in your Docker image / setup script) with:\n",
    "  Rscript -e 'install.packages(c(",
    paste0('"', missing_packages, '"', collapse = ", "),
    "))'\n",
    "Live install.packages() inside a request is disabled because it can hang ",
    "indefinitely (e.g. no outbound internet, no CRAN binary) with no visible error."
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

# ============================================================
# MODELS
# ============================================================
MAPPING_MODELS <- list(
  "SIM"     = list(label = "Simple Interval Mapping (SIM)",    desc = "Single-marker regression; LOD = (n/2)log10(RSS0/RSS1)"),
  "CIM"     = list(label = "Composite Interval Mapping (CIM)", desc = "Regression with stepwise-selected background cofactors"),
  "MQM"     = list(label = "Multiple QTL Model (MQM)",         desc = "Regression with a larger multi-QTL background model"),
  "RF"      = list(label = "Random Forest QTL",                desc = "Non-parametric; permutation importance"),
  "SVR"     = list(label = "Support Vector Regression",        desc = "Kernel-based; permutation importance"),
  "XGBOOST" = list(label = "XGBoost QTL",                     desc = "Gradient boosting; gain importance")
)

LINEAR_METHODS <- c("SIM", "CIM", "MQM")

MODEL_COLORS <- c(
  "SIM"     = "#CC0000",
  "CIM"     = "#FF3333",
  "MQM"     = "#990000",
  "RF"      = "#FF6600",
  "SVR"     = "#CC3300",
  "XGBOOST" = "#FF0066"
)

# ============================================================
# HELPERS
# ============================================================
cn <- function(x) trimws(gsub('[\r\n\"]', '', as.character(x)))

# Chromosome names in map order: 1, 2, ..., 10, 11 (not "1","10","11","2"),
# followed by any non-numeric names (X, Y, Un) alphabetically.
chr_order <- function(x) {
  u <- unique(as.character(x))
  n <- suppressWarnings(as.numeric(u))
  u[order(is.na(n), n, u)]
}

impute_mean <- function(M) {
  M <- as.matrix(M)
  for (j in seq_len(ncol(M))) {
    na <- is.na(M[, j])
    if (any(na)) {
      mu <- if (all(na)) 0 else mean(M[!na, j])
      M[na, j] <- mu
    }
  }
  M
}

impute_mode <- function(M) {
  M <- as.matrix(M)
  for (j in seq_len(ncol(M))) {
    na <- is.na(M[, j])
    if (any(na)) {
      if (all(na)) { M[, j] <- 0; next }
      tab <- table(M[!na, j])
      M[na, j] <- as.numeric(names(tab)[which.max(tab)])
    }
  }
  M
}

# ============================================================
# GENOTYPE LOADER
#   rows 1-3: marker names, chromosome, position (cM); then one
#   row per individual.  Any negative code (-1, -9 ...) or a
#   blank / non-numeric cell is a MISSING genotype.
# ============================================================
load_genotype_file <- function(filepath) {
  lines <- readLines(filepath, warn = FALSE)
  lines <- gsub("\r", "", lines)
  lines <- lines[nchar(trimws(lines)) > 0]
  if (length(lines) < 4) stop("Need >= 4 rows: header + chr + pos + >= 1 sample")

  header     <- cn(strsplit(lines[1], ",")[[1]])
  chr_row    <- cn(strsplit(lines[2], ",")[[1]])
  pos_row    <- cn(strsplit(lines[3], ",")[[1]])
  data_lines <- lines[-(1:3)]

  gd <- read.csv(text = paste(data_lines, collapse = "\n"),
                 header = FALSE, stringsAsFactors = FALSE, quote = '"', fill = TRUE)
  if (ncol(gd) < length(header)) {
    for (k in (ncol(gd)+1):length(header)) gd[[paste0("V",k)]] <- NA
  }
  if (ncol(gd) > length(header)) gd <- gd[, 1:length(header)]
  colnames(gd) <- header
  rownames(gd) <- cn(as.character(gd[, 1]))
  gd <- gd[, -1, drop = FALSE]
  for (col in colnames(gd)) {
    v <- suppressWarnings(as.numeric(as.character(gd[[col]])))
    v[!is.na(v) & v < 0] <- NA            # -1 / -9 style missing codes
    gd[[col]] <- v
  }

  marker_info <- data.frame(
    Marker     = colnames(gd),
    Chromosome = cn(chr_row[-1][seq_along(colnames(gd))]),
    Position   = suppressWarnings(as.numeric(cn(pos_row[-1][seq_along(colnames(gd))]))),
    stringsAsFactors = FALSE
  )
  marker_info <- marker_info[!is.na(marker_info$Position) & marker_info$Chromosome != "", ]
  # map order: chromosome (numeric aware), then cM position
  marker_info <- marker_info[order(match(marker_info$Chromosome, chr_order(marker_info$Chromosome)),
                                   marker_info$Position), ]
  valid_cols  <- colnames(gd)[colnames(gd) %in% marker_info$Marker]
  gd          <- gd[, valid_cols, drop = FALSE]

  list(geno = gd, marker_info = marker_info,
       n_samples = nrow(gd), n_markers = ncol(gd),
       missing_fraction = mean(is.na(as.matrix(gd))),
       chromosomes = chr_order(marker_info$Chromosome))
}

# ============================================================
# MISSING GENOTYPES
#   A missing genotype is filled with its most probable class given the nearest
#   observed markers on each side of the same individual (Kosambi distances,
#   backcross / RIL 2-state or intercross 3-state transition matrix), the
#   same flanking-marker information R/qtl uses (Haley-Knott).  Without this each
#   marker would be tested on a different subset of individuals, which makes
#   the LOD curve jagged and can split one QTL into several fake peaks.
# ============================================================
kosambi_r <- function(d_cM) 0.5 * tanh(d_cM / 50)

impute_flanking <- function(X, marker_info) {
  X <- as.matrix(X)
  cls <- sort(unique(as.vector(X[!is.na(X)])))
  K <- length(cls)
  if (K < 2 || K > 3 || !anyNA(X)) return(list(X = X, n_imputed = 0L))
  mi <- marker_info[match(colnames(X), marker_info$Marker), , drop = FALSE]
  trans <- function(r) {
    if (K == 2) matrix(c(1 - r, r, r, 1 - r), 2, byrow = TRUE)
    else matrix(c((1 - r)^2, 2 * r * (1 - r), r^2,
                  r * (1 - r), (1 - r)^2 + r^2, r * (1 - r),
                  r^2, 2 * r * (1 - r), (1 - r)^2), 3, byrow = TRUE)
  }
  n_imp <- 0L
  for (chr in unique(as.character(mi$Chromosome))) {
    cols <- which(as.character(mi$Chromosome) == chr)
    cols <- cols[order(mi$Position[cols])]
    pos  <- mi$Position[cols]
    sub  <- X[, cols, drop = FALSE]
    for (i in which(rowSums(is.na(sub)) > 0)) {
      g <- sub[i, ]
      obs <- which(!is.na(g))
      if (!length(obs)) next
      for (j in which(is.na(g))) {
        l  <- obs[obs < j]; l <- if (length(l)) l[length(l)] else NA
        rr <- obs[obs > j]; rr <- if (length(rr)) rr[1] else NA
        post <- rep(1, K)
        if (!is.na(l))  post <- post * trans(kosambi_r(pos[j] - pos[l]))[match(g[l], cls), ]
        if (!is.na(rr)) post <- post * trans(kosambi_r(pos[rr] - pos[j]))[, match(g[rr], cls)]
        sub[i, j] <- cls[which.max(post)]
        n_imp <- n_imp + 1L
      }
    }
    X[, cols] <- sub
  }
  list(X = X, n_imputed = n_imp)
}

# ============================================================
# REGRESSION CORE
#   LOD  = (n/2) * log10(RSS0 / RSS1)            (Haley & Knott 1992)
#   P    = F-test of the extra genotype terms     (nested models)
# ============================================================
geno_design <- function(g) {
  lv <- sort(unique(g))
  if (length(lv) > 5) return(matrix(g, ncol = 1))          # continuous dosage
  if (length(lv) == 2) return(matrix(as.numeric(g == lv[2]), ncol = 1))
  vapply(lv[-1], function(l) as.numeric(g == l), numeric(length(g)))   # a + d for an F2
}

lr_test <- function(y, X0, X1) {
  n  <- length(y)
  f0 <- lm.fit(X0, y)
  f1 <- lm.fit(X1, y)
  rss0 <- sum(f0$residuals^2)
  rss1 <- sum(f1$residuals^2)
  df1  <- f1$rank - f0$rank
  df2  <- n - f1$rank
  if (!is.finite(rss0) || !is.finite(rss1) || rss0 <= 0 || df2 < 1) return(NULL)
  if (df1 < 1 || rss1 <= 0) return(list(lod = 0, p = 1, rss0 = rss0, rss1 = rss0, df = 0, n = n))
  lod <- max(0, n / 2 * log10(rss0 / rss1))
  Fv  <- ((rss0 - rss1) / df1) / (rss1 / df2)
  list(lod = lod, p = pf(Fv, df1, df2, lower.tail = FALSE),
       rss0 = rss0, rss1 = rss1, df = df1, n = n)
}

# Additive / dominance effects of a marker.  a = half the difference between
# the two homozygote classes, d = heterozygote deviation from their mean.
#   pve        = marginal R^2 of the marker (single-marker model)
#   a, se, d   = effects estimated in the model that also contains `cov` (the
#                CIM / MQM background cofactors), so a QTL partly masked by a
#                linked one gets its conditional effect, as in composite
#                interval mapping.  Positive a = the higher-coded parent's
#                allele raises the trait.
effect_stats <- function(yv, gv, cov = NULL) {
  res <- list(pve = NA_real_, a = NA_real_, se = NA_real_, d = NA_real_,
              ratio = NA_real_, cls = "Additive")
  lv <- sort(unique(gv)); nl <- length(lv)
  if (nl < 2) return(res)
  d_code <- NULL
  if (nl > 3) {
    a_code <- gv
  } else {
    lo <- lv[1]; hi <- lv[nl]
    a_code <- (gv - (lo + hi) / 2) / ((hi - lo) / 2)
    if (nl == 3) d_code <- as.numeric(gv == lv[2])
  }
  Xm <- cbind(rep(1, length(yv)), a_code, d_code)
  tss <- sum((yv - mean(yv))^2)
  fm <- lm.fit(Xm, yv)
  if (tss > 0) res$pve <- max(0, (1 - sum(fm$residuals^2) / tss) * 100)
  Xf <- if (is.null(cov)) Xm else cbind(Xm, cov)
  ff <- lm.fit(Xf, yv)
  p  <- ncol(Xf)
  res$a <- unname(ff$coefficients[2])
  if (ff$rank == p && length(yv) > p) {
    s2 <- sum(ff$residuals^2) / (length(yv) - p)
    V  <- s2 * chol2inv(ff$qr$qr[1:p, 1:p, drop = FALSE])
    res$se <- sqrt(V[2, 2])
  }
  if (!is.null(d_code)) {
    res$d <- unname(ff$coefficients[3])
    if (is.finite(res$a) && res$a != 0 && is.finite(res$d)) {
      res$ratio <- abs(res$d / res$a)
      res$cls <- if (res$ratio > 1) "Overdominance"
                 else if (res$ratio > 0.5) "Dominance"
                 else if (res$ratio > 0.2) "Partial Dom." else "Additive"
    }
  }
  res
}

rows_to_df <- function(rows) {
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  cols <- names(rows[[1]])
  out <- lapply(cols, function(k) unlist(lapply(rows, function(r) r[[k]]), use.names = FALSE))
  names(out) <- cols
  as.data.frame(out, stringsAsFactors = FALSE)
}

# ============================================================
# LINEAR SCAN (SIM / CIM / MQM share this)
#   cofactors = background markers; those within `window_cM` of the
#   tested position (same chromosome) are dropped from the model so
#   they cannot absorb the QTL being tested (Zeng 1994).
# ============================================================
qtl_scan_linear <- function(X, y, marker_info, cofactors = character(0),
                            window_cM = 10, method = "SIM") {
  mi <- marker_info[match(colnames(X), marker_info$Marker), , drop = FALSE]
  Xc <- NULL; cof_chr <- cof_pos <- NULL
  if (length(cofactors)) {
    Xc <- impute_mean(X[, cofactors, drop = FALSE])
    cm <- marker_info[match(cofactors, marker_info$Marker), , drop = FALSE]
    cof_chr <- as.character(cm$Chromosome); cof_pos <- cm$Position
  }
  rows <- vector("list", ncol(X))
  for (j in seq_len(ncol(X))) {
    g  <- X[, j]
    ok <- !is.na(g) & !is.na(y)
    if (sum(ok) < 8 || length(unique(g[ok])) < 2) next
    chr <- as.character(mi$Chromosome[j]); pos <- mi$Position[j]
    if (is.na(pos)) next
    use <- integer(0)
    if (length(cofactors)) {
      use <- which(cofactors != colnames(X)[j] &
                   !(cof_chr == chr & abs(cof_pos - pos) <= window_cM))
    }
    yv <- y[ok]; gv <- g[ok]
    X0 <- cbind(rep(1, length(yv)), if (length(use)) Xc[ok, use, drop = FALSE])
    tst <- lr_test(yv, X0, cbind(X0, geno_design(gv)))
    if (is.null(tst)) next
    ef <- effect_stats(yv, gv, if (length(use)) Xc[ok, use, drop = FALSE] else NULL)
    tss <- sum((yv - mean(yv))^2)
    rows[[j]] <- list(
      Marker = colnames(X)[j], Chromosome = chr, Position = pos,
      LOD = tst$lod, P_Value = tst$p,
      PVE_Marker_Percent = ef$pve,
      # variance the marker explains over and above the cofactors (% of total)
      PVE_CIM_Percent = if (length(use) && tss > 0) max(0, (tst$rss0 - tst$rss1) / tss * 100) else ef$pve,
      Additive_Effect = ef$a, SE_Additive = ef$se,
      CI_Lower = ef$a - 1.96 * ef$se, CI_Upper = ef$a + 1.96 * ef$se,
      Dominance_Effect = ef$d, Dominance_Ratio = ef$ratio, Dominance_Class = ef$cls,
      N_Cofactors = length(use), N_Samples = length(yv), Method = method
    )
  }
  res <- rows_to_df(rows)
  if (is.null(res)) return(data.frame(Message = paste("No results from", method)))
  res <- res[order(match(res$Chromosome, chr_order(res$Chromosome)), res$Position), , drop = FALSE]
  rownames(res) <- NULL
  res
}

# Forward selection with backward elimination of background markers.
# Candidates are screened by their correlation with the residuals of the
# current model, then tested exactly with a partial F-test.
select_cofactors <- function(X, y, marker_info, max_cof = 5, p_enter = 0.01,
                             p_remove = 0.05, min_dist = 10) {
  if (max_cof < 1) return(character(0))
  keep <- !is.na(y)
  Xi <- impute_mean(X[keep, , drop = FALSE]); yy <- y[keep]
  n <- length(yy)
  mi  <- marker_info[match(colnames(Xi), marker_info$Marker), , drop = FALSE]
  chr <- as.character(mi$Chromosome); pos <- mi$Position
  varying <- apply(Xi, 2, sd) > 0
  one <- rep(1, n)
  sel <- integer(0)

  for (step in seq_len(max_cof)) {
    X0  <- cbind(one, Xi[, sel, drop = FALSE])
    res <- lm.fit(X0, yy)$residuals
    sc  <- suppressWarnings(abs(as.numeric(cor(Xi, res))))
    sc[is.na(sc) | !varying] <- 0
    for (s in sel) sc[chr == chr[s] & abs(pos - pos[s]) <= min_dist] <- 0
    sc[sel] <- 0
    if (!any(sc > 0)) break
    cand <- order(sc, decreasing = TRUE)[seq_len(min(30, sum(sc > 0)))]
    best <- NA_integer_; best_p <- 1
    for (cj in cand) {
      t <- lr_test(yy, X0, cbind(X0, Xi[, cj]))
      if (!is.null(t) && t$p < best_p) { best_p <- t$p; best <- cj }
    }
    if (is.na(best) || best_p >= p_enter) break
    sel <- c(sel, best)
  }

  while (length(sel) > 1) {
    ps <- vapply(seq_along(sel), function(k) {
      X0 <- cbind(one, Xi[, sel[-k], drop = FALSE])
      t  <- lr_test(yy, X0, cbind(X0, Xi[, sel[k]]))
      if (is.null(t)) 1 else t$p
    }, numeric(1))
    if (max(ps) > p_remove) sel <- sel[-which.max(ps)] else break
  }
  colnames(Xi)[sel]
}

qtl_scan_sim <- function(X, y, marker_info) {
  qtl_scan_linear(X, y, marker_info, character(0), 0, "SIM")
}

qtl_scan_cim <- function(X, y, marker_info, window_cM = 10, n_cofactors = 5) {
  # cofactors closer together than the window would be dropped together anyway
  cof <- select_cofactors(X, y, marker_info, max_cof = n_cofactors, min_dist = window_cM)
  cat("CIM cofactors (", length(cof), "): ", paste(cof, collapse = ", "), "\n", sep = "")
  res <- qtl_scan_linear(X, y, marker_info, cof, window_cM, "CIM")
  attr(res, "cofactors") <- cof
  res
}

# MQM-style: a larger background model (up to 15 cofactors) and the wider
# 25 cM window that R/qtl's mqmscan uses by default.
qtl_scan_mqm <- function(X, y, marker_info, window_cM = 25, max_cof = 15) {
  cof <- select_cofactors(X, y, marker_info, max_cof = max_cof, p_enter = 0.01, p_remove = 0.02,
                           min_dist = window_cM)
  cat("MQM cofactors (", length(cof), "): ", paste(cof, collapse = ", "), "\n", sep = "")
  res <- qtl_scan_linear(X, y, marker_info, cof, window_cM, "MQM")
  attr(res, "cofactors") <- cof
  res
}

# ============================================================
# MACHINE-LEARNING SCANS
#   Marker importance is rescaled to 0..scale and used as a LOD-like score.
#   It is NOT a likelihood LOD; P_Value = 10^-score is only a monotone
#   transform so that the same peak / interval code can be used.
#   Every marker gets a row (score 0 when unimportant) so the profile is
#   continuous along the chromosome.
# ============================================================
ml_result <- function(imp, X, y, marker_info, method, scale) {
  imp$Importance[!is.finite(imp$Importance)] <- 0
  imp$Importance <- pmax(imp$Importance, 0)
  mx <- max(imp$Importance)
  if (!is.finite(mx) || mx <= 0) return(data.frame(Message = "No important markers"))
  imp$LOD <- imp$Importance / mx * scale
  imp <- merge(imp, marker_info, by = "Marker", all.x = TRUE)
  imp <- imp[!is.na(imp$Position), , drop = FALSE]
  rows <- lapply(seq_len(nrow(imp)), function(i) {
    m  <- imp$Marker[i]
    g  <- X[, m]; v <- !is.na(g)
    ef <- if (sum(v) >= 5) effect_stats(y[v], g[v]) else list(pve = NA_real_, a = NA_real_, se = NA_real_, d = NA_real_, ratio = NA_real_, cls = "Additive")
    list(Marker = m, Chromosome = as.character(imp$Chromosome[i]), Position = imp$Position[i],
         LOD = imp$LOD[i], P_Value = 10^(-imp$LOD[i]),
         PVE_Marker_Percent = ef$pve, PVE_CIM_Percent = ef$pve,
         Additive_Effect = ef$a, SE_Additive = ef$se,
         CI_Lower = ef$a - 1.96 * ef$se, CI_Upper = ef$a + 1.96 * ef$se,
         Dominance_Effect = ef$d, Dominance_Ratio = ef$ratio, Dominance_Class = ef$cls,
         N_Cofactors = 0, N_Samples = length(y), Method = method)
  })
  res <- rows_to_df(rows)
  res <- res[order(match(res$Chromosome, chr_order(res$Chromosome)), res$Position), , drop = FALSE]
  rownames(res) <- NULL
  res
}

qtl_scan_rf <- function(X, y, marker_info) {
  X0 <- X                       # keep NAs for the effect estimates
  X <- impute_mode(X)
  rf_fit <- tryCatch(randomForest(x = X, y = y, ntree = 500, importance = TRUE), error = function(e) NULL)
  if (is.null(rf_fit)) return(data.frame(Message = "Random Forest failed"))
  imp <- importance(rf_fit, type = 1)[, 1]
  ml_result(data.frame(Marker = names(imp), Importance = as.numeric(imp)), X0, y, marker_info, "RF", 10)
}

qtl_scan_svr <- function(X, y, marker_info) {
  X0 <- X
  X <- impute_mode(X)
  svr_fit <- tryCatch(svm(x = X, y = y, kernel = "radial", scale = TRUE), error = function(e) NULL)
  if (is.null(svr_fit)) return(data.frame(Message = "SVR failed"))
  baseline_error <- mean((y - predict(svr_fit, X))^2)
  importance <- vapply(colnames(X), function(m) {
    Xp <- X; Xp[, m] <- sample(Xp[, m])
    (mean((y - predict(svr_fit, Xp))^2) - baseline_error) / baseline_error
  }, numeric(1))
  ml_result(data.frame(Marker = names(importance), Importance = as.numeric(importance)), X0, y, marker_info, "SVR", 10)
}

qtl_scan_xgboost <- function(X, y, marker_info) {
  X0 <- X
  X <- impute_mode(X)
  dtrain <- xgb.DMatrix(data = as.matrix(X), label = y)
  xgb_fit <- tryCatch(xgb.train(
    data = dtrain, nrounds = 100, objective = "reg:squarederror", verbose = 0,
    params = list(max_depth = 3, eta = 0.1, subsample = 0.8, colsample_bytree = 0.8)),
    error = function(e) NULL)
  if (is.null(xgb_fit)) return(data.frame(Message = "XGBoost failed"))
  im  <- xgb.importance(model = xgb_fit)
  imp <- data.frame(Marker = colnames(X), Importance = 0, stringsAsFactors = FALSE)
  imp$Importance[match(im$Feature, imp$Marker)] <- im$Gain
  ml_result(imp, X0, y, marker_info, "XGBOOST", 12)
}

# ============================================================
# MAIN DISPATCHER
# ============================================================
qtl_scan <- function(genotype_matrix, phenotype_vector, marker_info,
                     method = "CIM", n_cofactors = 5, window_size = 10) {
  cat("\n=== QTL SCAN === Method:", method, "\n")
  valid_rows <- !is.na(phenotype_vector)
  X <- as.matrix(genotype_matrix[valid_rows, , drop = FALSE])
  y <- phenotype_vector[valid_rows]
  nzv <- apply(X, 2, function(x) length(unique(na.omit(x))) < 2)
  X   <- X[, !nzv, drop = FALSE]
  cat("Individuals:", nrow(X), " Markers after NZV removal:", ncol(X), "\n")
  if (method == "SIM")     return(qtl_scan_sim(X, y, marker_info))
  if (method == "CIM")     return(qtl_scan_cim(X, y, marker_info, window_size, n_cofactors))
  if (method == "MQM")     return(qtl_scan_mqm(X, y, marker_info))
  if (method == "RF")      return(qtl_scan_rf(X, y, marker_info))
  if (method == "SVR")     return(qtl_scan_svr(X, y, marker_info))
  if (method == "XGBOOST") return(qtl_scan_xgboost(X, y, marker_info))
  data.frame(Message = paste("Unknown method:", method))
}

# ============================================================
# GENOME-WIDE PERMUTATION THRESHOLD  (Churchill & Doerge 1994)
#   Shuffle the phenotype, record the largest single-marker-regression LOD
#   anywhere in the genome, repeat; the (1 - alpha) quantile of those maxima
#   is the genome-wide LOD threshold.  Computed with matrix algebra
#   (RSS from class sums), so 1000 permutations take seconds.
# ============================================================
perm_max_lod <- function(X, y, n_perm = 1000, seed = 42) {
  ok <- !is.na(y)
  X <- as.matrix(X[ok, , drop = FALSE]); y <- y[ok]
  n <- length(y)
  vals <- sort(unique(as.vector(X[!is.na(X)])))
  categorical <- length(vals) <= 3
  M  <- (!is.na(X)) * 1
  Nm <- colSums(M)
  Nsafe <- pmax(Nm, 1)
  if (categorical) {
    Ig <- lapply(vals, function(v) { m <- (X == v); m[is.na(m)] <- FALSE; m * 1 })
    Ng <- lapply(Ig, colSums)
  } else {
    Xz <- X; Xz[is.na(Xz)] <- 0
    Sx <- colSums(Xz); Sxx <- colSums(Xz^2)
    Sxx_c <- Sxx - Sx^2 / Nsafe
  }
  usable <- Nm >= 8
  set.seed(seed)
  out <- numeric(n_perm); done <- 0
  while (done < n_perm) {
    B  <- min(200, n_perm - done)
    Y  <- vapply(seq_len(B), function(b) sample(y), numeric(n))
    SY <- crossprod(M, Y); SY2 <- crossprod(M, Y^2)
    rss0 <- SY2 - SY^2 / Nsafe
    if (categorical) {
      ssb <- 0
      for (k in seq_along(Ig)) {
        ng <- Ng[[k]]; ng[ng == 0] <- Inf
        ssb <- ssb + crossprod(Ig[[k]], Y)^2 / ng
      }
      rss1 <- SY2 - ssb
    } else {
      Sxy <- crossprod(Xz, Y) - Sx * SY / Nsafe
      den <- Sxx_c; den[den <= 0] <- Inf
      rss1 <- rss0 - Sxy^2 / den
    }
    lod <- Nm / 2 * log10(rss0 / rss1)
    lod[!is.finite(lod)] <- 0
    lod[!usable, ] <- 0
    out[done + seq_len(B)] <- apply(lod, 2, max)
    done <- done + B
  }
  out
}

# ------------------------------------------------------------
# PERMUTATION THRESHOLD FOR CIM / MQM  (Doerge & Churchill 1996)
#   The whole procedure is repeated on every shuffled phenotype: the background
#   cofactors are selected again and the composite scan is run again; the largest
#   LOD of the scan is kept.  The single-marker null (above) is NOT used for them,
#   because with cofactors the maxima are larger (the cofactors absorb chance
#   associations and leave less residual variance), which made the cut-off too low.
#
#   Only the maximum is needed, so the scan is computed in blocks: the markers that
#   share the same cofactor set (all cofactors outside their window) are regressed
#   together on the residuals of the cofactor model.  Markers that cannot be handled
#   that way (missing calls, more than 3 genotype classes) use the same lr_test() as
#   qtl_scan_linear(), so the maximum equals max(qtl_scan_*()$LOD) exactly.
# ------------------------------------------------------------
# Genotype classes of every marker, worked out once (they do not change between permutations):
# the 0/1 indicator of the 2nd class (and of the 3rd class for 3-class markers).
cim_prep <- function(X) {
  n <- nrow(X); m <- ncol(X)
  nlev <- integer(m); D1 <- D2 <- matrix(0, n, m)
  for (j in seq_len(m)) {
    g <- X[, j]; lv <- sort(unique(g[!is.na(g)])); nlev[j] <- length(lv)
    if (anyNA(g) || !(length(lv) %in% 2:3)) next
    D1[, j] <- as.numeric(g == lv[2])
    if (length(lv) == 3) D2[, j] <- as.numeric(g == lv[3])
  }
  list(nlev = nlev, fast = colSums(is.na(X)) == 0 & nlev %in% 2:3, D1 = D1, D2 = D2)
}

cim_max_lod <- function(X, y, marker_info, cof, window_cM, prep = cim_prep(X)) {
  n <- length(y)
  mi  <- marker_info[match(colnames(X), marker_info$Marker), , drop = FALSE]
  chr <- as.character(mi$Chromosome); pos <- mi$Position
  k   <- length(cof)
  Xc  <- if (k) impute_mean(X[, cof, drop = FALSE]) else NULL
  # code[j]: which cofactors are left out for marker j (itself, or on its chromosome within the window)
  code <- rep(0, ncol(X))
  if (k) {
    cm <- marker_info[match(cof, marker_info$Marker), , drop = FALSE]
    cof_chr <- as.character(cm$Chromosome); cof_pos <- cm$Position
    excl <- outer(seq_len(k), seq_len(ncol(X)), function(a, j) {
      cof[a] == colnames(X)[j] | (cof_chr[a] == chr[j] & abs(cof_pos[a] - pos[j]) <= window_cM) })
    excl[is.na(excl)] <- FALSE
    code <- colSums(excl * 2^(seq_len(k) - 1))
  }
  best <- 0
  for (cd in unique(code)) {
    cols <- which(code == cd & !is.na(pos))
    if (!length(cols)) next
    use <- if (k) which(bitwAnd(as.integer(cd), as.integer(2^(seq_len(k) - 1))) == 0) else integer(0)
    X0  <- cbind(rep(1, n), if (length(use)) Xc[, use, drop = FALSE])
    fast <- cols[prep$fast[cols]]
    if (length(fast)) {
      qrx <- qr(X0); r0 <- qr.resid(qrx, y); rss0 <- sum(r0^2)
      for (nl in 2:3) {
        idx <- fast[prep$nlev[fast] == nl]
        if (!length(idx)) next
        D1 <- qr.resid(qrx, prep$D1[, idx, drop = FALSE])
        u  <- colSums(D1 * r0); A <- colSums(D1^2)
        if (nl == 2) {
          expl <- ifelse(A > 1e-8 * n, u^2 / A, 0)
        } else {
          D2 <- qr.resid(qrx, prep$D2[, idx, drop = FALSE])
          v <- colSums(D2 * r0); C <- colSums(D2^2); B <- colSums(D1 * D2)
          det <- A * C - B^2
          expl <- ifelse(det > 1e-8 * n^2, (u^2 * C - 2 * u * v * B + v^2 * A) / det, 0)
        }
        rss1 <- pmax(rss0 - expl, 0)
        ok   <- rss1 > 0 & expl > 0
        if (any(ok)) best <- max(best, max(n / 2 * log10(rss0 / rss1[ok])))
      }
    }
    for (j in setdiff(cols, fast)) {                     # slow path: the exact code of qtl_scan_linear()
      g <- X[, j]; ok <- !is.na(g) & !is.na(y)
      if (sum(ok) < 8 || length(unique(g[ok])) < 2) next
      X0j <- cbind(rep(1, sum(ok)), if (length(use)) Xc[ok, use, drop = FALSE])
      tst <- lr_test(y[ok], X0j, cbind(X0j, geno_design(g[ok])))
      if (!is.null(tst)) best <- max(best, tst$lod)
    }
  }
  best
}

perm_max_lod_cim <- function(X, y, marker_info, method = "CIM", n_perm = 1000, window_cM = 10,
                             n_cofactors = 5, seed = 42, progress = TRUE) {
  ok <- !is.na(y)
  X <- as.matrix(X[ok, , drop = FALSE]); y <- y[ok]
  X <- X[, apply(X, 2, function(x) length(unique(na.omit(x))) >= 2), drop = FALSE]     # as in qtl_scan()
  if (method == "MQM") { window_cM <- 25; max_cof <- 15; p_remove <- 0.02 } else { max_cof <- n_cofactors; p_remove <- 0.05 }
  prep <- cim_prep(X)
  set.seed(seed)
  out <- numeric(n_perm); t0 <- Sys.time()
  for (b in seq_len(n_perm)) {
    yp  <- sample(y)
    cof <- select_cofactors(X, yp, marker_info, max_cof = max_cof, p_enter = 0.01, p_remove = p_remove, min_dist = window_cM)
    out[b] <- cim_max_lod(X, yp, marker_info, cof, window_cM, prep)
    if (progress && (b %% 50 == 0 || b == n_perm)) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      cat(sprintf("  permutation %d / %d  (%.0f s elapsed, about %.0f s left)\n", b, n_perm, el, el / b * (n_perm - b)))
      flush.console()
    }
  }
  out
}

# ============================================================
# QTL REGIONS: peaks, flanking markers, support intervals
#
# 1. Per chromosome find local LOD maxima that pass the threshold
#    (LOD >= thr and nominal P <= p_max).
# 2. Two maxima are separate QTL only if the LOD curve falls at least
#    `peak_drop` below the lower one between them (R/qtl2 find_peaks);
#    otherwise the lower one is a shoulder of the higher peak.
# 3. Support interval (R/qtl lodint): walk left and right from the peak
#    while LOD >= peak - lod_drop.  The first marker outside on each side is
#    the flanking marker (expandtomarkers = TRUE), so the region always
#    covers the interval.  CI_Start/End are the interpolated positions where
#    the curve crosses peak - lod_drop.
# 4. Bayesian credible interval (R/qtl bayesint): posterior ~ 10^LOD x marker
#    spacing, restricted to this peak's domain.
# ============================================================
REGION_COLUMNS <- c(
  "QTL_ID", "Chromosome",
  "Left_Marker", "Left_Position", "Peak_Marker", "Peak_Position", "Right_Marker", "Right_Position",
  "Interval_Size_cM", "CI_Start_cM", "CI_End_cM", "Bayes_Start_cM", "Bayes_End_cM",
  "N_Markers_In_Region", "LOD", "P_Value", "Genomewide_P",
  "PVE_Marker_Percent", "PVE_CIM_Percent", "Additive_Effect", "SE_Additive",
  "CI_Lower", "CI_Upper", "Dominance_Effect", "Dominance_Ratio", "Dominance_Class",
  "N_Cofactors", "N_Samples", "Method", "Threshold_LOD", "LOD_Drop"
)

empty_regions <- function() {
  chr <- c("QTL_ID", "Chromosome", "Left_Marker", "Peak_Marker", "Right_Marker", "Dominance_Class", "Method")
  out <- lapply(REGION_COLUMNS, function(k) if (k %in% chr) character(0) else numeric(0))
  names(out) <- REGION_COLUMNS
  as.data.frame(out, stringsAsFactors = FALSE)
}

interp_cross <- function(x_out, y_out, x_in, y_in, target) {
  if (!is.finite(y_in) || !is.finite(y_out) || y_in == y_out) return(x_in)
  x <- x_out + (target - y_out) * (x_in - x_out) / (y_in - y_out)
  min(max(x, min(x_out, x_in)), max(x_out, x_in))
}

call_qtl_regions <- function(scan_all, thr, lod_drop = 1.5, peak_drop = 2, p_max = 1,
                             bayes_prob = 0.95, bayes = TRUE, perm_null = NULL) {
  peak_drop <- max(peak_drop, lod_drop)
  sc <- scan_all[!is.na(scan_all$LOD) & !is.na(scan_all$Position), , drop = FALSE]
  out <- list()
  keep_cols <- c("Marker", "Position", "LOD", "P_Value", "PVE_Marker_Percent", "PVE_CIM_Percent",
                 "Additive_Effect", "SE_Additive", "CI_Lower", "CI_Upper", "Dominance_Effect",
                 "Dominance_Ratio", "Dominance_Class", "N_Cofactors", "N_Samples", "Method")

  for (ch in chr_order(sc$Chromosome)) {
    d <- sc[as.character(sc$Chromosome) == ch, , drop = FALSE]
    d <- d[order(d$Position), , drop = FALSE]
    n <- nrow(d)
    lod <- d$LOD; pos <- d$Position
    left <- c(-Inf, lod[-n]); right <- c(lod[-1], -Inf)
    p_ok <- is.na(d$P_Value) | d$P_Value <= p_max
    cand <- which(lod > left & lod >= right & lod >= thr & p_ok)
    if (!length(cand)) next
    cand <- cand[order(-lod[cand])]

    peaks <- integer(0)
    for (j in cand) {
      separate <- TRUE
      for (k in peaks) {
        if (min(lod[min(j, k):max(j, k)]) > lod[j] - peak_drop) { separate <- FALSE; break }
      }
      if (separate) peaks <- c(peaks, j)
    }
    peaks <- sort(peaks); k <- length(peaks)

    dstart <- rep(1L, k); dend <- rep(n, k)
    if (k > 1) for (a in 1:(k - 1)) {
      seg <- peaks[a]:peaks[a + 1]
      v <- seg[which.min(lod[seg])]
      v <- min(max(v, peaks[a] + 1L), peaks[a + 1] - 1L)
      dend[a] <- v - 1L; dstart[a + 1] <- v + 1L
    }

    dx <- if (n > 1) diff(pos) else 1
    w  <- if (n > 1) (c(0, dx) + c(dx, 0)) / 2 else 1
    w[w <= 0] <- if (any(w > 0)) min(w[w > 0]) / 10 else 1

    for (a in seq_len(k)) {
      i <- peaks[a]; target <- lod[i] - lod_drop
      li <- i; while (li > dstart[a] && lod[li - 1] >= target) li <- li - 1L
      ri <- i; while (ri < dend[a]   && lod[ri + 1] >= target) ri <- ri + 1L
      lf <- if (li > 1) li - 1L else 1L
      rf <- if (ri < n) ri + 1L else n
      ci_s <- if (li > 1) interp_cross(pos[li - 1], lod[li - 1], pos[li], lod[li], target) else pos[1]
      ci_e <- if (ri < n) interp_cross(pos[ri + 1], lod[ri + 1], pos[ri], lod[ri], target) else pos[n]

      b_s <- b_e <- NA_real_
      if (bayes) {
        idx  <- dstart[a]:dend[a]
        post <- 10^(lod[idx] - lod[i]) * w[idx]
        post <- post / sum(post)
        o    <- order(post, decreasing = TRUE)
        m    <- which(cumsum(post[o]) >= bayes_prob)[1]
        sel  <- c(idx[o[seq_len(m)]], i)
        b_s  <- min(pos[sel]); b_e <- max(pos[sel])
      }

      row <- as.list(d[i, keep_cols, drop = FALSE])
      row <- lapply(row, function(v) if (is.factor(v)) as.character(v) else v)
      gw  <- if (!is.null(perm_null)) (sum(perm_null >= lod[i]) + 1) / (length(perm_null) + 1) else NA_real_
      out[[length(out) + 1]] <- c(
        list(QTL_ID = NA_character_, Chromosome = ch,
             Left_Marker = d$Marker[lf], Left_Position = pos[lf],
             Peak_Marker = row$Marker, Peak_Position = row$Position,
             Right_Marker = d$Marker[rf], Right_Position = pos[rf],
             Interval_Size_cM = pos[rf] - pos[lf],
             CI_Start_cM = ci_s, CI_End_cM = ci_e, Bayes_Start_cM = b_s, Bayes_End_cM = b_e,
             N_Markers_In_Region = rf - lf + 1L),
        row[c("LOD", "P_Value")], list(Genomewide_P = gw),
        row[c("PVE_Marker_Percent", "PVE_CIM_Percent", "Additive_Effect", "SE_Additive", "CI_Lower",
              "CI_Upper", "Dominance_Effect", "Dominance_Ratio", "Dominance_Class", "N_Cofactors",
              "N_Samples", "Method")],
        list(Threshold_LOD = thr, LOD_Drop = lod_drop))
    }
  }
  res <- rows_to_df(out)
  if (is.null(res)) return(empty_regions())
  res <- res[order(match(res$Chromosome, chr_order(res$Chromosome)), res$Peak_Position), , drop = FALSE]
  res$QTL_ID <- ave(seq_len(nrow(res)), res$Chromosome, FUN = seq_along)
  res$QTL_ID <- paste0("QTL_", res$Chromosome, ".", res$QTL_ID)
  rownames(res) <- NULL
  res[, REGION_COLUMNS, drop = FALSE]
}

# Label every scanned marker with the QTL region it lies in.
tag_scan_regions <- function(scan_all, regions) {
  scan_all$QTL_ID <- NA_character_
  for (i in seq_len(nrow(regions))) {
    hit <- as.character(scan_all$Chromosome) == regions$Chromosome[i] &
           scan_all$Position >= regions$Left_Position[i] & scan_all$Position <= regions$Right_Position[i]
    scan_all$QTL_ID[hit] <- regions$QTL_ID[i]
  }
  scan_all
}

# ============================================================
# MANHATTAN (one panel per chromosome, QTL regions shaded)
# ============================================================
draw_manhattan <- function(scan_all, regions, thr, plot_cfg, title = NULL) {
  pd <- scan_all[!is.na(scan_all$LOD) & !is.na(scan_all$Position) & !is.na(scan_all$Chromosome), , drop = FALSE]
  if (nrow(pd) == 0) { plot.new(); text(0.5, 0.5, "No scan results"); return(invisible()) }
  pd$Chromosome <- as.character(pd$Chromosome)
  chroms <- chr_order(pd$Chromosome)
  n_cols <- ceiling(sqrt(length(chroms))); n_rows <- ceiling(length(chroms) / n_cols)
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1), bg = "white", oma = c(0, 0, if (is.null(title)) 0 else 2, 0))
  lty <- switch(plot_cfg$threshold_style, solid = 1, dash = 2, dot = 3, dashdot = 4, 2)
  shade <- adjustcolor(plot_cfg$sig_color, alpha.f = 0.14)

  for (chr in chroms) {
    d <- pd[pd$Chromosome == chr, , drop = FALSE]
    d <- d[order(d$Position), , drop = FALSE]
    ymax <- max(d$LOD, thr, na.rm = TRUE) * 1.12
    plot(d$Position, d$LOD, type = "n", xlab = paste("Chr", chr, "position (cM)"), ylab = "LOD",
         main = paste("Chr", chr), ylim = c(0, ymax), frame.plot = FALSE)
    rg <- regions[as.character(regions$Chromosome) == chr, , drop = FALSE]
    for (i in seq_len(nrow(rg))) {
      rect(rg$Left_Position[i], 0, rg$Right_Position[i], ymax, col = shade, border = NA)
    }
    lines(d$Position, d$LOD, col = plot_cfg$nonsig_color, lwd = 1)
    points(d$Position, d$LOD, pch = 16, cex = as.numeric(plot_cfg$nonsig_point_size) / 7, col = plot_cfg$nonsig_color)
    abline(h = thr, lty = lty, col = plot_cfg$sig_color, lwd = as.numeric(plot_cfg$threshold_width))
    for (i in seq_len(nrow(rg))) {
      points(rg$Peak_Position[i], rg$LOD[i], pch = 18, cex = as.numeric(plot_cfg$sig_point_size) / 5, col = plot_cfg$sig_color)
      text(rg$Peak_Position[i], rg$LOD[i], rg$Peak_Marker[i], pos = 3, cex = 0.5, col = "#333333", xpd = NA)
    }
  }
  if (!is.null(title)) mtext(title, outer = TRUE, font = 2, cex = 0.9)
}

# ============================================================
# CHROMOSOME IDEOGRAM: each QTL is drawn as its flanking-marker interval
# (bar) with the peak marked (dot), one lane per model.
# ============================================================
make_chromosome_ideogram_reference <- function(qtl_list, all_markers_df,
                                               lod_threshold = 2.0,
                                               model_colors  = MODEL_COLORS,
                                               title_str     = NULL,
                                               band_dark_col  = "#2E6B35",
                                               band_light_col = "#6AAB6E",
                                               border_col     = "#1A4A1F",
                                               qtl_size_min   = 0.7,
                                               qtl_size_max   = 2.0,
                                               x_gap          = 1.8,
                                               lane_step      = 0.26) {
  if (is.null(all_markers_df) || nrow(all_markers_df) == 0) return(invisible(NULL))

  chroms    <- chr_order(all_markers_df$Chromosome)
  n_chr     <- length(chroms)
  bar_w     <- 0.55
  max_pos   <- max(all_markers_df$Position, na.rm = TRUE)
  model_names <- names(qtl_list)

  par(mar = c(5, 4, 4, 8), bg = "white", xpd = TRUE)
  plot(0, 0, type = "n", xlim = c(0.2, n_chr * x_gap + 0.3), ylim = c(-8, max_pos * 1.12),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "", frame.plot = FALSE, bty = "n")
  y_ticks <- pretty(c(0, max_pos), n = 7)
  axis(2, at = y_ticks, labels = y_ticks, cex.axis = 0.72, las = 1, col = "#555555", col.axis = "#333333", tcl = -0.3, lwd = 0.8)
  mtext("Genetic Position (cM)", side = 2, line = 2.5, cex = 0.82, col = "#222222", font = 2)

  for (ci in seq_along(chroms)) {
    chr <- chroms[ci]
    x_center <- (ci - 1) * x_gap + 0.9
    chr_m <- all_markers_df[as.character(all_markers_df$Chromosome) == chr, ]
    if (nrow(chr_m) == 0) next
    chr_len <- max(chr_m$Position, na.rm = TRUE)
    chr_min <- min(chr_m$Position, na.rm = TRUE)

    n_bands <- max(8, round((chr_len - chr_min) / 8))
    band_h  <- (chr_len - chr_min) / n_bands
    for (bi in 1:n_bands) {
      bcol <- if (bi %% 2 == 1) band_dark_col else band_light_col
      rect(x_center - bar_w/2, chr_min + (bi-1)*band_h, x_center + bar_w/2, chr_min + bi*band_h, col = bcol, border = NA)
    }
    rect(x_center - bar_w/2, chr_min, x_center + bar_w/2, chr_len, col = NA, border = border_col, lwd = 1.8)
    symbols(x_center, chr_min, circles = bar_w/2, inches = FALSE, add = TRUE, fg = border_col, bg = band_light_col, lwd = 1.2)
    symbols(x_center, chr_len, circles = bar_w/2, inches = FALSE, add = TRUE, fg = border_col, bg = band_dark_col, lwd = 1.2)

    for (ti in seq_along(model_names)) {
      trait <- model_names[ti]
      q <- qtl_list[[trait]]
      if (is.null(q) || nrow(q) == 0) next
      q <- q[as.character(q$Chromosome) == chr, , drop = FALSE]
      if (nrow(q) == 0) next
      col_m <- model_colors[trait]; if (is.na(col_m)) col_m <- "#CC0000"
      x_off <- x_center + bar_w/2 + 0.12 + (ti - 1) * lane_step
      for (ri in seq_len(nrow(q))) {
        segments(x_off, q$Left_Position[ri], x_off, q$Right_Position[ri], col = col_m, lwd = 4, lend = 1)
        segments(x_center + bar_w/2, q$Peak_Position[ri], x_off, q$Peak_Position[ri], col = col_m, lwd = 1)
        points(x_off, q$Peak_Position[ri], pch = 16, col = col_m,
               cex = pmax(qtl_size_min, pmin(qtl_size_max, 0.6 + q$LOD[ri] / 8)))
        text(x_off + 0.12, q$Peak_Position[ri], round(q$LOD[ri], 1), cex = 0.55, col = col_m, font = 2, adj = 0)
      }
    }
    text(x_center, chr_min - max_pos * 0.04, paste0("Chr", chr), cex = 0.72, font = 2, col = "#1A1A1A", adj = 0.5)
  }

  legend_x <- n_chr * x_gap + 0.5
  legend_y <- max_pos * 0.95
  text(legend_x, legend_y + max_pos * 0.07, "Models", cex = 0.75, font = 2, col = "#1A1A1A", adj = 0)
  for (ti in seq_along(model_names)) {
    tr <- model_names[ti]
    col_m <- model_colors[tr]; if (is.na(col_m)) col_m <- "#CC0000"
    y_leg <- legend_y - (ti - 1) * max_pos * 0.07
    q <- qtl_list[[tr]]
    n_sig <- if (is.null(q)) 0 else nrow(q)
    points(legend_x, y_leg, pch = 16, col = col_m, cex = 1.2)
    text(legend_x + 0.18, y_leg, paste0(tr, " (", n_sig, " QTL)"), cex = 0.65, col = col_m, adj = 0)
  }
  text(legend_x, legend_y - length(model_names) * max_pos * 0.07 - max_pos * 0.03,
       "bar = flanking-marker interval\ndot = QTL peak (LOD)", cex = 0.55, col = "#555555", adj = c(0, 1))

  ttl <- if (!is.null(title_str)) title_str else paste0("QTL regions  |  LOD >= ", signif(lod_threshold, 3))
  title(main = ttl, cex.main = 0.95, font.main = 2, col.main = "#1A1A1A", line = 2)
}
