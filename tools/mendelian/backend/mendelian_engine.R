# mendelian_engine.R — Mendelian segregation test (RIL 1:1 and F2 1:2:1)
# ---------------------------------------------------------------
# Same statistics as the Mendelian Segregation Checker (app.R), restructured from
# Shiny reactives into a plain batch function so the desktop app can call it
# (see mendelian_cli.R). Input: Marker, Chromosome, genotypes coded 0/1/2.
#
#   RIL : chi-square against 1:1 on the 0 and 2 classes (heterozygotes ignored), df = 1
#   F2  : chi-square against 1:2:1, df = 2
#   low_power  : an expected class count is below 5 (the chi-square is unreliable)
#   p_adj      : p corrected for the number of markers tested (none / Bonferroni / FDR-BH)
#   distorted  : p_adj < alpha
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(jsonlite)
  library(openxlsx)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

MTC_LABELS <- c(none = "uncorrected", bonferroni = "Bonferroni", fdr = "FDR-BH")

# Chi-square test of one marker (a vector of 0/1/2/NA calls)
test_marker <- function(vec, pop) {
  v <- vec[!is.na(vec)]
  empty <- list(n = 0, c0 = NA, c1 = NA, c2 = NA, chisq = NA, df = NA, p = NA, het_rate = NA, low_power = NA)
  if (length(v) == 0) return(empty)

  c0 <- sum(v == 0); c1 <- sum(v == 1); c2 <- sum(v == 2)

  if (pop == "F2") {
    n <- c0 + c1 + c2
    if (n == 0) return(modifyList(empty, list(n = 0, c0 = c0, c1 = c1, c2 = c2)))
    exp <- c(0.25, 0.50, 0.25) * n
    obs <- c(c0, c1, c2)
    low_power <- any(exp < 5)                       # validity of the chi-square rests on the EXPECTED counts
    keep <- exp > 0
    chisq <- sum((obs[keep] - exp[keep])^2 / pmax(exp[keep], .Machine$double.eps))
    list(n = n, c0 = c0, c1 = c1, c2 = c2, chisq = chisq, df = 2,
         p = pchisq(chisq, df = 2, lower.tail = FALSE), het_rate = c1 / n, low_power = low_power)
  } else {                                          # RIL: only the 0 and 2 classes are tested
    n <- c0 + c2
    if (n == 0) return(modifyList(empty, list(n = 0, c0 = c0, c1 = c1, c2 = c2)))
    exp <- c(0.5, 0.5) * n
    obs <- c(c0, c2)
    low_power <- any(exp < 5)
    chisq <- sum((obs - exp)^2 / pmax(exp, .Machine$double.eps))
    list(n = n, c0 = c0, c1 = c1, c2 = c2, chisq = chisq, df = 1,
         p = pchisq(chisq, df = 1, lower.tail = FALSE),
         het_rate = ifelse((c0 + c1 + c2) > 0, c1 / (c0 + c1 + c2), NA_real_), low_power = low_power)
  }
}

# Reads the CSV; codes in `missing_codes` become NA, genotypes become numeric
read_and_validate <- function(input_file, missing_codes) {
  df <- readr::read_csv(input_file, show_col_types = FALSE, progress = FALSE)
  if (ncol(df) < 3) stop("The file needs at least 3 columns: Marker, Chromosome, then the genotypes.")
  names(df)[1:2] <- c("Marker", "Chromosome")
  geno <- df[, -(1:2), drop = FALSE]
  geno[] <- lapply(geno, function(x) {
    x <- as.character(x)
    x[x %in% missing_codes] <- NA_character_
    suppressWarnings(as.numeric(x))
  })
  dplyr::bind_cols(df[1:2], geno)
}

# Chromosome labels in natural order (1, 2, ... 10, then names such as X or Un)
chr_levels <- function(x) {
  u <- unique(as.character(x))
  num <- suppressWarnings(as.numeric(u))
  c(u[!is.na(num)][order(num[!is.na(num)])], sort(u[is.na(num)]))
}

theme_mend <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), plot.background = element_rect(fill = "white", colour = NA),
          legend.position = "bottom")
}

save_png <- function(path, p, width = 7.5, height = 4.6) {
  ggsave(path, p, width = width, height = height, dpi = 200, bg = "white",
         device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png")
}

run_mendelian <- function(input_file, pop, pthr, missing_codes, drop_unplaced, outdir, mtc_method = "fdr") {

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  pop <- if (identical(toupper(pop), "F2")) "F2" else "RIL"
  mtc_method <- if (mtc_method %in% names(MTC_LABELS)) mtc_method else "fdr"
  mtc_name <- MTC_LABELS[[mtc_method]]

  df <- read_and_validate(input_file, missing_codes)
  geno_mat <- as.matrix(df[, -(1:2), drop = FALSE])

  res_df <- dplyr::bind_rows(apply(geno_mat, 1, test_marker, pop = pop)) |> dplyr::mutate(row = dplyr::row_number())

  # Genotype calls stay in the full table (after the statistics) so the file can feed the next step
  geno_cols <- df[, -(1:2), drop = FALSE]
  out <- dplyr::bind_cols(df[, 1:2], res_df, geno_cols) |> dplyr::mutate(`-log10p` = -log10(p))

  # Correction over every tested marker (before the unplaced-marker filter, so it reflects the true number of tests)
  out <- out |>
    dplyr::mutate(p_adj = if (mtc_method == "none") p else p.adjust(p, method = if (mtc_method == "bonferroni") "bonferroni" else "BH"),
                  mtc_method_used = mtc_method,
                  distorted = !is.na(p_adj) & p_adj < pthr)

  if (isTRUE(drop_unplaced)) out <- out |> dplyr::filter(!is.na(Chromosome) & Chromosome != "" & !is.na(p))

  sample_cols <- colnames(geno_cols)
  stat_cols <- setdiff(colnames(out), sample_cols)

  # ── Files ────────────────────────────────────────────────────────────────
  pass <- out |> dplyr::filter(!is.na(p_adj) & !distorted)        # markers that pass the quality check
  full_csv <- file.path(outdir, "results_full.csv")
  filt_csv <- file.path(outdir, "results_filtered.csv")
  lm_xlsx  <- file.path(outdir, "markers_LMready.xlsx")
  readr::write_csv(out, full_csv)
  readr::write_csv(pass, filt_csv)
  # Linkage Mapping reads column 1 = Marker, column 2 = Chromosome and everything after as genotype calls,
  # so the LM-ready file holds exactly those columns (0/1/2/NA coding is detected there).
  openxlsx::write.xlsx(pass[, c("Marker", "Chromosome", sample_cols), drop = FALSE], lm_xlsx)

  # ── Figures ──────────────────────────────────────────────────────────────
  dd <- out |> dplyr::filter(!is.na(p))
  n_tests <- nrow(dd)
  manhattan_png <- file.path(outdir, "manhattan.png")
  rateplot_png  <- file.path(outdir, "rateplot.png")

  if (n_tests > 0) {
    pm <- ggplot(dd, aes(x = row, y = `-log10p`, color = distorted)) +
      geom_point(alpha = 0.9, size = 1.6) +
      scale_color_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"), drop = FALSE,
                         name = paste0("Distorted (", mtc_name, ")"))
    # Bonferroni and "none" have one raw-p cut-off to draw; FDR depends on the rank, so the colour carries the call
    if (mtc_method == "bonferroni") pm <- pm + geom_hline(yintercept = -log10(pthr / n_tests), linetype = 2)
    if (mtc_method == "none")       pm <- pm + geom_hline(yintercept = -log10(pthr), linetype = 2)
    pm <- pm + labs(x = "Marker index", y = expression(-log[10](p)),
                    caption = if (mtc_method == "fdr") "FDR (Benjamini-Hochberg) has no single reference line; the point colour shows the corrected call." else NULL) +
      theme_mend()
    save_png(manhattan_png, pm)

    lv <- chr_levels(dd$Chromosome)
    rate <- dd |> dplyr::group_by(Chromosome) |>
      dplyr::summarise(n = dplyr::n(), distorted = sum(distorted), rate = distorted / n, .groups = "drop")
    rate$Chromosome <- factor(as.character(rate$Chromosome), levels = lv)
    pr <- ggplot(rate, aes(x = Chromosome, y = rate, fill = Chromosome)) +
      geom_col() +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(x = "Chromosome", y = paste0("% distorted (", mtc_name, ")")) +
      theme_mend() + theme(legend.position = "none")
    save_png(rateplot_png, pr)
  }

  # ── Table for the window: the statistics only (the genotype columns stay in the downloads) ──
  results_json <- file.path(outdir, "results.json")
  jsonlite::write_json(out[, stat_cols, drop = FALSE], results_json, na = "null", auto_unbox = TRUE, digits = 6)

  list(
    status = "success",
    n_markers_total = nrow(df),
    n_markers_kept = nrow(out),
    n_tested = sum(!is.na(out$p)),
    n_distorted = sum(out$distorted, na.rm = TRUE),
    n_low_power = sum(out$low_power, na.rm = TRUE),
    n_pass = nrow(pass),
    n_filtered_nondistorted = nrow(pass),
    pop = pop,
    pthr = pthr,
    mtc_method = mtc_method,
    mtc_label = mtc_name,
    results_full_csv = full_csv,
    results_filtered_csv = filt_csv,
    lm_ready_xlsx = lm_xlsx,
    manhattan_png = if (file.exists(manhattan_png)) manhattan_png else NULL,
    rateplot_png = if (file.exists(rateplot_png)) rateplot_png else NULL,
    results_json = results_json
  )
}
