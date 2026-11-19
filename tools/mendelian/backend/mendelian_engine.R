# mendelian_engine.R - Mendelian segregation test for explicitly defined populations
# ---------------------------------------------------------------
# Input: Marker, Chromosome, genotypes coded 0/1/2 (see mendelian_cli.R).
#
# Every population type is an explicit configuration (POPULATIONS below): the genotype classes that are
# expected and their ratio. A type that is not in the table is rejected, never tested as something else.
#
#   RIL, DH : 1:1 between the two homozygous classes (0 and 2); heterozygous calls are not part of the test
#   BC1P1   : backcross to parent 1, 1:1 between classes 0 and 1 (class 2 cannot occur)
#   BC1P2   : backcross to parent 2, 1:1 between classes 1 and 2 (class 0 cannot occur)
#   F2      : 1:2:1     F3 : 3:2:3     F4 : 7:2:7   (selfed generations, heterozygosity halves each generation)
#   F1      : rejected - every plant is heterozygous, there is no segregation to test
#
# Tests
#   two classes   : exact two-sided binomial test (p = 0.5), for any number of calls
#   three classes : exact multinomial test (every possible outcome is enumerated; the p-value is the total
#                   probability of the outcomes that are not more likely than the observed one) up to
#                   EXACT_MAX calls; above that the expected counts are large and the chi-square test is used
#   inconclusive  : the marker has no calls, or so few that even the most extreme outcome could not give a
#                   p-value below alpha, so "not distorted" would say nothing. Such markers are neither called
#                   distorted nor passed, and are left out of the multiple-testing correction.
#   p_adj         : p corrected over the conclusive markers (none / Bonferroni / FDR-BH)
#   distorted     : conclusive and p_adj < alpha
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

# -- Population configurations --------------------------------------------------------------
# classes: genotype codes expected; ratio: their expected ratio; het: expected share of heterozygous calls
POPULATIONS <- list(
  RIL   = list(label = "RIL",                  classes = c(0, 2),    ratio = c(1, 1),    het = 0),
  DH    = list(label = "Doubled haploid (DH)", classes = c(0, 2),    ratio = c(1, 1),    het = 0),
  BC1P1 = list(label = "BC1 to parent 1",      classes = c(0, 1),    ratio = c(1, 1),    het = 0.5),
  BC1P2 = list(label = "BC1 to parent 2",      classes = c(1, 2),    ratio = c(1, 1),    het = 0.5),
  F2    = list(label = "F2",                   classes = c(0, 1, 2), ratio = c(1, 2, 1), het = 0.5),
  F3    = list(label = "F3",                   classes = c(0, 1, 2), ratio = c(3, 2, 3), het = 0.25),
  F4    = list(label = "F4",                   classes = c(0, 1, 2), ratio = c(7, 2, 7), het = 0.125)
)

population_config <- function(pop) {
  key <- toupper(trimws(as.character(pop %||% "")))
  if (length(key) != 1 || !nzchar(key) || !key %in% names(POPULATIONS)) {
    why <- if (identical(key, "F1")) " F1 plants are all heterozygous, so there is no segregation to test." else
           if (identical(key, "BC1")) " Say which parent the backcross is to: BC1P1 or BC1P2." else ""
    stop(sprintf("Unsupported population type '%s'. Supported types: %s.%s",
                 paste(as.character(pop %||% ""), collapse = ""), paste(names(POPULATIONS), collapse = ", "), why),
         call. = FALSE)
  }
  cfg <- POPULATIONS[[key]]
  cfg$key <- key
  cfg
}

EXACT_MAX <- 200   # largest number of calls for which the multinomial test is enumerated exactly

# All outcomes (a, b, c) of n calls in three classes with their probabilities, sorted; cached per (ratio, n)
.mn_cache <- new.env(parent = emptyenv())
multinomial_table <- function(n, prob) {
  key <- paste(n, paste(prob, collapse = ","))
  if (!is.null(.mn_cache[[key]])) return(.mn_cache[[key]])
  a <- rep(0:n, times = n + 1 - (0:n))
  b <- unlist(lapply(0:n, function(i) 0:(n - i)))
  cc <- n - a - b
  lp <- lgamma(n + 1) - lgamma(a + 1) - lgamma(b + 1) - lgamma(cc + 1) +
        a * log(prob[1]) + b * log(prob[2]) + cc * log(prob[3])
  pr <- sort(exp(lp))
  tab <- list(p = pr, cum = cumsum(pr))
  assign(key, tab, envir = .mn_cache)
  tab
}

# p-value of an observed outcome: the total probability of the outcomes that are not more likely than it
# p_min is the p-value of the most extreme outcome, the smallest p this n can ever give
exact_multinomial <- function(obs, prob) {
  n <- sum(obs)
  tab <- multinomial_table(n, prob)
  lp <- lgamma(n + 1) - sum(lgamma(obs + 1)) + sum(obs * log(prob))
  tol <- 1 + 1e-7
  idx <- findInterval(exp(lp) * tol, tab$p)
  p_min <- tab$cum[findInterval(tab$p[1] * tol, tab$p)]
  list(p = min(1, tab$cum[max(idx, 1)]), p_min = min(1, p_min))
}

# Test of one marker (a vector of 0/1/2/NA calls) in the population `cfg`
test_marker <- function(vec, cfg, pthr) {
  v <- vec[!is.na(vec)]
  counts <- c(sum(v == 0), sum(v == 1), sum(v == 2))
  cls <- counts[cfg$classes + 1]
  n <- sum(cls)
  out <- list(n = n, n_unexpected = sum(counts) - n, c0 = counts[1], c1 = counts[2], c2 = counts[3],
              chisq = NA_real_, df = NA_real_, test = "none", p = NA_real_, p_min = NA_real_,
              inconclusive = TRUE, het_rate = if (sum(counts) > 0) counts[2] / sum(counts) else NA_real_)
  if (n == 0) return(out)
  prob <- cfg$ratio / sum(cfg$ratio)
  exp_n <- prob * n
  out$chisq <- sum((cls - exp_n)^2 / exp_n)
  out$df <- length(cls) - 1
  if (length(cls) == 2) {
    out$test <- "exact binomial"
    out$p <- stats::binom.test(cls[1], n, prob[1])$p.value
    out$p_min <- min(1, 2 * 0.5^n)
  } else if (n <= EXACT_MAX) {
    out$test <- "exact multinomial"
    r <- exact_multinomial(cls, prob)
    out$p <- r$p; out$p_min <- r$p_min
  } else {
    out$test <- "chi-square (large n)"
    out$p <- stats::pchisq(out$chisq, df = out$df, lower.tail = FALSE)
    out$p_min <- 0
  }
  # unresolved: even the most extreme outcome could not reach alpha, so "not distorted" carries no information
  out$inconclusive <- out$p_min >= pthr
  out
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
  cfg <- population_config(pop)          # rejects a population type that is not defined
  pop <- cfg$key
  mtc_method <- if (mtc_method %in% names(MTC_LABELS)) mtc_method else "fdr"
  mtc_name <- MTC_LABELS[[mtc_method]]

  df <- read_and_validate(input_file, missing_codes)
  geno_mat <- as.matrix(df[, -(1:2), drop = FALSE])

  res_df <- dplyr::bind_rows(lapply(seq_len(nrow(geno_mat)), function(i) test_marker(geno_mat[i, ], cfg, pthr))) |>
    dplyr::mutate(row = dplyr::row_number())

  # Genotype calls stay in the full table (after the statistics) so the file can feed the next step
  geno_cols <- df[, -(1:2), drop = FALSE]
  out <- dplyr::bind_cols(df[, 1:2], res_df, geno_cols) |> dplyr::mutate(`-log10p` = -log10(p))

  # Correction over the conclusive markers (before the unplaced-marker filter, so it reflects the true number of tests).
  # Inconclusive markers cannot reach significance at this alpha: they are left out of the correction and never called.
  conclusive <- !is.na(out$p) & !out$inconclusive
  p_adj <- rep(NA_real_, nrow(out))
  p_adj[conclusive] <- if (mtc_method == "none") out$p[conclusive] else
    p.adjust(out$p[conclusive], method = if (mtc_method == "bonferroni") "bonferroni" else "BH")
  out <- out |>
    dplyr::mutate(p_adj = p_adj,
                  mtc_method_used = mtc_method,
                  distorted = !is.na(p_adj) & p_adj < pthr,
                  call = dplyr::case_when(is.na(p) ~ "no data", inconclusive ~ "inconclusive",
                                          distorted ~ "distorted", TRUE ~ "not distorted")) |>
    dplyr::select(-p_min)

  if (isTRUE(drop_unplaced)) out <- out |> dplyr::filter(!is.na(Chromosome) & Chromosome != "" & !is.na(p))

  sample_cols <- colnames(geno_cols)
  stat_cols <- setdiff(colnames(out), sample_cols)

  # ── Files ────────────────────────────────────────────────────────────────
  pass <- out |> dplyr::filter(call == "not distorted")          # markers that pass the quality check
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
  n_tests <- sum(!dd$inconclusive)                 # markers that were really tested (conclusive)
  manhattan_png <- file.path(outdir, "manhattan.png")
  rateplot_png  <- file.path(outdir, "rateplot.png")
  hetplot_png   <- file.path(outdir, "heterozygosity.png")

  if (n_tests > 0) {
    # Manhattan plot: one colour per chromosome, in chromosome order along the marker index
    lv_chr <- chr_levels(dd$Chromosome)
    dd$ChrF <- factor(as.character(dd$Chromosome), levels = lv_chr)
    # The dashed line is the cut-off of the chosen correction, expressed as a raw p-value:
    #   Bonferroni: alpha / number of tests;  none: alpha;
    #   FDR (Benjamini-Hochberg): the largest raw p among the markers called distorted (nothing is drawn when none is)
    cut_p <- switch(mtc_method,
      bonferroni = pthr / n_tests,
      none = pthr,
      fdr = if (any(dd$distorted)) max(dd$p[dd$distorted], na.rm = TRUE) else NA_real_)
    pm <- ggplot(dd, aes(x = row, y = `-log10p`, color = ChrF)) +
      geom_point(alpha = 0.9, size = 1.6)
    if (is.finite(cut_p)) pm <- pm + geom_hline(yintercept = -log10(cut_p), linetype = 2, colour = "grey20")
    pm <- pm +
      guides(color = guide_legend(title = "Chr", nrow = min(3, length(lv_chr)), byrow = FALSE, override.aes = list(size = 2.6, alpha = 1))) +
      labs(x = "Marker index", y = expression(-log[10](p)),
           caption = paste0("Dashed line: ", mtc_name, if (is.finite(cut_p)) paste0(" cut-off (raw p = ", signif(cut_p, 2), ")")
                            else if (mtc_method == "fdr") " - no marker is distorted, so there is no cut-off to draw" else "",
                            if (any(dd$inconclusive)) paste0("\n", sum(dd$inconclusive), " inconclusive marker(s) (too few calls to reach alpha) are drawn but never called.") else "")) +
      theme_mend()
    save_png(manhattan_png, pm)

    # Heterozygosity plot: % heterozygous calls of every marker, coloured by chromosome like the segregation plot.
    # The dashed line is the heterozygosity expected in the population (F2 50 %, F3 25 %, BC1 50 %, ...); RILs and DHs are
    # homozygous, so heterozygous calls there are residual or errors.
    hd <- dd[!is.na(dd$het_rate), , drop = FALSE]
    hd$het_pct <- 100 * hd$het_rate
    ph <- ggplot(hd, aes(x = row, y = het_pct, color = ChrF)) +
      geom_point(alpha = 0.9, size = 1.6)
    if (cfg$het > 0) ph <- ph + geom_hline(yintercept = 100 * cfg$het, linetype = 2, colour = "grey20")
    ph <- ph +
      # a floor on the axis keeps a file without heterozygous calls (an inbred RIL set) from being stretched to nothing
      scale_y_continuous(labels = function(v) paste0(v, "%")) +
      coord_cartesian(ylim = c(0, max(if (cfg$het > 0) 100 * cfg$het * 1.2 else 5, max(hd$het_pct, na.rm = TRUE) * 1.1))) +
      guides(color = guide_legend(title = "Chr", nrow = min(3, length(lv_chr)), byrow = FALSE, override.aes = list(size = 2.6, alpha = 1))) +
      labs(x = "Marker index", y = "Heterozygous calls (%)",
           caption = if (cfg$het > 0) sprintf("Dashed line: %s %% heterozygous, as expected in %s.", format(100 * cfg$het), cfg$label) else
             paste0(cfg$label, " plants are expected to be homozygous: markers with many heterozygous calls deserve a check.")) +
      theme_mend()
    save_png(hetplot_png, ph)

    lv <- chr_levels(dd$Chromosome)
    rate <- dd |> dplyr::filter(!inconclusive) |> dplyr::group_by(Chromosome) |>
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
    pop_label = cfg$label,
    n_tested = sum(!is.na(out$p) & !out$inconclusive),
    n_inconclusive = sum(out$call == "inconclusive"),
    n_no_data = sum(out$call == "no data"),
    n_unexpected_calls = sum(out$n_unexpected, na.rm = TRUE),
    n_distorted = sum(out$distorted, na.rm = TRUE),
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
    hetplot_png = if (file.exists(hetplot_png)) hetplot_png else NULL,
    results_json = results_json
  )
}
