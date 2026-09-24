# mendelian_engine.R — Core Mendelian segregation logic (RIL & F2)
# ---------------------------------------------------------------
# This is the same statistical logic as the Hugging Face app.R Shiny app,
# restructured from Shiny reactives into a plain batch function so it can
# be called from the command line (see mendelian_cli.R) and driven by the
# desktop app's Python backend. No math has been changed.
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(purrr)
  library(scales)
  library(jsonlite)
})

# Same per-marker chi-square test as app.R's test_marker()
test_marker <- function(vec, pop) {
  v <- vec
  v <- v[!is.na(v)]
  if (length(v) < 5) {
    return(list(n = length(v), c0 = NA, c1 = NA, c2 = NA, chisq = NA, df = NA, p = NA, het_rate = NA))
  }

  c0 <- sum(v == 0)
  c1 <- sum(v == 1)
  c2 <- sum(v == 2)

  if (pop == "F2") {
    n <- c0 + c1 + c2
    if (n == 0) return(list(n = 0, c0 = c0, c1 = c1, c2 = c2, chisq = NA, df = NA, p = NA, het_rate = NA))
    exp <- c(0.25, 0.50, 0.25) * n
    obs <- c(c0, c1, c2)
    keep <- exp > 0
    chisq <- sum((obs[keep] - exp[keep])^2 / pmax(exp[keep], .Machine$double.eps))
    pval <- pchisq(chisq, df = 2, lower.tail = FALSE)
    het_rate <- c1 / n
    list(n = n, c0 = c0, c1 = c1, c2 = c2, chisq = chisq, df = 2, p = pval, het_rate = het_rate)

  } else { # RIL
    n <- c0 + c2
    if (n == 0) return(list(n = 0, c0 = c0, c1 = c1, c2 = c2, chisq = NA, df = NA, p = NA, het_rate = NA))
    exp <- c(0.5, 0.5) * n
    obs <- c(c0, c2)
    chisq <- sum((obs - exp)^2 / pmax(exp, .Machine$double.eps))
    pval <- pchisq(chisq, df = 1, lower.tail = FALSE)
    het_rate <- ifelse((c0 + c1 + c2) > 0, c1 / (c0 + c1 + c2), NA_real_)
    list(n = n, c0 = c0, c1 = c1, c2 = c2, chisq = chisq, df = 1, p = pval, het_rate = het_rate)
  }
}

# Reads + validates the CSV exactly as app.R's raw_data()/valid_data() reactives did
read_and_validate <- function(input_file, missing_codes) {
  df <- readr::read_csv(input_file, show_col_types = FALSE, progress = FALSE)

  if (ncol(df) < 3) {
    stop("CSV must have >=3 columns: Marker, Chromosome, and genotypes from 3rd col.")
  }
  names(df)[1:2] <- c("Marker", "Chromosome")

  geno <- df[, -(1:2)]
  geno[] <- lapply(geno, function(x) {
    x <- as.character(x)
    x[x %in% missing_codes] <- NA_character_
    suppressWarnings(as.numeric(x))
  })

  dplyr::bind_cols(df[1:2], geno)
}

# Runs the full segregation test across all markers, exactly as app.R's
# results eventReactive() did, and writes outputs to disk instead of
# rendering them in a Shiny UI.
#
# Returns a list (converted to JSON by the caller in mendelian_cli.R).
run_mendelian <- function(input_file, pop, pthr, missing_codes, drop_unplaced, outdir) {

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  df <- read_and_validate(input_file, missing_codes)
  geno_mat <- as.matrix(df[, -(1:2)])

  res_list <- apply(geno_mat, 1, test_marker, pop = pop)
  res_df <- dplyr::bind_rows(res_list) |> dplyr::mutate(row = dplyr::row_number())

  out <- dplyr::bind_cols(df[, 1:2], res_df) |>
    dplyr::mutate(`-log10p` = -log10(p))

  if (isTRUE(drop_unplaced)) {
    out <- out |> dplyr::filter(!is.na(Chromosome) & Chromosome != "" & !is.na(p))
  }

  # --- Downloads: same two CSVs app.R offered ---------------------------
  filtered <- out |> dplyr::filter(!is.na(p) & p >= pthr)

  results_full_path <- file.path(outdir, "results_full.csv")
  results_filtered_path <- file.path(outdir, "results_filtered.csv")
  readr::write_csv(out, results_full_path)
  readr::write_csv(filtered, results_filtered_path)

  # --- Plots: same two plots app.R rendered ------------------------------
  manhattan_path <- file.path(outdir, "manhattan.png")
  rateplot_path <- file.path(outdir, "rateplot.png")

  dd <- out |> dplyr::filter(!is.na(p))
  dd <- dd |> dplyr::group_by(Chromosome) |> dplyr::mutate(idx = dplyr::row_number()) |> dplyr::ungroup()

  p_manhattan <- ggplot2::ggplot(dd, ggplot2::aes(x = row, y = `-log10p`, group = Chromosome, color = as.factor(Chromosome))) +
    ggplot2::geom_point(alpha = 0.9, size = 1.6) +
    ggplot2::geom_hline(yintercept = -log10(pthr), linetype = 2) +
    ggplot2::labs(x = "Marker index", y = expression(-log[10](p)), color = "Chr") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
  ggplot2::ggsave(manhattan_path, p_manhattan, width = 10, height = 5.2, dpi = 130)

  dd2 <- out |> dplyr::filter(!is.na(p)) |> dplyr::mutate(distorted = p < pthr)
  rate <- dd2 |> dplyr::group_by(Chromosome) |>
    dplyr::summarise(n = dplyr::n(), distorted = sum(distorted), rate = distorted / n, .groups = "drop") |>
    dplyr::arrange(as.numeric(as.character(Chromosome)))

  p_rate <- ggplot2::ggplot(rate, ggplot2::aes(x = reorder(Chromosome, as.numeric(as.character(Chromosome))), y = rate, fill = Chromosome)) +
    ggplot2::geom_col() +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Chromosome", y = "% distorted (p < threshold)") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none")
  ggplot2::ggsave(rateplot_path, p_rate, width = 10, height = 5.2, dpi = 130)

  # --- Table data for the frontend, and summary stats --------------------
  results_json_path <- file.path(outdir, "results.json")
  jsonlite::write_json(out, results_json_path, na = "null", auto_unbox = TRUE, digits = 6)

  list(
    status = "success",
    n_markers_total = nrow(df),
    n_markers_kept = nrow(out),
    n_distorted = sum(out$p < pthr, na.rm = TRUE),
    n_filtered_nondistorted = nrow(filtered),
    pop = pop,
    pthr = pthr,
    results_full_csv = results_full_path,
    results_filtered_csv = results_filtered_path,
    manhattan_png = manhattan_path,
    rateplot_png = rateplot_path,
    results_json = results_json_path
  )
}
