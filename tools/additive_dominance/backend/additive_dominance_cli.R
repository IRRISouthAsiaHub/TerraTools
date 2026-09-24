#!/usr/bin/env Rscript

# ============================================================
# Genomic Relationship Matrix Builder — headless CLI
#
# Non-interactive port of app__13_.R: computes the additive
# (rrBLUP::A.mat) and/or dominance (sommer::D.mat) genomic
# relationship matrices from a SNP dosage (0/1/2) file. Every
# computational step below (delimiter auto-detection,
# imputation, monomorphic-marker filtering, MAF, -1/0/1
# recoding) is copied verbatim from the original app's
# observeEvent(input$run_btn, ...) block.
#
# Usage:
#   Rscript additive_dominance_cli.R <config.json> <output_dir>
#
# config.json fields:
#   snp_file       - path to .csv/.txt/.tsv (rows = samples,
#                    cols = SNP markers, first column = sample
#                    IDs, values in {0,1,2,NA})
#   impute_method  - "mean" | "median" | "zero"
#   matrix_types   - array subset of ["additive","dominance"]
# ============================================================

suppressMessages({
    library(rrBLUP)
    library(sommer)
    library(ggplot2)
    library(reshape2)
    library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    cat(toJSON(list(
        status = "error",
        message = "Usage: Rscript additive_dominance_cli.R <config.json> <output_dir>"
    ), auto_unbox = TRUE))
    quit(status = 2)
}

config_file <- args[[1]]
output_dir  <- normalizePath(args[[2]], mustWork = FALSE)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

config <- jsonlite::fromJSON(config_file)

fail <- function(message) {
    cat(toJSON(list(status = "error", message = message), auto_unbox = TRUE))
    quit(status = 1)
}


# ============================================================
# HELPERS
# ============================================================

dark_theme <- theme(
    plot.background  = element_rect(fill = "#13161e", color = NA),
    panel.background = element_rect(fill = "#13161e", color = NA),
    panel.grid       = element_blank(),
    axis.text        = element_text(color = "#6b7280", size = 7),
    axis.ticks       = element_line(color = "#2a2f3d"),
    axis.title       = element_text(color = "#6b7280", size = 10),
    plot.title       = element_text(color = "#00e5a0", size = 12,
                                     family = "mono", hjust = 0),
    legend.background = element_rect(fill = "#13161e", color = NA),
    legend.text  = element_text(color = "#6b7280", size = 8),
    legend.title = element_text(color = "#6b7280", size = 8)
)

save_dark_plot <- function(p, name, height_in = 4.2, width_in = 5.4) {
    path <- file.path(output_dir, paste0(name, ".png"))
    ggsave(path, plot = p, height = height_in, width = width_in,
           dpi = 110, units = "in", bg = "#13161e", limitsize = FALSE)
    path
}

heatmap_plot <- function(mat, title, col1 = "#0a2340", col2 = "#00e5a0") {
    n   <- min(nrow(mat), 60)
    sub <- mat[1:n, 1:n]
    df  <- reshape2::melt(sub)
    ggplot(df, aes(Var1, Var2, fill = value)) +
        geom_tile() +
        scale_fill_gradient2(low = col1, mid = "#1c2030",
                              high = col2, midpoint = median(df$value)) +
        labs(title = title, x = NULL, y = NULL, fill = "Value") +
        theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
              axis.ticks  = element_blank()) +
        dark_theme
}

hist_plot <- function(mat, title, col) {
    vals <- mat[upper.tri(mat)]
    df   <- data.frame(v = vals)
    ggplot(df, aes(v)) +
        geom_histogram(bins = 60, fill = col, alpha = .85, color = NA) +
        geom_vline(xintercept = mean(vals), color = "#fff",
                   linetype = "dashed", linewidth = .5) +
        labs(title = title, x = "Relationship coefficient", y = "Count") +
        dark_theme
}

matrix_preview <- function(mat, n_max = 100) {
    n  <- min(nrow(mat), n_max)
    df <- as.data.frame(round(mat[1:n, 1:n], 4))
    df <- cbind(Sample = rownames(df), df)
    rownames(df) <- NULL
    df
}

to_records <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(list())
    jsonlite::fromJSON(
        jsonlite::toJSON(df, na = "null", digits = NA),
        simplifyVector = FALSE
    )
}

# Smart file reader: auto-detects CSV vs TSV, same as the app.
read_snp_file <- function(path) {
    first_line <- readLines(path, n = 1)
    delim <- if (grepl("\t", first_line)) "\t" else ","
    read.table(path, header = TRUE, sep = delim,
               stringsAsFactors = FALSE, check.names = FALSE,
               quote = "", fill = TRUE)
}


# ============================================================
# LOAD + VALIDATE
# ============================================================

snp_file <- config$snp_file

if (is.null(snp_file) || !file.exists(snp_file)) {
    fail(paste("SNP file not found:", snp_file))
}

impute_method <- config$impute_method
if (is.null(impute_method)) impute_method <- "mean"

matrix_types <- config$matrix_types
if (is.null(matrix_types) || length(matrix_types) == 0) {
    matrix_types <- c("additive", "dominance")
}

result <- tryCatch({

    # ---- 1. Read file ----
    geno_data <- read_snp_file(snp_file)

    if (ncol(geno_data) < 2) {
        stop("File must have at least 2 columns (sample ID + markers).")
    }

    sample_ids  <- as.character(geno_data[, 1])
    geno_matrix <- as.matrix(geno_data[, -1, drop = FALSE])
    storage.mode(geno_matrix) <- "numeric"

    if (nrow(geno_matrix) < 2) {
        stop("Need at least 2 samples to build a relationship matrix.")
    }

    pct_miss <- round(100 * mean(is.na(geno_matrix)), 2)

    # ---- 2. Imputation ----
    geno_imputed <- apply(geno_matrix, 2, function(x) {
        if (any(is.na(x))) {
            fill <- switch(impute_method,
                            mean   = mean(x, na.rm = TRUE),
                            median = median(x, na.rm = TRUE),
                            zero   = 0)
            x[is.na(x)] <- fill
        }
        x
    })

    # ---- 3. Remove monomorphic markers ----
    marker_var <- apply(geno_imputed, 2, var, na.rm = TRUE)
    mono_count <- sum(marker_var == 0)
    geno_imputed <- geno_imputed[, marker_var > 0, drop = FALSE]

    if (ncol(geno_imputed) == 0) {
        stop("All markers are monomorphic after filtering. Cannot compute matrices.")
    }

    # ---- 4. MAF (informational) ----
    maf <- apply(geno_imputed, 2, function(x) {
        p <- mean(x, na.rm = TRUE) / 2
        min(p, 1 - p)
    })
    maf_pass <- sum(maf > 0.05)

    # ---- 5. Recode 0/1/2 -> -1/0/1 ----
    geno_coded <- geno_imputed - 1

    # ---- 6. Build matrices ----
    A_mat <- NULL
    D_mat <- NULL

    if ("additive" %in% matrix_types) {
        A_mat <- rrBLUP::A.mat(geno_coded)
        rownames(A_mat) <- colnames(A_mat) <- sample_ids
    }

    if ("dominance" %in% matrix_types) {
        D_mat <- sommer::D.mat(geno_coded)
        rownames(D_mat) <- colnames(D_mat) <- sample_ids
    }

    out <- list(
        status = "success",
        summary = list(
            n_samples = nrow(geno_imputed),
            n_snps = ncol(geno_imputed),
            pct_missing = pct_miss,
            maf_pass = maf_pass,
            mono_removed = mono_count,
            matrix_types = as.list(matrix_types)
        )
    )

    # ---- Additive outputs ----
    if (!is.null(A_mat)) {

        out$heatmap_additive <- save_dark_plot(
            heatmap_plot(A_mat, "ADDITIVE (A)"), "heatmap_additive"
        )

        out$hist_additive <- save_dark_plot(
            hist_plot(A_mat, "ADDITIVE \u2014 off-diagonal distribution", "#00e5a0"),
            "hist_additive"
        )

        out$table_additive <- to_records(matrix_preview(A_mat))

        additive_csv <- file.path(output_dir, "Additive_Matrix.csv")
        write.csv(A_mat, additive_csv, row.names = TRUE)
        out$additive_csv <- additive_csv
    }

    # ---- Dominance outputs ----
    if (!is.null(D_mat)) {

        out$heatmap_dominance <- save_dark_plot(
            heatmap_plot(D_mat, "DOMINANCE (D)", col2 = "#6c63ff"), "heatmap_dominance"
        )

        out$hist_dominance <- save_dark_plot(
            hist_plot(D_mat, "DOMINANCE \u2014 off-diagonal distribution", "#6c63ff"),
            "hist_dominance"
        )

        out$table_dominance <- to_records(matrix_preview(D_mat))

        dominance_csv <- file.path(output_dir, "Dominance_Matrix.csv")
        write.csv(D_mat, dominance_csv, row.names = TRUE)
        out$dominance_csv <- dominance_csv
    }

    out

}, error = function(e) {
    list(status = "error", message = conditionMessage(e))
})


# ============================================================
# WRITE RESULT
# ============================================================

result_json_file <- file.path(output_dir, "additive_dominance_results.json")

json_text <- jsonlite::toJSON(result, auto_unbox = TRUE, na = "null", digits = NA)

writeLines(json_text, result_json_file, useBytes = TRUE)

message("Result JSON written to: ", result_json_file)

cat(json_text, "\n")

quit(status = if (identical(result$status, "error")) 1 else 0)
