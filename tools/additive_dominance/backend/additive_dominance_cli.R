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
#   train_file,    - INSTEAD of snp_file: the Training and Testing files
#   test_file        written by the Clustering tool. Markers are centred on
#                    the TRAINING lines' means and the additive matrix is
#                    X Xt / number of markers, cut into the blocks
#                    training x training and testing x training (the testing
#                    lines do not influence the training block). Imputation
#                    and the monomorphic filter also use the training lines.
#                    Columns named phenotypes / Cluster are not markers.
#   impute_method  - "mean" | "median" | "zero"
#   matrix_types   - array subset of ["additive","dominance"]
#
# -1 is read as a missing call (the coding of the Clustering tool).
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

train_file <- config$train_file
test_file  <- config$test_file
split_mode <- !is.null(train_file) && nzchar(train_file)

if (split_mode) {
    if (!file.exists(train_file)) fail(paste("Training file not found:", train_file))
    if (is.null(test_file) || !nzchar(test_file) || !file.exists(test_file)) {
        fail(paste("Testing file not found:", test_file))
    }
} else {
    snp_file <- config$snp_file
    if (is.null(snp_file) || !file.exists(snp_file)) {
        fail(paste("SNP file not found:", snp_file))
    }
}

impute_method <- config$impute_method
if (is.null(impute_method)) impute_method <- "mean"

matrix_types <- config$matrix_types
if (is.null(matrix_types) || length(matrix_types) == 0) {
    matrix_types <- c("additive", "dominance")
}

result <- tryCatch({

    # ---- 1. Read file(s) ----
    is_non_marker <- function(nm) grepl("^(phenotypes?|cluster)$", nm, ignore.case = TRUE)

    read_markers <- function(path) {
        d <- read_snp_file(path)
        if (ncol(d) < 2) {
            stop("File must have at least 2 columns (sample ID + markers).")
        }
        g <- d[, -1, drop = FALSE]
        list(ids = as.character(d[, 1]), geno = g[, !is_non_marker(colnames(g)), drop = FALSE])
    }

    split_info <- NULL
    if (split_mode) {
        tr <- read_markers(train_file)
        te <- read_markers(test_file)
        shared_ids <- intersect(tr$ids, te$ids)
        if (length(shared_ids) > 0) {
            stop(paste0(length(shared_ids), " line(s) are in both the training and the testing file (e.g. ",
                        paste(head(shared_ids, 3), collapse = ", "), "). A line must be in only one of them."))
        }
        common <- intersect(colnames(tr$geno), colnames(te$geno))
        if (length(common) == 0) {
            stop("The training and testing files have no marker columns in common.")
        }
        geno_matrix <- as.matrix(rbind(tr$geno[, common, drop = FALSE], te$geno[, common, drop = FALSE]))
        sample_ids  <- c(tr$ids, te$ids)
        split_info  <- list(n_train = length(tr$ids), n_test = length(te$ids),
                            n_markers_train = ncol(tr$geno), n_markers_test = ncol(te$geno),
                            n_markers_common = length(common))
    } else {
        d <- read_markers(snp_file)
        geno_matrix <- as.matrix(d$geno)
        sample_ids  <- d$ids
    }
    storage.mode(geno_matrix) <- "numeric"

    # -1 is a missing call (as in the Clustering tool), not a genotype
    n_minus1 <- sum(geno_matrix == -1, na.rm = TRUE)
    geno_matrix[geno_matrix == -1] <- NA

    if (nrow(geno_matrix) < 2) {
        stop("Need at least 2 samples to build a relationship matrix.")
    }

    pct_miss <- round(100 * mean(is.na(geno_matrix)), 2)

    # Rows that define the marker statistics: the training lines in split mode, every line otherwise
    ref_rows <- if (split_mode) seq_len(split_info$n_train) else seq_len(nrow(geno_matrix))

    # ---- 2. Imputation ----
    geno_imputed <- apply(geno_matrix, 2, function(x) {
        if (any(is.na(x))) {
            xr <- x[ref_rows]
            if (all(is.na(xr))) xr <- x
            fill <- switch(impute_method,
                            mean   = mean(xr, na.rm = TRUE),
                            median = median(xr, na.rm = TRUE),
                            zero   = 0)
            x[is.na(x)] <- fill
        }
        x
    })

    # ---- 3. Remove monomorphic markers ----
    marker_var <- apply(geno_imputed[ref_rows, , drop = FALSE], 2, var, na.rm = TRUE)
    marker_var[is.na(marker_var)] <- 0
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
        if (split_mode) {
            # centre every marker on the training lines' mean, then X Xt / number of markers
            Xc <- sweep(geno_imputed, 2, colMeans(geno_imputed[ref_rows, , drop = FALSE]))
            A_mat <- tcrossprod(Xc) / ncol(Xc)
        } else {
            A_mat <- rrBLUP::A.mat(geno_coded)
        }
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
            minus1_as_missing = n_minus1,
            matrix_types = as.list(matrix_types),
            mode = if (split_mode) "split" else "single",
            additive_method = if (split_mode) "centred on training means, X Xt / markers" else "rrBLUP::A.mat"
        )
    )
    if (split_mode) {
        out$summary$split <- split_info
        n_tr    <- split_info$n_train
        tr_ids  <- sample_ids[seq_len(n_tr)]
        te_ids  <- sample_ids[(n_tr + 1):length(sample_ids)]
        # old block files must not survive a run that did not write them
        unlink(file.path(output_dir, c("training_additive.csv", "test_additive.csv",
                                       "training_dominance.csv", "test_dominance.csv")))
    }

    # Rows x columns of M as a CSV whose first column is `sample_id` and whose other columns are named by `cols`
    write_block <- function(M, rows, cols, path) {
        B  <- signif(M[rows, cols, drop = FALSE], 7)
        df <- data.frame(sample_id = rows, B, check.names = FALSE, stringsAsFactors = FALSE)
        write.csv(df, path, row.names = FALSE)
        path
    }

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

        if (split_mode) {
            # training x training  /  testing x training
            out$train_additive_csv <- write_block(A_mat, tr_ids, tr_ids, file.path(output_dir, "training_additive.csv"))
            out$test_additive_csv  <- write_block(A_mat, te_ids, tr_ids, file.path(output_dir, "test_additive.csv"))
        } else {
            additive_csv <- file.path(output_dir, "Additive_Matrix.csv")
            write.csv(A_mat, additive_csv, row.names = TRUE)
            out$additive_csv <- additive_csv
        }
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

        if (split_mode) {
            out$train_dominance_csv <- write_block(D_mat, tr_ids, tr_ids, file.path(output_dir, "training_dominance.csv"))
            out$test_dominance_csv  <- write_block(D_mat, te_ids, tr_ids, file.path(output_dir, "test_dominance.csv"))
        } else {
            dominance_csv <- file.path(output_dir, "Dominance_Matrix.csv")
            write.csv(D_mat, dominance_csv, row.names = TRUE)
            out$dominance_csv <- dominance_csv
        }
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
