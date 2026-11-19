# =============================================================================
# 6. TABLES: readable column names, sensible rounding, one row per accession
# =============================================================================
# The analysis functions above return tidy "long" tables that are easy to plot. The functions below
# turn them into the tables a person wants to read: plain column names, units, rounding, wide layout.
rename_cols <- function(df, map) {
  nm <- names(df); i <- match(nm, names(map)); nm[!is.na(i)] <- unname(map[i[!is.na(i)]]); names(df) <- nm
  as.data.frame(df, check.names = FALSE)
}
H2 <- "H²"; GXE <- "G×E"

# Layout of each environment; the standard errors of differences (augmented design) are added as columns
tidy_design_check <- function(df, sed = NULL) {
  out <- rename_cols(df, c(
    Genotypes = "Accessions", Reps = "Replications", Block_Size = "Plots per block",
    Min_Replication = "Min. plots per accession", Max_Replication = "Max. plots per accession",
    Repeated_Genotypes = "Accessions repeated (checks)", Single_Plot_Genotypes = "Accessions with one plot"))
  if (!is.null(sed) && nrow(sed)) {
    for (cmp in unique(sed$Comparison)) {
      s <- sed[sed$Comparison == cmp, ]
      out[[paste0("SED: ", tolower(substr(cmp, 1, 1)), substring(cmp, 2))]] <- round(s$SED_mean[match(out$Environment, s$Environment)], 3)
    }
  }
  out
}

# Descriptive statistics of every trait: all environments together ("All") and in each environment
describe_traits <- function(wd, traits, by_env = TRUE) {
  long <- wd %>% mutate(Environment = as.character(Environment)) %>% select(Environment, all_of(traits)) %>%
    pivot_longer(-Environment, names_to = "Trait", values_to = "Value")
  rows <- describe_values(long, "Trait", quartiles = TRUE) %>% mutate(Environment = "All")
  if (by_env && n_distinct(long$Environment) > 1) rows <- bind_rows(rows, describe_values(long, c("Trait", "Environment"), quartiles = TRUE))
  rows %>% filter(N > 0) %>%                                   # a trait that was not measured in an environment has no row
    mutate(Trait = factor(Trait, levels = traits)) %>% arrange(Trait, Environment != "All", Environment) %>%
    mutate(Trait = as.character(Trait)) %>% select(Trait, Environment, everything())
}

tidy_desc <- function(df) rename_cols(df, c(CV = "CV (%)", Kurtosis = "Excess kurtosis", Q1 = "Q1 (25%)", Q3 = "Q3 (75%)"))

# One row per accession: mean and rank in every environment, then overall mean and rank
tidy_means_ranks <- function(env_means) {
  em <- env_means %>% mutate(Genotype = as.character(Genotype), Environment = as.character(Environment))
  out <- data.frame(Accession = sort(unique(em$Genotype)), stringsAsFactors = FALSE)
  for (e in unique(em$Environment)) {
    s <- em[em$Environment == e, ]; i <- match(out$Accession, s$Genotype)
    out[[paste(e, "mean")]] <- s$Mean[i]
    out[[paste(e, "rank")]] <- rank(-s$Mean, ties.method = "min", na.last = "keep")[i]
  }
  ov <- tapply(em$Mean, em$Genotype, mean)
  out[["Overall mean"]] <- round(as.numeric(ov[out$Accession]), 3)
  out[["Overall rank"]] <- rank(-out[["Overall mean"]], ties.method = "min")
  out[order(out[["Overall rank"]]), ]
}

# BLUEs: one row per accession; BLUE and SE per environment when there are several
tidy_blues <- function(df) {
  df <- df %>% mutate(Genotype = as.character(Genotype), Environment = as.character(Environment))
  envs <- unique(df$Environment)
  if (length(envs) == 1) {
    out <- data.frame(Accession = df$Genotype, BLUE = round(df$BLUE, 3), SE = round(df$SE, 3), stringsAsFactors = FALSE)
    out$Rank <- rank(-out$BLUE, ties.method = "min")
    return(out[order(out$Rank), ])
  }
  out <- data.frame(Accession = sort(unique(df$Genotype)), stringsAsFactors = FALSE)
  for (e in envs) {
    s <- df[df$Environment == e, ]; i <- match(out$Accession, s$Genotype)
    out[[paste(e, "BLUE")]] <- round(s$BLUE[i], 3); out[[paste(e, "SE")]] <- round(s$SE[i], 3)
  }
  mb <- rowMeans(as.matrix(out[paste(envs, "BLUE")]), na.rm = TRUE)
  out[["Mean BLUE"]] <- round(mb, 3); out$Rank <- rank(-mb, ties.method = "min")
  out[order(out$Rank), ]
}

# One BLUE table: within-environment BLUEs, then the mixed-model BLUE over all environments, ranked on the latter
tidy_blues_all <- function(env_blues, mm) {
  if (is.null(env_blues)) return(tidy_blues(mm))
  out <- tidy_blues(env_blues)
  if (is.null(mm) || n_distinct(env_blues$Environment) < 2) return(out)
  mm <- mm %>% mutate(Genotype = as.character(Genotype))
  across <- if (n_distinct(mm$Environment) == 1) mm[, c("Genotype", "BLUE", "SE")] else
    mm %>% group_by(Genotype) %>% summarise(BLUE = mean(BLUE), SE = NA_real_, .groups = "drop")
  i <- match(out$Accession, across$Genotype)
  out[["Mean BLUE"]] <- NULL; out$Rank <- NULL
  out[["All-environment BLUE"]] <- round(across$BLUE[i], 3)
  if (any(!is.na(across$SE))) out[["All-environment SE"]] <- round(across$SE[i], 3)
  out$Rank <- rank(-out[["All-environment BLUE"]], ties.method = "min", na.last = "keep")
  out[order(out$Rank), ]
}

tidy_env_cor <- function(cm) data.frame(Environment = rownames(cm), round(cm, 2), check.names = FALSE, row.names = NULL)


tidy_herit <- function(h) rename_cols(h, setNames(
  c("Accessions", "Genotypic variance (Vg)", paste(GXE, "variance (Vge)"), "Block variance", "Error variance (Ve)",
    paste(H2, "plot basis"), paste(H2, "entry-mean")),
  c("N_Genotypes", "Vg", "Vge", "Vblk", "Ve", "H2_Plot", "H2_Mean")))

# BLUPs: one row per trait x environment x accession (best first within each); Pooled = all environments together
tidy_blups <- function(b) {
  out <- data.frame(Trait = b$Trait, Environment = b$Environment, Accession = b$Accession, Plots = b$Plots,
                    Environments = b$Environments, `Raw mean` = signif(b$Raw_Mean, 6), BLUP = signif(b$BLUP, 5),
                    `SE of BLUP` = signif(b$SE, 4), Reliability = round(b$Reliability, 3),
                    `Predicted value` = signif(b$Predicted, 6), Rank = b$Rank, check.names = FALSE, stringsAsFactors = FALSE)
  out[order(match(out$Trait, unique(out$Trait)), match(out$Environment, unique(out$Environment)), out$Rank), ]
}

# Variance components and heritability of the BLUP models, one row per trait x environment
tidy_blup_summary <- function(s) {
  out <- data.frame(Trait = s$Trait, Environment = s$Environment, Accessions = s$Accessions,
                    `Genotypic variance (Vg)` = signif(s$Vg, 4), `G×E variance (Vge)` = signif(s$Vge, 4),
                    `Error variance (Ve)` = signif(s$Ve, 4),
                    `H² Cullis (BLUP)` = round(s$H2_Cullis, 3), `Mean reliability` = round(s$Mean_Reliability, 3),
                    `H² entry-mean (classical)` = round(s$H2_Mean, 3), Note = s$Note, check.names = FALSE, stringsAsFactors = FALSE)
  out
}

# One PCA table: trait loadings, then the eigenvalue and variance explained of each component (first 10 components)
tidy_pca <- function(pca) {
  k <- min(10, nrow(pca$variance)); pcs <- paste0("PC", seq_len(k))
  load <- data.frame(`Trait / statistic` = pca$loadings$Trait, round(pca$loadings[pcs], 3) + 0, check.names = FALSE)   # + 0 turns -0 into 0
  v <- pca$variance[seq_len(k), ]
  stat <- data.frame(`Trait / statistic` = c("Eigenvalue", "Variance explained (%)", "Cumulative variance (%)"),
                     rbind(round(v$Variance, 3), round(v$Prop_Variance * 100, 2), round(v$Cumulative * 100, 2)),
                     check.names = FALSE)
  names(stat)[-1] <- pcs
  rbind(load, stat)
}

# One row per accession: its cluster and its scores on the first principal components
# One clustering table: the mean of each cluster (first rows), then every accession with its cluster, its scores on
# the first principal components and its trait means.
tidy_cluster_members <- function(clust, pca, mv) {
  acc <- as.character(clust$clusters$Genotype)
  out <- data.frame(Accession = acc, Cluster = as.character(clust$clusters$Cluster), stringsAsFactors = FALSE)
  pcs <- intersect(paste0("PC", 1:5), names(pca$scores)); i <- match(acc, as.character(pca$scores$Genotype))
  for (p in pcs) out[[p]] <- pca$scores[[p]][i]
  j <- match(acc, as.character(mv$data$Genotype))
  for (tr in mv$traits) out[[tr]] <- mv$data[[tr]][j]
  num <- c(pcs, mv$traits)
  means <- do.call(rbind, lapply(sort(unique(out$Cluster)), function(k) {
    s <- out[out$Cluster == k, ]
    cbind(data.frame(Accession = sprintf("Mean of cluster %s (n = %d)", k, nrow(s)), Cluster = k, stringsAsFactors = FALSE),
          data.frame(lapply(s[num], mean), check.names = FALSE))
  }))
  out <- out[order(out$Cluster, out$Accession), ]
  out <- rbind(means, out); out[num] <- round(out[num], 3) + 0
  rownames(out) <- NULL; out
}

# One stability table, one row per accession: stability parameters, AMMI scores (IPCA, ASV) and the G x E effect
# in every environment. Rows "[Environment] ..." at the bottom give the AMMI scores of the environments.
tidy_stability <- function(stab) {
  p <- stab$parameters[order(stab$parameters$Stability_Rank), ]
  out <- data.frame(Accession = p$Genotype, Mean = round(p$Overall_Mean, 3), `CV (%)` = p$CV, Range = round(p$Range, 3),
                    `Wricke ecovalence (Wi)` = p$Ecovalence, `Regression slope (b)` = p$Slope_b,
                    `R² of regression` = p$Reg_R2, `Stability rank` = p$Stability_Rank,
                    stringsAsFactors = FALSE, check.names = FALSE)
  sc <- stab$ammi_scores
  if (!is.null(sc)) {
    g <- sc[sc$Type == "Genotype", ]; i <- match(out$Accession, g$Name)
    out$IPCA1 <- round(g$IPCA1[i], 3); out$IPCA2 <- round(g$IPCA2[i], 3); out$ASV <- round(g$ASV[i], 3)
    out[["ASV rank"]] <- rank(out$ASV, ties.method = "min", na.last = "keep")
  }
  a <- stab$ammi_data
  for (e in unique(as.character(a$Environment))) {
    s <- a[as.character(a$Environment) == e, ]
    out[[paste(e, GXE, "effect")]] <- round(s$GE_Interaction[match(out$Accession, as.character(s$Genotype))], 3)
  }
  if (!is.null(sc) && any(sc$Type == "Environment")) {
    e <- sc[sc$Type == "Environment", ]
    env_rows <- out[rep(NA_integer_, nrow(e)), ]
    env_rows$Accession <- paste("[Environment]", e$Name); env_rows$Mean <- round(e$Mean, 3)
    env_rows$IPCA1 <- round(e$IPCA1, 3); env_rows$IPCA2 <- round(e$IPCA2, 3)
    out <- rbind(out, env_rows)
  }
  rownames(out) <- NULL
  out
}

# =============================================================================
# 7. EXCEL WORKBOOKS (styled sheets, a contents page, one ANOVA sheet per trait)
# =============================================================================
HEADER_STYLE <- function() createStyle(fontColour = "#FFFFFF", fgFill = "#1E3A8A", textDecoration = "bold",
                                       halign = "center", valign = "center", wrapText = TRUE)

# column widths that fit the header and the first rows, never absurdly wide
auto_widths <- function(df) vapply(seq_along(df), function(i) {
  x <- c(names(df)[i], utils::head(as.character(df[[i]]), 200))
  min(45, max(10, max(nchar(x), na.rm = TRUE) + 2))
}, numeric(1))

# One "ANOVA" sheet: for every trait a consolidated table (sources down, seasons and pooled as grouped columns)
add_anova_sheet <- function(wb, ta, sh = "ANOVA") {
  if (is.null(ta) || is.null(ta$wide)) return(NULL)
  addWorksheet(wb, sh)
  title_st <- createStyle(textDecoration = "bold", fontSize = 12, fgFill = "#DBEAFE")
  hdr <- HEADER_STYLE()
  total_st <- createStyle(textDecoration = "bold", border = "top")
  foot_st <- createStyle(fontColour = "#374151", fgFill = "#F3F4F6", halign = "center")
  row <- 1; maxcol <- 2
  for (tr in names(ta$wide)) {
    w <- ta$wide[[tr]]; tb <- w$table; nc <- ncol(tb); nf <- length(ANOVA_FOOTER); nd <- nrow(tb) - nf
    maxcol <- max(maxcol, nc)
    writeData(wb, sh, paste(tr, "— ANOVA"), startRow = row, startCol = 1)
    addStyle(wb, sh, title_st, rows = row, cols = 1:nc, gridExpand = TRUE)
    labels <- "Source"
    for (g in names(w$groups)) {
      idx <- match(w$groups[[g]], names(tb))
      writeData(wb, sh, if (g == "Pooled") "Pooled over environments" else g, startRow = row + 1, startCol = idx[1])
      mergeCells(wb, sh, cols = idx, rows = row + 1)
      labels <- c(labels, substring(w$groups[[g]], nchar(g) + 3))
    }
    addStyle(wb, sh, hdr, rows = row + 1:2, cols = 1:nc, gridExpand = TRUE)
    writeData(wb, sh, as.data.frame(t(labels)), startRow = row + 2, colNames = FALSE)
    writeData(wb, sh, tb[seq_len(nd), ], startRow = row + 3, colNames = FALSE)
    addStyle(wb, sh, total_st, rows = row + 2 + which(tb$Source[seq_len(nd)] == "Total"), cols = 1:nc, gridExpand = TRUE)
    for (k in seq_len(nf)) {
      fr <- row + 2 + nd + k
      writeData(wb, sh, tb$Source[nd + k], startRow = fr, startCol = 1)
      for (g in names(w$groups)) {
        idx <- match(w$groups[[g]], names(tb))
        writeData(wb, sh, tb[[idx[1]]][nd + k], startRow = fr, startCol = idx[1])
        mergeCells(wb, sh, cols = idx, rows = fr)
      }
      addStyle(wb, sh, foot_st, rows = fr, cols = 2:nc, gridExpand = TRUE)
    }
    row <- row + 3 + nrow(tb) + 1
  }
  setColWidths(wb, sh, 1, 44); setColWidths(wb, sh, 2:maxcol, 10)
  sh
}

write_anova_workbook <- function(ta, path) {
  wb <- createWorkbook(); add_anova_sheet(wb, ta); saveWorkbook(wb, path, overwrite = TRUE); path
}

# items: named list  sheet -> list(df, label).  Writes a "Contents" sheet first.
write_styled_workbook <- function(items, path, ta = NULL) {
  wb <- createWorkbook(); addWorksheet(wb, "Contents")
  contents <- data.frame(Sheet = character(0), `What it contains` = character(0), check.names = FALSE)
  for (sh in names(items)) {
    it <- items[[sh]]; df <- as.data.frame(it$df, check.names = FALSE)
    addWorksheet(wb, sh); writeData(wb, sh, df, headerStyle = HEADER_STYLE())
    freezePane(wb, sh, firstRow = TRUE); setColWidths(wb, sh, seq_along(df), widths = auto_widths(df))
    contents[nrow(contents) + 1, ] <- c(sh, it$label)
  }
  if (!is.null(add_anova_sheet(wb, ta)))
    contents[nrow(contents) + 1, ] <- c("ANOVA", "ANOVA of every trait: one table per trait, seasons and pooled analysis side by side")
  writeData(wb, "Contents", contents, headerStyle = HEADER_STYLE())
  setColWidths(wb, "Contents", 1:2, widths = c(28, 90))
  saveWorkbook(wb, path, overwrite = TRUE)
  path
}

# Order of the sheets in the workbook
SHEET_ORDER <- c("Design_Check", "Descriptive_Statistics", "Env_Means_Ranks", "Env_Correlations", "BLUEs",
                 "Heritability", "BLUPs", "BLUP_Summary", "Trait_Correlations", "PCA", "Clusters", "Stability_GxE")

# Writes the workbooks for everything run so far in `outdir` (from the state saved by run_phenotypic):
# Phenotypic_Analysis.xlsx (all tables) and ANOVA_tables.xlsx (the ANOVA tables, when an ANOVA was run).
export_workbook <- function(outdir) {
  state_file <- file.path(outdir, "state.rds")
  if (!file.exists(state_file)) stop("Nothing to export yet: run an analysis first.", call. = FALSE)
  st <- readRDS(state_file)
  if (is.null(st$ta)) unlink(file.path(outdir, "ANOVA_tables.xlsx"))     # no ANOVA (any more): no stale ANOVA file
  sheets <- st$sheets[order(match(names(st$sheets), SHEET_ORDER), na.last = TRUE)]
  blup_sheets <- sheets[intersect(c("BLUPs", "BLUP_Summary"), names(sheets))]
  if (!length(blup_sheets)) unlink(file.path(outdir, "BLUPs.xlsx"))         # no BLUPs (any more): no stale BLUP file
  list(all = write_styled_workbook(c(list(Data_Summary = st$summary), sheets), file.path(outdir, "Phenotypic_Analysis.xlsx"), st$ta),
       anova = if (!is.null(st$ta)) write_anova_workbook(st$ta, file.path(outdir, "ANOVA_tables.xlsx")),
       blups = if (length(blup_sheets)) write_styled_workbook(blup_sheets, file.path(outdir, "BLUPs.xlsx")))
}
