# LINKAGE RUNNER R BACKEND
# Generated from the supplied app (13)(1).R; Shiny UI/server removed.

load_and_convert_data <- function(filepath) {
  wb_raw <- tryCatch(loadWorkbook(filepath),
                     error = function(e) stop("Cannot open workbook: ", e$message))
  df_raw <- read.xlsx(filepath, sheet = 1, colNames = TRUE)
  if (!("Chromosome" %in% colnames(df_raw))) stop("Column 'Chromosome' not found.")
  if (!("Marker"     %in% colnames(df_raw))) stop("Column 'Marker' not found.")

  chr_col     <- as.character(df_raw$Chromosome)
  has_formula <- any(grepl("^=", chr_col, ignore.case = FALSE), na.rm = TRUE)
  if (has_formula) {
    inferred   <- infer_chromosome_from_marker(df_raw$Marker)
    n_inferred <- sum(!is.na(inferred))
    if (n_inferred / nrow(df_raw) >= 0.5) {
      df_raw$Chromosome <- inferred
    } else {
      df_raw$Chromosome <- seq_len(nrow(df_raw))
      warning("Could not infer chromosomes from marker names.")
    }
  }

  if (ncol(df_raw) < 3) {
    df_raw$.input_order <- seq_len(nrow(df_raw))
    return(df_raw)
  }

  geno_cols   <- df_raw[, 3:ncol(df_raw), drop = FALSE]
  sample_vals <- unique(unlist(geno_cols[1:min(50, nrow(geno_cols)), ]))
  sample_vals <- sample_vals[!is.na(sample_vals)]

  numeric_encoding <- length(sample_vals) > 0 && all(sample_vals %in% c(-1, 0, 1, 2, NA))
  text_encoding    <- length(sample_vals) > 0 && all(sample_vals %in% c("a", "b", "ab", "-", NA, "NA"))

  if (numeric_encoding) {
    df_raw[, 3:ncol(df_raw)] <- lapply(
      df_raw[, 3:ncol(df_raw), drop = FALSE],
      function(col) {
        col <- as.character(col)
        # Standardized to match Mendelian Segregation / QTL Mapping:
        # 0 = AA (homozygous Parent 1), 1 = Aa (heterozygous),
        # 2 = aa (homozygous Parent 2). Previously this was 2/1/0
        # respectively — a different, non-standard mapping that forced
        # users to manually recode the same file per tool.
        col[col == "0"]  <- "a"
        col[col == "1"]  <- "ab"
        col[col == "2"]  <- "b"
        col[col == "-1"] <- "-"
        col[is.na(col) | col == "NA"] <- "-"
        col
      })
  } else if (!text_encoding) {
    df_raw[, 3:ncol(df_raw)] <- lapply(
      df_raw[, 3:ncol(df_raw), drop = FALSE],
      function(col) {
        col <- as.character(col)
        col[!(col %in% c("a", "b", "ab"))] <- "-"
        col
      })
  }

  df_raw$Chromosome <- as.character(df_raw$Chromosome)
  df_raw$Chromosome[is.na(df_raw$Chromosome) | df_raw$Chromosome == "NA"] <- "Unknown"
  df_raw$.input_order <- seq_len(nrow(df_raw))
  df_raw
}

infer_chromosome_from_marker <- function(marker_names) {
  markers <- as.character(marker_names)
  chrom   <- rep(NA_character_, length(markers))
  patterns <- list(
    list(regex = "(?i)chr(\\d{1,2})[_\\-]",                 group = 1),
    list(regex = "(?i)[_\\-]C(\\d{1,2})[_\\-]",            group = 1),
    list(regex = "(?i)^MSU\\d+_(\\d{1,2})_",               group = 1),
    list(regex = "(?i)^DTY(\\d{1,2})[\\-_]",               group = 1),
    list(regex = "(?i)^[A-Za-z]{1,6}(\\d{1,2})[a-zA-Z_]", group = 1),
    list(regex = "(?i)^Saltol",                             group = NA),
    list(regex = "(?i)^SCT(\\d{1,2})[_\\-]",               group = 1),
    list(regex = "(?i)^R(\\d{1,2})",                        group = 1),
    list(regex = "(?i)^Pi(\\d{1,2})",                       group = 1)
  )
  for (i in seq_along(markers)) {
    m <- markers[i]
    for (pat in patterns) {
      if (is.na(pat$group)) {
        if (grepl(pat$regex, m, perl = TRUE)) { chrom[i] <- "1"; break }
      } else {
        mt     <- regmatches(m, regexpr(pat$regex, m, perl = TRUE))
        if (length(mt) > 0) {
          sub_mt <- regmatches(mt, regexec(pat$regex, mt, perl = TRUE))[[1]]
          if (length(sub_mt) >= pat$group + 1) {
            chr_num <- as.integer(sub_mt[pat$group + 1])
            if (!is.na(chr_num) && chr_num >= 1 && chr_num <= 30) {
              chrom[i] <- as.character(chr_num); break
            }
          }
        }
      }
    }
  }
  chrom
}

# ═══════════════════════════════════════════════════════════════════════════════
#  CORE GENETIC FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

calculate_recombination <- function(population_type, genotype1, genotype2) {
  g1   <- as.character(genotype1)
  g2   <- as.character(genotype2)
  keep <- !(g1 %in% c("-", "NA", NA) | g2 %in% c("-", "NA", NA))
  g1   <- g1[keep]; g2 <- g2[keep]
  n    <- length(g1)
  if (n < 4) return(NA_real_)

  if (population_type == "RIL") {
    hom <- !(g1 == "ab" | g2 == "ab")
    g1h <- g1[hom]; g2h <- g2[hom]
    if (length(g1h) < 4) return(NA_real_)
    r_obs <- mean(g1h != g2h)
    rf    <- r_obs / (2 * (1 - r_obs))
    return(min(rf, 0.4999))
  } else if (population_type == "F2") {
    AA_AA <- sum(g1=="a"  & g2=="a");  AA_Aa <- sum(g1=="a"  & g2=="ab")
    AA_BB <- sum(g1=="a"  & g2=="b");  Aa_AA <- sum(g1=="ab" & g2=="a")
    Aa_Aa <- sum(g1=="ab" & g2=="ab"); Aa_BB <- sum(g1=="ab" & g2=="b")
    BB_AA <- sum(g1=="b"  & g2=="a");  BB_Aa <- sum(g1=="b"  & g2=="ab")
    BB_BB <- sum(g1=="b"  & g2=="b")
    r <- 0.25
    for (iter in 1:200) {
      r_old <- r; q <- 1 - r
      p_parental <- q^2 + r^2; p_recomb <- 2 * r * q
      norm <- p_parental + p_recomb
      n_recomb_AaAa <- if (norm > 0) Aa_Aa * p_recomb / norm else 0
      recomb_count  <- 2*(AA_BB+BB_AA) + (AA_Aa+Aa_AA+Aa_BB+BB_Aa) + 2*n_recomb_AaAa
      r_new <- max(1e-6, min(recomb_count / (2*n), 0.4999))
      r <- r_new
      if (abs(r - r_old) < 1e-8) break
    }
    return(r)
  } else stop("population_type must be 'RIL' or 'F2'")
}

calculate_cM_kosambi <- function(r) {
  if (is.na(r)) return(NA_real_)
  r <- max(0, min(r, 0.4999))
  if (r == 0) return(0)
  25 * log((1 + 2*r) / (1 - 2*r))
}

detect_monomorphic_markers <- function(df) {
  geno_cols_idx <- setdiff(3:ncol(df), which(colnames(df) == ".input_order"))
  if (length(geno_cols_idx) == 0) return(setNames(rep(FALSE, nrow(df)), df$Marker))
  is_mono <- apply(df[, geno_cols_idx, drop = FALSE], 1, function(row) {
    vals        <- row[!is.na(row) & !(row %in% c("-", "NA", ""))]
    if (length(vals) == 0) return(FALSE)
    unique_vals <- unique(vals[vals %in% c("a", "b", "ab")])
    length(unique_vals) <= 1
  })
  names(is_mono) <- df$Marker
  is_mono
}

perform_chi_square_test <- function(df, population_type) {
  geno_cols_idx <- setdiff(3:ncol(df), which(colnames(df) == ".input_order"))
  if (population_type == "F2") {
    expected <- c(a = 0.25, ab = 0.50, b = 0.25)
    valid_classes <- c("a", "ab", "b")
  } else {
    expected <- c(a = 0.50, b = 0.50)
    valid_classes <- c("a", "b")
  }
  p_values <- apply(df[, geno_cols_idx, drop = FALSE], 1, function(row) {
    vals <- row[!is.na(row) & row %in% valid_classes]
    if (length(vals) < 10) return(1.0)
    obs    <- table(factor(vals, levels = valid_classes))
    result <- tryCatch(suppressWarnings(chisq.test(obs, p = expected)),
                       error = function(e) list(p.value = 1.0))
    result$p.value
  })
  names(p_values) <- df$Marker
  p_values
}

correct_errors_and_impute_missing <- function(df) {
  valid_geno    <- c("a", "b", "ab")
  geno_start    <- which(colnames(df) == "Chromosome") + 1
  if (length(geno_start) == 0 || geno_start > ncol(df)) return(df)
  geno_cols_idx <- setdiff(geno_start:ncol(df), which(colnames(df) == ".input_order"))
  for (j in geno_cols_idx) {
    x <- as.character(df[[j]])
    x[x %in% c("-", "NA", "", " ", "na")] <- NA
    x[!is.na(x) & !(x %in% valid_geno)]   <- NA
    df[[j]] <- x
  }
  if (length(geno_cols_idx) > 0) {
    na_col   <- colMeans(is.na(df[, geno_cols_idx, drop = FALSE]))
    keep_idx <- c(seq_len(geno_start - 1),
                  geno_cols_idx[na_col < 1],
                  which(colnames(df) == ".input_order"))
    keep_idx <- keep_idx[keep_idx <= ncol(df)]
    df       <- df[, sort(unique(keep_idx)), drop = FALSE]
  }
  geno_cols_idx2 <- setdiff(geno_start:ncol(df), which(colnames(df) == ".input_order"))
  if (length(geno_cols_idx2) > 0) {
    df[, geno_cols_idx2] <- t(apply(df[, geno_cols_idx2, drop = FALSE], 1, function(row) {
      nas <- is.na(row)
      if (all(nas) || !any(nas)) return(row)
      tbl      <- table(row[!nas])
      mode_val <- names(which.max(tbl))
      row[nas] <- mode_val
      row
    }))
  }
  df
}

quality_control <- function(df, missing_threshold = 0.1) {
  if (ncol(df) < 3 || nrow(df) == 0) return(df[0, ])
  geno_cols_idx <- setdiff(3:(ncol(df)), which(colnames(df) == ".input_order"))
  if (length(geno_cols_idx) == 0) return(df)
  geno      <- df[, geno_cols_idx, drop = FALSE]
  miss_frac <- rowMeans(is.na(geno))
  miss_frac[is.na(miss_frac)] <- 1
  df[miss_frac <= missing_threshold, ]
}

# ═══════════════════════════════════════════════════════════════════════════════
#  POSITION BUILDING
# ═══════════════════════════════════════════════════════════════════════════════

build_positions_from_input_order <- function(df, population_type, recombination_threshold) {
  chromosomes   <- unique(df$Chromosome[order(df$.input_order)])
  geno_cols_idx <- setdiff(3:ncol(df), which(colnames(df) == ".input_order"))
  positions_list     <- list()
  recombination_rows <- list()

  for (chr in chromosomes) {
    sub <- df[df$Chromosome == chr, ]
    sub <- sub[order(sub$.input_order), ]
    mks <- sub$Marker
    if (length(mks) == 0) next
    pos <- setNames(numeric(length(mks)), mks)
    if (length(mks) > 1) {
      for (i in 2:length(mks)) {
        m_prev <- mks[i-1]; m_curr <- mks[i]
        g1 <- unlist(sub[sub$Marker == m_prev, geno_cols_idx])
        g2 <- unlist(sub[sub$Marker == m_curr, geno_cols_idx])
        rf <- calculate_recombination(population_type, g1, g2)
        if (!is.na(rf)) {
          cM <- calculate_cM_kosambi(rf)
          recombination_rows[[length(recombination_rows)+1]] <-
            data.frame(Marker1=m_prev, Marker2=m_curr,
                       RecombinationFraction=rf, MapDistance_cM=cM,
                       Chromosome=chr, stringsAsFactors=FALSE)
          pos[i] <- if (rf < recombination_threshold) pos[i-1] + abs(cM)
                    else                               pos[i-1] + 50
        } else {
          pos[i] <- pos[i-1] + 5
        }
      }
    }
    for (mk in names(pos)) positions_list[[mk]] <- pos[mk]
  }
  list(positions_list=positions_list, recombination_rows=recombination_rows)
}

scale_linkage_map <- function(linkage_map, max_total_cm, method = "proportional") {
  pos_map <- linkage_map[!is.na(linkage_map$Position), ]
  if (nrow(pos_map) == 0)
    return(list(scaled_map=linkage_map, chrom_lengths_before=data.frame(),
                chrom_lengths_after=data.frame()))
  before <- pos_map %>%
    group_by(Chromosome) %>%
    summarise(Original_Length = max(Position, na.rm=TRUE), .groups="drop")
  total <- sum(before$Original_Length, na.rm=TRUE)
  if (total <= max_total_cm || total == 0)
    return(list(scaled_map=linkage_map, chrom_lengths_before=before,
                chrom_lengths_after=before))
  scaled <- linkage_map
  if (method == "proportional") {
    before$sf <- max_total_cm / total
  } else {
    n_chr <- nrow(before)
    before$sf <- ifelse(before$Original_Length > 0,
                        (max_total_cm/n_chr)/before$Original_Length, 1)
  }
  for (chr in before$Chromosome) {
    sf  <- before$sf[before$Chromosome == chr]
    idx <- which(scaled$Chromosome == chr & !is.na(scaled$Position))
    if (length(idx) > 0) {
      mn <- min(scaled$Position[idx])
      scaled$Position[idx] <- (scaled$Position[idx] - mn) * sf
    }
  }
  after <- scaled %>%
    filter(!is.na(Position)) %>%
    group_by(Chromosome) %>%
    summarise(Original_Length = max(Position, na.rm=TRUE), .groups="drop")
  list(scaled_map=scaled, chrom_lengths_before=before, chrom_lengths_after=after)
}

# ═══════════════════════════════════════════════════════════════════════════════
#  MASTER ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════

create_linkage_map <- function(df,
                               population_type         = "F2",
                               missing_threshold       = 0.1,
                               recombination_threshold = 0.5,
                               max_total_cm            = 1500,
                               scaling_method          = "proportional",
                               enable_chi_square       = FALSE,
                               chi_square_alpha        = 0.05,
                               use_bonferroni          = FALSE) {

  empty_result <- function(df_in = df) {
    list(original_data        = df_in,
         recombination_data   = data.frame(),
         linkage_map          = data.frame(Marker=character(0),
                                           Chromosome=character(0),
                                           Position=numeric(0)),
         all_markers_positions = data.frame(),
         retained_markers_list = data.frame(),
         chrom_lengths_before = data.frame(),
         chrom_lengths_after  = data.frame(),
         removal_log          = data.frame(),
         qc_report            = "No data passed QC.")
  }

  df_orig <- df
  if (!(".input_order" %in% colnames(df))) df$.input_order <- seq_len(nrow(df))

  removal_log <- data.frame(
    Marker     = character(0), Chromosome = character(0),
    Reason     = character(0), Detail     = character(0),
    stringsAsFactors = FALSE
  )

  n_input <- nrow(df)

  df <- correct_errors_and_impute_missing(df)

  geno_cols_idx <- setdiff(3:ncol(df), which(colnames(df) == ".input_order"))
  geno_mat      <- df[, geno_cols_idx, drop = FALSE]
  miss_frac     <- if (ncol(geno_mat) > 0) rowMeans(is.na(geno_mat)) else rep(0, nrow(df))
  miss_frac[is.na(miss_frac)] <- 1
  failed_miss   <- df$Marker[miss_frac > missing_threshold]

  if (length(failed_miss) > 0) {
    removal_log <- rbind(removal_log, data.frame(
      Marker     = failed_miss,
      Chromosome = df$Chromosome[miss_frac > missing_threshold],
      Reason     = "Excess missing data",
      Detail     = paste0(round(miss_frac[miss_frac > missing_threshold]*100,1),
                          "% missing (threshold: ", round(missing_threshold*100,1), "%)"),
      stringsAsFactors = FALSE
    ))
  }
  df <- df[miss_frac <= missing_threshold, ]
  if (nrow(df) == 0 || ncol(df) < 3) {
    res <- empty_result(df_orig); res$removal_log <- removal_log; return(res)
  }

  is_mono      <- detect_monomorphic_markers(df)
  mono_markers <- names(is_mono)[is_mono]
  if (length(mono_markers) > 0) {
    mono_df       <- df[df$Marker %in% mono_markers, ]
    geno_cols_idx2 <- setdiff(3:ncol(df), which(colnames(df) == ".input_order"))
    details <- sapply(mono_markers, function(mk) {
      row  <- df[df$Marker == mk, geno_cols_idx2, drop = FALSE]
      vals <- unlist(row)
      vals <- vals[!is.na(vals) & vals %in% c("a","b","ab")]
      if (length(vals) == 0) return("No informative genotypes")
      tbl      <- table(vals)
      dominant <- names(which.max(tbl))
      allele_name <- switch(dominant,
                            a  = "Homozygous Parent 1 (AA)",
                            b  = "Homozygous Parent 2 (BB)",
                            ab = "Heterozygous only (ab)",
                            dominant)
      paste0("Fixed as: ", allele_name,
             " (", round(max(tbl)/sum(tbl)*100, 1), "% of genotyped individuals)")
    })
    removal_log <- rbind(removal_log, data.frame(
      Marker     = mono_markers,
      Chromosome = mono_df$Chromosome,
      Reason     = "Monomorphic marker",
      Detail     = details,
      stringsAsFactors = FALSE
    ))
    df <- df[!df$Marker %in% mono_markers, ]
  }

  if (nrow(df) == 0) {
    res <- empty_result(df_orig); res$removal_log <- removal_log; return(res)
  }

  n_after_mono  <- nrow(df)
  n_failed_chi  <- 0
  chi_note       <- "Chi-square filter : DISABLED"

  if (enable_chi_square) {
    p_values   <- perform_chi_square_test(df, population_type)
    alpha_used <- if (use_bonferroni) chi_square_alpha / nrow(df) else chi_square_alpha
    failed_chi_markers <- names(p_values)[p_values < alpha_used]
    if (length(failed_chi_markers) > 0) {
      chi_df      <- df[df$Marker %in% failed_chi_markers, ]
      chi_details <- sapply(failed_chi_markers, function(mk) {
        pv <- p_values[mk]
        paste0("p = ", formatC(pv, format="e", digits=3),
               " (alpha = ", formatC(alpha_used, format="e", digits=3), ")")
      })
      removal_log <- rbind(removal_log, data.frame(
        Marker     = failed_chi_markers,
        Chromosome = chi_df$Chromosome,
        Reason     = "Segregation distortion (chi-square)",
        Detail     = chi_details,
        stringsAsFactors = FALSE
      ))
      df <- df[!df$Marker %in% failed_chi_markers, ]
      n_failed_chi <- length(failed_chi_markers)
    }
    bonf_label <- if (use_bonferroni)
      paste0("Bonferroni-corrected (alpha/", n_after_mono, ")")
    else "Uncorrected"
    chi_note <- paste0(
      "Chi-square filter : ENABLED\n",
      "  Significance    : alpha = ", chi_square_alpha, " (", bonf_label, ")\n",
      "  Effective alpha : ", formatC(alpha_used, format="e", digits=3), "\n",
      "  Markers removed : ", n_failed_chi
    )
  }

  if (nrow(df) == 0) {
    res <- empty_result(df_orig); res$removal_log <- removal_log; return(res)
  }

  n_after_qc <- nrow(df)

  qc_msg <- paste0(
    "QC Report (v6.0)\n",
    "================\n",
    "  Markers input              : ", n_input, "\n",
    "  Removed (missing data)     : ", length(failed_miss), "\n",
    "  Removed (monomorphic)      : ", length(mono_markers), "\n",
    "  Removed (chi-square)       : ", n_failed_chi, "\n",
    "  Markers retained           : ", n_after_qc, "\n\n",
    chi_note, "\n\n",
    "  Marker ORDER preserved from input file.\n",
    "  First marker of each chromosome = 0 cM.\n"
  )

  pos_result         <- build_positions_from_input_order(df, population_type, recombination_threshold)
  positions_list     <- pos_result$positions_list
  recombination_rows <- pos_result$recombination_rows
  rec_df <- if (length(recombination_rows) > 0) do.call(rbind, recombination_rows)
            else                                 data.frame()

  if (length(positions_list) == 0) {
    lm_empty <- df[, c("Marker","Chromosome")]
    lm_empty$Position <- NA_real_
    retained_list <- data.frame(Marker=df$Marker, Chromosome=df$Chromosome,
                                Status="Retained", stringsAsFactors=FALSE)
    all_markers   <- build_all_markers_positions(df_orig, removal_log, lm_empty)
    return(list(original_data=df_orig, recombination_data=rec_df,
                linkage_map=lm_empty,
                all_markers_positions=all_markers,
                retained_markers_list=retained_list,
                chrom_lengths_before=data.frame(),
                chrom_lengths_after=data.frame(),
                removal_log=removal_log, qc_report=qc_msg))
  }

  pos_df <- data.frame(Marker   = names(positions_list),
                       Position = unlist(positions_list),
                       stringsAsFactors = FALSE)

  linkage_map <- merge(df[, c("Marker","Chromosome",".input_order")],
                       pos_df, by="Marker", all.x=TRUE)
  linkage_map <- linkage_map[order(linkage_map$Chromosome, linkage_map$.input_order), ]
  linkage_map$.input_order <- NULL

  scaled <- scale_linkage_map(linkage_map, max_total_cm, scaling_method)

  all_markers_positions <- build_all_markers_positions(df_orig, removal_log, scaled$scaled_map)

  retained_markers_list <- scaled$scaled_map %>%
    select(Marker, Chromosome, Position) %>%
    mutate(Status = "Retained") %>%
    arrange(Chromosome, Position)

  list(original_data         = df_orig,
       recombination_data    = rec_df,
       linkage_map           = scaled$scaled_map,
       all_markers_positions = all_markers_positions,
       retained_markers_list = retained_markers_list,
       chrom_lengths_before  = scaled$chrom_lengths_before,
       chrom_lengths_after   = scaled$chrom_lengths_after,
       removal_log           = removal_log,
       qc_report             = qc_msg)
}

build_all_markers_positions <- function(df_orig, removal_log, linkage_map) {
  df_all <- df_orig[, c("Marker","Chromosome"), drop=FALSE]
  df_all$Chromosome <- as.character(df_all$Chromosome)
  lm_sub <- linkage_map[, c("Marker","Position"), drop=FALSE]
  df_all <- merge(df_all, lm_sub, by="Marker", all.x=TRUE)
  removed_markers <- unique(removal_log$Marker)
  df_all$Status <- ifelse(df_all$Marker %in% removed_markers, "Removed", "Retained")
  if (nrow(removal_log) > 0) {
    reason_df <- removal_log[, c("Marker","Reason","Detail"), drop=FALSE]
    df_all <- merge(df_all, reason_df, by="Marker", all.x=TRUE)
  } else {
    df_all$Reason <- NA_character_
    df_all$Detail <- NA_character_
  }
  df_all$Position[df_all$Status == "Removed"] <- NA_real_
  df_all <- df_all[order(df_all$Chromosome, df_all$Marker), ]
  df_all[, c("Marker","Chromosome","Status","Position","Reason","Detail")]
}

# ═══════════════════════════════════════════════════════════════════════════════
#  CYTOGENETIC PLOT DATA BUILDER (FIXED)
# ═══════════════════════════════════════════════════════════════════════════════

build_cyto_plot_data <- function(linkage_map, df_retained, selected_chromosomes) {
  if (is.null(selected_chromosomes) || length(selected_chromosomes) == 0) return(list())
  
  # Ensure linkage_map has Position column
  if (is.null(linkage_map$Position)) linkage_map$Position <- NA_real_
  
  # Get genotype columns (all columns except Marker, Chromosome, .input_order, Position, Status, Reason, Detail)
  exclude_cols <- c("Marker", "Chromosome", ".input_order", "Position", "Status", "Reason", "Detail")
  geno_cols <- setdiff(colnames(df_retained), exclude_cols)
  
  # Sort chromosomes numerically if possible
  chr_sorted <- tryCatch(
    selected_chromosomes[order(as.numeric(selected_chromosomes))],
    warning = function(w) sort(selected_chromosomes),
    error   = function(e) sort(selected_chromosomes)
  )
  
  result <- list()
  for (chr in chr_sorted) {
    # Get markers for this chromosome with valid positions
    sub_map <- linkage_map[linkage_map$Chromosome == chr & !is.na(linkage_map$Position), ]
    if (nrow(sub_map) == 0) next
    sub_map <- sub_map[order(sub_map$Position), ]
    
    # Get genotype data for these markers
    sub_geno <- df_retained[df_retained$Marker %in% sub_map$Marker, ]
    
    marker_data <- list()
    for (i in seq_len(nrow(sub_map))) {
      mk <- sub_map$Marker[i]
      pos <- sub_map$Position[i]
      
      # Get genotype row for this marker
      geno_row <- sub_geno[sub_geno$Marker == mk, geno_cols, drop = FALSE]
      
      if (nrow(geno_row) == 0 || length(geno_cols) == 0) {
        marker_data[[i]] <- list(
          marker = mk,
          position = round(pos, 3),
          allele_a = 0,
          allele_b = 0,
          allele_ab = 0,
          allele_miss = 1
        )
        next
      }
      
      vals <- unlist(geno_row)
      vals <- vals[!is.na(vals)]
      n_total <- length(vals)
      if (n_total == 0) n_total <- 1
      
      marker_data[[i]] <- list(
        marker = mk,
        position = round(pos, 3),
        allele_a = round(sum(vals == "a", na.rm = TRUE) / n_total, 4),
        allele_b = round(sum(vals == "b", na.rm = TRUE) / n_total, 4),
        allele_ab = round(sum(vals == "ab", na.rm = TRUE) / n_total, 4),
        allele_miss = round(sum(!(vals %in% c("a","b","ab")), na.rm = TRUE) / n_total, 4)
      )
    }
    
    if (length(marker_data) > 0) {
      result[[chr]] <- list(
        chromosome = chr,
        max_pos = round(max(sub_map$Position, na.rm = TRUE), 2),
        n_markers = nrow(sub_map),
        markers = marker_data
      )
    }
  }
  result
}

create_ggplot_linkage_map <- function(linkage_map,
                                      selected_chromosomes,
                                      show_labels     = FALSE,
                                      label_size      = 2.5,
                                      color_palette   = "Set2",
                                      chrom_width     = 0.4,
                                      point_size      = 2.5,
                                      chrom_spacing   = 2.5,
                                      plot_theme      = "dark") {

  if (is.null(selected_chromosomes) || length(selected_chromosomes) == 0)
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No chromosomes selected.", size=6) + theme_void())

  chr_sorted <- tryCatch(
    selected_chromosomes[order(as.numeric(selected_chromosomes))],
    warning = function(w) sort(selected_chromosomes),
    error   = function(e) sort(selected_chromosomes)
  )

  plot_data <- linkage_map %>%
    filter(Chromosome %in% selected_chromosomes, !is.na(Position)) %>%
    mutate(Chromosome = factor(Chromosome, levels = chr_sorted)) %>%
    arrange(Chromosome, Position)

  if (nrow(plot_data) == 0) {
    return(ggplot() +
             annotate("text", x=0.5, y=0.5,
                      label="No valid marker positions.\nTry relaxing QC parameters.",
                      size=5, color="#E74C3C", hjust=0.5) + theme_void())
  }

  n_chr <- length(chr_sorted)

  chrom_stats <- plot_data %>%
    group_by(Chromosome) %>%
    summarise(max_pos    = max(Position, na.rm=TRUE),
              min_pos    = min(Position, na.rm=TRUE),
              num_markers = n(), .groups="drop") %>%
    mutate(Chromosome = factor(Chromosome, levels = chr_sorted)) %>%
    arrange(Chromosome) %>%
    mutate(x_center = as.numeric(Chromosome) * chrom_spacing)

  plot_data <- plot_data %>%
    left_join(chrom_stats %>% select(Chromosome, x_center), by="Chromosome")

  # Color palette
  n_needed <- max(3, n_chr)
  chrom_colors <- tryCatch({
    if (n_chr <= 8) brewer.pal(n_needed, color_palette)[seq_len(n_chr)]
    else colorRampPalette(brewer.pal(8, color_palette))(n_chr)
  }, error = function(e) {
    colorRampPalette(brewer.pal(8, "Set2"))(n_chr)
  })
  names(chrom_colors) <- chr_sorted

  # Theme setup
  if (plot_theme == "dark") {
    bg_col    <- "#1A1A2E"
    panel_col <- "#16213E"
    grid_col  <- "#2D3561"
    text_col  <- "#E8E8F0"
    axis_col  <- "#8A8AB0"
    strip_col <- "#0F3460"
  } else if (plot_theme == "minimal_white") {
    bg_col    <- "#FFFFFF"
    panel_col <- "#FAFAFA"
    grid_col  <- "#EEEEEE"
    text_col  <- "#222222"
    axis_col  <- "#666666"
    strip_col <- "#F0F0F0"
  } else {
    bg_col    <- "#F0F4FF"
    panel_col <- "#FFFFFF"
    grid_col  <- "#D0DCF0"
    text_col  <- "#1E3A8A"
    axis_col  <- "#4A6FA5"
    strip_col <- "#DBEAFE"
  }

  # Half-width for chromosome bars
  hw <- chrom_width / 2

  p <- ggplot() +

    # Chromosome body (thick segment with rounded aesthetic via geom_rect)
    geom_rect(data = chrom_stats,
              aes(xmin = x_center - hw,
                  xmax = x_center + hw,
                  ymin = min_pos,
                  ymax = max_pos,
                  fill = as.factor(Chromosome)),
              alpha = 0.18, color = NA) +

    # Chromosome axis line
    geom_segment(data = chrom_stats,
                 aes(x = x_center, xend = x_center,
                     y = min_pos,  yend = max_pos,
                     color = as.factor(Chromosome)),
                 linewidth = chrom_width * 12,
                 alpha = 0.5, lineend = "round") +

    # Marker tick marks (horizontal lines)
    geom_segment(data = plot_data,
                 aes(x    = x_center - hw * 1.6,
                     xend = x_center + hw * 1.6,
                     y    = Position,
                     yend = Position,
                     color = as.factor(Chromosome)),
                 linewidth = 0.6, alpha = 0.85) +

    # Marker points
    geom_point(data = plot_data,
               aes(x = x_center, y = Position,
                   color = as.factor(Chromosome)),
               size = point_size, alpha = 0.92, shape = 21,
               fill = "white", stroke = 0.8) +

    # Chromosome name label (bottom)
    geom_text(data = chrom_stats,
              aes(x = x_center,
                  y = max_pos + max(chrom_stats$max_pos) * 0.025,
                  label = paste0("Chr ", Chromosome),
                  color = as.factor(Chromosome)),
              size  = 3.8, fontface = "bold", vjust = 0) +

    # Marker count label
    geom_text(data = chrom_stats,
              aes(x = x_center,
                  y = max_pos + max(chrom_stats$max_pos) * 0.07,
                  label = paste0("n=", num_markers)),
              color = axis_col, size = 3.0, vjust = 0) +

    # Length label (top)
    geom_text(data = chrom_stats,
              aes(x = x_center,
                  y = min_pos - max(chrom_stats$max_pos) * 0.025,
                  label = paste0(round(max_pos,1), " cM"),
                  color = as.factor(Chromosome)),
              size = 3.0, vjust = 1, fontface = "italic") +

    scale_color_manual(values = chrom_colors) +
    scale_fill_manual(values  = chrom_colors) +
    scale_y_reverse(name = "Genetic Position (cM)",
                    expand = expansion(mult = c(0.04, 0.15))) +
    scale_x_continuous(breaks = chrom_stats$x_center,
                        labels = chrom_stats$Chromosome,
                        name   = "Linkage Group",
                        expand = expansion(mult = 0.08)) +
    guides(color = "none", fill = "none") +
    labs(title    = "Genetic Linkage Map",
         subtitle = paste0(nrow(plot_data), " markers across ",
                           n_chr, " linkage groups"),
         caption  = "Map distances in Kosambi cM  |  Linkage Map Creator v6.0") +
    theme_minimal(base_size = 13) +
    theme(
      plot.background    = element_rect(fill = bg_col, color = NA),
      panel.background   = element_rect(fill = panel_col, color = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = grid_col, linewidth = 0.4, linetype = "dashed"),
      panel.grid.minor.y = element_line(color = grid_col, linewidth = 0.2),
      axis.text.x        = element_text(color = text_col, size = 11, face = "bold"),
      axis.text.y        = element_text(color = axis_col, size = 10),
      axis.title.x       = element_text(color = text_col, size = 13, face = "bold",
                                         margin = margin(t = 12)),
      axis.title.y       = element_text(color = text_col, size = 13, face = "bold",
                                         margin = margin(r = 12)),
      plot.title         = element_text(color = text_col, size = 18, face = "bold",
                                         hjust = 0.5, margin = margin(b = 4)),
      plot.subtitle      = element_text(color = axis_col, size = 12,
                                         hjust = 0.5, margin = margin(b = 14)),
      plot.caption       = element_text(color = axis_col, size = 9, hjust = 1),
      plot.margin        = margin(20, 30, 20, 20)
    )

  # Optional marker labels
  if (show_labels && nrow(plot_data) <= 600) {
    p <- p + ggrepel::geom_text_repel(
      data        = plot_data,
      aes(x = x_center, y = Position,
          label = Marker, color = as.factor(Chromosome)),
      size          = label_size,
      direction     = "y",
      nudge_x       = hw * 2.2,
      segment.size  = 0.25,
      segment.alpha = 0.55,
      max.overlaps  = 25,
      show.legend   = FALSE,
      fontface      = "plain"
    )
  }

  p
}

# ═══════════════════════════════════════════════════════════════════════════════
#  RECOMBINATION HEATMAP
# ═══════════════════════════════════════════════════════════════════════════════

create_recombination_heatmap <- function(rec_data, selected_chr, top_n = 60) {
  if (is.null(rec_data) || nrow(rec_data) == 0)
    return(ggplot() + annotate("text",x=.5,y=.5,label="No recombination data.",size=5) + theme_void())

  df <- rec_data
  if (!is.null(selected_chr) && length(selected_chr) > 0)
    df <- df[df$Chromosome %in% selected_chr, ]
  if (nrow(df) == 0)
    return(ggplot() + annotate("text",x=.5,y=.5,label="No data for selected chromosomes.",size=5) + theme_void())

  # Limit pairs for readability
  if (nrow(df) > top_n) df <- df[1:top_n, ]

  df$pair <- paste0(df$Marker1, "\n", df$Marker2)
  df$pair <- factor(df$pair, levels = rev(unique(df$pair)))

  ggplot(df, aes(x = 1, y = pair, fill = RecombinationFraction)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = paste0(round(MapDistance_cM, 1), " cM\nr=",
                                  round(RecombinationFraction, 3))),
              size = 2.8, color = "white", fontface = "bold") +
    scale_fill_gradientn(
      colors = c("#0D3B66","#1A6FA5","#4ECDC4","#FFE66D","#FF6B6B","#C0392B"),
      limits = c(0, 0.5),
      name   = "Recombination\nFraction (r)",
      breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
      labels = c("0.0","0.1","0.2","0.3","0.4","0.5")
    ) +
    facet_wrap(~Chromosome, scales = "free_y", ncol = 3) +
    labs(title    = "Adjacent Marker Recombination Fractions",
         subtitle = paste0("Showing ", nrow(df), " marker pairs"),
         x = NULL, y = "Marker Pair") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      axis.text.y     = element_text(size = 8),
      strip.text      = element_text(face = "bold", size = 10),
      legend.position = "right",
      plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle   = element_text(size = 10, hjust = 0.5, color = "gray50"),
      panel.grid      = element_blank(),
      plot.background = element_rect(fill = "#FAFAFA", color = NA)
    )
}

# ═══════════════════════════════════════════════════════════════════════════════
#  MARKER DENSITY PLOT
# ═══════════════════════════════════════════════════════════════════════════════

create_marker_density_plot <- function(linkage_map, selected_chr, bin_width = 10) {
  if (is.null(linkage_map) || nrow(linkage_map) == 0)
    return(ggplot() + annotate("text",x=.5,y=.5,label="No data.",size=5) + theme_void())

  df <- linkage_map %>%
    filter(!is.na(Position))

  if (!is.null(selected_chr) && length(selected_chr) > 0)
    df <- df %>% filter(Chromosome %in% selected_chr)

  if (nrow(df) == 0)
    return(ggplot() + annotate("text",x=.5,y=.5,label="No positioned markers.",size=5) + theme_void())

  chr_sorted <- tryCatch(
    sort(unique(df$Chromosome), decreasing = FALSE),
    warning = function(w) sort(unique(df$Chromosome))
  )
  df$Chromosome <- factor(df$Chromosome, levels = chr_sorted)

  n_chr <- length(chr_sorted)
  pal <- tryCatch(
    colorRampPalette(brewer.pal(min(n_chr, 9), "Set1"))(n_chr),
    error = function(e) colorRampPalette(c("#E74C3C","#3498DB","#2ECC71","#F39C12","#9B59B6"))(n_chr)
  )

  ggplot(df, aes(x = Position, fill = Chromosome, color = Chromosome)) +
    geom_histogram(binwidth = bin_width, alpha = 0.75, linewidth = 0.3) +
    geom_density(aes(y = after_stat(count) * bin_width),
                 alpha = 0, linewidth = 1.0, linetype = "dashed") +
    facet_wrap(~Chromosome, scales = "free", ncol = ceiling(sqrt(n_chr))) +
    scale_fill_manual(values  = pal) +
    scale_color_manual(values = pal) +
    labs(title    = "Marker Density Distribution",
         subtitle = paste0("Bin width: ", bin_width, " cM"),
         x        = "Genetic Position (cM)",
         y        = "Number of Markers",
         fill     = "Chromosome",
         color    = "Chromosome") +
    guides(fill = "none", color = "none") +
    theme_minimal(base_size = 12) +
    theme(
      strip.text      = element_text(face = "bold", size = 11,
                                     color = "#1E3A8A", margin = margin(b=4)),
      strip.background = element_rect(fill = "#DBEAFE", color = NA),
      panel.grid.minor = element_blank(),
      plot.title      = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle   = element_text(size = 10, hjust = 0.5, color = "gray50"),
      plot.background = element_rect(fill = "#F8FAFF", color = NA),
      panel.background = element_rect(fill = "#FFFFFF", color = NA)
    )
}

# ═══════════════════════════════════════════════════════════════════════════════


# =============================================================================
# Linkage Map Creator - R backend runner
# Ported directly from the supplied Shiny app (app (13)(1).R).
#
# This file is intentionally Shiny-free. It exposes:
#   run_linkage(config)
#   render_linkage_plots(result, config, output_dir)
#   write_linkage_outputs(result, config, output_dir)
#
# CLI mode:
#   Rscript linkage_runner.R '<json config>'
#
# Expected config:
# {
#   "input_path": "...xlsx",
#   "population_type": "F2" | "RIL",
#   "missing_threshold": 0.30,
#   "recombination_threshold": 0.45,
#   "max_total_cm": 1500,
#   "scaling_method": "proportional" | "equal",
#   "enable_chi_square": false,
#   "chi_square_alpha": 0.05,
#   "use_bonferroni": false,
#   "output_dir": "..."
# }
#
# JSON result keys are aligned to the desktop shell contract:
# linkage_map, removal_log, retained_markers, all_markers,
# linkage_map_columns, removal_log_columns, retained_markers_columns,
# all_markers_columns, chromosomes, qc_report, scaling_summary,
# recombination_data, recombination_data_columns, cyto_data,
# map_plot, density_plot, heatmap_plot, result_xlsx, all_markers_xlsx,
# retained_markers_xlsx, removal_log_xlsx.
# =============================================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(ggrepel)
  library(scales)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# -----------------------------------------------------------------------------
# Helpers used by the desktop bridge
# -----------------------------------------------------------------------------

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  if (is.logical(x)) return(isTRUE(x))
  if (is.numeric(x)) return(x != 0)
  tolower(as.character(x)) %in% c("true", "1", "yes", "y", "on")
}

as_num <- function(x, default) {
  if (is.null(x) || length(x) == 0 || is.na(suppressWarnings(as.numeric(x)))) {
    return(default)
  }
  as.numeric(x)
}

as_chr_vec <- function(x, default = character(0)) {
  if (is.null(x)) return(default)
  as.character(unlist(x, use.names = FALSE))
}

sort_chromosomes <- function(x) {
  x <- unique(as.character(x))
  if (!length(x)) return(x)
  n <- suppressWarnings(as.numeric(x))
  if (all(!is.na(n))) x[order(n)] else sort(x)
}

empty_df <- function(columns) {
  out <- as.data.frame(setNames(vector("list", length(columns)), columns),
                       stringsAsFactors = FALSE)
  for (nm in columns) {
    if (nm %in% c("Position", "RecombinationFraction", "MapDistance_cM",
                  "Original_Length", "Scaled_Length")) {
      out[[nm]] <- numeric(0)
    } else {
      out[[nm]] <- character(0)
    }
  }
  out
}

# -----------------------------------------------------------------------------
# Desktop-facing analysis function
# -----------------------------------------------------------------------------

run_linkage <- function(config) {
  input_path <- config$input_path %||% config$filepath
  if (is.null(input_path) || !nzchar(input_path)) {
    stop("config.input_path is required.")
  }
  if (!file.exists(input_path)) stop("Input file not found: ", input_path)

  population_type <- as.character(config$population_type %||% "F2")
  if (!population_type %in% c("F2", "RIL")) {
    stop("population_type must be 'F2' or 'RIL'.")
  }

  df <- load_and_convert_data(input_path)

  result <- create_linkage_map(
    df                      = df,
    population_type         = population_type,
    missing_threshold       = as_num(config$missing_threshold, 0.30),
    recombination_threshold = as_num(config$recombination_threshold, 0.45),
    max_total_cm            = as_num(config$max_total_cm, 1500),
    scaling_method          = as.character(config$scaling_method %||% "proportional"),
    enable_chi_square       = as_bool(config$enable_chi_square, FALSE),
    chi_square_alpha        = as_num(config$chi_square_alpha, 0.05),
    use_bonferroni           = as_bool(config$use_bonferroni, FALSE)
  )

  output_dir <- as.character(config$output_dir %||%
                             file.path(dirname(normalizePath(input_path)),
                                       "linkage_results"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  rendered <- render_linkage_plots(result, config, output_dir)
  written  <- write_linkage_outputs(result, output_dir)

  lm <- result$linkage_map
  rec <- result$recombination_data

  chromosomes <- sort_chromosomes(lm$Chromosome)
  if (!length(chromosomes)) {
    chromosomes <- sort_chromosomes(df$Chromosome)
  }

  scaling_summary <- list(
    method = as.character(config$scaling_method %||% "proportional"),
    max_total_cm = as_num(config$max_total_cm, 1500),
    chromosomes_before = result$chrom_lengths_before,
    chromosomes_after  = result$chrom_lengths_after
  )

  # The cytogenetic renderer expects retained genotype rows.
  cyto_chroms <- as_chr_vec(config$selected_chromosomes, chromosomes)
  cyto_data <- build_cyto_plot_data(
    result$linkage_map,
    result$original_data[result$original_data$Marker %in%
                           result$retained_markers_list$Marker, , drop = FALSE],
    cyto_chroms
  )

  list(
    ok = TRUE,
    linkage_map = result$linkage_map,
    linkage_map_columns = colnames(result$linkage_map),

    removal_log = result$removal_log,
    removal_log_columns = colnames(result$removal_log),

    retained_markers = result$retained_markers_list,
    retained_markers_columns = colnames(result$retained_markers_list),

    all_markers = result$all_markers_positions,
    all_markers_columns = colnames(result$all_markers_positions),

    recombination_data = rec,
    recombination_data_columns = if (ncol(rec)) colnames(rec) else character(0),

    chromosomes = chromosomes,
    qc_report = result$qc_report,

    scaling_summary = scaling_summary,
    chrom_lengths_before = result$chrom_lengths_before,
    chrom_lengths_after = result$chrom_lengths_after,

    cyto_data = cyto_data,

    map_plot = rendered$map_plot,
    density_plot = rendered$density_plot,
    heatmap_plot = rendered$heatmap_plot,

    result_xlsx = written$result_xlsx,
    all_markers_xlsx = written$all_markers_xlsx,
    retained_markers_xlsx = written$retained_markers_xlsx,
    removal_log_xlsx = written$removal_log_xlsx
  )
}

# -----------------------------------------------------------------------------
# Re-rendering from an already-computed result
# -----------------------------------------------------------------------------

render_linkage_plots <- function(result, config, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  chromosomes <- sort_chromosomes(result$linkage_map$Chromosome)
  selected <- as_chr_vec(config$selected_chromosomes, chromosomes)
  if (!length(selected)) selected <- chromosomes

  show_labels <- as_bool(config$show_labels, FALSE)
  label_size <- as_num(config$label_size, 2.5)
  color_palette <- as.character(config$color_palette %||% "Set2")
  chrom_width <- as_num(config$chrom_width, 0.4)
  point_size <- as_num(config$point_size, 2.5)
  chrom_spacing <- as_num(config$chrom_spacing, 2.5)
  plot_theme <- as.character(config$plot_theme %||% "dark")
  bin_width <- as_num(config$bin_width, 10)
  top_n <- as.integer(as_num(config$heatmap_top_n, 60))

  map_png <- file.path(output_dir, "linkage_map.png")
  density_png <- file.path(output_dir, "marker_density.png")
  heatmap_png <- file.path(output_dir, "recombination_heatmap.png")

  p_map <- create_ggplot_linkage_map(
    result$linkage_map, selected,
    show_labels = show_labels,
    label_size = label_size,
    color_palette = color_palette,
    chrom_width = chrom_width,
    point_size = point_size,
    chrom_spacing = chrom_spacing,
    plot_theme = plot_theme
  )

  p_density <- create_marker_density_plot(
    result$linkage_map, selected, bin_width = bin_width
  )

  p_heat <- create_recombination_heatmap(
    result$recombination_data, selected, top_n = top_n
  )

  ggsave(map_png, p_map, width = as_num(config$map_width, 12),
         height = as_num(config$map_height, 10), dpi = 300, bg = "white")
  ggsave(density_png, p_density, width = as_num(config$density_width, 12),
         height = as_num(config$density_height, 9), dpi = 300, bg = "white")
  ggsave(heatmap_png, p_heat, width = as_num(config$heatmap_width, 12),
         height = as_num(config$heatmap_height, 10), dpi = 300, bg = "white")

  list(
    map_plot = normalizePath(map_png, winslash = "/", mustWork = FALSE),
    density_plot = normalizePath(density_png, winslash = "/", mustWork = FALSE),
    heatmap_plot = normalizePath(heatmap_png, winslash = "/", mustWork = FALSE)
  )
}

# -----------------------------------------------------------------------------
# Excel exports
# -----------------------------------------------------------------------------

write_linkage_outputs <- function(result, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  result_xlsx <- file.path(output_dir, "linkage_results.xlsx")
  all_xlsx <- file.path(output_dir, "all_markers.xlsx")
  retained_xlsx <- file.path(output_dir, "retained_markers.xlsx")
  removal_xlsx <- file.path(output_dir, "removal_log.xlsx")

  # Full workbook
  wb <- createWorkbook()
  addWorksheet(wb, "Linkage Map")
  writeData(wb, "Linkage Map", result$linkage_map)

  addWorksheet(wb, "Recombination")
  writeData(wb, "Recombination", result$recombination_data)

  addWorksheet(wb, "Retained Markers")
  writeData(wb, "Retained Markers", result$retained_markers_list)

  addWorksheet(wb, "All Markers")
  writeData(wb, "All Markers", result$all_markers_positions)

  addWorksheet(wb, "Removal Log")
  writeData(wb, "Removal Log", result$removal_log)

  addWorksheet(wb, "QC Report")
  writeData(wb, "QC Report",
            data.frame(QC_Report = result$qc_report,
                       stringsAsFactors = FALSE))

  addWorksheet(wb, "Chromosome Lengths")
  lengths_df <- result$chrom_lengths_after
  if (nrow(lengths_df) == 0) lengths_df <- result$chrom_lengths_before
  writeData(wb, "Chromosome Lengths", lengths_df)

  saveWorkbook(wb, result_xlsx, overwrite = TRUE)

  wb_all <- createWorkbook()
  addWorksheet(wb_all, "All Markers")
  writeData(wb_all, "All Markers", result$all_markers_positions)
  saveWorkbook(wb_all, all_xlsx, overwrite = TRUE)

  wb_ret <- createWorkbook()
  addWorksheet(wb_ret, "Retained Markers")
  writeData(wb_ret, "Retained Markers", result$retained_markers_list)
  saveWorkbook(wb_ret, retained_xlsx, overwrite = TRUE)

  wb_rem <- createWorkbook()
  addWorksheet(wb_rem, "Removal Log")
  writeData(wb_rem, "Removal Log", result$removal_log)
  saveWorkbook(wb_rem, removal_xlsx, overwrite = TRUE)

  list(
    result_xlsx = normalizePath(result_xlsx, winslash = "/", mustWork = FALSE),
    all_markers_xlsx = normalizePath(all_xlsx, winslash = "/", mustWork = FALSE),
    retained_markers_xlsx = normalizePath(retained_xlsx, winslash = "/", mustWork = FALSE),
    removal_log_xlsx = normalizePath(removal_xlsx, winslash = "/", mustWork = FALSE)
  )
}

# -----------------------------------------------------------------------------
# JSON serialization / CLI
# -----------------------------------------------------------------------------

sanitize_for_json <- function(x) {
  if (is.data.frame(x)) {
    # Preserve empty data-frame column schemas for the Python bridge.
    return(unname(x))
  }
  if (is.list(x)) return(lapply(x, sanitize_for_json))
  x
}

run_linkage_cli <- function(json_text) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for CLI mode.")
  }

  config <- jsonlite::fromJSON(json_text, simplifyVector = TRUE)

  ans <- tryCatch(
    {
      res <- run_linkage(config)
      jsonlite::toJSON(sanitize_for_json(res),
                       dataframe = "rows",
                       auto_unbox = TRUE,
                       na = "null",
                       null = "null",
                       pretty = FALSE)
    },
    error = function(e) {
      jsonlite::toJSON(
        list(ok = FALSE, error = conditionMessage(e)),
        auto_unbox = TRUE,
        na = "null"
      )
    }
  )
  cat(ans, "\n")
  invisible(ans)
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    cat(jsonlite::toJSON(
      list(ok = FALSE, error = "Missing JSON configuration argument."),
      auto_unbox = TRUE
    ), "\n")
    quit(status = 2)
  }
  run_linkage_cli(args[[1]])
}

