# =============================================================================
# 4. ANALYSIS FUNCTIONS (each returns NULL when the data cannot support it)
# =============================================================================

# Mean of each trait per genotype: mean of the genotype x environment means, so every
# environment carries equal weight. Used by correlation, PCA and clustering.
genotype_means <- function(d, traits) {
  d %>%
    mutate(Genotype = as.character(Genotype), Environment = as.character(Environment)) %>%
    group_by(Genotype, Environment) %>%
    summarise(across(all_of(traits), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    group_by(Genotype) %>%
    summarise(across(all_of(traits), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    mutate(across(all_of(traits), ~ ifelse(is.nan(.x), NA_real_, .x)))
}

# Descriptive statistics of `Value` by one grouping column (Environment or Trait).
# Kurtosis is EXCESS kurtosis (normal distribution = 0).
describe_values <- function(df, by, quartiles = FALSE) {
  res <- df %>% group_by(across(all_of(by))) %>%
    summarise(N = sum(!is.na(Value)), Mean = round(mean(Value, na.rm = TRUE), 3),
              Median = round(median(Value, na.rm = TRUE), 3), SD = round(sd(Value, na.rm = TRUE), 3),
              CV = round(sd(Value, na.rm = TRUE) / mean(Value, na.rm = TRUE) * 100, 2),
              Min = round(suppressWarnings(min(Value, na.rm = TRUE)), 3),
              Max = round(suppressWarnings(max(Value, na.rm = TRUE)), 3),
              Q1 = round(quantile(Value, .25, na.rm = TRUE, names = FALSE), 3),
              Q3 = round(quantile(Value, .75, na.rm = TRUE, names = FALSE), 3),
              Skewness = round(moments::skewness(Value, na.rm = TRUE), 3),
              Kurtosis = round(moments::kurtosis(Value, na.rm = TRUE) - 3, 3),
              Missing = sum(is.na(Value)), .groups = "drop")
  if (!quartiles) res <- select(res, -Q1, -Q3)
  res[] <- lapply(res, function(x) if (is.numeric(x)) replace(x, !is.finite(x), NA) else x)   # no data -> NA, not Inf
  res
}

trait_long <- function(d, trait) {
  d %>% select(Environment, Genotype, all_of(trait)) %>% rename(Value = !!trait)
}

env_means_table <- function(wd, trait) {
  trait_long(wd, trait) %>% filter(!is.na(Value)) %>%
    group_by(Genotype, Environment) %>%
    summarise(Mean = round(mean(Value), 3), SD = round(sd(Value), 3), N = n(), .groups = "drop")
}

env_rank_tables <- function(env_means) {
  env_means %>% group_by(Environment) %>%
    mutate(Rank = rank(-Mean, ties.method = "min")) %>% ungroup() %>%
    select(Genotype, Environment, Mean, Rank)
}

# Pearson correlation between environments, on genotype means
env_correlation <- function(wd, trait) {
  wide <- trait_long(wd, trait) %>% filter(!is.na(Value)) %>%
    group_by(Genotype, Environment) %>% summarise(Mean = mean(Value), .groups = "drop") %>%
    pivot_wider(names_from = Environment, values_from = Mean) %>% column_to_rownames("Genotype")
  if (ncol(wide) < 2) return(NULL)
  cm <- cor(wide, use = "pairwise.complete.obs")
  # helper names that cannot clash with an environment called "Env1" / "Env2"
  tab <- cm %>% as.data.frame() %>% rownames_to_column(".env_a") %>%
    pivot_longer(-.env_a, names_to = ".env_b", values_to = "Correlation") %>%
    filter(.env_a != .env_b) %>% mutate(Correlation = round(Correlation, 3)) %>%
    rename(Env1 = .env_a, Env2 = .env_b)
  list(matrix = cm, table = tab)
}

# ANOVA of one trait in ONE environment (data `d_env`). Augmented RCBD: Federer's classical table;
# otherwise a sequential aov. Always the columns Part, term, df, sumsq, meansq, statistic, p.value.
anova_one_env <- function(d_env, trait, design) {
  df <- model_df(d_env, trait, design)
  if (nrow(df) < 6 || length(unique(df$Genotype)) < 2) return(NULL)
  tryCatch({
    prep <- if (design == "augmented") augmented_prepare(df)
    if (!is.null(prep)) {
      augmented_anova(df, prep)
    } else {
      m <- aov(as.formula(design_formula(design, "aov", df = df)), data = df)
      sm <- broom::tidy(m)
      sm$Part <- if (design == "augmented") "Sequential (Genotype first, blocks not estimable from checks)" else "Sequential ANOVA"
      sm[, c("Part", "term", "df", "sumsq", "meansq", "statistic", "p.value")]
    }
  }, error = function(e) NULL)
}

# Adjusted means of a mixed model: the estimates and standard errors are exactly those of emmeans(m, specs).
# emmeans spends nearly all its time on Kenward-Roger degrees of freedom, which are never used here, so the
# Kenward-Roger covariance matrix is computed once and handed to emmeans (same estimates, same SEs, ~40x faster).
# Falls back to plain emmeans() when pbkrtest is missing or the data are above emmeans' own size limit for it.
lmer_emmeans <- function(m, specs) {
  V <- NULL
  if (requireNamespace("pbkrtest", quietly = TRUE) && nrow(stats::model.frame(m)) <= emmeans::get_emm_option("pbkrtest.limit"))
    V <- tryCatch(as.matrix(pbkrtest::vcovAdj(m)), error = function(e) NULL)
  em <- if (is.null(V)) emmeans(m, specs) else emmeans(m, specs, lmer.df = "asymptotic", vcov. = V)
  as.data.frame(em)
}

# Within-environment BLUEs. Augmented RCBD: fixed blocks (classical analysis); otherwise a mixed
# model with random blocks. Returns the table plus the messages produced on the way.
env_blues_table <- function(wd, trait, design) {
  notes <- character(0); warns <- character(0); errs <- list(); sed_list <- list()
  res <- lapply(unique(as.character(wd$Environment)), function(env) {
    df <- model_df(wd[as.character(wd$Environment) == env, , drop = FALSE], trait, design)
    if (nrow(df) < 6 || length(unique(df$Genotype)) < 2) return(NULL)
    tryCatch({
      prep <- if (design == "augmented") augmented_prepare(df)
      if (!is.null(prep)) {
        bn <- augmented_block_note(df, env)
        if (!is.null(bn)) notes <<- c(notes, bn)
        am <- augmented_means(df, prep)
        sdt <- am$sed; sdt$Environment <- env; sed_list[[env]] <<- sdt
        em <- am$blues
      } else {
        if (design == "augmented")
          warns <<- c(warns, sprintf(
            "%s: genotypes cannot be separated from blocks (no check is repeated across the blocks, or every block holds the same genotypes), so the classical fixed-block analysis is not possible; blocks were treated as random instead.", env))
        m <- lmer(as.formula(design_formula(design, "lmer_geno", df = df)), data = df)
        em <- lmer_emmeans(m, ~Genotype)
      }
      em$Environment <- env
      em$Genotype <- as.character(em$Genotype)
      em[, c("Environment", "Genotype", "emmean", "SE")] %>% rename(BLUE = emmean)
    }, error = function(e) { errs[[paste0("env_blues_", env)]] <<- conditionMessage(e); NULL })
  })
  res <- Filter(Negate(is.null), res)
  sed <- NULL
  if (length(sed_list)) {
    sed <- dplyr::bind_rows(sed_list)
    sed <- sed[, c("Environment", "Comparison", "Pairs", "SED_mean", "SED_min", "SED_max")]
  }
  list(blues = if (length(res)) dplyr::bind_rows(res) else NULL, sed = sed,
       notes = notes, warnings = warns, errors = errs)
}

# Pooled ANOVA across environments. Blocks are nested in environments. For RCBD and alpha-lattice
# Genotype x Environment is a separate term (so Genotype is tested against pure error) and the table
# also gives F_vs_GxE (Genotype MS / G x E MS), the test to use when environments are a sample of
# target environments. Augmented RCBD keeps G x E in the residual. Always a true ANOVA table.
pooled_anova <- function(wd, trait, design) {
  td <- model_df(wd, trait, design)
  if (nrow(td) < 10) return(NULL)
  fit <- function(gxe) lm(as.formula(design_formula(design, "lm_pooled", interaction = gxe, df = td)), data = td)
  sm <- NULL; used_gxe <- FALSE
  single_env <- nlevels(td$Environment) < 2     # trait measured in one season only: plain one-environment model
  # R labels the interaction in order of appearance ("Environment:Genotype"), so match either spelling
  is_gxe <- function(x) x %in% c("Genotype:Environment", "Environment:Genotype")
  # Augmented RCBD: new entries have one plot per environment, so Genotype x Environment would
  # absorb every plot of a new entry (residual df only from the checks) and is left in the residual.
  m_gxe <- if (design == "augmented" || single_env) NULL else tryCatch(fit(TRUE), error = function(e) NULL)
  a <- NULL
  if (!is.null(m_gxe) && df.residual(m_gxe) > 0) {
    a_gxe <- anova(m_gxe)
    if (any(is_gxe(rownames(a_gxe)) & a_gxe$Df > 0)) { sm <- m_gxe; used_gxe <- TRUE; a <- a_gxe }
  }
  if (is.null(sm)) {
    sm <- if (single_env) lm(as.formula(design_formula(design, "aov", df = td)), data = td) else fit(FALSE)
  }
  if (is.null(a)) a <- anova(sm)

  tbl <- broom::tidy(a)
  tbl$term[is_gxe(tbl$term)] <- "Genotype:Environment"
  tbl <- tbl[tbl$df > 0 | tbl$term == "Residuals", ]
  tbl$F_vs_GxE <- NA_real_; tbl$p_vs_GxE <- NA_real_
  if (used_gxe && design != "augmented") {
    g <- tbl[tbl$term == "Genotype", ]; ge <- tbl[tbl$term == "Genotype:Environment", ]
    if (nrow(g) == 1 && nrow(ge) == 1) {
      f <- g$meansq / ge$meansq
      tbl$F_vs_GxE[tbl$term == "Genotype"] <- f
      tbl$p_vs_GxE[tbl$term == "Genotype"] <- pf(f, g$df, ge$df, lower.tail = FALSE)
    }
  }
  list(model = sm, data = td, table = tbl, used_gxe = used_gxe)
}

# Residuals of a pooled model (only needed for the diagnostic plots of one trait, so they are computed on request)
pooled_residuals <- function(pa) {
  rd <- data.frame(Fitted = fitted(pa$model), Residuals = residuals(pa$model),
                   Standardized_Residuals = rstandard(pa$model), Genotype = pa$data$Genotype)
  rd[is.finite(rd$Standardized_Residuals), ]    # plots that are fitted exactly (e.g. unreplicated new entries) have none
}

# ── Single-trait ANOVA: every trait, in every environment (season) and pooled ─
# Maps an ANOVA row to the source of variation it stands for
anova_source <- function(term) {
  t <- trimws(term)
  ifelse(grepl("^Residual", t), "Residual",
  ifelse(t %in% c("Genotype:Environment", "Environment:Genotype"), "G x E",
  ifelse(t %in% c("Genotype", "Treatment (adjusted for blocks)"), "Genotype",
  ifelse(t == "Environment", "Environment", "Blocks/Reps"))))
}

# The rows of an ANOVA table that add up to the total sum of squares (no sub-partitions, no "ignoring blocks" half)
anova_main_rows <- function(tbl) {
  keep <- !startsWith(tbl$term, " ")
  if ("Part" %in% names(tbl)) keep <- keep & (is.na(tbl$Part) | tbl$Part != "Treatments ignoring blocks")
  tbl[keep, , drop = FALSE]
}

# Plain-language names of the ANOVA sources, and significance stars
FRIENDLY_SOURCES <- c(
  "Block (ignoring treatments)" = "Blocks",
  "Treatment (adjusted for blocks)" = "Genotypes (adjusted for blocks)",
  "Among checks" = "   Checks",
  "New entries, and new entries vs checks" = "   New entries (and entries vs checks)",
  "Residuals" = "Error",
  "Rep" = "Replications", "Block" = "Blocks", "RepBlock" = "Blocks (within replications)",
  "EnvRep" = "Replications (within environments)", "EnvBlock" = "Blocks (within environments)",
  "EnvRepBlock" = "Blocks (within replications and environments)",
  "Genotype" = "Genotypes", "Environment" = "Environments", "Genotype:Environment" = "Genotype × Environment")

friendly_source <- function(term) {
  t <- trimws(term)
  nm <- unname(ifelse(t %in% names(FRIENDLY_SOURCES), FRIENDLY_SOURCES[t], t))
  ifelse(startsWith(nm, " "), nm, paste0(sub("^(\\s*).*$", "\\1", term), nm))   # sub-rows keep their indentation
}

sig_stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

# A classic ANOVA table (Source, DF, SumSq, MeanSq, F, P, Sig, with a Total row) from the raw rows of
# one trait in one environment. Only the block-adjusted half of the augmented table is shown.
anova_classic <- function(tb) {
  if ("Part" %in% names(tb)) tb <- tb[is.na(tb$Part) | tb$Part != "Treatments ignoring blocks", , drop = FALSE]
  main <- !startsWith(tb$term, " ")
  out <- data.frame(Source = friendly_source(tb$term), DF = tb$df, `Sum Sq` = round(tb$sumsq, 3),
                    `Mean Sq` = round(tb$meansq, 3), F = round(tb$statistic, 3), P = signif(tb$p.value, 4),
                    `Sig.` = sig_stars(tb$p.value), stringsAsFactors = FALSE, check.names = FALSE)
  out <- rbind(out, data.frame(Source = "Total", DF = sum(tb$df[main]), `Sum Sq` = round(sum(tb$sumsq[main]), 3),
                               `Mean Sq` = NA_real_, F = NA_real_, P = NA_real_, `Sig.` = "", stringsAsFactors = FALSE,
                               check.names = FALSE))
  # pooled table of an RCBD / alpha-lattice: Genotype is also tested against the G x E mean square
  if (all(c("F_vs_GxE", "p_vs_GxE") %in% names(tb)) && any(!is.na(tb$F_vs_GxE))) {
    out[[paste("F vs", "G×E")]] <- c(round(tb$F_vs_GxE, 3), NA)
    out[[paste("P vs", "G×E")]] <- c(signif(tb$p_vs_GxE, 4), NA)
  }
  out
}

# ── One consolidated ANOVA table per trait ─────────────────────────────────────
# Rows are the sources of variation; each environment (season) and the pooled analysis is a column group
# (DF, SS, MS, F, P, Sig.), so seasons can be compared side by side. Below the ANOVA come plots, mean
# and CV of each group. Flat column names are "<group>: <column>" (e.g. "Env1: DF").
ANOVA_ROW_NAMES <- c("Blocks (within environments)" = "Blocks", "Replications (within environments)" = "Replications",
                     "Blocks (within replications and environments)" = "Blocks (within replications)")
ANOVA_FOOTER <- c("Plots", "Mean", "CV (%)")
ANOVA_COLS <- c(DF = "DF", SS = "Sum Sq", MS = "Mean Sq", F = "F", P = "P", Sig = "Sig.")

# an: long ANOVA of one trait (Environment, Source, DF, Sum Sq, ...), sm: its summary rows (N, Mean, CV_pct)
anova_wide <- function(an, sm, env_order = unique(an$Environment)) {
  an <- as.data.frame(an, check.names = FALSE); an$Source <- as.character(an$Source)
  an$Source <- ifelse(an$Source %in% names(ANOVA_ROW_NAMES), unname(ANOVA_ROW_NAMES[an$Source]), an$Source)
  # pooled Genotype (entered after the blocks) is the same test as "adjusted for blocks" of the augmented design
  if (any(an$Source == "Genotypes (adjusted for blocks)"))
    an$Source[an$Environment == "Pooled" & an$Source == "Genotypes"] <- "Genotypes (adjusted for blocks)"
  groups <- c(intersect(setdiff(env_order, "Pooled"), an$Environment), intersect("Pooled", an$Environment))
  canon <- c("Environments", "Replications", "Blocks", "Blocks (within replications)", "Genotypes",
             "Genotypes (adjusted for blocks)", "   Checks", "   New entries (and entries vs checks)", "Genotype × Environment")
  src <- unique(an$Source)
  rows <- c(intersect(canon, src), setdiff(src, c(canon, "Error", "Total")), intersect("Error", src), "Total")
  out <- data.frame(Source = rows, stringsAsFactors = FALSE)
  foot <- data.frame(Source = ANOVA_FOOTER, stringsAsFactors = FALSE)
  for (g in groups) {
    s <- an[an$Environment == g, , drop = FALSE]; i <- match(rows, s$Source)
    spec <- ANOVA_COLS
    if ("F vs G×E" %in% names(s) && any(!is.na(s[["F vs G×E"]]))) spec <- c(spec, `F vs G×E` = "F vs G×E", `P vs G×E` = "P vs G×E")
    for (lab in names(spec)) {
      nm <- paste0(g, ": ", if (lab == "Sig") "Sig." else lab)
      out[[nm]] <- s[[spec[[lab]]]][i]
      foot[[nm]] <- NA
    }
    first <- paste0(g, ": DF"); m <- sm[sm$Environment == g, ][1, ]
    foot[[first]] <- c(m$N, m$Mean, m$CV_pct)
  }
  cols_of <- lapply(groups, function(g) names(out)[startsWith(names(out), paste0(g, ": "))]); names(cols_of) <- groups
  list(table = rbind(out, foot), groups = cols_of)
}

# All traits in one CSV-ready table: Trait, Source, then the column groups in environment order, Pooled last
anova_wide_all <- function(wide, env_order) {
  tabs <- lapply(names(wide), function(tr) cbind(Trait = tr, wide[[tr]]$table, stringsAsFactors = FALSE))
  all <- dplyr::bind_rows(tabs)
  keep <- setdiff(names(all), c("Trait", "Source"))
  g <- sub(": .*$", "", keep)
  gl <- c(setdiff(env_order, "Pooled"), "Pooled")
  all[, c("Trait", "Source", keep[order(match(g, gl), seq_along(keep))])]
}

# The same results as a tidy table for download: one row per trait x environment x source, nothing but the cells that
# exist (no empty season columns, no NA), and the plots, mean and CV of each analysis on its Total row.
anova_long_all <- function(wide) {
  nf <- length(ANOVA_FOOTER)
  parts <- unlist(lapply(names(wide), function(tr) {
    tb <- wide[[tr]]$table; grp <- wide[[tr]]$groups; nd <- nrow(tb) - nf
    lapply(names(grp), function(g) {
      cols <- grp[[g]]; lab <- substring(cols, nchar(g) + 3)
      d <- tb[seq_len(nd), cols, drop = FALSE]; names(d) <- lab
      keep <- rowSums(!is.na(d[setdiff(lab, "Sig.")])) > 0           # sources that belong to this analysis
      d <- d[keep, , drop = FALSE]
      src <- tb$Source[seq_len(nd)][keep]
      src <- ifelse(startsWith(src, " "), paste("of which:", trimws(src)), src)
      foot <- tb[nd + seq_len(nf), cols[1]]                          # plots, mean, CV
      res <- data.frame(Trait = tr, Environment = if (g == "Pooled") "Pooled over environments" else g, Source = src,
                        d, check.names = FALSE, stringsAsFactors = FALSE)
      tot <- src == "Total"
      res$Plots <- NA_real_; res$Mean <- NA_real_; res[["CV (%)"]] <- NA_real_
      res$Plots[tot] <- foot[1]; res$Mean[tot] <- foot[2]; res[["CV (%)"]][tot] <- foot[3]
      res
    })
  }), recursive = FALSE)
  all <- dplyr::bind_rows(parts)
  for (pc in intersect(c("P", "P vs G×E"), names(all)))          # re-parsed from 4 significant digits: no 4.109999999999999e-94 in the file
    all[[pc]] <- as.numeric(sprintf("%.4g", all[[pc]]))
  lead <- c("Trait", "Environment", "Source", "DF", "SS", "MS", "F", "P", "Sig.")
  all[, c(lead, setdiff(names(all), c(lead, "Plots", "Mean", "CV (%)")), "Plots", "Mean", "CV (%)")]
}

# The wide ANOVA table as the app shows it, rebuilt from the state of a run folder (NULL when no ANOVA was run)
anova_wide_table <- function(outdir) {
  f <- file.path(outdir, "state.rds")
  w <- if (file.exists(f)) readRDS(f)$ta$wide
  if (is.null(w)) NULL else anova_wide_all(w, unique(unlist(lapply(w, function(x) names(x$groups)))))
}
# Returns list(anova, summary, partition, pooled):
#   anova     classic ANOVA tables: every trait in every environment (season) and pooled over environments
#   summary   one line per trait x environment: genotype F, P and stars, error df, CV%, mean squares
#             (heritability is added by the pipeline from heritability_all)
#   partition share of the total sum of squares by source of variation (for the stacked bar plot)
trait_anova_all <- function(wd, traits, design, keep_pooled = NULL) {
  envs <- unique(as.character(wd$Environment))
  tabs <- list(); means <- list(); pooled <- list()
  add <- function(tbl, tr, env, d) {
    tbl$Trait <- tr; tbl$Environment <- env
    tabs[[length(tabs) + 1]] <<- tbl
    v <- d[[tr]]
    means[[length(means) + 1]] <<- data.frame(Trait = tr, Environment = env, N = sum(!is.na(v)), Mean = mean(v, na.rm = TRUE))
  }
  for (tr in traits) {
    for (env in envs) {
      d <- wd[as.character(wd$Environment) == env, , drop = FALSE]
      sm <- anova_one_env(d, tr, design)
      if (!is.null(sm)) add(sm, tr, env, d)
    }
    d_tr <- wd[!is.na(wd[[tr]]), , drop = FALSE]
    if (n_distinct(d_tr$Environment) > 1) {
      pa <- tryCatch(pooled_anova(d_tr, tr, design), error = function(e) NULL)
      if (!is.null(pa)) { if (identical(tr, keep_pooled)) pooled[[tr]] <- pa; tb <- as.data.frame(pa$table); tb$Part <- "Pooled over environments"; add(tb, tr, "Pooled", d_tr) }
    }
  }
  if (!length(tabs)) return(NULL)
  anova_tbl <- dplyr::bind_rows(tabs)
  for (cc in c("F_vs_GxE", "p_vs_GxE")) if (!cc %in% names(anova_tbl)) anova_tbl[[cc]] <- NA_real_
  mean_tbl <- dplyr::bind_rows(means)
  keys <- unique(anova_tbl[, c("Trait", "Environment")])

  per_key <- lapply(seq_len(nrow(keys)), function(i) {
    tr <- keys$Trait[i]; env <- keys$Environment[i]
    tb <- anova_tbl[anova_tbl$Trait == tr & anova_tbl$Environment == env, ]
    mn <- anova_main_rows(tb); mn$Source <- anova_source(mn$term)
    pick <- function(s) mn[mn$Source == s, , drop = FALSE][1, ]      # an all-NA row when the source is absent
    g <- pick("Genotype"); e <- pick("Residual"); ge <- pick("G x E")
    mm <- mean_tbl[mean_tbl$Trait == tr & mean_tbl$Environment == env, ][1, ]
    ss <- tapply(mn$sumsq, mn$Source, sum)
    list(anova = cbind(Trait = tr, Environment = env, anova_classic(tb)),
         summary = data.frame(
           Trait = tr, Environment = env, N = mm$N, Mean = signif(mm$Mean, 5), DF_Error = e$df,
           F_Genotype = round(g$statistic, 3), P_Genotype = signif(g$p.value, 4), Sig = sig_stars(g$p.value),
           CV_pct = round(sqrt(e$meansq) / mm$Mean * 100, 2),
           MS_Genotype = signif(g$meansq, 5), MS_Error = signif(e$meansq, 5),
           MS_GxE = signif(ge$meansq, 5), P_Genotype_vs_GxE = signif(g$p_vs_GxE, 4), stringsAsFactors = FALSE),
         partition = data.frame(Trait = tr, Environment = env, Source = names(ss),
                                SS = as.numeric(ss), Percent = round(as.numeric(ss) / sum(ss) * 100, 1),
                                stringsAsFactors = FALSE))
  })
  list(anova = dplyr::bind_rows(lapply(per_key, `[[`, "anova")),
       summary = dplyr::bind_rows(lapply(per_key, `[[`, "summary")),
       partition = dplyr::bind_rows(lapply(per_key, `[[`, "partition")),
       pooled = pooled)          # the pooled fit of `keep_pooled` only (for its residual plots; others are not kept in memory)
}

# BLUEs across environments (lmer, blocks nested in environments).
#   include_interaction = TRUE : Genotype | Environment BLUEs (Environment:Genotype in the model)
#   otherwise                   : ONE BLUE per genotype across environments. Without the interaction
#                                 term the model gives the same genotype differences in every
#                                 environment, so per-environment BLUEs would be meaningless.
across_env_blues <- function(wd, trait, design, include_interaction) {
  td <- model_df(wd, trait, design)
  if (nrow(td) < 10) return(NULL)
  multi <- nlevels(td$Environment) > 1
  note <- NULL
  if (!multi) {
    m <- lmer(as.formula(design_formula(design, "lmer_geno", df = td)), data = td)
    em <- lmer_emmeans(m, ~Genotype); em$Environment <- as.character(td$Environment[1])
  } else if (isTRUE(include_interaction)) {
    m <- lmer(as.formula(design_formula(design, "lmer_env", across_env = TRUE, interaction = TRUE, df = td)), data = td)
    em <- lmer_emmeans(m, ~ Genotype | Environment)
  } else {
    m <- lmer(as.formula(design_formula(design, "lmer_env", across_env = TRUE, interaction = FALSE, df = td)), data = td)
    em <- lmer_emmeans(m, ~Genotype); em$Environment <- "All environments"
    note <- "Mixed-model BLUEs are across all environments (one value per genotype), because the Genotype x Environment interaction is switched off. Tick 'include interaction' for BLUEs per environment."
  }
  em <- em[!is.na(em$emmean), ]
  em$Genotype <- as.character(em$Genotype); em$Environment <- as.character(em$Environment)
  list(table = em[, c("Genotype", "Environment", "emmean", "SE")] %>% rename(BLUE = emmean), note = note)
}

# ── Variance components and heritability, per trait and per environment (season) ──
# Genotypes that are checks of an augmented design: repeated within an environment
check_genotypes <- function(td) {
  cnt <- table(as.character(td$Genotype), as.character(td$Environment))
  rownames(cnt)[rowSums(cnt > 1) > 0]
}

# One mixed model for one trait.
#   scope "pooled": all environments together   Value ~ Environment [* Check] + (1|Genotype) + (1|Genotype:Environment) + (1|blocks)
#   scope "single": one environment             Value ~ [Check] + (1|Genotype) + (1|blocks)
# Augmented RCBD: every check gets its own FIXED effect (Check). If the checks were random genotypes, the
# differences between check varieties would be counted as genotypic variance and heritability would be too
# high (in a simulation with known Vg = 16 this gave 21.9, against 16.3 with fixed checks).
# Entry-mean heritability uses the harmonic mean number of environments (ne) and plots (np) per ENTRY:
#   pooled  H2 = Vg / (Vg + Vge/ne + Ve/np)        single  H2 = Vg / (Vg + Ve/np)
herit_model <- function(td, design, scope) {
  td <- droplevels(td)
  chk <- if (design == "augmented") check_genotypes(td) else character(0)
  td$Check <- factor(ifelse(as.character(td$Genotype) %in% chk, as.character(td$Genotype), "entry"))
  has_chk <- length(chk) > 0 && nlevels(td$Check) > 1
  fixed <- if (scope == "pooled") c("Environment", if (has_chk) c("Check", "Environment:Check")) else if (has_chk) "Check"
  kind <- if (scope == "pooled") "lmer_var" else "lmer_var_single"
  fml <- design_formula(design, kind, across_env = (scope == "pooled"), df = td, fixed = fixed)
  m <- suppressMessages(lmer(as.formula(fml), data = td))
  list(m = m, td = td, chk = chk, has_chk = has_chk)
}

herit_one <- function(td, design, scope) {
  hm <- herit_model(td, design, scope)
  m <- hm$m; td <- hm$td; chk <- hm$chk; has_chk <- hm$has_chk
  vc <- as.data.frame(VarCorr(m))
  get_vc <- function(g) { v <- vc$vcov[vc$grp == g]; if (length(v)) v[1] else 0 }
  vg <- get_vc("Genotype"); vge <- get_vc("Genotype:Environment"); ve <- get_vc("Residual")
  vblk <- sum(vc$vcov[!vc$grp %in% c("Genotype", "Genotype:Environment", "Residual")])

  ent <- td[!(as.character(td$Genotype) %in% chk), , drop = FALSE]
  per_g <- ent %>% group_by(Genotype) %>% summarise(np = n(), ne = n_distinct(Environment), .groups = "drop")
  np <- harmonic_mean(per_g$np); ne <- harmonic_mean(per_g$ne)
  h2_plot <- vg / (vg + vge + ve)
  h2_mean <- if (scope == "pooled") vg / (vg + vge / ne + ve / np) else vg / (vg + ve / np)

  note <- character(0)
  if (vg < 1e-8) note <- c(note, "No genotypic variance detected (H² = 0)")
  if (design == "augmented" && !has_chk) note <- c(note, "No checks found: all genotypes treated as random")
  if (design == "augmented" && has_chk)
    note <- c(note, sprintf("Error from %d check plots", sum(as.character(td$Genotype) %in% chk)))
  data.frame(N_Genotypes = nlevels(ent$Genotype), Vg = vg, Vge = if (scope == "pooled") vge else NA_real_,
                     Vblk = vblk, Ve = ve, H2_Plot = h2_plot, H2_Mean = h2_mean,
                     Note = paste(note, collapse = "; "), stringsAsFactors = FALSE)
}

# Every trait: pooled over environments (needs >= 2) and in each environment.
# Returns one row per Trait x Environment ("Pooled" = all environments together), or NULL.
heritability_all <- function(wd, traits, design) {
  rows <- list()
  for (tr in traits) {
    d_tr <- wd[!is.na(wd[[tr]]), , drop = FALSE]
    envs <- unique(as.character(d_tr$Environment))
    jobs <- c(if (length(envs) > 1) list(list(env = "Pooled", d = d_tr, scope = "pooled")),
              lapply(envs, function(e) list(env = e, d = d_tr[as.character(d_tr$Environment) == e, , drop = FALSE], scope = "single")))
    for (j in jobs) {
      td <- model_df(j$d, tr, design)
      if (nrow(td) < 10 || nlevels(td$Genotype) < 5) next
      r <- tryCatch(herit_one(td, design, j$scope), error = function(e) NULL)
      if (is.null(r)) next
      rows[[length(rows) + 1]] <- cbind(Trait = tr, Environment = j$env, r)
    }
  }
  if (length(rows) == 0) return(NULL)
  h2 <- dplyr::bind_rows(rows)
  num <- c("Vg", "Vge", "Vblk", "Ve")
  h2[num] <- lapply(h2[num], signif, 4)
  h2[c("H2_Plot", "H2_Mean")] <- lapply(h2[c("H2_Plot", "H2_Mean")], round, 3)
  h2
}

# ── BLUPs: best linear unbiased predictions of the accessions ────────────────────────────────
# The model is the one of the heritability analysis (herit_model): the checks of an augmented design are fixed effects
# (they get BLUEs, not BLUPs), accessions, Genotype x Environment and blocks are random, REML. The BLUPs are the
# accession effects of Henderson's mixed model equations, which are solved here directly with the variances at their
# REML estimates, because the same equations give the prediction error variance (PEV) of every BLUP INCLUDING the
# uncertainty of the fixed effects; the condVar of lme4 leaves that part out and so reports a PEV that is too small.
#   PEV_i        = Ve * [C^-1]_ii            SE_i = sqrt(PEV_i)          Reliability_i = 1 - PEV_i / Vg
#   H2 Cullis    = 1 - vbar / (2 Vg)         vbar = mean variance of a DIFFERENCE between two BLUPs (Cullis et al. 2006)
#                = 1 - (n tr(P) - 1'P1) / (n (n - 1) Vg),   P = prediction error (co)variance matrix of the n BLUPs
# 1 - mean(PEV)/Vg is the mean reliability, a different (slightly higher) number; both are reported.
blup_mme <- function(m, group = "Genotype", exclude = character(0), chunk = 300) {
  X <- getME(m, "X"); y <- getME(m, "y"); s2 <- sigma(m)^2
  vc <- as.data.frame(VarCorr(m)); vc <- vc[!is.na(vc$grp) & vc$grp != "Residual" & is.na(vc$var2), , drop = FALSE]
  Zl <- getME(m, "Ztlist")
  blocks <- list(); lam <- numeric(0); is_g <- logical(0); gnames <- character(0)
  for (nm in names(Zl)) {
    g <- sub("\\.\\(Intercept\\)$", "", nm)
    v <- vc$vcov[vc$grp == g][1]
    if (is.na(v) || v < 1e-8 * s2) next                      # a term without variance: all its effects are 0
    zt <- Zl[[nm]]
    if (g == group) zt <- zt[!(rownames(zt) %in% exclude), , drop = FALSE]
    if (!nrow(zt)) next
    blocks[[length(blocks) + 1]] <- zt
    lam <- c(lam, rep(s2 / v, nrow(zt)))
    is_g <- c(is_g, rep(g == group, nrow(zt)))
    if (g == group) gnames <- rownames(zt)
  }
  if (!any(is_g)) return(NULL)
  Zt <- methods::as(do.call(rbind, blocks), "CsparseMatrix")            # random-effect design, effects x plots
  p <- ncol(X); gi <- p + which(is_g)
  ZtX <- as.matrix(Zt %*% X)
  C <- Matrix::forceSymmetric(Matrix::rbind2(
    Matrix::cbind2(Matrix::Matrix(crossprod(X), sparse = TRUE), Matrix::Matrix(t(ZtX), sparse = TRUE)),
    Matrix::cbind2(Matrix::Matrix(ZtX, sparse = TRUE), Matrix::tcrossprod(Zt) + Matrix::Diagonal(x = lam))))
  ch <- Matrix::Cholesky(C, LDL = FALSE, perm = TRUE)
  sol <- as.numeric(Matrix::solve(ch, c(as.numeric(crossprod(X, y)), as.numeric(Zt %*% y)), system = "A"))
  u <- setNames(sol[gi], gnames)
  # prediction error variances: the diagonal of the genotype block of C^-1 and the sum of all its entries, in chunks
  n <- length(gi); pev <- numeric(n); tot <- 0
  for (s in seq(1, n, by = chunk)) {
    j <- s:min(n, s + chunk - 1)
    B <- Matrix::sparseMatrix(i = gi[j], j = seq_along(j), x = 1, dims = c(nrow(C), length(j)))
    S <- as.matrix(Matrix::solve(ch, B, system = "A"))
    pev[j] <- s2 * S[cbind(gi[j], seq_along(j))]
    tot <- tot + sum(S[gi, , drop = FALSE])
  }
  names(pev) <- gnames
  list(u = u, pev = pev, trace = sum(pev), total = s2 * tot, n = n, fixed = as.numeric(X %*% lme4::fixef(m)))
}

# BLUPs of one trait in one scope: "pooled" (all environments: the across-environment value of an accession) or "single"
# (one environment). Returns list(table, summary).
blup_one <- function(td, design, scope) {
  hm <- herit_model(td, design, scope)
  m <- hm$m; chk <- hm$chk
  vc <- as.data.frame(VarCorr(m))
  get_vc <- function(g) { v <- vc$vcov[vc$grp == g]; if (length(v)) v[1] else 0 }
  vg <- get_vc("Genotype"); vge <- get_vc("Genotype:Environment"); ve <- get_vc("Residual")
  used <- hm$td[match(rownames(m@frame), rownames(hm$td)), , drop = FALSE]
  isent <- !(as.character(used$Genotype) %in% chk)
  ent <- used[isent, , drop = FALSE]
  per_g <- ent %>% group_by(Genotype) %>% summarise(np = n(), ne = n_distinct(Environment), .groups = "drop")
  np <- harmonic_mean(per_g$np); ne <- harmonic_mean(per_g$ne)
  h2_mean <- if (scope == "pooled") vg / (vg + vge / ne + ve / np) else vg / (vg + ve / np)

  note <- character(0)
  zero_vc <- vc$grp[vc$grp != "Residual" & vc$vcov < 1e-6 * max(ve, 1e-12)]
  if (length(zero_vc))
    note <- c(note, paste0("Variance estimated as 0 (singular fit): ", paste(sub("Genotype:Environment", GXE, zero_vc, fixed = TRUE), collapse = ", ")))
  if (design == "augmented" && length(chk)) note <- c(note, sprintf("%d checks are fixed effects and have no BLUP", length(chk)))
  r <- if (vg > 1e-8 * (vg + ve)) blup_mme(m, "Genotype", exclude = chk) else NULL
  tab <- NULL
  if (is.null(r)) {
    note <- c(note, "No genotypic variance detected: every BLUP is 0, heritability is 0")
    h2c <- 0; rel <- 0
  } else {
    n <- r$n
    vbar <- 2 * (n * r$trace - r$total) / (n * (n - 1))             # mean variance of a difference between two BLUPs
    h2c <- if (n > 1) max(0, 1 - vbar / (2 * vg)) else NA_real_
    rel <- 1 - mean(r$pev) / vg
    base <- mean(r$fixed[isent])            # the fixed part (environments, blocks) averaged over the plots of the accessions
    acc <- names(r$u)
    cnt <- table(as.character(ent$Genotype))[acc]
    nenv <- tapply(as.character(ent$Environment), as.character(ent$Genotype), function(z) length(unique(z)))[acc]
    rawm <- tapply(ent$Value, as.character(ent$Genotype), mean)[acc]
    tab <- data.frame(Accession = acc, Plots = as.integer(cnt), Environments = as.integer(nenv), Raw_Mean = as.numeric(rawm),
                      BLUP = as.numeric(r$u), SE = sqrt(pmax(r$pev, 0)), Reliability = 1 - r$pev / vg,
                      Predicted = base + as.numeric(r$u), stringsAsFactors = FALSE)
    tab$Rank <- rank(-tab$BLUP, ties.method = "min")
    tab <- tab[order(tab$Rank), ]
  }
  sm <- data.frame(Accessions = nlevels(droplevels(ent$Genotype)), Vg = vg, Vge = if (scope == "pooled") vge else NA_real_,
                   Ve = ve, H2_Cullis = h2c, Mean_Reliability = rel, H2_Mean = h2_mean,
                   Note = paste(note, collapse = "; "), stringsAsFactors = FALSE)
  list(table = tab, summary = sm)
}

# Every trait: pooled over environments (needs >= 2) and, when by_env, in each environment.
# Returns list(table, summary, errors): one table row per Trait x Environment x Accession, one summary row per Trait x Environment
# (table and summary are NULL when no model could be fitted).
blups_all <- function(wd, traits, design, by_env = TRUE) {
  tabs <- list(); sums <- list(); errs <- character(0)
  for (tr in traits) {
    d_tr <- wd[!is.na(wd[[tr]]), , drop = FALSE]
    envs <- unique(as.character(d_tr$Environment))
    jobs <- c(if (length(envs) > 1) list(list(env = "Pooled", d = d_tr, scope = "pooled")),
              if (by_env || length(envs) == 1) lapply(envs, function(e)
                list(env = e, d = d_tr[as.character(d_tr$Environment) == e, , drop = FALSE], scope = "single")))
    for (j in jobs) {
      td <- model_df(j$d, tr, design)
      if (nrow(td) < 10 || nlevels(td$Genotype) < 5) next
      r <- tryCatch(blup_one(td, design, j$scope), error = function(e) { errs <<- c(errs, sprintf("%s, %s: %s", tr, j$env, conditionMessage(e))); NULL })
      if (is.null(r)) next
      if (!is.null(r$table)) tabs[[length(tabs) + 1]] <- cbind(Trait = tr, Environment = j$env, r$table, stringsAsFactors = FALSE)
      sums[[length(sums) + 1]] <- cbind(Trait = tr, Environment = j$env, r$summary, stringsAsFactors = FALSE)
    }
  }
  list(table = if (length(tabs)) dplyr::bind_rows(tabs) else NULL, summary = if (length(sums)) dplyr::bind_rows(sums) else NULL, errors = errs)
}

# ── Multivariate analyses: one row per genotype ──────────────────────────────
# Traits with no value at all in these data (for example a trait that was not measured in the chosen environments)
# are left out first; otherwise they would leave no accession with complete data.
multivariate_data <- function(wd, traits) {
  measured <- traits[vapply(traits, function(t) any(!is.na(wd[[t]])), logical(1))]
  if (length(measured) < 2) return(NULL)
  gm <- genotype_means(wd, measured)
  gm <- gm[stats::complete.cases(gm[, measured]), , drop = FALSE]
  keep <- measured[vapply(measured, function(t) isTRUE(stats::sd(gm[[t]]) > 0), logical(1))]
  if (nrow(gm) < 5 || length(keep) < 2) return(NULL)
  list(data = gm, traits = keep, dropped = setdiff(traits, keep))
}

correlation_analysis <- function(mv, method = "pearson") {
  method <- if (method %in% c("pearson", "spearman", "kendall")) method else "pearson"
  cm <- cor(mv$data[, mv$traits], use = "complete.obs", method = method)
  # every pair once, with its P value (test on the accession means), strongest first
  pr <- t(utils::combn(mv$traits, 2))
  tab <- dplyr::bind_rows(lapply(seq_len(nrow(pr)), function(i) {
    ct <- suppressWarnings(stats::cor.test(mv$data[[pr[i, 1]]], mv$data[[pr[i, 2]]], method = method))
    data.frame(`Trait 1` = pr[i, 1], `Trait 2` = pr[i, 2], r = round(unname(ct$estimate), 3),
               Accessions = nrow(mv$data), P = signif(max(ct$p.value, 2.2e-16), 4), `Sig.` = sig_stars(ct$p.value),
               stringsAsFactors = FALSE, check.names = FALSE)
  }))
  tab <- tab[order(-abs(tab$r)), ]
  dup <- tab[abs(tab$r) >= 0.999, ]                  # traits that carry the same information
  if (method != "pearson") names(tab)[names(tab) == "r"] <- if (method == "spearman") "rho" else "tau"
  list(matrix = cm, table = tab, duplicates = dup, method = method)
}

# scale = TRUE: traits are standardised first (correlation-matrix PCA); FALSE: covariance-matrix PCA
pca_analysis <- function(mv, scale = TRUE) {
  pr <- prcomp(mv$data[, mv$traits], center = TRUE, scale. = isTRUE(scale))
  v <- pr$sdev^2
  list(scores = as.data.frame(pr$x) %>% mutate(Genotype = mv$data$Genotype),
       variance = data.frame(PC = paste0("PC", seq_along(v)), Variance = v,
                             Prop_Variance = v / sum(v), Cumulative = cumsum(v / sum(v))),
       loadings = as.data.frame(pr$rotation) %>% rownames_to_column("Trait"))
}

# k-means (fixed seed, so results are reproducible) or hierarchical (Ward.D2) clustering of genotypes
cluster_analysis <- function(mv, k, method) {
  trs <- mv$traits
  cd <- scale(as.matrix(mv$data[, trs])); rownames(cd) <- mv$data$Genotype
  if (is.na(k) || k < 2 || nrow(cd) < k + 2) return(NULL)
  hc <- NULL
  if (identical(method, "kmeans")) {
    set.seed(123)
    clusters <- kmeans(cd, centers = k, nstart = 25)$cluster
  } else {
    hc <- hclust(dist(cd), method = "ward.D2")
    clusters <- cutree(hc, k = k)
  }
  cl <- factor(unname(clusters))
  list(clusters = data.frame(Genotype = mv$data$Genotype, Cluster = cl),
       summary = mv$data %>% mutate(Cluster = cl) %>% group_by(Cluster) %>%
         summarise(across(all_of(trs), ~ mean(.x)), Count = n(), .groups = "drop"),
       data = cd, method = method, hclust = hc)
}

# ── Stability and AMMI ───────────────────────────────────────────────────────
# Uses genotypes present in EVERY environment (the formulas need a complete genotype x environment table).
stability_analysis <- function(data, trait) {
  cells <- trait_long(data, trait) %>% filter(!is.na(Value)) %>%
    mutate(Genotype = as.character(Genotype), Environment = as.character(Environment)) %>%
    group_by(Genotype, Environment) %>% summarise(Mean_Value = mean(Value), .groups = "drop")
  n_env <- n_distinct(cells$Environment)
  if (n_env < 2) return(NULL)
  full <- cells %>% count(Genotype) %>% filter(n == n_env) %>% pull(Genotype)
  if (length(full) < 3) return(NULL)
  n_excluded <- n_distinct(cells$Genotype) - length(full)
  cells <- filter(cells, Genotype %in% full)

  mat <- cells %>% pivot_wider(names_from = Environment, values_from = Mean_Value) %>%
    column_to_rownames("Genotype") %>% as.matrix()
  mu  <- mean(mat); gm <- rowMeans(mat); em <- colMeans(mat); idx <- em - mu    # idx = environmental index I_j
  ge  <- sweep(sweep(mat, 1, gm), 2, em) + mu                                    # GE_ij = x_ij - x_i. - x_.j + x_..
  wi  <- rowSums(ge^2)                                                           # Wricke's ecovalence
  dev <- mat - gm
  slope <- as.vector(dev %*% idx) / sum(idx^2)                                   # Finlay-Wilkinson slope
  r2 <- 1 - rowSums((dev - outer(slope, idx))^2) / rowSums(dev^2)

  parameters <- tibble(
    Genotype = rownames(mat), Overall_Mean = gm,
    CV = round(apply(mat, 1, sd) / gm * 100, 2), Ecovalence = round(wi, 3),
    Min_Value = apply(mat, 1, min), Max_Value = apply(mat, 1, max),
    Range = apply(mat, 1, max) - apply(mat, 1, min),
    Slope_b = round(slope, 3), Reg_R2 = round(r2, 3), Stability_Rank = rank(wi, ties.method = "min"))

  ge_long <- as.data.frame(ge) %>% rownames_to_column("Genotype") %>%
    pivot_longer(-Genotype, names_to = "Environment", values_to = "GE_Interaction")
  ammi_data <- cells %>%
    left_join(tibble(Environment = names(em), Env_Mean = em, Env_Index = idx), by = "Environment") %>%
    left_join(ge_long, by = c("Genotype", "Environment")) %>%
    left_join(tibble(Genotype = names(gm), Overall_Mean = gm), by = "Genotype")

  # AMMI: singular value decomposition of the GE matrix
  sv <- svd(ge); k <- min(2, min(dim(ge)) - 1)
  ammi_scores <- NULL; ipca_var <- NA
  if (k >= 1) {
    ipca_var <- sv$d^2 / sum(sv$d^2)
    gs <- sweep(sv$u[, 1:k, drop = FALSE], 2, sqrt(sv$d[1:k]), "*")
    es <- sweep(sv$v[, 1:k, drop = FALSE], 2, sqrt(sv$d[1:k]), "*")
    asv <- if (k == 2 && sv$d[2] > 0) sqrt((sv$d[1]^2 / sv$d[2]^2 * gs[, 1])^2 + gs[, 2]^2) else rep(NA_real_, nrow(gs))
    ammi_scores <- bind_rows(
      tibble(Type = "Genotype", Name = rownames(mat), Mean = gm, IPCA1 = gs[, 1],
             IPCA2 = if (k == 2) gs[, 2] else NA_real_, ASV = asv),
      tibble(Type = "Environment", Name = colnames(mat), Mean = em, IPCA1 = es[, 1],
             IPCA2 = if (k == 2) es[, 2] else NA_real_, ASV = NA_real_))
  }
  list(parameters = parameters, ammi_data = ammi_data, ammi_scores = ammi_scores,
       ipca_var = ipca_var, n_excluded = n_excluded)
}
