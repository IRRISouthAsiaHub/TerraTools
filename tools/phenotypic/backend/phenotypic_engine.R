# phenotypic_engine.R — Advanced Phenotypic Analysis engine
# ---------------------------------------------------------------
# Same statistical logic as the Hugging Face app.R Shiny app, restructured
# from Shiny reactives into a single batch pipeline (run_phenotypic()) that
# runs every analysis once per request and writes outputs to disk, instead
# of recomputing reactively per-tab-click. No math has been changed.
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(lme4)
  library(emmeans)
  library(plotly)
  library(htmlwidgets)
  library(GGally)
  library(corrplot)
  library(agricolae)
  library(car)
  library(broom)
  library(cluster)
  library(factoextra)
  library(pheatmap)
  library(RColorBrewer)
  library(moments)
  library(reshape2)
  library(gridExtra)
  library(openxlsx)
  library(stringr)
  library(scales)
  library(jsonlite)
})

# ── Color palette (unchanged) ────────────────────────────────────────────
get_color_palette <- function(scheme, n) {
  n <- max(n, 1)
  if (scheme == "blue") {
    return(colorRampPalette(c("#E6F0FA", "#1E3A8A"))(n))
  } else if (scheme %in% c("viridis", "plasma", "magma") &&
             requireNamespace("viridisLite", quietly = TRUE)) {
    fn <- switch(scheme, viridis = viridisLite::viridis,
                 plasma = viridisLite::plasma, magma = viridisLite::magma)
    return(fn(n))
  }
  colorRampPalette(c("#1E3A8A", "#6B8E23", "#FF8C00", "#DC143C", "#4169E1"))(n)
}

# ── Data validation (unchanged) ──────────────────────────────────────────
validate_data <- function(data) {
  validation <- list(errors = character(0), warnings = character(0))
  required_cols <- c("Genotype", "Environment", "Rep", "Block")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0)
    validation$errors <- c(validation$errors, paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  for (col in c("Genotype", "Environment", "Rep", "Block")) {
    if (col %in% names(data) && !is.factor(data[[col]])) {
      data[[col]] <- as.factor(data[[col]])
      validation$warnings <- c(validation$warnings, paste(col, "converted to factor"))
    }
  }
  if ("Genotype" %in% names(data) && length(unique(data$Genotype)) < 3)
    validation$errors <- c(validation$errors, "Need at least 3 genotypes")
  if ("Environment" %in% names(data) && length(unique(data$Environment)) < 2)
    validation$errors <- c(validation$errors, "Need at least 2 environments")
  if (ncol(data) >= 10) {
    tc <- names(data)[10:ncol(data)]
    nt <- tc[sapply(data[tc], is.numeric)]
    if (length(nt) == 0)
      validation$errors <- c(validation$errors, "No numeric trait columns from col 10 onward")
  } else {
    validation$errors <- c(validation$errors, "Data has fewer than 10 columns")
  }
  list(data = data, validation = validation)
}

# ── Helper: save a plotly widget as a standalone (non-selfcontained) HTML
#    file. selfcontained=FALSE avoids a pandoc dependency while still being
#    fully offline (dependency files sit in a sibling "_files" folder). ────
save_plotly <- function(p, path) {
  htmlwidgets::saveWidget(plotly::as_widget(p), path, selfcontained = FALSE)
}

# ── Helper: consistent x-axis label rotation logic (unchanged pattern used
#    throughout app.R for environment-labeled axes) ─────────────────────
label_geom <- function(n_env) {
  angle <- if (n_env > 6) 90 else if (n_env > 3) 45 else 0
  hjust <- if (angle == 90) 1 else if (angle == 45) 1 else 0.5
  vjust <- if (angle == 90) 0.5 else 1
  bmargin <- if (angle >= 45) 60 else 20
  list(angle = angle, hjust = hjust, vjust = vjust, bmargin = bmargin)
}

# =============================================================================
# MAIN PIPELINE
# =============================================================================
run_phenotypic <- function(input_file, trait, include_interaction,
                            analysis_mode, selected_environments,
                            n_clusters, clustering_method, color_scheme,
                            outdir) {

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  plots_dir <- file.path(outdir, "plots")
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

  # ── Read + validate ──────────────────────────────────────────────────
  ext <- tolower(tools::file_ext(input_file))
  data <- if (ext == "csv") {
    read.csv(input_file, stringsAsFactors = FALSE)
  } else {
    as.data.frame(readxl::read_excel(input_file))
  }
  vres <- validate_data(data)
  data <- vres$data
  if (length(vres$validation$errors) > 0) {
    stop(paste(vres$validation$errors, collapse = "; "))
  }

  trait_cols <- if (ncol(data) >= 10) names(data)[10:ncol(data)] else character(0)
  traits <- trait_cols[sapply(data[trait_cols], is.numeric)]
  if (!(trait %in% traits)) {
    if (length(traits) == 0) stop("No numeric trait columns found.")
    trait <- traits[1]
  }

  working_data <- data
  active_envs <- NULL
  if (identical(analysis_mode, "specific") && length(selected_environments) > 0) {
    working_data <- data %>% filter(as.character(Environment) %in% selected_environments)
    active_envs <- selected_environments
  }

  out <- list(
    status = "success",
    trait = trait,
    traits_all = traits,
    n_obs = nrow(data), n_vars = ncol(data),
    n_geno = length(unique(data$Genotype)),
    n_env = length(unique(data$Environment)),
    active_envs = active_envs,
    tables = list(), plots = list(), errors = list()
  )

  safe <- function(name, expr) {
    tryCatch(expr, error = function(e) {
      out$errors[[name]] <<- conditionMessage(e)
      NULL
    })
  }

  write_csv_out <- function(df, name) {
    p <- file.path(outdir, paste0(name, ".csv"))
    readr::write_csv(df, p)
    out$tables[[name]] <<- p
  }
  write_plot_out <- function(p_gg, name, height_px = 600, width_px = 900, tooltip = NULL) {
    d <- file.path(plots_dir, name)
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
    html_fp <- file.path(d, paste0(name, ".html"))
    png_fp <- file.path(d, paste0(name, ".png"))

    # Interactive version (unchanged behavior)
    widget <- if (!is.null(tooltip)) {
      ggplotly(p_gg, tooltip = tooltip, height = height_px, width = width_px)
    } else {
      ggplotly(p_gg, height = height_px, width = width_px)
    }
    save_plotly(widget, html_fp)

    # High-quality static companion — rendered directly from the ggplot
    # object at 300 DPI (a real vector-quality re-render, not a screenshot
    # of the interactive widget), so "Download PNG" gives a print-quality
    # image regardless of how big the on-screen iframe is.
    tryCatch(
      ggsave(png_fp, plot = p_gg, width = 11, height = 7, dpi = 300, units = "in", bg = "white"),
      error = function(e) NULL
    )

    out$plots[[name]] <<- list(html = html_fp, png = if (file.exists(png_fp)) png_fp else NULL)
  }
  write_png_out <- function(expr_fn, name, width = 1600, height = 1200) {
    fp <- file.path(plots_dir, paste0(name, ".png"))
    png(fp, width = width, height = height, res = 200)
    tryCatch(expr_fn(), finally = dev.off())
    out$plots[[name]] <<- list(html = NULL, png = fp)
  }

  # =========================================================================
  # ENVIRONMENT-LEVEL ANALYSIS (uses working_data, i.e. env-filtered)
  # =========================================================================

  # -- Per-env descriptive stats --
  env_desc <- safe("env_desc_stats", {
    working_data %>%
      select(Environment, Genotype, all_of(trait)) %>%
      rename(Value = !!trait) %>% filter(!is.na(Value)) %>%
      group_by(Environment) %>%
      summarise(
        N = n(), Mean = round(mean(Value, na.rm = TRUE), 3),
        Median = round(median(Value, na.rm = TRUE), 3),
        SD = round(sd(Value, na.rm = TRUE), 3),
        CV = round(sd(Value, na.rm = TRUE) / mean(Value, na.rm = TRUE) * 100, 2),
        Min = round(min(Value, na.rm = TRUE), 3), Max = round(max(Value, na.rm = TRUE), 3),
        Skewness = round(moments::skewness(Value, na.rm = TRUE), 3),
        Kurtosis = round(moments::kurtosis(Value, na.rm = TRUE), 3),
        Missing = sum(is.na(Value)), .groups = "drop")
  })
  if (!is.null(env_desc)) write_csv_out(env_desc, "env_desc_stats")

  # -- Env boxplot + violin --
  env_df <- safe("env_df", {
    working_data %>% select(Environment, Genotype, all_of(trait)) %>%
      rename(Value = !!trait) %>% filter(!is.na(Value)) %>%
      mutate(Environment = factor(Environment))
  })
  if (!is.null(env_df)) {
    n_env <- length(unique(env_df$Environment))
    lg <- label_geom(n_env)
    safe("env_boxplot", {
      p <- ggplot(env_df, aes(x = Environment, y = Value, fill = Environment)) +
        geom_boxplot(alpha = 0.7, outlier.colour = "red", outlier.size = 2) +
        geom_jitter(width = 0.2, alpha = 0.35, size = 0.8, color = "#444444") +
        labs(title = paste("Boxplot of", trait, "per Environment"), x = NULL, y = trait) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = lg$angle, hjust = lg$hjust, vjust = lg$vjust, size = 11, color = "#000"),
              axis.text.y = element_text(size = 11, color = "#000"), axis.title = element_text(size = 12, color = "#000"),
              plot.title = element_text(size = 13, color = "#000", face = "bold"),
              plot.margin = margin(t = 10, r = 10, b = lg$bmargin, l = 10), legend.position = "none") +
        scale_fill_manual(values = get_color_palette(color_scheme, n_env)) +
        scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 12))
      write_plot_out(p, "env_boxplot", height_px = 600, width_px = 900)
    })
    safe("env_violin_plot", {
      p <- ggplot(env_df, aes(x = Environment, y = Value, fill = Environment)) +
        geom_violin(alpha = 0.7, trim = FALSE) +
        geom_boxplot(width = 0.08, fill = "white", alpha = 0.85, outlier.size = 1.5) +
        stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "#1E3A8A") +
        labs(title = paste("Violin Plot of", trait, "per Environment"), x = NULL, y = trait) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = lg$angle, hjust = lg$hjust, vjust = lg$vjust, size = 11, color = "#000"),
              axis.text.y = element_text(size = 11, color = "#000"), axis.title = element_text(size = 12, color = "#000"),
              plot.title = element_text(size = 13, color = "#000", face = "bold"),
              plot.margin = margin(t = 10, r = 10, b = lg$bmargin, l = 10), legend.position = "none") +
        scale_fill_manual(values = get_color_palette(color_scheme, n_env)) +
        scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 12))
      write_plot_out(p, "env_violin_plot", height_px = 600, width_px = 900)
    })
  }

  # -- Per-env ANOVA --
  env_anova <- safe("env_anova", {
    envs <- unique(as.character(working_data$Environment))
    results <- lapply(envs, function(env) {
      df <- working_data %>% filter(as.character(Environment) == env) %>%
        select(Genotype, Rep, Block, all_of(trait)) %>% rename(Value = !!trait) %>%
        filter(!is.na(Value)) %>%
        mutate(Genotype = as.factor(Genotype), Rep = as.factor(Rep), Block = as.factor(Block))
      if (nrow(df) < 6 || length(unique(df$Genotype)) < 2) return(NULL)
      tryCatch({
        m <- aov(Value ~ Genotype + Rep + Block, data = df)
        sm <- broom::tidy(m); sm$Environment <- env; sm
      }, error = function(e) NULL)
    })
    do.call(rbind, Filter(Negate(is.null), results))
  })
  if (!is.null(env_anova)) write_csv_out(env_anova, "env_anova")

  # -- Genotype means per env + rank --
  env_means <- safe("env_means", {
    working_data %>% select(Genotype, Environment, all_of(trait)) %>%
      rename(Value = !!trait) %>% filter(!is.na(Value)) %>%
      group_by(Genotype, Environment) %>%
      summarise(Mean = round(mean(Value, na.rm = TRUE), 3), SD = round(sd(Value, na.rm = TRUE), 3),
                N = n(), .groups = "drop")
  })
  if (!is.null(env_means)) {
    write_csv_out(env_means, "env_means")
    safe("env_means_plot", {
      n_geno <- length(unique(env_means$Genotype))
      txt_size <- if (n_geno > 40) 6 else if (n_geno > 20) 8 else 10
      dyn_h <- max(600, n_geno * 18 + 100)
      p <- ggplot(env_means, aes(x = reorder(Genotype, Mean), y = Mean, fill = Environment)) +
        geom_col(position = "dodge", alpha = 0.85) + coord_flip() +
        labs(title = paste("Genotype Means per Environment —", trait), x = NULL, y = "Mean Value") +
        theme_minimal() +
        theme(axis.text.y = element_text(size = txt_size, color = "#000", lineheight = 0.9),
              axis.text.x = element_text(size = 10, color = "#000"), axis.title = element_text(size = 12, color = "#000"),
              plot.title = element_text(size = 13, color = "#000", face = "bold")) +
        scale_fill_manual(values = get_color_palette(color_scheme, length(unique(env_means$Environment)))) +
        scale_x_discrete(labels = function(x) stringr::str_trunc(x, 20))
      write_plot_out(p, "env_means_plot", height_px = dyn_h, width_px = 900)
    })

    env_rank <- safe("env_rank", {
      ranks <- env_means %>% group_by(Environment) %>%
        mutate(Rank = rank(-Mean, ties.method = "average")) %>% ungroup() %>%
        select(Genotype, Environment, Mean, Rank)
      wide_rank <- ranks %>% select(Genotype, Environment, Rank) %>%
        pivot_wider(names_from = Environment, values_from = Rank, names_prefix = "Rank_")
      list(long = ranks, wide_rank = wide_rank)
    })
    if (!is.null(env_rank)) {
      write_csv_out(env_rank$wide_rank, "env_rank")
      safe("env_rank_plot", {
        df <- env_rank$long
        n_env2 <- length(unique(df$Environment)); n_geno2 <- length(unique(df$Genotype))
        lg2 <- label_geom(n_env2); dyn_h <- max(600, n_geno2 * 14 + 120)
        p <- ggplot(df, aes(x = Environment, y = Rank, color = Genotype, group = Genotype,
                             text = paste0("Genotype: ", Genotype, "<br>Rank: ", round(Rank, 1), "<br>Mean: ", round(Mean, 3)))) +
          geom_line(linewidth = 0.9, alpha = 0.75) + geom_point(size = 2.5, alpha = 0.9) +
          scale_y_reverse(breaks = seq(1, n_geno2, by = max(1, floor(n_geno2 / 10)))) +
          labs(title = paste("Genotype Rank Across Environments —", trait), x = NULL, y = "Rank (1 = Best)") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = lg2$angle, hjust = lg2$hjust, vjust = lg2$vjust, size = 11, color = "#000"),
                axis.text.y = element_text(size = 10, color = "#000"), axis.title.y = element_text(size = 12, color = "#000"),
                plot.title = element_text(size = 13, color = "#000", face = "bold"),
                legend.position = if (n_geno2 > 25) "none" else "right") +
          scale_color_manual(values = get_color_palette(color_scheme, n_geno2)) +
          scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))
        write_plot_out(p, "env_rank_plot", height_px = dyn_h, width_px = 900, tooltip = "text")
      })
    }
  }

  # -- Env-env correlation --
  env_cor <- safe("env_cor", {
    wide <- working_data %>% select(Genotype, Environment, all_of(trait)) %>%
      rename(Value = !!trait) %>% filter(!is.na(Value)) %>%
      group_by(Genotype, Environment) %>% summarise(Mean = mean(Value, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Environment, values_from = Mean) %>% column_to_rownames("Genotype")
    if (ncol(wide) < 2) return(NULL)
    cor_mat <- cor(wide, use = "pairwise.complete.obs")
    cor_long <- cor_mat %>% as.data.frame() %>% rownames_to_column("Env1") %>%
      pivot_longer(-Env1, names_to = "Env2", values_to = "Correlation") %>%
      filter(Env1 != Env2) %>% mutate(Correlation = round(Correlation, 3))
    list(matrix = cor_mat, table = cor_long)
  })
  if (!is.null(env_cor)) {
    write_csv_out(env_cor$table, "env_cor")
    safe("env_cor_heatmap", {
      write_png_out(function() {
        corrplot::corrplot(env_cor$matrix, method = "color", type = "upper", order = "hclust",
                            tl.cex = 0.9, tl.col = "#000000", tl.srt = 45, addCoef.col = "#000000",
                            number.cex = 0.8, title = "Environment Correlation Matrix", mar = c(0, 0, 2, 0),
                            col = colorRampPalette(get_color_palette(color_scheme, 3))(200))
      }, "env_cor_heatmap")
    })
  }

  # -- Per-env BLUEs (lmer + emmeans) --
  env_blues <- safe("env_blues", {
    envs <- unique(as.character(working_data$Environment))
    results <- lapply(envs, function(env) {
      df <- working_data %>% filter(as.character(Environment) == env) %>%
        select(Genotype, Rep, Block, all_of(trait)) %>% rename(Value = !!trait) %>%
        filter(!is.na(Value)) %>%
        mutate(Genotype = as.factor(Genotype), Rep = as.factor(Rep), Block = as.factor(Block))
      if (nrow(df) < 6 || length(unique(df$Genotype)) < 2) return(NULL)
      tryCatch({
        m <- lmer(Value ~ Genotype + (1 | Rep) + (1 | Block), data = df)
        em <- as.data.frame(emmeans(m, ~Genotype))
        em$Environment <- env
        em %>% select(Environment, Genotype, emmean, SE) %>% rename(BLUE = emmean)
      }, error = function(e) NULL)
    })
    do.call(rbind, Filter(Negate(is.null), results))
  })
  if (!is.null(env_blues)) {
    write_csv_out(env_blues, "env_blues")
    safe("env_blues_plot", {
      n_geno3 <- length(unique(env_blues$Genotype)); n_env3 <- length(unique(env_blues$Environment))
      txt_size <- if (n_geno3 > 40) 6 else if (n_geno3 > 20) 8 else 10
      dyn_h <- max(600, n_geno3 * 16 * max(1, ceiling(n_env3 / 3)) + 120)
      p <- ggplot(env_blues, aes(x = reorder(Genotype, BLUE), y = BLUE, color = Environment)) +
        geom_point(size = 2, alpha = 0.85) +
        geom_errorbar(aes(ymin = BLUE - SE, ymax = BLUE + SE), width = 0.4, alpha = 0.65) +
        facet_wrap(~Environment, scales = "free_x", ncol = min(n_env3, 3)) + coord_flip() +
        labs(title = paste("Within-Environment BLUEs —", trait), x = NULL, y = "BLUE") +
        theme_minimal() +
        theme(axis.text.y = element_text(size = txt_size, color = "#000", lineheight = 0.9),
              strip.text = element_text(size = 11, face = "bold", color = "#1E3A8A"),
              plot.title = element_text(size = 13, color = "#000", face = "bold"), legend.position = "none") +
        scale_color_manual(values = get_color_palette(color_scheme, n_env3))
      write_plot_out(p, "env_blues_plot", height_px = dyn_h, width_px = 900)
    })
  }

  # =========================================================================
  # ACROSS-ENVIRONMENT (GLOBAL) ANALYSES
  # =========================================================================

  # -- Descriptive stats (all traits) --
  desc_stats <- safe("descriptive_stats", {
    working_data %>% select(all_of(traits)) %>%
      pivot_longer(everything(), names_to = "Trait", values_to = "Value") %>%
      group_by(Trait) %>%
      summarise(N = sum(!is.na(Value)), Mean = round(mean(Value, na.rm = TRUE), 3),
                Median = round(median(Value, na.rm = TRUE), 3), SD = round(sd(Value, na.rm = TRUE), 3),
                CV = round(sd(Value, na.rm = TRUE) / mean(Value, na.rm = TRUE) * 100, 2),
                Min = round(min(Value, na.rm = TRUE), 3), Max = round(max(Value, na.rm = TRUE), 3),
                Q1 = round(quantile(Value, .25, na.rm = TRUE), 3), Q3 = round(quantile(Value, .75, na.rm = TRUE), 3),
                Skewness = round(moments::skewness(Value, na.rm = TRUE), 3),
                Kurtosis = round(moments::kurtosis(Value, na.rm = TRUE), 3),
                Missing = sum(is.na(Value)), .groups = "drop")
  })
  if (!is.null(desc_stats)) write_csv_out(desc_stats, "descriptive_stats")

  safe("distribution_plot", {
    td <- working_data[[trait]]; td <- td[!is.na(td)]
    if (length(td) >= 2) {
      pal <- get_color_palette(color_scheme, 2)
      p <- data.frame(Value = td) %>% ggplot(aes(x = Value)) +
        geom_histogram(aes(y = after_stat(density)), bins = 30, fill = pal[1], alpha = 0.7) +
        geom_density(color = pal[2], linewidth = 1) +
        labs(title = paste("Distribution of", trait), x = trait, y = "Density") +
        theme_minimal()
      write_plot_out(p, "distribution_plot", height_px = 600, width_px = 900)
    }
  })

  # -- ANOVA --
  anova_res <- safe("anova_analysis", {
    td <- working_data %>% select(Genotype, Environment, Rep, Block, all_of(trait)) %>%
      rename(Value = !!trait) %>% filter(!is.na(Value)) %>%
      mutate(across(c(Genotype, Environment, Rep, Block), as.factor))
    if (nrow(td) < 10) return(NULL)
    sm <- lm(Value ~ Environment + Genotype, data = td)
    rd <- data.frame(Fitted = fitted(sm), Residuals = residuals(sm),
                      Standardized_Residuals = rstandard(sm), Genotype = td$Genotype)
    list(simple_model = sm, residuals = rd)
  })
  if (!is.null(anova_res)) {
    write_csv_out(broom::tidy(anova_res$simple_model), "anova_table")
    safe("residual_plot", {
      rd <- anova_res$residuals
      pal <- get_color_palette(color_scheme, 2)
      p <- ggplot(rd, aes(x = Fitted, y = Residuals, color = Genotype)) +
        geom_point(size = 1, alpha = 0.6) +
        geom_smooth(method = "loess", se = FALSE, color = pal[2]) +
        geom_hline(yintercept = 0, linetype = "dashed", color = pal[2]) +
        labs(title = "Residuals vs Fitted", x = "Fitted", y = "Residuals") + theme_minimal() +
        scale_color_manual(values = get_color_palette(color_scheme, length(unique(rd$Genotype))))
      write_plot_out(p, "residual_plot", height_px = 600, width_px = 900)
    })
    safe("qq_plot", {
      rd <- anova_res$residuals
      pal <- get_color_palette(color_scheme, 2)
      p <- ggplot(rd, aes(sample = Standardized_Residuals, color = Genotype)) +
        stat_qq(size = 1) + stat_qq_line(color = pal[2]) +
        labs(title = "Normal Q-Q Plot", x = "Theoretical Quantiles", y = "Sample Quantiles") + theme_minimal() +
        scale_color_manual(values = get_color_palette(color_scheme, length(unique(rd$Genotype))))
      write_plot_out(p, "qq_plot", height_px = 600, width_px = 900)
    })
  }

  # -- Mixed Model / BLUEs (across environments) --
  mm_res <- safe("mixed_model_analysis", {
    td <- working_data %>% select(Genotype, Environment, Rep, Block, all_of(trait)) %>%
      rename(Value = !!trait) %>% filter(!is.na(Value)) %>%
      mutate(across(c(Genotype, Environment, Rep, Block), as.factor))
    if (nrow(td) < 10) return(NULL)
    model <- if (include_interaction)
      lmer(Value ~ Environment + Genotype + Environment:Genotype + (1 | Rep) + (1 | Block), data = td)
    else
      lmer(Value ~ Environment + Genotype + (1 | Rep) + (1 | Block), data = td)
    blues_df <- as.data.frame(emmeans(model, ~ Genotype | Environment)) %>%
      select(Genotype, Environment, emmean, SE) %>% rename(BLUE = emmean)
    blues_df
  })
  if (!is.null(mm_res)) {
    write_csv_out(mm_res, "mixed_model_blues")
    safe("blues_plot", {
      n_env4 <- length(unique(mm_res$Environment)); n_gen4 <- length(unique(mm_res$Genotype))
      lg4 <- label_geom(n_env4)
      p <- ggplot(mm_res, aes(x = Environment, y = BLUE, color = Genotype, group = Genotype)) +
        geom_line(linewidth = 1) + geom_point(size = 2) +
        labs(title = paste("BLUEs for", trait), x = NULL, y = "BLUE") + theme_minimal() +
        theme(axis.text.x = element_text(angle = lg4$angle, hjust = lg4$hjust, vjust = lg4$vjust),
              legend.position = if (n_gen4 > 25) "none" else "right") +
        scale_color_manual(values = get_color_palette(color_scheme, n_gen4))
      write_plot_out(p, "blues_plot", height_px = 600, width_px = 900)
    })
  }

  # -- Heritability (all traits) --
  herit_res <- safe("heritability_analysis", {
    results <- list()
    for (tr in traits) {
      tryCatch({
        td <- working_data %>% select(Genotype, Environment, Rep, Block, all_of(tr)) %>%
          rename(Value = !!tr) %>% filter(!is.na(Value)) %>%
          mutate(across(c(Genotype, Environment, Rep, Block), as.factor))
        if (nrow(td) < 10 || length(unique(td$Environment)) < 2) next
        m <- lmer(Value ~ (1 | Genotype) + (1 | Environment) + (1 | Genotype:Environment) +
                     (1 | Rep) + (1 | Block), data = td)
        vc <- as.data.frame(VarCorr(m)); vc$Trait <- tr
        vg <- vc$vcov[vc$grp == "Genotype"]; vg <- if (length(vg) == 0) 0 else vg
        vge <- vc$vcov[vc$grp == "Genotype:Environment"]; vge <- if (length(vge) == 0) 0 else vge
        ve <- vc$vcov[vc$grp == "Residual"]; ve <- if (length(ve) == 0) 1 else ve
        ne <- length(unique(td$Environment)); nr <- length(unique(td$Rep))
        h2b <- vg / (vg + vge + ve)
        h2e <- vg / (vg + (vge / ne) + (ve / (ne * nr)))
        results[[tr]] <- list(trait = tr, variance_components = vc, h2_broad = h2b, h2_entry = h2e)
      }, error = function(e) NULL)
    }
    if (length(results) == 0) NULL else results
  })
  if (!is.null(herit_res)) {
    vc_df <- do.call(rbind, lapply(herit_res, `[[`, "variance_components"))
    write_csv_out(vc_df, "variance_components")
    h2_df <- data.frame(Trait = names(herit_res),
                         H2_Broad = sapply(herit_res, function(x) round(x$h2_broad, 3)),
                         H2_Entry = sapply(herit_res, function(x) round(x$h2_entry, 3)))
    write_csv_out(h2_df, "heritability")
    safe("heritability_plot", {
      dfp <- h2_df %>% pivot_longer(cols = c(H2_Broad, H2_Entry), names_to = "Type", values_to = "Heritability")
      pal <- get_color_palette(color_scheme, 2)
      p <- ggplot(dfp, aes(x = Trait, y = Heritability, fill = Type)) +
        geom_col(position = "dodge", alpha = 0.8) +
        geom_hline(yintercept = c(.3, .6), linetype = "dashed", color = pal[2]) +
        labs(title = "Heritability by Trait", x = "Trait", y = "Heritability") + theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_fill_manual(values = pal)
      write_plot_out(p, "heritability_plot", height_px = 600, width_px = 900)
    })
  }

  # -- Correlations (all traits) --
  cor_res <- safe("correlation_analysis", {
    td <- working_data %>% select(all_of(traits)) %>% na.omit()
    if (nrow(td) < 10 || length(traits) < 2) return(NULL)
    cm <- cor(td, use = "complete.obs")
    ct <- cm %>% as.data.frame() %>% rownames_to_column("Trait1") %>%
      pivot_longer(-Trait1, names_to = "Trait2", values_to = "Correlation") %>%
      filter(Trait1 != Trait2) %>%
      mutate(Correlation = round(Correlation, 3), Abs_Correlation = abs(Correlation)) %>%
      arrange(desc(Abs_Correlation))
    list(matrix = cm, table = ct)
  })
  if (!is.null(cor_res)) {
    write_csv_out(cor_res$table, "correlation_table")
    safe("correlation_heatmap", {
      write_png_out(function() {
        corrplot::corrplot(cor_res$matrix, method = "color", type = "upper", order = "hclust",
                            tl.cex = 0.8, tl.col = "#000000", tl.srt = 45, addCoef.col = "#000000",
                            number.cex = 0.7, title = "Trait Correlation Matrix", mar = c(0, 0, 2, 0),
                            col = colorRampPalette(get_color_palette(color_scheme, 3))(200))
      }, "correlation_heatmap")
    })
    safe("correlation_network", {
      write_png_out(function() {
        corrplot::corrplot(cor_res$matrix, method = "circle", type = "full", order = "hclust",
                            tl.cex = 0.8, tl.col = "#000000", title = "Correlation Network", mar = c(0, 0, 2, 0),
                            col = colorRampPalette(get_color_palette(color_scheme, 3))(200))
      }, "correlation_network")
    })
  }

  # -- PCA --
  pca_res <- safe("pca_analysis", {
    td <- working_data %>% select(Genotype, all_of(traits)) %>% na.omit()
    if (nrow(td) < 10 || length(traits) < 2) return(NULL)
    pd <- scale(td %>% select(-Genotype))
    pr <- prcomp(pd, center = FALSE, scale. = FALSE)
    scores <- as.data.frame(pr$x) %>% mutate(Genotype = td$Genotype)
    vexp <- data.frame(PC = paste0("PC", 1:length(pr$sdev)), Variance = pr$sdev^2,
                        Prop_Variance = pr$sdev^2 / sum(pr$sdev^2),
                        Cumulative = cumsum(pr$sdev^2 / sum(pr$sdev^2)))
    loadings <- as.data.frame(pr$rotation) %>% rownames_to_column("Trait")
    list(scores = scores, variance = vexp, loadings = loadings)
  })
  if (!is.null(pca_res)) {
    write_csv_out(pca_res$scores, "pca_scores")
    write_csv_out(pca_res$variance, "pca_variance")
    write_csv_out(pca_res$loadings, "pca_loadings")
    safe("pca_plot", {
      sc <- pca_res$scores
      p <- ggplot(sc, aes(x = PC1, y = PC2, color = Genotype, text = Genotype)) +
        geom_point(size = 2, alpha = 0.8) +
        labs(title = "PCA Biplot",
             x = paste0("PC1 (", round(pca_res$variance$Prop_Variance[1] * 100, 1), "%)"),
             y = paste0("PC2 (", round(pca_res$variance$Prop_Variance[2] * 100, 1), "%)")) +
        theme_minimal() +
        scale_color_manual(values = get_color_palette(color_scheme, length(unique(sc$Genotype))))
      write_plot_out(p, "pca_plot", height_px = 600, width_px = 900, tooltip = c("text", "x", "y"))
    })
    safe("pca_variance_plot", {
      vd <- pca_res$variance[1:min(10, nrow(pca_res$variance)), ]
      pal <- get_color_palette(color_scheme, 2)
      p <- ggplot(vd, aes(x = PC, y = Prop_Variance)) +
        geom_col(fill = pal[1], alpha = 0.8) +
        geom_line(aes(y = Cumulative, group = 1), color = pal[2], linewidth = 1) +
        geom_point(aes(y = Cumulative), color = pal[2], size = 2) +
        labs(title = "PCA Variance Explained", x = "PC", y = "Proportion of Variance") +
        theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      write_plot_out(p, "pca_variance_plot", height_px = 600, width_px = 900)
    })
  }

  # -- Clustering --
  clust_res <- safe("clustering_analysis", {
    td <- working_data %>% select(Genotype, all_of(traits)) %>% na.omit()
    if (nrow(td) < n_clusters + 2 || length(traits) < 2) return(NULL)
    cd <- scale(td %>% select(-Genotype)); rownames(cd) <- td$Genotype
    if (clustering_method == "kmeans") {
      kr <- kmeans(cd, centers = n_clusters, nstart = 25)
      clusters <- kr$cluster; hclust_res <- NULL
    } else {
      hr <- hclust(dist(cd), method = "ward.D2")
      clusters <- cutree(hr, k = n_clusters); hclust_res <- hr
    }
    cdf <- data.frame(Genotype = td$Genotype, Cluster = as.factor(clusters))
    cs <- td %>% mutate(Cluster = as.factor(clusters)) %>% group_by(Cluster) %>%
      summarise_at(vars(all_of(traits)), mean, na.rm = TRUE) %>%
      mutate(Count = as.numeric(table(clusters)))
    list(clusters = cdf, summary = cs, data = cd, method = clustering_method, hclust = hclust_res)
  })
  if (!is.null(clust_res)) {
    write_csv_out(clust_res$clusters, "cluster_assignments")
    write_csv_out(clust_res$summary, "cluster_summary")
    safe("cluster_plot", {
      cd <- clust_res$data; clu <- clust_res$clusters
      if (ncol(cd) > 2) {
        pf <- prcomp(cd)
        pdd <- data.frame(PC1 = pf$x[, 1], PC2 = pf$x[, 2], Genotype = clu$Genotype, Cluster = clu$Cluster)
        p <- ggplot(pdd, aes(x = PC1, y = PC2, color = Cluster, text = Genotype)) +
          geom_point(size = 3, alpha = 0.8) + labs(title = "Clustering (PCA Space)", x = "PC1", y = "PC2")
      } else {
        pdd <- data.frame(X = cd[, 1], Y = cd[, 2], Genotype = clu$Genotype, Cluster = clu$Cluster)
        p <- ggplot(pdd, aes(x = X, y = Y, color = Cluster, text = Genotype)) +
          geom_point(size = 3, alpha = 0.8) + labs(title = "Genotype Clustering", x = traits[1], y = traits[2])
      }
      p <- p + theme_minimal() + scale_color_manual(values = get_color_palette(color_scheme, n_clusters))
      write_plot_out(p, "cluster_plot", height_px = 600, width_px = 900, tooltip = c("text", "color"))
    })
    if (identical(clust_res$method, "hierarchical") && !is.null(clust_res$hclust)) {
      safe("dendrogram", {
        write_png_out(function() {
          plot(clust_res$hclust, main = "Hierarchical Clustering Dendrogram",
               xlab = "Genotypes", ylab = "Distance", cex = 0.8, col = "#000000")
          rect.hclust(clust_res$hclust, k = n_clusters, border = get_color_palette(color_scheme, 2)[2])
        }, "dendrogram")
      })
    }
  }

  # -- Stability / AMMI (always uses FULL data, not env-filtered — matches app.R) --
  stab_res <- safe("stability_analysis", {
    sd2 <- data %>% select(Genotype, Environment, all_of(trait)) %>%
      rename(Value = !!trait) %>% filter(!is.na(Value)) %>%
      group_by(Genotype, Environment) %>% summarise(Mean_Value = mean(Value, na.rm = TRUE), .groups = "drop")
    if (nrow(sd2) < 10) return(NULL)
    om <- mean(sd2$Mean_Value, na.rm = TRUE)
    sp <- sd2 %>% group_by(Genotype) %>%
      summarise(Overall_Mean = mean(Mean_Value, na.rm = TRUE),
                CV = sd(Mean_Value, na.rm = TRUE) / mean(Mean_Value, na.rm = TRUE) * 100,
                Ecovalence = sum((Mean_Value - om)^2, na.rm = TRUE),
                Min_Value = min(Mean_Value, na.rm = TRUE), Max_Value = max(Mean_Value, na.rm = TRUE),
                Range = Max_Value - Min_Value, .groups = "drop") %>%
      mutate(CV = round(CV, 2), Ecovalence = round(Ecovalence, 3), Stability_Rank = rank(Ecovalence))
    em <- sd2 %>% group_by(Environment) %>% summarise(Env_Mean = mean(Mean_Value, na.rm = TRUE), .groups = "drop")
    ad <- sd2 %>% left_join(em, by = "Environment") %>%
      mutate(Env_Index = Env_Mean - om, GE_Interaction = Mean_Value - om - Env_Index) %>%
      left_join(sp %>% select(Genotype, Overall_Mean), by = "Genotype")
    list(parameters = sp, ammi_data = ad)
  })
  if (!is.null(stab_res)) {
    write_csv_out(stab_res$parameters, "stability_parameters")
    write_csv_out(stab_res$ammi_data, "ammi_data")
    safe("stability_plot", {
      sp <- stab_res$parameters
      pal <- get_color_palette(color_scheme, 2)
      p <- ggplot(sp, aes(x = Overall_Mean, y = Ecovalence, color = Genotype, text = Genotype)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_hline(yintercept = median(sp$Ecovalence), linetype = "dashed", color = pal[2]) +
        geom_vline(xintercept = median(sp$Overall_Mean), linetype = "dashed", color = pal[2]) +
        labs(title = paste("Stability —", trait), x = "Mean Performance", y = "Ecovalence (Wi\u00b2)") +
        theme_minimal() + scale_color_manual(values = get_color_palette(color_scheme, length(unique(sp$Genotype))))
      write_plot_out(p, "stability_plot", height_px = 600, width_px = 900, tooltip = c("text", "x", "y"))
    })
    safe("ammi_plot", {
      ad <- stab_res$ammi_data
      n_env5 <- length(unique(ad$Environment)); n_gen5 <- length(unique(ad$Genotype))
      lg5 <- label_geom(n_env5)
      p <- ggplot(ad, aes(x = Environment, y = Mean_Value, color = Genotype, group = Genotype)) +
        geom_line(linewidth = 1) + geom_point(size = 2) +
        labs(title = paste("G\u00d7E Interaction —", trait), x = NULL, y = "Mean Value") + theme_minimal() +
        theme(axis.text.x = element_text(angle = lg5$angle, hjust = lg5$hjust, vjust = lg5$vjust),
              legend.position = if (n_gen5 > 25) "none" else "right") +
        scale_color_manual(values = get_color_palette(color_scheme, n_gen5))
      write_plot_out(p, "ammi_plot", height_px = 600, width_px = 900, tooltip = c("color", "x", "y"))
    })
  }

  # =========================================================================
  # EXCEL WORKBOOKS — mirrors the two download buttons in app.R
  # =========================================================================
  wb_all <- createWorkbook()
  add_sheet <- function(wb, name, df) {
    if (is.null(df)) return(invisible(NULL))
    addWorksheet(wb, name); writeData(wb, name, df)
  }
  add_sheet(wb_all, "Data_Summary", data.frame(
    Metric = c("Observations", "Variables", "Genotypes", "Environments", "Traits", "Active Environments", "Date"),
    Value = c(nrow(working_data), ncol(working_data), length(unique(working_data$Genotype)),
              length(unique(working_data$Environment)), length(traits),
              if (!is.null(active_envs)) paste(active_envs, collapse = ", ") else "All",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"))))
  add_sheet(wb_all, "Descriptive_Statistics", desc_stats)
  if (!is.null(anova_res)) add_sheet(wb_all, "ANOVA_Results", broom::tidy(anova_res$simple_model))
  add_sheet(wb_all, "BLUEs", mm_res)
  if (!is.null(herit_res)) {
    add_sheet(wb_all, "Variance_Components", vc_df)
    add_sheet(wb_all, "Heritability_Estimates", h2_df)
  }
  if (!is.null(cor_res)) {
    add_sheet(wb_all, "Correlation_Matrix", rownames_to_column(as.data.frame(cor_res$matrix), "Trait"))
    add_sheet(wb_all, "Pairwise_Correlations", cor_res$table)
  }
  if (!is.null(pca_res)) {
    add_sheet(wb_all, "PCA_Scores", pca_res$scores)
    add_sheet(wb_all, "PCA_Variance", pca_res$variance)
    add_sheet(wb_all, "PCA_Loadings", pca_res$loadings)
  }
  if (!is.null(clust_res)) {
    add_sheet(wb_all, "Cluster_Assignments", clust_res$clusters)
    add_sheet(wb_all, "Cluster_Summary", clust_res$summary)
  }
  if (!is.null(stab_res)) {
    add_sheet(wb_all, "Stability_Parameters", stab_res$parameters)
    add_sheet(wb_all, "GxE_Interaction", stab_res$ammi_data)
  }
  all_xlsx_path <- file.path(outdir, "Phenotypic_Analysis.xlsx")
  saveWorkbook(wb_all, all_xlsx_path, overwrite = TRUE)
  out$excel_all <- all_xlsx_path

  wb_env <- createWorkbook()
  add_sheet(wb_env, "Env_Descriptive_Stats", env_desc)
  add_sheet(wb_env, "Env_Genotype_Means", env_means)
  if (!is.null(env_rank)) add_sheet(wb_env, "Env_Rank_Table", env_rank$wide_rank)
  add_sheet(wb_env, "Env_ANOVA_Results", env_anova)
  if (!is.null(env_cor)) add_sheet(wb_env, "Env_Correlations", env_cor$table)
  add_sheet(wb_env, "Env_BLUEs", env_blues)
  env_xlsx_path <- file.path(outdir, "Environment_Analysis.xlsx")
  saveWorkbook(wb_env, env_xlsx_path, overwrite = TRUE)
  out$excel_env <- env_xlsx_path

  # Write the master summary.json
  summary_path <- file.path(outdir, "summary.json")
  jsonlite::write_json(out, summary_path, auto_unbox = TRUE, na = "null", null = "null")

  out
}
