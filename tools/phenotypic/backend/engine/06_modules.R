# =============================================================================
# ANALYSIS MODULES — one function per analysis tab
# =============================================================================
# Every module is a function without arguments. run_phenotypic() gives each one the environment of the running
# request (the data, `settings`, helpers such as safe(), cfg(), write_table(), fig() ...), so a module only
# contains its own analysis: it reads its traits and settings, computes, and writes its tables and figures.
# A module touches nothing that belongs to another one. Anything it must hand on to later runs (`ta`, `sed_keep`,
# `h2_keep`) is assigned with <<-.

# =============================================================================
# DESCRIPTIVE STATISTICS
# =============================================================================
mod_descriptive <- function() {
  trs <- mod_traits(); ptr <- mod_trait("plot_trait", trs)
  desc <- safe("descriptive_stats", describe_traits(working_data, trs, by_env = isTRUE(cfg("by_env", TRUE))))
  if (!is.null(desc))
    write_table(tidy_desc(desc), "descriptive_stats", "Descriptive_Statistics",
                "Mean, SD, CV, range and shape of every trait: all environments together and in each environment")
  env_df <- safe("env_df", trait_long(working_data, ptr) %>% filter(!is.na(Value)) %>% mutate(Environment = factor(Environment)))
  if (!is.null(env_df)) {
    n_env_box <- nlevels(env_df$Environment); box_size <- c(max(5.5, 1.5 + 1.1 * n_env_box), 5)
    fig("env_boxplot",     write_plot_out(plot_env_box(env_df, ptr, color_scheme), "env_boxplot", size_in = box_size))
    fig("env_violin_plot", write_plot_out(plot_env_violin(env_df, ptr, color_scheme), "env_violin_plot", size_in = box_size))
  }
  fig("distribution_plot", {
    td <- working_data[[ptr]]; td <- td[!is.na(td)]
    if (length(td) >= 2) write_plot_out(plot_distribution(td, ptr, color_scheme), "distribution_plot", size_in = FIG_SMALL)
  })
}

# =============================================================================
# ENVIRONMENTS: means, ranks and correlation between environments
# =============================================================================
mod_environments <- function() {
  trait <- mod_trait()
  env_means <- safe("env_means", env_means_table(working_data, trait))
  if (!is.null(env_means)) {
    write_table(tidy_means_ranks(env_means), "env_means", "Env_Means_Ranks",
                paste("Mean and rank of every accession for", trait, "in each environment"))
    fig("env_means_plot", write_plot_out(plot_env_means(env_means, trait, color_scheme), "env_means_plot",
                                         height_px = max(600, n_distinct(env_means$Genotype) * 18 + 100),
                                         size_in = fig_tall(n_distinct(env_means$Genotype))))
    env_rank <- safe("env_rank", env_rank_tables(env_means))
    if (!is.null(env_rank) && n_distinct(env_rank$Environment) > 1)       # ranks across environments need at least two
      fig("env_rank_plot", write_plot_out(plot_env_rank(env_rank, trait, color_scheme), "env_rank_plot",
                                          tooltip = "text", size_in = FIG_WIDE))
  }
  env_cor <- safe("env_cor", env_correlation(working_data, trait))
  if (!is.null(env_cor)) {
    write_table(tidy_env_cor(env_cor$matrix), "env_cor", "Env_Correlations",
                paste("Correlation of accession means of", trait, "between environments"))
    fig("env_cor_heatmap", write_png_out(function() {
      corrplot::corrplot(env_cor$matrix, method = "color", type = "upper", order = "hclust",
                         tl.cex = 0.9, tl.col = "#000000", tl.srt = 45, addCoef.col = "#000000",
                         number.cex = 0.8, title = "Environment Correlation Matrix", mar = c(0, 0, 2, 0),
                         col = heat_colors)
    }, "env_cor_heatmap", size_in = c(5.5, 5.5)))
  } else add_note("Correlation between environments needs the trait in at least two environments.")
}

# =============================================================================
# BLUEs: within each environment and over all environments
# =============================================================================
mod_blues <- function() {
  trait <- mod_trait()
  eb <- safe("env_blues", env_blues_table(working_data, trait, design))
  env_blues <- NULL; env_sed <- NULL
  if (!is.null(eb)) {
    env_blues <- eb$blues; env_sed <- eb$sed
    add_note(eb$notes); add_warning(eb$warnings)
    for (nm in names(eb$errors)) add_error(nm, eb$errors[[nm]])
  }
  sed_keep <<- env_sed
  if (!is.null(dres))      # the standard errors of differences belong to the design table
    write_table(tidy_design_check(dres$table, env_sed), "design_check", "Design_Check",
                "Layout of every environment, design warnings and standard errors of differences (augmented design)")
  if (!is.null(env_blues)) {
    n_g <- n_distinct(env_blues$Genotype); n_e <- n_distinct(env_blues$Environment)
    fig("env_blues_plot", write_plot_out(plot_env_blues(env_blues, trait, color_scheme), "env_blues_plot",
                                         height_px = max(600, n_g * 16 * max(1, ceiling(n_e / 3)) + 120),
                                         size_in = c(min(11, 1.5 + 3.2 * min(n_e, 3)),
                                                     min(12, max(5, 1.5 + 0.16 * n_g * ceiling(n_e / 3))))))
  }
  mm <- safe("mixed_model_analysis", across_env_blues(working_data, trait, design, isTRUE(cfg("interaction", TRUE))))
  if (!is.null(mm)) {
    mm_res <- mm$table
    if (!is.null(mm$note)) add_note(mm$note)
    fig("blues_plot", {
      if (n_distinct(mm_res$Environment) > 1) {
        write_plot_out(plot_lines(mm_res, "BLUE", paste("BLUEs for", trait), "BLUE", color_scheme), "blues_plot",
                       tooltip = "text", size_in = FIG_WIDE)
      } else {
        write_plot_out(plot_blues_ranked(mm_res, trait, color_scheme), "blues_plot",
                       height_px = max(600, n_distinct(mm_res$Genotype) * 16 + 120),
                       size_in = fig_tall(n_distinct(mm_res$Genotype)))
      }
    })
  }
  if (!is.null(env_blues) || !is.null(mm))
    write_table(tidy_blues_all(env_blues, if (!is.null(mm)) mm$table else NULL), "blues", "BLUEs",
                paste("BLUEs of", trait, ": within each environment (augmented-design or mixed model) and over all environments (mixed model)"))
}

# =============================================================================
# ANOVA: every selected trait in every season and pooled; residual diagnostics of one trait
# =============================================================================
mod_anova <- function() {
  trs <- mod_traits()
  dtr <- mod_trait("diag_trait", trs)                          # trait of the residual plots
  new_ta <- safe("trait_anova", trait_anova_all(working_data, trs, design, keep_pooled = dtr))
  pooled_fits <- list()
  if (!is.null(new_ta)) {
    pooled_fits <- new_ta$pooled; new_ta$pooled <- NULL       # fitted models are not kept in the saved state
    # one consolidated ANOVA table per trait (seasons and pooled side by side); all traits in one CSV / Excel sheet
    env_order <- unique(as.character(working_data$Environment))
    new_ta$wide <- setNames(lapply(unique(new_ta$anova$Trait), function(tr)
      anova_wide(new_ta$anova[new_ta$anova$Trait == tr, ], new_ta$summary[new_ta$summary$Trait == tr, ], env_order)),
      unique(new_ta$anova$Trait))
    write_table(anova_long_all(new_ta$wide), "anova", "ANOVA", "ANOVA of every trait", workbook = FALSE)
    ta <<- new_ta                                              # kept for the workbooks
    n_tr <- n_distinct(new_ta$partition$Trait); nc <- min(4, n_tr)
    fig("trait_variance_plot", write_plot_out(plot_variance_partition(new_ta$partition), "trait_variance_plot",
                                              height_px = 650,
                                              size_in = c(min(9.5, 1.5 + 1.9 * nc), 1.6 + 2.6 * ceiling(n_tr / 4))))
  }
  # residual diagnostics of the pooled model of one trait
  if (draw()) {
    anova_res <- pooled_fits[[dtr]] %||% safe("anova_diagnostics", pooled_anova(working_data, dtr, design))
    if (!is.null(anova_res)) {
      rd <- pooled_residuals(anova_res)
      if (!anova_res$used_gxe && n_distinct(working_data$Environment[!is.na(working_data[[dtr]])]) > 1)
        add_note("Pooled ANOVA: Genotype x Environment is not a separate term (entries are not replicated within an environment, or the layout cannot estimate it), so it is part of the residual.")
      safe("residual_plot", write_plot_out(plot_residuals(rd, color_scheme), "residual_plot", size_in = c(6.5, 5)))
      safe("qq_plot",       write_plot_out(plot_qq(rd, color_scheme), "qq_plot", size_in = c(5.5, 5.5)))
    }
  }
}

# =============================================================================
# HERITABILITY: variance components and heritability, per trait, season and pooled
# =============================================================================
mod_herit <- function() {
  trs <- mod_traits()
  h2_df <- h2_cached(trs)
  if (!is.null(h2_df)) {
    h2_df <- h2_df[order(match(h2_df$Trait, trs)), ]
    write_table(tidy_herit(h2_df), "heritability", "Heritability", "Variance components and heritability per trait and environment")
    fig("heritability_plot", write_plot_out(plot_heritability(h2_df, color_scheme), "heritability_plot", tooltip = "text",
                                            size_in = c(min(12, max(7, 1.5 + 1.0 * n_distinct(h2_df$Trait))), 5)))
    one_season <- setdiff(trs, h2_df$Trait[h2_df$Environment == "Pooled"])
    if (length(one_season))
      add_note(paste("Heritability across environments needs a trait measured in two or more environments; per-environment values only for:",
                     paste(one_season, collapse = ", ")))
    if (design == "augmented")
      add_note("Heritability: checks are fixed effects, so differences between check varieties do not inflate the genotypic variance; the error variance comes from the check plots.")
  }
}

# =============================================================================
# CORRELATION between traits (on accession means)
# =============================================================================
mod_correlation <- function() {
  mv <- get_mv()
  meth <- cfg("method", "pearson")
  cor_res <- safe("correlation_analysis", correlation_analysis(mv, meth))
  if (!is.null(cor_res)) {
    mname <- c(pearson = "Pearson", spearman = "Spearman", kendall = "Kendall")[[cor_res$method]]
    write_table(cor_res$table, "correlation_table", "Trait_Correlations",
                paste(mname, "correlation between every pair of traits, with P value"))
    if (nrow(cor_res$duplicates))
      add_note(paste0("These traits are perfectly correlated, so they carry the same information: ",
                      paste0(cor_res$duplicates[["Trait 1"]], " & ", cor_res$duplicates[["Trait 2"]], collapse = "; "),
                      ". Consider keeping only one in the multi-trait analyses."))
    fig("correlation_heatmap", write_png_out(function() {
      corrplot::corrplot(cor_res$matrix, method = "color", type = "upper", order = "hclust",
                         tl.cex = 0.8, tl.col = "#000000", tl.srt = 45, addCoef.col = "#000000",
                         number.cex = 0.7, title = paste("Trait Correlation Matrix (", mname, ")", sep = ""),
                         mar = c(0, 0, 2, 0), col = heat_colors)
    }, "correlation_heatmap", size_in = c(7.5, 7.5)))
    fig("correlation_network", write_png_out(function() {
      corrplot::corrplot(cor_res$matrix, method = "circle", type = "full", order = "hclust",
                         tl.cex = 0.8, tl.col = "#000000", title = "Correlation Circles", mar = c(0, 0, 2, 0),
                         col = heat_colors)
    }, "correlation_network", size_in = c(7.5, 7.5)))
  }
}

# =============================================================================
# PCA
# =============================================================================
mod_pca <- function() {
  mv <- get_mv()
  pca_res <- safe("pca_analysis", pca_analysis(mv, isTRUE(cfg("scale", TRUE))))
  if (!is.null(pca_res)) {
    write_table(tidy_pca(pca_res), "pca", "PCA", "Trait loadings, eigenvalues and variance explained of the principal components")
    fig("pca_plot", write_plot_out(plot_pca_biplot(pca_res, color_scheme), "pca_plot", tooltip = c("text", "x", "y"),
                                   size_in = pca_fig_size(pca_res)))
    fig("pca_variance_plot", write_plot_out(plot_pca_variance(pca_res$variance, color_scheme), "pca_variance_plot",
                                            tooltip = "text", size_in = c(7, 5)))
  }
}

# =============================================================================
# CLUSTERING of accessions
# =============================================================================
mod_cluster <- function() {
  mv <- get_mv()
  k <- suppressWarnings(as.integer(cfg("k", 3))); method <- cfg("method", "kmeans")
  clust_res <- safe("clustering_analysis", cluster_analysis(mv, k, method))
  if (is.null(clust_res)) add_note(sprintf("%d clusters are too many for %d accessions; lower the number of clusters.", k, nrow(mv$data)))
  else {
    pca_res <- safe("pca_for_clusters", pca_analysis(mv, TRUE))     # the table lists the PC scores of every accession
    if (!is.null(pca_res))
      write_table(tidy_cluster_members(clust_res, pca_res, mv), "clusters", "Clusters",
                  "Cluster means, then every accession with its cluster, principal component scores and trait means")
    fig("cluster_plot", write_plot_out(plot_clusters(clust_res, k, color_scheme), "cluster_plot",
                                       tooltip = c("text", "color"), size_in = FIG_SQUARE))
    if (!is.null(clust_res$hclust))
      fig("dendrogram", write_png_out(function() {
        plot(clust_res$hclust, main = "Hierarchical Clustering Dendrogram", xlab = "Accessions", ylab = "Distance",
             cex = 0.7, col = "#000000")
        rect.hclust(clust_res$hclust, k = k, border = get_color_palette(color_scheme, 2)[2])
      }, "dendrogram", size_in = c(max(8, min(14, 0.14 * nrow(clust_res$data))), 5.5)))
  }
}

# =============================================================================
# STABILITY and AMMI (always the FULL data, not environment-filtered)
# =============================================================================
mod_stability <- function() {
  trait <- mod_trait()
  stab_res <- safe("stability_analysis", stability_analysis(data, trait))
  if (is.null(stab_res)) add_note(paste("Stability needs", trait, "in at least two environments and at least three accessions present in all of them."))
  else {
    if (stab_res$n_excluded > 0)
      add_note(sprintf("Stability and AMMI use the %d genotype(s) present in every environment; %d genotype(s) missing from at least one environment were left out.",
                       nrow(stab_res$parameters), stab_res$n_excluded))
    write_table(tidy_stability(stab_res), "stability", "Stability_GxE",
                paste("Stability of", trait, ": Wricke ecovalence, regression slope, AMMI scores and", GXE, "effects; sorted by stability rank"))
    fig("stability_plot", write_plot_out(plot_stability(stab_res$parameters, trait, color_scheme), "stability_plot",
                                         tooltip = c("text", "x", "y"), size_in = c(7, 5.5)))
    fig("ammi_plot", write_plot_out(
      plot_lines(stab_res$ammi_data, "Mean_Value", paste("G×E Interaction —", trait), "Mean Value", color_scheme),
      "ammi_plot", tooltip = "text", size_in = FIG_WIDE))
    if (!is.null(stab_res$ammi_scores))
      fig("ammi_biplot", write_plot_out(plot_ammi_biplot(stab_res$ammi_scores, stab_res$ipca_var, trait, color_scheme),
                                        "ammi_biplot", tooltip = c("text", "x", "y"), size_in = c(6.5, 5.5)))
  }
}

MODULE_FUNS <- list(descriptive = mod_descriptive, environments = mod_environments, blues = mod_blues, anova = mod_anova,
                    herit = mod_herit, correlation = mod_correlation, pca = mod_pca, cluster = mod_cluster,
                    stability = mod_stability)
ALL_MODULES <- names(MODULE_FUNS)

# What each module writes (table files, workbook sheets, figures). A module that runs again first removes its old
# output, so a folder never holds files that belong to settings that are no longer the current ones.
MODULE_OUTPUTS <- list(
  descriptive  = list(tables = "descriptive_stats", sheets = "Descriptive_Statistics",
                      plots = c("env_boxplot", "env_violin_plot", "distribution_plot")),
  environments = list(tables = c("env_means", "env_cor"), sheets = c("Env_Means_Ranks", "Env_Correlations"),
                      plots = c("env_means_plot", "env_rank_plot", "env_cor_heatmap")),
  blues        = list(tables = "blues", sheets = "BLUEs", plots = c("env_blues_plot", "blues_plot")),
  anova        = list(tables = "anova", sheets = character(0), plots = c("trait_variance_plot", "residual_plot", "qq_plot")),
  herit        = list(tables = "heritability", sheets = "Heritability", plots = "heritability_plot"),
  correlation  = list(tables = "correlation_table", sheets = "Trait_Correlations",
                      plots = c("correlation_heatmap", "correlation_network")),
  pca          = list(tables = "pca", sheets = "PCA", plots = c("pca_plot", "pca_variance_plot")),
  cluster      = list(tables = "clusters", sheets = "Clusters", plots = c("cluster_plot", "dendrogram")),
  stability    = list(tables = "stability", sheets = "Stability_GxE", plots = c("stability_plot", "ammi_plot", "ammi_biplot")))