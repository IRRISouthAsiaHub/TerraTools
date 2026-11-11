# =============================================================================
# 5. PLOT BUILDERS (ggplot objects; run_phenotypic turns them into HTML + PNG)
# =============================================================================
ACCENT <- c("#1E3A8A", "#DC143C", "#FF8C00", "#2E8B57", "#8A2BE2")   # colours for highlighted genotypes
GREY   <- "#B8BEC6"                                                   # everything that is not highlighted

# Blue-white-red palette for correlations (the same in every colour scheme, so 0 is always neutral)
CORR_COLORS <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(200)

SOURCE_COLORS <- c("Environment" = "#F59E0B", "Blocks/Reps" = "#94A3B8", "Genotype" = "#1E3A8A",
                   "G x E" = "#10B981", "Residual" = "#E5E7EB")

# ── Distributions per environment ────────────────────────────────────────────
plot_env_box <- function(df, trait, scheme) {
  n <- nlevels(droplevels(df$Environment))
  ggplot(df, aes(x = Environment, y = Value, fill = Environment)) +
    geom_boxplot(alpha = 0.55, width = 0.5, outlier.shape = NA) +          # outliers are drawn by the points
    geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1), alpha = 0.5, size = 1.4, color = "#333333") +   # sideways only, fixed seed: identical figure every run
    labs(title = paste("Boxplot of", trait, "per Environment"), x = NULL, y = trait) +
    theme_pheno(label_geom(n), "none") +
    scale_fill_manual(values = get_color_palette(scheme, n)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 12))
}

plot_env_violin <- function(df, trait, scheme) {
  n <- nlevels(droplevels(df$Environment))
  ggplot(df, aes(x = Environment, y = Value, fill = Environment)) +
    geom_violin(alpha = 0.55, trim = TRUE, adjust = 1.5) +                  # smoother, no tails beyond the data
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.9, outlier.shape = NA) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3.5, color = "#B22222") +
    labs(title = paste("Violin Plot of", trait, "per Environment (red diamond = mean)"), x = NULL, y = trait) +
    theme_pheno(label_geom(n), "none") +
    scale_fill_manual(values = get_color_palette(scheme, n)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 12))
}

# Genotypes sorted by overall mean; one dot per environment joined by a grey bar
plot_env_means <- function(env_means, trait, scheme) {
  ord <- env_means %>% group_by(Genotype) %>% summarise(m = mean(Mean), .groups = "drop") %>% arrange(m)
  d <- env_means %>% mutate(Genotype = factor(Genotype, levels = ord$Genotype))
  seg <- d %>% group_by(Genotype) %>% summarise(lo = min(Mean), hi = max(Mean), .groups = "drop")
  ggplot() +
    geom_segment(data = seg, aes(x = lo, xend = hi, y = Genotype, yend = Genotype), color = GREY, linewidth = 0.7) +
    geom_point(data = d, aes(x = Mean, y = Genotype, color = Environment), size = 2.4) +
    labs(title = paste("Genotype Means per Environment —", trait), x = paste("Mean", trait), y = NULL) +
    theme_pheno() +
    theme(axis.text.y = element_text(size = text_size_for(n_distinct(d$Genotype)))) +
    scale_color_manual(values = get_color_palette(scheme, n_distinct(d$Environment))) +
    scale_y_discrete(labels = function(x) stringr::str_trunc(x, 22))
}

# Genotype lines across environments. With more than 12 genotypes only the `top_n` highest (or best
# ranked) are coloured and labelled; the others stay grey, so the picture stays readable.
plot_lines <- function(df, y, title, ylab, scheme, reverse_y = FALSE, top_n = 5) {
  envs <- unique(as.character(df$Environment))
  df$Environment <- factor(df$Environment, levels = envs)
  n_geno <- n_distinct(df$Genotype)
  df$.tip <- paste0("Genotype: ", df$Genotype, "<br>", ylab, ": ", round(df[[y]], 2))
  th <- theme_pheno(label_geom(length(envs)), if (n_geno <= 12) "right" else "bottom")
  if (n_geno <= 12) {
    p <- ggplot(df, aes(x = Environment, y = .data[[y]], group = Genotype, color = Genotype, text = .tip)) +
      geom_line(linewidth = 1) + geom_point(size = 2.5) +
      scale_color_manual(values = get_color_palette(scheme, n_geno))
  } else {
    m <- tapply(df[[y]], df$Genotype, mean, na.rm = TRUE)
    top <- names(sort(m, decreasing = !reverse_y))[seq_len(min(top_n, length(m)))]
    df$.grp <- factor(ifelse(df$Genotype %in% top, as.character(df$Genotype), "Other"), levels = c(top, "Other"))
    hi <- df[df$.grp != "Other", ]; lo <- df[df$.grp == "Other", ]
    last <- hi[hi$Environment == tail(envs, 1), ]
    p <- ggplot(df, aes(x = Environment, y = .data[[y]], group = Genotype, text = .tip)) +
      geom_line(data = lo, color = GREY, alpha = 0.6, linewidth = 0.5) +
      geom_point(data = lo, color = GREY, alpha = 0.6, size = 1.2) +
      geom_line(data = hi, aes(color = .grp), linewidth = 1.2) +
      geom_point(data = hi, aes(color = .grp), size = 2.8) +
      geom_text(data = last, aes(label = Genotype, color = .grp), hjust = -0.12, size = 3.3,
                check_overlap = TRUE, show.legend = FALSE) +
      scale_color_manual(values = setNames(ACCENT[seq_along(top)], top),
                         name = if (reverse_y) paste("Best", length(top)) else paste("Highest", length(top))) +
      scale_x_discrete(expand = expansion(add = c(0.3, 1.2))) +
      guides(color = guide_legend(nrow = 1))
  }
  if (reverse_y) p <- p + scale_y_reverse()
  p + labs(title = title, x = NULL, y = ylab) + th
}

plot_env_rank <- function(df, trait, scheme)
  plot_lines(df, "Rank", paste("Genotype Rank Across Environments —", trait), "Rank (1 = best)", scheme, reverse_y = TRUE)

plot_env_blues <- function(env_blues, trait, scheme) {
  n_geno <- n_distinct(env_blues$Genotype); n_env <- n_distinct(env_blues$Environment)
  ggplot(env_blues, aes(x = reorder(Genotype, BLUE, FUN = mean), y = BLUE, color = Environment)) +
    geom_errorbar(aes(ymin = BLUE - SE, ymax = BLUE + SE), width = 0.4, alpha = 0.6) +
    geom_point(size = 2, alpha = 0.9) +
    facet_wrap(~Environment, scales = "free_x", ncol = min(n_env, 3)) + coord_flip() +
    labs(title = paste("Within-Environment BLUEs —", trait, "(bars = ±1 SE)"), x = NULL, y = "BLUE") +
    theme_pheno(legend = "none") +
    theme(axis.text.y = element_text(size = text_size_for(n_geno), lineheight = 0.9),
          strip.text = element_text(size = 11, face = "bold", color = "#1E3A8A")) +
    scale_color_manual(values = get_color_palette(scheme, n_env)) +
    scale_x_discrete(labels = function(x) stringr::str_trunc(x, 22))
}

# ── Across-environment analyses ──────────────────────────────────────────────
plot_distribution <- function(values, trait, scheme) {
  pal <- get_color_palette(scheme, 2)
  ggplot(data.frame(Value = values), aes(x = Value)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = pal[1], alpha = 0.7) +
    geom_density(color = pal[2], linewidth = 1) +
    labs(title = paste("Distribution of", trait), x = trait, y = "Density") + theme_minimal()
}

plot_residuals <- function(rd, scheme) {
  pal <- get_color_palette(scheme, 2)
  # a robust trend line: mean residual in groups of fitted values (loess breaks down on tied fitted values)
  brk <- unique(quantile(rd$Fitted, seq(0, 1, length.out = 9), names = FALSE))
  grp <- cut(rd$Fitted, breaks = brk, include.lowest = TRUE)
  trend <- data.frame(Fitted = tapply(rd$Fitted, grp, mean), Residuals = tapply(rd$Residuals, grp, mean))
  trend <- trend[stats::complete.cases(trend), ]
  ggplot(rd, aes(x = Fitted, y = Residuals, text = Genotype)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(size = 1.8, alpha = 0.65, color = pal[2]) +
    geom_line(data = trend, aes(x = Fitted, y = Residuals), inherit.aes = FALSE, color = "#DC143C", linewidth = 0.9) +
    geom_point(data = trend, aes(x = Fitted, y = Residuals), inherit.aes = FALSE, color = "#DC143C", size = 2.4) +
    labs(title = "Residuals vs Fitted (red = mean residual in 8 groups of fitted values)", x = "Fitted", y = "Residuals") +
    theme_pheno(legend = "none")
}

plot_qq <- function(rd, scheme) {
  pal <- get_color_palette(scheme, 2)
  z <- rd$Standardized_Residuals
  # quantiles are computed here: ggplotly cannot convert stat_qq() layers
  rd$Theoretical <- qnorm(ppoints(length(z)))[rank(z, ties.method = "first")]
  qs <- quantile(z, c(.25, .75), na.rm = TRUE); qt <- qnorm(c(.25, .75))
  sl <- diff(qs) / diff(qt)
  ggplot(rd, aes(x = Theoretical, y = Standardized_Residuals, text = Genotype)) +
    geom_abline(slope = sl, intercept = qs[[1]] - sl * qt[[1]], color = "#DC143C") +
    geom_point(size = 1.8, alpha = 0.7, color = pal[2]) +
    labs(title = "Normal Q-Q Plot", x = "Theoretical Quantiles", y = "Standardised residuals") + theme_pheno(legend = "none")
}

# BLUEs across environments (single group): genotypes ranked, with standard-error bars
plot_blues_ranked <- function(df, trait, scheme) {
  pal <- get_color_palette(scheme, 2)
  ggplot(df, aes(x = reorder(Genotype, BLUE), y = BLUE)) +
    geom_errorbar(aes(ymin = BLUE - SE, ymax = BLUE + SE), width = 0.4, alpha = 0.6, color = pal[2]) +
    geom_point(size = 2, color = pal[2], alpha = 0.9) +
    coord_flip() + labs(title = paste("BLUEs across environments —", trait, "(bars = ±1 SE)"), x = NULL, y = "BLUE") +
    theme_pheno() + theme(axis.text.y = element_text(size = text_size_for(n_distinct(df$Genotype)))) +
    scale_x_discrete(labels = function(x) stringr::str_trunc(x, 22))
}

# Entry-mean heritability of every trait: one bar per season, and one for all seasons pooled
plot_heritability <- function(h2_df, scheme) {
  envs <- unique(h2_df$Environment); envs <- c(setdiff(envs, "Pooled"), intersect(envs, "Pooled"))
  h2_df$Environment <- factor(h2_df$Environment, levels = envs)
  h2_df$Trait <- factor(h2_df$Trait, levels = unique(h2_df$Trait))
  seasons <- setdiff(envs, "Pooled")
  cols <- c(setNames(get_color_palette(scheme, max(length(seasons), 2))[seq_along(seasons)], seasons), Pooled = "#1E3A8A")
  ggplot(h2_df, aes(x = Trait, y = H2_Mean, fill = Environment,
                    text = paste0(Trait, " (", Environment, "): H2 = ", sprintf("%.2f", H2_Mean)))) +
    geom_col(position = position_dodge(0.85), width = 0.75, alpha = 0.9) +
    geom_text(aes(label = sprintf("%.2f", H2_Mean)), position = position_dodge(0.85), vjust = -0.4, size = 3) +
    geom_hline(yintercept = c(.3, .6), linetype = "dashed", color = "grey50") +
    scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = cols, name = NULL) +
    labs(title = "Entry-mean heritability per trait (dashed lines: 0.3 low / 0.6 high)", x = NULL, y = "Heritability (H2)") +
    theme_pheno(legend = "bottom")
}

# Where the variance comes from: share of the total sum of squares, per trait and season
plot_variance_partition <- function(part) {
  envs <- unique(part$Environment); envs <- c(setdiff(envs, "Pooled"), intersect(envs, "Pooled"))
  part$Environment <- factor(part$Environment, levels = envs)
  part$Source <- factor(part$Source, levels = names(SOURCE_COLORS))
  part$Trait <- factor(part$Trait, levels = unique(part$Trait))
  txt <- setNames(rep("#111111", length(SOURCE_COLORS)), names(SOURCE_COLORS)); txt["Genotype"] <- "white"
  ggplot(part, aes(x = Environment, y = Percent, fill = Source)) +
    geom_col(width = 0.75, color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(Percent >= 8, paste0(round(Percent), "%"), ""), color = Source),
              position = position_stack(vjust = 0.5), size = 3) +
    facet_wrap(~Trait, ncol = min(4, n_distinct(part$Trait)), scales = "free_x") +
    scale_fill_manual(values = SOURCE_COLORS, drop = TRUE, name = "Source of variation") +
    scale_color_manual(values = txt, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(title = "Share of total variation (sum of squares) per trait and season", x = NULL, y = "% of total SS") +
    theme_pheno(legend = "bottom") +
    theme(strip.text = element_text(face = "bold", size = 11, color = "#1E3A8A"),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

# ── Multivariate ─────────────────────────────────────────────────────────────
# PCA biplot: genotype scores plus trait loading arrows (traits that coincide share one label)
plot_pca_biplot <- function(pca, scheme) {
  pal <- get_color_palette(scheme, 2)
  sc <- pca$scores; ld <- pca$loadings
  mult <- 0.8 * min(max(abs(sc$PC1)) / max(abs(ld$PC1)), max(abs(sc$PC2)) / max(abs(ld$PC2)))
  ld <- ld %>% mutate(x = PC1 * mult, y = PC2 * mult)
  # traits whose arrows coincide share one label; the label sits just beyond the arrow tip, on the side it points to
  lab <- ld %>% group_by(kx = round(x / max(abs(x)), 1), ky = round(y / max(abs(y)), 1)) %>%
    summarise(x = mean(x), y = mean(y), label = paste(Trait, collapse = " / "), .groups = "drop") %>%
    mutate(len = sqrt(x^2 + y^2), ux = x / len, uy = y / len,
           lx = x + 0.04 * max(abs(sc$PC1)) * ux, ly = y + 0.04 * max(abs(sc$PC2)) * uy,
           hj = 0.5 - 0.6 * ux, vj = 0.5 - 0.6 * uy)
  # name only the accessions farthest from the centre (the interesting ones); the rest stay as dots.
  # A name that would run into a trait label is left out (widths are estimated from the number of characters).
  xspan <- diff(range(c(sc$PC1, ld$x))) * 1.36; yspan <- diff(range(c(sc$PC2, ld$y))) * 1.24
  far <- sc[order(-(sc$PC1^2 + sc$PC2^2)), ][seq_len(min(6, nrow(sc))), ]
  hits <- vapply(seq_len(nrow(far)), function(i) {
    w <- nchar(as.character(far$Genotype[i])) * 0.0055 * xspan; y0 <- far$PC2[i] + 0.04 * yspan
    any(abs(far$PC1[i] - (lab$lx + (0.5 - lab$hj) * nchar(lab$label) * 0.016 * xspan)) < w + nchar(lab$label) * 0.008 * xspan &
        abs(y0 - (lab$ly + (0.5 - lab$vj) * 0.05 * yspan)) < 0.06 * yspan)
  }, logical(1))
  far <- far[!hits, , drop = FALSE]
  pct <- round(pca$variance$Prop_Variance[1:2] * 100, 1)
  ggplot() +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
    suppressWarnings(geom_point(data = sc, aes(x = PC1, y = PC2, text = Genotype),    # `text` is read by plotly
                                size = 2.2, alpha = 0.7, color = pal[2])) +
    geom_text(data = far, aes(x = PC1, y = PC2, label = Genotype), size = 2.6, color = "grey30",
              vjust = -0.9, check_overlap = TRUE) +
    geom_segment(data = ld, aes(x = 0, y = 0, xend = x, yend = y),
                 arrow = arrow(length = unit(0.2, "cm")), color = "#DC143C", linewidth = 0.7) +
    geom_text(data = lab, aes(x = lx, y = ly, label = label, hjust = hj, vjust = vj),
              color = "#DC143C", size = 3.6, fontface = "bold") +
    scale_x_continuous(expand = expansion(mult = 0.18)) + scale_y_continuous(expand = expansion(mult = 0.12)) +
    coord_fixed() +                                            # same scale on both axes, so angles and distances are true
    labs(title = "PCA biplot (accession means)", x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)")) +
    theme_pheno() + theme(panel.grid.minor = element_blank())
}

# Figure size (inches) that follows the shape of the biplot, so the equal-scale axes leave no empty bands
pca_fig_size <- function(pca) {
  sc <- pca$scores; ld <- pca$loadings
  aspect <- diff(range(sc$PC2)) * 1.24 / (diff(range(sc$PC1)) * 1.36)
  c(7, min(7.5, max(4.5, 1.6 + 5.3 * aspect)))
}

plot_pca_variance <- function(variance, scheme) {
  # Bars: variance explained by each component (dark = eigenvalue >= 1, the Kaiser rule); line: cumulative variance.
  # Dashed guides: 80 % of the variance, and the share an average trait carries (eigenvalue = 1).
  p_all <- nrow(variance)
  vd <- variance[1:min(10, p_all), ]
  vd$PC <- factor(vd$PC, levels = vd$PC)
  vd$Pct <- vd$Prop_Variance * 100; vd$Cum <- vd$Cumulative * 100
  vd$Rule <- factor(ifelse(vd$Variance >= 1, "Eigenvalue ≥ 1", "Eigenvalue < 1"), levels = c("Eigenvalue ≥ 1", "Eigenvalue < 1"))
  vd$Tip <- sprintf("%s: eigenvalue %.2f, %.1f%% (cumulative %.1f%%)", vd$PC, vd$Variance, vd$Pct, vd$Cum)
  pal <- get_color_palette(scheme, 2)
  n80 <- which(vd$Cum >= 80)[1]; n_keep <- sum(variance$Variance >= 1)
  sub <- sprintf("%d component%s with eigenvalue ≥ 1 (together %.0f%% of the variance); %s reach 80%%",
                 n_keep, if (n_keep == 1) "" else "s", sum(variance$Prop_Variance[variance$Variance >= 1]) * 100,
                 if (is.na(n80)) "more than 10 components are needed to" else sprintf("%d component%s", n80, if (n80 == 1) "" else "s"))
  nx <- nlevels(vd$PC)
  ggplot(vd, aes(x = PC)) +
    geom_hline(yintercept = 80, linetype = "dashed", color = "grey55") +
    geom_hline(yintercept = 100 / p_all, linetype = "dotted", color = "grey55") +
    suppressWarnings(geom_col(aes(y = Pct, fill = Rule, text = Tip), width = 0.68)) +      # `text` is read by plotly
    annotate("text", x = 0.55, y = 80, label = "80% of the variance", hjust = 0, vjust = -0.6, size = 3.1, color = "grey35") +
    annotate("text", x = nx + 0.45, y = 100 / p_all, label = "average trait (eigenvalue = 1)", hjust = 1, vjust = -0.6,
             size = 3.1, color = "grey35") +
    geom_text(data = vd[vd$Pct >= 10, ], aes(y = Pct, label = sprintf("%.1f%%", Pct),
                                             colour = ifelse(Rule == "Eigenvalue ≥ 1", "white", "#222222")),
              vjust = 1.6, size = 3.4, fontface = "bold") + scale_colour_identity() +  # tall bars: value inside
    geom_text(data = vd[vd$Pct < 10, ], aes(y = Pct, label = sprintf("%.1f%%", Pct)), vjust = -0.55, size = 3.3,
              color = "#222222") +                                                   # short bars: value above
    geom_line(aes(y = Cum, group = 1), color = pal[2], linewidth = 0.9) +
    suppressWarnings(geom_point(aes(y = Cum, text = Tip), color = pal[2], size = 2.4)) +
    geom_text(data = vd[-1, ], aes(y = Cum, label = sprintf("%.0f%%", Cum)), vjust = -1.1, size = 3.1, color = pal[2], fontface = "bold") +
    scale_fill_manual(values = c("Eigenvalue ≥ 1" = pal[1], "Eigenvalue < 1" = "#C5CCD8"), name = NULL) +
    scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 20), labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
    labs(title = "Variance explained by the principal components", subtitle = sub,
         x = NULL, y = "Variance explained (bars) and cumulative (line)") +
    theme_pheno(legend = "bottom") +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          plot.subtitle = element_text(size = 9.5, color = "grey30"))
}

plot_clusters <- function(clust, k, scheme) {
  pf <- prcomp(clust$data)
  pct <- round(pf$sdev^2 / sum(pf$sdev^2) * 100, 1)
  d <- data.frame(PC1 = pf$x[, 1], PC2 = pf$x[, 2], Genotype = clust$clusters$Genotype, Cluster = clust$clusters$Cluster)
  ggplot(d, aes(x = PC1, y = PC2, color = Cluster, text = Genotype)) +
    geom_point(size = 3, alpha = 0.85) +
    labs(title = "Clustering (PCA space)", x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)")) +
    theme_minimal() + scale_color_manual(values = get_color_palette(scheme, k))
}

# ── Stability ────────────────────────────────────────────────────────────────
# Genotypes at or below the median ecovalence are "stable"; the five most stable are labelled
plot_stability <- function(sp, trait, scheme) {
  my <- median(sp$Ecovalence); mx <- median(sp$Overall_Mean)
  sp$Group <- factor(ifelse(sp$Ecovalence <= my, "Stable (Wi at or below median)", "Less stable"),
                     levels = c("Stable (Wi at or below median)", "Less stable"))
  lab <- sp %>% filter(Group == levels(Group)[1]) %>% arrange(desc(Overall_Mean)) %>% head(6)   # spread along the x axis
  ggplot(sp, aes(x = Overall_Mean, y = Ecovalence, color = Group, text = Genotype)) +
    geom_hline(yintercept = my, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = mx, linetype = "dashed", color = "grey60") +
    geom_point(size = 2.8, alpha = 0.85) +
    geom_text(data = lab, aes(label = Genotype), vjust = -0.9, size = 3.2, check_overlap = TRUE, show.legend = FALSE) +
    scale_color_manual(values = c("Stable (Wi at or below median)" = ACCENT[1], "Less stable" = GREY), name = NULL) +
    labs(title = paste("Stability —", trait, "(lower Wi = more stable; labels: stable genotypes with the highest means)"),
         x = "Mean performance", y = "Ecovalence (Wricke Wi)") +
    theme_pheno(legend = "bottom")
}

# AMMI1 biplot: mean performance (x) against the first interaction axis (y)
plot_ammi_biplot <- function(sc, ipca_var, trait, scheme) {
  pal <- get_color_palette(scheme, 2)
  gen <- sc[sc$Type == "Genotype", ]; env <- sc[sc$Type == "Environment", ]
  ggplot(sc, aes(x = Mean, y = IPCA1, color = Type, text = Name)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = mean(gen$Mean), linetype = "dashed", color = "grey60") +
    geom_point(data = gen, size = 2.5, alpha = 0.8) +
    geom_point(data = env, size = 4, shape = 17) +
    geom_text(data = env, aes(label = Name), vjust = -1, size = 3.5, show.legend = FALSE) +
    labs(title = paste("AMMI1 Biplot —", trait), x = "Mean performance",
         y = paste0("IPCA1 (", round(ipca_var[1] * 100, 1), "% of G×E)")) +
    scale_x_continuous(expand = expansion(mult = 0.08)) + scale_y_continuous(expand = expansion(mult = 0.15)) +
    theme_minimal() + scale_color_manual(values = c(Genotype = pal[2], Environment = "#DC143C"), name = NULL)
}
