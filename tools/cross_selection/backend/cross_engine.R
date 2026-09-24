# ============================================================
# TerraTools — Cross Selection ENGINE
# Ported from the standalone Shiny "Half-Diallel Cross
# Optimizer" prototype into the TerraTools CLI pattern
# (see qtl_engine.R for the sibling QTL module).
#
# Notes on the port:
#  - reshape2 was imported by the original Shiny app but never
#    called (no melt()/dcast() anywhere) — dropped here.
#  - DT/shiny were UI-only and are dropped for the same reason
#    qtl_engine.R dropped shiny/DT/plotly/colourpicker/shinyjs.
#  - ggplot2 output now renders to PNG via the cli driver
#    instead of Shiny's renderPlot(), so plot dimensions and
#    colors are threaded through as explicit parameters instead
#    of being fixed in a reactive dark_theme().
# ============================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

required_packages <- c(
  "ggplot2", "dplyr", "scales", "jsonlite"
)

missing_packages <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them once ahead of time (e.g. in your Docker image / setup script) with:\n",
    "  Rscript -e 'install.packages(c(",
    paste0('"', missing_packages, '"', collapse = ", "),
    "))'\n",
    "Live install.packages() inside a request is disabled because it can hang ",
    "indefinitely (e.g. no outbound internet, no CRAN binary) with no visible error."
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

# ============================================================
# RELATEDNESS / PERFORMANCE CATEGORY COLORS
# ============================================================
RELATEDNESS_COLORS <- c(
  "Unrelated"          = "#3A7D44",
  "Distantly_Related"  = "#C8760A",
  "Moderately_Related" = "#E8834A",
  "Closely_Related"    = "#990000"
)

# ============================================================
# HELPERS
# ============================================================
cn <- function(x) trimws(gsub('[\r\n"]', "", as.character(x)))

get_color <- function(input_val, fallback) {
  if (is.null(input_val) || nchar(trimws(input_val)) == 0) return(fallback)
  input_val
}

# ============================================================
# DATA LOADERS
# ============================================================
load_data <- function(file_path, expected_cols = NULL) {
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  if (!is.null(expected_cols) && !all(expected_cols %in% colnames(data))) {
    stop(paste(
      "Missing required columns. Expected:",
      paste(expected_cols, collapse = ", ")
    ))
  }
  data
}

process_kmat <- function(kmat_raw) {
  Kmat <- as.matrix(kmat_raw[, -1])
  rownames(Kmat) <- trimws(kmat_raw[, 1])
  colnames(Kmat) <- rownames(Kmat)
  Kmat[Kmat < 0] <- 1e-6
  Kmat <- (Kmat + t(Kmat)) / 2
  diag(Kmat) <- 0
  Kmat
}

# ============================================================
# RELATEDNESS THINNING
# Greedily drops the lowest-criterion members of any cluster of
# genotypes that are more related to each other than `threshold`,
# keeping at most `max_per_cluster` members per cluster.
# ============================================================
relate_thinning <- function(K, criterion, threshold = 0.15, max_per_cluster = 5) {
  K[K < 0] <- 1e-6
  diag(K) <- 0
  K_binary <- K
  K_binary[K_binary > threshold] <- 1
  K_binary[K_binary <= threshold] <- 0

  for (i in 1:nrow(K_binary)) {
    index <- colSums(K_binary)
    index <- sort(index, decreasing = TRUE)
    index <- index[which(index > max_per_cluster)]
    if (length(index) == 0) break
    cluster_members <- names(which(K_binary[names(index)[1], ] == 1))
    tmp <- criterion[match(cluster_members, criterion[, 1]), ]
    tmp <- tmp[order(tmp$Criterion, decreasing = TRUE), ]
    if (nrow(tmp) > max_per_cluster) {
      tmp_rm <- match(tmp[-c(1:max_per_cluster), 1], rownames(K_binary))
      tmp_rm <- tmp_rm[!is.na(tmp_rm)]
      if (length(tmp_rm) > 0) K_binary <- K_binary[-tmp_rm, -tmp_rm]
    }
  }
  rownames(K_binary)
}

# ============================================================
# HALF-DIALLEL CROSS GENERATION
# One row per unordered parent pair. Predicted_Offspring_Value
# is the simple mid-parent value; Kinship_Coefficient comes
# straight from the (thinned) kinship matrix.
# ============================================================
generate_half_diallel <- function(parents, criterion, K) {
  n <- length(parents)
  if (n < 2) stop("Need at least 2 parents")
  total_crosses <- n * (n - 1) / 2

  P1 <- character(total_crosses)
  P2 <- character(total_crosses)
  P1_Value <- numeric(total_crosses)
  P2_Value <- numeric(total_crosses)
  Mid_Parent_Value <- numeric(total_crosses)
  Kinship <- numeric(total_crosses)

  idx <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      P1[idx] <- parents[i]
      P2[idx] <- parents[j]
      P1_Value[idx] <- criterion$Criterion[criterion$Genotype == parents[i]]
      P2_Value[idx] <- criterion$Criterion[criterion$Genotype == parents[j]]
      Mid_Parent_Value[idx] <- (P1_Value[idx] + P2_Value[idx]) / 2
      Kinship[idx] <- K[parents[i], parents[j]]
      idx <- idx + 1
    }
  }

  crosses <- data.frame(
    Cross_ID = 1:total_crosses,
    Parent1 = P1,
    Parent2 = P2,
    P1_Criterion_Value = round(P1_Value, 4),
    P2_Criterion_Value = round(P2_Value, 4),
    Predicted_Offspring_Value = round(Mid_Parent_Value, 4),
    Kinship_Coefficient = round(Kinship, 6),
    stringsAsFactors = FALSE
  )

  crosses$Relatedness_Category <- as.character(cut(
    crosses$Kinship_Coefficient,
    breaks = c(-Inf, 0.05, 0.125, 0.25, Inf),
    labels = c("Unrelated", "Distantly_Related", "Moderately_Related", "Closely_Related")
  ))

  q_breaks <- quantile(crosses$Predicted_Offspring_Value, probs = c(0, 0.25, 0.5, 0.75, 1))
  if (length(unique(q_breaks)) < 5) {
    crosses$Performance_Category <- "Medium"
  } else {
    crosses$Performance_Category <- as.character(cut(
      crosses$Predicted_Offspring_Value,
      breaks = q_breaks,
      labels = c("Low", "Medium", "High", "Very_High"),
      include.lowest = TRUE
    ))
  }

  crosses <- crosses[order(crosses$Predicted_Offspring_Value, decreasing = TRUE), ]
  rownames(crosses) <- NULL
  crosses
}

build_parent_usage <- function(crosses, criterion) {
  all_parents <- sort(unique(c(crosses$Parent1, crosses$Parent2)))
  usage <- data.frame(
    Genotype = all_parents,
    Times_Used = sapply(
      all_parents,
      function(x) sum(crosses$Parent1 == x | crosses$Parent2 == x)
    ),
    Criterion_Value = sapply(
      all_parents,
      function(x) criterion$Criterion[criterion$Genotype == x][1]
    ),
    stringsAsFactors = FALSE
  )
  usage[order(usage$Criterion_Value, decreasing = TRUE), ]
}

# ============================================================
# PLOT THEME — TerraTools palette, PNG-oriented (no reactive
# Shiny bindings; every visual knob comes in through plot_cfg)
# ============================================================
terratools_theme <- function(bg_color = "#FFFFFF", panel_color = "#F7F7F5",
                              grid_color = "#E3EAE2", text_color = "#4A534A",
                              title_color = "#1A1A1A") {
  theme_minimal(base_family = "sans") +
    theme(
      plot.background   = element_rect(fill = bg_color, color = NA),
      panel.background  = element_rect(fill = panel_color, color = NA),
      panel.grid.major  = element_line(color = grid_color, linewidth = 0.4),
      panel.grid.minor  = element_blank(),
      axis.text         = element_text(color = text_color, size = 9),
      axis.title        = element_text(color = text_color, size = 10),
      plot.title        = element_text(color = title_color, size = 12, face = "bold"),
      legend.background = element_rect(fill = bg_color, color = NA),
      legend.text       = element_text(color = text_color, size = 8),
      legend.title      = element_text(color = text_color, size = 8),
      strip.background  = element_rect(fill = panel_color, color = NA),
      strip.text        = element_text(color = text_color)
    )
}

# ============================================================
# PLOTS
# ============================================================
make_pov_histogram <- function(crosses, bar_color = "#3A7D44", mean_line_color = "#1A4A1F", bins = 40) {
  ggplot(crosses, aes(x = Predicted_Offspring_Value)) +
    geom_histogram(bins = bins, fill = bar_color, alpha = 0.85, color = "#FFFFFF") +
    geom_vline(xintercept = mean(crosses$Predicted_Offspring_Value),
               color = mean_line_color, linetype = "dashed", linewidth = 0.8) +
    labs(x = "Predicted Offspring Value", y = "Count", title = "Predicted Offspring Value Distribution") +
    terratools_theme()
}

make_kinship_histogram <- function(crosses, bar_color = "#C8760A", mod_line_color = "#E8834A",
                                    close_line_color = "#990000", bins = 40) {
  ggplot(crosses, aes(x = Kinship_Coefficient)) +
    geom_histogram(bins = bins, fill = bar_color, alpha = 0.85, color = "#FFFFFF") +
    geom_vline(xintercept = 0.125, color = mod_line_color, linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = 0.25, color = close_line_color, linetype = "dashed", linewidth = 0.8) +
    labs(x = "Kinship Coefficient (K)", y = "Count", title = "Kinship Coefficient Distribution") +
    terratools_theme()
}

make_relatedness_bar <- function(crosses, colors = RELATEDNESS_COLORS) {
  cnt <- as.data.frame(table(crosses$Relatedness_Category))
  colnames(cnt) <- c("Category", "Count")
  ggplot(cnt, aes(x = reorder(Category, -Count), y = Count, fill = Category)) +
    geom_col(alpha = 0.9, width = 0.6) +
    scale_fill_manual(values = colors, guide = "none") +
    labs(x = "", y = "Number of Crosses", title = "Crosses by Relatedness Category") +
    terratools_theme() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}

make_value_kinship_scatter <- function(crosses, colors = RELATEDNESS_COLORS) {
  ggplot(crosses, aes(x = Kinship_Coefficient, y = Predicted_Offspring_Value, color = Relatedness_Category)) +
    geom_point(alpha = 0.4, size = 1.1) +
    geom_vline(xintercept = 0.125, color = "#C8760A", linetype = "dashed", linewidth = 0.6, alpha = 0.7) +
    scale_color_manual(values = colors, name = "Relatedness") +
    labs(x = "Kinship (K)", y = "Predicted Offspring Value", title = "Value vs Kinship") +
    terratools_theme() +
    theme(legend.position = "bottom")
}

make_parent_bar <- function(parent_usage, bar_color = "#3A7D44") {
  parent_usage$Genotype <- factor(
    parent_usage$Genotype,
    levels = parent_usage$Genotype[order(parent_usage$Criterion_Value)]
  )
  ggplot(parent_usage, aes(x = Criterion_Value, y = Genotype)) +
    geom_col(fill = bar_color, alpha = 0.85, width = 0.7) +
    labs(x = "Criterion Value", y = "", title = "Parent Criterion Values (Ranked)") +
    terratools_theme() +
    theme(axis.text.y = element_text(size = 7))
}
