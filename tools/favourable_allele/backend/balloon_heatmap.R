# =====================================================================
# Favourable-allele BALLOON HEATMAP  (TerraTools · Favourable Allele Analysis)
#
# This is the script "rice_fixed (2).R" supplied with the Heatmap_Final folder:
# the data preparation, the quality control and the balloon heatmap (section 6)
# are unchanged. Only these things differ, so that TerraTools can run it:
#   * the input, the output folder and the figure size come from the command line
#     (Rscript balloon_heatmap.R <input.csv> <outdir> <width_mm> <height_mm> <dpi>)
#   * the input is a CSV with the same two header rows as the Excel sheet
#     (TerraTools writes it from the genes and lines chosen on the page)
#   * the bar chart and the clustered heatmap (sections 7 and 8) are not drawn
#   * library(ggnewscale) is not loaded (it was never used in the balloon plot)
#   * trait categories the original palette did not know ("Hybrid fertility",
#     "Plant growth and development", ...) get their own colour instead of grey
#
# Everything below is the original code.
# =====================================================================

args <- commandArgs(trailingOnly = TRUE)
INPUT   <- args[1]
OUTDIR  <- args[2]
FIG_W   <- as.numeric(args[3])
FIG_H   <- as.numeric(args[4])
FIG_DPI <- as.numeric(args[5])
# name of the block of non-landrace lines (the original calls it "Improved varieties")
LAB_VAR <- if (length(args) >= 6 && nzchar(args[6])) args[6] else "Improved varieties"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
setwd(OUTDIR)

library(ggplot2)
library(reshape2)

# ---------------------------------------------------------------------
# 0) READ — the sheet has TWO header rows:
#       row 1 = group tag (Variety / Landrace) over the sample columns
#       row 2 = accession name
#    Metadata columns occupy the first 12 columns.
# ---------------------------------------------------------------------
all_rows <- as.data.frame(read.csv(INPUT, header = FALSE, colClasses = "character",
                                   na.strings = c("", "NA"), check.names = FALSE,
                                   stringsAsFactors = FALSE, fileEncoding = "UTF-8"))
hdr  <- all_rows[1:2, , drop = FALSE]
body <- all_rows[-(1:2), , drop = FALSE]
rownames(body) <- NULL

meta_cols <- c("Gene_ID","Gene name","Gene_label","Broad trait category",
               "Sub category","Marker_in_8k","CHROM","POS","REF","ALT",
               "Unfavourable allele","Favourable allele")
n_meta <- length(meta_cols)                       # = 12

group_tag    <- trimws(as.character(unlist(hdr[1, (n_meta + 1):ncol(hdr)])))
sample_names <- trimws(as.character(unlist(hdr[2, (n_meta + 1):ncol(hdr)])))
keep_s       <- !(is.na(sample_names) | sample_names == "" | sample_names == "NA")
sample_names <- sample_names[keep_s]
group_tag    <- group_tag[keep_s]

colnames(body) <- c(meta_cols, sample_names)
raw <- body

# drop blank / junk trailing rows that have no Gene name
gn <- trimws(as.character(raw[["Gene name"]]))
raw <- raw[!(is.na(gn) | gn == "" | gn == "NA"), , drop = FALSE]
rownames(raw) <- NULL

for (m in meta_cols) raw[[m]] <- trimws(as.character(raw[[m]]))

samples   <- sample_names
landraces <- samples[group_tag == "Landrace"]
checks    <- samples[group_tag != "Landrace"]      # "Variety" = improved varieties

# pretty axis labels (no underscores)
sample_label <- setNames(gsub("_", " ", samples), samples)

# coerce genotype columns to integer codes
for (s in samples) raw[[s]] <- suppressWarnings(as.integer(raw[[s]]))

raw$REF        <- trimws(as.character(raw$REF))
raw$ALT        <- trimws(as.character(raw$ALT))
fav_allele     <- trimws(as.character(raw[["Favourable allele"]]))
unfav_allele   <- trimws(as.character(raw[["Unfavourable allele"]]))

# ---------------------------------------------------------------------
# 1) CODE  <->  GENOTYPE STRING   (explicit requested mapping)
# ---------------------------------------------------------------------
code_label <- function(code) {
  code <- suppressWarnings(as.integer(code))
  if (is.na(code) || code == -1) return("./.")
  if (code == 0) return("0/0")
  if (code == 1) return("0/1")
  if (code == 2) return("1/1")
  "./."
}
code_to_geno <- function(code, ref, alt) {
  code <- suppressWarnings(as.integer(code))
  if (is.na(code) || code == -1) return("./.")
  if (code == 0) return(paste0(ref, "/", ref))   # 0/0 homozygous REF
  if (code == 1) return(paste0(ref, "/", alt))   # 0/1 heterozygous
  if (code == 2) return(paste0(alt, "/", alt))   # 1/1 homozygous ALT
  "./."
}

# ---------------------------------------------------------------------
# 2) QUALITY CONTROL
# ---------------------------------------------------------------------
MISS_THRESHOLD <- 0.80   # drop a marker if missing in >= 80% of accessions

in_refalt <- function(a, ref, alt) (a == ref) | (a == alt)
contradict <- !in_refalt(fav_allele,   raw$REF, raw$ALT) |
              !in_refalt(unfav_allele, raw$REF, raw$ALT)

geno_mat_all <- as.matrix(raw[, samples])
miss_frac    <- apply(geno_mat_all, 1,
                      function(r) mean(is.na(r) | r == -1))
high_missing <- miss_frac >= MISS_THRESHOLD

drop_marker <- contradict | high_missing

qc_report <- data.frame(
  Marker_in_8k     = raw[["Marker_in_8k"]],
  Gene             = raw[["Gene name"]],
  REF = raw$REF, ALT = raw$ALT,
  Favourable = fav_allele, Unfavourable = unfav_allele,
  Missing_fraction = round(miss_frac, 3),
  Contradicting    = contradict,
  Missing_ge_80pct = high_missing,
  Dropped          = drop_marker,
  stringsAsFactors = FALSE)
write.csv(qc_report, "qc_marker_filter_report.csv", row.names = FALSE)

message(sprintf("QC: %d markers in, %d dropped (%d contradicting, %d >= %.0f%% missing), %d retained across %d genes.",
                nrow(raw), sum(drop_marker), sum(contradict),
                sum(high_missing), 100 * MISS_THRESHOLD,
                sum(!drop_marker),
                length(unique(raw[["Gene name"]][!drop_marker]))))

raw          <- raw[!drop_marker, , drop = FALSE]
fav_allele   <- fav_allele[!drop_marker]
unfav_allele <- unfav_allele[!drop_marker]
rownames(raw) <- NULL

# ---------------------------------------------------------------------
# 3) Favourable-allele code per marker (homozygous favourable only)
# ---------------------------------------------------------------------
fav_code <- ifelse(fav_allele == raw$ALT, 2L,
            ifelse(fav_allele == raw$REF, 0L, NA_integer_))

is_fav <- function(code, favcode) {
  code <- suppressWarnings(as.integer(code))
  if (is.na(code) || code == -1) return(NA_integer_)   # missing
  if (code == 1) return(0L)                            # het -> not favourable
  if (is.na(favcode)) return(NA_integer_)
  as.integer(code == favcode)                          # hom fav -> 1, hom unfav -> 0
}

# ---------------------------------------------------------------------
# 4) Readable genotype + canonical-code export
# ---------------------------------------------------------------------
nuc  <- raw[, meta_cols]
codes <- raw[, meta_cols]
for (s in samples) {
  nuc[[s]]   <- mapply(code_to_geno, raw[[s]], raw$REF, raw$ALT, USE.NAMES = FALSE)
  codes[[s]] <- vapply(raw[[s]], code_label, character(1))
}
write.csv(nuc[, c("Marker_in_8k", samples)],
          "nucleotide_genotypes.csv", row.names = FALSE)
write.csv(codes[, c("Marker_in_8k", samples)],
          "canonical_codes.csv", row.names = FALSE)

# ---------------------------------------------------------------------
# 5) Gene-level favourable fraction + per-accession favourable count
# ---------------------------------------------------------------------
gl <- data.frame(g = raw[["Gene name"]], lab = trimws(raw[["Gene_label"]]),
                 stringsAsFactors = FALSE)
gl$lab[is.na(gl$lab) | gl$lab == "" | gl$lab == "NA"] <-
  gl$g[is.na(gl$lab) | gl$lab == "" | gl$lab == "NA"]
gl <- gl[!duplicated(gl$g), ]
genes      <- gl$g
gene_label <- setNames(gl$lab, gl$g)

tc <- data.frame(g = raw[["Gene name"]],
                 cat = trimws(raw[["Broad trait category"]]),
                 stringsAsFactors = FALSE)
tc <- tc[!duplicated(tc$g), ]
gene_cat   <- setNames(tc$cat, tc$g)
cat_levels <- unname(unique(gene_cat[genes]))

n_markers <- nrow(raw)
gene_mat  <- matrix(NA_real_, nrow = length(samples), ncol = length(genes),
                    dimnames = list(samples, genes))
fav_count <- setNames(integer(length(samples)), samples)

for (s in samples) {
  fav_flags <- vapply(seq_len(n_markers),
                      function(i) is_fav(raw[[s]][i], fav_code[i]),
                      integer(1))
  fav_count[s] <- sum(fav_flags == 1, na.rm = TRUE)
  for (gi in seq_along(genes)) {
    vals <- fav_flags[raw[["Gene name"]] == genes[gi]]
    vals <- vals[!is.na(vals)]
    gene_mat[s, gi] <- if (length(vals) == 0) NA_real_ else mean(vals)
  }
}

write.csv(data.frame(Variety = samples, gene_mat, check.names = FALSE),
          "gene_freq_matrix.csv", row.names = FALSE)
write.csv(data.frame(Variety = samples,
                     Favourable = as.integer(fav_count[samples]),
                     Total = n_markers),
          "locus_counts.csv", row.names = FALSE)

# ---------------------------------------------------------------------
# trait-category palette  (keys MATCH the data strings exactly)
#   Verified categories present in Gene_mapping_final_4.xlsx:
#     Abiotic stress, Biotic stress, DSR trait, Floral,
#     Grain quality, Hybrid Fertility, Nutrient Use Efficiency, Yield
# ---------------------------------------------------------------------
cat_palette <- c(
  "Abiotic stress"          = "#59A14F",
  "Biotic stress"           = "#E15759",
  "DSR trait"               = "#9C755F",
  "Floral"                  = "#F28E2B",
  "Grain quality"           = "#B07AA1",
  "Hybrid Fertility"        = "#76B7B2",
  "Nutrient Use Efficiency" = "#4E79A7",
  "Nutrient use Efficiency" = "#4E79A7",
  "Multiple stress"         = "#8C61C9",
  "Yield"                   = "#EDC948")
# safety net for any unexpected category
miss <- setdiff(cat_levels, names(cat_palette))
if (length(miss)) cat_palette[miss] <- "#BAB0AC"

# shared group colours
land_col <- "#2f6f63"   # teal       - landraces
var_col  <- "#b5532f"   # terracotta - improved varieties
box_col  <- "#9c4a2f"

# =====================================================================
# 6) BALLOON / BUBBLE HEATMAP  -> gene_balloon.png   (produced FIRST)
# =====================================================================
balloon <- melt(gene_mat, varnames = c("Variety", "Gene"), value.name = "Fav")
balloon$GeneLab <- factor(gene_label[as.character(balloon$Gene)],
                          levels = unname(gene_label[genes]))
balloon$Cat     <- factor(gene_cat[as.character(balloon$Gene)], levels = cat_levels)
balloon$GroupF  <- factor(ifelse(balloon$Variety %in% landraces,
                                 "Landraces", LAB_VAR),
                          levels = c("Landraces", LAB_VAR))
balloon$VarF    <- factor(balloon$Variety, levels = rev(samples))

bp_pts <- balloon[!is.na(balloon$Fav), ]
bp_na  <- balloon[ is.na(balloon$Fav), ]

blue_low  <- "#cfe3f5"; blue_mid <- "#2a7fc0"; blue_high <- "#04203f"

# Legend breaks for the favourable-allele fraction, and the EXACT fill colour
# for each break taken from the same blue ramp used for the dots. Defining
# these once guarantees the size-legend key dots darken in step with their
# size (no recycling / mis-ordering of the override fills).
fav_breaks <- c(0, 0.25, 0.5, 0.75, 1)
fav_ramp   <- scales::gradient_n_pal(c(blue_low, blue_mid, blue_high),
                                     values = c(0, 0.5, 1))(fav_breaks)

# strip palette — keys MATCH the data strings exactly (this is the one the
# balloon strip actually uses; previously it caused the grey tiles)
cat_palette_nat <- c(
  "Abiotic stress"          = "#009E73",
  "Biotic stress"           = "#D55E00",
  "DSR trait"               = "#8C6D31",
  "Floral"                  = "#E69F00",
  "Grain quality"           = "#CC79A7",
  "Hybrid Fertility"        = "#56B4E9",
  "Nutrient use Efficiency" = "#0072B2",
  "Nutrient Use Efficiency" = "#0072B2",
  "Multiple stress"         = "#7B4FA3",
  "Yield"                   = "#F0E442",
  # TerraTools: categories of the built-in panel that the original palette did not list
  "Hybrid fertility"             = "#56B4E9",
  "Plant growth and development" = "#999999",
  "Unspecified"                  = "#BAB0AC")
miss2 <- setdiff(cat_levels, names(cat_palette_nat))
if (length(miss2)) cat_palette_nat[miss2] <- rep(c("#44AA99", "#AA4499", "#117733", "#882255"), length.out = length(miss2))

catstrip_df <- data.frame(
  GeneLab = factor(unname(gene_label[genes]), levels = unname(gene_label[genes])),
  Cat     = factor(unname(gene_cat[genes]),   levels = cat_levels), y = 1,
  stringsAsFactors = FALSE, row.names = NULL)

have_ggtext <- requireNamespace("ggtext", quietly = TRUE)
y_label_fun <- function(v) {
  disp <- sample_label[v]; col <- ifelse(v %in% landraces, land_col, var_col)
  if (have_ggtext) sprintf("<b style='color:%s'>%s</b>", col, disp) else disp
}

pA_main <- ggplot() +
  geom_point(data = bp_pts, aes(GeneLab, VarF, size = Fav, fill = Fav),
             shape = 21, colour = "grey15", stroke = 0.3) +
  geom_point(data = bp_na, aes(GeneLab, VarF),
             shape = 4, size = 1.6, colour = "grey60", stroke = 0.7) +
  # --- Colour + size both encode favourable-allele fraction ---
  # The legend is driven by the FILL scale (ggplot always colours fill-scale
  # keys correctly from the gradient — this is why the colourbar was right).
  # We turn that fill scale into discrete keys at fav_breaks, and let SIZE
  # ride along as an override so the keys also grow. This avoids the
  # unreliable "recolour a size legend via override.aes$fill" path that left
  # every key the same blue.
  scale_fill_gradientn(colours = c(blue_low, blue_mid, blue_high),
                       values = c(0, 0.5, 1), limits = c(0, 1),
                       breaks = fav_breaks,
                       labels = sprintf("%.2f", fav_breaks),
                       name = "Favourable\nallele fraction",
                       guide = guide_legend(
                         title = "Favourable\nallele fraction",
                         override.aes = list(
                           shape  = 21,
                           colour = "grey15",
                           stroke = 0.3,
                           # legend key sizes; floor (2.0) MUST match the plot
                           # floor below so the 0.00 key matches 0.00 dots.
                           size   = scales::rescale(fav_breaks, to = c(2.0, 7))
                         ))) +
  # range floor raised from 0.8 -> 2.0 so 0.00 dots are clearly visible
  # (empty cells still read as small dots, not near-invisible specks).
  scale_radius(range = c(2.0, 7), limits = c(0, 1), guide = "none") +
  facet_grid(GroupF ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_discrete(labels = y_label_fun) +
  labs(x = "Gene", y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    text             = element_text(colour = "grey20"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12,
                                    face = "bold.italic", colour = "grey20"),
    axis.title.x     = element_text(size = 13, face = "bold", margin = margin(t = 6)),
    axis.line        = element_line(linewidth = 0.4, colour = "grey30"),
    axis.ticks       = element_line(linewidth = 0.4, colour = "grey30"),
    panel.grid       = element_blank(),
    strip.background = element_blank(),
    strip.placement  = "outside",
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10),
    legend.key       = element_blank(),
    legend.key.height = unit(16, "pt"),
    panel.spacing    = unit(10, "pt"))

if (have_ggtext) {
  library(ggtext)
  strip_lab <- c(
    "Landraces"          = sprintf("<b style='color:%s'>Landraces</b>", land_col),
    setNames(sprintf("<b style='color:%s'>%s</b>", var_col, LAB_VAR), LAB_VAR))
  pA_main <- pA_main +
    theme(axis.text.y = ggtext::element_markdown(size = 12),
          strip.text.y.left = ggtext::element_markdown(angle = 90, size = 13)) +
    facet_grid(GroupF ~ ., scales = "free_y", space = "free_y", switch = "y",
               labeller = labeller(GroupF = strip_lab))
} else {
  pA_main <- pA_main +
    theme(axis.text.y = element_text(size = 12, face = "bold"),
          strip.text.y.left = element_text(angle = 90, size = 13, face = "bold"))
  message("Install 'ggtext' for bold, group-coloured labels.")
}

pA_strip <- ggplot(catstrip_df, aes(x = GeneLab, y = y, fill = Cat)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  scale_fill_manual(values = cat_palette_nat, name = "Trait category") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  guides(fill = guide_legend(byrow = TRUE, override.aes = list(colour = NA))) +
  theme_void() +
  theme(legend.title = element_text(size = 11),
        legend.text  = element_text(size = 10),
        legend.key.size = unit(12, "pt"), plot.margin = margin(0, 0, 2, 0))

fav_box_df <- data.frame(
  VarF   = factor(samples, levels = rev(samples)),
  GroupF = factor(ifelse(samples %in% landraces,
                         "Landraces", LAB_VAR),
                  levels = c("Landraces", LAB_VAR)),
  lab    = paste0(fav_count[samples], "/", n_markers),
  Count  = as.integer(fav_count[samples]), x = 1)

# Right-hand panel: horizontal BARS, faceted IDENTICALLY to pA_main
# (same GroupF facet, same free_y + space = "free_y") so every bar sits
# exactly parallel to its dot row. No in-panel title (a title would add
# height at the top and shift all rows down out of alignment); the column
# header is supplied by the patchwork design instead.
bar_fill <- ifelse(fav_box_df$VarF %in% landraces, land_col, var_col)

pA_box <- ggplot(fav_box_df, aes(x = Count, y = VarF)) +
  geom_col(aes(fill = GroupF), width = 0.62, show.legend = FALSE) +
  geom_text(aes(label = Count), hjust = -0.25,
            fontface = "bold", size = 3.6, colour = "grey20") +
  scale_fill_manual(values = c("Landraces" = land_col,
                               setNames(var_col, LAB_VAR))) +
  facet_grid(GroupF ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, n_markers * 1.15), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(strip.text = element_blank(),
        panel.spacing = unit(10, "pt"),     # MUST equal pA_main panel.spacing
        plot.margin = margin(0, 4, 0, 4))

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  # Layout grid:
  #   row 1  = trait strip (A) over the dot panel only;  D = "Fav /35" header
  #   rows 2+ = dot panel (B) parallel with the bar panel (C)
  # B and C share identical facet geometry, so their rows line up 1:1.
  hdr_df <- data.frame(x = 0, y = 0,
                       lab = paste0("Favourable\nalleles (/", n_markers, ")"))
  pA_hdr <- ggplot(hdr_df, aes(x, y, label = lab)) +
    geom_text(fontface = "bold", size = 3.6, lineheight = 0.95,
              hjust = 0) +                      # left-align so text starts at cell edge
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +             # let it overflow the cell, don't crop
    theme_void() +
    theme(plot.margin = margin(0, 2, 0, 0))

  design <- "
    AAAAAAAAAAAADD
    BBBBBBBBBBBBCC
    BBBBBBBBBBBBCC
    BBBBBBBBBBBBCC
    BBBBBBBBBBBBCC
    BBBBBBBBBBBBCC
    BBBBBBBBBBBBCC
  "
  pA <- pA_strip + pA_main + pA_box + pA_hdr +
    plot_layout(design = design, guides = "collect") &
    theme(plot.margin = margin(4, 14, 4, 4))    # extra right margin for the header
  ggsave("gene_balloon.png", pA, width = FIG_W, height = FIG_H,
         units = "mm", dpi = FIG_DPI, limitsize = FALSE)
} else {
  ggsave("gene_balloon.png", pA_main, width = FIG_W, height = FIG_H,
         units = "mm", dpi = FIG_DPI, limitsize = FALSE)
  message("Install 'patchwork' for the trait strip and Fav/total column.")
}

cat("Done. Files written:
",
    "  qc_marker_filter_report.csv
",
    "  nucleotide_genotypes.csv
  canonical_codes.csv
",
    "  gene_freq_matrix.csv
  locus_counts.csv
",
    "  gene_balloon.png
")
