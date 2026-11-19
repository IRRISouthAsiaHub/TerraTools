# GWAS figures drawn in R:
#   Manhattan (one per model)  base graphics: one pastel colour per chromosome, marker-density ideogram strip, solid / dashed
#                              significance lines, the name of each genome-wide / QTN marker above the plot on a vertical line.
#                              The y axis ends at the highest point or threshold of that model.
#   Circular Manhattan         CMplot: black chromosome band, one ring per model with its own y axis.
#   QQ plot                    CMplot: observed against expected -log10(p); the page takes the proportions of the data.
#
# usage: Rscript gwas_plots.R <config.json>
#   config: { input_csv, out_dir, models, thresholds: {model: {bonferroni_p, suggestive_p, permutation_p}}, hits: {model: [snp]},
#             lambdas: {model: x}, show_bonferroni, show_suggestive, max_labels, dpi, model_colors: [hex] }
# input_csv: SNP, Chromosome, Position and one column of p-values per model.

suppressPackageStartupMessages(library(jsonlite))

PASTEL <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
            "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

axis_ticks <- function(vmax, max_intervals = 4) {
  vmax <- max(vmax, 1e-9)
  mag <- 10^floor(log10(vmax))
  for (m in c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50)) {
    step <- m * mag
    if (ceiling(vmax / step - 1e-9) <= max_intervals) break
  }
  top <- ceiling(vmax / step - 1e-9) * step
  list(ticks = seq(0, top, by = step), top = top)
}

fmt <- function(v) format(v, trim = TRUE, scientific = FALSE, digits = 6)

# horizontal positions for labels that must stay `sep` apart inside [lo, hi]
spread <- function(xs, lo, hi, sep) {
  out <- xs
  if (length(out) > 1) for (k in 2:length(out)) out[k] <- max(out[k], out[k - 1] + sep)
  if (length(out) && out[length(out)] > hi) {
    out[length(out)] <- hi
    if (length(out) > 1) for (k in (length(out) - 1):1) out[k] <- min(out[k], out[k + 1] - sep)
  }
  pmax(out, lo)
}

layout_genome <- function(chr, pos) {
  chrs <- sort(unique(chr))
  gap <- max(pos) * 0.018 + 1
  x <- numeric(length(chr)); start <- numeric(length(chrs)); end <- numeric(length(chrs))
  off <- 0
  for (k in seq_along(chrs)) {
    m <- chr == chrs[k]
    x[m] <- pos[m] - min(pos[m]) + off
    start[k] <- min(x[m]); end[k] <- max(x[m])
    off <- off + (max(pos[m]) - min(pos[m])) + gap
  }
  list(x = x, chrs = chrs, start = start, end = end, total = max(end), gap = gap)
}

manhattan <- function(file, d, model, thr, hits, cfg) {
  p <- d[[model]]
  ok <- is.finite(p) & p > 0
  lp <- rep(NA_real_, length(p)); lp[ok] <- -log10(p[ok])
  L <- layout_genome(d$Chromosome, d$Position)
  chrs <- L$chrs; total <- L$total; gap <- L$gap
  lines_def <- list()
  if (isTRUE(cfg$show_suggestive) && !is.null(thr$suggestive_p))
    lines_def[[length(lines_def) + 1]] <- list(v = -log10(thr$suggestive_p), lty = 2, col = "#444444", name = "suggestive")
  if (isTRUE(cfg$show_bonferroni) && !is.null(thr$bonferroni_p))
    lines_def[[length(lines_def) + 1]] <- list(v = -log10(thr$bonferroni_p), lty = 1, col = "#111111", name = "Bonferroni")
  if (!is.null(thr$permutation_p) && is.finite(thr$permutation_p) && thr$permutation_p > 0)
    lines_def[[length(lines_def) + 1]] <- list(v = -log10(thr$permutation_p), lty = 1, col = "#3A8F3A", name = "simulated null")
  top_data <- if (any(is.finite(lp))) max(lp, na.rm = TRUE) else 1
  tk <- axis_ticks(max(c(top_data, vapply(lines_def, function(z) z$v, 0))) * 1.02)
  top <- tk$top

  # marker names above the plot: the strongest marker of each cluster of hits, at most max_labels
  hit_i <- which(d$SNP %in% hits & is.finite(lp))
  chosen <- integer(0)
  for (i in hit_i[order(-lp[hit_i])]) {
    if (all(abs(L$x[i] - L$x[chosen]) > total * 0.03)) chosen <- c(chosen, i)
    if (length(chosen) >= cfg$max_labels) break
  }
  chosen <- chosen[order(L$x[chosen])]
  label_in <- if (length(chosen)) max(nchar(as.character(d$SNP[chosen]))) * 0.105 + 0.35 else 0.35   # inches for the names

  W <- 14; H <- 4.6 + label_in
  png(file, width = W, height = H, units = "in", res = cfg$dpi, type = "cairo", bg = "white")
  par(mai = c(0.25, 1.15, label_in, 1.1), xaxs = "i", yaxs = "i", family = "sans")
  strip_h <- 0.085 * top; strip_y0 <- -0.045 * top
  plot(NA, xlim = c(-gap * 0.6, total + gap * 0.3), ylim = c(strip_y0 - strip_h - 0.12 * top, top), axes = FALSE, xlab = "", ylab = "")
  n <- sum(is.finite(lp))
  cx <- max(0.4, min(1.5, 2.3 - 0.36 * log10(max(n, 10))))
  col_of <- PASTEL[(match(d$Chromosome, chrs) - 1) %% length(PASTEL) + 1]
  for (k in seq_along(chrs)) {
    m <- which(d$Chromosome == chrs[k] & is.finite(lp))
    points(L$x[m], lp[m], pch = 16, cex = cx, col = PASTEL[(k - 1) %% length(PASTEL) + 1])
  }
  for (z in lines_def) {
    segments(-gap * 0.6, z$v, total + gap * 0.3, z$v, lty = z$lty, col = z$col, lwd = 1)
    text(total + gap * 0.45, z$v, z$name, adj = c(0, 0.5), cex = 0.7, col = z$col, xpd = NA)
  }
  hi <- hit_i
  if (length(hi)) points(L$x[hi], lp[hi], pch = 21, cex = cx * 1.25, col = "#222222", bg = col_of[hi], lwd = 0.6)
  if (!length(hit_i) && !is.null(hits))
    text(0, top, "No genome-wide significant associations detected", adj = c(0, 1.2), col = "#8A1F1F", font = 2, cex = 0.9)

  # y axis, left of the first chromosome only
  x0 <- -gap * 0.6
  segments(x0, 0, x0, top, lwd = 1.2)
  segments(x0 - gap * 0.18, tk$ticks, x0, tk$ticks, lwd = 1.2)
  text(x0 - gap * 0.3, tk$ticks, vapply(tk$ticks, fmt, ""), adj = c(1, 0.5), font = 2, cex = 1, xpd = NA)
  mtext(bquote(.(model) ~~ -log[10](italic(p))), side = 2, line = 3.1, cex = 1.35)

  # names above the plot, on vertical lines
  if (length(chosen)) {
    lx <- spread(L$x[chosen], 0, total, total * 0.021)
    yb <- top * 1.045
    for (j in seq_along(chosen)) {
      i <- chosen[j]
      segments(L$x[i], 0, L$x[i], top * 1.0, lwd = 0.7, col = "#333333", xpd = NA)
      segments(L$x[i], top, L$x[i], top * 1.02, lwd = 0.7, col = "#333333", xpd = NA)
      segments(L$x[i], top * 1.02, lx[j], yb, lwd = 0.7, col = "#333333", xpd = NA)
      text(lx[j], yb, as.character(d$SNP[i]), srt = 90, adj = c(0, 0.5), cex = 1.05, xpd = NA)
    }
  }

  # marker-density ideogram and chromosome names
  for (k in seq_along(chrs)) {
    m <- d$Chromosome == chrs[k]
    xs <- L$x[m]
    segments(xs, strip_y0 - strip_h * 0.06, xs, strip_y0 - strip_h * 0.94, lwd = 0.35, col = "#111111")
    rect(L$start[k], strip_y0, max(L$end[k], L$start[k] + total * 0.002), strip_y0 - strip_h, border = "#111111", lwd = 0.9)
    text((L$start[k] + L$end[k]) / 2, strip_y0 - strip_h - 0.06 * top, paste0("chr", chrs[k]),
         cex = if (length(chrs) <= 16) 0.95 else 0.7, srt = if (length(chrs) <= 16) 0 else 90, adj = if (length(chrs) <= 16) c(0.5, 1) else c(1, 0.5), xpd = NA)
  }
  dev.off()
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  cfg <- fromJSON(args[1], simplifyVector = TRUE, simplifyDataFrame = FALSE)
  out <- cfg$out_dir
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  d <- read.csv(cfg$input_csv, check.names = FALSE, stringsAsFactors = FALSE)
  models <- as.character(cfg$models)
  if (is.null(cfg$dpi)) cfg$dpi <- 300
  if (is.null(cfg$max_labels)) cfg$max_labels <- 12
  made <- list()

  for (m in models) {
    f <- file.path(out, if (length(models) == 1) "01_Manhattan.png" else paste0("01_Manhattan_", m, ".png"))
    res <- tryCatch({ manhattan(f, d, m, cfg$thresholds[[m]], as.character(unlist(cfg$hits[[m]])), cfg); TRUE },
                    error = function(e) { message("Manhattan ", m, " failed: ", conditionMessage(e)); try(dev.off(), silent = TRUE); FALSE })
    if (res) made[[basename(f)]] <- f
  }

  suppressPackageStartupMessages(library(CMplot))
  work <- tempfile("cmplot_"); dir.create(work); old <- setwd(work); on.exit(setwd(old), add = TRUE)
  chrs <- sort(unique(d$Chromosome))
  chr_lab <- paste0("Chr", chrs)

  # QQ: observed against expected, each axis ends at its own largest value and the page takes the proportions of the data
  tryCatch({
    n_ok <- sapply(models, function(m) sum(is.finite(d[[m]]) & d[[m]] > 0))
    xmax <- log10(max(n_ok)); ymax <- max(sapply(models, function(m) { p <- d[[m]]; -log10(min(p[is.finite(p) & p > 0])) }))
    ratio <- ymax / xmax
    base <- 5.2
    wp <- if (ratio >= 1) max(base / ratio, 2.6) else base
    hp <- if (ratio >= 1) base else max(base * ratio, 2.6)
    q <- d[, c("SNP", "Chromosome", "Position", models), drop = FALSE]
    lam <- unlist(cfg$lambdas)
    names(q)[4:ncol(q)] <- ifelse(!is.na(lam[models]), sprintf("%s  λ = %.3f", models, lam[models]), models)
    cols <- if (!is.null(cfg$model_colors)) as.character(unlist(cfg$model_colors))[seq_along(models)] else NULL
    CMplot(q, plot.type = "q", multraits = length(models) > 1, conf.int = FALSE, col = cols, file = "png", file.name = "qq",
           dpi = cfg$dpi, width = wp + 1.9, height = hp + 1.9, main = "", verbose = FALSE)
    f <- list.files(work, pattern = "QQplot.*png$", full.names = TRUE)
    if (length(f)) { file.copy(f[1], file.path(out, "02_QQ.png"), overwrite = TRUE); made[["02_QQ.png"]] <- file.path(out, "02_QQ.png") }
  }, error = function(e) message("QQ failed: ", conditionMessage(e)))

  # circular Manhattan: one ring per model, each with its own y axis
  tryCatch({
    q <- d[, c("SNP", "Chromosome", "Position", models), drop = FALSE]
    thr <- NULL
    bon <- unique(unlist(lapply(models, function(m) cfg$thresholds[[m]]$bonferroni_p)))
    if (length(bon) >= 1) thr <- if (length(bon) == 1) bon else list(bon)
    fc <- file.path(out, "03_Circos.png")
    png(fc, width = 10, height = 10, units = "in", res = cfg$dpi, type = "cairo", bg = "white")
    par(mar = c(0.5, 0.5, 0.5, 0.5))
    CMplot(q, plot.type = "c", r = 0.4, cir.chr.h = 1.3, cir.band = 1, chr.labels = chr_lab, threshold = if (length(bon) == 1) bon else NULL,
           amplify = length(bon) == 1, signal.col = "#E4572E", signal.cex = 1, signal.line = NULL, outward = FALSE,
           file.output = FALSE, verbose = FALSE)
    # which ring is which model: listed in the top-left corner, outer ring first
    u <- par("usr")
    par(xpd = NA)
    legend(u[1], u[4], legend = c("Rings, outer to inner", paste0(seq_along(models), "  ", models)), bty = "n", text.font = 2,
           cex = 1, text.col = "#222222", y.intersp = 1.1)
    dev.off()
    made[["03_Circos.png"]] <- fc
  }, error = function(e) message("Circular plot failed: ", conditionMessage(e)))

  cat(toJSON(list(ok = TRUE, files = unlist(made)), auto_unbox = TRUE), "\n")
}

main()
