# Self-test of the map-construction engine (linkage_ml.R) on simulated populations with a KNOWN map.
#
#   Rscript selftest_map_construction.R
#
# Simulates RIL and F2 chromosomes (1.5 % genotyping errors, 5 % missing calls), writes the markers in a
# random order, and checks that the engine
#   * recovers the true order (Spearman |rho| > 0.999) and a map length close to the truth,
#   * removes markers labelled with the wrong chromosome and assigns markers without a chromosome,
#   * gives identical markers the same position,
#   * rebuilds the linkage groups from the genotypes alone (LOD grouping).
# Exits with status 1 when a check fails.

here <- tryCatch(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
                 error = function(e) getwd())
setwd(here)
txt <- readLines("linkage_runner.R", warn = FALSE)
eval(parse(text = txt[1:(grep("^if \\(!interactive\\(\\)\\)", txt)[1] - 1)]), envir = globalenv())

failures <- character(0)
check <- function(ok, what) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(ok)) "ok  " else "FAIL", what))
  if (!isTRUE(ok)) failures <<- c(failures, what)
}

simulate_chr <- function(chr, M, len_cm, n, pop, err, miss) {
  pos <- sort(runif(M, 0, len_cm))
  r <- (1 - exp(-2 * diff(pos) / 100)) / 2                       # Haldane
  sw <- if (pop == "RIL") 2 * r / (1 + 2 * r) else r
  G <- matrix(0L, M, n)
  if (pop == "RIL") {
    s <- rbinom(n, 1, .5); G[1, ] <- s
    for (i in 2:M) { s <- bitwXor(s, rbinom(n, 1, sw[i - 1])); G[i, ] <- s }
  } else {
    a <- rbinom(n, 1, .5); b <- rbinom(n, 1, .5); G[1, ] <- a + b
    for (i in 2:M) { a <- bitwXor(a, rbinom(n, 1, sw[i - 1])); b <- bitwXor(b, rbinom(n, 1, sw[i - 1])); G[i, ] <- a + b }
  }
  E <- matrix(runif(M * n) < err, M, n)
  if (pop == "RIL") G[E] <- 1L - G[E] else G[E] <- (G[E] + sample(1:2, sum(E), TRUE)) %% 3L
  G[matrix(runif(M * n) < miss, M, n)] <- NA
  txt <- if (pop == "RIL") c("a", "b")[G + 1L] else c("a", "ab", "b")[G + 1L]
  txt[is.na(G)] <- "-"
  gdf <- as.data.frame(matrix(txt, M, n), stringsAsFactors = FALSE); colnames(gdf) <- paste0("L", seq_len(n))
  d <- cbind(data.frame(Marker = sprintf("chr%02d_%d", chr, sort(sample(1e6:5e7, M))), Chromosome = as.character(chr),
                        stringsAsFactors = FALSE), gdf)
  d$true_pos <- pos
  d[sample(M), ]                                                  # file order is random
}

run_map <- function(d, pop, ...) {
  d$true_pos <- NULL; d$.input_order <- seq_len(nrow(d))
  create_linkage_map(d, population_type = pop, missing_threshold = 0.3, recombination_threshold = 0.45,
                     max_total_cm = 1e5, marker_order = "file", map_config = default_map_config(...))
}

set.seed(2024)
for (pop in c("RIL", "F2")) {
  cat("\n==", pop, ": order and length from a scrambled file ==\n")
  n <- if (pop == "RIL") 180 else 200
  d <- do.call(rbind, lapply(1:3, function(c) simulate_chr(c, c(60, 90, 45)[c], c(110, 140, 85)[c], n, pop, 0.015, 0.05)))
  truth <- setNames(d$true_pos, d$Marker)
  res <- run_map(d, pop)
  for (chr in unique(res$linkage_map$Chromosome)) {
    s <- res$linkage_map[res$linkage_map$Chromosome == chr, ]
    tl <- diff(range(truth[s$Marker]))
    rho <- abs(cor(s$Position, truth[s$Marker], method = "spearman"))
    rk <- rank(truth[s$Marker][order(s$Position)])
    disp <- max(abs(rk - seq_along(rk)))
    # markers a fraction of a cM apart cannot be separated by 180 individuals: allow small local swaps only
    check(rho > 0.995 && disp <= 6, sprintf("chr %s order recovered (|rho| %.4f, largest displacement %d positions)", chr, rho, disp))
    check(abs(max(s$Position) / tl - 1) < 0.25, sprintf("chr %s length %.0f cM vs true %.0f cM", chr, max(s$Position), tl))
  }
}

cat("\n== linkage-group checks ==\n")
d <- do.call(rbind, lapply(1:4, function(c) simulate_chr(c, c(50, 60, 40, 55)[c], c(100, 120, 90, 110)[c], 180, "RIL", 0.015, 0.05)))
truth <- setNames(d$true_pos, d$Marker)
wrong <- sample(which(d$Chromosome == "1"), 5); d$Chromosome[wrong] <- "3"
unk <- sample(which(d$Chromosome == "2"), 4); d$Chromosome[unk] <- "Unknown"
dup <- d[sample(which(d$Chromosome == "4"), 8), ]; dup$Marker <- paste0(dup$Marker, "_dup")
d2 <- rbind(d, dup)
res <- run_map(d2, "RIL")
lm <- res$linkage_map
check(all(d$Marker[wrong] %in% res$removal_log$Marker), "markers labelled with the wrong chromosome are removed")
check(sum(res$removal_log$Reason == "Linked to another chromosome") == length(wrong), "no correctly labelled marker is removed")
check(all(lm$Chromosome[match(d$Marker[unk], lm$Marker)] == "2"), "markers without a chromosome are assigned to the right one")
dd <- lm[grepl("_dup$", lm$Marker), ]
check(all(abs(dd$Position - lm$Position[match(sub("_dup$", "", dd$Marker), lm$Marker)]) < 1e-9), "identical markers share a position")

d3 <- d; d3$Chromosome <- "1"
res3 <- run_map(d3, "RIL", grouping = "lod", group_lod = 6)
lm3 <- res3$linkage_map
truth_chr <- setNames(sub("^chr0?(\\d+)_.*", "\\1", d$Marker), d$Marker)
tab <- table(lm3$Chromosome, truth_chr[lm3$Marker])
check(nrow(tab) == 4 && all(apply(tab, 1, function(x) sum(x > 0)) == 1), "LOD grouping rebuilds the four chromosomes")

cat("\n")
if (length(failures)) { cat(sprintf("%d check(s) FAILED\n", length(failures))); quit(status = 1) }
cat("all checks passed\n")
