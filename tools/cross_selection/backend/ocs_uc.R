# =============================================================================
# ocs_uc.R - usefulness of a set of crosses with SimpleMating::getUsefA()
#
#   Rscript ocs_uc.R <cfg.json>
#
# cfg: dir (holds geno.csv, effects.csv, map.csv, kmat.csv, pairs.csv), prop_sel, type (DH | RIL), generation
# writes <dir>/usefulness.csv  (Parent1, Parent2, Mean, Variance, sd, Usefulness)
# Exit code 3 = SimpleMating is not installed.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
cfg <- jsonlite::fromJSON(args[1])
if (!requireNamespace("SimpleMating", quietly = TRUE)) {
  cat("R package SimpleMating is not installed\n", file = stderr())
  quit(status = 3)
}
fail <- function(msg) { cat(msg, "\n", file = stderr()); quit(status = 1) }
f <- function(x) file.path(cfg$dir, x)

g <- read.csv(f("geno.csv"), check.names = FALSE, colClasses = c(Genotype = "character"))
M <- as.matrix(g[, -1, drop = FALSE]); mode(M) <- "numeric"; rownames(M) <- g$Genotype
eff <- read.csv(f("effects.csv"), stringsAsFactors = FALSE)
map <- read.csv(f("map.csv"), stringsAsFactors = FALSE)
kk <- read.csv(f("kmat.csv"), check.names = FALSE, colClasses = c(Genotype = "character"))
K <- as.matrix(kk[, -1, drop = FALSE]); mode(K) <- "numeric"; rownames(K) <- kk$Genotype
pairs <- read.csv(f("pairs.csv"), stringsAsFactors = FALSE, colClasses = "character")

if (!identical(colnames(M), map$mkr) || !identical(eff$Marker, map$mkr)) fail("The marker order differs between the files.")

res <- tryCatch(
  suppressWarnings(SimpleMating::getUsefA(
    MatePlan = pairs, Markers = M, addEff = eff$Effect, K = K, Map.In = map,
    propSel = as.numeric(cfg$prop_sel), Type = cfg$type, Generation = as.integer(cfg$generation),
    display_progress = FALSE)),
  error = function(e) fail(paste("getUsefA failed:", conditionMessage(e))))

u <- res[[1]]
out <- data.frame(Parent1 = u$Parent1, Parent2 = u$Parent2, Mean = u$Mean, Variance = u$Variance, sd = u$sd,
                  Usefulness = u$Usefulness)
write.csv(out, f("usefulness.csv"), row.names = FALSE)
