# Run reference GWAS programs on the files written by make_inputs.py, so their p-values can be compared with this program's.
#
#   Rscript tools/gwas/benchmark/run_references.R benchmark_run
#
# Each block runs only if its package is installed, and otherwise prints how to install it. Nothing is faked.
# These calls follow the packages' documented interfaces but have NOT been run in the environment where this script was written
# (none of the reference programs were installed there). Check argument names against the version you install.
# Use the same number of principal components as pcs_used in truth.json.

args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[1] else "benchmark_run"
f <- function(x) file.path(dir, x)
truth <- jsonlite::fromJSON(f("truth.json"))
npc <- truth$pcs_used

geno <- read.csv(f("numeric_genotypes.csv"), check.names = FALSE)
map <- read.csv(f("map.csv"))
phe <- read.csv(f("phenotype.csv"))
pcs <- read.csv(f("pcs.csv"))
kin <- read.csv(f("kinship.csv"), check.names = FALSE)
stopifnot(identical(as.character(geno$Taxa), as.character(phe$Taxa)))
have <- function(pkg) requireNamespace(pkg, quietly = TRUE)
need <- function(pkg, how) message("[not run] ", pkg, " is not installed. Install with: ", how)

# ---- GAPIT (GLM, MLM, FarmCPU, Blink) ----
if (have("GAPIT") || file.exists("gapit_functions.txt")) {
  if (!have("GAPIT")) source("gapit_functions.txt") else library(GAPIT)
  GD <- data.frame(taxa = geno$Taxa, geno[, -1], check.names = FALSE)
  GM <- data.frame(SNP = map$SNP, Chromosome = map$Chromosome, Position = map$Position)
  CV <- if (npc > 0) pcs[, 1:(npc + 1)] else NULL
  for (m in c("GLM", "MLM", "FarmCPU", "Blink")) {
    out <- GAPIT(Y = phe, GD = GD, GM = GM, KI = kin, CV = CV, PCA.total = 0, model = m, file.output = FALSE)
    res <- out$GWAS
    write.csv(data.frame(SNP = res$SNP, P = res$P.value), f(paste0("gapit_", m, ".csv")), row.names = FALSE)
  }
} else need("GAPIT", "source('http://zzlab.net/GAPIT/gapit_functions.txt')  # or install from https://github.com/jiabowang/GAPIT")

# ---- rMVP (GLM, MLM, FarmCPU) ----
if (have("rMVP")) {
  library(rMVP)
  # rMVP reads markers x individuals matrices; write them in its own format, then call MVP() with the same K and covariates.
  message("rMVP is installed. Convert numeric_genotypes.csv with MVP.Data.Numeric2MVP(), then call MVP(phe=, geno=, map=, K=, ",
          "CV.GLM/CV.MLM/CV.FarmCPU=<first ", npc, " PCs>, nPC.*=0, method=c('GLM','MLM','FarmCPU')) and save the p-value columns as ",
          "rmvp_<model>.csv with columns SNP,P. The result layout differs between versions, so this step is left to the user.")
} else need("rMVP", "install.packages('rMVP')")

# ---- EMMAX and GEMMA are command-line programs: see README.md ----
