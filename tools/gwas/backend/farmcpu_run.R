# Runs FarmCPU with the rMVP R package on files written by farmcpu_package.py.
#
# Rscript farmcpu_run.R <folder> <n samples> <m markers> <max loops>
#
# <folder> holds geno.bin (float32, sample-major: for each sample all m markers), map.csv (SNP,CHROM,POS) and phenotype.csv
# (Taxa,trait). Writes farmcpu.csv (effect,se,p per marker, in map order). The package does all the statistics.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("expected 4 arguments")
dir <- args[1]; n <- as.integer(args[2]); m <- as.integer(args[3]); loops <- as.integer(args[4])
fail <- function(msg) { writeLines(msg, file.path(dir, "error.txt")); quit(status = 2) }

if (!requireNamespace("rMVP", quietly = TRUE)) fail("rMVP package is not installed")
suppressPackageStartupMessages({ library(rMVP); library(bigmemory) })

con <- file(file.path(dir, "geno.bin"), "rb")
v <- readBin(con, what = "numeric", n = n * m, size = 4)
close(con)
if (length(v) != n * m) fail("genotype file has the wrong size")
geno <- as.big.matrix(matrix(v, nrow = m, ncol = n), type = "double")       # markers x samples
rm(v)
map <- read.csv(file.path(dir, "map.csv"), stringsAsFactors = FALSE)
phe <- read.csv(file.path(dir, "phenotype.csv"), stringsAsFactors = FALSE)

res <- try(MVP(phe = phe, geno = geno, map = map, method = "FarmCPU", maxLoop = loops, method.bin = "static",
               file.output = FALSE, verbose = FALSE, ncpus = 1), silent = TRUE)
if (inherits(res, "try-error")) fail(as.character(res))
fc <- res$farmcpu.results
if (is.null(fc) || nrow(fc) != m) fail("rMVP returned no FarmCPU results")
write.csv(data.frame(effect = fc[, 1], se = fc[, 2], p = fc[, ncol(fc)]), file.path(dir, "farmcpu.csv"), row.names = FALSE)
writeLines(as.character(packageVersion("rMVP")), file.path(dir, "version.txt"))
