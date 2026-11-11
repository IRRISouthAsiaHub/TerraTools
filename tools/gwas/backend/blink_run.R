# Runs BLINK with the GAPIT R package on files written by blink_package.py.
#
# Rscript blink_run.R <folder> <n samples> <m markers>
#
# <folder> holds geno.bin (float32, sample-major), map.csv (SNP,CHROM,POS) and phenotype.csv (Taxa,trait1). Writes blink.csv
# (SNP,effect,p per marker). The package does all the statistics.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("expected 3 arguments")
dir <- args[1]; n <- as.integer(args[2]); m <- as.integer(args[3])
fail <- function(msg) { writeLines(msg, file.path(dir, "error.txt")); quit(status = 2) }

if (!requireNamespace("GAPIT", quietly = TRUE)) fail("GAPIT package is not installed")
suppressPackageStartupMessages(library(GAPIT))
setwd(dir)                                       # GAPIT may write scratch files to the working directory

v <- readBin(file.path(dir, "geno.bin"), what = "numeric", n = n * m, size = 4)
if (length(v) != n * m) fail("genotype file has the wrong size")
G <- matrix(v, nrow = m, ncol = n)               # markers x samples
rm(v)
map <- read.csv(file.path(dir, "map.csv"), stringsAsFactors = FALSE)
phe <- read.csv(file.path(dir, "phenotype.csv"), stringsAsFactors = FALSE)
GD <- data.frame(taxa = phe$Taxa, t(G), check.names = FALSE)
colnames(GD)[-1] <- map$SNP
GM <- data.frame(SNP = map$SNP, Chromosome = map$CHROM, Position = map$POS)

res <- try(GAPIT(Y = phe, GD = GD, GM = GM, model = "Blink", file.output = FALSE, SNP.MAF = 0), silent = TRUE)
if (inherits(res, "try-error")) fail(as.character(res))
g <- res$GWAS
if (is.null(g) || !all(c("SNP", "P.value") %in% colnames(g))) fail("GAPIT returned no BLINK results")
write.csv(data.frame(SNP = g$SNP, effect = g$effect, p = g$P.value), file.path(dir, "blink.csv"), row.names = FALSE)
writeLines(as.character(packageVersion("GAPIT")), file.path(dir, "version.txt"))
