# mendelian_cli.R — Command-line entrypoint for the Mendelian engine.
#
# Called by mendelian_runner.py via subprocess. Writes summary.json into
# --outdir on completion (success OR failure) so the Python side never has
# to parse stdout — it just reads that one file after the process exits.
#
# Usage:
#   Rscript mendelian_cli.R --input path.csv --pop RIL --pthr 0.05 \
#           --missing "NA,.,-9,-1" --mtc fdr --drop_unplaced true --outdir /tmp/xyz

suppressPackageStartupMessages(library(jsonlite))

this_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(this_dir, "mendelian_engine.R"))

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (startsWith(key, "--")) {
      key <- substring(key, 3)
      val <- if (i + 1 <= length(args)) args[[i + 1]] else ""
      out[[key]] <- val
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

outdir <- args[["outdir"]]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
summary_path <- file.path(outdir, "summary.json")

result <- tryCatch({
  run_mendelian(
    input_file     = args[["input"]],
    pop            = args[["pop"]],
    pthr           = as.numeric(args[["pthr"]]),
    missing_codes  = trimws(strsplit(args[["missing"]], ",")[[1]]),
    drop_unplaced  = tolower(args[["drop_unplaced"]]) %in% c("true", "1", "yes"),
    outdir         = outdir,
    mtc_method     = if (is.null(args[["mtc"]])) "fdr" else args[["mtc"]]
  )
}, error = function(e) {
  list(status = "error", message = conditionMessage(e))
})

jsonlite::write_json(result, summary_path, auto_unbox = TRUE, na = "null")

if (identical(result$status, "error")) {
  quit(status = 1)
} else {
  quit(status = 0)
}
