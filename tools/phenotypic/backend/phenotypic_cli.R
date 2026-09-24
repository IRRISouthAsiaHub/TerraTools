# phenotypic_cli.R — Command-line entrypoint for the phenotypic engine.
#
# Called by phenotypic_runner.py via subprocess. Writes summary.json into
# --outdir on completion (success OR failure) so Python never has to parse
# stdout — it just reads that one file after the process exits.

suppressPackageStartupMessages(library(jsonlite))

this_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(this_dir, "phenotypic_engine.R"))

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

selected_envs <- if (!is.null(args[["environments"]]) && nzchar(args[["environments"]]))
  trimws(strsplit(args[["environments"]], "\\|\\|")[[1]]) else character(0)

result <- tryCatch({
  run_phenotypic(
    input_file            = args[["input"]],
    trait                 = args[["trait"]],
    include_interaction   = tolower(args[["interaction"]]) %in% c("true", "1", "yes"),
    analysis_mode         = args[["mode"]],
    selected_environments = selected_envs,
    n_clusters            = as.integer(args[["n_clusters"]]),
    clustering_method     = args[["clustering_method"]],
    color_scheme          = args[["color_scheme"]],
    outdir                = outdir
  )
}, error = function(e) {
  list(status = "error", message = conditionMessage(e))
})

# run_phenotypic() already writes summary.json on success; only write it
# here if the top-level tryCatch caught something run_phenotypic() itself
# didn't handle (e.g. a failure before its own internal writer ran).
if (!file.exists(summary_path)) {
  jsonlite::write_json(result, summary_path, auto_unbox = TRUE, na = "null", null = "null")
}

if (identical(result$status, "error")) {
  quit(status = 1)
} else {
  quit(status = 0)
}
