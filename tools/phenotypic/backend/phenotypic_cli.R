# phenotypic_cli.R — command-line entrypoint of the phenotypic engine.
#
# Called by phenotypic_runner.py via subprocess. The result is written to <outdir>/result_<action>.json (success OR
# failure), so Python never has to parse stdout.
#
#   Rscript phenotypic_cli.R --action inspect --params params.json --outdir DIR
#       reads and validates the data file: columns, environments, usable traits, design check
#   Rscript phenotypic_cli.R --action run --params params.json --outdir DIR
#       runs ONE analysis module (params$module) with its settings; earlier results in DIR are kept
#   Rscript phenotypic_cli.R --action export --outdir DIR
#       writes the Excel workbooks for everything run so far in DIR
#
# params.json: input_file, design, analysis_mode, selected_environments, color_scheme, zero_as_na, module, settings

suppressPackageStartupMessages(library(jsonlite))

.tt_backend_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(.tt_backend_dir, "phenotypic_engine.R"), encoding = "UTF-8")

parse_args <- function(args) {
  out <- list(); i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (startsWith(key, "--")) {
      out[[substring(key, 3)]] <- if (i + 1 <= length(args)) args[[i + 1]] else ""
      i <- i + 2
    } else i <- i + 1
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
action <- if (is.null(args[["action"]])) "run" else args[["action"]]
outdir <- args[["outdir"]]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
result_path <- file.path(outdir, paste0("result_", action, ".json"))

chr <- function(x) if (is.null(x) || length(x) == 0) character(0) else as.character(unlist(x))
as_rows <- function(df) if (is.null(df)) list() else df

# ── inspect: what the data file holds, as the analyses will see it ───────────────────────
inspect_data <- function(p) {
  data <- read_data_file(p$input_file)
  # "auto": the design the columns suggest (a Block column means augmented / alpha-lattice, none means RCBD)
  design <- if (identical(p$design, "auto")) (if ("Block" %in% names(data)) "augmented" else "rcbd") else normalize_design(p$design)
  all_envs <- if ("Environment" %in% names(data)) unique(as.character(data$Environment)) else character(0)
  v <- validate_data(data, design)
  dd <- v$data
  sel <- chr(p$selected_environments)
  if (identical(p$analysis_mode, "specific") && length(sel) && "Environment" %in% names(dd))
    dd <- dd[as.character(dd$Environment) %in% sel, , drop = FALSE]
  all_tr <- find_traits(dd)
  z <- if (!isFALSE(p$zero_as_na)) zeros_to_missing(dd, all_tr) else list(data = dd, notes = character(0))
  dd <- z$data
  traits <- all_tr[vapply(all_tr, function(x) any(!is.na(dd[[x]])), logical(1))]
  empty <- setdiff(all_tr, traits)
  w <- v$validation$warnings
  w <- w[!grepl("converted to factor", w)]                                  # routine, not worth showing

  stats <- if (length(traits) && length(v$validation$errors) == 0) do.call(rbind, lapply(traits, function(tr) {
    x <- dd[[tr]]; xv <- x[!is.na(x)]
    data.frame(Trait = tr, Plots = length(x), `Missing (incl. zeros)` = sum(is.na(x)), Min = round(min(xv), 3),
               Mean = round(mean(xv), 3), Max = round(max(xv), 3), check.names = FALSE)
  }))
  dcheck <- NULL; dnotes <- character(0)
  if (length(v$validation$errors) == 0) {
    rep <- tryCatch(design_report(dd, design), error = function(e) NULL)
    if (!is.null(rep)) { dcheck <- tidy_design_check(rep$table); dnotes <- rep$notes }
  }
  list(status = "success", design = list(code = design, label = DESIGN_LABELS[[design]]),
       default_design = if ("Block" %in% names(data)) "augmented" else "rcbd",
       columns = names(data), environments = all_envs,
       n_plots = nrow(dd), n_genotypes = if ("Genotype" %in% names(dd)) length(unique(dd$Genotype)) else 0,
       n_environments = length(unique(as.character(dd$Environment))), n_variables = ncol(data),
       traits = traits, empty_traits = empty, usable_traits = usable_traits(data),
       errors = v$validation$errors, warnings = w, notes = z$notes, design_notes = dnotes,
       trait_stats = as_rows(stats), design_check = as_rows(dcheck))
}

# ── run: one analysis module ─────────────────────────────────────────────────────────────
run_module <- function(p) {
  m <- p$module
  if (is.null(m) || !(m %in% ALL_MODULES)) stop("Unknown analysis: ", m, call. = FALSE)
  s <- p$settings
  if (is.null(s)) s <- list()
  for (k in c("traits")) if (!is.null(s[[k]])) s[[k]] <- chr(s[[k]])
  res <- run_phenotypic(
    input_file = p$input_file, outdir = outdir, design = p$design,
    analysis_mode = if (is.null(p$analysis_mode)) "all" else p$analysis_mode,
    selected_environments = chr(p$selected_environments),
    color_scheme = if (is.null(p$color_scheme)) "default" else p$color_scheme,
    zero_as_na = !isFALSE(p$zero_as_na),
    modules = m, settings = setNames(list(s), m), accumulate = TRUE, workbook = FALSE)
  out <- list(status = "success", module = m, outdir = outdir, design = res$design, traits_all = res$traits_all,
              tables = res$tables, plots = res$plots, msgs = res$msgs)
  if (m == "anova") {                                       # the wide layout shown on screen (the CSV is the tidy version)
    wide <- tryCatch(anova_wide_table(outdir), error = function(e) NULL)
    if (!is.null(wide)) {
      wp <- file.path(outdir, "anova_wide.json")
      jsonlite::write_json(wide, wp, na = "null", digits = 10, dataframe = "rows")
      out$anova_wide <- wp
    }
  }
  out
}

result <- tryCatch({
  if (action == "export") {
    f <- export_workbook(outdir)
    list(status = "success", excel_all = f$all, excel_anova = f$anova)
  } else {
    p <- jsonlite::fromJSON(args[["params"]], simplifyVector = TRUE)
    if (action == "inspect") inspect_data(p) else run_module(p)
  }
}, error = function(e) list(status = "error", message = conditionMessage(e)))

jsonlite::write_json(result, result_path, auto_unbox = TRUE, na = "null", null = "null", digits = 10, dataframe = "rows")
quit(status = if (identical(result$status, "error")) 1 else 0)
