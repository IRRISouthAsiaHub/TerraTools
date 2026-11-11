# =============================================================================
# MAIN PIPELINE: run_phenotypic()
# =============================================================================

# Reads (or takes) the data, validates it and applies the missing-value rules:
#   * 0 in a trait column is a missing value (not measured) when zero_as_na is TRUE; non-finite values are missing;
#   * traits left with no value at all are dropped.
# Returns list(data, traits, notes).
prepare_data <- function(input_file, design, zero_as_na = TRUE) {
  data <- if (is.data.frame(input_file)) as.data.frame(input_file) else read_data_file(input_file)
  names(data) <- canonical_names(names(data))
  v <- validate_data(data, design)
  if (length(v$validation$errors) > 0) stop(paste(v$validation$errors, collapse = "; "), call. = FALSE)
  data <- v$data
  traits <- find_traits(data)
  for (tr in traits) data[[tr]][!is.finite(data[[tr]])] <- NA_real_

  notes <- character(0)
  if (isTRUE(zero_as_na)) {
    z <- zeros_to_missing(data, traits)
    data <- z$data; notes <- z$notes
  }
  empty <- traits[!vapply(traits, function(tr) any(!is.na(data[[tr]])), logical(1))]
  traits <- setdiff(traits, empty)
  if (length(traits) == 0) stop("Every trait is empty (all values missing or 0).", call. = FALSE)
  if (length(empty)) notes <- c(notes, paste("Left out because they have no data:", paste(empty, collapse = ", ")))
  list(data = data, traits = traits, notes = notes)
}

# The rows of the chosen environments ("all" keeps everything)
select_environments <- function(data, mode, selected) {
  if (!identical(mode, "specific")) return(data)
  keep <- intersect(as.character(selected), unique(as.character(data$Environment)))
  if (length(keep) == 0) stop("Choose at least one environment.", call. = FALSE)
  data[as.character(data$Environment) %in% keep, , drop = FALSE]
}

# run_phenotypic(): loads the data once and runs only the requested modules (see 06_modules.R).
#   input_file: path of a .csv / .xlsx / .xls file, or a data frame the caller has already loaded
#   modules   : which analyses to run (a subset of ALL_MODULES)
#   settings  : one list per module, e.g. list(cluster = list(traits = c("PH", "DTF"), k = 4, method = "kmeans"),
#               anova = list(traits = ..., diag_trait = "PH"), blues = list(trait = "PH", interaction = TRUE)).
#               Every module also takes plots = TRUE/FALSE. A module without settings uses every trait.
#   zero_as_na: 0 in a trait column is a missing value (not measured); on by default
#   accumulate: results of earlier runs in the same `outdir` are kept (the workbook lists every analysis done so far)
#   workbook  : write the Excel workbooks at the end of the run. The app passes FALSE and calls export_workbook() only
#               when the user asks for the file, so a run does not pay for it.
# Returns list(tables, plots, msgs, excel_all, excel_anova, ...). Every message is tagged with the module that
# produced it (module, type = note / warning / error, text).
run_phenotypic <- function(input_file, outdir, design = "augmented",
                           analysis_mode = "all", selected_environments = character(0),
                           color_scheme = "default", zero_as_na = TRUE,
                           modules = ALL_MODULES, settings = list(), accumulate = TRUE, workbook = TRUE) {

  design <- normalize_design(design)
  modules <- intersect(ALL_MODULES, modules)
  if (!length(modules)) stop("Select at least one analysis to run.", call. = FALSE)

  prep <- prepare_data(input_file, design, zero_as_na)
  data <- prep$data                                   # all environments (stability always uses these)
  traits_all <- prep$traits
  working_data <- select_environments(data, analysis_mode, selected_environments)   # what the other modules use

  plots_dir <- file.path(outdir, "plots")
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

  out <- list(status = "success", design = list(code = design, label = DESIGN_LABELS[[design]]),
              traits_all = traits_all, outdir = outdir, tables = list(), plots = list(), msgs = list())

  # ── Messages ───────────────────────────────────────────────────────────────
  cur <- "data"                                       # the module that is running; tags every message
  push <- function(type, text) for (s in text) out$msgs[[length(out$msgs) + 1]] <<- list(module = cur, type = type, text = s)
  add_note <- function(x) push("note", x)
  add_warning <- function(x) push("warning", x)
  add_error <- function(name, text) push("error", text)

  # Runs one analysis; an error is recorded and the rest of the module continues.
  # NB: `expr` must not use return() — pass a function call (see 03_analysis.R) instead.
  safe <- function(name, expr) tryCatch(expr, error = function(e) { add_error(name, conditionMessage(e)); NULL })
  # Runs a whole module; a problem with its settings (too few traits ...) stops only this module.
  run_mod <- function(m, expr) {
    cur <<- m
    tryCatch(expr, error = function(e) add_error(m, conditionMessage(e)))
  }

  # ── Settings of the running module ─────────────────────────────────────────
  cfg <- function(key, default = NULL) {
    v <- settings[[cur]][[key]]
    if (is.null(v) || length(v) == 0 || (length(v) == 1 && is.na(v))) default else v
  }
  mod_traits <- function() {                          # the traits ticked for this module (all when nothing was given)
    v <- settings[[cur]]$traits
    if (is.null(v)) return(traits_all)
    t <- intersect(as.character(v), traits_all)
    if (!length(t)) stop("Select at least one trait.", call. = FALSE)
    t
  }
  mod_trait <- function(key = "trait", from = traits_all) {   # one trait (the first one when none is given)
    v <- settings[[cur]][[key]]
    if (is.null(v) || !(v[1] %in% from)) from[1] else v[1]
  }
  draw <- function() isTRUE(cfg("plots", TRUE))
  fig <- function(name, expr) if (draw()) safe(name, expr)          # figures can be switched off per module

  # ── State of earlier runs in this folder: tables for the workbook, ANOVA tables, reusable results ──
  state_file <- file.path(outdir, "state.rds")
  st <- if (isTRUE(accumulate) && file.exists(state_file)) tryCatch(readRDS(state_file), error = function(e) NULL)
  sheets   <- st$sheets %||% list()
  ta       <- st$ta                                   # ANOVA tables of the last ANOVA run
  sed_keep <- st$sed                                  # standard errors of differences of the last BLUE run
  h2_keep  <- st$h2                                   # heritability already computed (same data settings)

  # Heritability of the requested traits; only traits that were not done before are fitted
  h2_cached <- function(trs) {
    have <- h2_keep
    miss <- setdiff(trs, if (is.null(have)) character(0) else unique(have$Trait))
    if (length(miss)) {
      hh <- safe("heritability_analysis", heritability_all(working_data, miss, design))
      if (!is.null(hh)) have <- dplyr::bind_rows(have, hh)
      h2_keep <<- have
    }
    if (is.null(have)) NULL else have[have$Trait %in% trs, , drop = FALSE]
  }

  # ── Output writers ─────────────────────────────────────────────────────────
  # A readable table as CSV (UTF-8 with BOM, so Excel shows the symbols); also kept for the Excel workbook
  write_table <- function(df, name, sheet, label, workbook = TRUE) {
    readr::write_excel_csv(df, file.path(outdir, paste0(name, ".csv")), na = "")        # empty cells, not "NA"
    out$tables[[name]] <<- file.path(outdir, paste0(name, ".csv"))
    if (workbook) sheets[[sheet]] <<- list(df = df, label = label)
  }

  # Interactive HTML (plotly), a 300-dpi PNG and a vector PDF, the last two at manuscript size (inches).
  # width_px = NULL: the interactive plot fills its frame, so it also fits a half-width card.
  write_plot_out <- function(p_gg, name, height_px = 600, width_px = NULL, tooltip = NULL, size_in = FIG_STD) {
    if (!inherits(p_gg, "ggplot")) stop("the figure ", name, " could not be built")     # never write an empty figure
    html_fp <- file.path(plots_dir, paste0(name, ".html")); png_fp <- file.path(plots_dir, paste0(name, ".png"))
    pdf_fp <- file.path(plots_dir, paste0(name, ".pdf"))
    widget <- if (!is.null(tooltip)) ggplotly(p_gg, tooltip = tooltip, height = height_px, width = width_px)
              else ggplotly(p_gg, height = height_px, width = width_px)
    save_plotly(widget, html_fp)
    save_figure(p_gg, png_fp, pdf_fp, size_in)
    out$plots[[name]] <<- list(html = html_fp, png = if (file.exists(png_fp)) png_fp else NULL,
                               pdf = if (file.exists(pdf_fp)) pdf_fp else NULL)
  }
  # Base-graphics figures (corrplot, dendrogram): same 300-dpi PNG and PDF
  write_png_out <- function(draw_fn, name, size_in = FIG_SQUARE) {
    png_fp <- file.path(plots_dir, paste0(name, ".png")); pdf_fp <- file.path(plots_dir, paste0(name, ".pdf"))
    png_device(png_fp, size_in[1], size_in[2])
    tryCatch(draw_fn(), finally = dev.off())
    pdf(pdf_fp, width = size_in[1], height = size_in[2])
    tryCatch(draw_fn(), finally = dev.off())
    out$plots[[name]] <<- list(html = NULL, png = png_fp, pdf = pdf_fp)
  }
  heat_colors <- CORR_COLORS

  # Genotype means of the module's traits (one row per accession) for correlation, PCA and clustering
  get_mv <- function() {
    trs <- mod_traits()
    if (length(trs) < 2) stop("Select at least two traits.", call. = FALSE)
    mv <- multivariate_data(working_data, trs)
    if (is.null(mv)) stop("Needs at least two traits that vary and five accessions with complete data for the selected traits.", call. = FALSE)
    if (length(mv$dropped)) add_note(paste("Left out (no data in these environments, or no variation):", paste(mv$dropped, collapse = ", ")))
    mv
  }

  # ── Data notes and design check (cheap, always) ────────────────────────────
  add_note(prep$notes)
  dres <- safe("design_check", design_report(working_data, design))
  if (!is.null(dres)) {
    add_warning(dres$notes)
    write_table(tidy_design_check(dres$table, sed_keep), "design_check", "Design_Check",
                "Layout of every environment, design warnings and standard errors of differences (augmented design)")
  }

  # Removes what an earlier run of the module wrote (files and workbook sheets) before it runs again
  clear_outputs <- function(m) {
    o <- MODULE_OUTPUTS[[m]]
    unlink(file.path(outdir, paste0(o$tables, ".csv")))
    unlink(file.path(plots_dir, c(outer(o$plots, c("html", "png", "pdf"), paste, sep = "."))))
    for (sh in o$sheets) sheets[[sh]] <<- NULL
    if (m == "anova") ta <<- NULL
    if (m == "blues") {                              # the design table loses the standard errors of the old BLUE run
      sed_keep <<- NULL
      if (!is.null(dres)) write_table(tidy_design_check(dres$table, NULL), "design_check", "Design_Check",
                                      "Layout of every environment, design warnings and standard errors of differences (augmented design)")
    }
  }

  # ── The requested analyses, in tab order. Each module function runs in this environment. ──
  for (m in modules) {
    clear_outputs(m)
    f <- MODULE_FUNS[[m]]
    environment(f) <- environment()
    run_mod(m, f())
  }

  # ── Save what later runs and the workbooks need ────────────────────────────
  summary_sheet <- list(df = data.frame(
    Item = c("Experimental design", "Observations", "Variables", "Accessions", "Environments", "Traits",
             "Active environments", "Date"),
    Value = c(DESIGN_LABELS[[design]], nrow(working_data), ncol(working_data), length(unique(working_data$Genotype)),
              length(unique(working_data$Environment)), length(traits_all),
              if (identical(analysis_mode, "specific")) paste(unique(as.character(working_data$Environment)), collapse = ", ") else "All",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"))),
    label = "Dataset, design and run settings")
  saveRDS(list(sheets = sheets, ta = ta, sed = sed_keep, h2 = h2_keep, summary = summary_sheet), state_file)
  if (isTRUE(workbook)) {
    files <- export_workbook(outdir)
    out$excel_all <- files$all; out$excel_anova <- files$anova
  }
  out
}
