# =============================================================================
# 2. DATA VALIDATION
# =============================================================================
# Reads a .csv / .xlsx / .xls file; column names of the design columns are normalised (case, spaces)
read_data_file <- function(path) {
  if (!is.character(path) || length(path) != 1 || !file.exists(path)) stop("Data file not found.", call. = FALSE)
  ext <- tolower(tools::file_ext(path))
  if (!ext %in% c("csv", "xlsx", "xls")) stop("Unsupported file type: use .csv, .xlsx or .xls.", call. = FALSE)
  data <- tryCatch(
    if (ext == "csv") read_csv_flexible(path) else as.data.frame(readxl::read_excel(path)),
    error = function(e) stop("Could not read the data file: ", conditionMessage(e), call. = FALSE))
  names(data) <- canonical_names(names(data))
  for (cl in intersect(DESIGN_COLS, names(data)))              # "G1 " and "G1" are the same genotype
    if (is.character(data[[cl]])) data[[cl]] <- trimws(data[[cl]])
  data
}

# CSV with a comma, semicolon or tab as separator (decimal comma with semicolons), with or without a byte-order mark,
# UTF-8 or Latin-1 text
read_csv_flexible <- function(path) {
  first <- readLines(path, n = 1, warn = FALSE, encoding = "UTF-8")
  seps <- c(",", ";", "\t")
  sep <- seps[which.max(vapply(seps, function(s) lengths(regmatches(first, gregexpr(s, first, fixed = TRUE))), 1L))]
  rd <- function(enc) read.csv(path, sep = sep, dec = if (sep == ";") "," else ".", stringsAsFactors = FALSE,
                               check.names = FALSE, strip.white = TRUE, fileEncoding = enc)
  tryCatch(suppressWarnings(rd("UTF-8-BOM")), error = function(e) rd("latin1"))
}

validate_data <- function(data, design = "augmented") {
  validation <- list(errors = character(0), warnings = character(0))
  missing_cols <- setdiff(design_required_cols(design), names(data))
  if (length(missing_cols) > 0)
    validation$errors <- c(validation$errors,
                           paste0("Missing columns for the ", DESIGN_LABELS[[design]], " design: ",
                                  paste(missing_cols, collapse = ", ")))
  for (col in intersect(DESIGN_COLS, names(data))) {
    if (!is.factor(data[[col]])) {
      data[[col]] <- as.factor(data[[col]])
      validation$warnings <- c(validation$warnings, paste(col, "converted to factor"))
    }
  }
  if ("Genotype" %in% names(data) && length(unique(data$Genotype)) < 3)
    validation$errors <- c(validation$errors, "Need at least 3 genotypes")
  if ("Environment" %in% names(data) && length(unique(data$Environment)) < 2)
    validation$errors <- c(validation$errors, "Need at least 2 environments")
  if (length(find_traits(data)) == 0)
    validation$errors <- c(validation$errors,
      "No numeric trait columns found: every column other than Genotype, Environment, Rep and Block must be a trait, and a trait column must be numeric")
  list(data = data, validation = validation)
}

# =============================================================================
# 3. EXPERIMENTAL-DESIGN MODELS
# =============================================================================

# Plot-level data for one trait with the design's blocking factors. Replicates and blocks are
# nested: block 1 of rep 1 is not block 1 of rep 2, and rep 1 of one environment is not rep 1
# of another (EnvRep / EnvBlock / RepBlock / EnvRepBlock).
model_df <- function(d, trait, design) {
  keep <- intersect(DESIGN_COLS, names(d))
  df <- as.data.frame(d[, c(keep, trait), drop = FALSE])
  names(df)[names(df) == trait] <- "Value"
  df <- df[!is.na(df$Value), , drop = FALSE]
  for (cl in keep) df[[cl]] <- droplevels(as.factor(df[[cl]]))   # drop levels that only had missing values
  df$EnvRep <- droplevels(interaction(df$Environment, df$Rep, drop = TRUE))
  if ("Block" %in% keep) {
    df$RepBlock    <- droplevels(interaction(df$Rep, df$Block, drop = TRUE))
    df$EnvBlock    <- droplevels(interaction(df$Environment, df$Block, drop = TRUE))
    df$EnvRepBlock <- droplevels(interaction(df$Environment, df$Rep, df$Block, drop = TRUE))
  }
  df
}

# Model formulas for each design.
#   kind: aov | lmer_geno | lmer_env | lmer_var | lmer_var_single | lm_pooled
#   across_env: blocking factors are nested in environments (all designs, including augmented)
#   fixed: extra fixed terms for the variance-component models (environments, checks)
#   Rep / Block terms are left out only when the data has a single Rep or Block level.
design_formula <- function(design, kind, across_env = FALSE, interaction = FALSE, df = NULL, fixed = character(0)) {
  many <- function(col) is.null(df) || !(col %in% names(df)) || nlevels(droplevels(df[[col]])) > 1
  nest <- function(col) if (across_env) paste0("Env", col) else col

  blk_fixed <- switch(design, rcbd = "Rep", alpha = c("Rep", "RepBlock"),
                      c(if (many("Rep")) "Rep", if (many("Block")) "Block"))
  blk_random <- switch(design,
    rcbd  = nest("Rep"),
    alpha = if (across_env) c("EnvRep", "EnvRepBlock") else c("Rep", "RepBlock"),
    c(if (many("Rep")) nest("Rep"), if (many("Block")) nest("Block")))
  rnd <- if (length(blk_random)) paste0("(1 | ", blk_random, ")") else character(0)

  pooled_blocks <- switch(design,
    rcbd  = c("Environment", "EnvRep"),
    alpha = c("Environment", "EnvRep", "EnvRepBlock"),
    c("Environment", if (many("Rep")) "EnvRep", if (many("Block")) "EnvBlock"))

  terms <- switch(kind,
    aov       = if (design == "augmented") c("Genotype", blk_fixed) else c(blk_fixed, "Genotype"),
    lmer_geno = c("Genotype", rnd),
    lmer_env  = c("Environment", "Genotype", if (interaction) "Environment:Genotype", rnd),
    lmer_var  = c(if (length(fixed)) fixed else "1", "(1 | Genotype)", "(1 | Genotype:Environment)", rnd),
    lmer_var_single = c(if (length(fixed)) fixed else "1", "(1 | Genotype)", rnd),
    lm_pooled = c(pooled_blocks, "Genotype", if (interaction) "Genotype:Environment"))
  paste("Value ~", paste(terms, collapse = " + "))
}

# ── Augmented RCBD: classical (Federer 1956) analysis, one environment ───────
# Blocks are fixed and their effects are estimated from the checks (genotypes repeated in the
# environment); new entries are adjusted for the block they sit in. A block is the Block column,
# or Rep + Block when the trial has several replicates. Same analysis as agricolae::DAU.test.
augmented_fit <- function(df) {
  blk <- df$Block
  if (nlevels(droplevels(df$Rep)) > 1 && "RepBlock" %in% names(df)) {
    # Rep + Block only when every such group holds a check; otherwise Rep is just a label
    # (for example the repeats of the checks) and the blocks are the Block column
    cnt <- table(as.character(df$Genotype))
    isc <- as.character(df$Genotype) %in% names(cnt)[cnt > 1]
    if (all(tapply(isc, df$RepBlock, any))) blk <- df$RepBlock
  }
  df$Blk <- droplevels(blk)
  df$Genotype <- droplevels(df$Genotype)
  df
}

augmented_block_note <- function(df, env) {
  if (nlevels(droplevels(df$Rep)) < 2 || !("RepBlock" %in% names(df))) return(NULL)
  f <- augmented_fit(df)
  if (nlevels(f$Blk) == nlevels(droplevels(df$RepBlock))) return(NULL)
  sprintf("%s: Rep was not used to define blocks, because some Rep + Block groups hold no check. The %d blocks are the values of the Block column.",
          env, nlevels(f$Blk))
}

# The classical augmented analysis is possible when genotypes can be told apart from blocks (checks repeated in
# >= 2 blocks, error df > 0). Returns list(df, fit): the data with its blocks and the fitted Value ~ Blk + Genotype
# (reused by the ANOVA and the BLUEs, so it is fitted once), or NULL when the analysis is not possible.
augmented_prepare <- function(df) {
  tryCatch({
    d <- augmented_fit(df)
    if (nlevels(d$Blk) < 2) {
      NULL
    } else {
      m <- lm(Value ~ Blk + Genotype, data = d)
      if (anyNA(coef(m)) || df.residual(m) <= 0) NULL else list(df = d, fit = m)
    }
  }, error = function(e) NULL)
}

augmented_anova <- function(df, prep = NULL) {
  df <- augmented_fit(df)
  y <- df$Value; n <- length(y)
  cnt <- table(df$Genotype); chk <- names(cnt)[cnt > 1]
  isc <- as.character(df$Genotype) %in% chk
  n_chk <- length(chk); n_new <- nlevels(df$Genotype) - n_chk
  a1 <- anova(if (is.null(prep)) lm(Value ~ Blk + Genotype, data = df) else prep$fit)     # treatments adjusted for blocks
  a2 <- anova(lm(Value ~ Genotype + Blk, data = df))     # treatments ignoring blocks
  edf <- a1["Residuals", "Df"]; mse <- a1["Residuals", "Mean Sq"]

  # Among-check sums of squares: adjusted for blocks (part A) and ignoring blocks (part B)
  ss_chk_adj <- 0; ss_chk_raw <- 0
  if (n_chk > 1) {
    cd <- droplevels(df[isc, ])
    ss_chk_raw <- anova(lm(Value ~ Genotype, data = cd))["Genotype", "Sum Sq"]
    ss_chk_adj <- if (nlevels(cd$Blk) > 1) anova(lm(Value ~ Blk + Genotype, data = cd))["Genotype", "Sum Sq"] else ss_chk_raw
  }

  rows <- list()
  add <- function(part, term, d, ss, test = TRUE) {
    f <- if (test && d > 0 && edf > 0) (ss / d) / mse else NA
    rows[[length(rows) + 1]] <<- data.frame(
      Part = part, term = term, df = d, sumsq = ss, meansq = ss / d, statistic = f,
      p.value = if (is.na(f)) NA else pf(f, d, edf, lower.tail = FALSE), stringsAsFactors = FALSE)
  }
  pa <- "Treatments adjusted for blocks"; pb <- "Treatments ignoring blocks"
  add(pa, "Block (ignoring treatments)", a1["Blk", "Df"], a1["Blk", "Sum Sq"], test = FALSE)
  add(pa, "Treatment (adjusted for blocks)", a1["Genotype", "Df"], a1["Genotype", "Sum Sq"])
  if (n_chk > 1 && n_new > 0) {
    add(pa, "   Among checks", n_chk - 1, ss_chk_adj)
    add(pa, "   New entries, and new entries vs checks", a1["Genotype", "Df"] - (n_chk - 1),
        a1["Genotype", "Sum Sq"] - ss_chk_adj)
  }
  add(pa, "Residuals", edf, a1["Residuals", "Sum Sq"], test = FALSE)
  add(pb, "Treatment (ignoring blocks)", a2["Genotype", "Df"], a2["Genotype", "Sum Sq"], test = FALSE)
  if (n_chk > 1 && n_new > 1) {
    ss_cvn <- sum(y[isc])^2 / sum(isc) + sum(y[!isc])^2 / sum(!isc) - sum(y)^2 / n
    add(pb, "   Among checks", n_chk - 1, ss_chk_raw)
    add(pb, "   Among new entries", n_new - 1, a2["Genotype", "Sum Sq"] - ss_cvn - ss_chk_raw)
    add(pb, "   Checks vs new entries", 1, ss_cvn)
  }
  add(pb, "Block (adjusted for treatments)", a2["Blk", "Df"], a2["Blk", "Sum Sq"])
  add(pb, "Residuals", edf, a1["Residuals", "Sum Sq"], test = FALSE)
  out <- do.call(rbind, rows)
  num <- c("sumsq", "meansq", "statistic", "p.value")
  out[, num] <- lapply(out[, num], function(v) signif(v, 7))
  out
}

# Adjusted means (BLUEs) of every genotype with blocks fixed, and the standard errors of the
# differences between them, grouped as in agricolae::DAU.test.
augmented_means <- function(df, prep = NULL) {
  df <- augmented_fit(df)
  m <- if (is.null(prep)) lm(Value ~ Blk + Genotype, data = df) else prep$fit
  # Blk is a nuisance factor: its levels are averaged over exactly as before, but without building the full
  # Blk x Genotype grid (which stops emmeans above ~10,000 rows, i.e. ~500 accessions in 20 blocks)
  em <- emmeans(m, ~Genotype, nuisance = "Blk")
  tab <- as.data.frame(em)
  V <- as.matrix(vcov(em))
  g <- as.character(tab$Genotype)
  cnt <- table(as.character(df$Genotype)); is_chk <- cnt[g] > 1
  blk <- vapply(g, function(x) as.character(df$Blk[df$Genotype == x][1]), "")
  idx <- which(upper.tri(V), arr.ind = TRUE)
  sed <- sqrt(diag(V)[idx[, 1]] + diag(V)[idx[, 2]] - 2 * V[idx])
  ci <- is_chk[idx[, 1]]; cj <- is_chk[idx[, 2]]
  grp <- ifelse(ci & cj, "Two checks",
         ifelse(!ci & !cj,
                ifelse(blk[idx[, 1]] == blk[idx[, 2]], "Two new entries, same block", "Two new entries, different blocks"),
                "A new entry and a check"))
  sed_tab <- do.call(rbind, lapply(unique(grp), function(k) {
    v <- sed[grp == k]
    data.frame(Comparison = k, Pairs = length(v), SED_mean = round(mean(v, na.rm = TRUE), 4),
               SED_min = round(min(v, na.rm = TRUE), 4), SED_max = round(max(v, na.rm = TRUE), 4),
               stringsAsFactors = FALSE)
  }))
  list(blues = tab, sed = sed_tab)
}

# Checks that the layout looks like the chosen design; returns a table (one row per
# environment) and a vector of plain-language notes.
design_report <- function(d, design) {
  notes <- character(0)
  rows <- lapply(unique(as.character(d$Environment)), function(env) {
    e <- d[as.character(d$Environment) == env, , drop = FALSE]
    nrep <- length(unique(e$Rep))
    has_block <- "Block" %in% names(e)
    cnt <- table(as.character(e$Genotype))
    blk_n <- NA; blk_size <- NA
    if (has_block) {
      key <- paste(e$Rep, e$Block, sep = "/")
      if (design == "augmented") {
        tmp <- data.frame(Genotype = factor(e$Genotype), Rep = factor(e$Rep), Block = factor(e$Block), Value = 1)
        tmp$RepBlock <- droplevels(interaction(tmp$Rep, tmp$Block, drop = TRUE))
        key <- as.character(augmented_fit(tmp)$Blk)
        has_chk <- tapply(as.character(e$Genotype) %in% names(cnt)[cnt > 1], key, any)
        if (any(!has_chk) && any(cnt > 1))
          notes <<- c(notes, sprintf("%s: %d block(s) hold no check (%s), so their effect cannot be estimated from checks.",
                                     env, sum(!has_chk), paste(names(has_chk)[!has_chk], collapse = ", ")))
      }
      blk_n <- length(unique(key))
      bs <- as.integer(table(key))
      blk_size <- if (min(bs) == max(bs)) as.character(min(bs)) else paste0(min(bs), "-", max(bs))
    }
    if (nrep < 2 && design != "augmented")
      notes <<- c(notes, sprintf("%s: only one replicate; the %s model needs at least two, so this environment is skipped in the per-environment models.",
                                 env, DESIGN_LABELS[[design]]))
    if (design == "rcbd" && any(cnt != nrep))
      notes <<- c(notes, sprintf("%s: %d genotype(s) are not present exactly once in every replicate, which a complete-block design expects.",
                                 env, sum(cnt != nrep)))
    if (design == "alpha") {
      bpr <- tapply(paste(e$Rep, e$Block), e$Rep, function(z) length(unique(z)))
      if (any(bpr < 2))
        notes <<- c(notes, sprintf("%s: some replicates have a single block; an alpha-lattice needs several incomplete blocks per replicate.", env))
      if (!is.na(blk_size) && grepl("-", blk_size))
        notes <<- c(notes, sprintf("%s: blocks differ in size (%s plots); alpha-lattice blocks normally have equal size.", env, blk_size))
      if (any(table(as.character(e$Genotype), as.character(e$Rep)) > 1))
        notes <<- c(notes, sprintf("%s: some genotypes occur more than once in a replicate.", env))
    }
    if (design == "augmented") {
      if (!any(cnt > 1))
        notes <<- c(notes, sprintf("%s: no genotype is repeated, so the checks that estimate the block effects are missing.", env))
      if (has_block && blk_n < 2)
        notes <<- c(notes, sprintf("%s: only one block; block effects cannot be estimated.", env))
    }
    data.frame(Environment = env, Plots = nrow(e), Genotypes = length(cnt), Reps = nrep, Blocks = blk_n,
               Block_Size = blk_size, Min_Replication = as.integer(min(cnt)), Max_Replication = as.integer(max(cnt)),
               Repeated_Genotypes = sum(cnt > 1), Single_Plot_Genotypes = sum(cnt == 1), stringsAsFactors = FALSE)
  })
  list(table = do.call(rbind, rows), notes = unique(notes))
}
