# ============================================================
# QTL-Mapper v7 ENHANCED — FIXED VERSION
# Fixes:
#  1. Added colourpicker to required_packages
#  2. Customization controls now fully wired to plots & downloads
#  3. make_chromosome_ideogram_reference accepts color/size params
#  4. Manhattan plot height/width driven by sliders (uiOutput)
#  5. Ideogram width now respected in renderPlot
#  6. PNG downloads use custom dimensions & styles
#  7. threshold_style mapped correctly to base R lty values
#  8. compare_bar PNG uses custom dimensions
#  9. All colourInput references safe-guarded with defaults
# ============================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

# Only the packages qtl_engine.R / qtl_cli.R actually call.
# (shiny, DT, plotly, shinyjs, colourpicker were Shiny-UI-only
#  leftovers from the original app and are never used here.)
required_packages <- c(
  "randomForest", "caret", "ggplot2",
  "dplyr", "glmnet", "scales", "e1071", "xgboost",
  "writexl", "jsonlite"
)

missing_packages <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them once ahead of time (e.g. in your Docker image / setup script) with:\n",
    "  Rscript -e 'install.packages(c(",
    paste0('"', missing_packages, '"', collapse = ", "),
    "))'\n",
    "Live install.packages() inside a request is disabled because it can hang ",
    "indefinitely (e.g. no outbound internet, no CRAN binary) with no visible error."
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

# ============================================================
# 6 MODELS
# ============================================================
MAPPING_MODELS <- list(
  "SIM"     = list(label = "Simple Interval Mapping (SIM)",    desc = "Single-marker regression; LOD = -log10(p)"),
  "CIM"     = list(label = "Composite Interval Mapping (CIM)", desc = "Cartographer-like CIM with cofactors"),
  "MQM"     = list(label = "Multiple QTL Model (MQM)",         desc = "Forward stepwise multi-QTL selection"),
  "RF"      = list(label = "Random Forest QTL",                desc = "Non-parametric; permutation importance"),
  "SVR"     = list(label = "Support Vector Regression",        desc = "Kernel-based; permutation importance"),
  "XGBOOST" = list(label = "XGBoost QTL",                     desc = "Gradient boosting; gain importance")
)

MODEL_COLORS <- c(
  "SIM"     = "#CC0000",
  "CIM"     = "#FF3333",
  "MQM"     = "#990000",
  "RF"      = "#FF6600",
  "SVR"     = "#CC3300",
  "XGBOOST" = "#FF0066"
)

# ============================================================
# HELPERS
# ============================================================
cn <- function(x) trimws(gsub('[\r\n\"]', '', as.character(x)))

calc_genetic_distance <- function(marker1, marker2, marker_info) {
  m1 <- marker_info[marker_info$Marker == marker1, ]
  m2 <- marker_info[marker_info$Marker == marker2, ]
  if (nrow(m1) == 0 || nrow(m2) == 0) return(Inf)
  if (m1$Chromosome[1] == m2$Chromosome[1]) return(abs(m1$Position[1] - m2$Position[1]))
  return(Inf)
}

# Safe colour getter with fallback
get_color <- function(input_val, fallback) {
  if (is.null(input_val) || nchar(trimws(input_val)) == 0) return(fallback)
  input_val
}

# ============================================================
# GENOTYPE LOADER
# ============================================================
load_genotype_file <- function(filepath) {
  lines <- readLines(filepath, warn = FALSE)
  lines <- gsub("\r", "", lines)
  lines <- lines[nchar(trimws(lines)) > 0]
  if (length(lines) < 4) stop("Need >= 4 rows: header + chr + pos + >= 1 sample")

  header     <- cn(strsplit(lines[1], ",")[[1]])
  chr_row    <- cn(strsplit(lines[2], ",")[[1]])
  pos_row    <- cn(strsplit(lines[3], ",")[[1]])
  data_lines <- lines[-(1:3)]

  gd <- read.csv(text = paste(data_lines, collapse = "\n"),
                 header = FALSE, stringsAsFactors = FALSE, quote = '"', fill = TRUE)
  if (ncol(gd) < length(header)) {
    for (k in (ncol(gd)+1):length(header)) gd[[paste0("V",k)]] <- NA
  }
  if (ncol(gd) > length(header)) gd <- gd[, 1:length(header)]
  colnames(gd) <- header
  rownames(gd) <- cn(as.character(gd[, 1]))
  gd <- gd[, -1, drop = FALSE]
  for (col in colnames(gd)) gd[[col]] <- suppressWarnings(as.numeric(as.character(gd[[col]])))

  marker_info <- data.frame(
    Marker     = colnames(gd),
    Chromosome = cn(chr_row[-1][seq_along(colnames(gd))]),
    Position   = suppressWarnings(as.numeric(cn(pos_row[-1][seq_along(colnames(gd))]))),
    stringsAsFactors = FALSE
  )
  marker_info <- marker_info[!is.na(marker_info$Position) & marker_info$Chromosome != "", ]
  valid_cols  <- colnames(gd)[colnames(gd) %in% marker_info$Marker]
  gd          <- gd[, valid_cols, drop = FALSE]

  list(geno = gd, marker_info = marker_info,
       n_samples = nrow(gd), n_markers = ncol(gd),
       chromosomes = sort(unique(marker_info$Chromosome)))
}

# ============================================================
# SIM
# ============================================================
qtl_scan_sim <- function(X_valid, y_valid, marker_info) {
  results_list <- list()
  for (m in colnames(X_valid)) {
    geno      <- X_valid[, m]
    valid_idx <- !is.na(geno) & !is.na(y_valid)
    geno_v    <- geno[valid_idx]; pheno_v <- y_valid[valid_idx]
    if (length(unique(geno_v)) < 2 || length(pheno_v) < 8) next
    m_info <- marker_info[marker_info$Marker == m, ]
    if (nrow(m_info) == 0) next
    m_chr <- as.character(m_info$Chromosome[1]); m_pos <- m_info$Position[1]
    if (is.na(m_pos)) next
    tryCatch({
      fit <- lm(pheno_v ~ geno_v)
      ct  <- coef(summary(fit))
      if (!"geno_v" %in% rownames(ct)) next
      p_val <- ct["geno_v", "Pr(>|t|)"]; lod <- -log10(max(p_val, 1e-300))
      a_code <- geno_v - 1; d_code <- as.numeric(geno_v == 1)
      has_het <- length(unique(d_code)) > 1
      lm_eff <- if (has_het) lm(pheno_v ~ a_code + d_code) else lm(pheno_v ~ a_code)
      pve_m  <- max(0, summary(lm_eff)$r.squared * 100)
      ce <- coef(lm_eff); se_t <- coef(summary(lm_eff))
      a_eff <- ce["a_code"]
      se_a  <- if ("a_code" %in% rownames(se_t)) se_t["a_code", "Std. Error"] else NA
      d_eff <- if (has_het) ce["d_code"] else NA
      d_rat <- if (!is.na(d_eff) && !is.na(a_eff) && a_eff != 0) abs(d_eff / a_eff) else NA
      dom_cl <- if (is.na(d_rat)) "Additive"
                else if (d_rat > 1) "Overdominance"
                else if (d_rat > 0.5) "Dominance"
                else if (d_rat > 0.2) "Partial Dom." else "Additive"
      results_list[[m]] <- data.frame(
        Marker=m, Chromosome=m_chr, Position=m_pos,
        LOD=lod, P_Value=p_val,
        PVE_Marker_Percent=pve_m, PVE_CIM_Percent=pve_m,
        Additive_Effect=a_eff, SE_Additive=se_a,
        CI_Lower=ifelse(is.na(se_a),NA,a_eff-1.96*se_a),
        CI_Upper=ifelse(is.na(se_a),NA,a_eff+1.96*se_a),
        Dominance_Effect=d_eff, Dominance_Ratio=d_rat, Dominance_Class=dom_cl,
        N_Cofactors=0, N_Samples=length(pheno_v), Method="SIM",
        stringsAsFactors=FALSE)
    }, error=function(e) {})
  }
  if (length(results_list) == 0) return(data.frame(Message="No results from SIM"))
  res <- do.call(rbind, results_list)
  res[order(res$Chromosome, res$Position), ]
}

# ============================================================
# CIM
# ============================================================
qtl_scan_cim_improved <- function(X_valid, y_valid, marker_info,
                                  window_size = 10, n_cofactors = 5) {
  cat("\n=== CIM ===\n")
  initial_pvals <- sapply(colnames(X_valid), function(m) {
    g <- X_valid[, m]; valid_idx <- !is.na(g) & !is.na(y_valid)
    if (sum(valid_idx) < 8 || length(unique(g[valid_idx])) < 2) return(1)
    tryCatch(coef(summary(lm(y_valid[valid_idx] ~ g[valid_idx])))[2, 4], error=function(e) 1)
  })
  sorted_markers <- names(sort(initial_pvals))
  cofactors <- character(0)
  min_distance <- 20
  for (candidate in sorted_markers) {
    if (length(cofactors) >= n_cofactors) break
    if (initial_pvals[candidate] > 0.01) break
    too_close <- FALSE
    if (length(cofactors) > 0) {
      for (existing_cof in cofactors) {
        if (calc_genetic_distance(candidate, existing_cof, marker_info) < min_distance) {
          too_close <- TRUE; break
        }
      }
    }
    if (!too_close) cofactors <- c(cofactors, candidate)
  }
  results_list <- list()
  for (m in colnames(X_valid)) {
    geno <- X_valid[, m]; valid_idx <- !is.na(geno) & !is.na(y_valid)
    geno_v <- geno[valid_idx]; pheno_v <- y_valid[valid_idx]
    if (length(unique(geno_v)) < 2 || length(pheno_v) < 8) next
    m_info <- marker_info[marker_info$Marker == m, ]
    if (nrow(m_info) == 0) next
    m_chr <- as.character(m_info$Chromosome[1]); m_pos <- m_info$Position[1]
    if (is.na(m_pos)) next
    chr_markers <- marker_info[marker_info$Chromosome == m_chr, ]
    chr_markers <- chr_markers[order(chr_markers$Position), ]
    focal_idx   <- which(chr_markers$Marker == m)
    exclude_markers <- character(0)
    if (length(focal_idx) > 0) {
      ws <- max(1, focal_idx - window_size); we <- min(nrow(chr_markers), focal_idx + window_size)
      exclude_markers <- chr_markers$Marker[ws:we]
    }
    active_cof <- setdiff(cofactors, c(m, exclude_markers))
    active_cof <- active_cof[active_cof %in% colnames(X_valid)]
    tryCatch({
      df <- data.frame(y = pheno_v, focal = geno_v)
      if (length(active_cof) > 0) {
        cof_matrix <- X_valid[valid_idx, active_cof, drop=FALSE]
        cof_matrix[is.na(cof_matrix)] <- 0
        df <- cbind(df, cof_matrix)
      }
      df <- df[complete.cases(df), ]
      if (nrow(df) < 8) next
      if (length(active_cof) > 0) {
        fit_full    <- lm(as.formula(paste("y ~ focal +", paste(active_cof, collapse=" + "))), data=df)
        fit_reduced <- lm(as.formula(paste("y ~", paste(active_cof, collapse=" + "))), data=df)
      } else {
        fit_full <- lm(y ~ focal, data=df); fit_reduced <- lm(y ~ 1, data=df)
      }
      anova_result <- anova(fit_reduced, fit_full)
      p_val <- anova_result[2, "Pr(>F)"]; lod <- -log10(max(p_val, 1e-300))
      ct <- coef(summary(fit_full))
      if (!"focal" %in% rownames(ct)) next
      RSS_full <- sum(residuals(fit_full)^2); RSS_reduced <- sum(residuals(fit_reduced)^2)
      pve_cim <- if (RSS_reduced > 0) max(0, (RSS_reduced - RSS_full) / RSS_reduced * 100) else NA
      a_code <- df$focal - 1; d_code <- as.numeric(df$focal == 1)
      has_het <- length(unique(d_code)) > 1
      lm_effects <- if (has_het) lm(df$y ~ a_code + d_code) else lm(df$y ~ a_code)
      pve_marker <- max(0, summary(lm_effects)$r.squared * 100)
      ce <- coef(lm_effects); se_table <- coef(summary(lm_effects))
      a_eff <- ce["a_code"]
      se_a  <- if ("a_code" %in% rownames(se_table)) se_table["a_code", "Std. Error"] else NA
      d_eff <- if (has_het) ce["d_code"] else NA
      d_ratio <- if (!is.na(d_eff) && !is.na(a_eff) && a_eff != 0) abs(d_eff / a_eff) else NA
      dom_class <- if (is.na(d_ratio)) "Additive"
                   else if (d_ratio > 1) "Overdominance"
                   else if (d_ratio > 0.5) "Dominance"
                   else if (d_ratio > 0.2) "Partial Dom." else "Additive"
      results_list[[m]] <- data.frame(
        Marker=m, Chromosome=m_chr, Position=m_pos,
        LOD=lod, P_Value=p_val,
        PVE_Marker_Percent=pve_marker, PVE_CIM_Percent=pve_cim,
        Additive_Effect=a_eff, SE_Additive=se_a,
        CI_Lower=ifelse(is.na(se_a),NA,a_eff-1.96*se_a),
        CI_Upper=ifelse(is.na(se_a),NA,a_eff+1.96*se_a),
        Dominance_Effect=d_eff, Dominance_Ratio=d_ratio, Dominance_Class=dom_class,
        N_Cofactors=length(active_cof), N_Samples=nrow(df), Method="CIM",
        stringsAsFactors=FALSE)
    }, error=function(e) cat("CIM error at", m, ":", e$message, "\n"))
  }
  if (length(results_list) == 0) return(data.frame(Message="CIM produced no results"))
  res <- do.call(rbind, results_list)
  res[order(res$Chromosome, res$Position), ]
}

# ============================================================
# MQM
# ============================================================
qtl_scan_mqm <- function(X, y, marker_info) {
  valid <- !is.na(y); X_v <- X[valid,,drop=FALSE]; y_v <- y[valid]
  p_vals <- sapply(colnames(X_v), function(m) {
    g <- X_v[,m]; v2 <- !is.na(g)
    if (length(unique(g[v2]))<2||sum(v2)<8) return(1)
    tryCatch(coef(summary(lm(y_v[v2]~g[v2])))[2,4],error=function(e) 1)
  })
  selected <- character(0); threshold <- 0.01
  for (cand in names(sort(p_vals))[1:min(20,length(p_vals))]) {
    test_model <- c(selected, cand)
    df_t <- data.frame(y=y_v, X_v[,test_model,drop=FALSE])
    df_t <- df_t[complete.cases(df_t),]
    if (nrow(df_t)<10) next
    tryCatch({
      fw <- lm(y~.,data=df_t)
      wo <- if(length(selected)>0) lm(as.formula(paste("y~",paste(selected,collapse="+"))),data=df_t) else lm(y~1,data=df_t)
      p  <- anova(wo,fw)[2,"Pr(>F)"]
      if (!is.na(p)&&p<threshold) selected <- test_model
    },error=function(e){})
  }
  if (length(selected)==0) selected <- names(sort(p_vals))[1:min(5,length(p_vals))]
  results_list <- list()
  for (m in colnames(X_v)) {
    g <- X_v[,m]; v2 <- !is.na(g)&!is.na(y_v)
    gv <- g[v2]; yv <- y_v[v2]
    if (length(unique(gv))<2||length(yv)<8) next
    mi <- marker_info[marker_info$Marker==m,]
    if (nrow(mi)==0) next
    cof_m <- setdiff(selected,m)
    tryCatch({
      df2 <- data.frame(y=yv,focal=gv)
      if (length(cof_m)>0){cm<-X_v[v2,cof_m,drop=FALSE];cm[is.na(cm)]<-0;df2<-cbind(df2,cm)}
      df2 <- df2[complete.cases(df2),]
      if (nrow(df2)<8) next
      ff <- if(length(cof_m)>0) lm(as.formula(paste("y~focal+",paste(cof_m,collapse="+"))),data=df2) else lm(y~focal,data=df2)
      fr <- if(length(cof_m)>0) lm(as.formula(paste("y~",paste(cof_m,collapse="+"))),data=df2) else lm(y~1,data=df2)
      ct2 <- coef(summary(ff))
      if (!"focal"%in%rownames(ct2)) next
      p2 <- ct2["focal","Pr(>|t|)"]; lod2 <- -log10(max(p2,1e-300))
      RF2 <- sum(residuals(ff)^2); RR2 <- sum(residuals(fr)^2)
      pvc <- if(RR2>0) max(0,(RR2-RF2)/RR2*100) else NA
      ac <- gv-1; dc2 <- as.numeric(gv==1); h2 <- length(unique(dc2))>1
      le <- if(h2) lm(yv~ac+dc2) else lm(yv~ac)
      pvm <- max(0,summary(le)$r.squared*100)
      ce2 <- coef(le); se2 <- coef(summary(le))
      ae <- ce2["ac"]; sea <- if("ac"%in%rownames(se2)) se2["ac","Std. Error"] else NA
      de <- if(h2) ce2["dc2"] else NA
      dr <- if(!is.na(de)&&!is.na(ae)&&ae!=0) abs(de/ae) else NA
      dc3 <- if(is.na(dr)) "Additive" else if(dr>1) "Overdominance" else if(dr>0.5) "Dominance" else if(dr>0.2) "Partial Dom." else "Additive"
      results_list[[m]] <- data.frame(
        Marker=m,Chromosome=as.character(mi$Chromosome[1]),Position=mi$Position[1],
        LOD=lod2,P_Value=p2,PVE_Marker_Percent=pvm,PVE_CIM_Percent=pvc,
        Additive_Effect=ae,SE_Additive=sea,
        CI_Lower=ifelse(is.na(sea),NA,ae-1.96*sea),CI_Upper=ifelse(is.na(sea),NA,ae+1.96*sea),
        Dominance_Effect=de,Dominance_Ratio=dr,Dominance_Class=dc3,
        N_Cofactors=length(cof_m),N_Samples=nrow(df2),Method="MQM",stringsAsFactors=FALSE)
    },error=function(e){})
  }
  if (length(results_list)==0) return(data.frame(Message="MQM no results"))
  res <- do.call(rbind,results_list); res[order(res$Chromosome,res$Position),]
}

# ============================================================
# RF
# ============================================================
qtl_scan_rf <- function(X, y, marker_info) {
  valid <- !is.na(y)&apply(X,1,function(r) !any(is.na(r)))
  X_v <- X[valid,,drop=FALSE]; y_v <- y[valid]
  rf_fit <- tryCatch(randomForest(x=X_v,y=y_v,ntree=500,importance=TRUE),error=function(e) NULL)
  if (is.null(rf_fit)) return(data.frame(Message="Random Forest failed"))
  imp <- importance(rf_fit,type=1)[,1]
  imp_df <- data.frame(Marker=names(imp),Importance=imp)
  imp_df <- imp_df[imp_df$Importance>0,]
  if (nrow(imp_df)==0) return(data.frame(Message="No important markers"))
  imp_df <- merge(imp_df,marker_info,by="Marker",all.x=TRUE)
  imp_df <- imp_df[!is.na(imp_df$Position),]
  max_imp <- max(imp_df$Importance); imp_df$LOD <- imp_df$Importance/max_imp*10
  r2_vec <- sapply(imp_df$Marker,function(m){
    g<-X_v[,m];v<-!is.na(g);if(sum(v)<5) return(NA)
    summary(lm(y_v[v]~g[v]))$r.squared*100})
  imp_df$PVE_Marker_Percent<-r2_vec; imp_df$PVE_CIM_Percent<-r2_vec
  imp_df$Additive_Effect<-NA; imp_df$P_Value<-10^(-imp_df$LOD)
  imp_df$SE_Additive<-NA; imp_df$CI_Lower<-NA; imp_df$CI_Upper<-NA
  imp_df$Dominance_Effect<-NA; imp_df$Dominance_Ratio<-NA
  imp_df$Dominance_Class<-"Additive"; imp_df$N_Cofactors<-0
  imp_df$N_Samples<-sum(valid); imp_df$Method<-"RF"
  imp_df <- imp_df[order(imp_df$Chromosome,imp_df$Position),]
  imp_df[,c("Marker","Chromosome","Position","LOD","P_Value","PVE_Marker_Percent","PVE_CIM_Percent",
            "Additive_Effect","SE_Additive","CI_Lower","CI_Upper","Dominance_Effect",
            "Dominance_Ratio","Dominance_Class","N_Cofactors","N_Samples","Method")]
}

# ============================================================
# SVR
# ============================================================
qtl_scan_svr <- function(X, y, marker_info) {
  valid <- !is.na(y)&apply(X,1,function(r) !any(is.na(r)))
  X_v <- X[valid,,drop=FALSE]; y_v <- y[valid]
  svr_fit <- tryCatch(svm(x=X_v,y=y_v,kernel="radial",scale=TRUE),error=function(e) NULL)
  if (is.null(svr_fit)) return(data.frame(Message="SVR failed"))
  baseline_pred  <- predict(svr_fit,X_v)
  baseline_error <- mean((y_v-baseline_pred)^2)
  importance <- sapply(colnames(X_v),function(m){
    X_perm <- X_v; X_perm[,m] <- sample(X_perm[,m])
    perm_pred <- predict(svr_fit,X_perm)
    perm_error <- mean((y_v-perm_pred)^2)
    (perm_error-baseline_error)/baseline_error})
  imp_df <- data.frame(Marker=names(importance),Importance=importance)
  imp_df <- imp_df[imp_df$Importance>0,]
  if (nrow(imp_df)==0) return(data.frame(Message="No important markers"))
  imp_df <- merge(imp_df,marker_info,by="Marker",all.x=TRUE)
  imp_df <- imp_df[!is.na(imp_df$Position),]
  max_imp <- max(imp_df$Importance); imp_df$LOD <- imp_df$Importance/max_imp*10
  r2_vec <- sapply(imp_df$Marker,function(m){
    g<-X_v[,m];v<-!is.na(g);if(sum(v)<5) return(NA)
    summary(lm(y_v[v]~g[v]))$r.squared*100})
  imp_df$PVE_Marker_Percent<-r2_vec; imp_df$PVE_CIM_Percent<-r2_vec
  imp_df$Additive_Effect<-NA; imp_df$P_Value<-10^(-imp_df$LOD)
  imp_df$SE_Additive<-NA; imp_df$CI_Lower<-NA; imp_df$CI_Upper<-NA
  imp_df$Dominance_Effect<-NA; imp_df$Dominance_Ratio<-NA
  imp_df$Dominance_Class<-"Additive"; imp_df$N_Cofactors<-0
  imp_df$N_Samples<-sum(valid); imp_df$Method<-"SVR"
  imp_df <- imp_df[order(imp_df$Chromosome,imp_df$Position),]
  imp_df[,c("Marker","Chromosome","Position","LOD","P_Value","PVE_Marker_Percent","PVE_CIM_Percent",
            "Additive_Effect","SE_Additive","CI_Lower","CI_Upper","Dominance_Effect",
            "Dominance_Ratio","Dominance_Class","N_Cofactors","N_Samples","Method")]
}

# ============================================================
# XGBOOST
# ============================================================
qtl_scan_xgboost <- function(X, y, marker_info) {
  valid <- !is.na(y)&apply(X,1,function(r) !any(is.na(r)))
  X_v <- X[valid,,drop=FALSE]; y_v <- y[valid]
  dtrain <- xgb.DMatrix(data=as.matrix(X_v),label=y_v)
  xgb_fit <- tryCatch(xgb.train(
    data=dtrain,nrounds=100,objective="reg:squarederror",verbose=0,
    params=list(max_depth=3,eta=0.1,subsample=0.8,colsample_bytree=0.8)),
    error=function(e) NULL)
  if (is.null(xgb_fit)) return(data.frame(Message="XGBoost failed"))
  importance_matrix <- xgb.importance(model=xgb_fit)
  imp_df <- data.frame(Marker=importance_matrix$Feature,Importance=importance_matrix$Gain)
  imp_df <- imp_df[imp_df$Importance>0,]
  if (nrow(imp_df)==0) return(data.frame(Message="No important markers"))
  imp_df <- merge(imp_df,marker_info,by="Marker",all.x=TRUE)
  imp_df <- imp_df[!is.na(imp_df$Position),]
  max_imp <- max(imp_df$Importance); imp_df$LOD <- imp_df$Importance/max_imp*12
  r2_vec <- sapply(imp_df$Marker,function(m){
    g<-X_v[,m];v<-!is.na(g);if(sum(v)<5) return(NA)
    summary(lm(y_v[v]~g[v]))$r.squared*100})
  imp_df$PVE_Marker_Percent<-r2_vec; imp_df$PVE_CIM_Percent<-r2_vec
  imp_df$Additive_Effect<-NA; imp_df$P_Value<-10^(-imp_df$LOD)
  imp_df$SE_Additive<-NA; imp_df$CI_Lower<-NA; imp_df$CI_Upper<-NA
  imp_df$Dominance_Effect<-NA; imp_df$Dominance_Ratio<-NA
  imp_df$Dominance_Class<-"Additive"; imp_df$N_Cofactors<-0
  imp_df$N_Samples<-sum(valid); imp_df$Method<-"XGBOOST"
  imp_df <- imp_df[order(imp_df$Chromosome,imp_df$Position),]
  imp_df[,c("Marker","Chromosome","Position","LOD","P_Value","PVE_Marker_Percent","PVE_CIM_Percent",
            "Additive_Effect","SE_Additive","CI_Lower","CI_Upper","Dominance_Effect",
            "Dominance_Ratio","Dominance_Class","N_Cofactors","N_Samples","Method")]
}

# ============================================================
# MAIN DISPATCHER
# ============================================================
qtl_scan <- function(genotype_matrix, phenotype_vector, marker_info,
                     method="CIM", n_cofactors=10, window_size=10) {
  cat("\n=== QTL SCAN === Method:", method, "\n")
  valid_rows <- !is.na(phenotype_vector)
  X_valid    <- genotype_matrix[valid_rows,,drop=FALSE]
  y_valid    <- phenotype_vector[valid_rows]
  nzv        <- apply(X_valid, 2, function(x) length(unique(na.omit(x))) < 2)
  X_valid    <- X_valid[, !nzv, drop=FALSE]
  cat("Markers after NZV removal:", ncol(X_valid), "\n")
  if (method == "SIM")     return(qtl_scan_sim(X_valid, y_valid, marker_info))
  if (method == "CIM")     return(qtl_scan_cim_improved(X_valid, y_valid, marker_info, window_size, n_cofactors))
  if (method == "MQM")     return(qtl_scan_mqm(X_valid, y_valid, marker_info))
  if (method == "RF")      return(qtl_scan_rf(X_valid, y_valid, marker_info))
  if (method == "SVR")     return(qtl_scan_svr(X_valid, y_valid, marker_info))
  if (method == "XGBOOST") return(qtl_scan_xgboost(X_valid, y_valid, marker_info))
  data.frame(Message=paste("Unknown method:", method))
}

# ============================================================
# LOD SUPPORT INTERVAL
# ============================================================
add_lod_support_interval <- function(sig_df, all_df, lod_drop=1.5) {
  sig_df$Low_Flank <- sig_df$High_Flank <- NA_character_
  sig_df$Interval_Start <- sig_df$Interval_End <- sig_df$Interval_Size_cM <- NA_real_
  for (i in 1:nrow(sig_df)) {
    m <- sig_df$Marker[i]; ch <- sig_df$Chromosome[i]
    peak_lod <- sig_df$LOD[i]; peak_pos <- sig_df$Position[i]
    chr_data <- all_df[as.character(all_df$Chromosome)==as.character(ch)&!is.na(all_df$LOD)&!is.na(all_df$Position),]
    chr_data <- chr_data[order(chr_data$Position),]
    if (nrow(chr_data)==0) next
    threshold <- peak_lod - lod_drop
    interval_markers <- chr_data[chr_data$LOD >= threshold,]
    if (nrow(interval_markers) > 0) {
      interval_start <- min(interval_markers$Position, peak_pos)
      interval_end   <- max(interval_markers$Position, peak_pos)
      sig_df$Interval_Start[i]   <- interval_start
      sig_df$Interval_End[i]     <- interval_end
      sig_df$Interval_Size_cM[i] <- interval_end - interval_start
      below_markers <- chr_data[chr_data$Position < interval_start,]
      above_markers <- chr_data[chr_data$Position > interval_end,]
      sig_df$Low_Flank[i]  <- if(nrow(below_markers)>0) below_markers$Marker[nrow(below_markers)] else m
      sig_df$High_Flank[i] <- if(nrow(above_markers)>0) above_markers$Marker[1] else m
    } else {
      peak_idx <- which(chr_data$Marker==m)
      if (length(peak_idx)>0) {
        sig_df$Low_Flank[i]  <- if(peak_idx>1) chr_data$Marker[peak_idx-1] else m
        sig_df$High_Flank[i] <- if(peak_idx<nrow(chr_data)) chr_data$Marker[peak_idx+1] else m
        sig_df$Interval_Start[i] <- sig_df$Position[i]
        sig_df$Interval_End[i]   <- sig_df$Position[i]
        sig_df$Interval_Size_cM[i] <- 0
      }
    }
  }
  sig_df
}

# ============================================================
# CHROMOSOME IDEOGRAM — NOW ACCEPTS COLOR & SIZE PARAMS
# ============================================================
make_chromosome_ideogram_reference <- function(qtl_list, all_markers_df,
                                               lod_threshold = 2.0,
                                               model_colors  = MODEL_COLORS,
                                               title_str     = NULL,
                                               band_dark_col  = "#2E6B35",
                                               band_light_col = "#6AAB6E",
                                               border_col     = "#1A4A1F",
                                               qtl_size_min   = 0.7,
                                               qtl_size_max   = 2.0) {
  if (is.null(all_markers_df) || nrow(all_markers_df) == 0) return(invisible(NULL))

  chroms    <- sort(unique(as.character(all_markers_df$Chromosome)))
  n_chr     <- length(chroms)
  x_gap     <- 1.8
  bar_w     <- 0.55
  max_pos   <- max(all_markers_df$Position, na.rm=TRUE)
  model_names <- names(qtl_list)

  par(mar=c(5, 4, 4, 8), bg="white", xpd=TRUE)

  plot(0, 0, type="n",
       xlim=c(0.2, n_chr * x_gap + 0.3),
       ylim=c(-8, max_pos * 1.12),
       xaxt="n", yaxt="n", xlab="", ylab="", main="", frame.plot=FALSE, bty="n")

  y_ticks <- pretty(c(0, max_pos), n=7)
  axis(2, at=y_ticks, labels=y_ticks, cex.axis=0.72, las=1,
       col="#555555", col.axis="#333333", tcl=-0.3, lwd=0.8)
  mtext("Genetic Position (cM)", side=2, line=2.5, cex=0.82, col="#222222", font=2)

  for (ci in seq_along(chroms)) {
    chr      <- chroms[ci]
    x_center <- (ci - 1) * x_gap + 0.9

    chr_m   <- all_markers_df[as.character(all_markers_df$Chromosome)==chr, ]
    if (nrow(chr_m)==0) next
    chr_len <- max(chr_m$Position, na.rm=TRUE)
    chr_min <- min(chr_m$Position, na.rm=TRUE)

    n_bands <- max(8, round((chr_len - chr_min) / 8))
    band_h  <- (chr_len - chr_min) / n_bands

    for (bi in 1:n_bands) {
      band_bot <- chr_min + (bi-1)*band_h
      band_top <- chr_min + bi*band_h
      bcol     <- if (bi %% 2 == 1) band_dark_col else band_light_col
      rect(x_center - bar_w/2, band_bot, x_center + bar_w/2, band_top,
           col=bcol, border=NA)
    }

    rect(x_center - bar_w/2, chr_min, x_center + bar_w/2, chr_len,
         col=NA, border=border_col, lwd=1.8)

    symbols(x_center, chr_min, circles=bar_w/2, inches=FALSE,
            add=TRUE, fg=border_col, bg=band_light_col, lwd=1.2)
    symbols(x_center, chr_len, circles=bar_w/2, inches=FALSE,
            add=TRUE, fg=border_col, bg=band_dark_col, lwd=1.2)

    cen_pos <- chr_min + (chr_len - chr_min)*0.42
    polygon(c(x_center-bar_w/2, x_center, x_center+bar_w/2, x_center),
            c(cen_pos - chr_len*0.025, cen_pos,
              cen_pos - chr_len*0.025, cen_pos - chr_len*0.05),
            col="white", border=border_col, lwd=1.2)
    polygon(c(x_center-bar_w/2, x_center, x_center+bar_w/2, x_center),
            c(cen_pos + chr_len*0.025, cen_pos,
              cen_pos + chr_len*0.025, cen_pos + chr_len*0.05),
            col="white", border=border_col, lwd=1.2)

    n_models <- length(model_names)
    for (ti in seq_along(model_names)) {
      trait    <- model_names[ti]
      qtl_data <- qtl_list[[trait]]
      if (is.null(qtl_data) || "Message" %in% colnames(qtl_data)) next

      sig_qtls <- qtl_data[as.character(qtl_data$Chromosome)==chr &
                           !is.na(qtl_data$LOD) &
                           qtl_data$LOD >= lod_threshold, ]
      if (nrow(sig_qtls)==0) next

      col_m <- model_colors[trait]
      if (is.na(col_m)) col_m <- "#CC0000"

      x_off <- x_center + bar_w/2 + 0.08 + (ti-1)*0.22

      for (ri in 1:nrow(sig_qtls)) {
        pos     <- sig_qtls$Position[ri]
        lod_val <- sig_qtls$LOD[ri]

        segments(x_center + bar_w/2, pos, x_off - 0.04, pos,
                 col=col_m, lwd=1.5)
        points(x_off, pos, pch=16, col=col_m,
               cex=pmax(qtl_size_min, pmin(qtl_size_max, 0.6 + lod_val/8)))
        text(x_off + 0.15, pos,
             paste0(round(lod_val,1)),
             cex=0.55, col=col_m, font=2, adj=0)
      }
    }

    text(x_center, chr_min - max_pos*0.04, paste0("Chr", chr),
         cex=0.72, font=2, col="#1A1A1A", adj=0.5)
  }

  legend_x <- n_chr * x_gap + 0.5
  legend_y  <- max_pos * 0.95
  text(legend_x, legend_y + max_pos*0.07, "Models", cex=0.75, font=2, col="#1A1A1A", adj=0)
  for (ti in seq_along(model_names)) {
    tr    <- model_names[ti]
    col_m <- model_colors[tr]; if (is.na(col_m)) col_m <- "#CC0000"
    y_leg <- legend_y - (ti-1)*max_pos*0.07
    qtl_d <- qtl_list[[tr]]
    n_sig <- if(!is.null(qtl_d)&&!"Message"%in%colnames(qtl_d))
               sum(!is.na(qtl_d$LOD)&qtl_d$LOD>=lod_threshold,na.rm=TRUE) else 0
    points(legend_x, y_leg, pch=16, col=col_m, cex=1.2)
    text(legend_x + 0.18, y_leg,
         paste0(tr, " (", n_sig, " QTL)"),
         cex=0.65, col=col_m, font=1, adj=0)
  }

  ttl <- if (!is.null(title_str)) title_str else paste0("QTL Chromosome Map  |  LOD ≥ ", lod_threshold)
  title(main=ttl, cex.main=0.95, font.main=2, col.main="#1A1A1A", line=2)
}

# ============================================================
# CHROMOSOME-WISE MANHATTAN PLOT (PLOTLY) — FULLY CUSTOMIZABLE
# ============================================================
make_plotly_manhattan_chromosomewise <- function(scan_all, lod_threshold=2.0,
                                                  title="QTL Manhattan Plot",
                                                  nonsig_color="#3A7D44",
                                                  sig_color="#CC0000",
                                                  point_size=6,
                                                  sig_point_size=10,
                                                  threshold_line_width=1.5,
                                                  threshold_line_style="dash") {
  blank <- plot_ly() %>% layout(title=title, plot_bgcolor="white", paper_bgcolor="white")
  if (is.null(scan_all) || "Message" %in% colnames(scan_all)) return(blank)

  pd <- scan_all[!is.na(scan_all$LOD) & !is.na(scan_all$Position) & !is.na(scan_all$Chromosome), ]
  if (nrow(pd) == 0) return(blank)

  pd$Chromosome <- as.character(pd$Chromosome)
  pd$sig <- pd$LOD >= lod_threshold

  chroms <- sort(unique(pd$Chromosome))
  subplot_list <- list()

  for (chr in chroms) {
    chr_data <- pd[pd$Chromosome == chr, ]
    chr_data <- chr_data[order(chr_data$Position), ]
    nonsig   <- chr_data[!chr_data$sig, ]
    sig_data <- chr_data[chr_data$sig, ]

    p <- plot_ly()

    if (nrow(nonsig) > 0) {
      p <- p %>% add_trace(
        data = nonsig, x = ~Position, y = ~LOD,
        type = "scatter", mode = "markers",
        marker = list(size = point_size, color = nonsig_color, opacity = 0.7),
        text = ~paste0("<b>", Marker, "</b><br>Chr: ", Chromosome,
                       "<br>Pos: ", round(Position, 1), " cM<br>LOD: ", round(LOD, 2),
                       "<br>PVE: ", round(PVE_Marker_Percent, 1), "%"),
        hoverinfo = "text", showlegend = FALSE, name = paste("Chr", chr)
      )
    }

    if (nrow(sig_data) > 0) {
      p <- p %>% add_trace(
        data = sig_data, x = ~Position, y = ~LOD,
        type = "scatter", mode = "markers",
        marker = list(size = sig_point_size, color = sig_color,
                      symbol = "diamond", line = list(width = 1, color = "#660000")),
        text = ~paste0("<b>", Marker, "</b><br>Chr: ", Chromosome,
                       "<br>Pos: ", round(Position, 1), " cM<br>LOD: ", round(LOD, 2),
                       "<br>PVE: ", round(PVE_Marker_Percent, 1), "%"),
        hoverinfo = "text", showlegend = FALSE, name = paste("Chr", chr, "Significant")
      )
    }

    p <- p %>% add_segments(
      x = min(chr_data$Position, na.rm=TRUE), xend = max(chr_data$Position, na.rm=TRUE),
      y = lod_threshold, yend = lod_threshold,
      line = list(color = sig_color, dash = threshold_line_style, width = threshold_line_width),
      showlegend = FALSE, hoverinfo = "none"
    )

    p <- p %>% layout(
      xaxis = list(title = paste("Chr", chr, "Position (cM)"), gridcolor = "#EEEEEE", zeroline = FALSE),
      yaxis = list(title = if(chr == chroms[1]) "LOD Score" else "", gridcolor = "#EEEEEE", zeroline = FALSE),
      plot_bgcolor = "white", paper_bgcolor = "white"
    )

    subplot_list[[chr]] <- p
  }

  n_chroms <- length(chroms)
  n_cols   <- ceiling(sqrt(n_chroms))
  n_rows   <- ceiling(n_chroms / n_cols)

  combined_plot <- subplot(
    subplot_list, nrows = n_rows, shareY = TRUE,
    titleX = TRUE, titleY = TRUE, margin = 0.05
  )

  combined_plot %>% layout(
    title = list(
      text = paste0(title, " | ", sum(pd$sig, na.rm=TRUE), " significant QTLs"),
      font = list(size = 14, color = "#111111")
    ),
    showlegend = FALSE,
    font = list(color = "#111111"),
    hoverlabel = list(bgcolor = "#FFFFFF", font = list(color = "#111111", size = 11))
  )
}

# ============================================================
# COMMON QTL FINDER
# ============================================================
find_common_qtl <- function(compare_res, lod_threshold, min_models=2, window_cM=5) {
  sig_list <- list()
  for (mth in names(compare_res)) {
    r <- compare_res[[mth]]
    if (is.null(r)||"Message"%in%colnames(r)) next
    sig <- r[!is.na(r$LOD)&r$LOD>=lod_threshold,
             c("Marker","Chromosome","Position","LOD","PVE_Marker_Percent","Additive_Effect"),drop=FALSE]
    if (nrow(sig)>0){sig$Model<-mth;sig_list[[mth]]<-sig}
  }
  if (length(sig_list)==0) return(NULL)
  all_sig <- do.call(rbind,sig_list)
  all_sig$Chromosome <- as.character(all_sig$Chromosome)
  marker_counts <- table(all_sig$Marker)
  all_sig$N_Models_Detected <- as.integer(marker_counts[all_sig$Marker])
  common <- all_sig[all_sig$N_Models_Detected>=min_models,]
  if (nrow(common)==0) return(NULL)
  markers_uniq <- unique(common$Marker)
  summary_rows <- lapply(markers_uniq, function(mk) {
    rows <- common[common$Marker==mk,]
    base <- rows[which.max(rows$LOD), c("Marker","Chromosome","Position","N_Models_Detected")]
    base$Models_List  <- paste(sort(unique(rows$Model)), collapse=", ")
    base$Mean_LOD     <- round(mean(rows$LOD, na.rm=TRUE), 3)
    base$Max_LOD      <- round(max(rows$LOD, na.rm=TRUE), 3)
    base$Mean_PVE     <- round(mean(rows$PVE_Marker_Percent, na.rm=TRUE), 3)
    base$Mean_Add_Eff <- round(mean(rows$Additive_Effect, na.rm=TRUE), 3)
    for (mth in names(compare_res)) {
      r_mth <- rows[rows$Model==mth,]
      base[[paste0("LOD_",mth)]] <- if(nrow(r_mth)>0) round(r_mth$LOD[1],3) else NA
    }
    base
  })
  result <- do.call(rbind, summary_rows)
  result[order(-result$N_Models_Detected, -result$Max_LOD),]
}
