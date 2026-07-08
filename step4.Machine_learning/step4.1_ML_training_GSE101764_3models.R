############################################################
## ML training using final 110 DMRs in GSE101764
## Models:
## 1. Elastic-net / LASSO logistic regression
## 2. Random Forest
## 3. Radial SVM
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
options(scipen = 20)
set.seed(123)

############################################################
## Step 0. Load required packages
############################################################

required_pkgs <- c(
  "data.table",
  "caret",
  "glmnet",
  "ranger",
  "kernlab",
  "pROC",
  "minfi",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19"
)

missing_pkgs <- required_pkgs[
  !sapply(required_pkgs, requireNamespace, quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "Please install these packages first: ",
      paste(missing_pkgs, collapse = ", ")
    )
  )
}

library(data.table)
library(caret)
library(glmnet)
library(ranger)
library(kernlab)
library(pROC)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

############################################################
## Step 1. Set input files
############################################################

## Final 110 DMRs after whole blood filtering
dmr_file <- "final_DMRs_after_GSE196696_whole_blood_filter.csv"

## GSE101764 ChAMP-filtered object
## This file should contain myLoad$beta and myLoad$pd
gse101764_file <- "step0.1-GSE101764_output.Rdata"

## Output directory
outdir <- "GSE101764_ML_training_110DMRs"

if (!dir.exists(outdir)) {
  dir.create(outdir)
}

if (!file.exists(dmr_file)) {
  stop(paste0("DMR file not found: ", dmr_file))
}

if (!file.exists(gse101764_file)) {
  stop(paste0("GSE101764 Rdata file not found: ", gse101764_file))
}

############################################################
## Step 2. Load final DMR table
############################################################
library(data.table)
dmr.df <- fread(dmr_file, data.table = FALSE)

cat("\nFinal DMR table:\n")
print(dim(dmr.df))
print(head(dmr.df))
print(colnames(dmr.df))

############################################################
## Step 3. Helper functions for DMR coordinates
############################################################

match_col <- function(df, candidates) {
  
  idx <- match(tolower(candidates), tolower(colnames(df)))
  idx <- idx[!is.na(idx)]
  
  if (length(idx) == 0) {
    return(NA_character_)
  } else {
    return(colnames(df)[idx[1]])
  }
}

standardize_chr <- function(x) {
  
  x <- as.character(x)
  x <- gsub("^CHR", "chr", x, ignore.case = TRUE)
  x <- ifelse(grepl("^chr", x, ignore.case = TRUE), x, paste0("chr", x))
  
  return(x)
}

parse_coord <- function(coord) {
  
  coord <- as.character(coord)
  coord <- gsub(",", "", coord)
  
  ok <- grepl("^chr[^:]+:[0-9]+-[0-9]+$", coord)
  
  chr <- rep(NA_character_, length(coord))
  start <- rep(NA_integer_, length(coord))
  end <- rep(NA_integer_, length(coord))
  
  chr[ok] <- sub(":.*", "", coord[ok])
  start[ok] <- as.integer(sub(".*:([0-9]+)-[0-9]+$", "\\1", coord[ok]))
  end[ok] <- as.integer(sub(".*:[0-9]+-([0-9]+)$", "\\1", coord[ok]))
  
  data.frame(
    chr = standardize_chr(chr),
    start = start,
    end = end,
    stringsAsFactors = FALSE
  )
}

find_coord_col <- function(df) {
  
  candidate_cols <- colnames(df)[
    grepl("coord|region|position|location|dmr", colnames(df), ignore.case = TRUE)
  ]
  
  if (length(candidate_cols) == 0) {
    candidate_cols <- colnames(df)
  }
  
  for (cn in candidate_cols) {
    
    vals <- as.character(df[[cn]])
    vals <- vals[!is.na(vals)]
    vals <- vals[seq_len(min(length(vals), 100))]
    
    if (any(grepl("^chr[^:]+:[0-9,]+-[0-9,]+$", vals))) {
      return(cn)
    }
  }
  
  return(NA_character_)
}

standardize_dmr_table <- function(dmr.df) {
  
  chr_col <- match_col(
    dmr.df,
    c("chr", "chrom", "chromosome", "seqnames")
  )
  
  start_col <- match_col(
    dmr.df,
    c("start", "start_position", "pos_start", "dmr_start")
  )
  
  end_col <- match_col(
    dmr.df,
    c("end", "end_position", "pos_end", "dmr_end")
  )
  
  id_col <- match_col(
    dmr.df,
    c("DMR_id", "dmr_id", "region", "coord", "coordinates")
  )
  
  coord_col <- find_coord_col(dmr.df)
  coord_parsed <- NULL
  
  if (!is.na(coord_col)) {
    coord_vals <- as.character(dmr.df[[coord_col]])
    
    if (any(grepl("^chr[^:]+:[0-9,]+-[0-9,]+$", coord_vals))) {
      coord_parsed <- parse_coord(coord_vals)
    }
  }
  
  if (!is.na(chr_col)) {
    chr_out <- standardize_chr(dmr.df[[chr_col]])
  } else if (!is.null(coord_parsed)) {
    chr_out <- coord_parsed$chr
  } else {
    chr_out <- rep(NA_character_, nrow(dmr.df))
  }
  
  if (!is.na(start_col)) {
    start_out <- as.integer(gsub(",", "", dmr.df[[start_col]]))
  } else if (!is.null(coord_parsed)) {
    start_out <- coord_parsed$start
  } else {
    start_out <- rep(NA_integer_, nrow(dmr.df))
  }
  
  if (!is.na(end_col)) {
    end_out <- as.integer(gsub(",", "", dmr.df[[end_col]]))
  } else if (!is.null(coord_parsed)) {
    end_out <- coord_parsed$end
  } else {
    end_out <- rep(NA_integer_, nrow(dmr.df))
  }
  
  if (!is.na(id_col)) {
    dmr_id <- as.character(dmr.df[[id_col]])
  } else {
    dmr_id <- paste0(chr_out, ":", start_out, "-", end_out)
  }
  
  dmr_id <- make.unique(dmr_id)
  
  out <- data.frame(
    DMR_id = dmr_id,
    chr = chr_out,
    start = start_out,
    end = end_out,
    stringsAsFactors = FALSE
  )
  
  if (any(is.na(out$chr)) || any(is.na(out$start)) || any(is.na(out$end))) {
    stop("Some DMRs do not have valid chr/start/end coordinates. Please check the DMR file columns.")
  }
  
  return(out)
}

dmr.tbl <- standardize_dmr_table(dmr.df)

cat("\nStandardized DMR table:\n")
print(dim(dmr.tbl))
print(head(dmr.tbl))

############################################################
## Step 4. Load GSE101764 beta matrix and phenotype table
############################################################

env <- new.env()
loaded_objects <- load(gse101764_file, envir = env)

if (!"myLoad" %in% loaded_objects) {
  stop("Object myLoad was not found in the GSE101764 Rdata file.")
}

myLoad <- env$myLoad

if (!"beta" %in% names(myLoad)) {
  stop("myLoad$beta was not found.")
}

if (!"pd" %in% names(myLoad)) {
  stop("myLoad$pd was not found.")
}

beta <- as.matrix(myLoad$beta)
mode(beta) <- "numeric"

pd <- as.data.frame(myLoad$pd)

cat("\nGSE101764 beta and pd:\n")
print(dim(beta))
print(dim(pd))
print(head(pd))

############################################################
## Step 5. Align phenotype table with beta matrix
############################################################

if ("ID" %in% colnames(pd) && setequal(as.character(pd$ID), colnames(beta))) {
  rownames(pd) <- as.character(pd$ID)
}

if (setequal(rownames(pd), colnames(beta))) {
  pd <- pd[colnames(beta), , drop = FALSE]
}

if (!identical(rownames(pd), colnames(beta))) {
  stop("Sample alignment failed between pd and beta.")
}

if (!"group" %in% colnames(pd)) {
  stop("Column 'group' was not found in pd.")
}

pd$group <- factor(pd$group, levels = c("Normal", "Tumor"))

cat("\nGroup distribution:\n")
print(table(pd$group, useNA = "ifany"))

############################################################
## Step 6. Map DMRs to 450K CpGs
############################################################

anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

anno.df <- data.frame(
  CpG = rownames(anno),
  chr = standardize_chr(anno$chr),
  pos = as.integer(anno$pos),
  stringsAsFactors = FALSE
)

anno.dt <- as.data.table(anno.df)

build_dmr_cpg_list <- function(dmr.tbl, anno.dt) {
  
  cpg.list <- vector("list", nrow(dmr.tbl))
  names(cpg.list) <- dmr.tbl$DMR_id
  
  for (i in seq_len(nrow(dmr.tbl))) {
    
    chr_i <- dmr.tbl$chr[i]
    start_i <- dmr.tbl$start[i]
    end_i <- dmr.tbl$end[i]
    
    cpgs <- anno.dt[
      chr == chr_i & pos >= start_i & pos <= end_i,
      CpG
    ]
    
    cpg.list[[i]] <- unique(cpgs)
  }
  
  return(cpg.list)
}

dmr.cpg.list <- build_dmr_cpg_list(dmr.tbl, anno.dt)

dmr.tbl$nCpG_450K_annotation <- lengths(dmr.cpg.list)
dmr.tbl$nCpG_present_GSE101764 <- sapply(
  dmr.cpg.list,
  function(x) length(intersect(x, rownames(beta)))
)

cat("\nCpG coverage summary:\n")
print(summary(dmr.tbl$nCpG_present_GSE101764))

fwrite(
  dmr.tbl,
  file = file.path(outdir, "DMR_CpG_mapping_summary_GSE101764.csv")
)

############################################################
## Step 7. Build DMR-level methylation matrix
## rows = samples
## columns = DMRs
############################################################

dmr_beta_mat <- matrix(
  NA_real_,
  nrow = ncol(beta),
  ncol = nrow(dmr.tbl)
)

rownames(dmr_beta_mat) <- colnames(beta)
colnames(dmr_beta_mat) <- dmr.tbl$DMR_id

for (i in seq_len(nrow(dmr.tbl))) {
  
  cpgs_present <- intersect(dmr.cpg.list[[i]], rownames(beta))
  
  if (length(cpgs_present) > 0) {
    dmr_beta_mat[, i] <- colMeans(
      beta[cpgs_present, , drop = FALSE],
      na.rm = TRUE
    )
  }
}

dmr_beta_mat[is.nan(dmr_beta_mat)] <- NA_real_

## Remove DMRs with all missing values
keep_dmr <- colSums(!is.na(dmr_beta_mat)) > 0

dmr_beta_mat <- dmr_beta_mat[, keep_dmr, drop = FALSE]
dmr.tbl.ml <- dmr.tbl[keep_dmr, , drop = FALSE]

cat("\nDMR-level beta matrix:\n")
print(dim(dmr_beta_mat))

## Median imputation for occasional missing values
for (j in seq_len(ncol(dmr_beta_mat))) {
  
  if (any(is.na(dmr_beta_mat[, j]))) {
    med_j <- median(dmr_beta_mat[, j], na.rm = TRUE)
    dmr_beta_mat[is.na(dmr_beta_mat[, j]), j] <- med_j
  }
}

## Make ML-safe column names
feature_name_map <- data.frame(
  DMR_id = colnames(dmr_beta_mat),
  Feature = make.names(colnames(dmr_beta_mat), unique = TRUE),
  stringsAsFactors = FALSE
)

colnames(dmr_beta_mat) <- feature_name_map$Feature

fwrite(
  feature_name_map,
  file = file.path(outdir, "DMR_feature_name_map.csv")
)

fwrite(
  data.frame(
    Sample = rownames(dmr_beta_mat),
    Group = pd$group,
    dmr_beta_mat,
    check.names = FALSE
  ),
  file = file.path(outdir, "GSE101764_DMR_level_beta_matrix_110DMRs.csv")
)

############################################################
## Step 8. Prepare ML data
############################################################

ml_df <- data.frame(
  group = pd$group,
  dmr_beta_mat,
  check.names = FALSE
)

## For caret twoClassSummary:
## The first level is treated as the positive class.
ml_df$group <- factor(
  ifelse(ml_df$group == "Tumor", "Tumor", "Normal"),
  levels = c("Tumor", "Normal")
)

cat("\nML dataset:\n")
print(dim(ml_df))
print(table(ml_df$group))

## Remove near-zero variance features
nzv <- nearZeroVar(ml_df[, -1, drop = FALSE])

if (length(nzv) > 0) {
  cat("\nRemoving near-zero variance features:\n")
  print(length(nzv))
  ml_df <- ml_df[, -c(nzv + 1), drop = FALSE]
}

cat("\nML dataset after NZV filtering:\n")
print(dim(ml_df))

############################################################
## Step 9. Cross-validation setup
############################################################

set.seed(123)

ctrl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  verboseIter = FALSE
)

############################################################
## Step 10. Model 1: Elastic-net / LASSO logistic regression
############################################################

set.seed(123)

glmnet_grid <- expand.grid(
  alpha = c(0, 0.25, 0.5, 0.75, 1),
  lambda = 10^seq(-4, 1, length = 100)
)

model_glmnet <- train(
  group ~ .,
  data = ml_df,
  method = "glmnet",
  metric = "ROC",
  trControl = ctrl,
  tuneGrid = glmnet_grid,
  preProcess = c("center", "scale")
)

cat("\nBest glmnet tuning parameters:\n")
print(model_glmnet$bestTune)
print(max(model_glmnet$results$ROC, na.rm = TRUE))

############################################################
## Step 11. Model 2: Random Forest
############################################################

set.seed(123)

p <- ncol(ml_df) - 1

rf_grid <- expand.grid(
  mtry = unique(pmax(1, round(c(sqrt(p), p / 5, p / 3)))),
  splitrule = "gini",
  min.node.size = c(1, 5, 10)
)

model_rf <- train(
  group ~ .,
  data = ml_df,
  method = "ranger",
  metric = "ROC",
  trControl = ctrl,
  tuneGrid = rf_grid,
  importance = "impurity",
  num.trees = 1000
)

cat("\nBest Random Forest tuning parameters:\n")
print(model_rf$bestTune)
print(max(model_rf$results$ROC, na.rm = TRUE))

############################################################
## Step 12. Model 3: Radial SVM
############################################################

set.seed(123)

model_svm <- train(
  group ~ .,
  data = ml_df,
  method = "svmRadial",
  metric = "ROC",
  trControl = ctrl,
  tuneLength = 10,
  preProcess = c("center", "scale")
)

cat("\nBest SVM tuning parameters:\n")
print(model_svm$bestTune)
print(max(model_svm$results$ROC, na.rm = TRUE))

############################################################
## Step 13. Compare model performance
############################################################

get_best_performance <- function(model, model_name) {
  
  best <- model$bestTune
  res <- model$results
  
  for (cn in colnames(best)) {
    res <- res[res[[cn]] == best[[cn]], , drop = FALSE]
  }
  
  data.frame(
    Model = model_name,
    ROC = res$ROC[1],
    Sens = res$Sens[1],
    Spec = res$Spec[1],
    stringsAsFactors = FALSE
  )
}

model_performance <- rbind(
  get_best_performance(model_glmnet, "Elastic-net / LASSO logistic regression"),
  get_best_performance(model_rf, "Random Forest"),
  get_best_performance(model_svm, "Radial SVM")
)

cat("\nModel performance based on repeated 5-fold CV:\n")
print(model_performance)

fwrite(
  model_performance,
  file = file.path(outdir, "GSE101764_ML_model_performance_repeatedCV.csv")
)

############################################################
## Step 14. Confusion matrices from cross-validated predictions
############################################################

get_best_predictions <- function(model) {
  
  pred <- model$pred
  best <- model$bestTune
  
  for (cn in colnames(best)) {
    pred <- pred[pred[[cn]] == best[[cn]], , drop = FALSE]
  }
  
  return(pred)
}

pred_glmnet <- get_best_predictions(model_glmnet)
pred_rf <- get_best_predictions(model_rf)
pred_svm <- get_best_predictions(model_svm)

cm_glmnet <- confusionMatrix(
  pred_glmnet$pred,
  pred_glmnet$obs,
  positive = "Tumor"
)

cm_rf <- confusionMatrix(
  pred_rf$pred,
  pred_rf$obs,
  positive = "Tumor"
)

cm_svm <- confusionMatrix(
  pred_svm$pred,
  pred_svm$obs,
  positive = "Tumor"
)

cat("\nConfusion matrix: glmnet\n")
print(cm_glmnet)

cat("\nConfusion matrix: Random Forest\n")
print(cm_rf)

cat("\nConfusion matrix: SVM\n")
print(cm_svm)

capture.output(
  cm_glmnet,
  file = file.path(outdir, "confusion_matrix_glmnet.txt")
)

capture.output(
  cm_rf,
  file = file.path(outdir, "confusion_matrix_random_forest.txt")
)

capture.output(
  cm_svm,
  file = file.path(outdir, "confusion_matrix_svm_radial.txt")
)

############################################################
## Step 15. Extract important DMRs / selected markers
############################################################

## 15.1 Selected DMRs from glmnet
best_lambda <- model_glmnet$bestTune$lambda

coef_glmnet <- coef(model_glmnet$finalModel, s = best_lambda)
coef_glmnet_mat <- as.matrix(coef_glmnet)

glmnet_selected <- data.frame(
  Feature = rownames(coef_glmnet_mat),
  Coefficient = coef_glmnet_mat[, 1],
  stringsAsFactors = FALSE
)

glmnet_selected <- glmnet_selected[
  glmnet_selected$Feature != "(Intercept)" &
    glmnet_selected$Coefficient != 0,
]

glmnet_selected <- merge(
  glmnet_selected,
  feature_name_map,
  by = "Feature",
  all.x = TRUE
)

glmnet_selected <- glmnet_selected[
  order(abs(glmnet_selected$Coefficient), decreasing = TRUE),
]

cat("\nNumber of DMRs selected by glmnet:\n")
print(nrow(glmnet_selected))

print(head(glmnet_selected, 20))

fwrite(
  glmnet_selected,
  file = file.path(outdir, "glmnet_selected_DMRs.csv")
)

## 15.2 Random Forest variable importance
rf_importance <- varImp(model_rf)$importance

rf_importance$Feature <- rownames(rf_importance)

rf_importance <- merge(
  rf_importance,
  feature_name_map,
  by = "Feature",
  all.x = TRUE
)

rf_importance <- rf_importance[
  order(rf_importance$Overall, decreasing = TRUE),
]

fwrite(
  rf_importance,
  file = file.path(outdir, "random_forest_variable_importance.csv")
)

## 15.3 SVM variable importance
svm_importance <- varImp(model_svm)$importance

svm_importance$Feature <- rownames(svm_importance)

## SVM varImp may return class-specific columns such as Tumor and Normal
## instead of an Overall column.
## Create an Overall score from the maximum class-specific importance.
if (!"Overall" %in% colnames(svm_importance)) {
  
  importance_cols <- setdiff(colnames(svm_importance), "Feature")
  
  svm_importance$Overall <- apply(
    svm_importance[, importance_cols, drop = FALSE],
    1,
    max,
    na.rm = TRUE
  )
}

svm_importance <- merge(
  svm_importance,
  feature_name_map,
  by = "Feature",
  all.x = TRUE
)

svm_importance <- svm_importance[
  order(svm_importance$Overall, decreasing = TRUE),
]

fwrite(
  svm_importance,
  file = file.path(outdir, "svm_variable_importance.csv")
)

head(svm_importance, 20)
############################################################
## Step 16. Save all ML objects for TCGA external testing
############################################################

save(
  dmr.tbl,
  dmr.tbl.ml,
  dmr.cpg.list,
  feature_name_map,
  ml_df,
  model_glmnet,
  model_rf,
  model_svm,
  model_performance,
  glmnet_selected,
  rf_importance,
  svm_importance,
  file = file.path(outdir, "GSE101764_ML_training_3models_110DMRs.Rdata")
)

## Step 17. Final output summary
cat("\nModel performance:\n")
print(model_performance)
