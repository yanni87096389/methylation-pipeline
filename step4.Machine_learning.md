---
title: "step4.Machine_learning"
author: "YC"
date: "2026-06-28"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## data preprocessing

### GSE101764
```{r }
###############################################
## GSE101764 450k methylation processed matrix pipeline
## Goal:1.get standard beta value matrix
##      2.Build phenotype table pD
##      3. Run ChAMP filtering.
############################################################

############################################################
## Module 0. Clear environment and load required packages
############################################################
rm(list = ls())
options(timeout = 100000) 
options(scipen = 20)
options(stringsAsFactors = F)
require(GEOquery)
require(Biobase)
library("impute")
library(ChAMP)

############################################################
## Module 1. Read GEO metadata
#Ensure that the methylation beta matrix and the phenotype information strictly correspond to each other
## Download/read GEO series matrix metadata.
eset <- getGEO("GSE101764",destdir = './',AnnotGPL = T,getGPL = F)
beta.m <- exprs(eset[[1]])
class(beta.m)
dim(beta.m)
############################################################
## Module 2. Build phenotype table pD

## Build a clean phenotype table.
## Modify this part when applying the pipeline to another dataset
pD.all <- pData(eset[[1]])
pD <- pD.all[,c(13,2)]
names(pD) <- c("group", "ID")
pD$group <- ifelse(grepl("mucosa",pD$group,ignore.case = T),
                   "Normal","Tumor")

#⭐Check if the extracted groupings are correct
check_df <- data.frame(Title = pD.all$characteristics_ch1.3, Group = pD$group)
head(check_df, 10)
tail(check_df,10)
table(pD$group)

############################################################
## Module 3. Align pD with beta.m
############################################################
## Check whether beta.m and pD contain the same sample set.
## setequal() ignores order.
identical(rownames(pD),colnames(beta.m))##This must be TRUE before champ.filter().
setequal(colnames(beta.m), pD$ID)
# Confirm that there are no duplicated sample names
sum(duplicated(colnames(beta.m)))
sum(duplicated(pD$ID))

# Set group order
pD$group <- factor(
  pD$group,
  levels = c("Normal", "Tumor"))

# Sort pD by group
sample_order <- order(pD$group)

pD <- pD[sample_order, , drop = FALSE]

# Reorder beta matrix columns according to pD
beta.m <- beta.m[, rownames(pD), drop = FALSE]

## Final strict alignment check.
identical(rownames(pD),colnames(beta.m))
dim(beta.m)
dim(pD)

############################################################
## Module 4. Run ChAMP filtering
############################################################
## Run ChAMP filtering using the cleaned beta matrix and aligned pD.
myLoad <- champ.filter(beta = beta.m, pd = pD, arraytype = "450k")
## Save the ChAMP object for downstream analysis.
save(myLoad, file = "step0.1-GSE101764_output.Rdata")
#myLoad <- load("step0.1-GSE101764_output.Rdata")
```

### TCGA-COAD 450k 
```{r }
###############################################
## TCGA-COAD 450k methylation processed matrix pipeline
## Goal:1.get standard beta value matrix
##      2.Build phenotype table pD
##      3.Run ChAMP filtering.
############################################################

############################################################
## Module 0. Clear environment and load required packages
############################################################

rm(list = ls())
options(timeout = 100000)
options(scipen = 20)
options(stringsAsFactors = FALSE)

require(data.table)
require(ChAMP)
require(stringr)

############################################################
## Module 1. Read TCGA methylation matrix
############################################################
## Input file:
##   TCGA-COAD.methylation450.tsv.gz
##
## Matrix format:
##   rows    = CpG probes
##   columns = TCGA samples
##   values  = beta values

beta.m <- fread(
  "TCGA-COAD.methylation450.tsv.gz",
  data.table = FALSE,
  check.names = FALSE
)

## Check matrix format
beta.m[1:5, 1:5]
class(beta.m)
dim(beta.m)
colnames(beta.m)[1:10]

############################################################
## Module 2. Build standard beta value matrix
############################################################

## Set probe IDs as row names
rownames(beta.m) <- beta.m[, 1]

## Remove the first column containing probe IDs
beta.m <- beta.m[, -1]

## Convert beta data frame to numeric matrix
beta.m <- as.matrix(beta.m)
storage.mode(beta.m) <- "numeric"

## Basic checks
class(beta.m)
dim(beta.m)
head(rownames(beta.m))
head(colnames(beta.m))

sum(is.na(rownames(beta.m)))
sum(is.na(colnames(beta.m)))

anyDuplicated(rownames(beta.m))
anyDuplicated(colnames(beta.m))

############################################################
## Module 3. Build phenotype table pD
############################################################

## Extract sample names from beta matrix
sample_names <- colnames(beta.m)

## Extract TCGA sample type code
## TCGA barcode example:
##   TCGA-XX-XXXX-01A-...
## The first two digits of the 4th field indicate sample type.
sample_type_code <- sapply(
  strsplit(sample_names, "-"),
  function(x) substr(x[4], 1, 2)
)

## Check sample type code distribution
table(sample_type_code)

## Define group:
##   01-09 = Tumor
##   10-19 = Normal
group <- ifelse(
  sample_type_code %in% c("01", "02", "03", "04", "05", "06", "07", "08", "09"),
  "Tumor",
  ifelse(
    sample_type_code %in% c("10", "11", "12", "13", "14", "15", "16", "17", "18", "19"),
    "Normal",
    "Other"
  )
)

## Build phenotype table
pD <- data.frame(
  ID = sample_names,
  sample_type_code = sample_type_code,
  group = group,
  project = "TCGA-COAD",
  stringsAsFactors = FALSE
)

## Check group assignment
table(pD$group)

check_df <- data.frame(
  SampleID = pD$ID,
  SampleTypeCode = pD$sample_type_code,
  Group = pD$group,
  Project = pD$project
)

head(check_df, 10)
tail(check_df, 10)

############################################################
## Module 4. Keep only Tumor and Normal samples
############################################################

keep_samples <- pD$group %in% c("Tumor", "Normal")

pD <- pD[keep_samples, , drop = FALSE]
beta.m <- beta.m[, pD$ID, drop = FALSE]

## Set pD row names as sample IDs
rownames(pD) <- pD$ID

## Final strict alignment check
identical(colnames(beta.m), rownames(pD))
stopifnot(identical(colnames(beta.m), rownames(pD)))

## Set group order
pD$group <- factor(
  pD$group,
  levels = c("Normal", "Tumor")
)

## Sort samples by group
sample_order <- order(pD$group)

pD <- pD[sample_order, , drop = FALSE]
beta.m <- beta.m[, rownames(pD), drop = FALSE]

## Final check before ChAMP filtering
stopifnot(identical(colnames(beta.m), rownames(pD)))

table(pD$group)
dim(beta.m)
dim(pD)

## Module 5. Run ChAMP filtering

## Run ChAMP filtering using the cleaned beta matrix and aligned pD.
myLoad <- champ.filter(
  beta = beta.m,
  pd = pD,
  arraytype = "450K")

## Module 6. Save the ChAMP object for downstream analysis

save(myLoad,file = "step0.2-TCGA-COAD_output.Rdata")

```

### TCGA-READ 450k

```{r }
###############################################
## TCGA-READ 450k methylation processed matrix pipeline
## Goal:1.get standard beta value matrix
##      2.Build phenotype table pD
##      3. Run ChAMP filtering.
############################################################

############################################################
## Module 0. Clear environment and load required packages
############################################################

rm(list = ls())
options(timeout = 100000)
options(scipen = 20)
options(stringsAsFactors = FALSE)

require(data.table)
require(ChAMP)
require(stringr)

############################################################
## Module 1. Read TCGA methylation matrix
############################################################
## Input file:
##   TCGA-READ.methylation450.tsv.gz
##
## Matrix format:
##   rows    = CpG probes
##   columns = TCGA samples
##   values  = beta values

beta.m <- fread(
  "TCGA-READ.methylation450.tsv.gz",
  data.table = FALSE,
  check.names = FALSE
)

## Check matrix format
beta.m[1:5, 1:5]
class(beta.m)
dim(beta.m)
colnames(beta.m)[1:10]

############################################################
## Module 2. Build standard beta value matrix
############################################################

## Set probe IDs as row names
rownames(beta.m) <- beta.m[, 1]

## Remove the first column containing probe IDs
beta.m <- beta.m[, -1]

## Convert beta data frame to numeric matrix
beta.m <- as.matrix(beta.m)
storage.mode(beta.m) <- "numeric"

## Basic checks
class(beta.m)
dim(beta.m)
head(rownames(beta.m))
head(colnames(beta.m))

sum(is.na(rownames(beta.m)))
sum(is.na(colnames(beta.m)))

anyDuplicated(rownames(beta.m))
anyDuplicated(colnames(beta.m))

############################################################
## Module 3. Build phenotype table pD
############################################################

## Extract sample names from beta matrix
sample_names <- colnames(beta.m)

## Extract TCGA sample type code
## TCGA barcode example:
##   TCGA-XX-XXXX-01A-...
## The first two digits of the 4th field indicate sample type.
sample_type_code <- sapply(
  strsplit(sample_names, "-"),
  function(x) substr(x[4], 1, 2)
)

## Check sample type code distribution
table(sample_type_code)

## Define group:
##   01-09 = Tumor
##   10-19 = Normal
group <- ifelse(
  sample_type_code %in% c("01", "02", "03", "04", "05", "06", "07", "08", "09"),
  "Tumor",
  ifelse(
    sample_type_code %in% c("10", "11", "12", "13", "14", "15", "16", "17", "18", "19"),
    "Normal",
    "Other"
  )
)

## Build phenotype table
pD <- data.frame(
  ID = sample_names,
  sample_type_code = sample_type_code,
  group = group,
  project = "TCGA-READ",
  stringsAsFactors = FALSE
)

## Check group assignment
table(pD$group)

check_df <- data.frame(
  SampleID = pD$ID,
  SampleTypeCode = pD$sample_type_code,
  Group = pD$group,
  Project = pD$project
)

head(check_df, 10)
tail(check_df, 10)

############################################################
## Module 4. Keep only Tumor and Normal samples
############################################################

keep_samples <- pD$group %in% c("Tumor", "Normal")

pD <- pD[keep_samples, , drop = FALSE]
beta.m <- beta.m[, pD$ID, drop = FALSE]

## Set pD row names as sample IDs
rownames(pD) <- pD$ID

## Final strict alignment check
identical(colnames(beta.m), rownames(pD))
stopifnot(identical(colnames(beta.m), rownames(pD)))

## Set group order
pD$group <- factor(
  pD$group,
  levels = c("Normal", "Tumor")
)

## Sort samples by group
sample_order <- order(pD$group)

pD <- pD[sample_order, , drop = FALSE]
beta.m <- beta.m[, rownames(pD), drop = FALSE]

## Final check before ChAMP filtering
stopifnot(identical(colnames(beta.m), rownames(pD)))

table(pD$group)
dim(beta.m)
dim(pD)

############################################################
## Module 5. Run ChAMP filtering
############################################################

## Run ChAMP filtering using the cleaned beta matrix and aligned pD.
myLoad <- champ.filter(
  beta = beta.m,
  pd = pD,
  arraytype = "450K"
)

############################################################
## Module 6. Save the ChAMP object for downstream analysis
############################################################

save(myLoad,file = "step0.3-TCGA-READ_output.Rdata")

```

## Machine learning training in GSE101764
We select 3 Models to do ML training:
1. Elastic-net / LASSO logistic regression
2. Random Forest
3. Radial SVM

```{r }
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

```


## external ML testing in TCGA

### TCGA full/balanced cohort preprocessing
```{r}
############################################################
## Create TCGA-COAD + TCGA-READ balanced cohort
## Input:
##   step0.2-TCGA-COAD_output.Rdata
##   step0.3-TCGA-READ_output.Rdata
##
## Output:
##   step0.TCGA_COAD_READ_full_output.Rdata
##   step0.TCGA_COAD_READ_balanced_60Tumor_40Normal_output.Rdata
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
options(scipen = 20)

library(data.table)

set.seed(123)

############################################################
## Module 1. Set input files
############################################################

coad_file <- "step0.2-TCGA-COAD_output.Rdata"
read_file <- "step0.3-TCGA-READ_output.Rdata"

target_tumor <- 60
target_normal <- 40

if (!file.exists(coad_file)) {
  stop(paste0("Cannot find file: ", coad_file))
}

if (!file.exists(read_file)) {
  stop(paste0("Cannot find file: ", read_file))
}

############################################################
## Module 2. Helper function to load myLoad object
############################################################

load_tcga_myload <- function(file, project_name) {
  
  env <- new.env()
  loaded_objects <- load(file, envir = env)
  
  if (!"myLoad" %in% loaded_objects) {
    stop(paste0(
      "Object myLoad was not found in ",
      file,
      ". Objects found: ",
      paste(loaded_objects, collapse = ", ")
    ))
  }
  
  obj <- env$myLoad
  
  if (!"beta" %in% names(obj)) {
    stop(paste0("myLoad$beta was not found in ", file))
  }
  
  if (!"pd" %in% names(obj)) {
    stop(paste0("myLoad$pd was not found in ", file))
  }
  
  beta <- as.matrix(obj$beta)
  storage.mode(beta) <- "numeric"
  
  pd <- as.data.frame(obj$pd)
  
  ## Make sure phenotype table has sample ID column
  if (!"ID" %in% colnames(pd)) {
    pd$ID <- rownames(pd)
  }
  
  ## Align pd to beta matrix
  if (!setequal(as.character(pd$ID), colnames(beta))) {
    
    if (setequal(rownames(pd), colnames(beta))) {
      pd$ID <- rownames(pd)
    } else {
      stop(paste0("Sample IDs in beta and pd do not match in ", file))
    }
  }
  
  pd <- pd[match(colnames(beta), as.character(pd$ID)), , drop = FALSE]
  rownames(pd) <- pd$ID
  
  stopifnot(identical(colnames(beta), rownames(pd)))
  
  ## Make sure group column exists
  if (!"group" %in% colnames(pd)) {
    
    possible_group_col <- grep(
      "group|sample_group|samplegroup|phenotype|status",
      colnames(pd),
      ignore.case = TRUE,
      value = TRUE
    )
    
    if (length(possible_group_col) >= 1) {
      colnames(pd)[colnames(pd) == possible_group_col[1]] <- "group"
    } else {
      stop(paste0("Cannot find group column in ", file))
    }
  }
  
  ## Make sure project column exists
  if (!"project" %in% colnames(pd)) {
    pd$project <- project_name
  }
  
  pd$project <- project_name
  
  ## Standardise group labels
  pd$group <- as.character(pd$group)
  pd$group <- ifelse(
    grepl("tumor|tumour|cancer|carcinoma", pd$group, ignore.case = TRUE),
    "Tumor",
    ifelse(
      grepl("normal|control|adjacent|solid tissue normal", pd$group, ignore.case = TRUE),
      "Normal",
      pd$group
    )
  )
  
  pd$group <- factor(pd$group, levels = c("Normal", "Tumor"))
  
  cat("\nLoaded:", project_name, "\n")
  cat("Beta dimension:\n")
  print(dim(beta))
  cat("pd dimension:\n")
  print(dim(pd))
  cat("Group distribution:\n")
  print(table(pd$group, useNA = "ifany"))
  
  return(list(beta = beta, pd = pd))
}

############################################################
## Module 3. Load COAD and READ
############################################################

tcga_coad <- load_tcga_myload(
  file = coad_file,
  project_name = "TCGA-COAD"
)

tcga_read <- load_tcga_myload(
  file = read_file,
  project_name = "TCGA-READ"
)

############################################################
## Module 4. Merge COAD and READ by common CpGs
############################################################

common_cpgs <- intersect(
  rownames(tcga_coad$beta),
  rownames(tcga_read$beta)
)

cat("\nNumber of common CpGs between COAD and READ:\n")
print(length(common_cpgs))

beta_full <- cbind(
  tcga_coad$beta[common_cpgs, , drop = FALSE],
  tcga_read$beta[common_cpgs, , drop = FALSE]
)

pd_full <- rbind(
  tcga_coad$pd,
  tcga_read$pd
)

## Check duplicated sample IDs
if (anyDuplicated(colnames(beta_full)) > 0) {
  duplicated_samples <- colnames(beta_full)[duplicated(colnames(beta_full))]
  stop(paste0(
    "Duplicated sample IDs found after merging: ",
    paste(duplicated_samples, collapse = ", ")
  ))
}

## Align pd with beta
pd_full <- pd_full[match(colnames(beta_full), pd_full$ID), , drop = FALSE]
rownames(pd_full) <- pd_full$ID

stopifnot(identical(colnames(beta_full), rownames(pd_full)))

############################################################
## Module 5. Keep only Tumor and Normal samples
############################################################

keep_samples <- pd_full$group %in% c("Normal", "Tumor")

pd_full <- pd_full[keep_samples, , drop = FALSE]
beta_full <- beta_full[, rownames(pd_full), drop = FALSE]

pd_full$group <- factor(
  as.character(pd_full$group),
  levels = c("Normal", "Tumor")
)

## Sort Normal first, Tumor second
sample_order <- order(pd_full$group)

pd_full <- pd_full[sample_order, , drop = FALSE]
beta_full <- beta_full[, rownames(pd_full), drop = FALSE]

stopifnot(identical(colnames(beta_full), rownames(pd_full)))

cat("\nMerged full TCGA COAD + READ cohort:\n")
print(table(pd_full$project, pd_full$group))
print(table(pd_full$group))
print(dim(beta_full))
print(dim(pd_full))

############################################################
## Module 6. Save full merged TCGA cohort
############################################################

myLoad <- list(
  beta = beta_full,
  pd = pd_full
)

save(
  myLoad,
  file = "step0.TCGA_COAD_READ_full_output.Rdata"
)

fwrite(
  pd_full,
  file = "TCGA_COAD_READ_full_sample_information.csv"
)

cat("\nSaved full TCGA cohort:\n")
print("step0.TCGA_COAD_READ_full_output.Rdata")

############################################################
## Module 7. Create balanced cohort: 60 Tumor + 40 Normal
############################################################

tumor_ids <- rownames(pd_full)[pd_full$group == "Tumor"]
normal_ids <- rownames(pd_full)[pd_full$group == "Normal"]

cat("\nAvailable samples before balancing:\n")
cat("Tumor:", length(tumor_ids), "\n")
cat("Normal:", length(normal_ids), "\n")

if (length(tumor_ids) < target_tumor) {
  stop(paste0(
    "Not enough Tumor samples. Available = ",
    length(tumor_ids),
    ", target = ",
    target_tumor
  ))
}

if (length(normal_ids) < target_normal) {
  stop(paste0(
    "Not enough Normal samples. Available = ",
    length(normal_ids),
    ", target = ",
    target_normal
  ))
}

selected_tumor_ids <- sample(
  tumor_ids,
  size = target_tumor,
  replace = FALSE
)

selected_normal_ids <- sample(
  normal_ids,
  size = target_normal,
  replace = FALSE
)

selected_ids <- c(selected_normal_ids, selected_tumor_ids)

pd_balanced <- pd_full[selected_ids, , drop = FALSE]
beta_balanced <- beta_full[, selected_ids, drop = FALSE]

pd_balanced$group <- factor(
  as.character(pd_balanced$group),
  levels = c("Normal", "Tumor")
)

## Sort Normal first, Tumor second
sample_order <- order(pd_balanced$group)

pd_balanced <- pd_balanced[sample_order, , drop = FALSE]
beta_balanced <- beta_balanced[, rownames(pd_balanced), drop = FALSE]

stopifnot(identical(colnames(beta_balanced), rownames(pd_balanced)))

cat("\nBalanced TCGA COAD + READ cohort:\n")
print(table(pd_balanced$project, pd_balanced$group))
print(table(pd_balanced$group))
print(dim(beta_balanced))
print(dim(pd_balanced))

############################################################
## Module 8. Save balanced TCGA cohort
############################################################

myLoad <- list(
  beta = beta_balanced,
  pd = pd_balanced
)

save(
  myLoad,
  file = "step0.TCGA_COAD_READ_balanced_60Tumor_40Normal_output.Rdata"
)

fwrite(
  pd_balanced,
  file = "TCGA_COAD_READ_balanced_60Tumor_40Normal_sample_information.csv"
)

fwrite(
  data.frame(
    Sample = rownames(pd_balanced),
    Project = pd_balanced$project,
    Group = pd_balanced$group,
    stringsAsFactors = FALSE
  ),
  file = "TCGA_COAD_READ_balanced_60Tumor_40Normal_selected_samples.csv"
)

cat("\nSaved balanced TCGA cohort:\n")
print("step0.TCGA_COAD_READ_balanced_60Tumor_40Normal_output.Rdata")

```

### TCGA full cohort external testing using 3 trained ML models
Current-directory version:
Input:
1. GSE101764_ML_training_3models_110DMRs.Rdata
2. step0.TCGA_COAD_READ_full_output.Rdata
Models:
1. Elastic-net logistic regression
2. Random Forest
3. Radial SVM
```{r}
rm(list = ls())
options(stringsAsFactors = FALSE)
options(scipen = 20)

set.seed(123)

############################################################
## Module 0. Load required packages
############################################################

library(data.table)
library(caret)
library(pROC)
library(ggplot2)

############################################################
## Module 1. Set input files in current directory
############################################################

training_model_file <- "GSE101764_ML_training_3models_110DMRs.Rdata"

tcga_full_file <- "step0.TCGA_COAD_READ_full_output.Rdata"

if (!file.exists(training_model_file)) {
  stop(paste0("Training model file not found: ", training_model_file))
}

if (!file.exists(tcga_full_file)) {
  stop(paste0("TCGA full cohort file not found: ", tcga_full_file))
}

############################################################
## Module 2. Load trained models and training information
############################################################

load(training_model_file)

required_objects <- c(
  "model_glmnet",
  "model_rf",
  "model_svm",
  "dmr.tbl",
  "dmr.cpg.list",
  "feature_name_map",
  "ml_df"
)

missing_objects <- required_objects[!sapply(required_objects, exists)]

if (length(missing_objects) > 0) {
  stop(paste0(
    "Missing objects in training Rdata: ",
    paste(missing_objects, collapse = ", ")
  ))
}

cat("\nTraining objects loaded successfully.\n")

cat("\nBest glmnet tuning parameters:\n")
print(model_glmnet$bestTune)

############################################################
## Module 3. Extract model-specific features
############################################################

get_model_features <- function(model, fallback_ml_df = ml_df) {
  
  if (!is.null(model$trainingData)) {
    features <- setdiff(colnames(model$trainingData), ".outcome")
  } else {
    features <- setdiff(colnames(fallback_ml_df), "group")
  }
  
  return(features)
}

features_glmnet <- get_model_features(model_glmnet)
features_rf <- get_model_features(model_rf)
features_svm <- get_model_features(model_svm)

all_training_features <- unique(c(
  features_glmnet,
  features_rf,
  features_svm
))

cat("\nNumber of model features:\n")
cat("Elastic-net:", length(features_glmnet), "\n")
cat("Random Forest:", length(features_rf), "\n")
cat("Radial SVM:", length(features_svm), "\n")
cat("All unique features:", length(all_training_features), "\n")

############################################################
## Module 4. Prepare training medians for missing-value imputation
############################################################

training_medians <- sapply(
  ml_df[, all_training_features, drop = FALSE],
  function(x) median(as.numeric(x), na.rm = TRUE)
)

training_medians[is.na(training_medians)] <- 0.5

############################################################
## Module 5. Load TCGA full myLoad object
############################################################

env <- new.env()
loaded_objects <- load(tcga_full_file, envir = env)

if (!"myLoad" %in% loaded_objects) {
  stop(paste0(
    "Object myLoad was not found in: ",
    tcga_full_file,
    ". Objects found: ",
    paste(loaded_objects, collapse = ", ")
  ))
}

myLoad_tcga <- env$myLoad

if (!"beta" %in% names(myLoad_tcga)) {
  stop("myLoad$beta was not found in TCGA full cohort file.")
}

if (!"pd" %in% names(myLoad_tcga)) {
  stop("myLoad$pd was not found in TCGA full cohort file.")
}

beta_tcga <- as.matrix(myLoad_tcga$beta)
storage.mode(beta_tcga) <- "numeric"

pd_tcga <- as.data.frame(myLoad_tcga$pd)

############################################################
## Module 6. Align TCGA phenotype table with beta matrix
############################################################

if ("ID" %in% colnames(pd_tcga) && setequal(as.character(pd_tcga$ID), colnames(beta_tcga))) {
  rownames(pd_tcga) <- as.character(pd_tcga$ID)
}

if (setequal(rownames(pd_tcga), colnames(beta_tcga))) {
  pd_tcga <- pd_tcga[colnames(beta_tcga), , drop = FALSE]
}

if (!identical(rownames(pd_tcga), colnames(beta_tcga))) {
  stop("Sample alignment failed between TCGA beta and pd.")
}

if (!"group" %in% colnames(pd_tcga)) {
  stop("Column group was not found in TCGA pd.")
}

pd_tcga$group <- factor(
  ifelse(as.character(pd_tcga$group) == "Tumor", "Tumor", "Normal"),
  levels = c("Tumor", "Normal")
)

cat("\nTCGA full cohort sample distribution:\n")
print(table(pd_tcga$group, useNA = "ifany"))

if ("project" %in% colnames(pd_tcga)) {
  print(table(pd_tcga$project, pd_tcga$group))
}

cat("\nTCGA beta dimension:\n")
print(dim(beta_tcga))

cat("\nTCGA pd dimension:\n")
print(dim(pd_tcga))

############################################################
## Module 7. Build TCGA DMR-level beta matrix
############################################################

dmr_beta_mat <- matrix(
  NA_real_,
  nrow = ncol(beta_tcga),
  ncol = nrow(dmr.tbl)
)

rownames(dmr_beta_mat) <- colnames(beta_tcga)
colnames(dmr_beta_mat) <- dmr.tbl$DMR_id

for (i in seq_len(nrow(dmr.tbl))) {
  
  cpgs_present <- intersect(dmr.cpg.list[[i]], rownames(beta_tcga))
  
  if (length(cpgs_present) > 0) {
    dmr_beta_mat[, i] <- colMeans(
      beta_tcga[cpgs_present, , drop = FALSE],
      na.rm = TRUE
    )
  }
}

dmr_beta_mat[is.nan(dmr_beta_mat)] <- NA_real_

############################################################
## Module 8. Check CpG coverage in TCGA full cohort
############################################################

cpg_coverage <- data.frame(
  DMR_id = dmr.tbl$DMR_id,
  nCpG_present_TCGA_full = sapply(
    dmr.cpg.list,
    function(x) length(intersect(x, rownames(beta_tcga)))
  ),
  stringsAsFactors = FALSE
)

fwrite(
  cpg_coverage,
  file = "TCGA_full_DMR_CpG_coverage.csv"
)

cat("\nCpG coverage summary in TCGA full cohort:\n")
print(summary(cpg_coverage$nCpG_present_TCGA_full))

cat("\nNumber of DMRs with zero CpG coverage in TCGA:\n")
print(sum(cpg_coverage$nCpG_present_TCGA_full == 0))

############################################################
## Module 9. Map DMR IDs to training feature names
############################################################

feature_names <- feature_name_map$Feature[
  match(colnames(dmr_beta_mat), feature_name_map$DMR_id)
]

if (any(is.na(feature_names))) {
  missing_dmr_ids <- colnames(dmr_beta_mat)[is.na(feature_names)]
  stop(paste0(
    "Some TCGA DMRs could not be mapped to training feature names: ",
    paste(head(missing_dmr_ids, 20), collapse = ", ")
  ))
}

colnames(dmr_beta_mat) <- feature_names

missing_features <- setdiff(all_training_features, colnames(dmr_beta_mat))

if (length(missing_features) > 0) {
  stop(paste0(
    "These training features are missing in TCGA matrix: ",
    paste(missing_features, collapse = ", ")
  ))
}

dmr_beta_mat <- dmr_beta_mat[, all_training_features, drop = FALSE]

############################################################
## Module 10. Impute missing values using training medians
############################################################

for (fn in all_training_features) {
  
  if (any(is.na(dmr_beta_mat[, fn]))) {
    dmr_beta_mat[is.na(dmr_beta_mat[, fn]), fn] <- training_medians[fn]
  }
}

if (any(is.na(dmr_beta_mat))) {
  stop("There are still missing values after imputation.")
}

############################################################
## Module 11. Build final TCGA full testing dataframe
############################################################

test_df_full <- data.frame(
  group = pd_tcga$group,
  dmr_beta_mat,
  check.names = FALSE
)

rownames(test_df_full) <- rownames(pd_tcga)

cat("\nTCGA full testing dataframe:\n")
print(dim(test_df_full))
print(table(test_df_full$group))

fwrite(
  data.frame(
    Sample = rownames(test_df_full),
    Group = test_df_full$group,
    dmr_beta_mat,
    check.names = FALSE
  ),
  file = "TCGA_full_DMR_level_beta_matrix_for_testing.csv"
)

############################################################
## Module 12. Helper function to evaluate one model
############################################################

evaluate_one_model <- function(model, model_name, test_df, model_features) {
  
  x_test <- test_df[, model_features, drop = FALSE]
  
  prob <- predict(
    model,
    newdata = x_test,
    type = "prob"
  )
  
  pred_class <- predict(
    model,
    newdata = x_test,
    type = "raw"
  )
  
  obs <- factor(
    as.character(test_df$group),
    levels = c("Tumor", "Normal")
  )
  
  pred_class <- factor(
    as.character(pred_class),
    levels = c("Tumor", "Normal")
  )
  
  roc_obj <- roc(
    response = obs,
    predictor = prob$Tumor,
    levels = c("Normal", "Tumor"),
    direction = "<",
    quiet = TRUE
  )
  
  auc_value <- as.numeric(auc(roc_obj))
  
  cm <- confusionMatrix(
    data = pred_class,
    reference = obs,
    positive = "Tumor"
  )
  
  performance <- data.frame(
    Model = model_name,
    AUC = auc_value,
    Accuracy = as.numeric(cm$overall["Accuracy"]),
    Sensitivity = as.numeric(cm$byClass["Sensitivity"]),
    Specificity = as.numeric(cm$byClass["Specificity"]),
    Balanced_Accuracy = as.numeric(cm$byClass["Balanced Accuracy"]),
    stringsAsFactors = FALSE
  )
  
  pred_df <- data.frame(
    Sample = rownames(test_df),
    Observed = obs,
    Predicted = pred_class,
    Prob_Tumor = prob$Tumor,
    Prob_Normal = prob$Normal,
    Model = model_name,
    stringsAsFactors = FALSE
  )
  
  return(list(
    performance = performance,
    predictions = pred_df,
    roc = roc_obj,
    confusion_matrix = cm
  ))
}

############################################################
## Module 13. Test three models in TCGA full cohort
############################################################

res_glmnet_full <- evaluate_one_model(
  model = model_glmnet,
  model_name = "Elastic-net",
  test_df = test_df_full,
  model_features = features_glmnet
)

res_rf_full <- evaluate_one_model(
  model = model_rf,
  model_name = "Random Forest",
  test_df = test_df_full,
  model_features = features_rf
)

res_svm_full <- evaluate_one_model(
  model = model_svm,
  model_name = "Radial SVM",
  test_df = test_df_full,
  model_features = features_svm
)

############################################################
## Module 14. Combine and save performance results
############################################################

performance_full <- rbind(
  res_glmnet_full$performance,
  res_rf_full$performance,
  res_svm_full$performance
)

predictions_full <- rbind(
  res_glmnet_full$predictions,
  res_rf_full$predictions,
  res_svm_full$predictions
)

cat("\nTCGA full model performance:\n")
print(performance_full)

fwrite(
  performance_full,
  file = "TCGA_full_model_performance.csv"
)

fwrite(
  predictions_full,
  file = "TCGA_full_model_predictions.csv"
)

capture.output(
  res_glmnet_full$confusion_matrix,
  file = "TCGA_full_confusion_matrix_elastic_net.txt"
)

capture.output(
  res_rf_full$confusion_matrix,
  file = "TCGA_full_confusion_matrix_random_forest.txt"
)

capture.output(
  res_svm_full$confusion_matrix,
  file = "TCGA_full_confusion_matrix_radial_svm.txt"
)

############################################################
## Module 15. Plot ROC curves for three models
############################################################

roc_to_df <- function(roc_obj, model_name, auc_value) {
  
  data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    Model = paste0(model_name, " (AUC = ", sprintf("%.3f", auc_value), ")"),
    stringsAsFactors = FALSE
  )
}

roc_df_full <- rbind(
  roc_to_df(
    res_glmnet_full$roc,
    "Elastic-net",
    res_glmnet_full$performance$AUC
  ),
  roc_to_df(
    res_rf_full$roc,
    "Random Forest",
    res_rf_full$performance$AUC
  ),
  roc_to_df(
    res_svm_full$roc,
    "Radial SVM",
    res_svm_full$performance$AUC
  )
)

p_roc_full <- ggplot(
  roc_df_full,
  aes(x = FPR, y = TPR, color = Model)
) +
  geom_line(linewidth = 1.2) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  coord_equal() +
  labs(
    title = "ROC curves in TCGA full cohort",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Model"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.65, 0.20),
    legend.background = element_rect(fill = "white", color = "black")
  )

p_roc_full

ggsave(
  filename = "TCGA_full_ROC_curves_3models.png",
  plot = p_roc_full,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "TCGA_full_ROC_curves_3models.pdf",
  plot = p_roc_full,
  width = 8,
  height = 6
)

############################################################
## Module 16. Plot model performance heatmap
############################################################

library(data.table)
library(ggplot2)

performance_heat <- performance_full[, c(
  "Model",
  "AUC",
  "Accuracy",
  "Sensitivity",
  "Specificity",
  "Balanced_Accuracy"
)]

performance_heat <- as.data.table(performance_heat)

performance_heat <- data.table::melt(
  performance_heat,
  id.vars = "Model",
  variable.name = "Metric",
  value.name = "Value"
)

performance_heat$Metric <- factor(
  performance_heat$Metric,
  levels = c("AUC", "Accuracy", "Sensitivity", "Specificity", "Balanced_Accuracy"),
  labels = c("AUC", "Accuracy", "Sensitivity", "Specificity", "Balanced accuracy")
)

performance_heat$Model <- factor(
  performance_heat$Model,
  levels = c("Elastic-net", "Random Forest", "Radial SVM")
)

p_heat_full <- ggplot(
  performance_heat,
  aes(x = Metric, y = Model, fill = Value)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = sprintf("%.3f", Value)),
    size = 5
  ) +
  scale_fill_gradient(
    low = "#fdc771",
    high = "#e14c48",
    limits = c(0.95, 1.00)
  ) +
  labs(
    title = "Model performance in TCGA full cohort",
    x = NULL,
    y = NULL,
    fill = "Performance"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    legend.position = "right"
  )

p_heat_full

ggsave(
  filename = "TCGA_full_model_performance_heatmap.png",
  plot = p_heat_full,
  width = 9,
  height = 4.5,
  dpi = 300
)

ggsave(
  filename = "TCGA_full_model_performance_heatmap.pdf",
  plot = p_heat_full,
  width = 9,
  height = 4.5
)
############################################################
## Module 17. Save all full TCGA testing objects
############################################################

save(
  test_df_full,
  performance_full,
  predictions_full,
  roc_df_full,
  p_roc_full,
  res_glmnet_full,
  res_rf_full,
  res_svm_full,
  file = "TCGA_full_testing_results_3models.Rdata"
)

############################################################
## Module 18. Final summary
############################################################

cat("\nFinal TCGA full model performance:\n")
print(performance_full)

```

### TCGA balanced cohort external testing using 3 trained ML models
Current-directory version:
60 Tumor vs 40 Normal
```{r}
rm(list = ls())
options(stringsAsFactors = FALSE)
options(scipen = 20)

set.seed(123)

############################################################
## Module 0. Load required packages
############################################################

library(data.table)
library(caret)
library(pROC)
library(ggplot2)

############################################################
## Module 1. Set input files in current directory
############################################################

training_model_file <- "GSE101764_ML_training_3models_110DMRs.Rdata"

tcga_balanced_file <- "step0.TCGA_COAD_READ_balanced_60Tumor_40Normal_output.Rdata"

if (!file.exists(training_model_file)) {
  stop(paste0("Training model file not found: ", training_model_file))
}

if (!file.exists(tcga_balanced_file)) {
  stop(paste0("TCGA balanced file not found: ", tcga_balanced_file))
}

############################################################
## Module 2. Load trained models and training information
############################################################

load(training_model_file)

required_objects <- c(
  "model_glmnet",
  "model_rf",
  "model_svm",
  "dmr.tbl",
  "dmr.cpg.list",
  "feature_name_map",
  "ml_df"
)

missing_objects <- required_objects[!sapply(required_objects, exists)]

if (length(missing_objects) > 0) {
  stop(paste0(
    "Missing objects in training Rdata: ",
    paste(missing_objects, collapse = ", ")
  ))
}

cat("\nTraining objects loaded successfully.\n")

############################################################
## Module 3. Extract model-specific features
############################################################

get_model_features <- function(model, fallback_ml_df = ml_df) {
  
  if (!is.null(model$trainingData)) {
    features <- setdiff(colnames(model$trainingData), ".outcome")
  } else {
    features <- setdiff(colnames(fallback_ml_df), "group")
  }
  
  return(features)
}

features_glmnet <- get_model_features(model_glmnet)
features_rf <- get_model_features(model_rf)
features_svm <- get_model_features(model_svm)

all_training_features <- unique(c(
  features_glmnet,
  features_rf,
  features_svm
))

cat("\nNumber of model features:\n")
cat("Elastic-net:", length(features_glmnet), "\n")
cat("Random Forest:", length(features_rf), "\n")
cat("Radial SVM:", length(features_svm), "\n")
cat("All unique features:", length(all_training_features), "\n")

############################################################
## Module 4. Prepare training medians for missing-value imputation
############################################################

training_medians <- sapply(
  ml_df[, all_training_features, drop = FALSE],
  function(x) median(as.numeric(x), na.rm = TRUE)
)

training_medians[is.na(training_medians)] <- 0.5

############################################################
## Module 5. Load TCGA balanced myLoad object
############################################################

env <- new.env()
loaded_objects <- load(tcga_balanced_file, envir = env)

if (!"myLoad" %in% loaded_objects) {
  stop(paste0(
    "Object myLoad was not found in: ",
    tcga_balanced_file,
    ". Objects found: ",
    paste(loaded_objects, collapse = ", ")
  ))
}

myLoad_tcga <- env$myLoad

if (!"beta" %in% names(myLoad_tcga)) {
  stop("myLoad$beta was not found in TCGA balanced file.")
}

if (!"pd" %in% names(myLoad_tcga)) {
  stop("myLoad$pd was not found in TCGA balanced file.")
}

beta_tcga <- as.matrix(myLoad_tcga$beta)
storage.mode(beta_tcga) <- "numeric"

pd_tcga <- as.data.frame(myLoad_tcga$pd)

############################################################
## Module 6. Align TCGA phenotype table with beta matrix
############################################################

if ("ID" %in% colnames(pd_tcga) && setequal(as.character(pd_tcga$ID), colnames(beta_tcga))) {
  rownames(pd_tcga) <- as.character(pd_tcga$ID)
}

if (setequal(rownames(pd_tcga), colnames(beta_tcga))) {
  pd_tcga <- pd_tcga[colnames(beta_tcga), , drop = FALSE]
}

if (!identical(rownames(pd_tcga), colnames(beta_tcga))) {
  stop("Sample alignment failed between TCGA beta and pd.")
}

if (!"group" %in% colnames(pd_tcga)) {
  stop("Column group was not found in TCGA pd.")
}

pd_tcga$group <- factor(
  ifelse(as.character(pd_tcga$group) == "Tumor", "Tumor", "Normal"),
  levels = c("Tumor", "Normal")
)

cat("\nTCGA balanced cohort sample distribution:\n")
print(table(pd_tcga$group, useNA = "ifany"))

if ("project" %in% colnames(pd_tcga)) {
  print(table(pd_tcga$project, pd_tcga$group))
}

cat("\nTCGA beta dimension:\n")
print(dim(beta_tcga))

cat("\nTCGA pd dimension:\n")
print(dim(pd_tcga))

############################################################
## Module 7. Build TCGA DMR-level beta matrix
############################################################

dmr_beta_mat <- matrix(
  NA_real_,
  nrow = ncol(beta_tcga),
  ncol = nrow(dmr.tbl)
)

rownames(dmr_beta_mat) <- colnames(beta_tcga)
colnames(dmr_beta_mat) <- dmr.tbl$DMR_id

for (i in seq_len(nrow(dmr.tbl))) {
  
  cpgs_present <- intersect(dmr.cpg.list[[i]], rownames(beta_tcga))
  
  if (length(cpgs_present) > 0) {
    dmr_beta_mat[, i] <- colMeans(
      beta_tcga[cpgs_present, , drop = FALSE],
      na.rm = TRUE
    )
  }
}

dmr_beta_mat[is.nan(dmr_beta_mat)] <- NA_real_

############################################################
## Module 8. Check CpG coverage in TCGA balanced cohort
############################################################

cpg_coverage <- data.frame(
  DMR_id = dmr.tbl$DMR_id,
  nCpG_present_TCGA_balanced = sapply(
    dmr.cpg.list,
    function(x) length(intersect(x, rownames(beta_tcga)))
  ),
  stringsAsFactors = FALSE
)

fwrite(
  cpg_coverage,
  file = "TCGA_balanced_DMR_CpG_coverage.csv"
)

cat("\nCpG coverage summary in TCGA balanced cohort:\n")
print(summary(cpg_coverage$nCpG_present_TCGA_balanced))

############################################################
## Module 9. Map DMR IDs to training feature names
############################################################

feature_names <- feature_name_map$Feature[
  match(colnames(dmr_beta_mat), feature_name_map$DMR_id)
]

if (any(is.na(feature_names))) {
  missing_dmr_ids <- colnames(dmr_beta_mat)[is.na(feature_names)]
  stop(paste0(
    "Some TCGA DMRs could not be mapped to training feature names: ",
    paste(head(missing_dmr_ids, 20), collapse = ", ")
  ))
}

colnames(dmr_beta_mat) <- feature_names

missing_features <- setdiff(all_training_features, colnames(dmr_beta_mat))

if (length(missing_features) > 0) {
  stop(paste0(
    "These training features are missing in TCGA matrix: ",
    paste(missing_features, collapse = ", ")
  ))
}

dmr_beta_mat <- dmr_beta_mat[, all_training_features, drop = FALSE]

############################################################
## Module 10. Impute missing values using training medians
############################################################

for (fn in all_training_features) {
  
  if (any(is.na(dmr_beta_mat[, fn]))) {
    dmr_beta_mat[is.na(dmr_beta_mat[, fn]), fn] <- training_medians[fn]
  }
}

if (any(is.na(dmr_beta_mat))) {
  stop("There are still missing values after imputation.")
}

############################################################
## Module 11. Build final TCGA balanced testing dataframe
############################################################

test_df_balanced <- data.frame(
  group = pd_tcga$group,
  dmr_beta_mat,
  check.names = FALSE
)

rownames(test_df_balanced) <- rownames(pd_tcga)

cat("\nTCGA balanced testing dataframe:\n")
print(dim(test_df_balanced))
print(table(test_df_balanced$group))

fwrite(
  data.frame(
    Sample = rownames(test_df_balanced),
    Group = test_df_balanced$group,
    dmr_beta_mat,
    check.names = FALSE
  ),
  file = "TCGA_balanced_DMR_level_beta_matrix_for_testing.csv"
)

############################################################
## Module 12. Helper function to evaluate one model
############################################################

evaluate_one_model <- function(model, model_name, test_df, model_features) {
  
  x_test <- test_df[, model_features, drop = FALSE]
  
  prob <- predict(
    model,
    newdata = x_test,
    type = "prob"
  )
  
  pred_class <- predict(
    model,
    newdata = x_test,
    type = "raw"
  )
  
  obs <- factor(
    as.character(test_df$group),
    levels = c("Tumor", "Normal")
  )
  
  pred_class <- factor(
    as.character(pred_class),
    levels = c("Tumor", "Normal")
  )
  
  roc_obj <- roc(
    response = obs,
    predictor = prob$Tumor,
    levels = c("Normal", "Tumor"),
    direction = "<",
    quiet = TRUE
  )
  
  auc_value <- as.numeric(auc(roc_obj))
  
  cm <- confusionMatrix(
    data = pred_class,
    reference = obs,
    positive = "Tumor"
  )
  
  performance <- data.frame(
    Model = model_name,
    AUC = auc_value,
    Accuracy = as.numeric(cm$overall["Accuracy"]),
    Sensitivity = as.numeric(cm$byClass["Sensitivity"]),
    Specificity = as.numeric(cm$byClass["Specificity"]),
    Balanced_Accuracy = as.numeric(cm$byClass["Balanced Accuracy"]),
    stringsAsFactors = FALSE
  )
  
  pred_df <- data.frame(
    Sample = rownames(test_df),
    Observed = obs,
    Predicted = pred_class,
    Prob_Tumor = prob$Tumor,
    Prob_Normal = prob$Normal,
    Model = model_name,
    stringsAsFactors = FALSE
  )
  
  return(list(
    performance = performance,
    predictions = pred_df,
    roc = roc_obj,
    confusion_matrix = cm
  ))
}

############################################################
## Module 13. Test three models in TCGA balanced cohort
############################################################

res_glmnet_balanced <- evaluate_one_model(
  model = model_glmnet,
  model_name = "Elastic-net",
  test_df = test_df_balanced,
  model_features = features_glmnet
)

res_rf_balanced <- evaluate_one_model(
  model = model_rf,
  model_name = "Random Forest",
  test_df = test_df_balanced,
  model_features = features_rf
)

res_svm_balanced <- evaluate_one_model(
  model = model_svm,
  model_name = "Radial SVM",
  test_df = test_df_balanced,
  model_features = features_svm
)

############################################################
## Module 14. Combine and save performance results
############################################################

performance_balanced <- rbind(
  res_glmnet_balanced$performance,
  res_rf_balanced$performance,
  res_svm_balanced$performance
)

predictions_balanced <- rbind(
  res_glmnet_balanced$predictions,
  res_rf_balanced$predictions,
  res_svm_balanced$predictions
)

cat("\nTCGA balanced model performance:\n")
print(performance_balanced)

fwrite(
  performance_balanced,
  file = "TCGA_balanced_model_performance.csv"
)

fwrite(
  predictions_balanced,
  file = "TCGA_balanced_model_predictions.csv"
)

capture.output(
  res_glmnet_balanced$confusion_matrix,
  file = "TCGA_balanced_confusion_matrix_elastic_net.txt"
)

capture.output(
  res_rf_balanced$confusion_matrix,
  file = "TCGA_balanced_confusion_matrix_random_forest.txt"
)

capture.output(
  res_svm_balanced$confusion_matrix,
  file = "TCGA_balanced_confusion_matrix_radial_svm.txt"
)

############################################################
## Module 15. Plot ROC curves for three models
############################################################

roc_to_df <- function(roc_obj, model_name, auc_value) {
  
  data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    Model = paste0(model_name, " (AUC = ", sprintf("%.3f", auc_value), ")"),
    stringsAsFactors = FALSE
  )
}

roc_df_balanced <- rbind(
  roc_to_df(
    res_glmnet_balanced$roc,
    "Elastic-net",
    res_glmnet_balanced$performance$AUC
  ),
  roc_to_df(
    res_rf_balanced$roc,
    "Random Forest",
    res_rf_balanced$performance$AUC
  ),
  roc_to_df(
    res_svm_balanced$roc,
    "Radial SVM",
    res_svm_balanced$performance$AUC
  )
)

p_roc_balanced <- ggplot(
  roc_df_balanced,
  aes(x = FPR, y = TPR, color = Model)
) +
  geom_line(linewidth = 1.2) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  coord_equal() +
  labs(
    title = "ROC curves in TCGA balanced cohort",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Model"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.65, 0.20),
    legend.background = element_rect(fill = "white", color = "black")
  )

p_roc_balanced

ggsave(
  filename = "TCGA_balanced_ROC_curves_3models.png",
  plot = p_roc_balanced,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "TCGA_balanced_ROC_curves_3models.pdf",
  plot = p_roc_balanced,
  width = 8,
  height = 6
)

############################################################
## Module 16. Plot model performance barplot
## Fixed version using data.table::melt
############################################################
if(F){performance_long <- performance_balanced[, c(
  "Model",
  "AUC",
  "Accuracy",
  "Sensitivity",
  "Specificity",
  "Balanced_Accuracy"
)]

## Convert data.frame to data.table before melt
performance_long <- as.data.table(performance_long)

performance_long <- data.table::melt(
  performance_long,
  id.vars = "Model",
  variable.name = "Metric",
  value.name = "Value"
)

performance_long$Metric <- factor(
  performance_long$Metric,
  levels = c("AUC", "Accuracy", "Sensitivity", "Specificity", "Balanced_Accuracy"),
  labels = c("AUC", "Accuracy", "Sensitivity", "Specificity", "Balanced accuracy")
)

p_perf_balanced <- ggplot(
  performance_long,
  aes(x = Model, y = Value, fill = Metric)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_text(
    aes(label = sprintf("%.3f", Value)),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 3.5
  ) +
  coord_cartesian(ylim = c(0.95, 1.00), clip = "off") +
  labs(
    title = "Model performance in TCGA balanced cohort",
    x = NULL,
    y = "Performance",
    fill = "Metric"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
    plot.margin = margin(10, 20, 10, 10)
  )

p_perf_balanced

ggsave(
  filename = "TCGA_balanced_model_performance_barplot.png",
  plot = p_perf_balanced,
  width = 10,
  height = 6,
  dpi = 300
)}
############################################################
## Alternative: Heatmap-style performance plot
## Recommended for very similar high-performance results
############################################################

library(data.table)
library(ggplot2)

performance_heat <- performance_balanced[, c(
  "Model",
  "AUC",
  "Accuracy",
  "Sensitivity",
  "Specificity",
  "Balanced_Accuracy"
)]

performance_heat <- as.data.table(performance_heat)

performance_heat <- data.table::melt(
  performance_heat,
  id.vars = "Model",
  variable.name = "Metric",
  value.name = "Value"
)

performance_heat$Metric <- factor(
  performance_heat$Metric,
  levels = c("AUC", "Accuracy", "Sensitivity", "Specificity", "Balanced_Accuracy"),
  labels = c("AUC", "Accuracy", "Sensitivity", "Specificity", "Balanced accuracy")
)

performance_heat$Model <- factor(
  performance_heat$Model,
  levels = c("Elastic-net", "Random Forest", "Radial SVM")
)

p_heat_balanced <- ggplot(
  performance_heat,
  aes(x = Metric, y = Model, fill = Value)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = sprintf("%.3f", Value)),
    size = 5
  ) +
  scale_fill_gradient(
    low = "#fdd0c4",
    high = "#ff2610",
    limits = c(0.95, 1.00)
  ) +
  labs(
    title = "Model performance in TCGA balanced cohort",
    x = NULL,
    y = NULL,
    fill = "Performance"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    legend.position = "right"
  )

p_heat_balanced

ggsave(
  filename = "TCGA_balanced_model_performance_heatmap.png",
  plot = p_heat_balanced,
  width = 9,
  height = 4.5,
  dpi = 300
)

ggsave(
  filename = "TCGA_balanced_model_performance_heatmap.pdf",
  plot = p_heat_balanced,
  width = 9,
  height = 4.5
)
############################################################
## Module 17. Save all balanced testing objects
############################################################

save(
  test_df_balanced,
  performance_balanced,
  predictions_balanced,
  roc_df_balanced,
  p_roc_balanced,
  p_perf_balanced,
  res_glmnet_balanced,
  res_rf_balanced,
  res_svm_balanced,
  file = "TCGA_balanced_testing_results_3models.Rdata"
)
cat("\nFinal TCGA balanced model performance:\n")
print(performance_balanced)
```