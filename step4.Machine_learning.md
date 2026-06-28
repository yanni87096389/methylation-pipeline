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

## Machine learning

```{r }
###############################################
## Machine learning using final 239 DMRs
## Goal:1. Use GSE101764 as training cohort
##      2. Use TCGA-COAD + TCGA-READ as external test cohort
##      3. Extract DMR-level beta features
##      4. Train LASSO logistic regression model
##      5. Evaluate model performance in TCGA
############################################################

############################################################
## Module 0. Clear environment and load required packages
############################################################

rm(list = ls())
options(timeout = 100000)
options(scipen = 20)
options(stringsAsFactors = FALSE)

library(data.table)
library(GenomicRanges)
library(IRanges)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(glmnet)
library(pROC)
library(ggplot2)

############################################################
## Module 1. Set input files
############################################################

dmr_file <- "final_DMRs_after_GSE196696_whole_blood_filter.csv"

train_dataset <- "GSE101764"
train_file <- "step0.1-GSE101764_output.Rdata"

test_files <- data.frame(
  dataset = c("TCGA-COAD", "TCGA-READ"),
  file = c(
    "step0.2-TCGA-COAD_output.Rdata",
    "step0.3-TCGA-READ_output.Rdata"
  ),
  stringsAsFactors = FALSE
)

output_dir <- "ML_239DMR_GSE101764_train_TCGA_test_outputs"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

min_feature_cpg <- 6
max_missing_rate <- 0.20

set.seed(2026)

############################################################
## Module 2. Check input files
############################################################

input_files <- c(
  dmr_file,
  train_file,
  test_files$file
)

print(file.exists(input_files))

if (!all(file.exists(input_files))) {
  stop("Some input files do not exist. Please check file names.")
}

############################################################
## Module 3. Load final DMRs
############################################################

dmrs <- fread(
  dmr_file,
  data.table = FALSE
)

dim(dmrs)
head(dmrs)

if (!"DMR_ID" %in% colnames(dmrs)) {
  if ("region_id" %in% colnames(dmrs)) {
    dmrs$DMR_ID <- paste0("DMR_", dmrs$region_id)
  } else {
    dmrs$DMR_ID <- paste0("DMR_", seq_len(nrow(dmrs)))
  }
}

stopifnot(all(c("DMR_ID", "chr", "start", "end") %in% colnames(dmrs)))

if (anyDuplicated(dmrs$DMR_ID) > 0) {
  stop("Duplicated DMR_ID found. Please check the DMR file.")
}

cat("Number of input DMRs for ML:", nrow(dmrs), "\n")

############################################################
## Module 4. Define helper functions
############################################################

fix_chr <- function(x) {
  x <- as.character(x)
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  x[x == "23"] <- "X"
  x[x == "24"] <- "Y"
  x[x %in% c("M", "MT")] <- "M"
  paste0("chr", x)
}

load_myLoad_Rdata <- function(file) {
  
  tmp_env <- new.env()
  
  load(file, envir = tmp_env)
  
  if (!"myLoad" %in% ls(tmp_env)) {
    stop(paste0(file, " does not contain an object named myLoad."))
  }
  
  return(tmp_env$myLoad)
}

standardize_group <- function(group_vector) {
  
  group_raw <- as.character(group_vector)
  group_low <- tolower(group_raw)
  
  group_new <- group_raw
  
  group_new[
    grepl("tumor|tumour|cancer|carcinoma|crc|coad|read|adenocarcinoma|primary", group_low)
  ] <- "Tumor"
  
  group_new[
    grepl("normal|control|healthy|adjacent|mucosa|cancer-free|non-cancer", group_low)
  ] <- "Normal"
  
  return(group_new)
}

safe_divide <- function(a, b) {
  ifelse(b == 0, NA_real_, a / b)
}

calculate_metrics <- function(true_label, pred_label, prob) {
  
  true_label <- factor(true_label, levels = c("Normal", "Tumor"))
  pred_label <- factor(pred_label, levels = c("Normal", "Tumor"))
  
  cm <- table(
    True = true_label,
    Predicted = pred_label
  )
  
  TN <- cm["Normal", "Normal"]
  FP <- cm["Normal", "Tumor"]
  FN <- cm["Tumor", "Normal"]
  TP <- cm["Tumor", "Tumor"]
  
  accuracy <- safe_divide(TP + TN, sum(cm))
  sensitivity <- safe_divide(TP, TP + FN)
  specificity <- safe_divide(TN, TN + FP)
  ppv <- safe_divide(TP, TP + FP)
  npv <- safe_divide(TN, TN + FN)
  
  auc_value <- as.numeric(
    auc(
      roc(
        response = ifelse(true_label == "Tumor", 1, 0),
        predictor = prob,
        quiet = TRUE
      )
    )
  )
  
  metrics <- data.frame(
    AUC = auc_value,
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    PPV = ppv,
    NPV = npv,
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      confusion_matrix = cm,
      metrics = metrics
    )
  )
}

make_stratified_foldid <- function(y, nfolds = 5) {
  
  y <- as.factor(y)
  foldid <- rep(NA_integer_, length(y))
  
  for (class_i in levels(y)) {
    
    idx <- which(y == class_i)
    idx <- sample(idx)
    
    foldid[idx] <- rep(
      seq_len(nfolds),
      length.out = length(idx)
    )
  }
  
  return(foldid)
}

############################################################
## Module 5. Prepare DMR coordinates
############################################################

dmrs$chr <- fix_chr(dmrs$chr)
dmrs$start <- as.numeric(dmrs$start)
dmrs$end <- as.numeric(dmrs$end)

dmrs <- dmrs[
  !is.na(dmrs$chr) &
    !is.na(dmrs$start) &
    !is.na(dmrs$end),
]

dmrs <- dmrs[
  dmrs$chr %in% paste0("chr", 1:22),
]

dmr_gr <- GRanges(
  seqnames = dmrs$chr,
  ranges = IRanges(
    start = dmrs$start,
    end = dmrs$end
  ),
  DMR_ID = dmrs$DMR_ID
)

cat("Number of final DMRs used as ML features:", nrow(dmrs), "\n")

############################################################
## Module 6. Load 450K annotation
############################################################

anno <- minfi::getAnnotation(
  IlluminaHumanMethylation450kanno.ilmn12.hg19
)

anno <- as.data.frame(anno)
anno$CpG <- rownames(anno)

anno$chr <- fix_chr(anno$chr)
anno$pos <- as.numeric(anno$pos)

anno <- anno[
  !is.na(anno$chr) &
    !is.na(anno$pos) &
    anno$chr %in% paste0("chr", 1:22),
]

############################################################
## Module 7. Function to extract DMR-level beta features
############################################################

extract_dmr_features <- function(dataset_name, rdata_file) {
  
  cat("\n############################################################\n")
  cat("Extracting DMR-level features from:", dataset_name, "\n")
  cat("############################################################\n")
  
  ############################################################
  ## 7.1 Load ChAMP-filtered object
  ############################################################
  
  myLoad <- load_myLoad_Rdata(rdata_file)
  
  beta <- myLoad$beta
  pd <- myLoad$pd
  
  stopifnot(identical(colnames(beta), rownames(pd)))
  
  pd$group_original <- as.character(pd$group)
  pd$group <- standardize_group(pd$group_original)
  
  cat("Original group table:\n")
  print(table(pd$group_original))
  
  cat("Standardized group table:\n")
  print(table(pd$group))
  
  ############################################################
  ## 7.2 Keep only Tumor and Normal samples
  ############################################################
  
  keep_samples <- pd$group %in% c("Normal", "Tumor")
  
  beta <- beta[, keep_samples, drop = FALSE]
  pd <- pd[keep_samples, , drop = FALSE]
  
  stopifnot(identical(colnames(beta), rownames(pd)))
  
  pd$group <- factor(
    pd$group,
    levels = c("Normal", "Tumor")
  )
  
  cat("Final group table used for ML:\n")
  print(table(pd$group))
  
  if (length(unique(pd$group)) < 2) {
    stop(paste0(dataset_name, " does not contain both Normal and Tumor samples."))
  }
  
  ############################################################
  ## 7.3 Map DMRs to CpGs available in this dataset
  ############################################################
  
  anno_use <- anno[
    anno$CpG %in% rownames(beta),
  ]
  
  common_cpgs <- intersect(
    rownames(beta),
    anno_use$CpG
  )
  
  beta <- beta[common_cpgs, , drop = FALSE]
  anno_use <- anno_use[match(common_cpgs, anno_use$CpG), ]
  
  stopifnot(identical(rownames(beta), anno_use$CpG))
  
  cpg_gr <- GRanges(
    seqnames = anno_use$chr,
    ranges = IRanges(
      start = anno_use$pos,
      end = anno_use$pos
    ),
    CpG = anno_use$CpG
  )
  
  hits <- findOverlaps(dmr_gr, cpg_gr)
  
  cat("Number of DMR-CpG overlaps:", length(hits), "\n")
  
  dmr_cpg_list <- split(
    mcols(cpg_gr)$CpG[subjectHits(hits)],
    queryHits(hits)
  )
  
  ############################################################
  ## 7.4 Calculate DMR-level beta values
  ############################################################
  
  feature_matrix <- matrix(
    NA_real_,
    nrow = ncol(beta),
    ncol = nrow(dmrs),
    dimnames = list(
      colnames(beta),
      dmrs$DMR_ID
    )
  )
  
  feature_info <- data.frame(
    DMR_ID = dmrs$DMR_ID,
    chr = dmrs$chr,
    start = dmrs$start,
    end = dmrs$end,
    dataset = dataset_name,
    nCpG_feature = 0,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(dmrs))) {
    
    cpgs <- dmr_cpg_list[[as.character(i)]]
    
    if (is.null(cpgs)) {
      next
    }
    
    cpgs <- unique(intersect(cpgs, rownames(beta)))
    
    feature_info$nCpG_feature[i] <- length(cpgs)
    
    if (length(cpgs) < min_feature_cpg) {
      next
    }
    
    dmr_value <- colMeans(
      beta[cpgs, , drop = FALSE],
      na.rm = TRUE
    )
    
    feature_matrix[, i] <- dmr_value
  }
  
  feature_df <- as.data.frame(feature_matrix)
  
  phenotype_df <- data.frame(
    Sample_ID = rownames(pd),
    Group = as.character(pd$group),
    Dataset = dataset_name,
    Group_original = as.character(pd$group_original),
    stringsAsFactors = FALSE
  )
  
  feature_df$Sample_ID <- rownames(feature_df)
  
  feature_df <- merge(
    phenotype_df,
    feature_df,
    by = "Sample_ID",
    all.x = TRUE
  )
  
  cat("Feature matrix dimension:\n")
  print(dim(feature_df))
  
  cat("DMRs with at least one CpG overlap:\n")
  print(sum(feature_info$nCpG_feature >= min_feature_cpg))
  
  return(
    list(
      feature_df = feature_df,
      feature_info = feature_info
    )
  )
}

############################################################
## Module 8. Extract training features from GSE101764
############################################################

train_output <- extract_dmr_features(
  dataset_name = train_dataset,
  rdata_file = train_file
)

train_feature_df <- train_output$feature_df
train_feature_info <- train_output$feature_info

dim(train_feature_df)
table(train_feature_df$Group)

write.csv(
  train_feature_df,
  file = file.path(output_dir, "GSE101764_239DMR_training_feature_matrix.csv"),
  row.names = FALSE
)

write.csv(
  train_feature_info,
  file = file.path(output_dir, "GSE101764_239DMR_feature_CpG_info.csv"),
  row.names = FALSE
)

############################################################
## Module 9. Extract external test features from TCGA
############################################################

test_outputs <- lapply(seq_len(nrow(test_files)), function(i) {
  
  extract_dmr_features(
    dataset_name = test_files$dataset[i],
    rdata_file = test_files$file[i]
  )
  
})

test_feature_df <- rbindlist(
  lapply(test_outputs, function(x) x$feature_df),
  fill = TRUE
)

test_feature_df <- as.data.frame(test_feature_df)

test_feature_info <- rbindlist(
  lapply(test_outputs, function(x) x$feature_info),
  fill = TRUE
)

test_feature_info <- as.data.frame(test_feature_info)

dim(test_feature_df)
table(test_feature_df$Group)
table(test_feature_df$Dataset)

write.csv(
  test_feature_df,
  file = file.path(output_dir, "TCGA_COAD_READ_239DMR_external_test_feature_matrix.csv"),
  row.names = FALSE
)

write.csv(
  test_feature_info,
  file = file.path(output_dir, "TCGA_COAD_READ_239DMR_feature_CpG_info.csv"),
  row.names = FALSE
)

############################################################
## Module 10. Prepare ML matrices
############################################################

meta_cols <- c(
  "Sample_ID",
  "Group",
  "Dataset",
  "Group_original"
)

feature_cols_train <- setdiff(colnames(train_feature_df), meta_cols)
feature_cols_test <- setdiff(colnames(test_feature_df), meta_cols)

common_features <- intersect(
  feature_cols_train,
  feature_cols_test
)

cat("Common DMR features between training and test:", length(common_features), "\n")

X_train_raw <- train_feature_df[, common_features, drop = FALSE]
X_test_raw <- test_feature_df[, common_features, drop = FALSE]

y_train <- factor(
  train_feature_df$Group,
  levels = c("Normal", "Tumor")
)

y_test <- factor(
  test_feature_df$Group,
  levels = c("Normal", "Tumor")
)

############################################################
## Module 11. Remove features with high missingness
############################################################

feature_missing_rate_train <- colMeans(is.na(X_train_raw))
feature_missing_rate_test <- colMeans(is.na(X_test_raw))

keep_features <- names(feature_missing_rate_train)[
  feature_missing_rate_train <= max_missing_rate &
    feature_missing_rate_test <= max_missing_rate
]

X_train_raw <- X_train_raw[, keep_features, drop = FALSE]
X_test_raw <- X_test_raw[, keep_features, drop = FALSE]

cat("Features retained after missingness filtering:", ncol(X_train_raw), "\n")

############################################################
## Module 12. Remove zero-variance features
############################################################

feature_sd <- apply(
  X_train_raw,
  2,
  sd,
  na.rm = TRUE
)

keep_nonzero_features <- names(feature_sd)[
  !is.na(feature_sd) &
    feature_sd > 0
]

X_train_raw <- X_train_raw[, keep_nonzero_features, drop = FALSE]
X_test_raw <- X_test_raw[, keep_nonzero_features, drop = FALSE]

cat("Features retained after zero-variance filtering:", ncol(X_train_raw), "\n")

############################################################
## Module 13. Median imputation using training-set medians
############################################################

train_median <- apply(
  X_train_raw,
  2,
  median,
  na.rm = TRUE
)

for (j in seq_len(ncol(X_train_raw))) {
  
  X_train_raw[is.na(X_train_raw[, j]), j] <- train_median[j]
  X_test_raw[is.na(X_test_raw[, j]), j] <- train_median[j]
}

############################################################
## Module 14. Standardize features using training-set mean and SD
############################################################

train_mean <- apply(
  X_train_raw,
  2,
  mean,
  na.rm = TRUE
)

train_sd <- apply(
  X_train_raw,
  2,
  sd,
  na.rm = TRUE
)

train_sd[train_sd == 0 | is.na(train_sd)] <- 1

X_train <- scale(
  X_train_raw,
  center = train_mean,
  scale = train_sd
)

X_test <- scale(
  X_test_raw,
  center = train_mean,
  scale = train_sd
)

X_train <- as.matrix(X_train)
X_test <- as.matrix(X_test)

############################################################
## Module 15. Encode labels
############################################################

y_train_binary <- ifelse(y_train == "Tumor", 1, 0)
y_test_binary <- ifelse(y_test == "Tumor", 1, 0)

table(y_train)
table(y_test)

############################################################
## Module 16. Train LASSO logistic regression model
############################################################

class_table <- table(y_train_binary)

sample_weights <- ifelse(
  y_train_binary == 1,
  length(y_train_binary) / (2 * class_table["1"]),
  length(y_train_binary) / (2 * class_table["0"])
)

sample_weights <- as.numeric(sample_weights)

foldid <- make_stratified_foldid(
  y = y_train_binary,
  nfolds = 5
)

cv_fit <- cv.glmnet(
  x = X_train,
  y = y_train_binary,
  family = "binomial",
  alpha = 1,
  type.measure = "auc",
  nfolds = 5,
  foldid = foldid,
  weights = sample_weights
)

best_lambda <- cv_fit$lambda.1se

final_model <- glmnet(
  x = X_train,
  y = y_train_binary,
  family = "binomial",
  alpha = 1,
  lambda = best_lambda,
  weights = sample_weights
)

############################################################
## Module 17. If lambda.1se selects no DMRs, use lambda.min
############################################################

coef_mat_tmp <- coef(final_model)

selected_tmp <- rownames(coef_mat_tmp)[
  as.numeric(coef_mat_tmp[, 1]) != 0
]

selected_tmp <- setdiff(selected_tmp, "(Intercept)")

if (length(selected_tmp) == 0) {
  
  cat("lambda.1se selected no DMR features. Switching to lambda.min.\n")
  
  best_lambda <- cv_fit$lambda.min
  
  final_model <- glmnet(
    x = X_train,
    y = y_train_binary,
    family = "binomial",
    alpha = 1,
    lambda = best_lambda,
    weights = sample_weights
  )
}

cat("Final lambda used:", best_lambda, "\n")

############################################################
## Module 18. Training-set prediction
############################################################

train_prob <- as.numeric(
  predict(
    final_model,
    newx = X_train,
    type = "response"
  )
)

train_roc <- roc(
  response = y_train_binary,
  predictor = train_prob,
  quiet = TRUE
)

train_auc <- as.numeric(auc(train_roc))

cat("Training AUC:", train_auc, "\n")

############################################################
## Module 19. Select threshold using Youden index in training set
############################################################

train_coords <- coords(
  train_roc,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity")
)

best_threshold <- as.numeric(train_coords["threshold"])

cat("Best threshold from training set:", best_threshold, "\n")

############################################################
## Module 20. TCGA external test prediction
############################################################

test_prob <- as.numeric(
  predict(
    final_model,
    newx = X_test,
    type = "response"
  )
)

test_pred <- ifelse(
  test_prob >= best_threshold,
  "Tumor",
  "Normal"
)

test_pred <- factor(
  test_pred,
  levels = c("Normal", "Tumor")
)

test_roc <- roc(
  response = y_test_binary,
  predictor = test_prob,
  quiet = TRUE
)

test_auc <- as.numeric(auc(test_roc))

cat("External TCGA test AUC:", test_auc, "\n")

############################################################
## Module 21. Calculate performance metrics
############################################################

train_pred <- ifelse(
  train_prob >= best_threshold,
  "Tumor",
  "Normal"
)

train_pred <- factor(
  train_pred,
  levels = c("Normal", "Tumor")
)

train_eval <- calculate_metrics(
  true_label = y_train,
  pred_label = train_pred,
  prob = train_prob
)

test_eval <- calculate_metrics(
  true_label = y_test,
  pred_label = test_pred,
  prob = test_prob
)

train_metrics <- train_eval$metrics
train_metrics$Cohort <- train_dataset

test_metrics <- test_eval$metrics
test_metrics$Cohort <- "TCGA-COAD_READ"

performance_summary <- rbind(
  train_metrics,
  test_metrics
)

cat("\nPerformance summary:\n")
print(performance_summary)

cat("\nExternal TCGA confusion matrix:\n")
print(test_eval$confusion_matrix)

############################################################
## Module 22. TCGA project-specific performance
############################################################

test_prediction_df <- data.frame(
  Sample_ID = test_feature_df$Sample_ID,
  Group = test_feature_df$Group,
  Dataset = test_feature_df$Dataset,
  Group_original = test_feature_df$Group_original,
  Predicted_prob_Tumor = test_prob,
  Predicted_label = as.character(test_pred),
  stringsAsFactors = FALSE
)

project_metrics_list <- list()

for (proj in unique(test_prediction_df$Dataset)) {
  
  tmp <- test_prediction_df[
    test_prediction_df$Dataset == proj,
  ]
  
  if (length(unique(tmp$Group)) < 2) {
    next
  }
  
  tmp_eval <- calculate_metrics(
    true_label = tmp$Group,
    pred_label = factor(
      tmp$Predicted_label,
      levels = c("Normal", "Tumor")
    ),
    prob = tmp$Predicted_prob_Tumor
  )
  
  tmp_metrics <- tmp_eval$metrics
  tmp_metrics$Cohort <- proj
  
  project_metrics_list[[proj]] <- tmp_metrics
}

project_metrics <- rbindlist(
  project_metrics_list,
  fill = TRUE
)

performance_summary <- rbindlist(
  list(
    as.data.table(performance_summary),
    project_metrics
  ),
  fill = TRUE
)

performance_summary <- as.data.frame(performance_summary)

print(performance_summary)

############################################################
## Module 23. Extract selected DMR features
############################################################

coef_mat <- coef(final_model)

selected_features <- data.frame(
  Feature = rownames(coef_mat),
  Coefficient = as.numeric(coef_mat[, 1]),
  stringsAsFactors = FALSE
)

selected_features <- selected_features[
  selected_features$Coefficient != 0,
]

selected_features <- selected_features[
  selected_features$Feature != "(Intercept)",
]

selected_features <- selected_features[
  order(abs(selected_features$Coefficient), decreasing = TRUE),
]

cat("Number of selected DMR features by LASSO:", nrow(selected_features), "\n")

selected_features <- merge(
  selected_features,
  dmrs,
  by.x = "Feature",
  by.y = "DMR_ID",
  all.x = TRUE
)

head(selected_features)

############################################################
## Module 24. Save ML outputs
############################################################

write.csv(
  performance_summary,
  file = file.path(output_dir, "ML_239DMR_performance_summary.csv"),
  row.names = FALSE
)

write.csv(
  test_prediction_df,
  file = file.path(output_dir, "ML_239DMR_TCGA_external_test_predictions.csv"),
  row.names = FALSE
)

write.csv(
  selected_features,
  file = file.path(output_dir, "ML_239DMR_LASSO_selected_features.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    Feature = colnames(X_train),
    stringsAsFactors = FALSE
  ),
  file = file.path(output_dir, "ML_239DMR_features_used_in_model.csv"),
  row.names = FALSE
)

save(
  dmrs,
  train_feature_df,
  test_feature_df,
  train_feature_info,
  test_feature_info,
  X_train_raw,
  X_test_raw,
  X_train,
  X_test,
  y_train,
  y_test,
  y_train_binary,
  y_test_binary,
  cv_fit,
  final_model,
  best_lambda,
  best_threshold,
  train_prob,
  test_prob,
  train_eval,
  test_eval,
  performance_summary,
  selected_features,
  test_prediction_df,
  train_mean,
  train_sd,
  train_median,
  file = file.path(output_dir, "step1-ML-239DMR-GSE101764-train-TCGA-test-output.Rdata")
)

############################################################
## Module 25. Visualization
############################################################

############################################################
## Figure 1. ROC curves
############################################################

roc_train_df <- data.frame(
  False_positive_rate = 1 - train_roc$specificities,
  Sensitivity = train_roc$sensitivities,
  Cohort = paste0(train_dataset, " training")
)

roc_test_df <- data.frame(
  False_positive_rate = 1 - test_roc$specificities,
  Sensitivity = test_roc$sensitivities,
  Cohort = "TCGA external test"
)

roc_df <- rbind(
  roc_train_df,
  roc_test_df
)

p_roc <- ggplot(
  roc_df,
  aes(
    x = False_positive_rate,
    y = Sensitivity,
    linetype = Cohort
  )
) +
  geom_line(linewidth = 1) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "ROC curves for 239-DMR LASSO model",
    subtitle = paste0(
      train_dataset,
      " training AUC = ",
      round(train_auc, 3),
      "; TCGA test AUC = ",
      round(test_auc, 3)
    ),
    x = "False positive rate",
    y = "True positive rate",
    linetype = NULL
  )

print(p_roc)

ggsave(
  filename = file.path(output_dir, "ML_239DMR_ROC_curve.png"),
  plot = p_roc,
  width = 6.5,
  height = 5.5,
  dpi = 300
)

############################################################
## Figure 2. Predicted tumor probability in TCGA
############################################################

p_prob <- ggplot(
  test_prediction_df,
  aes(
    x = Group,
    y = Predicted_prob_Tumor
  )
) +
  geom_boxplot(width = 0.45, outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.7, size = 1.8) +
  geom_hline(
    yintercept = best_threshold,
    linetype = "dashed"
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Predicted tumor probability in TCGA external test cohort",
    subtitle = paste0("Threshold selected from ", train_dataset),
    x = NULL,
    y = "Predicted probability of Tumor"
  )

print(p_prob)

ggsave(
  filename = file.path(output_dir, "ML_239DMR_TCGA_predicted_probability_boxplot.png"),
  plot = p_prob,
  width = 5.5,
  height = 5,
  dpi = 300
)

############################################################
## Figure 3. Predicted tumor probability by TCGA project
############################################################

p_prob_project <- ggplot(
  test_prediction_df,
  aes(
    x = Dataset,
    y = Predicted_prob_Tumor,
    shape = Group
  )
) +
  geom_boxplot(width = 0.45, outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.7, size = 1.8) +
  geom_hline(
    yintercept = best_threshold,
    linetype = "dashed"
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Predicted tumor probability by TCGA project",
    x = NULL,
    y = "Predicted probability of Tumor",
    shape = NULL
  )

print(p_prob_project)

ggsave(
  filename = file.path(output_dir, "ML_239DMR_TCGA_predicted_probability_by_project.png"),
  plot = p_prob_project,
  width = 6.5,
  height = 5,
  dpi = 300
)

############################################################
## Figure 4. Top selected DMR coefficients
############################################################

top_coef <- head(selected_features, 20)

if (nrow(top_coef) > 0) {
  
  top_coef$Feature <- factor(
    top_coef$Feature,
    levels = rev(top_coef$Feature)
  )
  
  p_coef <- ggplot(
    top_coef,
    aes(
      x = Feature,
      y = Coefficient
    )
  ) +
    geom_col(width = 0.7, fill = "grey70", color = "black") +
    coord_flip() +
    theme_bw(base_size = 13) +
    labs(
      title = "Top selected DMR features from LASSO model",
      x = NULL,
      y = "Model coefficient"
    )
  
  print(p_coef)
  
  ggsave(
    filename = file.path(output_dir, "ML_239DMR_top_LASSO_coefficients.png"),
    plot = p_coef,
    width = 7,
    height = 6,
    dpi = 300
  )
}

```


## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
