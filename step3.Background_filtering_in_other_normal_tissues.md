---
title: "step3.background_filtering_in_other_normal_tissues"
output: html_document
date: "2026-06-26"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## step0.data preprocessing
## GSE213478 850k dataset
```{r cars}
###############################################
## GSE213478 850k methylation processed matrix pipeline
## Goal: 1. Get standard beta value matrix
##       2. Build phenotype table pD
##       3. Run ChAMP filtering
############################################################

############################################################
## Module 0. Clear environment and load required packages
############################################################

rm(list = ls())
options(timeout = 100000) 
options(scipen = 20)
options(stringsAsFactors = FALSE)

require(GEOquery)
require(Biobase)
library(impute)
library(ChAMP)
library(data.table)

############################################################
## Module 1. Read GEO metadata and beta matrix
############################################################

## Download/read GEO series matrix metadata
eset <- getGEO("GSE213478", destdir = "./", AnnotGPL = TRUE, getGPL = FALSE)

## Read processed beta matrix
beta_213478 <- fread(
  "GSE213478_methylation_DNAm_noob_final_BMIQ_all_tissues_987.txt.gz",
  data.table = FALSE
)
## Check beta matrix
beta_213478[1:5, 1:5]
class(beta_213478)
dim(beta_213478)

## Set CpG probe IDs as row names
rownames(beta_213478) <- beta_213478[, 1]

## Remove the first CpG ID column
beta_213478 <- beta_213478[, -1]

## Check beta matrix after cleaning
class(beta_213478)
dim(beta_213478)
head(rownames(beta_213478))
head(colnames(beta_213478))

############################################################
## Module 2. Match beta matrix column names to GSM IDs
############################################################

## Extract phenotype table from GEO object
pD.all <- pData(eset[[1]])

## Check tissue distribution
table(pD.all$`tissue:ch1`)

## Beta matrix column names are GTEx sample IDs
beta_sample_ids <- colnames(beta_213478)

## Match beta matrix columns to pD.all$title
match_index <- match(beta_sample_ids, pD.all$title)

## Check whether all beta columns can be matched
sum(is.na(match_index))

## Show unmatched samples if any
beta_sample_ids[is.na(match_index)]

## Stop if any beta matrix column cannot be matched
if (any(is.na(match_index))) {
  stop("Some beta matrix columns cannot be matched to pD.all$title.")
}

## Build a mapping table between original beta column names and GSM accessions
colname_GSM_map <- data.frame(
  original_colname = beta_sample_ids,
  GSM = pD.all$geo_accession[match_index],
  stringsAsFactors = FALSE
)

## Check the mapping table
head(colname_GSM_map)
dim(colname_GSM_map)

## Rename beta matrix columns from GTEx IDs to GSM IDs
colnames(beta_213478) <- colname_GSM_map$GSM

## Check whether the new column names are all GSM IDs
all(grepl("^GSM", colnames(beta_213478)))

## Use the cleaned beta matrix as beta.m
beta.m <- beta_213478

############################################################
## Module 3. Build phenotype table pD
############################################################

## Reorder phenotype table according to beta matrix column order
pD.matched <- pD.all[match_index, ]

## Check the column index for tissue information
which(colnames(pD.matched) == "tissue:ch1")

## Build a clean phenotype table
## Column 57 = tissue type
## Column 2  = GSM accession
pD <- pD.matched[, c(57, 2)]

## Rename columns
names(pD) <- c("group", "ID")

## Standardize tissue names
pD$group <- ifelse(pD$group == "Breast - Mammary Tissue", "Breast",
                   ifelse(pD$group == "Colon - Transverse", "Colon",
                          ifelse(pD$group == "Kidney - Cortex", "Kidney",
                                 ifelse(pD$group == "Lung", "Lung",
                                        ifelse(pD$group == "Muscle - Skeletal", "Muscle",
                                               ifelse(pD$group == "Ovary", "Ovary",
                                                      ifelse(pD$group == "Prostate", "Prostate",
                                                             ifelse(pD$group == "Testis", "Testis",
                                                                    ifelse(pD$group == "Whole Blood", "Blood", NA)))))))))

## Convert ID and group to character
pD$ID <- as.character(pD$ID)
pD$group <- as.character(pD$group)

## Add sample type information
pD$sample_type <- "Normal"

## Set GSM IDs as row names
rownames(pD) <- pD$ID

## Check if the extracted groupings are correct
check_df <- data.frame(
  Original_colname = colname_GSM_map$original_colname,
  GSM = pD$ID,
  Tissue_raw = pD.matched$`tissue:ch1`,
  Group = pD$group,
  Sample_type = pD$sample_type,
  stringsAsFactors = FALSE
)

head(check_df, 10)
tail(check_df, 10)

## Check tissue group counts
table(pD$group)

## Check sample type counts
table(pD$sample_type)

## Check whether there are NA groups
table(is.na(pD$group))

if (any(is.na(pD$group))) {
  stop("Some tissue groups were converted to NA. Please check tissue name standardization.")
}

############################################################
## Module 4. Align pD with beta.m
############################################################

## Check whether beta.m and pD contain the same sample set
setequal(colnames(beta.m), pD$ID)

## Check duplicated sample names
sum(duplicated(colnames(beta.m)))
sum(duplicated(pD$ID))

## Strict check before sorting
all(colnames(beta.m) == pD$ID)

## Sort pD by group
sample_order <- order(pD$group)

pD <- pD[sample_order, , drop = FALSE]

## Reorder beta matrix columns according to pD
beta.m <- beta.m[, rownames(pD), drop = FALSE]

## Final strict alignment check
dim(beta.m)
dim(pD)

identical(colnames(beta.m), rownames(pD))
identical(colnames(beta.m), pD$ID)

## Convert beta.m to numeric matrix for ChAMP
beta.m <- as.matrix(beta.m)
mode(beta.m) <- "numeric"

## Final beta matrix check
class(beta.m)
dim(beta.m)

############################################################
## Module 5. Run ChAMP filtering
############################################################

## Run ChAMP filtering
## GSE213478 is EPIC/850K data, so arraytype should be "EPIC"
myLoad <- champ.filter(
  beta = beta.m,
  pd = pD,
  arraytype = "EPIC"
)
## Save the ChAMP object for downstream analysis.
saveRDS(myLoad, file = "step0.1-GSE213478_output.rds")
save(myLoad,file = "step0.1-GSE213478_output.Rdata")
#myLoad <- readRDS("step0.1-GSE213478_output.rds")
```

## GSE48472 450k dataset
```{r}
###############################################
## GSE48472 450k methylation processed matrix pipeline
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
eset <- getGEO("GSE48472",destdir = './',AnnotGPL = T,getGPL = F)
beta.m <- exprs(eset[[1]])
class(beta.m)
dim(beta.m)

############################################################
## Module 2. Build phenotype table pD
## Modify this part when applying the pipeline to another dataset
pD.all <- pData(eset[[1]])
colnames(pD.all)
table(pD.all$`tissue:ch1`)
which(colnames(pD.all)=="tissue:ch1")
pD <- pD.all[,c(33,2)]
names(pD) <- c("group", "ID")
## Standardize tissue names
pD$group <- ifelse(pD$group == "blood", "Blood",
                  ifelse(pD$group == "buccal", "Buccal",
                  ifelse(pD$group == "hair", "Hair",
                  ifelse(pD$group == "liver", "Liver",
                  ifelse(pD$group == "muscle", "Muscle",
                  ifelse(pD$group == "omentum", "Omentum",
                  ifelse(pD$group == "pancreas", "Pancreas",
                  ifelse(pD$group == "saliva", "Saliva",
                  ifelse(pD$group == "scfat", "SC_Fat",
                  ifelse(pD$group == "spleen", "Spleen", NA))))))))))
## Check if the extracted groupings are correct
check_df <- data.frame(
  Title = pD.all$title,
  Tissue_raw = pD.all$`tissue:ch1`,
  Group = pD$group,
  ID = pD$ID
)

head(check_df, 10)
tail(check_df, 10)
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

# Set pD row names as GSM IDs.
rownames(pD) <- pD$ID
## Final strict alignment check.
dim(beta.m)
dim(pD)

# Set group order
#pD$group <- factor(pD$group,levels = c("Normal", "Tumor"))

# Sort pD by group
sample_order <- order(pD$group)

pD <- pD[sample_order, , drop = FALSE]

# Reorder beta matrix columns according to pD
beta.m <- beta.m[, rownames(pD), drop = FALSE]

############################################################
## Module 4. Run ChAMP filtering
############################################################

## Run ChAMP filtering using the cleaned beta matrix and aligned pD.
myLoad <- champ.filter(beta = beta.m, pd = pD, arraytype = "450k")
## Save the ChAMP object for downstream analysis.
save(myLoad, file = "step0.2-GSE48472_output.Rdata")
saveRDS(myLoad, file = "step0.2-GSE48472_output.rds")
#myLoad <- readRDS("step0.2-GSE48472_output.rds")


```

## step1.Background filtering in other normal tissues
```{r}
## Background filtering of robust validated hyper-DMRs
## Goal:
##   Keep robust validated hyper-DMRs that show low methylation
##   in normal/background tissues.
## Main filtering criterion:
##   nCpG_background > 5
##   max_group_mean_beta < 0.20
## Input:
##   1. robust_validated_hyper_DMRs_at_least_3datasets.csv
##   2. step0.1-GSE213478_output.Rdata
##   3. step0.2-GSE48472_output.Rdata
##
## Output:
##   1. background_DMR_beta_summary_by_dataset.csv
##   2. background_DMR_beta_by_tissue_group.csv
##   3. background_cross_dataset_summary.csv
##   4. background_filtered_DMRs_beta_lt_0.2_nCpG_gt5.csv
##   5. step2-background-filtering-output.Rdata
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
options(scipen = 20)

############################################################
## 0. Load required packages
############################################################

library(data.table)
library(GenomicRanges)
library(IRanges)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

############################################################
## 1. Set filtering parameters
############################################################

beta_threshold <- 0.20

## Discovery DMR criterion was no.cpgs > 5.
## Therefore, for background filtering, we keep the same CpG coverage logic:
## nCpG_background > 5, equivalent to nCpG_background >= 6.
min_bg_cpg <- 6

############################################################
## 2. Load robust validated hyper-DMRs
############################################################

dmrs <- fread(
  "robust_validated_hyper_DMRs_at_least_3datasets.csv",
  data.table = FALSE
)

## Required coordinate columns
stopifnot(all(c("chr", "start", "end") %in% colnames(dmrs)))

## Add DMR ID if not already present
if (!"DMR_ID" %in% colnames(dmrs)) {
  if ("region_id" %in% colnames(dmrs)) {
    dmrs$DMR_ID <- paste0("DMR_", dmrs$region_id)
  } else {
    dmrs$DMR_ID <- paste0("DMR_", seq_len(nrow(dmrs)))
  }
}

## Keep discovery delta beta if available
if (!"discovery_deltaBeta" %in% colnames(dmrs)) {
  if ("mean_deltaBeta" %in% colnames(dmrs)) {
    dmrs$discovery_deltaBeta <- dmrs$mean_deltaBeta
  } else {
    dmrs$discovery_deltaBeta <- NA_real_
  }
}

############################################################
## 3. Standardize genomic coordinates
############################################################

fix_chr <- function(x) {
  x <- as.character(x)
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  x[x == "23"] <- "X"
  x[x == "24"] <- "Y"
  x[x %in% c("M", "MT")] <- "M"
  paste0("chr", x)
}

dmrs$chr <- fix_chr(dmrs$chr)
dmrs$start <- as.numeric(dmrs$start)
dmrs$end <- as.numeric(dmrs$end)

dmrs <- dmrs[
  !is.na(dmrs$chr) &
    !is.na(dmrs$start) &
    !is.na(dmrs$end),
]

dmrs <- dmrs[dmrs$chr %in% paste0("chr", 1:22), ]

dmr_gr <- GRanges(
  seqnames = dmrs$chr,
  ranges = IRanges(start = dmrs$start, end = dmrs$end),
  DMR_ID = dmrs$DMR_ID
)

############################################################
## 4. Function to load .Rdata file containing myLoad
############################################################

load_myLoad_Rdata <- function(file) {
  
  tmp_env <- new.env()
  
  load(file, envir = tmp_env)
  
  if (!"myLoad" %in% ls(tmp_env)) {
    stop(paste0(file, " does not contain an object named myLoad."))
  }
  
  return(tmp_env$myLoad)
}

############################################################
## 5. Function to load array annotation
############################################################

get_array_annotation <- function(arraytype) {
  
  if (arraytype == "450K") {
    
    anno <- minfi::getAnnotation(
      IlluminaHumanMethylation450kanno.ilmn12.hg19
    )
    
  } else if (arraytype == "EPIC") {
    
    anno <- minfi::getAnnotation(
      IlluminaHumanMethylationEPICanno.ilm10b4.hg19
    )
    
  } else {
    
    stop("arraytype must be either '450K' or 'EPIC'.")
  }
  
  anno <- as.data.frame(anno)
  anno$CpG <- rownames(anno)
  
  anno$chr <- fix_chr(anno$chr)
  anno$pos <- as.numeric(anno$pos)
  
  anno <- anno[
    !is.na(anno$chr) &
      !is.na(anno$pos) &
      anno$chr %in% paste0("chr", 1:22),
  ]
  
  return(anno)
}

############################################################
## 6. Function to background-filter one dataset
############################################################

background_filter_one_dataset <- function(dataset_name, rdata_file, arraytype) {
  
  cat("\n############################################################\n")
  cat("Background filtering dataset:", dataset_name, "\n")
  cat("############################################################\n")
  
  ############################################################
  ## 6.1 Load ChAMP-filtered object
  ############################################################
  
  myLoad <- load_myLoad_Rdata(rdata_file)
  
  beta_bg <- myLoad$beta
  pd_bg <- myLoad$pd
  
  ## Check sample alignment
  stopifnot(identical(colnames(beta_bg), rownames(pd_bg)))
  
  pd_bg$group <- as.character(pd_bg$group)
  
  ############################################################
  ## 6.2 Select normal/background samples
  ############################################################
  ## If the dataset has Normal/Tumor labels, keep only Normal.
  ## If the dataset contains normal tissue names, treat all groups as
  ## normal background tissue groups.
  
  if (all(c("Normal", "Tumor") %in% unique(pd_bg$group))) {
    
    keep_samples <- pd_bg$group == "Normal"
    
  } else {
    
    keep_samples <- !is.na(pd_bg$group)
  }
  
  beta_bg <- beta_bg[, keep_samples, drop = FALSE]
  pd_bg <- pd_bg[keep_samples, , drop = FALSE]
  
  stopifnot(identical(colnames(beta_bg), rownames(pd_bg)))
  
  cat("Normal/background sample groups:\n")
  print(table(pd_bg$group))
  
  ############################################################
  ## 6.3 Map DMRs to CpGs available in this dataset
  ############################################################
  
  anno <- get_array_annotation(arraytype)
  
  anno <- anno[anno$CpG %in% rownames(beta_bg), ]
  
  common_cpgs <- intersect(rownames(beta_bg), anno$CpG)
  
  beta_bg <- beta_bg[common_cpgs, , drop = FALSE]
  anno <- anno[match(common_cpgs, anno$CpG), ]
  
  stopifnot(identical(rownames(beta_bg), anno$CpG))
  
  cpg_gr <- GRanges(
    seqnames = anno$chr,
    ranges = IRanges(start = anno$pos, end = anno$pos),
    CpG = anno$CpG
  )
  
  hits <- findOverlaps(dmr_gr, cpg_gr)
  
  dmr_cpg_list <- split(
    mcols(cpg_gr)$CpG[subjectHits(hits)],
    queryHits(hits)
  )
  
  ############################################################
  ## 6.4 Calculate DMR-level beta in normal/background samples
  ############################################################
  
  dataset_summary_list <- list()
  tissue_mean_list <- list()
  
  for (i in seq_len(nrow(dmrs))) {
    
    cpgs <- dmr_cpg_list[[as.character(i)]]
    
    if (is.null(cpgs)) {
      
      dataset_summary_list[[i]] <- data.frame(
        DMR_ID = dmrs$DMR_ID[i],
        chr = dmrs$chr[i],
        start = dmrs$start[i],
        end = dmrs$end[i],
        dataset = dataset_name,
        nCpG_background = 0,
        mean_beta_all_background = NA_real_,
        max_sample_beta = NA_real_,
        max_group_mean_beta = NA_real_,
        background_pass = FALSE,
        stringsAsFactors = FALSE
      )
      
      next
    }
    
    cpgs <- unique(intersect(cpgs, rownames(beta_bg)))
    
    if (length(cpgs) < 1) {
      
      dataset_summary_list[[i]] <- data.frame(
        DMR_ID = dmrs$DMR_ID[i],
        chr = dmrs$chr[i],
        start = dmrs$start[i],
        end = dmrs$end[i],
        dataset = dataset_name,
        nCpG_background = 0,
        mean_beta_all_background = NA_real_,
        max_sample_beta = NA_real_,
        max_group_mean_beta = NA_real_,
        background_pass = FALSE,
        stringsAsFactors = FALSE
      )
      
      next
    }
    
    ## DMR-level beta for each sample
    dmr_value <- colMeans(
      beta_bg[cpgs, , drop = FALSE],
      na.rm = TRUE
    )
    
    ## Mean beta by normal tissue/group
    group_mean <- tapply(
      dmr_value,
      pd_bg$group,
      mean,
      na.rm = TRUE
    )
    
    group_mean_df <- data.frame(
      DMR_ID = dmrs$DMR_ID[i],
      chr = dmrs$chr[i],
      start = dmrs$start[i],
      end = dmrs$end[i],
      dataset = dataset_name,
      group = names(group_mean),
      mean_beta = as.numeric(group_mean),
      stringsAsFactors = FALSE
    )
    
    tissue_mean_list[[i]] <- group_mean_df
    
    max_group_mean_beta <- max(group_mean, na.rm = TRUE)
    
    ## Main background filtering criterion:
    ## 1. Enough CpG coverage: nCpG_background > 5
    ## 2. Low methylation in all normal/background tissue groups:
    ##    max_group_mean_beta < 0.20
    background_pass <- 
      length(cpgs) >= min_bg_cpg &&
      !is.na(max_group_mean_beta) &&
      max_group_mean_beta < beta_threshold
    
    dataset_summary_list[[i]] <- data.frame(
      DMR_ID = dmrs$DMR_ID[i],
      chr = dmrs$chr[i],
      start = dmrs$start[i],
      end = dmrs$end[i],
      dataset = dataset_name,
      nCpG_background = length(cpgs),
      mean_beta_all_background = mean(dmr_value, na.rm = TRUE),
      max_sample_beta = max(dmr_value, na.rm = TRUE),
      max_group_mean_beta = max_group_mean_beta,
      background_pass = background_pass,
      stringsAsFactors = FALSE
    )
  }
  
  dataset_summary <- rbindlist(dataset_summary_list, fill = TRUE)
  tissue_mean_df <- rbindlist(tissue_mean_list, fill = TRUE)
  
  return(
    list(
      dataset_summary = dataset_summary,
      tissue_mean_df = tissue_mean_df
    )
  )
}

############################################################
## 7. Background datasets
############################################################

background_files <- data.frame(
  dataset = c("GSE213478", "GSE48472"),
  file = c(
    "step0.1-GSE213478_output.Rdata",
    "step0.2-GSE48472_output.Rdata"
  ),
  arraytype = c("EPIC", "450K"),
  stringsAsFactors = FALSE
)

## Check whether files exist
print(file.exists(background_files$file))

if (!all(file.exists(background_files$file))) {
  stop("Some background .Rdata files do not exist. Please check file names.")
}

############################################################
## 8. Run background filtering
############################################################

background_outputs <- lapply(seq_len(nrow(background_files)), function(i) {
  
  background_filter_one_dataset(
    dataset_name = background_files$dataset[i],
    rdata_file = background_files$file[i],
    arraytype = background_files$arraytype[i]
  )
  
})

############################################################
## 9. Combine background filtering results
############################################################

background_summary_all <- rbindlist(
  lapply(background_outputs, function(x) x$dataset_summary),
  fill = TRUE
)

background_tissue_beta_all <- rbindlist(
  lapply(background_outputs, function(x) x$tissue_mean_df),
  fill = TRUE
)

write.csv(
  background_summary_all,
  file = "background_DMR_beta_summary_by_dataset.csv",
  row.names = FALSE
)

write.csv(
  background_tissue_beta_all,
  file = "background_DMR_beta_by_tissue_group.csv",
  row.names = FALSE
)

############################################################
## 10. Summarize background filtering across datasets
############################################################

safe_max <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  } else {
    return(max(x, na.rm = TRUE))
  }
}

safe_min <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  } else {
    return(min(x, na.rm = TRUE))
  }
}

background_cross_dataset_summary <- background_summary_all[
  ,
  .(
    n_background_datasets = .N,
    
    n_background_pass = sum(
      background_pass,
      na.rm = TRUE
    ),
    
    min_nCpG_background = safe_min(
      nCpG_background
    ),
    
    max_nCpG_background = safe_max(
      nCpG_background
    ),
    
    max_beta_across_background = safe_max(
      max_group_mean_beta
    ),
    
    mean_beta_across_background = mean(
      mean_beta_all_background,
      na.rm = TRUE
    ),
    
    background_pass_all = all(
      background_pass,
      na.rm = TRUE
    )
  ),
  by = .(
    DMR_ID,
    chr,
    start,
    end
  )
]

write.csv(
  background_cross_dataset_summary,
  file = "background_cross_dataset_summary.csv",
  row.names = FALSE
)

############################################################
## 11. Merge with original DMR information
############################################################

final_dmrs <- merge(
  as.data.table(dmrs),
  background_cross_dataset_summary,
  by = c("DMR_ID", "chr", "start", "end"),
  all.x = TRUE
)

############################################################
## 12. Keep final background-passed DMRs
############################################################
## Final filtering criterion:
##   1. This DMR passes background filtering in all background datasets.
##   2. The minimum CpG number across background datasets is >= 6.
##   3. The maximum normal/background tissue group beta is < 0.20.

background_filtered_dmrs <- final_dmrs[
  background_pass_all == TRUE &
    min_nCpG_background >= min_bg_cpg &
    max_beta_across_background < beta_threshold,
]

write.csv(
  background_filtered_dmrs,
  file = "background_filtered_DMRs_beta_lt_0.2_nCpG_gt5.csv",
  row.names = FALSE
)

############################################################
## 13. Save R object as .Rdata
############################################################

step2_background_output <- list(
  dmrs = dmrs,
  background_files = background_files,
  background_summary_all = background_summary_all,
  background_tissue_beta_all = background_tissue_beta_all,
  background_cross_dataset_summary = background_cross_dataset_summary,
  final_dmrs = final_dmrs,
  background_filtered_dmrs = background_filtered_dmrs,
  beta_threshold = beta_threshold,
  min_bg_cpg = min_bg_cpg
)

save(
  step2_background_output,
  file = "step2-background-filtering-output.Rdata"
)
## check
load("step2-background-filtering-output.Rdata")

dmrs <- step2_background_output$dmrs
background_summary_all <- step2_background_output$background_summary_all
background_tissue_beta_all <- step2_background_output$background_tissue_beta_all
background_cross_dataset_summary <- step2_background_output$background_cross_dataset_summary
background_filtered_dmrs <- step2_background_output$background_filtered_dmrs

beta_threshold <- step2_background_output$beta_threshold
min_bg_cpg <- step2_background_output$min_bg_cpg

dim(background_filtered_dmrs)

all(background_filtered_dmrs$background_pass_all == TRUE)
all(background_filtered_dmrs$min_nCpG_background >= min_bg_cpg)
all(background_filtered_dmrs$max_beta_across_background < beta_threshold)
#f <- read.csv("background_filtered_DMRs_beta_lt_0.2_nCpG_gt5.csv")
```


## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
