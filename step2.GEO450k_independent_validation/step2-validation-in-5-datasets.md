---
title: "GEO450k_independent_validation"
output: html_document
date: "2026-06-25"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## step0.data preprocessing

### GSE42752
```{r GSE42752}
## GSE42752 450k methylation processed matrix pipeline
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
eset <- getGEO("GSE42752",destdir = './',AnnotGPL = T,getGPL = F)
beta.m <- exprs(eset[[1]])
class(beta.m)
dim(beta.m)
############################################################
## Module 2. Build phenotype table pD

## Build a clean phenotype table.
## Modify this part when applying the pipeline to another dataset
pD.all <- pData(eset[[1]])
pD <- pD.all[,c(37,2)]
names(pD) <- c("group", "ID")
pD$group <- ifelse(grepl("normal",pD$group,ignore.case = T),
  "Normal","Tumor")

#⭐Check if the extracted groupings are correct
check_df <- data.frame(Title = pD.all$`tissue:ch1`, Group = pD$group)
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

# Set pD row names as GSM IDs.
rownames(pD) <- pD$ID
## Final strict alignment check.
dim(beta.m)
dim(pD)

# Set group order
pD$group <- factor(
  pD$group,
  levels = c("Normal", "Tumor"))

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
save(myLoad, file = "step0.1-GSE42752_output.Rdata")
#saveRDS(myLoad, file = "step0.1-GSE42752_output.rds")
#myLoad <- readRDS("step0.1-GSE42752_output.rds")

# For all downstream analyses, simply use this myLoad variable
#champ.filter() performs initial quality control:
#(1).Filters out probes with low detection p-value.
#(2).Removes probes with missing values.
#(3).Filters sex chromosomes or SNP-associated probes if needed.
#(4).The output myLoad is an object used for all downstream ChAMP analyses 
##(normalization, DMP, DMR, pathway analysis).


```

### GSE48684
```{r GSE48684}
###############################################
## GSE48684 450k methylation processed matrix pipeline
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
eset <- getGEO("GSE48684",destdir = './',AnnotGPL = T,getGPL = F)
beta.m <- exprs(eset[[1]])
class(beta.m)
dim(beta.m)
############################################################
## Module 2. Build phenotype table pD

## Build a clean phenotype table.
## Modify this part when applying the pipeline to another dataset
pD.all <- pData(eset[[1]])
pD <- pD.all[,c(22,2)]
pD_nt <- pD[pD$description %in% c("normal-C", "normal-H", "cancer"), ]
pD_nt1 <- pD_nt
pD_nt1$Group <- ifelse(pD_nt1$description %in% c("normal-C", "normal-H"),
  "Normal",
  "Tumor")
pD_nt1 <- pD_nt1[, c("geo_accession", "Group")]
colnames(pD_nt1)[colnames(pD_nt1) == "geo_accession"] <- "ID"
colnames(pD_nt1)[colnames(pD_nt1) == "Group"] <- "group"
pD <- pD_nt1

# Check the final group distribution
table(pD$group)

#⭐Check if the extracted groupings are correct
check_df <- data.frame(Title = pD_nt$description, Group = pD$group)
head(check_df, 10)
tail(check_df,10)

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
dim(pD)
dim(beta.m)
############################################################
## Module 4. Run ChAMP filtering
############################################################

## Run ChAMP filtering using the cleaned beta matrix and aligned pD.
myLoad <- champ.filter(beta = beta.m, pd = pD, arraytype = "450k")
## Save the ChAMP object for downstream analysis.
save(myLoad, file = "step0.2-GSE48684_output.Rdata")
#saveRDS(myLoad, file = "step0.2-GSE48684_output.rds")
#myLoad <- readRDS("step0.1-GSE42752_output.rds")

```

### GSE77718

```{r GSE77718}
###############################################
## GSE77718 450k methylation processed matrix pipeline
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
eset <- getGEO("GSE77718",destdir = './',AnnotGPL = T,getGPL = F)
beta.m <- exprs(eset[[1]])
class(beta.m)
dim(beta.m)
############################################################
## Module 2. Build phenotype table pD

## Build a clean phenotype table.
## Modify this part when applying the pipeline to another dataset
pD.all <- pData(eset[[1]])
pD <- pD.all[,c(33,2)]
names(pD) <- c("group", "ID")
pD$group <- ifelse(grepl("Normal",pD$group,ignore.case = T),
                   "Normal","Tumor")

#⭐Check if the extracted groupings are correct
check_df <- data.frame(Title = pD.all$`disease state:ch1`, Group = pD$group)
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
save(myLoad, file = "step0.3-GSE77718_output.Rdata")
#saveRDS(myLoad, file = "step0.3-GSE77718_output.rds")
#myLoad <- readRDS("step0.1-GSE42752_output.rds")


```

### GSE101764
```{r GSE101764}
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
save(myLoad, file = "step0.4-GSE101764_output.Rdata")
#saveRDS(myLoad, file = "step0.4-GSE101764_output.rds")
#myLoad <- readRDS("step0.4-GSE101764_output.rds")
```

### GSE131013
```{r GSE131013}
###############################################
## GSE131013 450k methylation processed matrix pipeline
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
library(data.table)

############################################################
## Module 1. Read GEO metadata
#Ensure that the methylation beta matrix and the phenotype information strictly correspond to each other
## Download/read GEO series matrix metadata.
eset <- getGEO("GSE131013",destdir = './',AnnotGPL = T,getGPL = F)
#beta.m <- exprs(eset[[1]])
#class(beta.m)
#dim(beta.m)
############################################################
## Module 2. Build phenotype table pD

## Build a clean phenotype table.
## Modify this part when applying the pipeline to another dataset
pD.all <- pData(eset[[1]])
pD <- pD.all[,c(10,2)]
names(pD) <- c("group", "ID")
pD$group <- ifelse(grepl("Normal",pD$group,ignore.case = T),
                   "Normal","Tumor")

#⭐Check if the extracted groupings are correct
check_df <- data.frame(Title = pD.all$characteristics_ch1, Group = pD$group)
head(check_df, 10)
tail(check_df,10)
table(pD$group)

## prepare standard beta value matrix
beta.df <- fread("GSE131013_normalized_matrix.txt.gz", data.table = FALSE)
dim(beta.df)
head(rownames(beta.df))
head(colnames(beta.df))
head(beta.df[, 1:5])
#extract columns including beta value
beta.cols <- colnames(beta.df)[
  colnames(beta.df) != "ID_REF" &
    !grepl("Detection\\.Pval$", colnames(beta.df))
]

length(beta.cols)
head(beta.cols)
tail(beta.cols)
## Make sure GEO accession is available as a column
if (!"geo_accession" %in% colnames(pD.all)) {
  pD.all$geo_accession <- rownames(pD.all)
}

## Set row names of phenotype metadata as GEO accession numbers
rownames(pD.all) <- pD.all$geo_accession

## Extract the internal sample ID from the GEO sample title
## Example:
## "Paired Adjacent Normal sample from A2004_N patient" -> "A2004_N"
## "Tumor sample from A2004_T patient" -> "A2004_T"
## "Mucosa sample from A6144_M healthy donor" -> "A6144_M"
pD.all$matrix_id <- sub(
  ".*from[[:space:]]+([^[:space:]]+).*",
  "\\1",
  pD.all$title
)

## Check whether matrix_id extraction worked correctly
head(pD.all[, c("title", "geo_accession", "matrix_id")], 20)

## Check how many extracted matrix IDs are found in beta column names
sum(pD.all$matrix_id %in% beta.cols)

## Check unmatched sample IDs between GEO metadata and beta matrix
## Sample IDs present in GEO metadata but not directly found in beta matrix
setdiff(pD.all$matrix_id, beta.cols)
## Sample IDs present in beta matrix but not directly found in GEO metadata
setdiff(beta.cols, pD.all$matrix_id)
## Check matching after removing extra "_0" suffix
beta.cols.clean <- sub("_0$", "", beta.cols)

sum(pD.all$matrix_id %in% beta.cols.clean)

setdiff(pD.all$matrix_id, beta.cols.clean)
setdiff(beta.cols.clean, pD.all$matrix_id)

## Build a mapping table for beta matrix sample IDs
beta.map <- data.frame(
  matrix_id_original = beta.cols,
  matrix_id_for_match = beta.cols,
  stringsAsFactors = FALSE
)

## Some sample IDs may have an extra "_0" suffix in the beta matrix
## Example: Z2015_T_0 may need to match Z2015_T in GEO metadata
beta.map$matrix_id_for_match_no0 <- sub("_0$", "", beta.map$matrix_id_for_match)

## First try exact matching
idx <- match(beta.map$matrix_id_for_match, pD.all$matrix_id)

## For unmatched samples, try matching after removing the "_0" suffix
unmatched <- is.na(idx)

idx[unmatched] <- match(
  beta.map$matrix_id_for_match_no0[unmatched],
  pD.all$matrix_id
)
## Build the final sample mapping table
sample.map <- data.frame(
  matrix_id_original = beta.map$matrix_id_original,
  matrix_id_for_match = ifelse(
    is.na(match(beta.map$matrix_id_for_match, pD.all$matrix_id)),
    beta.map$matrix_id_for_match_no0,
    beta.map$matrix_id_for_match
  ),
  GEO_accession = pD.all$geo_accession[idx],
  title = pD.all$title[idx],
  stringsAsFactors = FALSE
)
## Set row names of the sample mapping table as original matrix IDs
rownames(sample.map) <- sample.map$matrix_id_original

## Check unmatched samples
sample.map[is.na(sample.map$GEO_accession), ]
## Check the mapping table
head(sample.map, 20)
tail(sample.map, 20)

## Build the standard beta value matrix
## Extract beta values
beta.m <- beta.df[, sample.map$matrix_id_original]

## Set row names as CpG probe IDs
rownames(beta.m) <- beta.df$ID_REF

## Convert data frame to numeric matrix
beta.m <- as.matrix(beta.m)
mode(beta.m) <- "numeric"

## Set column names as GEO accession numbers
colnames(beta.m) <- sample.map$GEO_accession

## Check the standard beta matrix
dim(beta.m)
head(rownames(beta.m))
head(colnames(beta.m))
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
save(myLoad, file = "step0.5-GSE131013_output.Rdata")
#saveRDS(myLoad, file = "step0.5-GSE131013_output.rds")
#myLoad <- readRDS("step0.5-GSE131013_output.rds")


```

## step1.validation-in-5-datasets

```{r validation in 5 datasets}
## step1.validation-in-5-datasets

```{r validation in 5 datasets}
############################################################
## Independent validation of discovery hyper-DMRs
## in 5 external 450K datasets
##
## Discovery DMRs: read from CSV
## Validation datasets: read from .Rdata files containing myLoad$beta and myLoad$pd
############################################################

rm(list = ls())
options(timeout = 100000)
options(scipen = 20)
options(stringsAsFactors = FALSE)

############################################################
## Step 0. Load required packages
############################################################

library(data.table)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

############################################################
## Step 1. Set input files
############################################################

## Discovery hyper-DMR file from DMRcate meta-analysis
discovery_dmr_file <- "dmr_hyper_dmrcate_meta.csv"

## Validation datasets
## Each .Rdata file should contain an object named myLoad
## with myLoad$beta and myLoad$pd
dataset_files <- c(
  GSE42752  = "step0.1-GSE42752_output.Rdata",
  GSE48684  = "step0.2-GSE48684_output.Rdata",
  GSE77718  = "step0.3-GSE77718_output.Rdata",
  GSE101764 = "step0.4-GSE101764_output.Rdata",
  GSE131013 = "step0.5-GSE131013_output.Rdata"
)

## Output directory
outdir <- "independent_validation_5datasets"

if (!dir.exists(outdir)) {
  dir.create(outdir)
}

############################################################
## Step 2. Validation criteria
############################################################

## Important:
## The code below uses nCpG > min_cpg_required.
## Therefore min_cpg_required <- 5 means at least 6 CpGs are required.
## This matches your discovery criterion: no.cpgs >= 6.
##
## deltaBeta_validation is Tumor - Normal.
## For hyper-DMR validation, require Tumor - Normal > 0.20.

min_cpg_required <- 5
delta_cutoff <- 0.20
p_cutoff <- 0.05
fdr_cutoff <- 0.05

############################################################
## Step 3. Check input files
############################################################

if (!file.exists(discovery_dmr_file)) {
  stop(paste0("Discovery DMR CSV file is missing: ", discovery_dmr_file))
}

missing_files <- dataset_files[!file.exists(dataset_files)]

if (length(missing_files) > 0) {
  print(missing_files)
  stop("Some validation dataset files are missing. Please check dataset_files.")
}

cat("All input files exist.\n")

############################################################
## Step 4. Load discovery DMRs from CSV
############################################################

dmr_obj <- read.csv(
  discovery_dmr_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

## If the first column is an exported rowname column, use it as row names
if (ncol(dmr_obj) > 0 && colnames(dmr_obj)[1] %in% c("", "X", "V1", "...1")) {
  rownames(dmr_obj) <- dmr_obj[[1]]
  dmr_obj <- dmr_obj[, -1, drop = FALSE]
}

cat("\nDiscovery DMRs loaded from CSV:\n")
print(discovery_dmr_file)
print(dim(dmr_obj))
print(head(dmr_obj))

############################################################
## Step 5. Standardize discovery DMR table
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

find_cpg_col <- function(df) {
  
  candidate_cols <- colnames(df)[
    grepl("cpg|probe|probes|name|id", colnames(df), ignore.case = TRUE)
  ]
  
  if (length(candidate_cols) == 0) {
    return(NA_character_)
  }
  
  for (cn in candidate_cols) {
    
    vals <- as.character(df[[cn]])
    vals <- vals[!is.na(vals)]
    vals <- vals[seq_len(min(length(vals), 100))]
    
    if (any(grepl("cg[0-9]{8}", vals))) {
      return(cn)
    }
  }
  
  return(NA_character_)
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
  
  out <- data.frame(
    chr = standardize_chr(chr),
    start = start,
    end = end,
    stringsAsFactors = FALSE
  )
  
  return(out)
}

standardize_dmr_table <- function(dmr_obj) {
  
  dmr.df <- as.data.frame(dmr_obj)
  
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
  
  coord_col <- find_coord_col(dmr.df)
  cpg_col <- find_cpg_col(dmr.df)
  
  rowname_coord <- NULL
  
  if (!is.null(rownames(dmr.df))) {
    rn <- rownames(dmr.df)
    
    if (any(grepl("^chr[^:]+:[0-9,]+-[0-9,]+$", rn))) {
      rowname_coord <- parse_coord(rn)
    }
  }
  
  coord_parsed <- NULL
  
  if (!is.na(coord_col)) {
    coord_vals <- as.character(dmr.df[[coord_col]])
    
    if (any(grepl("^chr[^:]+:[0-9,]+-[0-9,]+$", coord_vals))) {
      coord_parsed <- parse_coord(coord_vals)
    }
  }
  
  ## chr
  if (!is.na(chr_col)) {
    chr_out <- standardize_chr(dmr.df[[chr_col]])
  } else if (!is.null(coord_parsed)) {
    chr_out <- coord_parsed$chr
  } else if (!is.null(rowname_coord)) {
    chr_out <- rowname_coord$chr
  } else {
    chr_out <- rep(NA_character_, nrow(dmr.df))
  }
  
  ## start
  if (!is.na(start_col)) {
    start_out <- as.integer(gsub(",", "", dmr.df[[start_col]]))
  } else if (!is.null(coord_parsed)) {
    start_out <- coord_parsed$start
  } else if (!is.null(rowname_coord)) {
    start_out <- rowname_coord$start
  } else {
    start_out <- rep(NA_integer_, nrow(dmr.df))
  }
  
  ## end
  if (!is.na(end_col)) {
    end_out <- as.integer(gsub(",", "", dmr.df[[end_col]]))
  } else if (!is.null(coord_parsed)) {
    end_out <- coord_parsed$end
  } else if (!is.null(rowname_coord)) {
    end_out <- rowname_coord$end
  } else {
    end_out <- rep(NA_integer_, nrow(dmr.df))
  }
  
  ## DMR ID
  if (!all(is.na(chr_out)) && !all(is.na(start_out)) && !all(is.na(end_out))) {
    dmr_id <- paste0(chr_out, ":", start_out, "-", end_out)
  } else if (!is.na(coord_col)) {
    dmr_id <- as.character(dmr.df[[coord_col]])
  } else if (!is.null(rownames(dmr.df)) &&
             !all(rownames(dmr.df) == as.character(seq_len(nrow(dmr.df))))) {
    dmr_id <- rownames(dmr.df)
  } else {
    dmr_id <- paste0("DMR_", seq_len(nrow(dmr.df)))
  }
  
  dmr_id <- make.unique(as.character(dmr_id))
  
  out <- data.frame(
    DMR_id = dmr_id,
    chr = chr_out,
    start = start_out,
    end = end_out,
    cpg_string = if (!is.na(cpg_col)) as.character(dmr.df[[cpg_col]]) else NA_character_,
    stringsAsFactors = FALSE
  )
  
  if (all(is.na(out$chr)) && all(is.na(out$cpg_string))) {
    stop("No chromosome coordinates or CpG probe IDs were found in the discovery DMR table.")
  }
  
  return(out)
}

dmr.tbl <- standardize_dmr_table(dmr_obj)

cat("\nStandardized DMR table:\n")
print(dim(dmr.tbl))
print(head(dmr.tbl))

cat("\nMissing coordinate counts:\n")
print(colSums(is.na(dmr.tbl[, c("chr", "start", "end")])))

############################################################
## Step 6. Load 450K annotation and map CpGs to discovery DMRs
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
    
    cpgs <- character(0)
    
    ## First, try to extract CpG IDs directly from the discovery DMR table
    if (!is.na(dmr.tbl$cpg_string[i]) &&
        grepl("cg[0-9]{8}", dmr.tbl$cpg_string[i])) {
      
      cpgs <- unique(unlist(regmatches(
        dmr.tbl$cpg_string[i],
        gregexpr("cg[0-9]{8}", dmr.tbl$cpg_string[i])
      )))
    }
    
    ## If no CpG IDs are available, map CpGs by chromosome coordinates
    if (length(cpgs) == 0 &&
        !is.na(dmr.tbl$chr[i]) &&
        !is.na(dmr.tbl$start[i]) &&
        !is.na(dmr.tbl$end[i])) {
      
      chr_i <- dmr.tbl$chr[i]
      start_i <- dmr.tbl$start[i]
      end_i <- dmr.tbl$end[i]
      
      cpgs <- anno.dt[
        chr == chr_i & pos >= start_i & pos <= end_i,
        CpG
      ]
    }
    
    cpg.list[[i]] <- unique(cpgs)
  }
  
  return(cpg.list)
}

dmr.cpg.list <- build_dmr_cpg_list(dmr.tbl, anno.dt)

dmr.tbl$nCpG_450K_annotation <- lengths(dmr.cpg.list)

cat("\nSummary of CpGs per DMR based on 450K annotation:\n")
print(summary(dmr.tbl$nCpG_450K_annotation))

cat("\nNumber of DMRs with 0 mapped CpGs:\n")
print(sum(dmr.tbl$nCpG_450K_annotation == 0))

############################################################
## Step 7. Helper functions for validation datasets
############################################################

load_validation_dataset <- function(file) {
  
  env <- new.env()
  loaded_objects <- load(file, envir = env)
  
  if (!"myLoad" %in% loaded_objects) {
    stop(paste0(
      "Object 'myLoad' was not found in file: ",
      file,
      ". Objects found: ",
      paste(loaded_objects, collapse = ", ")
    ))
  }
  
  obj <- get("myLoad", envir = env)
  
  if (!is.list(obj)) {
    stop(paste0("myLoad is not a list in file: ", file))
  }
  
  if (!"beta" %in% names(obj)) {
    stop(paste0("Cannot find myLoad$beta in file: ", file))
  }
  
  if (!"pd" %in% names(obj)) {
    stop(paste0("Cannot find myLoad$pd in file: ", file))
  }
  
  beta <- obj$beta
  pd <- obj$pd
  
  beta <- as.matrix(beta)
  mode(beta) <- "numeric"
  
  pd <- as.data.frame(pd)
  
  return(list(beta = beta, pd = pd))
}

normalize_group <- function(group_vector) {
  
  x <- as.character(group_vector)
  out <- rep(NA_character_, length(x))
  
  out[grepl("tumor|tumour|cancer|carcinoma|crc", x, ignore.case = TRUE)] <- "Tumor"
  out[grepl("normal|control|adjacent|mucosa|healthy", x, ignore.case = TRUE)] <- "Normal"
  
  return(factor(out, levels = c("Normal", "Tumor")))
}

prepare_validation_dataset <- function(dataset_name, file) {
  
  dat <- load_validation_dataset(file)
  
  beta <- dat$beta
  pd <- dat$pd
  
  ## Set row names of phenotype table using ID column if possible
  if ("ID" %in% colnames(pd) && setequal(as.character(pd$ID), colnames(beta))) {
    rownames(pd) <- as.character(pd$ID)
  }
  
  ## If sample names are the same but order is different, reorder phenotype table
  if (setequal(rownames(pd), colnames(beta))) {
    pd <- pd[colnames(beta), , drop = FALSE]
  }
  
  ## Stop if beta matrix and phenotype table cannot be aligned
  if (!identical(rownames(pd), colnames(beta))) {
    stop(paste0(
      "Sample alignment failed in ",
      dataset_name,
      ". Please check rownames(pd), pd$ID, and colnames(beta)."
    ))
  }
  
  ## Find group column
  group_col <- NA_character_
  
  for (cn in c("group", "Group", "condition", "Condition", "phenotype", "Phenotype")) {
    if (cn %in% colnames(pd)) {
      group_col <- cn
      break
    }
  }
  
  if (is.na(group_col)) {
    stop(paste0("No group column was found in phenotype table of ", dataset_name))
  }
  
  ## Standardize group labels to Normal and Tumor
  pd$group <- normalize_group(pd[[group_col]])
  
  ## Keep only Normal and Tumor samples
  keep <- !is.na(pd$group) & pd$group %in% c("Normal", "Tumor")
  
  beta <- beta[, keep, drop = FALSE]
  pd <- pd[keep, , drop = FALSE]
  
  ## Set clean ID column
  pd$ID <- rownames(pd)
  
  ## Check sample size
  if (sum(pd$group == "Normal") < 2 || sum(pd$group == "Tumor") < 2) {
    stop(paste0(
      "Too few samples after filtering in ",
      dataset_name,
      ". Need at least 2 Normal and 2 Tumor samples."
    ))
  }
  
  ## Sort samples by group
  sample_order <- order(pd$group)
  pd <- pd[sample_order, , drop = FALSE]
  beta <- beta[, rownames(pd), drop = FALSE]
  
  ## Final strict check
  stopifnot(identical(rownames(pd), colnames(beta)))
  
  cat("\nDataset:", dataset_name, "\n")
  cat("beta dimension:\n")
  print(dim(beta))
  cat("group distribution:\n")
  print(table(pd$group, useNA = "ifany"))
  cat("alignment check:\n")
  print(identical(rownames(pd), colnames(beta)))
  
  return(list(beta = beta, pd = pd))
}

############################################################
## Step 8. Validate discovery DMRs in one dataset
############################################################

validate_one_dataset <- function(dataset_name, file, dmr.tbl, dmr.cpg.list) {
  
  dat <- prepare_validation_dataset(dataset_name, file)
  
  beta <- dat$beta
  pd <- dat$pd
  
  normal_samples <- rownames(pd)[pd$group == "Normal"]
  tumor_samples <- rownames(pd)[pd$group == "Tumor"]
  
  res.list <- vector("list", nrow(dmr.tbl))
  
  for (i in seq_len(nrow(dmr.tbl))) {
    
    dmr_id <- dmr.tbl$DMR_id[i]
    
    ## Use list index instead of DMR_id name to avoid duplicated-name problems
    cpgs_all <- dmr.cpg.list[[i]]
    cpgs_present <- intersect(cpgs_all, rownames(beta))
    
    nCpG <- length(cpgs_present)
    
    if (nCpG <= min_cpg_required) {
      
      res.list[[i]] <- data.frame(
        dataset = dataset_name,
        DMR_id = dmr_id,
        chr = dmr.tbl$chr[i],
        start = dmr.tbl$start[i],
        end = dmr.tbl$end[i],
        nCpG = nCpG,
        mean_Normal = NA_real_,
        mean_Tumor = NA_real_,
        deltaBeta_validation = NA_real_,
        p.value = NA_real_,
        stringsAsFactors = FALSE
      )
      
      next
    }
    
    ## Calculate mean beta value of this DMR for each sample
    dmr_beta <- colMeans(beta[cpgs_present, , drop = FALSE], na.rm = TRUE)
    dmr_beta[is.nan(dmr_beta)] <- NA_real_
    
    normal_values <- dmr_beta[normal_samples]
    tumor_values <- dmr_beta[tumor_samples]
    
    mean_normal <- mean(normal_values, na.rm = TRUE)
    mean_tumor <- mean(tumor_values, na.rm = TRUE)
    
    delta_beta <- mean_tumor - mean_normal
    
    pval <- tryCatch(
      t.test(tumor_values, normal_values)$p.value,
      error = function(e) NA_real_
    )
    
    res.list[[i]] <- data.frame(
      dataset = dataset_name,
      DMR_id = dmr_id,
      chr = dmr.tbl$chr[i],
      start = dmr.tbl$start[i],
      end = dmr.tbl$end[i],
      nCpG = nCpG,
      mean_Normal = mean_normal,
      mean_Tumor = mean_tumor,
      deltaBeta_validation = delta_beta,
      p.value = pval,
      stringsAsFactors = FALSE
    )
  }
  
  res <- rbindlist(res.list, fill = TRUE)
  
  ## Adjust p values within each validation dataset
  res$FDR <- p.adjust(res$p.value, method = "BH")
  
  ## Main validation criterion for discovery hyper-DMRs
  res$validated_p005 <- with(
    res,
    nCpG > min_cpg_required &
      deltaBeta_validation > delta_cutoff &
      p.value < p_cutoff
  )
  
  ## More stringent FDR-based criterion
  res$validated_FDR005 <- with(
    res,
    nCpG > min_cpg_required &
      deltaBeta_validation > delta_cutoff &
      FDR < fdr_cutoff
  )
  
  return(res)
}

############################################################
## Step 9. Run independent validation across all datasets
############################################################

all.validation.list <- list()

for (dataset_name in names(dataset_files)) {
  
  cat("\n============================================================\n")
  cat("Running validation for:", dataset_name, "\n")
  cat("============================================================\n")
  
  all.validation.list[[dataset_name]] <- validate_one_dataset(
    dataset_name = dataset_name,
    file = dataset_files[[dataset_name]],
    dmr.tbl = dmr.tbl,
    dmr.cpg.list = dmr.cpg.list
  )
}

all.validation <- rbindlist(all.validation.list, fill = TRUE)

cat("\nAll validation result dimension:\n")
print(dim(all.validation))
print(head(all.validation))

############################################################
## Step 10. Summarize validation results across datasets
############################################################

safe_mean <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  return(mean(x))
}

safe_max <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  return(max(x))
}

safe_min <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  return(min(x))
}

summarize_one_dmr <- function(df) {
  
  pvals <- df$p.value
  pvals <- pvals[!is.na(pvals)]
  
  fisher_p <- NA_real_
  
  if (length(pvals) > 0) {
    pvals <- pmax(pvals, .Machine$double.xmin)
    fisher_p <- pchisq(
      -2 * sum(log(pvals)),
      df = 2 * length(pvals),
      lower.tail = FALSE
    )
  }
  
  data.frame(
    DMR_id = df$DMR_id[1],
    chr = df$chr[1],
    start = df$start[1],
    end = df$end[1],
    tested_datasets = sum(!is.na(df$p.value)),
    validated_datasets_p005 = sum(df$validated_p005, na.rm = TRUE),
    validated_datasets_FDR005 = sum(df$validated_FDR005, na.rm = TRUE),
    mean_deltaBeta = safe_mean(df$deltaBeta_validation),
    max_deltaBeta = safe_max(df$deltaBeta_validation),
    min_p.value = safe_min(pvals),
    fisher_p.value = fisher_p,
    stringsAsFactors = FALSE
  )
}

summary.list <- lapply(
  split(all.validation, all.validation$DMR_id),
  summarize_one_dmr
)

validation.summary <- rbindlist(summary.list, fill = TRUE)

validation.summary$meta_FDR <- p.adjust(
  validation.summary$fisher_p.value,
  method = "BH"
)

validation.summary$validated_any_p005 <- validation.summary$validated_datasets_p005 >= 1
validation.summary$validated_at_least2_p005 <- validation.summary$validated_datasets_p005 >= 2
validation.summary$validated_any_FDR005 <- validation.summary$validated_datasets_FDR005 >= 1

## Sort DMRs by validation strength
validation.summary <- validation.summary[
  order(
    -validated_datasets_p005,
    -mean_deltaBeta,
    min_p.value
  ),
]

cat("\nValidation summary dimension:\n")
print(dim(validation.summary))
print(head(validation.summary))

############################################################
## Step 11. Save validation results
############################################################

## Save all per-dataset results
fwrite(
  all.validation,
  file = file.path(outdir, "validation_results_all_DMRs_all_5datasets.csv")
)

## Save summary results across all validation datasets
fwrite(
  validation.summary,
  file = file.path(outdir, "validation_summary_across_5datasets.csv")
)

## Save DMRs validated in at least one dataset using p < 0.05 criterion
validated_p005 <- validation.summary[
  validation.summary$validated_any_p005 == TRUE,
]

fwrite(
  validated_p005,
  file = file.path(outdir, "validated_hyper_DMRs_at_least_one_dataset_p005.csv")
)

## Save DMRs validated in at least two datasets using p < 0.05 criterion
validated_at_least2_p005 <- validation.summary[
  validation.summary$validated_at_least2_p005 == TRUE,
]

fwrite(
  validated_at_least2_p005,
  file = file.path(outdir, "validated_hyper_DMRs_at_least_two_datasets_p005.csv")
)

## Save DMRs validated in at least one dataset using FDR < 0.05 criterion
validated_FDR005 <- validation.summary[
  validation.summary$validated_any_FDR005 == TRUE,
]

fwrite(
  validated_FDR005,
  file = file.path(outdir, "validated_hyper_DMRs_at_least_one_dataset_FDR005.csv")
)

## Save per-dataset results separately
for (dataset_name in names(dataset_files)) {
  
  tmp <- all.validation[all.validation$dataset == dataset_name, ]
  
  fwrite(
    tmp,
    file = file.path(
      outdir,
      paste0(dataset_name, "_validation_results_all_DMRs.csv")
    )
  )
}

############################################################
## Step 12. Dataset-level summary and barplot
############################################################

dataset.summary <- rbindlist(lapply(
  split(all.validation, all.validation$dataset),
  function(df) {
    data.frame(
      dataset = df$dataset[1],
      tested_DMRs = sum(!is.na(df$p.value)),
      validated_p005 = sum(df$validated_p005, na.rm = TRUE),
      validated_FDR005 = sum(df$validated_FDR005, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))

fwrite(
  dataset.summary,
  file = file.path(outdir, "validation_summary_by_dataset.csv")
)

png(
  filename = file.path(outdir, "validation_summary_barplot.png"),
  width = 1600,
  height = 900,
  res = 150
)

barplot(
  dataset.summary$validated_p005,
  names.arg = dataset.summary$dataset,
  las = 2,
  ylab = "Number of validated hyper-DMRs",
  main = "Independent validation of discovery hyper-DMRs"
)

dev.off()

cat("\nDataset-level validation summary:\n")
print(dataset.summary)

############################################################
## Step 13. Export robustly validated hyper-DMRs
############################################################

robust_validated_DMRs <- validation.summary[
  validation.summary$validated_datasets_p005 >= 3,
]

high_confidence_DMRs <- validation.summary[
  validation.summary$validated_datasets_p005 >= 4,
]

all5_validated_DMRs <- validation.summary[
  validation.summary$validated_datasets_p005 == 5,
]

cat("\nRobust validation results:\n")
cat("Validated in at least 3 datasets:\n")
print(dim(robust_validated_DMRs))

cat("Validated in at least 4 datasets:\n")
print(dim(high_confidence_DMRs))

cat("Validated in all 5 datasets:\n")
print(dim(all5_validated_DMRs))

fwrite(
  robust_validated_DMRs,
  file = file.path(outdir, "robust_validated_hyper_DMRs_at_least_3datasets.csv")
)

fwrite(
  high_confidence_DMRs,
  file = file.path(outdir, "high_confidence_hyper_DMRs_at_least_4datasets.csv")
)

fwrite(
  all5_validated_DMRs,
  file = file.path(outdir, "hyper_DMRs_validated_in_all_5datasets.csv")
)

############################################################
## Step 14. Final objects to inspect
############################################################

cat("\nFinished independent validation.\n")
cat("Main output directory:\n")
print(outdir)

cat("\nMain files to check:\n")
print(file.path(outdir, "validation_summary_across_5datasets.csv"))
print(file.path(outdir, "validation_summary_by_dataset.csv"))
print(file.path(outdir, "robust_validated_hyper_DMRs_at_least_3datasets.csv"))
print(file.path(outdir, "high_confidence_hyper_DMRs_at_least_4datasets.csv"))

head(validation.summary)
dataset.summary

```



## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
