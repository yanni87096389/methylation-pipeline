###############################################
## GSE196696 EPIC methylation processed matrix pipeline
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
## Module 1. Read GEO metadata and processed beta matrix
############################################################

## Download/read GEO series matrix metadata.
eset <- getGEO("GSE196696",destdir = './',AnnotGPL = T,getGPL = F)
#beta.m <- eset[[1]]
#dim(beta.m)
pD.all <- pData(eset[[1]])
table(pD.all$`tissue:ch1`)
which(colnames(pD.all)=="tissue:ch1")
pD <- pD.all[,c(39,2)]
names(pD) <- c("group", "ID")

############################################################
## Read processed methylation matrix
############################################################
## Make sure this file is in your current working directory:
## GSE196696_processed_data.tsv.gz
#nohup wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196696/suppl/GSE196696_processed_data.tsv.gz > download_GSE196696.log 2>&1 &

file_196696 <- "GSE196696_processed_data.tsv.gz"

## Check file format first
test_196696 <- fread(
  file_196696,
  nrows = 5,
  data.table = FALSE,
  check.names = FALSE
)

dim(test_196696)
colnames(test_196696)[1:10]
test_196696[1:5, 1:5]

## Read full processed data
raw_196696 <- fread(
  file_196696,
  data.table = FALSE,
  check.names = FALSE
)
raw_196696[1:5, 1:5]
class(raw_196696)
dim(raw_196696)
colnames(raw_196696)[1:10]

############################################################
## Extract beta value columns
############################################################

raw_df <- as.data.frame(raw_196696)
pd_df  <- as.data.frame(pD.all)

## Extract array ID from title
## This is used to match processed matrix columns to GSM IDs.
pd_df$array_id <- sub(".* -\\s*", "", pd_df$title)

## Check extracted array IDs
head(pd_df[, c("geo_accession", "title", "array_id")])

## Keep beta value columns only
## Remove ID_REF and Detection Pval columns.
beta_cols <- colnames(raw_df)[
  colnames(raw_df) != "ID_REF" &
    !grepl("^Detection Pval", colnames(raw_df))
]

length(beta_cols)
head(beta_cols)

############################################################
## Map processed matrix columns to GSM IDs
############################################################

## Build array_id to GSM mapping
array_to_gsm <- setNames(
  pd_df$geo_accession,
  pd_df$array_id
)

## Check unmatched beta columns
unmatched <- setdiff(beta_cols, names(array_to_gsm))

if (length(unmatched) > 0) {
  print(head(unmatched))
  stop("Some beta columns cannot be matched to GSM IDs. Please check title format or column names.")
}

## Convert beta matrix column names from array IDs to GSM IDs
gsm_names <- unname(array_to_gsm[beta_cols])

## Check mapping
head(
  data.frame(
    raw_column = beta_cols,
    GSM = gsm_names
  )
)

############################################################
## Build standard beta value matrix
############################################################

beta_df <- raw_df[, c("ID_REF", beta_cols), drop = FALSE]

## Check duplicated probes
if (anyDuplicated(beta_df$ID_REF) > 0) {
  stop("Duplicated ID_REF found. Please handle duplicated probes first.")
}

## Build beta matrix
beta.m <- beta_df[, beta_cols, drop = FALSE]

rownames(beta.m) <- beta_df$ID_REF
colnames(beta.m) <- gsm_names

## Convert to numeric matrix
beta.m <- as.matrix(beta.m)
storage.mode(beta.m) <- "numeric"

## Basic checks
class(beta.m)
dim(beta.m)
head(rownames(beta.m))
head(colnames(beta.m))

sum(is.na(colnames(beta.m)))
sum(is.na(rownames(beta.m)))
anyDuplicated(colnames(beta.m))
anyDuplicated(rownames(beta.m))

############################################################
## Module 2. Build phenotype table pD
############################################################

## GSE196696 is whole blood background dataset.
## For background filtering, all samples are treated as normal whole blood.

pD <- data.frame(
  group = "Whole_blood",
  ID = pd_df$geo_accession,
  stringsAsFactors = FALSE
)

## Check pD
head(pD)
tail(pD)
table(pD$group)
dim(pD)

## Check if all beta matrix samples are in pD
setequal(colnames(beta.m), pD$ID)

############################################################
## Module 3. Align pD with beta.m
############################################################

## Confirm that there are no duplicated sample names
sum(duplicated(colnames(beta.m)))
sum(duplicated(pD$ID))

## Set pD row names as GSM IDs
rownames(pD) <- pD$ID

## Align pD according to beta.m columns
pD <- pD[match(colnames(beta.m), pD$ID), , drop = FALSE]

## Final strict alignment check before sorting
stopifnot(identical(colnames(beta.m), rownames(pD)))

dim(beta.m)
dim(pD)

## Set group as factor
pD$group <- factor(pD$group)

## Sort pD by group
sample_order <- order(pD$group)

pD <- pD[sample_order, , drop = FALSE]

## Reorder beta matrix columns according to pD
beta.m <- beta.m[, rownames(pD), drop = FALSE]

## Final strict alignment check before champ.filter
identical(colnames(beta.m), rownames(pD))

table(pD$group)
dim(beta.m)
dim(pD)

############################################################
## Module 4. Run ChAMP filtering
############################################################

## Run ChAMP filtering using the cleaned beta matrix and aligned pD.
## GSE196696 is EPIC/850K array.
myLoad <- champ.filter(
  beta = beta.m,
  pd = pD,
  arraytype = "EPIC"
)
##  Save the ChAMP object for downstream analysis
save(myLoad,file = "step0.3-GSE196696_output.Rdata")
