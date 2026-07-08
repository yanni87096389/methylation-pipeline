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
