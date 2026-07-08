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