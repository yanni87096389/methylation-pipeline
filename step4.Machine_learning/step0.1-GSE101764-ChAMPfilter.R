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