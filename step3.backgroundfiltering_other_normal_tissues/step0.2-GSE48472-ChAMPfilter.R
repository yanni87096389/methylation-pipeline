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

