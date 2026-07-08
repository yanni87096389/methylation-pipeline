---
title: "Background_filtering_in_blood_cells"
output: html_document
date: "2026-06-27"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## data preprocessing
### GSE196696 EPIC methylation
```{r}
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

```

## background filtering in whole blood

```{r whole blood filtering}
############################################################
## Step 3.1 Background filtering in GSE196696 whole blood cells
##
## Input:
##   1. background_filtered_DMRs_beta_lt_0.2_nCpG_gt5.csv
##   2. step0.3-GSE196696_output.Rdata
##
## Goal:
##   Further filter candidate DMRs using GSE196696 whole blood cells
##
## Filtering criterion:
##   nCpG_GSE196696 > 5
##   max_group_mean_beta_GSE196696 < 0.20
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
options(scipen = 20)

############################################################
## Module 0. Load required packages
############################################################

library(data.table)
library(GenomicRanges)
library(IRanges)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(ggplot2)

############################################################
## Module 1. Set filtering parameters
############################################################

beta_threshold <- 0.20

## Keep the same CpG coverage criterion as before:
## nCpG > 5, equivalent to nCpG >= 6
min_bg_cpg <- 6

dataset_name <- "GSE196696"

############################################################
## Module 2. Load candidate DMRs from previous background filtering
############################################################

candidate_dmrs <- fread(
  "background_filtered_DMRs_beta_lt_0.2_nCpG_gt5.csv",
  data.table = FALSE
)

dim(candidate_dmrs)
head(candidate_dmrs)

## Check required columns
stopifnot(all(c("DMR_ID", "chr", "start", "end") %in% colnames(candidate_dmrs)))
cat("Input candidate DMRs:", nrow(candidate_dmrs), "\n")

############################################################
## Module 3. Standardize genomic coordinates
############################################################

fix_chr <- function(x) {
  x <- as.character(x)
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  x[x == "23"] <- "X"
  x[x == "24"] <- "Y"
  x[x %in% c("M", "MT")] <- "M"
  paste0("chr", x)
}

candidate_dmrs$chr <- fix_chr(candidate_dmrs$chr)
candidate_dmrs$start <- as.numeric(candidate_dmrs$start)
candidate_dmrs$end <- as.numeric(candidate_dmrs$end)

candidate_dmrs <- candidate_dmrs[
  !is.na(candidate_dmrs$chr) &
    !is.na(candidate_dmrs$start) &
    !is.na(candidate_dmrs$end),
]

candidate_dmrs <- candidate_dmrs[
  candidate_dmrs$chr %in% paste0("chr", 1:22),
]

dim(candidate_dmrs)

dmr_gr <- GRanges(
  seqnames = candidate_dmrs$chr,
  ranges = IRanges(
    start = candidate_dmrs$start,
    end = candidate_dmrs$end
  ),
  DMR_ID = candidate_dmrs$DMR_ID
)

############################################################
## Module 4. Load GSE196696 ChAMP-filtered object
############################################################

load("step0.3-GSE196696_output.Rdata")

## The object should be named myLoad
beta_bg <- myLoad$beta
pd_bg <- myLoad$pd

dim(beta_bg)
dim(pd_bg)

## Check sample alignment
stopifnot(identical(colnames(beta_bg), rownames(pd_bg)))
pd_bg$group <- as.character(pd_bg$group)
table(pd_bg$group)

############################################################
## Module 5. Load EPIC annotation and map CpGs
############################################################

anno <- minfi::getAnnotation(
  IlluminaHumanMethylationEPICanno.ilm10b4.hg19
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

## Keep CpGs available in GSE196696 beta matrix
anno <- anno[anno$CpG %in% rownames(beta_bg), ]

common_cpgs <- intersect(rownames(beta_bg), anno$CpG)

beta_bg <- beta_bg[common_cpgs, , drop = FALSE]
anno <- anno[match(common_cpgs, anno$CpG), ]

stopifnot(identical(rownames(beta_bg), anno$CpG))

dim(beta_bg)
dim(anno)

cpg_gr <- GRanges(
  seqnames = anno$chr,
  ranges = IRanges(
    start = anno$pos,
    end = anno$pos
  ),
  CpG = anno$CpG
)

############################################################
## Module 6. Find CpGs located within each DMR
############################################################

hits <- findOverlaps(dmr_gr, cpg_gr)
cat("Number of DMR-CpG overlaps in GSE196696:", length(hits), "\n")
dmr_cpg_list <- split(
  mcols(cpg_gr)$CpG[subjectHits(hits)],
  queryHits(hits)
)

############################################################
## Module 7. Calculate DMR-level beta in GSE196696
############################################################

dmr_beta_GSE196696 <- matrix(
  NA_real_,
  nrow = nrow(candidate_dmrs),
  ncol = ncol(beta_bg),
  dimnames = list(candidate_dmrs$DMR_ID, colnames(beta_bg))
)

background_results_GSE196696 <- candidate_dmrs

background_results_GSE196696$nCpG_GSE196696 <- 0
background_results_GSE196696$meanBeta_GSE196696 <- NA_real_
background_results_GSE196696$medianBeta_GSE196696 <- NA_real_
background_results_GSE196696$maxSampleBeta_GSE196696 <- NA_real_
background_results_GSE196696$p95Beta_GSE196696 <- NA_real_
background_results_GSE196696$max_group_mean_beta_GSE196696 <- NA_real_
background_results_GSE196696$propSamples_beta_gt_0.20_GSE196696 <- NA_real_
background_results_GSE196696$GSE196696_pass <- FALSE

for (i in seq_len(nrow(candidate_dmrs))) {
  
  cpgs <- dmr_cpg_list[[as.character(i)]]
  
  if (is.null(cpgs)) {
    next
  }
  
  cpgs <- unique(intersect(cpgs, rownames(beta_bg)))
  
  background_results_GSE196696$nCpG_GSE196696[i] <- length(cpgs)
  
  if (length(cpgs) < 1) {
    next
  }
  
  ## DMR-level beta for each sample
  dmr_value <- colMeans(
    beta_bg[cpgs, , drop = FALSE],
    na.rm = TRUE
  )
  
  dmr_beta_GSE196696[i, ] <- dmr_value
  
  ## Group-level mean beta
  ## For this dataset, group should be Whole_blood
  group_mean <- tapply(
    dmr_value,
    pd_bg$group,
    mean,
    na.rm = TRUE
  )
  
  max_group_mean_beta <- max(group_mean, na.rm = TRUE)
  
  background_results_GSE196696$meanBeta_GSE196696[i] <- mean(
    dmr_value,
    na.rm = TRUE
  )
  
  background_results_GSE196696$medianBeta_GSE196696[i] <- median(
    dmr_value,
    na.rm = TRUE
  )
  
  background_results_GSE196696$maxSampleBeta_GSE196696[i] <- max(
    dmr_value,
    na.rm = TRUE
  )
  
  background_results_GSE196696$p95Beta_GSE196696[i] <- as.numeric(
    quantile(
      dmr_value,
      probs = 0.95,
      na.rm = TRUE
    )
  )
  
  background_results_GSE196696$max_group_mean_beta_GSE196696[i] <-
    max_group_mean_beta
  
  background_results_GSE196696$propSamples_beta_gt_0.20_GSE196696[i] <- mean(
    dmr_value > beta_threshold,
    na.rm = TRUE
  )
  
  ## Final GSE196696 background filtering criterion
  background_results_GSE196696$GSE196696_pass[i] <-
    length(cpgs) >= min_bg_cpg &&
    !is.na(max_group_mean_beta) &&
    max_group_mean_beta < beta_threshold
}

############################################################
## Module 8. Filter final DMRs after GSE196696
############################################################

final_dmrs_after_GSE196696 <- background_results_GSE196696[
  background_results_GSE196696$GSE196696_pass == TRUE &
    background_results_GSE196696$nCpG_GSE196696 >= min_bg_cpg &
    background_results_GSE196696$max_group_mean_beta_GSE196696 < beta_threshold,
]

removed_dmrs_by_GSE196696 <- background_results_GSE196696[
  background_results_GSE196696$GSE196696_pass == FALSE |
    is.na(background_results_GSE196696$GSE196696_pass),
]

dim(final_dmrs_after_GSE196696)
dim(removed_dmrs_by_GSE196696)

############################################################
## Module 9. Summary table
############################################################

background_summary_GSE196696 <- data.frame(
  dataset = dataset_name,
  Input_DMRs = nrow(background_results_GSE196696),
  
  DMRs_with_CpG_overlap = sum(
    background_results_GSE196696$nCpG_GSE196696 > 0,
    na.rm = TRUE
  ),
  
  DMRs_with_nCpG_gt5 = sum(
    background_results_GSE196696$nCpG_GSE196696 >= min_bg_cpg,
    na.rm = TRUE
  ),
  
  DMRs_with_meanBeta_lt_0.20 = sum(
    background_results_GSE196696$max_group_mean_beta_GSE196696 < beta_threshold,
    na.rm = TRUE
  ),
  
  Final_passed_DMRs = nrow(final_dmrs_after_GSE196696),
  
  Final_retained_percent = round(
    100 * nrow(final_dmrs_after_GSE196696) /
      nrow(background_results_GSE196696),
    2
  )
)

print(background_summary_GSE196696)

############################################################
## Module 10. Save outputs
############################################################

write.csv(
  background_results_GSE196696,
  file = "GSE196696_background_filter_all_DMRs.csv",
  row.names = FALSE
)

write.csv(
  final_dmrs_after_GSE196696,
  file = "final_DMRs_after_GSE196696_whole_blood_filter.csv",
  row.names = FALSE
)

write.csv(
  removed_dmrs_by_GSE196696,
  file = "GSE196696_removed_DMRs_by_whole_blood_filter.csv",
  row.names = FALSE
)

write.csv(
  background_summary_GSE196696,
  file = "GSE196696_background_filter_summary.csv",
  row.names = FALSE
)

step3_GSE196696_background_output <- list(
  candidate_dmrs = candidate_dmrs,
  beta_bg = beta_bg,
  pd_bg = pd_bg,
  dmr_beta_GSE196696 = dmr_beta_GSE196696,
  background_results_GSE196696 = background_results_GSE196696,
  final_dmrs_after_GSE196696 = final_dmrs_after_GSE196696,
  removed_dmrs_by_GSE196696 = removed_dmrs_by_GSE196696,
  background_summary_GSE196696 = background_summary_GSE196696,
  beta_threshold = beta_threshold,
  min_bg_cpg = min_bg_cpg
)

save(
  step3_GSE196696_background_output,
  file = "step3.1-GSE196696-whole-blood-background-filter-output.Rdata"
)

############################################################
## Module 11. Simple visualisation
############################################################

plot_df <- background_results_GSE196696

plot_df$Pass_status <- ifelse(
  plot_df$GSE196696_pass == TRUE,
  "Passed",
  "Not passed"
)

############################################################
## Figure 1. Filtering flow
############################################################

flow_df <- data.frame(
  Step = factor(
    c(
      "Input DMRs",
      "CpG overlap",
      "nCpG > 5",
      "Beta < 0.20",
      "Final passed"
    ),
    levels = c(
      "Input DMRs",
      "CpG overlap",
      "nCpG > 5",
      "Beta < 0.20",
      "Final passed"
    )
  ),
  Count = c(
    nrow(background_results_GSE196696),
    sum(background_results_GSE196696$nCpG_GSE196696 > 0, na.rm = TRUE),
    sum(background_results_GSE196696$nCpG_GSE196696 >= min_bg_cpg, na.rm = TRUE),
    sum(background_results_GSE196696$max_group_mean_beta_GSE196696 < beta_threshold, na.rm = TRUE),
    nrow(final_dmrs_after_GSE196696)
  )
)

p_flow <- ggplot(flow_df, aes(x = Step, y = Count)) +
  geom_col(width = 0.65, fill = "grey70", color = "black") +
  geom_text(aes(label = Count), vjust = -0.4, size = 4) +
  theme_bw(base_size = 13) +
  labs(
    title = "GSE196696 whole blood background filtering",
    x = NULL,
    y = "Number of DMRs"
  ) +
  ylim(0, max(flow_df$Count) * 1.15) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

print(p_flow)

ggsave(
  "GSE196696_filtering_flow.png",
  p_flow,
  width = 7,
  height = 4.8,
  dpi = 300
)

############################################################
## Figure 2. Whole blood beta distribution
############################################################

p_beta <- ggplot(
  plot_df,
  aes(x = max_group_mean_beta_GSE196696)
) +
  geom_histogram(
    bins = 40,
    fill = "grey70",
    color = "black"
  ) +
  geom_vline(
    xintercept = beta_threshold,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "GSE196696 whole blood background methylation",
    x = "Maximum group mean beta in whole blood",
    y = "Number of DMRs"
  )

print(p_beta)

ggsave(
  "GSE196696_whole_blood_beta_distribution.png",
  p_beta,
  width = 7,
  height = 5,
  dpi = 300
)

############################################################
## Figure 3. CpG coverage vs whole blood beta
############################################################

p_scatter <- ggplot(
  plot_df,
  aes(
    x = nCpG_GSE196696,
    y = max_group_mean_beta_GSE196696,
    shape = Pass_status
  )
) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(
    yintercept = beta_threshold,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = min_bg_cpg,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "CpG coverage and whole blood methylation in GSE196696",
    x = "Number of overlapping CpGs",
    y = "Maximum group mean beta in whole blood",
    shape = NULL
  )

print(p_scatter)

ggsave(
  "GSE196696_CpG_vs_whole_blood_beta_scatter.png",
  p_scatter,
  width = 7,
  height = 5,
  dpi = 300
)

```






## Including Plots
You can also embed plots, for example:
```{r pressure, echo=FALSE}
plot(pressure)
```
Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
