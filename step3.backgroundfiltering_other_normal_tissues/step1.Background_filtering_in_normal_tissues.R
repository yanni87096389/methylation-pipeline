############################################################
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