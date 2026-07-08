library(GenomicRanges)
library(data.table)

############################################################
## 1. Load two DMR result tables
############################################################

load("dmr_hyper_custom.Rdata")
load("dmr_hyper_dmrcate_meta.Rdata")

############################################################
## 2. Convert custom DMRs to GRanges
############################################################
## custom table columns:
## chr, start, end

gr_custom <- GRanges(
  seqnames = as.character(dmr_hyper_custom$chr),
  ranges = IRanges(
    start = as.integer(dmr_hyper_custom$start),
    end = as.integer(dmr_hyper_custom$end)
  )
)

############################################################
## 3. Convert DMRcate DMRs to GRanges
############################################################
## DMRcate table columns:
## seqnames, start, end

gr_dmrcate <- GRanges(
  seqnames = as.character(dmr_hyper_dmrcate$seqnames),
  ranges = IRanges(
    start = as.integer(dmr_hyper_dmrcate$start),
    end = as.integer(dmr_hyper_dmrcate$end)
  )
)

############################################################
## 4. Find overlapping regions
############################################################

hits <- findOverlaps(
  gr_custom,
  gr_dmrcate,
  ignore.strand = TRUE
)

hits
length(hits)
############################################################
## 5. Count overlap numbers
############################################################

n_custom <- length(gr_custom)
n_dmrcate <- length(gr_dmrcate)

custom_overlap_ids <- unique(queryHits(hits))
dmrcate_overlap_ids <- unique(subjectHits(hits))

n_custom_overlap <- length(custom_overlap_ids)
n_dmrcate_overlap <- length(dmrcate_overlap_ids)

prop_custom_overlap <- n_custom_overlap / n_custom
prop_dmrcate_overlap <- n_dmrcate_overlap / n_dmrcate

overlap_summary <- data.frame(
  comparison = c(
    "Total custom hyper-DMRs",
    "Total DMRcate hyper-DMRs",
    "Custom DMRs overlapping DMRcate",
    "DMRcate DMRs overlapping custom",
    "Proportion of custom supported by DMRcate",
    "Proportion of DMRcate supported by custom"
  ),
  value = c(
    n_custom,
    n_dmrcate,
    n_custom_overlap,
    n_dmrcate_overlap,
    prop_custom_overlap,
    prop_dmrcate_overlap
  )
)

overlap_summary############################################################
## 6. Build overlap pair table
############################################################

custom_id <- queryHits(hits)
dmrcate_id <- subjectHits(hits)

overlap_pairs <- data.frame(
  custom_id = custom_id,
  custom_chr = as.character(seqnames(gr_custom))[custom_id],
  custom_start = start(gr_custom)[custom_id],
  custom_end = end(gr_custom)[custom_id],
  custom_width = width(gr_custom)[custom_id],
  
  dmrcate_id = dmrcate_id,
  dmrcate_chr = as.character(seqnames(gr_dmrcate))[dmrcate_id],
  dmrcate_start = start(gr_dmrcate)[dmrcate_id],
  dmrcate_end = end(gr_dmrcate)[dmrcate_id],
  dmrcate_width = width(gr_dmrcate)[dmrcate_id]
)

head(overlap_pairs)
dim(overlap_pairs)
############################################################
## 7. Calculate overlap width and reciprocal overlap percentage
############################################################

intersect_ranges <- pintersect(
  gr_custom[custom_id],
  gr_dmrcate[dmrcate_id]
)

overlap_pairs$overlap_bp <- width(intersect_ranges)

overlap_pairs$pct_custom_overlap <- 
  overlap_pairs$overlap_bp / overlap_pairs$custom_width

overlap_pairs$pct_dmrcate_overlap <- 
  overlap_pairs$overlap_bp / overlap_pairs$dmrcate_width

head(overlap_pairs)
