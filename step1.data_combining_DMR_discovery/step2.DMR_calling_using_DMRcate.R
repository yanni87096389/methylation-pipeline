rm(list = ls())
options(timeout = 100000) 
options(scipen = 20)


load("step0.1-GSE199057_output.Rdata")
g_GSE199057 <- myLoad$beta
g_GSE199057_group <- myLoad$pd
rm(myLoad)
g1 <- g_GSE199057
pd1 <- g_GSE199057_group

load("step0.2-GSE149282_output.Rdata")
g_GSE149282 <- myLoad$beta
g_GSE149282_group <- myLoad$pd
rm(myLoad)
g2 <-g_GSE149282
pd2 <- g_GSE149282_group

load("step0.3-GSE119526_output.Rdata")
g_GSE119526 <- myLoad$beta
g_GSE119526_group <- myLoad$pd
rm(myLoad)
g3 <-g_GSE119526
pd3 <- g_GSE119526_group

dim(g1)
dim(g2)
dim(g3)

dim(pd1)
dim(pd2)
dim(pd3)

colnames(pd1)
colnames(pd2)
colnames(pd3)

head(pd1);tail(pd1)
head(pd2);tail(pd2)
head(pd3);tail(pd3)

# unify CpG space
common_cpg <- Reduce(intersect, list(
  rownames(g1),
  rownames(g2),
  rownames(g3)
))

g1 <- g1[common_cpg, ]
g2 <- g2[common_cpg, ]
g3 <- g3[common_cpg, ]

all(colnames(g1) == pd1$ID)
all(colnames(g2) == pd2$ID)
all(colnames(g3) == pd3$ID)

# build phenotype
group1 <- factor(pd1$group, levels = c("Normal", "Tumor"))
group2 <- factor(pd2$group, levels = c("Normal", "Tumor"))
group3 <- factor(pd3$group, levels = c("Normal", "Tumor"))

## LIMMA per dataset

run_limma <- function(beta, group){
  
  beta <- as.matrix(beta)
  
  # Fixed effect direction: Tumor - Normal
  group <- factor(group, levels = c("Normal", "Tumor"))
  
  stopifnot(ncol(beta) == length(group))
  stopifnot(!any(is.na(group)))
  
  # Prevent beta values of 0 or 1 from producing infinite M values
  beta <- pmin(pmax(beta, 1e-6), 1 - 1e-6)
  
  # Convert beta values to M values
  mval <- log2(beta / (1 - beta))
  
  # Use a no-intercept design to explicitly define Tumor - Normal
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  contrast.matrix <- makeContrasts(
    Tumor_vs_Normal = Tumor - Normal,
    levels = design
  )
  
  fit <- lmFit(mval, design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit)
  
  res <- topTable(
    fit,
    coef = "Tumor_vs_Normal",
    number = Inf,
    sort.by = "none"
  )
  
  # CpG ID
  res$CpG <- rownames(res)
  
  # Manually calculate standard error
  se_vec <- sqrt(fit$s2.post) * fit$stdev.unscaled[, "Tumor_vs_Normal"]
  names(se_vec) <- rownames(fit$coefficients)
  res$SE <- se_vec[res$CpG]
  
  # Retain delta beta for biological interpretation
  res$deltaBeta <- rowMeans(beta[, group == "Tumor", drop = FALSE]) -
    rowMeans(beta[, group == "Normal", drop = FALSE])
  
  return(res)
}
res1 <- run_limma(g1, group1)
res2 <- run_limma(g2, group2)
res3 <- run_limma(g3, group3)

# check no significant probes
nrow(subset(res1,logFC>0 & adj.P.Val<0.05))
nrow(subset(res1,logFC<0 & adj.P.Val<0.05))

nrow(subset(res2,logFC>0 & adj.P.Val<0.05))
nrow(subset(res2,logFC<0 & adj.P.Val<0.05))

nrow(subset(res3,logFC>0 & adj.P.Val<0.05))
nrow(subset(res3,logFC<0 & adj.P.Val<0.05))

# check similarity between
cor(res1$logFC,res2$logFC,use="complete.obs")
cor(res1$logFC,res3$logFC,use="complete.obs")
cor(res2$logFC,res3$logFC,use="complete.obs")

# CpG meta-analysis
meta_df <- data.frame(
  CpG = res1$CpG,
  
  beta1 = res1$logFC,
  se1   = res1$SE,
  
  beta2 = res2$logFC,
  se2   = res2$SE,
  
  beta3 = res3$logFC,
  se3   = res3$SE
)
## meta-analysis per CpG
## fixed-effect meta-analysis per CpG

b1 <- meta_df$beta1
b2 <- meta_df$beta2
b3 <- meta_df$beta3

s1 <- meta_df$se1
s2 <- meta_df$se2
s3 <- meta_df$se3

## inverse-variance weights
w1 <- 1 / s1^2
w2 <- 1 / s2^2
w3 <- 1 / s3^2

## fixed-effect pooled beta
beta_meta <- (w1 * b1 + w2 * b2 + w3 * b3) / (w1 + w2 + w3)

## pooled standard error
se_meta <- sqrt(1 / (w1 + w2 + w3))

## z statistic
z_meta <- beta_meta / se_meta

## two-sided p value
p_meta <- 2 * pnorm(abs(z_meta), lower.tail = FALSE)

## output
meta_res_fixed <- data.frame(
  CpG = meta_df$CpG,
  beta = beta_meta,
  se = se_meta,
  z = z_meta,
  p = p_meta
)

## FDR correction
meta_res_fixed$FDR <- p.adjust(meta_res_fixed$p, method = "BH")

###############################################################################
## DMR calling using DMRcate from meta CpG signal
## Replace the whole previous custom DMR calling section with this chunk

############################################################
## Module 1. Load required packages
############################################################

library(DMRcate)
library(GenomicRanges)
library(data.table)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(kableExtra)

############################################################
## Module 2. Annotation
############################################################
## We keep the same EPIC hg19 annotation as in the previous custom method.

ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

ann_df <- data.frame(
  CpG = rownames(ann),
  chr = as.character(ann$chr),
  pos = as.integer(ann$pos),
  gene = as.character(ann$UCSC_RefGene_Name),
  relation_to_island = as.character(ann$Relation_to_Island),
  stringsAsFactors = FALSE
)

############################################################
## Module 3. Prepare mean delta beta
############################################################
## deltaBeta direction: Tumor - Normal
## Positive value means hypermethylated in Tumor.

delta_df <- Reduce(function(x, y) merge(x, y, by = "CpG", all = FALSE), list(
  data.frame(CpG = res1$CpG, deltaBeta1 = res1$deltaBeta),
  data.frame(CpG = res2$CpG, deltaBeta2 = res2$deltaBeta),
  data.frame(CpG = res3$CpG, deltaBeta3 = res3$deltaBeta)
))

delta_df$mean_deltaBeta <- rowMeans(
  delta_df[, c("deltaBeta1", "deltaBeta2", "deltaBeta3")],
  na.rm = TRUE
)

############################################################
## Module 4. Merge meta result + annotation + delta beta
############################################################
## meta_res_fixed contains:
## CpG, beta, se, z, p, FDR
##
## beta here is the meta-analysis effect size in M-value/logFC scale.
## mean_deltaBeta is retained for biological interpretation and filtering.

dat <- merge(meta_res_fixed, ann_df, by = "CpG")
dat <- merge(dat, delta_df[, c("CpG", "mean_deltaBeta")], by = "CpG")

dat <- dat[
  is.finite(dat$beta) &
    is.finite(dat$se) &
    is.finite(dat$z) &
    is.finite(dat$p) &
    is.finite(dat$FDR) &
    is.finite(dat$mean_deltaBeta) &
    !is.na(dat$chr) &
    !is.na(dat$pos),
]

## Keep autosomes only, same as your previous method
dat <- dat[dat$chr %in% paste0("chr", 1:22), ]

## Sort by chromosome and position
dat$chr_order <- match(dat$chr, paste0("chr", 1:22))
dat <- dat[order(dat$chr_order, dat$pos), ]

head(dat, 10) |>
  kbl(caption = "CpG data used for DMRcate") |>
  kable_paper("hover", full_width = FALSE)

############################################################
## Module 5. Build CpGannotated object for DMRcate
############################################################
## DMRcate normally uses cpg.annotate() from M-values + design.
## Here, because you already performed CpG-level meta-analysis,
## we manually construct a CpGannotated object using meta statistics.
##
## Required fields:
## stat    = meta z statistic
## rawpval = meta p value
## diff    = mean delta beta, Tumor - Normal
## ind.fdr = meta FDR
## is.sig  = whether the CpG is significant at CpG level

gr <- GRanges(
  seqnames = dat$chr,
  ranges = IRanges(start = dat$pos, end = dat$pos),
  stat = dat$z,
  rawpval = dat$p,
  diff = dat$mean_deltaBeta,
  ind.fdr = dat$FDR,
  is.sig = dat$FDR < 0.05
)

names(gr) <- dat$CpG

myannotation <- new("CpGannotated", ranges = gr)

myannotation

############################################################
## Module 6. Run DMRcate
############################################################
## lambda = 1000:
##   CpGs separated by >= 1000 bp tend to be split into different regions.
##
## C = 2:
##   Recommended commonly used setting for array data.
##
## min.cpgs = 6:
##   Your previous criterion was no.cpgs > 5, so here this becomes at least 6 CpGs.
##
## pcutoff = "fdr":
##   Use the FDR indexing stored in the CpGannotated object.

dmrcoutput <- dmrcate(
  myannotation,
  lambda = 1000,
  C = 2,
  pcutoff = "fdr",
  min.cpgs = 6
)

dmrcoutput

############################################################
## Module 7. Extract DMR ranges
############################################################

dmr_ranges <- extractRanges(dmrcoutput, genome = "hg19")

dmr_ranges

############################################################
## Module 8. Convert DMRcate result to data frame
############################################################

dmr_dmrcate <- as.data.frame(dmr_ranges)

## Make coordinates easier to read
dmr_dmrcate$coord <- paste0(
  dmr_dmrcate$seqnames, ":",
  dmr_dmrcate$start, "-",
  dmr_dmrcate$end
)

## Reorder columns
first_cols <- c("coord", "seqnames", "start", "end", "width", "no.cpgs")
other_cols <- setdiff(colnames(dmr_dmrcate), first_cols)
dmr_dmrcate <- dmr_dmrcate[, c(first_cols, other_cols)]

head(dmr_dmrcate, 10) |>
  kbl(caption = "DMRcate DMR result") |>
  kable_paper("hover", full_width = FALSE)

dim(dmr_dmrcate)

############################################################
## Module 9. Final filter for hypermethylated candidate DMRs
############################################################
## Because we set diff = mean_deltaBeta above,
## DMRcate's meandiff should represent mean Tumor - Normal delta beta
## across CpGs in each DMR.
##
## Keep the same biological filter as before:
## 1. at least 6 CpGs
## 2. mean delta beta > 0.2
## 3. hypermethylated in Tumor

dmr_hyper_dmrcate <- dmr_dmrcate[
  dmr_dmrcate$no.cpgs >= 6 &
    dmr_dmrcate$meandiff > 0.2,
]

## Sort by strongest regional evidence, then largest methylation difference
## DMRcate output usually contains min_smoothed_fdr.
if ("min_smoothed_fdr" %in% colnames(dmr_hyper_dmrcate)) {
  dmr_hyper_dmrcate <- dmr_hyper_dmrcate[
    order(dmr_hyper_dmrcate$min_smoothed_fdr, -dmr_hyper_dmrcate$meandiff),
  ]
} else {
  dmr_hyper_dmrcate <- dmr_hyper_dmrcate[
    order(-dmr_hyper_dmrcate$meandiff),
  ]
}

head(dmr_hyper_dmrcate, 10) |>
  kbl(caption = "Hypermethylated DMRs called by DMRcate") |>
  kable_paper("hover", full_width = FALSE)

dim(dmr_hyper_dmrcate)

############################################################
## Module 10. Save output
############################################################

save(
  dmr_hyper_dmrcate,
  dmr_dmrcate,
  dmr_ranges,
  dmrcoutput,
  myannotation,
  file = "dmr_hyper_dmrcate_meta.Rdata"
)

write.csv(
  dmr_hyper_dmrcate,
  file = "dmr_hyper_dmrcate_meta.csv",
  row.names = FALSE
)
