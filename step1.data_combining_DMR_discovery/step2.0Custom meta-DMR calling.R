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
## Custom meta-DMR calling from meta CpG signal
## Criteria:
## 1. CpG-level meta FDR < 0.05
## 2. Adjacent significant CpGs within 1000 bp are merged
## 3. DMR CpGs > 5
## 4. DMR mean Δβ > 0.2


## 1. Annotation

ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

ann_df <- data.frame(
  CpG = rownames(ann),
  chr = as.character(ann$chr),
  pos = as.integer(ann$pos),
  gene = as.character(ann$UCSC_RefGene_Name),
  relation_to_island = as.character(ann$Relation_to_Island),
  stringsAsFactors = FALSE
)

## 2. Prepare mean delta beta
## deltaBeta direction: Tumor - Normal

delta_df <- Reduce(function(x, y) merge(x, y, by = "CpG", all = FALSE), list(
  data.frame(CpG = res1$CpG, deltaBeta1 = res1$deltaBeta),
  data.frame(CpG = res2$CpG, deltaBeta2 = res2$deltaBeta),
  data.frame(CpG = res3$CpG, deltaBeta3 = res3$deltaBeta)
))

delta_df$mean_deltaBeta <- rowMeans(
  delta_df[, c("deltaBeta1", "deltaBeta2", "deltaBeta3")],
  na.rm = TRUE
)

## 3. Merge meta result + annotation + delta beta

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

dat <- dat[dat$chr %in% paste0("chr", 1:22), ]

dat <- as.data.table(dat)

head(dat,10) |>
  kbl(caption="CpG data") |>
  kable_paper("hover", full_width = F)

## 4. Define significant CpGs for DMR seed
## Hyper-DMR: Tumor - Normal > 0
sig <- dat[
  FDR < 0.05 &
    mean_deltaBeta > 0
]

sig[, chr_order := match(chr, paste0("chr", 1:22))]
setorder(sig, chr_order, pos)

## 5. Merge nearby significant CpGs into regions
## gap <= 1000 bp belongs to the same region
max_gap <- 1000

sig[, gap := pos - data.table::shift(pos, n = 1L, type = "lag"), by = chr]
sig[, new_region := is.na(gap) | gap > max_gap, by = chr]
sig[, region_id := cumsum(new_region)]

## 6. Summarize candidate DMRs
dmr_custom <- sig[
  ,
  .(
    chr = chr[1],
    start = min(pos),
    end = max(pos),
    width = max(pos) - min(pos) + 1,
    no.cpgs = .N,
    mean_deltaBeta = mean(mean_deltaBeta, na.rm = TRUE),
    max_deltaBeta = max(mean_deltaBeta, na.rm = TRUE),
    min_deltaBeta = min(mean_deltaBeta, na.rm = TRUE),
    min_FDR = min(FDR, na.rm = TRUE),
    mean_FDR = mean(FDR, na.rm = TRUE),
    min_p = min(p, na.rm = TRUE),
    mean_z = mean(z, na.rm = TRUE),
    CpGs = paste(CpG, collapse = ";"),
    genes = paste(
      base::unique(
        base::unlist(
          base::strsplit(paste(gene, collapse = ";"), ";")
        )
      ),
      collapse = ";"
    ),
    relation_to_island = paste(base::unique(relation_to_island), collapse = ";")
  ),
  by = region_id
]

head(dmr_custom,10) |>
  kbl(caption="DMR data") |>
  kable_paper("hover", full_width = F)

dim(dmr_custom)

## 7. Final filter: CpGs > 5 and mean Δβ > 0.2
dmr_hyper_custom <- dmr_custom[
  no.cpgs > 5 &
    mean_deltaBeta > 0.2
]

setorder(dmr_hyper_custom, min_FDR, -mean_deltaBeta)

head(dmr_hyper_custom,10) |>
  kbl(caption="DMR hyper") |>
  kable_paper("hover", full_width = F)

dim(dmr_hyper_custom)

#dmr_hyper_custom
save(dmr_hyper_custom,file = "dmr_hyper_custom.Rdata")
