# =============================================================
# Multi-factor differential expression analysis
# Arabidopsis ionizing radiation dataset: OSD-498, 502, 508, 510, 658, 782
#
# Design notes:
#   - 6 independent studies pooled together -> "Study" must be modeled
#     as a blocking factor (batch), or study-to-study technical variation
#     will dominate true radiation/genotype effects.
#   - OSD-658 and OSD-782 use different radiation sources/doses than the
#     other 4 studies (simulated GCR, and Cs-137 vs Co-60), so a single
#     "irradiated vs not" comparison across ALL samples is not always
#     biologically meaningful -- consider subsetting by study/radiation
#     type for specific comparisons, using the full model only for
#     overall variance partitioning / exploratory analysis (PCA, clustering).
#   - sog1-1 appears in OSD-508 and OSD-510 only; myb3r135 in OSD-502 only.
#     WT is the only genotype common across ALL 6 studies.
# =============================================================

# ---- 1. Install/load packages -------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c("DESeq2", "limma", "edgeR", "pheatmap", "ggplot2", "RColorBrewer")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p %in% c("DESeq2", "limma", "edgeR")) {
      BiocManager::install(p, update = FALSE, ask = FALSE)
    } else {
      install.packages(p)
    }
  }
}

library(DESeq2)
library(limma)
library(edgeR)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)

# ---- 2. Load data ---------------------------------------------------------
# Update these paths to wherever you saved the two files
counts_file  <- "renamed_counts.csv"
factors_file <- "factors_matrix.csv"

counts_raw <- read.csv(counts_file, row.names = 1, check.names = FALSE)
factors    <- read.csv(factors_file, row.names = "SampleName")

# Sanity check: sample names must match exactly between the two files
stopifnot(all(colnames(counts_raw) == rownames(factors)))

# Make factor columns actual factors, with WT / none as the reference level
factors$Study     <- factor(factors$Study)
factors$Genotype  <- factor(factors$Genotype, levels = c("WT", setdiff(unique(factors$Genotype), "WT")))
factors$Radiation <- factor(factors$Radiation, levels = c("none", setdiff(unique(factors$Radiation), "none")))
factors$Timepoint <- factor(factors$Timepoint)
factors$Replicate <- factor(factors$Replicate)

cat("Samples:", nrow(factors), "\n")
cat("Studies:", levels(factors$Study), "\n")
cat("Genotypes:", levels(factors$Genotype), "\n")
cat("Radiation levels:", levels(factors$Radiation), "\n")

# ---- 3. Round normalized counts for DESeq2 --------------------------------
# DESeq2 expects integer-like counts; these are already-normalized values
# from GeneLab, so round them and treat as offset-corrected counts.
# (If you have access to the RAW/unnormalized counts tables instead,
#  prefer those for DESeq2 -- they're statistically more appropriate.
#  This script uses what's available: the normalized counts table.)
counts_int <- round(as.matrix(counts_raw))
counts_int[counts_int < 0] <- 0

# ---- 4. Exploratory analysis: does Study dominate the variance? ----------
dds_explore <- DESeqDataSetFromMatrix(
  countData = counts_int,
  colData   = factors,
  design    = ~ Study
)
dds_explore <- dds_explore[rowSums(counts(dds_explore)) > 10, ]
vsd <- vst(dds_explore, blind = TRUE)

pca_data <- plotPCA(vsd, intgroup = "Study", returnData = TRUE)
p <- ggplot(pca_data, aes(PC1, PC2, color = Study)) +
  geom_point(size = 3) +
  labs(title = "PCA colored by Study (batch check)") +
  theme_minimal()
ggsave("pca_by_study.png", p, width = 7, height = 5)
cat("Saved pca_by_study.png -- check whether samples cluster by Study\n")
cat("(If yes, Study MUST be included as a covariate in every model below.)\n")

# ---- 5. Subset analyses (recommended approach) ----------------------------
# Because radiation source/dose differs by study, the cleanest comparisons
# are WITHIN study or across studies that share the same radiation type.

## 5a. OSD-498 + OSD-510: same design family (WT/sog1-1 x Co-60 100Gy x timepoints)
## Note: OSD-498 has no sog1-1 samples, so this combined model only works
## for the WT-only time course unless you restrict to WT.
sub_a <- factors[factors$Study %in% c("OSD498", "OSD510"), ]
counts_a <- counts_int[, rownames(sub_a)]

dds_a <- DESeqDataSetFromMatrix(
  countData = counts_a,
  colData   = sub_a,
  design    = ~ Study + Genotype + Radiation
)
dds_a <- dds_a[rowSums(counts(dds_a)) > 10, ]
dds_a <- DESeq(dds_a)

res_a <- results(dds_a, contrast = c("Radiation", "gammaCo60_100Gy", "none"))
res_a <- res_a[order(res_a$padj), ]
write.csv(as.data.frame(res_a), "DEG_OSD498_510_radiation_effect.csv")
cat("Saved DEG_OSD498_510_radiation_effect.csv --",
    sum(res_a$padj < 0.05, na.rm = TRUE), "genes at padj < 0.05\n")

## 5b. OSD-508 + OSD-510: full sog1-1 vs WT x timepoint x radiation comparison
sub_b <- factors[factors$Study %in% c("OSD508", "OSD510"), ]
counts_b <- counts_int[, rownames(sub_b)]

dds_b <- DESeqDataSetFromMatrix(
  countData = counts_b,
  colData   = sub_b,
  design    = ~ Study + Genotype + Radiation + Genotype:Radiation
)
dds_b <- dds_b[rowSums(counts(dds_b)) > 10, ]
dds_b <- DESeq(dds_b)

# Does the radiation response differ between WT and sog1-1?
res_b_interaction <- results(dds_b, name = "Genotypesog1_1.RadiationgammaCo60_100Gy")
res_b_interaction <- res_b_interaction[order(res_b_interaction$padj), ]
write.csv(as.data.frame(res_b_interaction), "DEG_OSD508_510_genotype_x_radiation_interaction.csv")
cat("Saved interaction results --",
    sum(res_b_interaction$padj < 0.05, na.rm = TRUE),
    "genes show genotype-dependent radiation response\n")

## 5c. OSD-782: dose-response at each timepoint (Cs-137)
sub_c <- factors[factors$Study == "OSD782", ]
counts_c <- counts_int[, rownames(sub_c)]

dds_c <- DESeqDataSetFromMatrix(
  countData = counts_c,
  colData   = sub_c,
  design    = ~ Timepoint + Radiation
)
dds_c <- dds_c[rowSums(counts(dds_c)) > 10, ]
dds_c <- DESeq(dds_c)

res_c_low  <- results(dds_c, contrast = c("Radiation", "gammaCs137_10cGy", "none"))
res_c_high <- results(dds_c, contrast = c("Radiation", "gammaCs137_100cGy", "none"))
write.csv(as.data.frame(res_c_low[order(res_c_low$padj), ]),  "DEG_OSD782_10cGy_vs_none.csv")
write.csv(as.data.frame(res_c_high[order(res_c_high$padj), ]), "DEG_OSD782_100cGy_vs_none.csv")
cat("Saved OSD782 dose-response results (10cGy and 100cGy vs none)\n")

## 5d. OSD-658: simulated GCR dose-response (40 / 80 cGy)
sub_d <- factors[factors$Study == "OSD658", ]
counts_d <- counts_int[, rownames(sub_d)]

dds_d <- DESeqDataSetFromMatrix(
  countData = counts_d,
  colData   = sub_d,
  design    = ~ Radiation
)
dds_d <- dds_d[rowSums(counts(dds_d)) > 10, ]
dds_d <- DESeq(dds_d)

res_d_40 <- results(dds_d, contrast = c("Radiation", "GCR_40cGy", "none"))
res_d_80 <- results(dds_d, contrast = c("Radiation", "GCR_80cGy", "none"))
write.csv(as.data.frame(res_d_40[order(res_d_40$padj), ]), "DEG_OSD658_GCR40_vs_none.csv")
write.csv(as.data.frame(res_d_80[order(res_d_80$padj), ]), "DEG_OSD658_GCR80_vs_none.csv")
cat("Saved OSD658 GCR dose-response results\n")

# ---- 6. (Optional) full combined model for overall exploration only ------
# Use this ONLY for exploratory clustering/heatmaps across all 6 studies,
# not for formal hypothesis testing across radiation types, since the
# "Radiation" factor isn't biologically comparable across all studies
# (different isotopes, particle types, and doses).
dds_full <- DESeqDataSetFromMatrix(
  countData = counts_int,
  colData   = factors,
  design    = ~ Study + Genotype
)
dds_full <- dds_full[rowSums(counts(dds_full)) > 10, ]
vsd_full <- vst(dds_full, blind = FALSE)

# Heatmap of top 50 most variable genes, annotated by Study + Radiation
topvar <- head(order(rowVars(assay(vsd_full)), decreasing = TRUE), 50)
annotation_col <- factors[, c("Study", "Genotype", "Radiation")]

pheatmap(assay(vsd_full)[topvar, ],
         annotation_col = annotation_col,
         show_rownames  = FALSE,
         show_colnames  = TRUE,
         fontsize_col   = 5,
         main = "Top 50 variable genes across all 6 studies")

cat("\nDone. Review pca_by_study.png first -- if Study is a major source of\n")
cat("variance (likely, given different platforms/years/labs), trust the\n")
cat("per-study/subset comparisons (5a-5d) over any pooled cross-study test.\n")
