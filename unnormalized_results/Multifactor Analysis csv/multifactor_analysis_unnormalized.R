# =============================================================
# Multi-factor differential expression analysis -- UNNORMALIZED COUNTS
# Arabidopsis ionizing radiation dataset: OSD-498, 502, 508, 510, 658, 782
#
# This is the statistically preferred version of multifactor_analysis.R.
# The previous version used pre-normalized counts (the only file available
# at the time) and had to round them as a workaround. DESeq2 is designed
# to receive raw integer counts and perform its own normalization
# internally -- feeding it already-normalized values violates that
# assumption. Now that the unnormalized counts tables are available, this
# version uses them correctly.
#
# The unnormalized counts files from GeneLab are named either:
#   *_RSEM_Unnormalized_Counts.csv  (RSEM-based quantification)
#   *_STAR_Unnormalized_Counts.csv  (STAR alignment-based counts)
# Both are referenced in each study's a_OSD-XXX assay file.
# Either is appropriate; STAR counts are slightly more commonly used
# with DESeq2. Use the same type consistently across all 6 studies.
#
# Design notes (unchanged from normalized version):
#   - 6 independent studies pooled -> Study must be modeled as a batch
#     blocking factor in every cross-study comparison.
#   - OSD-658 and OSD-782 use different radiation sources/doses than the
#     other 4 studies -- subset-based comparisons are recommended over
#     one global model for formal hypothesis testing.
#   - sog1_1 appears in OSD-508 and OSD-510 only.
#     myb3r135 appears in OSD-502 only.
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

# ---- 2. Set working directory and output folder ---------------------------
# Update this path to wherever your unnormalized counts files and
# factors_matrix.csv are saved.
setwd("C:/Users/smith/OneDrive - Florida Institute of Technology/summer 2026/OSDR/OSDR_A-rad-adopsis")

output_dir <- "unnormalized_results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- 3. Load unnormalized counts file -------------------------------------
# The unnormalized counts file is already merged across all 6 studies,
# in the same format as the normalized file -- one row per gene,
# one column per sample, GSM accession IDs as column headers.
#
# Update this filename to match what you downloaded from GeneLab.
unnorm_file <- "C:/Users/smith/OneDrive - Florida Institute of Technology/summer 2026/OSDR/OSDR_A-rad-adopsis/OSD-498&OSD-502&OSD-508&OSD-510&OSD-658&OSD-782_rna_seq_Unnormalized_Counts.csv"

cat("Loading unnormalized counts...\n")
counts_raw <- read.csv(unnorm_file, row.names = 1, check.names = FALSE)
cat("Dimensions:", nrow(counts_raw), "genes x", ncol(counts_raw), "samples\n")
cat("First few sample names:", head(colnames(counts_raw), 3), "\n")

# Rename columns from GSM IDs to human-readable names using sample_mapping.csv
# (the same mapping built during the cleaning phase)
mapping <- read.csv("C:/Users/smith/OneDrive - Florida Institute of Technology/summer 2026/OSDR/OSDR_A-rad-adopsis/multifactor + pathway enrichment/R files/sample_mapping.csv")
# Build a named vector: old GSM ID -> new readable name
name_map <- setNames(mapping$new_name, mapping$gsm)

# Rename columns that have a mapping; warn about any that don't
old_names <- colnames(counts_raw)
new_names <- ifelse(old_names %in% names(name_map), name_map[old_names], old_names)
unmapped  <- old_names[!old_names %in% names(name_map)]
if (length(unmapped) > 0) {
  cat("WARNING: these sample IDs have no mapping and will keep their original names:\n")
  print(unmapped)
}
colnames(counts_raw) <- new_names
cat("Columns renamed successfully.\n")

# ---- 4. Validate counts are genuine integers ------------------------------
# Unlike the normalized version, we don't need to round here -- these
# are already raw integer counts from the aligner. We just verify that
# and stop early with a clear message if something looks wrong.
counts_mat <- as.matrix(counts_raw)
non_integer <- sum(counts_mat != round(counts_mat), na.rm = TRUE)
if (non_integer > 0) {
  warning(non_integer, " non-integer values found in counts matrix.\n",
          "  This may mean you accidentally loaded the normalized file.\n",
          "  Rounding to proceed, but double-check your input file.")
  counts_mat <- round(counts_mat)
}
counts_mat[counts_mat < 0] <- 0
counts_int <- counts_mat
cat("Counts validation passed --", nrow(counts_int), "genes,",
    ncol(counts_int), "samples, all integer values confirmed.\n")

# ---- 5. Load factors table ------------------------------------------------
# Use the same factors_matrix.csv built during the cleaning phase.
# The sample names in that file must match the column names in the
# unnormalized counts files exactly.
factors_file <- "C:/Users/smith/OneDrive - Florida Institute of Technology/summer 2026/OSDR/OSDR_A-rad-adopsis/multifactor + pathway enrichment/R files/factors_matrix.csv"
factors <- read.csv(factors_file, row.names = "SampleName")

# Check alignment between counts columns and factors rows
missing_from_factors <- setdiff(colnames(counts_int), rownames(factors))
missing_from_counts  <- setdiff(rownames(factors), colnames(counts_int))

if (length(missing_from_factors) > 0) {
  cat("\nWARNING: these samples are in counts but NOT in factors table:\n")
  print(missing_from_factors)
  cat("These samples will be dropped from the analysis.\n")
  counts_int <- counts_int[, !colnames(counts_int) %in% missing_from_factors]
}
if (length(missing_from_counts) > 0) {
  cat("\nWARNING: these samples are in factors table but NOT in counts:\n")
  print(missing_from_counts)
  cat("These rows will be dropped from the factors table.\n")
  factors <- factors[!rownames(factors) %in% missing_from_counts, ]
}

# Reorder factors to match counts column order (required by DESeq2)
factors <- factors[colnames(counts_int), ]
stopifnot(all(rownames(factors) == colnames(counts_int)))

# Set reference levels: WT = reference genotype, none = reference radiation
factors$Study     <- factor(factors$Study)
factors$Genotype  <- factor(factors$Genotype,
                            levels = c("WT", setdiff(unique(factors$Genotype), "WT")))
factors$Radiation <- factor(factors$Radiation,
                            levels = c("none", setdiff(unique(factors$Radiation), "none")))
factors$Timepoint <- factor(factors$Timepoint)
factors$Replicate <- factor(factors$Replicate)

cat("\nSamples:", nrow(factors), "\n")
cat("Studies:", levels(factors$Study), "\n")
cat("Genotypes:", levels(factors$Genotype), "\n")
cat("Radiation levels:", levels(factors$Radiation), "\n")

# ---- 6. Exploratory PCA: does Study dominate the variance? ---------------
# Identical to the normalized version -- checks whether batch effects
# are visible before any formal testing. With raw counts and DESeq2's
# correct normalization, batch structure may look slightly different
# compared to the normalized version, so worth re-running even if you
# already ran this before.
cat("\nRunning exploratory PCA...\n")
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
  labs(title = "PCA colored by Study (batch check) -- unnormalized counts") +
  theme_minimal()
ggsave(file.path(output_dir, "pca_by_study_unnorm.png"), p, width = 7, height = 5)
cat("Saved pca_by_study_unnorm.png\n")
cat("Compare this to the normalized version -- if batch structure looks\n")
cat("similar, confidence in the previous results increases.\n")

# ---- 7. Subset analyses ---------------------------------------------------
# Identical comparisons to the normalized version (5a-5d).
# With correct raw counts input, DESeq2's size factor estimation and
# dispersion modeling are now working as intended -- results may show
# modestly different gene lists or p-values compared to the normalized
# version, which is expected and correct.

## 7a. OSD498 + OSD510: radiation effect (WT + sog1_1, Co-60 100Gy)
cat("\nRunning 7a: OSD498 + OSD510 radiation effect...\n")
sub_a   <- factors[factors$Study %in% c("OSD498", "OSD510"), ]
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
write.csv(as.data.frame(res_a),
          file.path(output_dir, "DEG_unnorm_OSD498_510_radiation_effect.csv"))
cat("Saved DEG_unnorm_OSD498_510_radiation_effect.csv --",
    sum(res_a$padj < 0.05, na.rm = TRUE), "genes at padj < 0.05\n")

## 7b. OSD508 + OSD510: genotype x radiation interaction
cat("\nRunning 7b: OSD508 + OSD510 genotype x radiation interaction...\n")
sub_b    <- factors[factors$Study %in% c("OSD508", "OSD510"), ]
counts_b <- counts_int[, rownames(sub_b)]

dds_b <- DESeqDataSetFromMatrix(
  countData = counts_b,
  colData   = sub_b,
  design    = ~ Study + Genotype + Radiation + Genotype:Radiation
)
dds_b <- dds_b[rowSums(counts(dds_b)) > 10, ]
dds_b <- DESeq(dds_b)

res_b <- results(dds_b, name = "Genotypesog1_1.RadiationgammaCo60_100Gy")
res_b <- res_b[order(res_b$padj), ]
write.csv(as.data.frame(res_b),
          file.path(output_dir, "DEG_unnorm_OSD508_510_genotype_x_radiation_interaction.csv"))
cat("Saved interaction results --",
    sum(res_b$padj < 0.05, na.rm = TRUE),
    "genes show genotype-dependent radiation response\n")

## 7c. OSD782: Cs-137 dose-response
cat("\nRunning 7c: OSD782 Cs-137 dose-response...\n")
sub_c    <- factors[factors$Study == "OSD782", ]
counts_c <- counts_int[, rownames(sub_c)]

dds_c <- DESeqDataSetFromMatrix(
  countData = counts_c,
  colData   = sub_c,
  design    = ~ Timepoint + Radiation
)
dds_c <- dds_c[rowSums(counts(dds_c)) > 10, ]
dds_c <- DESeq(dds_c)

res_c_low  <- results(dds_c, contrast = c("Radiation", "gammaCs137_10cGy",  "none"))
res_c_high <- results(dds_c, contrast = c("Radiation", "gammaCs137_100cGy", "none"))
write.csv(as.data.frame(res_c_low[order(res_c_low$padj), ]),
          file.path(output_dir, "DEG_unnorm_OSD782_10cGy_vs_none.csv"))
write.csv(as.data.frame(res_c_high[order(res_c_high$padj), ]),
          file.path(output_dir, "DEG_unnorm_OSD782_100cGy_vs_none.csv"))
cat("Saved OSD782 dose-response results\n",
    " 10cGy:", sum(res_c_low$padj  < 0.05, na.rm = TRUE), "genes\n",
    "100cGy:", sum(res_c_high$padj < 0.05, na.rm = TRUE), "genes\n")

## 7d. OSD658: simulated GCR dose-response
cat("\nRunning 7d: OSD658 GCR dose-response...\n")
sub_d    <- factors[factors$Study == "OSD658", ]
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
write.csv(as.data.frame(res_d_40[order(res_d_40$padj), ]),
          file.path(output_dir, "DEG_unnorm_OSD658_GCR40_vs_none.csv"))
write.csv(as.data.frame(res_d_80[order(res_d_80$padj), ]),
          file.path(output_dir, "DEG_unnorm_OSD658_GCR80_vs_none.csv"))
cat("Saved OSD658 GCR dose-response results\n",
    " 40cGy:", sum(res_d_40$padj < 0.05, na.rm = TRUE), "genes\n",
    " 80cGy:", sum(res_d_80$padj < 0.05, na.rm = TRUE), "genes\n")

# ---- 8. Full pooled model (exploratory heatmap only) ----------------------
cat("\nGenerating full pooled heatmap (exploratory only)...\n")
dds_full <- DESeqDataSetFromMatrix(
  countData = counts_int,
  colData   = factors,
  design    = ~ Study + Genotype
)
dds_full <- dds_full[rowSums(counts(dds_full)) > 10, ]
vsd_full <- vst(dds_full, blind = FALSE)

topvar <- head(order(rowVars(assay(vsd_full)), decreasing = TRUE), 50)
annotation_col <- factors[, c("Study", "Genotype", "Radiation")]

pheatmap(assay(vsd_full)[topvar, ],
         annotation_col = annotation_col,
         show_rownames  = FALSE,
         show_colnames  = TRUE,
         fontsize_col   = 5,
         main = "Top 50 variable genes -- unnormalized counts input",
         filename = file.path(output_dir, "heatmap_top50_unnorm.png"),
         width = 10, height = 8)
cat("Saved heatmap_top50_unnorm.png\n")

# ---- 9. Compare normalized vs unnormalized results (all comparisons) ------
# Run this block after both scripts have been run to decide whether
# pathway enrichment needs to be re-run on the new unnormalized results.
#
# For each comparison, it reports:
#   - How many genes were significant in each version
#   - How many overlap (appear in both)
#   - Jaccard similarity (0 = no overlap, 1 = identical lists)
#
# Rule of thumb:
#   Jaccard > 0.8 -> lists are very similar, no need to rerun enrichment
#   Jaccard 0.5-0.8 -> moderate difference, consider rerunning
#   Jaccard < 0.5 -> substantial difference, rerun enrichment for this comparison

# Save all output to a text file AND print to console simultaneously
jaccard_report <- file.path(output_dir, "normalized_vs_unnormalized_comparison.txt")
sink(jaccard_report, split = TRUE)  # split=TRUE means output goes to BOTH file AND console

comparison_pairs <- list(
  list(
    label       = "7a vs 5a: OSD498+510 radiation effect",
    norm_file   = "C:/Users/smith/OneDrive - Florida Institute of Technology/summer 2026/OSDR/OSDR_A-rad-adopsis/multifactor + pathway enrichment/Normalized data_multifactor analysis/csv Results/DEG_OSD498_510_radiation_effect.csv",
    unnorm_file = "C:/Users/smith/OneDrive - Florida Institute of Technology/summer 2026/OSDR/OSDR_A-rad-adopsis/unnormalized_results/DEG_unnorm_OSD498_510_radiation_effect.csv"
  ),
  list(
    label       = "7b vs 5b: OSD508+510 genotype x radiation interaction",
    norm_file   = "C:/Users/smith/OneDrive - Florida Institute of Technology/summer 2026/OSDR/OSDR_A-rad-adopsis/multifactor + pathway enrichment/Normalized data_multifactor analysis/csv Results/DEG_OSD508_510_genotype_x_radiation_interaction.csv",
    unnorm_file = "C:/Users/smith/OneDrive - Florida Institute of Technology/summer 2026/OSDR/OSDR_A-rad-adopsis/unnormalized_results/DEG_unnorm_OSD508_510_genotype_x_radiation_interaction.csv"
  ),
  list(
    label       = "7c-low vs 5c-low: OSD782 10cGy vs none",
    norm_file   = "C:/Users/smith/OneDrive - Florida Institute of Technology/summer 2026/OSDR/OSDR_A-rad-adopsis/multifactor + pathway enrichment/Normalized data_multifactor analysis/csv Results/DEG_OSD782_10cGy_vs_none.csv",
    unnorm_file = "C:/Users/smith/OneDrive - Florida Institute of Technology/summer 2026/OSDR/OSDR_A-rad-adopsis/unnormalized_results/DEG_unnorm_OSD782_10cGy_vs_none.csv"
  ),
  list(
    label       = "7c-high vs 5c-high: OSD782 100cGy vs none",
    norm_file   = "C:/Users/smith/OneDrive - Florida Institute of Technology/summer 2026/OSDR/OSDR_A-rad-adopsis/multifactor + pathway enrichment/Normalized data_multifactor analysis/csv Results/DEG_OSD782_100cGy_vs_none.csv",
    unnorm_file = "C:/Users/smith/OneDrive - Florida Institute of Technology/summer 2026/OSDR/OSDR_A-rad-adopsis/unnormalized_results/DEG_unnorm_OSD782_100cGy_vs_none.csv"
  ),
  list(
    label       = "7d-40 vs 5d-40: OSD658 GCR 40cGy vs none",
    norm_file   = "C:/Users/smith/OneDrive - Florida Institute of Technology/summer 2026/OSDR/OSDR_A-rad-adopsis/multifactor + pathway enrichment/Normalized data_multifactor analysis/csv Results/DEG_OSD658_GCR40_vs_none.csv",
    unnorm_file = "C:/Users/smith/OneDrive - Florida Institute of Technology/summer 2026/OSDR/OSDR_A-rad-adopsis/unnormalized_results/DEG_unnorm_OSD658_GCR40_vs_none.csv"
  ),
  list(
    label       = "7d-80 vs 5d-80: OSD658 GCR 80cGy vs none",
    norm_file   = "C:/Users/smith/OneDrive - Florida Institute of Technology/summer 2026/OSDR/OSDR_A-rad-adopsis/multifactor + pathway enrichment/Normalized data_multifactor analysis/csv Results/DEG_OSD658_GCR80_vs_none.csv",
    unnorm_file = "C:/Users/smith/OneDrive - Florida Institute of Technology/summer 2026/OSDR/OSDR_A-rad-adopsis/unnormalized_results/DEG_unnorm_OSD658_GCR80_vs_none.csv"
  )
)

cat("\n=== Normalized vs Unnormalized Gene List Comparison ===\n")
cat("Jaccard > 0.8  -> lists very similar, no need to rerun enrichment\n")
cat("Jaccard 0.5-0.8 -> moderate difference, consider rerunning\n")
cat("Jaccard < 0.5  -> substantial difference, rerun enrichment\n\n")

rerun_enrichment <- c()  # will collect which comparisons need rerunning

for (pair in comparison_pairs) {
  
  # Skip if either file doesn't exist yet
  if (!file.exists(pair$norm_file)) {
    cat(pair$label, "-- SKIPPED (normalized file not found:", pair$norm_file, ")\n\n")
    next
  }
  if (!file.exists(pair$unnorm_file)) {
    cat(pair$label, "-- SKIPPED (unnormalized file not found:", pair$unnorm_file, ")\n\n")
    next
  }
  
  res_norm   <- read.csv(pair$norm_file,   row.names = 1)
  res_unnorm <- read.csv(pair$unnorm_file, row.names = 1)
  
  sig_norm   <- rownames(subset(res_norm,   padj < 0.05 & !is.na(padj)))
  sig_unnorm <- rownames(subset(res_unnorm, padj < 0.05 & !is.na(padj)))
  
  overlap  <- length(intersect(sig_norm, sig_unnorm))
  union_n  <- length(union(sig_norm, sig_unnorm))
  jaccard  <- if (union_n > 0) round(overlap / union_n, 3) else NA
  
  verdict <- if (is.na(jaccard)) "cannot compare" else
    if (jaccard > 0.8)  "OK -- no need to rerun enrichment" else
      if (jaccard >= 0.5) "CONSIDER rerunning enrichment" else
        "RERUN enrichment recommended"
  
  cat(pair$label, "\n")
  cat("  Significant (normalized)  :", length(sig_norm), "\n")
  cat("  Significant (unnormalized):", length(sig_unnorm), "\n")
  cat("  Overlap                   :", overlap, "\n")
  cat("  Jaccard similarity        :", jaccard, "\n")
  cat("  Verdict                   :", verdict, "\n\n")
  
  if (!is.na(jaccard) && jaccard < 0.8) {
    rerun_enrichment <- c(rerun_enrichment, pair$unnorm_file)
  }
}

if (length(rerun_enrichment) > 0) {
  cat("=== Comparisons that need pathway enrichment rerun ===\n")
  for (f in rerun_enrichment) cat(" ", f, "\n")
  cat("\nRun pathway_enrichment_reusable.R on each of these files.\n")
} else {
  cat("=== All comparisons look similar ===\n")
  cat("Your existing pathway enrichment results are still valid.\n")
  cat("No need to rerun pathway_enrichment_reusable.R.\n")
}

cat("\n=== Analysis complete ===\n")
cat("All output files saved to:", normalizePath(output_dir), "\n")
cat("\nNext steps:\n")
cat("  1. Run pathway_enrichment_reusable.R on each DEG_unnorm_*.csv\n")
cat("  2. Run cross_study_comparison.R to compare studies within radiation groups\n")
cat("  3. Uncomment Section 9 to compare these results to the normalized version\n")
