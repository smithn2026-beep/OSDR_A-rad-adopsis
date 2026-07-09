# =============================================================================
# Tissue-Specific Separation of Arabidopsis Radiation DEGs
# =============================================================================
#
# Dataset: DEG_OSD498_510_radiation_effect.csv
#   - Arabidopsis thaliana exposed to ionizing radiation
#   - NASA GeneLab studies OSD-498 and OSD-510
#   - DESeq2 differential expression results (whole seedlings)
#
# Approach:
#   The DEG file has no tissue column — it comes from whole-seedling experiments.
#   To separate DEGs into tissue-specific responses, we annotate each gene by
#   its tissue-specific expression pattern using the AtGenExpress developmental
#   expression atlas (Schmid et al. 2005, Nature Genetics), then compute the
#   Tau tissue-specificity index to assign each DEG to its predominant tissue.
#
# Outputs:
#   1. Hierarchical CSV files (broad organ / sub-tissue)
#   2. Summary table (tissue_deg_summary.csv)
#   3. Visualizations: bar chart, volcano plots, heatmap, UpSet diagram
#
# Author: Biomni (Phylo)
# Date: 2026-07-09
# =============================================================================

# ── 0. SETUP: Install and load packages ───────────────────────────────────────
#    This section ensures all required packages are available.
#    If you already have them installed, this step is fast.

# List of required packages
required_cran <- c("ggplot2", "ggrepel", "UpSetR", "gridExtra")
required_bioc <- c("GEOquery", "ath1121501.db", "AnnotationDbi",
                   "ComplexHeatmap", "circlize", "BiocManager")

# Install missing CRAN packages
for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Install missing Bioconductor packages
for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

# Load all packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(UpSetR)
  library(gridExtra)
  library(grid)
  library(GEOquery)
  library(ath1121501.db)
  library(AnnotationDbi)
  library(ComplexHeatmap)
  library(circlize)
  library(Biobase)
})

# Color palette (colorblind-friendly)
UP_COLOR   <- "#0072B2"  # blue for upregulated
DOWN_COLOR <- "#D55E00"  # vermillion for downregulated

# Tau specificity threshold (standard from the literature)
# Genes with Tau >= 0.6 are considered tissue-specific
TAU_THRESHOLD <- 0.6

# ── 1. LOAD AND VALIDATE THE DEG DATASET ──────────────────────────────────────
#    Read the DESeq2 results CSV and inspect its structure.
#
#    The file has 8 columns:
#      - Column 1 (unnamed): AGI gene ID (e.g., AT5G60250)
#      - baseMean:           Average expression across samples
#      - log2FoldChange:     Log2 fold change (irradiated vs control)
#      - lfcSE:               Standard error of log2FoldChange
#      - stat:                Wald statistic
#      - pvalue:              Raw p-value
#      - padj:                Adjusted p-value (FDR)
#      - Column 8 (unnamed):  DEG flag ("yes" or "no")

# >>> SET YOUR INPUT FILE PATH HERE <<<
deg_file <- "DEG_OSD498_510_radiation_effect.csv"

# Read the CSV
deg <- read.csv(deg_file, header = TRUE, stringsAsFactors = FALSE)

# Rename the unnamed columns
colnames(deg)[1] <- "gene_id"       # First column = gene ID
colnames(deg)[ncol(deg)] <- "deg_flag"  # Last column = yes/no flag

cat("Dataset dimensions:", nrow(deg), "genes x", ncol(deg), "columns\n")
cat("DEGs (flag='yes'):", sum(deg$deg_flag == "yes"), "\n")
cat("Non-DEGs (flag='no'):", sum(deg$deg_flag == "no"), "\n")

# Validate AGI ID format (AT{1-5}Gnnnnn for nuclear genes)
# Some genes are organellar: ATCG* (chloroplast), ATMG* (mitochondria)
# These won't be in the AtGenExpress nuclear atlas and will be flagged "unannotated"
valid_nuclear <- grepl("^AT[1-5]G[0-9]{5}$", deg$gene_id)
cat("Nuclear AGI IDs:", sum(valid_nuclear), "\n")
cat("Organellar IDs (ATCG/ATMG):", sum(!valid_nuclear), "\n")

# Add regulation direction (up/down/non-DEG)
deg$regulation <- ifelse(deg$deg_flag == "yes",
                         ifelse(deg$log2FoldChange > 0, "up", "down"),
                         "non_DEG")

# ── 2. DOWNLOAD THE ATGENEXPRESS DEVELOPMENTAL ATLAS ──────────────────────────
#    The AtGenExpress atlas (Schmid et al. 2005) is the canonical Arabidopsis
#    tissue expression reference. It contains Affymetrix ATH1 microarray data
#    from 79+ diverse tissue samples covering the entire life cycle.
#
#    We download 6 GEO series:
#      GSE5629: Seedlings and whole plants
#      GSE5630: Leaves
#      GSE5631: Roots
#      GSE5632: Flowers and pollen
#      GSE5633: Shoots and stems
#      GSE5634: Siliques and seeds
#
#    >>> This step requires internet access and may take several minutes. <<<

geodir <- "geo_cache"  # Local cache directory for GEO downloads
dir.create(geodir, showWarnings = FALSE, recursive = TRUE)

series_ids <- c("GSE5629", "GSE5630", "GSE5631", "GSE5632", "GSE5633", "GSE5634")
geo_data <- list()

for (sid in series_ids) {
  cat("Downloading", sid, "...\n")
  geo_data[[sid]] <- tryCatch({
    gse <- getGEO(sid, destdir = geodir, getGPL = TRUE)
    if (is.list(gse)) gse <- gse[[1]]
    cat("  ", ncol(gse), "samples,", nrow(gse), "features\n")
    gse
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  })
}

# ── 3. EXTRACT TISSUE LABELS AND BUILD HIERARCHICAL MAPPING ────────────────────
#    Each AtGenExpress sample has a tissue annotation in its metadata.
#    We map these to a hierarchical structure:
#
#      Broad Organ          Sub-Tissues
#      ---------            -----------
#      Root                 root
#      Seedling             seedling_green_parts, whole_plant_pre_bolting
#      Leaf_Shoot           cotyledon, hypocotyl, rosette_leaf, cauline_leaf,
#                           senescing_leaf, rosette_vegetative, stem,
#                           shoot_apex_vegetative, shoot_apex_transition, ...
#      Flower               sepal, petal, stamen, carpel, pollen, pedicel,
#                           flower_stage_9, flower_stage_12, flower_stage_15, ...
#      Seed_Silique         silique_stage_3/4/5, seed_stage_6/7/8/9/10

# Extract tissue from characteristics_ch1 columns
extract_tissue <- function(pdata) {
  char_cols <- grep("characteristics_ch1", colnames(pdata), value = TRUE)
  tissue <- rep(NA_character_, nrow(pdata))
  for (cc in char_cols) {
    vals <- pdata[[cc]]
    tissue_matches <- grepl("Tissue:", vals)
    if (any(tissue_matches)) {
      extracted <- sub(".*Tissue:\\s*", "", vals)
      extracted <- sub("\\s*;.*", "", extracted)
      extracted[!tissue_matches] <- NA
      need <- is.na(tissue)
      tissue[need] <- extracted[need]
    }
  }
  return(tissue)
}

# Build combined sample metadata
all_samples <- data.frame(
  sample_id = character(), series = character(),
  tissue_raw = character(), stringsAsFactors = FALSE
)

for (sid in names(geo_data)) {
  gse <- geo_data[[sid]]
  if (is.null(gse)) next
  pdata <- pData(gse)
  df <- data.frame(
    sample_id = pdata$title, series = sid,
    tissue_raw = extract_tissue(pdata), stringsAsFactors = FALSE
  )
  all_samples <- rbind(all_samples, df)
}

# Map raw tissue labels to (broad_organ, sub_tissue)
build_tissue_map <- function(raw_tissue) {
  t <- tolower(trimws(raw_tissue))

  # ROOT
  if (grepl("root", t)) return(list(broad = "Root", sub = "root"))

  # SEEDLING
  if (grepl("seedling", t)) return(list(broad = "Seedling", sub = "seedling_green_parts"))
  if (grepl("developmental drift", t)) return(list(broad = "Seedling", sub = "whole_plant_pre_bolting"))

  # LEAF / SHOOT
  if (grepl("cotyledon", t)) return(list(broad = "Leaf_Shoot", sub = "cotyledon"))
  if (grepl("hypocotyl", t)) return(list(broad = "Leaf_Shoot", sub = "hypocotyl"))
  if (grepl("senescing", t)) return(list(broad = "Leaf_Shoot", sub = "senescing_leaf"))
  if (grepl("rosette leaf", t)) return(list(broad = "Leaf_Shoot", sub = "rosette_leaf"))
  if (grepl("^leaf", t) || grepl("leaves 1", t)) return(list(broad = "Leaf_Shoot", sub = "rosette_leaf"))
  if (grepl("cauline", t)) return(list(broad = "Leaf_Shoot", sub = "cauline_leaf"))
  if (grepl("veg rosette", t)) return(list(broad = "Leaf_Shoot", sub = "rosette_vegetative"))
  if (grepl("shoot apex, vegetative", t) && !grepl("inflorescence", t)) {
    if (grepl("young leaves", t)) return(list(broad = "Leaf_Shoot", sub = "shoot_apex_vegetative_with_leaves"))
    return(list(broad = "Leaf_Shoot", sub = "shoot_apex_vegetative"))
  }
  if (grepl("shoot apex, transition", t)) return(list(broad = "Leaf_Shoot", sub = "shoot_apex_transition"))
  if (grepl("stem", t) || grepl("1st node", t)) return(list(broad = "Leaf_Shoot", sub = "stem"))

  # FLOWER
  if (grepl("pollen", t)) return(list(broad = "Flower", sub = "mature_pollen"))
  if (grepl("carpel", t)) return(list(broad = "Flower", sub = "carpel"))
  if (grepl("petal", t)) return(list(broad = "Flower", sub = "petal"))
  if (grepl("sepal", t)) return(list(broad = "Flower", sub = "sepal"))
  if (grepl("stamen", t)) return(list(broad = "Flower", sub = "stamen"))
  if (grepl("pedicel", t)) return(list(broad = "Flower", sub = "pedicel"))
  if (grepl("flower stage 12 equivalent", t) || grepl("^flower$", t))
    return(list(broad = "Flower", sub = "flower_whole"))
  if (grepl("flowers stage 9", t)) return(list(broad = "Flower", sub = "flower_stage_9"))
  if (grepl("flowers stage 10", t)) return(list(broad = "Flower", sub = "flower_stage_10_11"))
  if (grepl("flowers stage 12", t)) return(list(broad = "Flower", sub = "flower_stage_12"))
  if (grepl("flowers stage 15", t)) return(list(broad = "Flower", sub = "flower_stage_15"))
  if (grepl("shoot apex, inflorescence", t)) return(list(broad = "Flower", sub = "inflorescence_apex"))

  # SEED / SILIQUE — check "seeds," before "silique" (label contains both words)
  if (grepl("^seeds?,", t)) {
    if (grepl("stage 6", t)) return(list(broad = "Seed_Silique", sub = "seed_stage_6"))
    if (grepl("stage 7", t)) return(list(broad = "Seed_Silique", sub = "seed_stage_7"))
    if (grepl("stage 8", t)) return(list(broad = "Seed_Silique", sub = "seed_stage_8"))
    if (grepl("stage 9", t)) return(list(broad = "Seed_Silique", sub = "seed_stage_9"))
    if (grepl("stage 10", t)) return(list(broad = "Seed_Silique", sub = "seed_stage_10"))
    return(list(broad = "Seed_Silique", sub = "seed"))
  }
  if (grepl("silique", t)) {
    if (grepl("stage 3", t)) return(list(broad = "Seed_Silique", sub = "silique_stage_3"))
    if (grepl("stage 4", t)) return(list(broad = "Seed_Silique", sub = "silique_stage_4"))
    if (grepl("stage 5", t)) return(list(broad = "Seed_Silique", sub = "silique_stage_5"))
    return(list(broad = "Seed_Silique", sub = "silique"))
  }

  return(list(broad = "Other", sub = "unknown"))
}

# Apply mapping
all_samples$broad_organ <- sapply(all_samples$tissue_raw, function(x) build_tissue_map(x)$broad)
all_samples$sub_tissue <- sapply(all_samples$tissue_raw, function(x) build_tissue_map(x)$sub)

cat("\nTissue mapping complete:\n")
print(table(all_samples$broad_organ))

# ── 4. COMBINE EXPRESSION MATRICES AND MAP PROBES TO AGI GENE IDs ──────────────
#    The ATH1 microarray uses probe set IDs (e.g., 244901_at).
#    We map these to Arabidopsis gene IDs (AGI format: AT{1-5}Gnnnnn)
#    using the ath1121501.db annotation package.

# Combine expression from all series (same platform = same probe order)
combined_expr <- NULL
for (sid in names(geo_data)) {
  gse <- geo_data[[sid]]
  if (is.null(gse)) next
  expr <- exprs(gse)
  if (is.null(combined_expr)) {
    combined_expr <- expr
  } else {
    combined_expr <- cbind(combined_expr, expr)
  }
}
colnames(combined_expr) <- all_samples$sample_id
cat("\nCombined expression:", nrow(combined_expr), "probes x", ncol(combined_expr), "samples\n")

# Map probes to AGI gene IDs via the TAIR column
probe_ids <- rownames(combined_expr)
agi_map <- AnnotationDbi::select(ath1121501.db,
  keys = probe_ids, columns = c("PROBEID", "TAIR"), keytype = "PROBEID")

# Keep only nuclear AGI IDs, deduplicate probes per gene
agi_map <- agi_map[!is.na(agi_map$TAIR), ]
agi_map <- agi_map[grepl("^AT[1-5]G[0-9]{5}", agi_map$TAIR), ]
agi_map <- agi_map[!duplicated(agi_map$PROBEID), ]

# Subset and rename expression matrix
expr_mapped <- combined_expr[agi_map$PROBEID, , drop = FALSE]
rownames(expr_mapped) <- agi_map$TAIR

# For genes with multiple probes, keep the one with highest mean expression
gene_ids <- rownames(expr_mapped)
dup_genes <- unique(gene_ids[duplicated(gene_ids)])
if (length(dup_genes) > 0) {
  mean_expr <- rowMeans(expr_mapped, na.rm = TRUE)
  keep_rows <- seq_len(nrow(expr_mapped))
  for (g in dup_genes) {
    idx <- which(gene_ids == g)
    best <- idx[which.max(mean_expr[idx])]
    keep_rows <- setdiff(keep_rows, setdiff(idx, best))
  }
  expr_mapped <- expr_mapped[keep_rows, , drop = FALSE]
}
cat("Gene-level expression:", nrow(expr_mapped), "genes\n")

# ── 5. BUILD TISSUE-LEVEL EXPRESSION MATRICES ─────────────────────────────────
#    Average expression across all samples within each tissue group.
#    This gives us one expression value per gene per tissue.

# Sub-tissue level (mean per sub_tissue)
sub_tissues <- unique(all_samples$sub_tissue)
sub_tissue_expr <- sapply(sub_tissues, function(st) {
  samples_st <- intersect(all_samples$sample_id[all_samples$sub_tissue == st],
                          colnames(expr_mapped))
  if (length(samples_st) > 0) rowMeans(expr_mapped[, samples_st, drop = FALSE], na.rm = TRUE)
  else rep(NA, nrow(expr_mapped))
})
rownames(sub_tissue_expr) <- rownames(expr_mapped)

# Broad-organ level (mean per broad_organ)
broad_organs <- unique(all_samples$broad_organ)
broad_organ_expr <- sapply(broad_organs, function(bo) {
  samples_bo <- intersect(all_samples$sample_id[all_samples$broad_organ == bo],
                          colnames(expr_mapped))
  if (length(samples_bo) > 0) rowMeans(expr_mapped[, samples_bo, drop = FALSE], na.rm = TRUE)
  else rep(NA, nrow(expr_mapped))
})
rownames(broad_organ_expr) <- rownames(expr_mapped)

cat("Sub-tissue matrix:", nrow(sub_tissue_expr), "genes x", ncol(sub_tissue_expr), "tissues\n")
cat("Broad-organ matrix:", nrow(broad_organ_expr), "genes x", ncol(broad_organ_expr), "organs\n")

# ── 6. COMPUTE THE TAU TISSUE-SPECIFICITY INDEX ───────────────────────────────
#    The Tau index quantifies how tissue-specific a gene's expression is.
#
#    Formula:  Tau = sum(1 - x_i / max(x)) / (n - 1)
#
#    Where:
#      x_i = expression of the gene in tissue i
#      n   = number of tissues
#
#    Interpretation:
#      Tau = 0  → gene is expressed equally in all tissues (ubiquitous/constitutive)
#      Tau = 1  → gene is expressed in only one tissue (perfectly tissue-specific)
#      Tau >= 0.6 → considered tissue-specific (standard threshold)
#      Tau < 0.6  → considered constitutive (broadly expressed)

compute_tau <- function(expr_matrix) {
  expr_matrix[expr_matrix < 0] <- 0  # clamp negatives to zero
  apply(expr_matrix, 1, function(x) {
    mx <- max(x, na.rm = TRUE)
    if (mx == 0 || is.na(mx)) return(NA)
    n <- sum(!is.na(x))
    if (n <= 1) return(NA)
    sum(1 - (x / mx), na.rm = TRUE) / (n - 1)
  })
}

tau_subtissue <- compute_tau(sub_tissue_expr)
tau_broad <- compute_tau(broad_organ_expr)

# Assign each gene to its predominant tissue (highest expression)
predominant_subtissue <- colnames(sub_tissue_expr)[apply(sub_tissue_expr, 1, function(x) {
  if (all(is.na(x)) || max(x, na.rm = TRUE) == 0) return(NA)
  which.max(x)
})]

# Classify: tissue-specific vs constitutive vs unannotated
specificity <- ifelse(tau_subtissue >= TAU_THRESHOLD, "tissue_specific", "constitutive")
specificity[is.na(tau_subtissue)] <- "unannotated"

# Map sub-tissue to broad organ
sub_to_broad <- unique(all_samples[, c("sub_tissue", "broad_organ")])
predominant_broad <- sub_to_broad$broad_organ[
  match(predominant_subtissue, sub_to_broad$sub_tissue)]

# Build gene annotation table
gene_annotation <- data.frame(
  gene_id = rownames(sub_tissue_expr),
  tau_subtissue = tau_subtissue[rownames(sub_tissue_expr)],
  tau_broad = tau_broad[rownames(sub_tissue_expr)],
  predominant_subtissue = predominant_subtissue,
  predominant_broad_organ = predominant_broad,
  specificity = specificity,
  stringsAsFactors = FALSE
)

# ── 7. MERGE WITH DEG DATA AND ANNOTATE ────────────────────────────────────────
#    Combine the DESeq2 results with the tissue annotations.

deg_annotated <- merge(deg, gene_annotation, by = "gene_id", all.x = TRUE)

# Flag genes not in the atlas as "unannotated"
deg_annotated$specificity[is.na(deg_annotated$specificity)] <- "unannotated"
deg_annotated$predominant_broad_organ[is.na(deg_annotated$predominant_broad_organ)] <- "unannotated"
deg_annotated$predominant_subtissue[is.na(deg_annotated$predominant_subtissue)] <- "unannotated"

# Reorder columns
col_order <- c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat",
               "pvalue", "padj", "deg_flag", "regulation",
               "tau_subtissue", "tau_broad", "specificity",
               "predominant_broad_organ", "predominant_subtissue")
deg_annotated <- deg_annotated[, col_order]

degs_only <- deg_annotated[deg_annotated$deg_flag == "yes", ]
cat("\n=== DEG Tissue Assignment ===\n")
cat("Tissue-specific:", sum(degs_only$specificity == "tissue_specific"), "\n")
cat("Constitutive:", sum(degs_only$specificity == "constitutive"), "\n")
cat("Unannotated:", sum(degs_only$specificity == "unannotated"), "\n")

# ── 8. EXPORT TISSUE-SPECIFIC CSV FILES (HIERARCHICAL) ─────────────────────────
#    Creates a folder structure:
#      tissue_specific_degs/
#        ├── Flower/
#        │   ├── Flower_all_DEGs.csv
#        │   ├── carpel_DEGs.csv
#        │   ├── petal_DEGs.csv
#        │   └── ...
#        ├── Leaf_Shoot/
#        ├── Root/
#        ├── Seed_Silique/
#        ├── Seedling/
#        ├── constitutive_DEGs.csv
#        ├── unannotated_DEGs.csv
#        └── all_degs_with_tissue_annotation.csv

out_base <- "tissue_specific_degs"
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

tissue_specific <- degs_only[degs_only$specificity == "tissue_specific", ]

for (bo in sort(unique(tissue_specific$predominant_broad_organ))) {
  organ_dir <- file.path(out_base, bo)
  dir.create(organ_dir, showWarnings = FALSE, recursive = TRUE)
  organ_degs <- tissue_specific[tissue_specific$predominant_broad_organ == bo, ]
  write.csv(organ_degs, file.path(organ_dir, paste0(bo, "_all_DEGs.csv")), row.names = FALSE)
  for (st in sort(unique(organ_degs$predominant_subtissue))) {
    st_degs <- organ_degs[organ_degs$predominant_subtissue == st, ]
    write.csv(st_degs, file.path(organ_dir, paste0(st, "_DEGs.csv")), row.names = FALSE)
  }
}

write.csv(degs_only[degs_only$specificity == "constitutive", ],
          file.path(out_base, "constitutive_DEGs.csv"), row.names = FALSE)
write.csv(degs_only[degs_only$specificity == "unannotated", ],
          file.path(out_base, "unannotated_DEGs.csv"), row.names = FALSE)
write.csv(deg_annotated, file.path(out_base, "all_degs_with_tissue_annotation.csv"), row.names = FALSE)

cat("CSV files exported to", out_base, "/\n")

# ── 9. GENERATE SUMMARY TABLE ─────────────────────────────────────────────────
#    tissue_deg_summary.csv: DEG counts, up/down breakdown, and statistics
#    per tissue in a hierarchical layout.

hierarchical_summary <- data.frame(
  broad_organ = character(), sub_tissue = character(),
  total_DEGs = integer(), upregulated = integer(), downregulated = integer(),
  median_log2FC = numeric(), mean_log2FC = numeric(), mean_tau_subtissue = numeric(),
  stringsAsFactors = FALSE
)

for (bo in sort(unique(tissue_specific$predominant_broad_organ))) {
  for (st in sort(unique(tissue_specific$predominant_subtissue[
    tissue_specific$predominant_broad_organ == bo]))) {
    subset_df <- tissue_specific[
      tissue_specific$predominant_broad_organ == bo &
      tissue_specific$predominant_subtissue == st, ]
    hierarchical_summary <- rbind(hierarchical_summary, data.frame(
      broad_organ = bo, sub_tissue = st,
      total_DEGs = nrow(subset_df),
      upregulated = sum(subset_df$regulation == "up"),
      downregulated = sum(subset_df$regulation == "down"),
      median_log2FC = round(median(subset_df$log2FoldChange, na.rm = TRUE), 4),
      mean_log2FC = round(mean(subset_df$log2FoldChange, na.rm = TRUE), 4),
      mean_tau_subtissue = round(mean(subset_df$tau_subtissue, na.rm = TRUE), 4),
      stringsAsFactors = FALSE
    ))
  }
}

# Add constitutive, unannotated, and total rows
for (spec in c("constitutive", "unannotated")) {
  subset_df <- degs_only[degs_only$specificity == spec, ]
  hierarchical_summary <- rbind(hierarchical_summary, data.frame(
    broad_organ = spec, sub_tissue = "(all)",
    total_DEGs = nrow(subset_df),
    upregulated = sum(subset_df$regulation == "up"),
    downregulated = sum(subset_df$regulation == "down"),
    median_log2FC = round(median(subset_df$log2FoldChange, na.rm = TRUE), 4),
    mean_log2FC = round(mean(subset_df$log2FoldChange, na.rm = TRUE), 4),
    mean_tau_subtissue = round(mean(subset_df$tau_subtissue, na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  ))
}

hierarchical_summary <- rbind(hierarchical_summary, data.frame(
  broad_organ = "TOTAL", sub_tissue = "(all)",
  total_DEGs = nrow(degs_only),
  upregulated = sum(degs_only$regulation == "up"),
  downregulated = sum(degs_only$regulation == "down"),
  median_log2FC = round(median(degs_only$log2FoldChange, na.rm = TRUE), 4),
  mean_log2FC = round(mean(degs_only$log2FoldChange, na.rm = TRUE), 4),
  mean_tau_subtissue = round(mean(degs_only$tau_subtissue, na.rm = TRUE), 4),
  stringsAsFactors = FALSE
))

write.csv(hierarchical_summary, "tissue_deg_summary.csv", row.names = FALSE)
cat("Summary table saved to tissue_deg_summary.csv\n")

# ── 10. VISUALIZATION 1: BAR CHART (DEG counts per broad organ) ────────────────
#    Stacked bar chart showing upregulated vs downregulated DEGs per organ.

tissue_degs <- degs_only[degs_only$specificity == "tissue_specific", ]
organ_order <- names(sort(table(tissue_degs$predominant_broad_organ), decreasing = TRUE))
tissue_degs$predominant_broad_organ <- factor(tissue_degs$predominant_broad_organ, levels = organ_order)

p_bar <- ggplot(tissue_degs, aes(x = predominant_broad_organ, fill = regulation)) +
  geom_bar(position = "stack", width = 0.7) +
  scale_fill_manual(values = c("up" = UP_COLOR, "down" = DOWN_COLOR),
                    labels = c("Upregulated", "Downregulated"), name = "Regulation") +
  labs(title = "Radiation-Responsive DEGs by Tissue (Broad Organ)",
       subtitle = "Arabidopsis thaliana — OSD-498/510 radiation effect",
       x = "Broad Organ", y = "Number of DEGs") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 11, color = "grey40"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
        legend.position = "top", panel.grid.minor = element_blank()) +
  geom_text(stat = "count", aes(label = after_stat(count)),
            position = position_stack(vjust = 0.5), size = 3.5, color = "white")

ggsave("fig1_deg_counts_by_organ.png", p_bar, width = 8, height = 6, dpi = 300, bg = "white")
cat("Saved fig1_deg_counts_by_organ.png\n")

# ── 11. VISUALIZATION 2: VOLCANO PLOTS (one per broad organ) ───────────────────
#    Each panel shows log2FoldChange vs -log10(padj) for genes assigned to that organ.
#    Constitutive genes are shown as grey background for context.

all_genes <- deg_annotated
all_genes$neg_log10_padj <- -log10(all_genes$padj)
all_genes$neg_log10_padj[all_genes$neg_log10_padj > 50] <- 50  # cap for visualization

volcano_list <- list()
for (bo in organ_order) {
  organ_genes <- all_genes[all_genes$predominant_broad_organ == bo, ]
  constitutive <- all_genes[all_genes$specificity == "constitutive", ]

  plot_data <- rbind(
    data.frame(constitutive[, c("log2FoldChange", "neg_log10_padj", "deg_flag", "regulation")],
               category = "constitutive"),
    data.frame(organ_genes[, c("log2FoldChange", "neg_log10_padj", "deg_flag", "regulation")],
               category = bo)
  )

  p <- ggplot(plot_data, aes(x = log2FoldChange, y = neg_log10_padj)) +
    geom_point(data = subset(plot_data, category == "constitutive"),
               color = "grey80", alpha = 0.3, size = 0.8) +
    geom_point(data = subset(plot_data, category == bo & deg_flag == "no"),
               color = "grey60", alpha = 0.4, size = 0.8) +
    geom_point(data = subset(plot_data, category == bo & deg_flag == "yes" & regulation == "up"),
               color = UP_COLOR, alpha = 0.7, size = 1.2) +
    geom_point(data = subset(plot_data, category == bo & deg_flag == "yes" & regulation == "down"),
               color = DOWN_COLOR, alpha = 0.7, size = 1.2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey50", linewidth = 0.4) +
    labs(title = bo, x = "log2 Fold Change", y = "-log10(padj)") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          legend.position = "none", panel.grid.minor = element_blank()) +
    coord_cartesian(xlim = c(-4, 7), ylim = c(0, 52))

  volcano_list[[bo]] <- p
}

ggsave("fig2_volcano_plots_by_organ.png",
       arrangeGrob(grobs = volcano_list, ncol = 3,
         top = textGrob("Volcano Plots by Tissue: Radiation-Responsive Genes",
                        gp = gpar(fontface = "bold", fontsize = 14))),
       width = 14, height = 10, dpi = 300, bg = "white")
cat("Saved fig2_volcano_plots_by_organ.png\n")

# ── 12. VISUALIZATION 3: HEATMAP (top DEGs x tissues) ─────────────────────────
#    Shows the tissue expression specificity of the 50 most significant DEGs.
#    Rows = genes, columns = sub-tissues, color = z-scored expression.

top_degs <- degs_only[order(degs_only$padj), ][1:50, ]
genes_in_atlas <- intersect(top_degs$gene_id, rownames(sub_tissue_expr))
expr_top <- sub_tissue_expr[genes_in_atlas, , drop = FALSE]
expr_z <- t(scale(t(expr_top)))  # z-score per gene

# Order columns by broad organ
col_order <- sub_to_broad[order(sub_to_broad$broad_organ), "sub_tissue"]
col_order <- intersect(col_order, colnames(expr_z))
expr_z <- expr_z[, col_order]

# Column annotation (broad organ)
col_anno_data <- data.frame(
  broad_organ = sub_to_broad$broad_organ[match(col_order, sub_to_broad$sub_tissue)],
  row.names = col_order
)
organ_colors <- c("Root" = "#0072B2", "Seedling" = "#009E73", "Leaf_Shoot" = "#E69F00",
                  "Flower" = "#CC79A7", "Seed_Silique" = "#D55E00")
col_anno <- HeatmapAnnotation(Broad_Organ = col_anno_data$broad_organ,
  col = list(Broad_Organ = organ_colors))

# Row annotation (regulation direction)
row_anno_data <- data.frame(
  regulation = top_degs$regulation[match(genes_in_atlas, top_degs$gene_id)],
  row.names = genes_in_atlas
)
row_anno <- rowAnnotation(Regulation = row_anno_data$regulation,
  col = list(Regulation = c("up" = UP_COLOR, "down" = DOWN_COLOR)))

png("fig3_heatmap_top_degs.png", width = 12, height = 8, units = "in", res = 300, bg = "white")
draw(Heatmap(expr_z,
  name = "Z-score\nExpression",
  col = colorRamp2(c(-2, 0, 2), c(DOWN_COLOR, "white", UP_COLOR)),
  top_annotation = col_anno, right_annotation = row_anno,
  cluster_columns = FALSE, cluster_rows = TRUE,
  show_row_names = TRUE, show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
  column_names_rot = 45,
  column_title = "Top 50 DEGs by Significance — Tissue Expression Specificity",
  column_title_gp = gpar(fontface = "bold", fontsize = 13)))
dev.off()
cat("Saved fig3_heatmap_top_degs.png\n")

# ── 13. VISUALIZATION 4: UPSET DIAGRAM (organ overlaps) ────────────────────────
#    Shows how many DEGs are expressed above median in each combination of organs.
#    This reveals which DEGs are active across multiple tissue types.

degs_in_atlas <- degs_only[degs_only$gene_id %in% rownames(broad_organ_expr), ]
expr_degs <- broad_organ_expr[degs_in_atlas$gene_id, , drop = FALSE]
gene_medians <- apply(expr_degs, 1, median, na.rm = TRUE)

binary_matrix <- matrix(0, nrow = nrow(expr_degs), ncol = ncol(expr_degs))
rownames(binary_matrix) <- rownames(expr_degs)
colnames(binary_matrix) <- colnames(expr_degs)
for (i in seq_len(nrow(expr_degs))) {
  binary_matrix[i, ] <- as.integer(expr_degs[i, ] > gene_medians[i])
}
binary_df <- as.data.frame(binary_matrix)

png("fig4_upset_organ_overlaps.png", width = 10, height = 6, units = "in", res = 300, bg = "white")
print(upset(binary_df,
  sets = colnames(expr_degs), order.by = "freq", nsets = 5, nintersects = 20,
  sets.bar.color = c("#0072B2", "#009E73", "#E69F00", "#CC79A7", "#D55E00"),
  main.bar.color = "grey30", matrix.color = "#0279EE",
  point.size = 3, line.size = 0.8, text.scale = 1.2,
  mb.ratio = c(0.6, 0.4),
  mainbar.y.label = "DEG Intersection Size",
  sets.x.label = "DEGs Expressed Above Median"))
dev.off()
cat("Saved fig4_upset_organ_overlaps.png\n")

# ── 14. STATISTICAL ANALYSIS: Is the tissue-specific response significant? ────
#    Tests whether DEG status depends on tissue assignment.
#
#    Framework: For each broad organ, we know:
#      - N_total = number of genes tissue-specific to that organ (background)
#      - N_deg   = number of those genes that became DEGs (observed)
#
#    Tests:
#      1. Chi-square test of independence (overall tissue effect)
#      2. Fisher's exact test per tissue (each tissue vs all others)
#      3. Pairwise Fisher's exact tests (organ vs organ, e.g., root vs shoot)
#      4. Chi-square test on up/down direction (does response direction differ?)

# Build contingency table: broad_organ × deg_flag (tissue-specific genes only)
atlas_genes <- rownames(broad_organ_expr)
deg_in_atlas <- deg_annotated[deg_annotated$gene_id %in% atlas_genes, ]
tissue_specific_all <- deg_in_atlas[deg_in_atlas$specificity == "tissue_specific", ]
contingency <- table(tissue_specific_all$predominant_broad_organ,
                     tissue_specific_all$deg_flag)

cat("\n=== Contingency table: Organ × DEG status ===\n")
print(contingency)

# TEST 1: Chi-square test of independence
cat("\n=== TEST 1: Chi-square test of independence ===\n")
chi_result <- chisq.test(contingency)
print(chi_result)
cat("Residuals (positive = enriched, negative = depleted):\n")
print(round(chi_result$residuals[, "yes"], 3))

# TEST 2: Per-tissue Fisher's exact test (organ vs rest)
cat("\n=== TEST 2: Per-tissue Fisher's exact test ===\n")
organs <- rownames(contingency)
fisher_results <- data.frame(
  organ = character(), deg_in_organ = integer(), total_in_organ = integer(),
  deg_rate_organ = numeric(), deg_rate_rest = numeric(), odds_ratio = numeric(),
  p_value = numeric(), direction = character(), stringsAsFactors = FALSE
)

for (org in organs) {
  org_yes <- contingency[org, "yes"]
  org_no  <- contingency[org, "no"]
  rest_yes <- sum(contingency[, "yes"]) - org_yes
  rest_no  <- sum(contingency[, "no"]) - org_no
  mat <- matrix(c(org_yes, org_no, rest_yes, rest_no), nrow = 2)
  ft <- fisher.test(mat, alternative = "two.sided")
  fisher_results <- rbind(fisher_results, data.frame(
    organ = org, deg_in_organ = org_yes, total_in_organ = org_yes + org_no,
    deg_rate_organ = round(org_yes / (org_yes + org_no), 4),
    deg_rate_rest = round(rest_yes / (rest_yes + rest_no), 4),
    odds_ratio = round(ft$estimate, 3), p_value = ft$p.value,
    direction = ifelse(org_yes / (org_yes + org_no) > rest_yes / (rest_yes + rest_no),
                       "enriched", "depleted"), stringsAsFactors = FALSE
  ))
}
fisher_results$padj <- p.adjust(fisher_results$p_value, method = "bonferroni")
fisher_results$significant <- ifelse(fisher_results$padj < 0.05, "YES", "no")
print(fisher_results, row.names = FALSE)

write.csv(fisher_results, "tissue_enrichment_statistics.csv", row.names = FALSE)

# TEST 3: Pairwise Fisher's exact tests (organ vs organ)
cat("\n=== TEST 3: Pairwise Fisher's exact tests ===\n")
pairwise_results <- data.frame(
  organ1 = character(), organ2 = character(), rate1 = numeric(), rate2 = numeric(),
  odds_ratio = numeric(), p_value = numeric(), stringsAsFactors = FALSE
)

for (i in 1:(length(organs) - 1)) {
  for (j in (i + 1):length(organs)) {
    mat <- matrix(c(
      contingency[organs[i], "yes"], contingency[organs[i], "no"],
      contingency[organs[j], "yes"], contingency[organs[j], "no"]
    ), nrow = 2)
    ft <- fisher.test(mat, alternative = "two.sided")
    pairwise_results <- rbind(pairwise_results, data.frame(
      organ1 = organs[i], organ2 = organs[j],
      rate1 = round(contingency[organs[i], "yes"] / sum(contingency[organs[i], ]), 4),
      rate2 = round(contingency[organs[j], "yes"] / sum(contingency[organs[j], ]), 4),
      odds_ratio = round(ft$estimate, 3), p_value = ft$p.value,
      stringsAsFactors = FALSE
    ))
  }
}
pairwise_results$padj <- p.adjust(pairwise_results$p_value, method = "bonferroni")
pairwise_results$higher <- ifelse(pairwise_results$rate1 > pairwise_results$rate2,
  paste0(pairwise_results$organ1, " > ", pairwise_results$organ2),
  paste0(pairwise_results$organ2, " > ", pairwise_results$organ1))
pairwise_results$significant <- ifelse(pairwise_results$padj < 0.05, "YES", "no")
pairwise_results <- pairwise_results[order(pairwise_results$p_value), ]
print(pairwise_results, row.names = FALSE)

write.csv(pairwise_results, "pairwise_tissue_comparisons.csv", row.names = FALSE)

# TEST 4: Direction of response (up/down ratio across tissues)
cat("\n=== TEST 4: Up/down direction differs by tissue? ===\n")
ts_degs <- tissue_specific_all[tissue_specific_all$deg_flag == "yes", ]
direction_table <- table(ts_degs$predominant_broad_organ, ts_degs$regulation)
print(direction_table)
chi_dir <- chisq.test(direction_table)
print(chi_dir)

# ── 15. VISUALIZATION 5: Forest plot of tissue enrichment ─────────────────────
#    Shows DEG rate per organ with 95% CI and significance markers.

overall_rate <- sum(contingency[, "yes"]) / sum(contingency)

plot_data <- data.frame(
  organ = factor(fisher_results$organ, levels = fisher_results$organ[order(fisher_results$deg_rate_organ)]),
  deg_rate = fisher_results$deg_rate_organ,
  direction = fisher_results$direction,
  ci_lower = sapply(1:nrow(fisher_results), function(i)
    binom.test(fisher_results$deg_in_organ[i], fisher_results$total_in_organ[i])$conf.int[1]),
  ci_upper = sapply(1:nrow(fisher_results), function(i)
    binom.test(fisher_results$deg_in_organ[i], fisher_results$total_in_organ[i])$conf.int[2]),
  sig_label = ifelse(fisher_results$padj < 0.001, "***",
              ifelse(fisher_results$padj < 0.01, "**",
              ifelse(fisher_results$padj < 0.05, "*", "ns"))),
  stringsAsFactors = FALSE
)

p_forest <- ggplot(plot_data, aes(x = deg_rate, y = organ, color = direction)) +
  geom_vline(xintercept = overall_rate, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, linewidth = 0.8) +
  geom_point(size = 4) +
  geom_text(aes(x = deg_rate, y = organ, label = sprintf("%.1f%% %s", deg_rate * 100, sig_label)),
            hjust = -0.2, vjust = 0.5, size = 3.5, color = "black") +
  scale_color_manual(values = c("enriched" = "#0072B2", "depleted" = "#D55E00"),
                     labels = c("Enriched for DEGs", "Depleted for DEGs"), name = "") +
  scale_x_continuous(limits = c(0.15, 0.70)) +
  labs(title = "Tissue-Specific Radiation Response: DEG Rate by Organ",
       subtitle = "Fisher's exact test (organ vs rest), Bonferroni-corrected  |  *** p<0.001, ** p<0.01, * p<0.05",
       x = "Proportion of tissue-specific genes that are DEGs", y = "Broad Organ") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 10, color = "grey40"),
        legend.position = "top", panel.grid.minor = element_blank(),
        axis.text.y = element_text(face = "bold"))

ggsave("fig5_tissue_enrichment_test.png", p_forest, width = 9, height = 5, dpi = 300, bg = "white")
cat("Saved fig5_tissue_enrichment_test.png\n")

# ── DONE ───────────────────────────────────────────────────────────────────────
cat("\n=== Analysis complete! ===\n")
cat("Outputs:\n")
cat("  1. tissue_specific_degs/  (hierarchical CSV files)\n")
cat("  2. tissue_deg_summary.csv (summary table)\n")
cat("  3. fig1_deg_counts_by_organ.png\n")
cat("  4. fig2_volcano_plots_by_organ.png\n")
cat("  5. fig3_heatmap_top_degs.png\n")
cat("  6. fig4_upset_organ_overlaps.png\n")
cat("  7. fig5_tissue_enrichment_test.png (statistical test)\n")
cat("  8. tissue_enrichment_statistics.csv (per-tissue test results)\n")
cat("  9. pairwise_tissue_comparisons.csv (pairwise test results)\n")
