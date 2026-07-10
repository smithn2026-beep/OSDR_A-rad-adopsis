# ============================================================================
# Tissue-Specific DEG Analysis v2 — Corrected (Unnormalized) DESeq2 Data
# Arabidopsis Ionizing Radiation: 6 Comparisons from NASA GeneLab OSD Studies
#
# FULLY SELF-CONTAINED: Downloads AtGenExpress atlas from GEO, builds tissue
# expression reference, maps probes to genes, computes Tau index, then runs
# the complete tissue-specific DEG pipeline on 6 DESeq2 output tables.
#
# Input: 6 DESeq2 output CSVs (re-run with raw unnormalized counts)
# Method: AtGenExpress developmental atlas (Schmid et al. 2005) + Tau index
# Threshold: padj < 0.05 AND |log2FC| >= 1
#
# Comparisons:
#   1. radiation_effect       — OSD498+510, Co-60 100Gy vs none (WT)
#   2. genotype_interaction   — OSD508+510, sog1-1 x radiation interaction
#   3. GCR40                  — OSD658, simulated GCR 40cGy vs none (WT)
#   4. GCR80                  — OSD658, simulated GCR 80cGy vs none (WT)
#   5. Cs137_100cGy           — OSD782, Cs-137 100cGy vs none (WT)
#   6. Cs137_10cGy            — OSD782, Cs-137 10cGy vs none (WT)
#
# Requirements: R >= 4.0, with packages:
#   GEOquery, ath1121501.db, AnnotationDbi, Biobase (Bioconductor)
#   data.table, ggplot2, ComplexHeatmap, circlize, UpSetR, grid, gridExtra,
#   pheatmap (CRAN/Bioconductor)
# ============================================================================

# ==== Section 0: Package Installation & Configuration ======================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_pkgs <- c("GEOquery", "ath1121501.db", "AnnotationDbi", "Biobase",
               "ComplexHeatmap", "circlize")
cran_pkgs <- c("data.table", "ggplot2", "UpSetR", "gridExtra", "pheatmap")

for (p in c(bioc_pkgs, cran_pkgs)) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p %in% bioc_pkgs) BiocManager::install(p, update = FALSE, ask = FALSE)
    else install.packages(p)
  }
}

suppressPackageStartupMessages({
  library(GEOquery)
  library(ath1121501.db)
  library(AnnotationDbi)
  library(Biobase)
  library(data.table)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(UpSetR)
  library(grid)
  library(gridExtra)
  library(pheatmap)
})

# ---- Configuration ----
TAU_THRESHOLD  <- 0.6
PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD  <- 1.0

UP_COLOR   <- "#0072B2"
DOWN_COLOR <- "#D55E00"
ORGAN_COLORS <- c(
  Root = "#0072B2", Seedling = "#009E73", Leaf_Shoot = "#E69F00",
  Flower = "#CC79A7", Seed_Silique = "#D55E00"
)
BROAD_ORGANS <- c("Root", "Seedling", "Leaf_Shoot", "Flower", "Seed_Silique")

# ---- Paths (adjust to your local setup) ----
RESULTS_DIR <- "/mnt/results"
GEO_CACHE   <- "/workspace/geo_cache"
ATLAS_PATH  <- "/workspace/atgenexpress_atlas.rds"
INPUT_DIR   <- "/mnt/user-uploads"
OLD_FILE    <- file.path(INPUT_DIR, "DEG_OSD498_510_radiation_effect.csv")

COMPARISON_FILES <- list(
  list(label = "radiation_effect",     file = file.path(INPUT_DIR, "DEG_unnorm_OSD498_510_radiation_effect.csv")),
  list(label = "genotype_interaction", file = file.path(INPUT_DIR, "DEG_unnorm_OSD508_510_genotype_x_radiation_interaction.csv")),
  list(label = "GCR40",                file = file.path(INPUT_DIR, "DEG_unnorm_OSD658_GCR40_vs_none.csv")),
  list(label = "GCR80",                file = file.path(INPUT_DIR, "DEG_unnorm_OSD658_GCR80_vs_none.csv")),
  list(label = "Cs137_100cGy",         file = file.path(INPUT_DIR, "DEG_unnorm_OSD782_100cGy_vs_none.csv")),
  list(label = "Cs137_10cGy",          file = file.path(INPUT_DIR, "DEG_unnorm_OSD782_10cGy_vs_none.csv"))
)

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(GEO_CACHE,   showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Section 1: Download AtGenExpress Developmental Atlas from GEO
# ============================================================================
# Schmid et al. 2005 (Nat Genet 37:501-6). 237 samples, 33 sub-tissues,
# ATH1 Affymetrix platform (GPL198), 22,810 probes.
#   GSE5629 = seedlings/whole plants (24), GSE5630 = leaves (60),
#   GSE5631 = roots (21), GSE5632 = flowers/pollen (66),
#   GSE5633 = shoots/stems (42), GSE5634 = siliques/seeds (24)

GEO_SERIES <- c("GSE5629", "GSE5630", "GSE5631", "GSE5632", "GSE5633", "GSE5634")

# Tissue mapping: raw GEO labels -> broad organ + sub-tissue
# CRITICAL: check "^seeds?," before "silique" so "seeds, stage 6, w/o siliques"
# is not misclassified as silique (it contains the word "siliques")
build_tissue_map <- function(raw_tissue) {
  t <- tolower(trimws(raw_tissue))
  if (grepl("root", t))               return(list(broad = "Root", sub = "root"))
  if (grepl("seedling", t))           return(list(broad = "Seedling", sub = "seedling_green_parts"))
  if (grepl("developmental drift", t))return(list(broad = "Seedling", sub = "whole_plant_pre_bolting"))
  # Seeds before siliques (see comment above)
  if (grepl("^seeds?,", t)) {
    if (grepl("stage 6", t))  return(list(broad = "Seed_Silique", sub = "seed_stage_6"))
    if (grepl("stage 7", t))  return(list(broad = "Seed_Silique", sub = "seed_stage_7"))
    if (grepl("stage 8", t))  return(list(broad = "Seed_Silique", sub = "seed_stage_8"))
    if (grepl("stage 9", t))  return(list(broad = "Seed_Silique", sub = "seed_stage_9"))
    if (grepl("stage 10", t)) return(list(broad = "Seed_Silique", sub = "seed_stage_10"))
  }
  if (grepl("silique", t)) {
    if (grepl("stage 3", t)) return(list(broad = "Seed_Silique", sub = "silique_stage_3"))
    if (grepl("stage 4", t)) return(list(broad = "Seed_Silique", sub = "silique_stage_4"))
    if (grepl("stage 5", t)) return(list(broad = "Seed_Silique", sub = "silique_stage_5"))
  }
  if (grepl("cotyledon", t))      return(list(broad = "Leaf_Shoot", sub = "cotyledon"))
  if (grepl("hypocotyl", t))      return(list(broad = "Leaf_Shoot", sub = "hypocotyl"))
  if (grepl("cauline", t))        return(list(broad = "Leaf_Shoot", sub = "cauline_leaf"))
  if (grepl("senescing", t))      return(list(broad = "Leaf_Shoot", sub = "senescing_leaf"))
  if (grepl("rosette leaf", t) || grepl("^leaf", t) || grepl("leaves 1", t) || grepl("leaf 7", t))
                                  return(list(broad = "Leaf_Shoot", sub = "rosette_leaf"))
  if (grepl("veg rosette", t))    return(list(broad = "Leaf_Shoot", sub = "rosette_vegetative"))
  if (grepl("shoot apex.*transition", t))
                                  return(list(broad = "Leaf_Shoot", sub = "shoot_apex_transition"))
  if (grepl("shoot apex.*vegetative.*young leaves", t))
                                  return(list(broad = "Leaf_Shoot", sub = "shoot_apex_vegetative_with_leaves"))
  if (grepl("shoot apex.*vegetative", t))
                                  return(list(broad = "Leaf_Shoot", sub = "shoot_apex_vegetative"))
  if (grepl("stem", t) || grepl("1st node", t))
                                  return(list(broad = "Leaf_Shoot", sub = "stem"))
  if (grepl("mature pollen", t))  return(list(broad = "Flower", sub = "mature_pollen"))
  if (grepl("inflorescence", t))  return(list(broad = "Flower", sub = "inflorescence_apex"))
  if (grepl("pedicel", t))        return(list(broad = "Flower", sub = "pedicel"))
  if (grepl("carpel", t))         return(list(broad = "Flower", sub = "carpel"))
  if (grepl("sepal", t))          return(list(broad = "Flower", sub = "sepal"))
  if (grepl("petal", t))          return(list(broad = "Flower", sub = "petal"))
  if (grepl("stamen", t))         return(list(broad = "Flower", sub = "stamen"))
  if (grepl("stage 10", t))       return(list(broad = "Flower", sub = "flower_stage_10_11"))
  if (grepl("stage 12", t) && !grepl("equivalent", t))
                                  return(list(broad = "Flower", sub = "flower_stage_12"))
  if (grepl("stage 15", t))       return(list(broad = "Flower", sub = "flower_stage_15"))
  if (grepl("stage 9", t))        return(list(broad = "Flower", sub = "flower_stage_9"))
  if (grepl("flower", t))         return(list(broad = "Flower", sub = "flower_whole"))
  return(list(broad = NA, sub = NA))
}

cat("=== Section 1: Downloading AtGenExpress atlas from GEO ===\n")
geo_data <- list()
all_samples <- data.frame()

for (gse_id in GEO_SERIES) {
  cat(sprintf("  %s ...", gse_id))
  gse <- getGEO(gse_id, destdir = GEO_CACHE, getGPL = TRUE)[[1]]
  geo_data[[gse_id]] <- gse
  pdata <- pData(gse)
  # Find tissue column in characteristics_ch1.* or source_name_ch1
  tissue_col <- "source_name_ch1"
  for (col in colnames(pdata)) {
    if (grepl("characteristics_ch1", col)) {
      vals <- pdata[[col]]
      if (any(grepl("tissue|root|leaf|flower|seed|silique|stem|shoot|cotyledon|hypocotyl|pollen|petal|sepal|stamen|carpel|pedicel|senescing|rosette|cauline|inflorescence|seedling|developmental", vals, ignore.case = TRUE))) {
        tissue_col <- col; break
      }
    }
  }
  tissue_labels <- gsub("^Tissue:\\s*", "", pdata[[tissue_col]], ignore.case = TRUE)
  tissue_labels <- trimws(tissue_labels)
  for (i in seq_along(sampleNames(gse))) {
    m <- build_tissue_map(tissue_labels[i])
    all_samples <- rbind(all_samples, data.frame(
      sample_id = sampleNames(gse)[i], series = gse_id,
      tissue_raw = tissue_labels[i], broad_organ = m$broad, sub_tissue = m$sub,
      stringsAsFactors = FALSE))
  }
  cat(sprintf(" %d samples\n", ncol(exprs(gse))))
}
cat(sprintf("Total: %d samples, %d sub-tissues, %d unmapped\n",
  nrow(all_samples), length(unique(all_samples$sub_tissue)), sum(is.na(all_samples$broad_organ))))

# ============================================================================
# Section 2: Combine Expression & Map Probes to AGI Gene IDs
# ============================================================================
cat("\n=== Section 2: Combining expression and mapping probes ===\n")
combined_expr <- NULL
for (sid in names(geo_data)) {
  e <- exprs(geo_data[[sid]])
  combined_expr <- if (is.null(combined_expr)) e else cbind(combined_expr, e)
}
cat(sprintf("Combined: %d probes x %d samples\n", nrow(combined_expr), ncol(combined_expr)))

agi_map <- AnnotationDbi::select(ath1121501.db,
  keys = rownames(combined_expr), columns = c("PROBEID", "TAIR", "SYMBOL", "GENENAME"),
  keytype = "PROBEID")
agi_map_clean <- agi_map[!is.na(agi_map$TAIR), ]
agi_map_clean <- agi_map_clean[grepl("^AT[1-5]G[0-9]{5}", agi_map_clean$TAIR), ]
agi_map_clean <- agi_map_clean[!duplicated(agi_map_clean$PROBEID), ]
cat(sprintf("Probes with nuclear AGI: %d, Unique genes: %d\n",
  nrow(agi_map_clean), length(unique(agi_map_clean$TAIR))))

# ============================================================================
# Section 3: Build Tissue-Level Expression Matrices
# ============================================================================
cat("\n=== Section 3: Building tissue-level expression matrices ===\n")
expr_mapped <- combined_expr[agi_map_clean$PROBEID, ]
# Deduplicate: keep probe with highest mean expression per gene
agi_map_clean$mean_expr <- rowMeans(expr_mapped, na.rm = TRUE)
agi_map_clean <- agi_map_clean[order(agi_map_clean$TAIR, -agi_map_clean$mean_expr), ]
agi_map_clean <- agi_map_clean[!duplicated(agi_map_clean$TAIR), ]
expr_dedup <- expr_mapped[agi_map_clean$PROBEID, ]
rownames(expr_dedup) <- agi_map_clean$TAIR
cat(sprintf("After dedup: %d genes\n", nrow(agi_map_clean)))

sub_tissues <- sort(unique(all_samples$sub_tissue))
sub_tissue_expr <- sapply(sub_tissues, function(st) {
  ss <- intersect(all_samples$sample_id[all_samples$sub_tissue == st], colnames(expr_dedup))
  if (length(ss) == 0) rep(NA, nrow(expr_dedup)) else rowMeans(expr_dedup[, ss, drop = FALSE], na.rm = TRUE)
})
rownames(sub_tissue_expr) <- rownames(expr_dedup)

broad_organs <- sort(unique(all_samples$broad_organ))
broad_organ_expr <- sapply(broad_organs, function(bo) {
  ss <- intersect(all_samples$sample_id[all_samples$broad_organ == bo], colnames(expr_dedup))
  if (length(ss) == 0) rep(NA, nrow(expr_dedup)) else rowMeans(expr_dedup[, ss, drop = FALSE], na.rm = TRUE)
})
rownames(broad_organ_expr) <- rownames(expr_dedup)

cat(sprintf("Sub-tissue: %d x %d, Broad-organ: %d x %d\n",
  nrow(sub_tissue_expr), ncol(sub_tissue_expr), nrow(broad_organ_expr), ncol(broad_organ_expr)))

atlas <- list(sub_tissue_expr = sub_tissue_expr, broad_organ_expr = broad_organ_expr,
              all_samples = all_samples, agi_map_clean = agi_map_clean)
saveRDS(atlas, ATLAS_PATH)
cat(sprintf("Saved atlas to %s\n", ATLAS_PATH))

# ============================================================================
# Section 4: Compute Tau Tissue-Specificity Index
# ============================================================================
# tau = sum(1 - x_i / max(x)) / (n - 1); range [0,1]; threshold 0.6
compute_tau <- function(expr_matrix) {
  apply(expr_matrix, 1, function(x) {
    x <- as.numeric(x)
    if (max(x, na.rm = TRUE) == 0) return(NA)
    sum(1 - x / max(x, na.rm = TRUE), na.rm = TRUE) / (sum(!is.na(x)) - 1)
  })
}

cat("\n=== Section 4: Computing Tau index ===\n")
tau_subtissue <- compute_tau(sub_tissue_expr)
tau_broad     <- compute_tau(broad_organ_expr)
cat(sprintf("Sub-tissue Tau: mean=%.4f median=%.4f\n",
  mean(tau_subtissue, na.rm = TRUE), median(tau_subtissue, na.rm = TRUE)))

gene_ids <- rownames(sub_tissue_expr)
tissue_ann <- data.table(gene_id = gene_ids,
  tau_subtissue = tau_subtissue[gene_ids], tau_broad = tau_broad[gene_ids],
  is_tissue_specific = tau_subtissue[gene_ids] >= TAU_THRESHOLD)
for (g in gene_ids) {
  if (tissue_ann[gene_id == g, is_tissue_specific]) {
    tissue_ann[gene_id == g, predominant_subtissue :=
      colnames(sub_tissue_expr)[which.max(sub_tissue_expr[g, ])]]
    tissue_ann[gene_id == g, predominant_broad_organ :=
      colnames(broad_organ_expr)[which.max(broad_organ_expr[g, ])]]
  }
}
cat(sprintf("Tissue-specific: %d, Constitutive: %d\n",
  sum(tissue_ann$is_tissue_specific, na.rm = TRUE),
  sum(!tissue_ann$is_tissue_specific, na.rm = TRUE)))
st_bo_map <- unique(all_samples[, c("sub_tissue", "broad_organ")])

# ============================================================================
# Section 5: Load & Validate DESeq2 CSVs
# ============================================================================
load_deg_data <- function(file_path) {
  dt <- fread(file_path)
  colnames(dt)[1] <- "gene_id"
  dt[, deg_flag := ifelse(!is.na(padj) & padj < PADJ_THRESHOLD & abs(log2FoldChange) >= LFC_THRESHOLD, "yes", "no")]
  dt[, direction := ifelse(deg_flag == "yes" & log2FoldChange > 0, "up",
                    ifelse(deg_flag == "yes" & log2FoldChange < 0, "down", "none"))]
  dt[, is_nuclear := grepl("^AT[1-5]G[0-9]{5}$", gene_id)]
  dt
}

cat("\n=== Section 5: Loading DESeq2 output tables ===\n")
deg_list <- lapply(COMPARISON_FILES, function(cf) load_deg_data(cf$file))
names(deg_list) <- sapply(COMPARISON_FILES, `[[`, "label")
for (label in names(deg_list)) {
  dt <- deg_list[[label]]
  cat(sprintf("[%s] %d genes, %d DEGs (%d up, %d down)\n", label, nrow(dt),
    sum(dt$deg_flag == "yes"), sum(dt$direction == "up"), sum(dt$direction == "down")))
}

# ============================================================================
# Section 6: Merge DEGs with Tissue Annotations
# ============================================================================
cat("\n=== Section 6: Merging with tissue annotations ===\n")
deg_annotated_list <- lapply(names(deg_list), function(label) {
  dt <- merge(deg_list[[label]], tissue_ann, by = "gene_id", all.x = TRUE)
  dt[, tissue_class := fifelse(is.na(tau_subtissue), "unannotated",
    fifelse(is_tissue_specific == TRUE, "tissue_specific", "constitutive"))]
  dt[is.na(tau_subtissue), tissue_class := fifelse(
    grepl("^AT[CM]G[0-9]{5}$", gene_id), "unannotated_organellar", "unannotated")]
  dt
})
names(deg_annotated_list) <- names(deg_list)

# ============================================================================
# Section 7: Export Hierarchical CSVs
# ============================================================================
export_cols <- c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat",
  "pvalue", "padj", "deg_flag", "direction", "tau_subtissue", "tau_broad",
  "predominant_subtissue", "predominant_broad_organ", "tissue_class")

cat("\n=== Section 7: Exporting CSVs ===\n")
for (label in names(deg_annotated_list)) {
  dt <- deg_annotated_list[[label]]
  comp_dir <- file.path(RESULTS_DIR, "tissue_specific_degs_v2", label)
  dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)
  degs_only <- dt[deg_flag == "yes"]
  fwrite(dt[, ..export_cols], file.path(comp_dir, "all_degs_with_tissue_annotation.csv"))
  fwrite(degs_only[tissue_class == "constitutive", ..export_cols], file.path(comp_dir, "constitutive_DEGs.csv"))
  fwrite(degs_only[tissue_class %like% "unannotated", ..export_cols], file.path(comp_dir, "unannotated_DEGs.csv"))
  for (bo in BROAD_ORGANS) {
    organ_dir <- file.path(comp_dir, bo)
    dir.create(organ_dir, showWarnings = FALSE, recursive = TRUE)
    organ_degs <- degs_only[tissue_class == "tissue_specific" & predominant_broad_organ == bo, ..export_cols]
    fwrite(organ_degs, file.path(organ_dir, "all_DEGs.csv"))
    for (st in as.character(st_bo_map[st_bo_map$broad_organ == bo, "sub_tissue"])) {
      fwrite(organ_degs[predominant_subtissue == st, ..export_cols], file.path(organ_dir, paste0(st, ".csv")))
    }
  }
}

# ============================================================================
# Section 8: Summary Tables
# ============================================================================
cat("\n=== Section 8: Summary tables ===\n")
for (label in names(deg_annotated_list)) {
  degs <- deg_annotated_list[[label]][deg_flag == "yes" & tissue_class == "tissue_specific"]
  rows <- list()
  for (bo in BROAD_ORGANS) {
    for (st in as.character(st_bo_map[st_bo_map$broad_organ == bo, "sub_tissue"])) {
      sd <- degs[predominant_subtissue == st]
      if (nrow(sd) > 0)
        rows[[length(rows)+1]] <- data.table(broad_organ=bo, sub_tissue=st, total_DEGs=nrow(sd),
          upregulated=sum(sd$direction=="up"), downregulated=sum(sd$direction=="down"),
          median_log2FC=round(median(sd$log2FoldChange,na.rm=TRUE),3),
          mean_log2FC=round(mean(sd$log2FoldChange,na.rm=TRUE),3),
          mean_tau_subtissue=round(mean(sd$tau_subtissue,na.rm=TRUE),4))
    }
  }
  fwrite(rbindlist(rows), file.path(RESULTS_DIR, paste0("tissue_deg_summary_", label, ".csv")))
}
cross_rows <- list()
for (label in names(deg_annotated_list)) {
  degs <- deg_annotated_list[[label]][deg_flag == "yes" & tissue_class == "tissue_specific"]
  for (bo in BROAD_ORGANS) {
    n <- tissue_ann[is_tissue_specific == TRUE & predominant_broad_organ == bo, .N]
    cross_rows[[length(cross_rows)+1]] <- data.table(broad_organ=bo, comparison=label,
      DEGs=degs[predominant_broad_organ==bo,.N],
      upregulated=degs[predominant_broad_organ==bo&direction=="up",.N],
      downregulated=degs[predominant_broad_organ==bo&direction=="down",.N],
      total_atlas_genes=n, deg_rate=round(degs[predominant_broad_organ==bo,.N]/n*100,2))
  }
}
cross_summary <- rbindlist(cross_rows)
fwrite(cross_summary, file.path(RESULTS_DIR, "cross_comparison_tissue_summary.csv"))
fwrite(dcast(cross_summary, broad_organ ~ comparison, value.var = "deg_rate"),
  file.path(RESULTS_DIR, "cross_comparison_direction.csv"))

# ============================================================================
# Section 9: Per-Comparison Statistical Tests
# ============================================================================
run_stat_tests <- function(dt, label) {
  ts_genes <- dt[tissue_class == "tissue_specific"]
  contig <- table(ts_genes$predominant_broad_organ, ts_genes$deg_flag)
  for (bo in BROAD_ORGANS) {
    if (!(bo %in% rownames(contig))) {
      contig <- rbind(contig, setNames(c(0, 0), colnames(contig)))
      rownames(contig)[nrow(contig)] <- bo
    }
  }
  contig <- contig[BROAD_ORGANS, ]
  td <- sum(contig[, "yes"]); tn <- sum(contig[, "no"])
  chi_test <- if (td >= 5) chisq.test(contig) else NULL
  fr <- list()
  for (bo in BROAD_ORGANS) {
    fm <- matrix(c(contig[bo,"yes"], contig[bo,"no"], td-contig[bo,"yes"], tn-contig[bo,"no"]), nrow=2)
    ft <- fisher.test(fm)
    fr[[bo]] <- data.table(comparison=label, broad_organ=bo, organ_DEGs=contig[bo,"yes"],
      organ_total=sum(contig[bo,]), organ_rate=round(contig[bo,"yes"]/sum(contig[bo,])*100,2),
      odds_ratio=round(unname(ft$estimate),3), ci95_low=round(ft$conf.int[1],3),
      ci95_high=round(ft$conf.int[2],3), pvalue=ft$p.value,
      padj_bonf=p.adjust(ft$p.value,method="bonferroni",n=5),
      enrichment=ifelse(unname(ft$estimate)>1&ft$p.value<0.05,"enriched",
        ifelse(unname(ft$estimate)<1&ft$p.value<0.05,"depleted","ns")))
  }
  fisher_dt <- rbindlist(fr)
  pairwise_dt <- NULL
  if (td >= 5) {
    pr <- list()
    for (p in combn(BROAD_ORGANS, 2, simplify=FALSE)) {
      m <- matrix(c(contig[p[1],"yes"],contig[p[1],"no"],contig[p[2],"yes"],contig[p[2],"no"]),nrow=2)
      ft <- fisher.test(m)
      pr[[paste(p[1],p[2],sep="_vs_")]] <- data.table(comparison=label,organ_1=p[1],organ_2=p[2],
        rate_1=round(contig[p[1],"yes"]/sum(contig[p[1],])*100,2),
        rate_2=round(contig[p[2],"yes"]/sum(contig[p[2],])*100,2),
        odds_ratio=round(unname(ft$estimate),3), pvalue=ft$p.value,
        padj_bonf=p.adjust(ft$p.value,method="bonferroni",n=10))
    }
    pairwise_dt <- rbindlist(pr)
  }
  dir_result <- NULL
  dg <- ts_genes[deg_flag == "yes"]
  if (nrow(dg) >= 5 && sum(dg$direction=="up") > 0 && sum(dg$direction=="down") > 0) {
    dc <- table(dg$predominant_broad_organ, dg$direction)
    dc <- dc[rowSums(dc) > 0, , drop=FALSE]
    if (nrow(dc) >= 2) {
      dchi <- chisq.test(dc)
      dir_result <- list(chi=dchi, pcts=data.table(comparison=label, broad_organ=rownames(dc),
        n_DEGs=rowSums(dc), pct_up=round(dc[,"up"]/rowSums(dc)*100,1),
        pct_down=round(dc[,"down"]/rowSums(dc)*100,1)))
    }
  }
  list(chi=chi_test, fisher=fisher_dt, pairwise=pairwise_dt, direction=dir_result)
}

cat("\n=== Section 9: Per-comparison statistical tests ===\n")
all_stats <- lapply(names(deg_annotated_list), function(l) run_stat_tests(deg_annotated_list[[l]], l))
names(all_stats) <- names(deg_annotated_list)
for (label in names(all_stats)) {
  fwrite(all_stats[[label]]$fisher, file.path(RESULTS_DIR, paste0("tissue_enrichment_stats_", label, ".csv")))
  if (!is.null(all_stats[[label]]$pairwise))
    fwrite(all_stats[[label]]$pairwise, file.path(RESULTS_DIR, paste0("pairwise_comparisons_", label, ".csv")))
}
fwrite(rbindlist(lapply(all_stats, function(x) x$fisher)), file.path(RESULTS_DIR, "cross_comparison_statistics.csv"))

# ============================================================================
# Section 10: Cross-Comparison Statistical Tests
# ============================================================================
cat("\n=== Section 10: Cross-comparison tests ===\n")
co60  <- deg_annotated_list[["radiation_effect"]][tissue_class == "tissue_specific"]
cs137 <- deg_annotated_list[["Cs137_100cGy"]][tissue_class == "tissue_specific"]
rad_array <- array(0, dim=c(2,2,5), dimnames=list(c("Co60","Cs137"),c("no","yes"),BROAD_ORGANS))
for (bo in BROAD_ORGANS) {
  rad_array["Co60","yes",bo] <- sum(co60$predominant_broad_organ==bo & co60$deg_flag=="yes")
  rad_array["Co60","no",bo]  <- sum(co60$predominant_broad_organ==bo & co60$deg_flag=="no")
  rad_array["Cs137","yes",bo]<- sum(cs137$predominant_broad_organ==bo & cs137$deg_flag=="yes")
  rad_array["Cs137","no",bo] <- sum(cs137$predominant_broad_organ==bo & cs137$deg_flag=="no")
}
cmh_rad <- mantelhaen.test(rad_array)

cs10  <- deg_annotated_list[["Cs137_10cGy"]][tissue_class == "tissue_specific"]
cs100 <- deg_annotated_list[["Cs137_100cGy"]][tissue_class == "tissue_specific"]
dose_array <- array(0, dim=c(2,2,5), dimnames=list(c("10cGy","100cGy"),c("no","yes"),BROAD_ORGANS))
for (bo in BROAD_ORGANS) {
  dose_array["10cGy","yes",bo] <- sum(cs10$predominant_broad_organ==bo & cs10$deg_flag=="yes")
  dose_array["10cGy","no",bo]  <- sum(cs10$predominant_broad_organ==bo & cs10$deg_flag=="no")
  dose_array["100cGy","yes",bo]<- sum(cs100$predominant_broad_organ==bo & cs100$deg_flag=="yes")
  dose_array["100cGy","no",bo] <- sum(cs100$predominant_broad_organ==bo & cs100$deg_flag=="no")
}
cmh_dose <- mantelhaen.test(dose_array)

rad_eff <- deg_annotated_list[["radiation_effect"]][tissue_class=="tissue_specific" & deg_flag=="yes"]
gen_int <- deg_annotated_list[["genotype_interaction"]][tissue_class=="tissue_specific" & deg_flag=="yes"]
gen_chi <- chisq.test(rbind(
  table(factor(rad_eff$predominant_broad_organ, levels=BROAD_ORGANS)),
  table(factor(gen_int$predominant_broad_organ, levels=BROAD_ORGANS))))
dir_chi <- chisq.test(rbind(
  table(factor(rad_eff$direction, levels=c("up","down"))),
  table(factor(gen_int$direction, levels=c("up","down")))))
all4_chi <- chisq.test(rbind(
  table(factor(rad_eff$predominant_broad_organ, levels=BROAD_ORGANS)),
  table(factor(gen_int$predominant_broad_organ, levels=BROAD_ORGANS)),
  table(factor(cs100[tissue_class=="tissue_specific"&deg_flag=="yes"]$predominant_broad_organ, levels=BROAD_ORGANS)),
  table(factor(cs10[tissue_class=="tissue_specific"&deg_flag=="yes"]$predominant_broad_organ, levels=BROAD_ORGANS))))

cross_tests <- rbindlist(list(
  data.table(test="CMH (Co-60 vs Cs-137)", statistic=round(unname(cmh_rad$statistic),2),
    df=unname(cmh_rad$parameter), pvalue=cmh_rad$p.value,
    common_OR=round(as.numeric(cmh_rad$estimate),3),
    ci95_low=round(cmh_rad$conf.int[1],3), ci95_high=round(cmh_rad$conf.int[2],3)),
  data.table(test="CMH (Cs-137 10cGy vs 100cGy)", statistic=round(unname(cmh_dose$statistic),2),
    df=unname(cmh_dose$parameter), pvalue=cmh_dose$p.value,
    common_OR=round(as.numeric(cmh_dose$estimate),3),
    ci95_low=round(cmh_dose$conf.int[1],3), ci95_high=round(cmh_dose$conf.int[2],3)),
  data.table(test="Chi-sq (genotype organ distribution)", statistic=round(unname(gen_chi$statistic),2),
    df=unname(gen_chi$parameter), pvalue=gen_chi$p.value,
    common_OR=NA_real_, ci95_low=NA_real_, ci95_high=NA_real_),
  data.table(test="Chi-sq (direction: radiation vs interaction)", statistic=round(unname(dir_chi$statistic),2),
    df=unname(dir_chi$parameter), pvalue=dir_chi$p.value,
    common_OR=NA_real_, ci95_low=NA_real_, ci95_high=NA_real_),
  data.table(test="Chi-sq (overall profile, 4 comparisons)", statistic=round(unname(all4_chi$statistic),2),
    df=unname(all4_chi$parameter), pvalue=all4_chi$p.value,
    common_OR=NA_real_, ci95_low=NA_real_, ci95_high=NA_real_)
), fill=TRUE)
fwrite(cross_tests, file.path(RESULTS_DIR, "cross_comparison_statistical_tests.csv"))

# ============================================================================
# Section 11: Per-Comparison Figures
# ============================================================================
make_fig1_bar <- function(label, outfile) {
  degs <- deg_annotated_list[[label]][deg_flag=="yes" & tissue_class=="tissue_specific"]
  pd <- degs[, .N, by=.(predominant_broad_organ, direction)][direction %in% c("up","down")]
  pd[, predominant_broad_organ := factor(predominant_broad_organ, levels=BROAD_ORGANS)]
  p <- ggplot(pd, aes(x=predominant_broad_organ, y=N, fill=direction)) +
    geom_bar(stat="identity", position="stack") +
    scale_fill_manual(values=c("up"=UP_COLOR,"down"=DOWN_COLOR)) +
    labs(title=paste0("DEGs by Organ: ", label), x="Organ", y="Count", fill="Direction") +
    theme_minimal(base_size=14) +
    theme(text=element_text(family="Liberation Sans"), plot.title=element_text(face="bold"),
          axis.text.x=element_text(angle=30, hjust=1))
  ggsave(outfile, p, width=8, height=6, dpi=300, bg="white")
}

make_fig2_volcano <- function(label, outfile) {
  ts <- deg_annotated_list[[label]][tissue_class=="tissue_specific"]
  plots <- list()
  for (bo in BROAD_ORGANS) {
    od <- ts[predominant_broad_organ==bo]
    if (nrow(od)==0) next
    od[, sig := ifelse(deg_flag=="yes","DEG","NS")]
    p <- ggplot(od, aes(x=log2FoldChange, y=-log10(padj), color=sig)) +
      geom_point(alpha=0.6, size=1.5) +
      scale_color_manual(values=c("DEG"=ORGAN_COLORS[bo],"NS"="grey80")) +
      geom_vline(xintercept=c(-1,1), linetype="dashed", alpha=0.5) +
      geom_hline(yintercept=-log10(0.05), linetype="dashed", alpha=0.5) +
      labs(title=bo, x="log2FC", y="-log10(padj)") +
      theme_minimal(base_size=12) +
      theme(legend.position="none", plot.title=element_text(face="bold", color=ORGAN_COLORS[bo]))
    plots[[bo]] <- p
  }
  if (length(plots) > 0) {
    nc <- min(length(plots), 3)
    comb <- arrangeGrob(grobs=plots, ncol=nc,
      top=textGrob(paste0("Volcano: ", label), gp=gpar(fontsize=18, fontface="bold")))
    ggsave(outfile, comb, width=5*nc, height=4*ceiling(length(plots)/nc), dpi=300, bg="white")
  }
}

make_fig3_heatmap <- function(label, outfile) {
  degs <- deg_annotated_list[[label]][deg_flag=="yes" & tissue_class=="tissue_specific"]
  top50 <- degs[order(padj)][1:min(50, nrow(degs))]
  gi <- intersect(top50$gene_id, rownames(sub_tissue_expr))
  if (length(gi) < 5) return(FALSE)
  ez <- t(scale(t(sub_tissue_expr[gi, ])))
  dv <- top50[match(rownames(ez), gene_id), direction]
  ra <- rowAnnotation(Direction=dv, col=list(Direction=c("up"=UP_COLOR,"down"=DOWN_COLOR)))
  ht <- Heatmap(ez, name="Z", col=colorRamp2(c(-2,0,2),c("#0072B2","white","#D55E00")),
    show_row_names=FALSE, left_annotation=ra,
    column_title=paste0("Top DEGs: ", label), cluster_columns=FALSE)
  png(outfile, width=3600, height=2400, res=300, bg="white")
  draw(ht); dev.off()
}

make_fig4_upset <- function(label, outfile) {
  degs <- deg_annotated_list[[label]][deg_flag=="yes" & tissue_class=="tissue_specific"]
  ol <- lapply(BROAD_ORGANS, function(bo) degs[predominant_broad_organ==bo, gene_id])
  names(ol) <- BROAD_ORGANS
  png(outfile, width=3000, height=1800, res=300, bg="white")
  print(upset(fromList(ol), order.by="freq", sets.bar.color=unname(ORGAN_COLORS[BROAD_ORGANS]), text.scale=1.2))
  dev.off()
}

make_fig5_forest <- function(label, outfile) {
  sd <- all_stats[[label]]$fisher
  if (is.null(sd)) return(FALSE)
  pd <- copy(sd)
  pd[, sig_label := ifelse(padj_bonf<0.001,"***",ifelse(padj_bonf<0.01,"**",ifelse(padj_bonf<0.05,"*","ns")))]
  pd[, broad_organ := factor(broad_organ, levels=rev(BROAD_ORGANS))]
  p <- ggplot(pd, aes(x=organ_rate, y=broad_organ)) +
    geom_vline(xintercept=mean(pd$organ_rate), linetype="dashed", alpha=0.5) +
    geom_point(aes(color=enrichment), size=4) +
    geom_text(aes(label=sig_label), hjust=-0.5, size=5) +
    scale_color_manual(values=c("enriched"="#009E73","depleted"="#D55E00","ns"="grey60")) +
    labs(title=paste0("Tissue Enrichment: ", label), x="DEG Rate (%)", y="Organ", color="Result") +
    theme_minimal(base_size=14) +
    theme(text=element_text(family="Liberation Sans"), plot.title=element_text(face="bold"))
  ggsave(outfile, p, width=9, height=5, dpi=300, bg="white")
}

cat("\n=== Section 11: Per-comparison figures ===\n")
fig_names <- c("deg_counts_by_organ","volcano_plots_by_organ","heatmap_top_degs","upset_organ_overlaps","tissue_enrichment_test")
fig_funcs <- list(make_fig1_bar, make_fig2_volcano, make_fig3_heatmap, make_fig4_upset, make_fig5_forest)
for (label in names(deg_annotated_list)) {
  nd <- sum(deg_annotated_list[[label]]$deg_flag=="yes" & deg_annotated_list[[label]]$tissue_class=="tissue_specific")
  for (fn in 1:5) {
    if (nd < ifelse(fn %in% c(3,4), 5, 3)) next
    outfile <- file.path(RESULTS_DIR, sprintf("fig%d_%s_%s.png", fn, fig_names[fn], label))
    tryCatch({ fig_funcs[[fn]](label, outfile); cat(sprintf("  [%s] fig%d OK\n", label, fn)) },
      error = function(e) cat(sprintf("  [%s] fig%d ERROR: %s\n", label, fn, e$message)))
  }
}

# ============================================================================
# Section 12: Cross-Comparison Figures
# ============================================================================
cat("\n=== Section 12: Cross-comparison figures ===\n")
rate_mat <- as.matrix(dcast(cross_summary, broad_organ ~ comparison, value.var="deg_rate")[, -1])
rownames(rate_mat) <- BROAD_ORGANS; rate_mat[is.na(rate_mat)] <- 0
rad_types <- c("Gamma Co-60","Co-60 int.","Gamma Cs-137","Gamma Cs-137","Sim. GCR","Sim. GCR")
pheatmap(rate_mat, color=colorRampPalette(c("white","#E69F00","#D55E00"))(100),
  cluster_rows=FALSE, cluster_cols=FALSE,
  annotation_col=data.frame(Radiation=rad_types, row.names=colnames(rate_mat)),
  display_numbers=TRUE, number_format="%.1f",
  main="Tissue-Specific DEG Rate (%) by Organ x Comparison",
  filename=file.path(RESULTS_DIR, "fig_cross1_tissue_response_heatmap.png"), width=10, height=5)

comp_order <- c("radiation_effect","genotype_interaction","Cs137_10cGy","Cs137_100cGy","GCR40","GCR80")
dm <- melt(cross_summary[, .(broad_organ,comparison,upregulated,downregulated)],
  id.vars=c("broad_organ","comparison"), variable.name="direction", value.name="count")
dm[, broad_organ := factor(broad_organ, levels=BROAD_ORGANS)]
dm[, comparison := factor(comparison, levels=comp_order)]
p2 <- ggplot(dm, aes(x=broad_organ, y=count, fill=direction)) +
  geom_bar(stat="identity", position="dodge") + facet_wrap(~comparison, ncol=3) +
  scale_fill_manual(values=c("upregulated"=UP_COLOR,"downregulated"=DOWN_COLOR)) +
  labs(title="DEG Direction by Organ Across Comparisons", x="Organ", y="Count", fill="Direction") +
  theme_minimal(base_size=12) + theme(axis.text.x=element_text(angle=30, hjust=1))
ggsave(file.path(RESULTS_DIR, "fig_cross2_direction_comparison.png"), p2, width=12, height=8, dpi=300, bg="white")

dd <- rbind(
  cross_summary[comparison=="Cs137_10cGy", .(broad_organ, dose="10 cGy", source="Cs-137", DEGs)],
  cross_summary[comparison=="Cs137_100cGy", .(broad_organ, dose="100 cGy", source="Cs-137", DEGs)],
  cross_summary[comparison=="GCR40", .(broad_organ, dose="40 cGy", source="GCR", DEGs)],
  cross_summary[comparison=="GCR80", .(broad_organ, dose="80 cGy", source="GCR", DEGs)])
dd[, broad_organ := factor(broad_organ, levels=BROAD_ORGANS)]
dd[, dose := factor(dose, levels=c("10 cGy","40 cGy","80 cGy","100 cGy"))]
p3 <- ggplot(dd, aes(x=dose, y=DEGs, group=broad_organ, color=broad_organ)) +
  geom_line(linewidth=1.2) + geom_point(size=3) + facet_wrap(~source, scales="free_x") +
  scale_color_manual(values=ORGAN_COLORS) +
  labs(title="Dose-Response: Tissue-Specific DEGs by Organ", x="Dose", y="DEGs", color="Organ") +
  theme_minimal(base_size=14)
ggsave(file.path(RESULTS_DIR, "fig_cross3_dose_response.png"), p3, width=10, height=6, dpi=300, bg="white")

rt <- rbind(
  cross_summary[comparison=="radiation_effect", .(broad_organ, condition="Co-60 (100 Gy)", deg_rate)],
  cross_summary[comparison=="Cs137_100cGy", .(broad_organ, condition="Cs-137 (100 cGy)", deg_rate)])
rt[, broad_organ := factor(broad_organ, levels=BROAD_ORGANS)]
p4 <- ggplot(rt, aes(x=broad_organ, y=deg_rate, fill=condition)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=c("Co-60 (100 Gy)"="#0072B2","Cs-137 (100 cGy)"="#009E73")) +
  labs(title="Radiation Type: Co-60 vs Cs-137", x="Organ", y="DEG Rate (%)", fill="Radiation") +
  theme_minimal(base_size=14) + theme(axis.text.x=element_text(angle=30, hjust=1))
ggsave(file.path(RESULTS_DIR, "fig_cross4_radiation_type_comparison.png"), p4, width=10, height=6, dpi=300, bg="white")

# ============================================================================
# Section 13: Old vs New Comparison (OSD498_510)
# ============================================================================
if (file.exists(OLD_FILE)) {
  cat("\n=== Section 13: Old vs New comparison ===\n")
  old_dt <- fread(OLD_FILE)
  colnames(old_dt)[1] <- "gene_id"
  colnames(old_dt)[ncol(old_dt)] <- "deg_flag_old"
  old_dt[, deg_flag_old := trimws(gsub("\r", "", deg_flag_old))]
  new_dt <- deg_annotated_list[["radiation_effect"]]
  comp <- merge(old_dt[, .(gene_id, log2FC_old=log2FoldChange, padj_old=padj, deg_flag_old)],
    new_dt[, .(gene_id, log2FC_new=log2FoldChange, padj_new=padj, deg_flag_new=deg_flag)],
    by="gene_id", all=TRUE)
  comp[, status := fifelse(deg_flag_old=="yes"&deg_flag_new=="yes","both_DEG",
    fifelse(deg_flag_old=="yes"&deg_flag_new=="no","lost_DEG",
    fifelse(deg_flag_old=="no"&deg_flag_new=="yes","gained_DEG","both_non_DEG")))]
  fwrite(comp, file.path(RESULTS_DIR, "old_vs_new_comparison_OSD498_510.csv"))
  cat(sprintf("  Old: %d, New: %d, Shared: %d, Lost: %d, Gained: %d\n",
    comp[deg_flag_old=="yes",.N], comp[deg_flag_new=="yes",.N],
    comp[status=="both_DEG",.N], comp[status=="lost_DEG",.N], comp[status=="gained_DEG",.N]))
}

# ============================================================================
# Section 14: Session Info
# ============================================================================
cat("\n=== Section 14: Session Info ===\n")
print(sessionInfo())
cat("\n=== Analysis complete ===\n")
