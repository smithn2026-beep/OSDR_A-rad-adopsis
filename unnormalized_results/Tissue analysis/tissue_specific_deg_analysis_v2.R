# ============================================================================
# Tissue-Specific DEG Analysis v2 — Corrected (Unnormalized) DESeq2 Data
# Arabidopsis Ionizing Radiation: 6 Comparisons from NASA GeneLab OSD Studies
#
# Input: 6 DESeq2 output CSVs (re-run with raw unnormalized counts)
# Method: AtGenExpress developmental atlas + Tau tissue-specificity index
# Threshold: padj < 0.05 AND |log2FC| >= 1
#
# Comparisons:
#   1. radiation_effect       — OSD498+510, Co-60 100Gy vs none (WT)
#   2. genotype_interaction   — OSD508+510, sog1-1 x radiation interaction
#   3. GCR40                  — OSD658, simulated GCR 40cGy vs none (WT)
#   4. GCR80                  — OSD658, simulated GCR 80cGy vs none (WT)
#   5. Cs137_100cGy           — OSD782, Cs-137 100cGy vs none (WT)
#   6. Cs137_10cGy            — OSD782, Cs-137 10cGy vs none (WT)
# ============================================================================

# ==== Section 1: Setup & Configuration =====================================

library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)
library(grid)
library(gridExtra)
library(pheatmap)

# Configuration
TAU_THRESHOLD <- 0.6
PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD <- 1.0  # |log2FC| >= 1 (2-fold change)

# Color palette
UP_COLOR <- "#0072B2"
DOWN_COLOR <- "#D55E00"
ORGAN_COLORS <- c(
  Root = "#0072B2", Seedling = "#009E73", Leaf_Shoot = "#E69F00",
  Flower = "#CC79A7", Seed_Silique = "#D55E00"
)
BROAD_ORGANS <- c("Root", "Seedling", "Leaf_Shoot", "Flower", "Seed_Silique")

# Define the 6 comparison files
COMPARISON_FILES <- list(
  list(label = "radiation_effect",     file = "DEG_unnorm_OSD498_510_radiation_effect.csv",
       study = "OSD498+510", radiation = "Gamma Co-60 100Gy", design = "~ Study + Genotype + Radiation"),
  list(label = "genotype_interaction", file = "DEG_unnorm_OSD508_510_genotype_x_radiation_interaction.csv",
       study = "OSD508+510", radiation = "Gamma Co-60 100Gy (interaction)", design = "~ Study + Genotype + Radiation + Genotype:Radiation"),
  list(label = "GCR40",                file = "DEG_unnorm_OSD658_GCR40_vs_none.csv",
       study = "OSD658", radiation = "Simulated GCR 40cGy", design = "~ Radiation"),
  list(label = "GCR80",                file = "DEG_unnorm_OSD658_GCR80_vs_none.csv",
       study = "OSD658", radiation = "Simulated GCR 80cGy", design = "~ Radiation"),
  list(label = "Cs137_100cGy",         file = "DEG_unnorm_OSD782_100cGy_vs_none.csv",
       study = "OSD782", radiation = "Gamma Cs-137 100cGy", design = "~ Timepoint + Radiation"),
  list(label = "Cs137_10cGy",          file = "DEG_unnorm_OSD782_10cGy_vs_none.csv",
       study = "OSD782", radiation = "Gamma Cs-137 10cGy", design = "~ Timepoint + Radiation")
)

RESULTS_DIR <- "/mnt/results"
ATLAS_PATH <- "/workspace/atgenexpress_atlas.rds"

# ==== Section 2: Load & Validate DESeq2 CSVs ===============================

load_deg_data <- function(file_path) {
  dt <- fread(file_path)
  colnames(dt)[1] <- "gene_id"
  
  # Apply DEG threshold
  dt[, deg_flag := ifelse(!is.na(padj) & padj < PADJ_THRESHOLD & abs(log2FoldChange) >= LFC_THRESHOLD, "yes", "no")]
  dt[, direction := ifelse(deg_flag == "yes" & log2FoldChange > 0, "up",
                    ifelse(deg_flag == "yes" & log2FoldChange < 0, "down", "none"))]
  
  # Validate AGI IDs
  dt[, is_nuclear := grepl("^AT[1-5]G[0-9]{5}$", gene_id)]
  dt[, is_organellar := grepl("^AT[CM]G[0-9]{5}$", gene_id)]
  
  return(dt)
}

deg_list <- lapply(COMPARISON_FILES, function(cf) load_deg_data(cf$file))
names(deg_list) <- sapply(COMPARISON_FILES, `[[`, "label")

# Print overview
for (label in names(deg_list)) {
  dt <- deg_list[[label]]
  cat(sprintf("[%s] %d genes, %d DEGs (%d up, %d down)\n",
              label, nrow(dt), sum(dt$deg_flag == "yes"),
              sum(dt$direction == "up"), sum(dt$direction == "down")))
}

# ==== Section 3: Load AtGenExpress Atlas & Compute Tau =====================

compute_tau <- function(expr_matrix) {
  tau <- apply(expr_matrix, 1, function(x) {
    x <- as.numeric(x)
    if (max(x, na.rm = TRUE) == 0) return(NA)
    sum(1 - x / max(x, na.rm = TRUE), na.rm = TRUE) / (sum(!is.na(x)) - 1)
  })
  return(tau)
}

atlas <- readRDS(ATLAS_PATH)
tau_subtissue <- compute_tau(atlas$sub_tissue_expr)
tau_broad <- compute_tau(atlas$broad_organ_expr)

# Build tissue annotation table
gene_ids <- rownames(atlas$sub_tissue_expr)
tissue_ann <- data.table(gene_id = gene_ids,
  tau_subtissue = tau_subtissue[gene_ids],
  tau_broad = tau_broad[gene_ids],
  is_tissue_specific = tau_subtissue[gene_ids] >= TAU_THRESHOLD)

for (g in gene_ids) {
  if (tissue_ann[gene_id == g, is_tissue_specific]) {
    tissue_ann[gene_id == g, predominant_subtissue := colnames(atlas$sub_tissue_expr)[which.max(atlas$sub_tissue_expr[g, ])]]
    tissue_ann[gene_id == g, predominant_broad_organ := colnames(atlas$broad_organ_expr)[which.max(atlas$broad_organ_expr[g, ])]]
  }
}

# Sub-tissue to broad-organ mapping
st_bo_map <- unique(atlas$all_samples[, c("sub_tissue", "broad_organ")])

# ==== Section 4: Merge DEGs with Tissue Annotations ========================

deg_annotated_list <- lapply(names(deg_list), function(label) {
  dt <- merge(deg_list[[label]], tissue_ann, by = "gene_id", all.x = TRUE)
  dt[, tissue_class := fifelse(
    is.na(tau_subtissue), "unannotated",
    fifelse(is_tissue_specific == TRUE, "tissue_specific", "constitutive"))]
  dt[is.na(tau_subtissue), tissue_class := fifelse(
    grepl("^AT[CM]G[0-9]{5}$", gene_id), "unannotated_organellar", "unannotated")]
  return(dt)
})
names(deg_annotated_list) <- names(deg_list)

# ==== Section 5: Export Hierarchical CSVs ==================================

export_cols <- c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat",
                 "pvalue", "padj", "deg_flag", "direction",
                 "tau_subtissue", "tau_broad", "predominant_subtissue",
                 "predominant_broad_organ", "tissue_class")

for (label in names(deg_annotated_list)) {
  dt <- deg_annotated_list[[label]]
  comp_dir <- file.path(RESULTS_DIR, "tissue_specific_degs_v2", label)
  degs_only <- dt[deg_flag == "yes"]
  
  fwrite(dt[, ..export_cols], file.path(comp_dir, "all_degs_with_tissue_annotation.csv"))
  fwrite(degs_only[tissue_class == "constitutive", ..export_cols], file.path(comp_dir, "constitutive_DEGs.csv"))
  fwrite(degs_only[tissue_class %like% "unannotated", ..export_cols], file.path(comp_dir, "unannotated_DEGs.csv"))
  
  for (bo in BROAD_ORGANS) {
    organ_dir <- file.path(comp_dir, bo)
    dir.create(organ_dir, showWarnings = FALSE, recursive = TRUE)
    organ_degs <- degs_only[tissue_class == "tissue_specific" & predominant_broad_organ == bo, ..export_cols]
    fwrite(organ_degs, file.path(organ_dir, "all_DEGs.csv"))
    organ_subtissues <- as.character(st_bo_map[st_bo_map$broad_organ == bo, "sub_tissue"])
    for (st in organ_subtissues) {
      fwrite(organ_degs[predominant_subtissue == st, ..export_cols], file.path(organ_dir, paste0(st, ".csv")))
    }
  }
}

# ==== Section 6: Summary Tables ============================================

for (label in names(deg_annotated_list)) {
  degs <- deg_annotated_list[[label]][deg_flag == "yes" & tissue_class == "tissue_specific"]
  summary_rows <- list()
  for (bo in BROAD_ORGANS) {
    organ_subtissues <- as.character(st_bo_map[st_bo_map$broad_organ == bo, "sub_tissue"])
    for (st in organ_subtissues) {
      st_degs <- degs[predominant_subtissue == st]
      if (nrow(st_degs) > 0) {
        summary_rows[[length(summary_rows)+1]] <- data.table(
          broad_organ = bo, sub_tissue = st, total_DEGs = nrow(st_degs),
          upregulated = sum(st_degs$direction == "up"), downregulated = sum(st_degs$direction == "down"),
          median_log2FC = round(median(st_degs$log2FoldChange, na.rm = TRUE), 3),
          mean_log2FC = round(mean(st_degs$log2FoldChange, na.rm = TRUE), 3),
          mean_tau_subtissue = round(mean(st_degs$tau_subtissue, na.rm = TRUE), 4))
      }
    }
  }
  fwrite(rbindlist(summary_rows), file.path(RESULTS_DIR, paste0("tissue_deg_summary_", label, ".csv")))
}

# Cross-comparison summary
cross_rows <- list()
for (label in names(deg_annotated_list)) {
  degs <- deg_annotated_list[[label]][deg_flag == "yes" & tissue_class == "tissue_specific"]
  for (bo in BROAD_ORGANS) {
    total_in_organ <- tissue_ann[is_tissue_specific == TRUE & predominant_broad_organ == bo, .N]
    cross_rows[[length(cross_rows)+1]] <- data.table(
      broad_organ = bo, comparison = label,
      DEGs = degs[predominant_broad_organ == bo, .N],
      upregulated = degs[predominant_broad_organ == bo & direction == "up", .N],
      downregulated = degs[predominant_broad_organ == bo & direction == "down", .N],
      total_atlas_genes = total_in_organ,
      deg_rate = round(degs[predominant_broad_organ == bo, .N] / total_in_organ * 100, 2))
  }
}
cross_summary <- rbindlist(cross_rows)
fwrite(cross_summary, file.path(RESULTS_DIR, "cross_comparison_tissue_summary.csv"))

# ==== Section 7: Per-Comparison Statistical Tests ==========================

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
  
  total_degs <- sum(contig[, "yes"])
  total_non <- sum(contig[, "no"])
  
  # Test 1: Chi-square
  chi_test <- if (total_degs >= 5) chisq.test(contig) else NULL
  
  # Test 2: Per-tissue Fisher's exact
  fisher_rows <- list()
  for (bo in BROAD_ORGANS) {
    fisher_mat <- matrix(c(contig[bo, "yes"], contig[bo, "no"],
                           total_degs - contig[bo, "yes"], total_non - contig[bo, "no"]), nrow = 2)
    ft <- fisher.test(fisher_mat)
    fisher_rows[[bo]] <- data.table(comparison = label, broad_organ = bo,
      organ_DEGs = contig[bo, "yes"], organ_total = sum(contig[bo, ]),
      organ_rate = round(contig[bo, "yes"] / sum(contig[bo, ]) * 100, 2),
      odds_ratio = round(unname(ft$estimate), 3),
      ci95_low = round(ft$conf.int[1], 3), ci95_high = round(ft$conf.int[2], 3),
      pvalue = ft$p.value,
      padj_bonf = p.adjust(ft$p.value, method = "bonferroni", n = 5),
      enrichment = ifelse(unname(ft$estimate) > 1 & ft$p.value < 0.05, "enriched",
                   ifelse(unname(ft$estimate) < 1 & ft$p.value < 0.05, "depleted", "ns")))
  }
  fisher_dt <- rbindlist(fisher_rows)
  
  # Test 3: Pairwise Fisher's
  pairwise_dt <- NULL
  if (total_degs >= 5) {
    pairwise_rows <- list()
    for (pr in combn(BROAD_ORGANS, 2, simplify = FALSE)) {
      mat <- matrix(c(contig[pr[1], "yes"], contig[pr[1], "no"],
                      contig[pr[2], "yes"], contig[pr[2], "no"]), nrow = 2)
      ft <- fisher.test(mat)
      pairwise_rows[[paste(pr[1], pr[2], sep = "_vs_")]] <- data.table(
        comparison = label, organ_1 = pr[1], organ_2 = pr[2],
        rate_1 = round(contig[pr[1], "yes"] / sum(contig[pr[1], ]) * 100, 2),
        rate_2 = round(contig[pr[2], "yes"] / sum(contig[pr[2], ]) * 100, 2),
        odds_ratio = round(unname(ft$estimate), 3), pvalue = ft$p.value,
        padj_bonf = p.adjust(ft$p.value, method = "bonferroni", n = 10))
    }
    pairwise_dt <- rbindlist(pairwise_rows)
  }
  
  # Test 4: Direction chi-square
  dir_result <- NULL
  degs_only <- ts_genes[deg_flag == "yes"]
  if (nrow(degs_only) >= 5 && sum(degs_only$direction == "up") > 0 && sum(degs_only$direction == "down") > 0) {
    dir_contig <- table(degs_only$predominant_broad_organ, degs_only$direction)
    dir_contig <- dir_contig[rowSums(dir_contig) > 0, , drop = FALSE]
    if (nrow(dir_contig) >= 2) {
      dir_chi <- chisq.test(dir_contig)
      dir_result <- list(chi = dir_chi, pcts = data.table(comparison = label,
        broad_organ = rownames(dir_contig), n_DEGs = rowSums(dir_contig),
        pct_up = round(dir_contig[, "up"] / rowSums(dir_contig) * 100, 1),
        pct_down = round(dir_contig[, "down"] / rowSums(dir_contig) * 100, 1)))
    }
  }
  
  return(list(chi = chi_test, fisher = fisher_dt, pairwise = pairwise_dt, direction = dir_result))
}

all_stats <- lapply(names(deg_annotated_list), function(label) {
  run_stat_tests(deg_annotated_list[[label]], label)
})
names(all_stats) <- names(deg_annotated_list)

# Save per-comparison stats
for (label in names(all_stats)) {
  fwrite(all_stats[[label]]$fisher, file.path(RESULTS_DIR, paste0("tissue_enrichment_stats_", label, ".csv")))
  if (!is.null(all_stats[[label]]$pairwise)) {
    fwrite(all_stats[[label]]$pairwise, file.path(RESULTS_DIR, paste0("pairwise_comparisons_", label, ".csv")))
  }
}
enrichment_all <- rbindlist(lapply(all_stats, function(x) x$fisher))
fwrite(enrichment_all, file.path(RESULTS_DIR, "cross_comparison_statistics.csv"))

# ==== Section 8: Cross-Comparison Statistical Tests ========================

# CMH test: radiation type (Co-60 vs Cs-137)
co60 <- deg_annotated_list[["radiation_effect"]][tissue_class == "tissue_specific"]
cs137 <- deg_annotated_list[["Cs137_100cGy"]][tissue_class == "tissue_specific"]
rad_array <- array(0, dim = c(2, 2, 5), dimnames = list(c("Co60","Cs137"), c("no","yes"), BROAD_ORGANS))
for (bo in BROAD_ORGANS) {
  rad_array["Co60","yes",bo] <- sum(co60$predominant_broad_organ == bo & co60$deg_flag == "yes")
  rad_array["Co60","no",bo]  <- sum(co60$predominant_broad_organ == bo & co60$deg_flag == "no")
  rad_array["Cs137","yes",bo] <- sum(cs137$predominant_broad_organ == bo & cs137$deg_flag == "yes")
  rad_array["Cs137","no",bo]  <- sum(cs137$predominant_broad_organ == bo & cs137$deg_flag == "no")
}
cmh_rad <- mantelhaen.test(rad_array)

# CMH test: Cs-137 dose-response
cs10 <- deg_annotated_list[["Cs137_10cGy"]][tissue_class == "tissue_specific"]
cs100 <- deg_annotated_list[["Cs137_100cGy"]][tissue_class == "tissue_specific"]
dose_array <- array(0, dim = c(2, 2, 5), dimnames = list(c("10cGy","100cGy"), c("no","yes"), BROAD_ORGANS))
for (bo in BROAD_ORGANS) {
  dose_array["10cGy","yes",bo] <- sum(cs10$predominant_broad_organ == bo & cs10$deg_flag == "yes")
  dose_array["10cGy","no",bo]  <- sum(cs10$predominant_broad_organ == bo & cs10$deg_flag == "no")
  dose_array["100cGy","yes",bo] <- sum(cs100$predominant_broad_organ == bo & cs100$deg_flag == "yes")
  dose_array["100cGy","no",bo]  <- sum(cs100$predominant_broad_organ == bo & cs100$deg_flag == "no")
}
cmh_dose <- mantelhaen.test(dose_array)

# Chi-square: genotype effect (organ distribution)
rad_eff <- deg_annotated_list[["radiation_effect"]][tissue_class == "tissue_specific" & deg_flag == "yes"]
gen_int <- deg_annotated_list[["genotype_interaction"]][tissue_class == "tissue_specific" & deg_flag == "yes"]
gen_chi <- chisq.test(rbind(
  table(factor(rad_eff$predominant_broad_organ, levels = BROAD_ORGANS)),
  table(factor(gen_int$predominant_broad_organ, levels = BROAD_ORGANS))))

# Chi-square: direction comparison
dir_chi <- chisq.test(rbind(
  table(factor(rad_eff$direction, levels = c("up","down"))),
  table(factor(gen_int$direction, levels = c("up","down")))))

# Chi-square: overall profile
all4_chi <- chisq.test(rbind(
  table(factor(rad_eff$predominant_broad_organ, levels = BROAD_ORGANS)),
  table(factor(gen_int$predominant_broad_organ, levels = BROAD_ORGANS)),
  table(factor(cs100[tissue_class == "tissue_specific" & deg_flag == "yes"]$predominant_broad_organ, levels = BROAD_ORGANS)),
  table(factor(cs10[tissue_class == "tissue_specific" & deg_flag == "yes"]$predominant_broad_organ, levels = BROAD_ORGANS))))

cross_tests <- rbindlist(list(
  data.table(test = "CMH (Co-60 vs Cs-137)", statistic = round(unname(cmh_rad$statistic),2),
    df = unname(cmh_rad$parameter), pvalue = cmh_rad$p.value,
    common_OR = round(as.numeric(cmh_rad$estimate),3),
    ci95_low = round(cmh_rad$conf.int[1],3), ci95_high = round(cmh_rad$conf.int[2],3)),
  data.table(test = "CMH (Cs-137 10cGy vs 100cGy)", statistic = round(unname(cmh_dose$statistic),2),
    df = unname(cmh_dose$parameter), pvalue = cmh_dose$p.value,
    common_OR = round(as.numeric(cmh_dose$estimate),3),
    ci95_low = round(cmh_dose$conf.int[1],3), ci95_high = round(cmh_dose$conf.int[2],3)),
  data.table(test = "Chi-sq (genotype organ distribution)", statistic = round(unname(gen_chi$statistic),2),
    df = unname(gen_chi$parameter), pvalue = gen_chi$p.value,
    common_OR = NA_real_, ci95_low = NA_real_, ci95_high = NA_real_),
  data.table(test = "Chi-sq (direction: radiation vs interaction)", statistic = round(unname(dir_chi$statistic),2),
    df = unname(dir_chi$parameter), pvalue = dir_chi$p.value,
    common_OR = NA_real_, ci95_low = NA_real_, ci95_high = NA_real_),
  data.table(test = "Chi-sq (overall profile, 4 comparisons)", statistic = round(unname(all4_chi$statistic),2),
    df = unname(all4_chi$parameter), pvalue = all4_chi$p.value,
    common_OR = NA_real_, ci95_low = NA_real_, ci95_high = NA_real_)
), fill = TRUE)
fwrite(cross_tests, file.path(RESULTS_DIR, "cross_comparison_statistical_tests.csv"))

# ==== Section 9: Per-Comparison Figures ====================================

make_fig1_bar <- function(label, outfile) {
  degs <- deg_annotated_list[[label]][deg_flag == "yes" & tissue_class == "tissue_specific"]
  plot_data <- degs[, .N, by = .(predominant_broad_organ, direction)][direction %in% c("up","down")]
  plot_data[, predominant_broad_organ := factor(predominant_broad_organ, levels = BROAD_ORGANS)]
  p <- ggplot(plot_data, aes(x = predominant_broad_organ, y = N, fill = direction)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("up" = UP_COLOR, "down" = DOWN_COLOR)) +
    labs(title = paste0("DEGs by Organ: ", label), x = "Organ", y = "Count", fill = "Direction") +
    theme_minimal(base_size = 14) +
    theme(text = element_text(family = "Liberation Sans"), plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(outfile, p, width = 8, height = 6, dpi = 300, bg = "white")
}

make_fig2_volcano <- function(label, outfile) {
  ts <- deg_annotated_list[[label]][tissue_class == "tissue_specific"]
  plots <- list()
  for (bo in BROAD_ORGANS) {
    od <- ts[predominant_broad_organ == bo]
    if (nrow(od) == 0) next
    od[, sig := ifelse(deg_flag == "yes", "DEG", "NS")]
    p <- ggplot(od, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("DEG" = ORGAN_COLORS[bo], "NS" = "grey80")) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
      labs(title = bo, x = "log2FC", y = "-log10(padj)") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none", plot.title = element_text(face = "bold", color = ORGAN_COLORS[bo]))
    plots[[bo]] <- p
  }
  if (length(plots) > 0) {
    ncol_g <- min(length(plots), 3)
    combined <- arrangeGrob(grobs = plots, ncol = ncol_g,
                            top = textGrob(paste0("Volcano: ", label), gp = gpar(fontsize = 18, fontface = "bold")))
    ggsave(outfile, combined, width = 5 * ncol_g, height = 4 * ceiling(length(plots)/ncol_g), dpi = 300, bg = "white")
  }
}

make_fig3_heatmap <- function(label, outfile) {
  degs <- deg_annotated_list[[label]][deg_flag == "yes" & tissue_class == "tissue_specific"]
  top50 <- degs[order(padj)][1:min(50, nrow(degs))]
  genes_in_atlas <- intersect(top50$gene_id, rownames(atlas$sub_tissue_expr))
  if (length(genes_in_atlas) < 5) return(FALSE)
  expr_z <- t(scale(t(atlas$sub_tissue_expr[genes_in_atlas, ])))
  dir_vec <- top50[match(rownames(expr_z), gene_id), direction]
  ra <- rowAnnotation(Direction = dir_vec, col = list(Direction = c("up" = UP_COLOR, "down" = DOWN_COLOR)))
  ht <- Heatmap(expr_z, name = "Z", col = colorRamp2(c(-2,0,2), c("#0072B2","white","#D55E00")),
                show_row_names = FALSE, left_annotation = ra,
                column_title = paste0("Top DEGs: ", label), cluster_columns = FALSE)
  png(outfile, width = 3600, height = 2400, res = 300, bg = "white")
  draw(ht)
  dev.off()
}

make_fig4_upset <- function(label, outfile) {
  degs <- deg_annotated_list[[label]][deg_flag == "yes" & tissue_class == "tissue_specific"]
  organ_lists <- lapply(BROAD_ORGANS, function(bo) degs[predominant_broad_organ == bo, gene_id])
  names(organ_lists) <- BROAD_ORGANS
  png(outfile, width = 3000, height = 1800, res = 300, bg = "white")
  print(upset(fromList(organ_lists), order.by = "freq",
              sets.bar.color = unname(ORGAN_COLORS[BROAD_ORGANS]), text.scale = 1.2))
  dev.off()
}

make_fig5_forest <- function(label, outfile) {
  stats_dt <- all_stats[[label]]$fisher
  if (is.null(stats_dt)) return(FALSE)
  plot_data <- copy(stats_dt)
  plot_data[, sig_label := ifelse(padj_bonf < 0.001, "***", ifelse(padj_bonf < 0.01, "**", ifelse(padj_bonf < 0.05, "*", "ns")))]
  plot_data[, broad_organ := factor(broad_organ, levels = rev(BROAD_ORGANS))]
  p <- ggplot(plot_data, aes(x = organ_rate, y = broad_organ)) +
    geom_vline(xintercept = mean(plot_data$organ_rate), linetype = "dashed", alpha = 0.5) +
    geom_point(aes(color = enrichment), size = 4) +
    geom_text(aes(label = sig_label), hjust = -0.5, size = 5) +
    scale_color_manual(values = c("enriched" = "#009E73", "depleted" = "#D55E00", "ns" = "grey60")) +
    labs(title = paste0("Tissue Enrichment: ", label), x = "DEG Rate (%)", y = "Organ", color = "Result") +
    theme_minimal(base_size = 14) +
    theme(text = element_text(family = "Liberation Sans"), plot.title = element_text(face = "bold"))
  ggsave(outfile, p, width = 9, height = 5, dpi = 300, bg = "white")
}

fig_names <- c("deg_counts_by_organ", "volcano_plots_by_organ", "heatmap_top_degs",
               "upset_organ_overlaps", "tissue_enrichment_test")
fig_funcs <- list(make_fig1_bar, make_fig2_volcano, make_fig3_heatmap, make_fig4_upset, make_fig5_forest)

for (label in names(deg_annotated_list)) {
  n_degs <- sum(deg_annotated_list[[label]]$deg_flag == "yes" & deg_annotated_list[[label]]$tissue_class == "tissue_specific")
  for (fig_num in 1:5) {
    min_degs <- if (fig_num %in% c(3, 4)) 5 else 3
    if (n_degs < min_degs) next
    outfile <- file.path(RESULTS_DIR, sprintf("fig%d_%s_%s.png", fig_num, fig_names[fig_num], label))
    tryCatch({ fig_funcs[[fig_num]](label, outfile) }, error = function(e) cat(sprintf("  [ERROR] %s fig%d: %s\n", label, fig_num, e$message)))
  }
}

# ==== Section 10: Cross-Comparison Figures =================================

# Fig cross1: Tissue response heatmap
rate_mat <- as.matrix(dcast(cross_summary, broad_organ ~ comparison, value.var = "deg_rate")[, -1])
rownames(rate_mat) <- BROAD_ORGANS
rate_mat[is.na(rate_mat)] <- 0
rad_types <- c("Gamma Co-60", "Co-60 interaction", "Gamma Cs-137", "Gamma Cs-137", "Sim. GCR", "Sim. GCR")
pheatmap(rate_mat, color = colorRampPalette(c("white","#E69F00","#D55E00"))(100),
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col = data.frame(Radiation = rad_types, row.names = colnames(rate_mat)),
         display_numbers = TRUE, number_format = "%.1f",
         main = "Tissue-Specific DEG Rate (%) by Organ x Comparison",
         filename = file.path(RESULTS_DIR, "fig_cross1_tissue_response_heatmap.png"),
         width = 10, height = 5)

# Fig cross2: Direction comparison
dir_melt <- melt(cross_summary[, .(broad_organ, comparison, upregulated, downregulated)],
                 id.vars = c("broad_organ","comparison"), variable.name = "direction", value.name = "count")
p2 <- ggplot(dir_melt, aes(x = broad_organ, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") + facet_wrap(~comparison, ncol = 3) +
  scale_fill_manual(values = c("upregulated" = UP_COLOR, "downregulated" = DOWN_COLOR)) +
  labs(title = "DEG Direction by Organ Across Comparisons", x = "Organ", y = "Count", fill = "Direction") +
  theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(RESULTS_DIR, "fig_cross2_direction_comparison.png"), p2, width = 12, height = 8, dpi = 300, bg = "white")

# Fig cross3: Dose-response
dose_data <- rbind(
  cross_summary[comparison == "Cs137_10cGy", .(broad_organ, dose = "10 cGy", source = "Cs-137", DEGs)],
  cross_summary[comparison == "Cs137_100cGy", .(broad_organ, dose = "100 cGy", source = "Cs-137", DEGs)],
  cross_summary[comparison == "GCR40", .(broad_organ, dose = "40 cGy", source = "GCR", DEGs)],
  cross_summary[comparison == "GCR80", .(broad_organ, dose = "80 cGy", source = "GCR", DEGs)])
p3 <- ggplot(dose_data, aes(x = dose, y = DEGs, group = broad_organ, color = broad_organ)) +
  geom_line(linewidth = 1.2) + geom_point(size = 3) + facet_wrap(~source, scales = "free_x") +
  scale_color_manual(values = ORGAN_COLORS) +
  labs(title = "Dose-Response: Tissue-Specific DEGs by Organ", x = "Dose", y = "DEGs", color = "Organ") +
  theme_minimal(base_size = 14)
ggsave(file.path(RESULTS_DIR, "fig_cross3_dose_response.png"), p3, width = 10, height = 6, dpi = 300, bg = "white")

# Fig cross4: Radiation type comparison
rad_type_data <- rbind(
  cross_summary[comparison == "radiation_effect", .(broad_organ, condition = "Co-60 (100 Gy)", deg_rate)],
  cross_summary[comparison == "Cs137_100cGy", .(broad_organ, condition = "Cs-137 (100 cGy)", deg_rate)])
p4 <- ggplot(rad_type_data, aes(x = broad_organ, y = deg_rate, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Co-60 (100 Gy)" = "#0072B2", "Cs-137 (100 cGy)" = "#009E73")) +
  labs(title = "Radiation Type: Co-60 vs Cs-137", x = "Organ", y = "DEG Rate (%)", fill = "Radiation") +
  theme_minimal(base_size = 14) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(RESULTS_DIR, "fig_cross4_radiation_type_comparison.png"), p4, width = 10, height = 6, dpi = 300, bg = "white")

# ==== Section 11: Old vs New Comparison (OSD498_510) =======================

old_dt <- fread("DEG_OSD498_510_radiation_effect.csv")
colnames(old_dt)[1] <- "gene_id"
colnames(old_dt)[ncol(old_dt)] <- "deg_flag_old"
old_dt[, deg_flag_old := trimws(gsub("\r", "", deg_flag_old))]

new_dt <- deg_annotated_list[["radiation_effect"]]
comparison_dt <- merge(
  old_dt[, .(gene_id, log2FC_old = log2FoldChange, padj_old = padj, deg_flag_old)],
  new_dt[, .(gene_id, log2FC_new = log2FoldChange, padj_new = padj, deg_flag_new = deg_flag)],
  by = "gene_id", all = TRUE)

comparison_dt[, status := fifelse(deg_flag_old == "yes" & deg_flag_new == "yes", "both_DEG",
  fifelse(deg_flag_old == "yes" & deg_flag_new == "no", "lost_DEG",
  fifelse(deg_flag_old == "no" & deg_flag_new == "yes", "gained_DEG", "both_non_DEG")))]

fwrite(comparison_dt, file.path(RESULTS_DIR, "old_vs_new_comparison_OSD498_510.csv"))

# ==== Section 12: Session Info =============================================

sessionInfo()
