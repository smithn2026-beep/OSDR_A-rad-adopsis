# =============================================================
# Sub-Tissue GO + KEGG + KEGG Module Enrichment Analysis
# Run this AFTER tissue_specific_deg_analysis_v2.R
#
# This script performs pathway enrichment (GO Biological Process,
# KEGG pathways, and KEGG modules) for each significant sub-tissue's
# DEG list. It extends pathway_enrichment_498v510.R by adding KEGG
# and KEGG module enrichment alongside the original GO analysis.
#
# Requirements:
#   - The annotated DEG list from tissue_specific_deg_analysis_v2.R
#     (saved as deg_annotated_list_v2.rds)
#   - The sub-tissue statistics (saved as subtissue_stats.rds)
#   - Internet access for KEGG API (enrichKEGG downloads pathway lists)
# =============================================================

# ---- 1. Install/load packages ---------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c("clusterProfiler", "org.At.tair.db", "enrichplot", "ggplot2",
          "AnnotationDbi", "data.table")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, update = FALSE, ask = FALSE)
  }
}

library(clusterProfiler)
library(org.At.tair.db)
library(enrichplot)
library(ggplot2)
library(AnnotationDbi)
library(data.table)

# ---- 2. Configuration -----------------------------------------------------
# Path to the annotated DEG list (from tissue_specific_deg_analysis_v2.R)
DEG_LIST_PATH <- "deg_annotated_list_v2.rds"

# Path to the sub-tissue statistics (from tissue_specific_deg_analysis_v2.R)
SUBTISSUE_STATS_PATH <- "subtissue_stats.rds"

# Output directory for enrichment results
OUTPUT_DIR <- "enrichment_results"

# Significance thresholds
PADJ_CUTOFF <- 0.05
QVALUE_CUTOFF <- 0.2
MIN_GENES <- 5  # minimum gene list size for enrichment

# KEGG organism code for Arabidopsis thaliana
KEGG_ORG <- "ath"

# ---- 3. The reusable enrichment function ----------------------------------
# Arguments:
#   gene_list : vector of TAIR gene IDs (e.g., AT1G01010)
#   universe  : vector of all TAIR gene IDs tested by DESeq2
#   label     : used for output filenames (e.g., "radiation_effect_senescing_leaf_up")
#
# Runs three enrichment tests:
#   1. GO Biological Process (enrichGO)
#   2. KEGG pathway (enrichKEGG)
#   3. KEGG module (enrichMKEGG)
#
# Returns a list with all three result objects.
run_subtissue_enrichment <- function(gene_list, universe, label) {

  cat("\n=========================================\n")
  cat("Running enrichment for:", label, "\n")
  cat("Gene list size:", length(gene_list), "\n")
  cat("=========================================\n")

  if (length(gene_list) < MIN_GENES) {
    cat("Too few genes (<", MIN_GENES, ") for meaningful enrichment -- skipping.\n")
    return(invisible(list(label = label, ego = NULL, ekg = NULL, emk = NULL)))
  }

  results <- list(label = label)

  # --- 1. GO enrichment (Biological Process) ---
  # enrichGO: clusterProfiler function that performs hypergeometric test
  # for over-representation of GO terms in the gene list vs universe.
  # org.At.tair.db: Arabidopsis annotation database mapping TAIR IDs to GO terms.
  # ont = "BP": Biological Process ontology (one of three GO domains).
  # pAdjustMethod = "BH": Benjamini-Hochberg false discovery rate correction.
  cat("\n--- GO BP enrichment ---\n")
  ego <- tryCatch({
    enrichGO(gene = gene_list, universe = universe, OrgDb = org.At.tair.db,
             keyType = "TAIR", ont = "BP", pAdjustMethod = "BH",
             pvalueCutoff = PADJ_CUTOFF, qvalueCutoff = QVALUE_CUTOFF)
  }, error = function(e) {
    cat("  GO error:", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    out_csv <- file.path(OUTPUT_DIR, paste0("GO_enrichment_", label, ".csv"))
    write.csv(as.data.frame(ego), out_csv, row.names = FALSE)
    cat("  Saved:", basename(out_csv), "(", nrow(as.data.frame(ego)), "terms )\n")

    # dotplot: enrichplot function creating a dot plot of enriched terms
    # ggtitle: ggplot2 function adding a title
    # ggsave: ggplot2 function saving the plot to PNG
    p <- dotplot(ego, showCategory = 15) + ggtitle(paste0(label, " -- GO BP"))
    ggsave(file.path(OUTPUT_DIR, paste0("GO_dotplot_", label, ".png")), p, width = 9, height = 7)
    cat("  Saved dotplot\n")
  } else {
    cat("  No enriched GO terms found.\n")
  }
  results$ego <- ego

  # --- 2. KEGG pathway enrichment ---
  # enrichKEGG: clusterProfiler function that tests for over-representation
  # of KEGG pathways. Downloads pathway-gene mappings from the KEGG REST API
  # (requires internet access). organism = "ath" for Arabidopsis thaliana.
  cat("\n--- KEGG pathway enrichment ---\n")
  ekg <- tryCatch({
    enrichKEGG(gene = gene_list, organism = KEGG_ORG, universe = universe,
               pAdjustMethod = "BH", pvalueCutoff = PADJ_CUTOFF,
               keyType = "kegg")
  }, error = function(e) {
    cat("  KEGG error:", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(ekg) && nrow(as.data.frame(ekg)) > 0) {
    out_csv <- file.path(OUTPUT_DIR, paste0("KEGG_enrichment_", label, ".csv"))
    write.csv(as.data.frame(ekg), out_csv, row.names = FALSE)
    cat("  Saved:", basename(out_csv), "(", nrow(as.data.frame(ekg)), "terms )\n")

    p <- dotplot(ekg, showCategory = 15) + ggtitle(paste0(label, " -- KEGG"))
    ggsave(file.path(OUTPUT_DIR, paste0("KEGG_dotplot_", label, ".png")), p, width = 9, height = 7)
    cat("  Saved dotplot\n")
  } else {
    cat("  No enriched KEGG pathways found.\n")
  }
  results$ekg <- ekg

  # --- 3. KEGG module enrichment ---
  # enrichMKEGG: clusterProfiler function that tests for over-representation
  # of KEGG modules (functional units within pathways). Also downloads from
  # KEGG REST API. Some small gene lists may return "No gene can be mapped"
  # if none of the genes are in KEGG module annotations.
  cat("\n--- KEGG module enrichment ---\n")
  emk <- tryCatch({
    enrichMKEGG(gene = gene_list, organism = KEGG_ORG, universe = universe,
                pAdjustMethod = "BH", pvalueCutoff = PADJ_CUTOFF)
  }, error = function(e) {
    cat("  KEGG module error:", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(emk) && nrow(as.data.frame(emk)) > 0) {
    out_csv <- file.path(OUTPUT_DIR, paste0("KEGG_module_enrichment_", label, ".csv"))
    write.csv(as.data.frame(emk), out_csv, row.names = FALSE)
    cat("  Saved:", basename(out_csv), "(", nrow(as.data.frame(emk)), "terms )\n")

    p <- dotplot(emk, showCategory = 15) + ggtitle(paste0(label, " -- KEGG Module"))
    ggsave(file.path(OUTPUT_DIR, paste0("KEGG_module_dotplot_", label, ".png")), p, width = 9, height = 7)
    cat("  Saved dotplot\n")
  } else {
    cat("  No enriched KEGG modules found.\n")
  }
  results$emk <- emk

  invisible(results)
}

# ---- 4. Load data and run enrichment --------------------------------------
# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Load the annotated DEG list (from tissue_specific_deg_analysis_v2.R)
deg_list <- readRDS(DEG_LIST_PATH)
for (label in names(deg_list)) {
  if (!is.data.table(deg_list[[label]])) setDT(deg_list[[label]])
}

# Load sub-tissue statistics (from tissue_specific_deg_analysis_v2.R)
subtissue_stats <- readRDS(SUBTISSUE_STATS_PATH)

# Comparisons with significant sub-tissues
comparisons_to_run <- c("radiation_effect", "Cs137_100cGy", "Cs137_10cGy")

# Summary table
enrichment_summary <- data.table(
  comparison = character(), sub_tissue = character(), direction = character(),
  n_genes = integer(), n_go_terms = integer(), n_kegg_pathways = integer(),
  n_kegg_modules = integer()
)

# Run enrichment for each significant sub-tissue
for (comp in comparisons_to_run) {
  cat(sprintf("\n\n========== %s ==========\n", comp))

  dt <- deg_list[[comp]]
  universe <- dt[!is.na(padj), gene_id]  # all genes tested by DESeq2

  # Get significant sub-tissues for this comparison
  sig_st <- subtissue_stats[[comp]]$fisher[padj_bonf < 0.05, sub_tissue]
  cat("Significant sub-tissues:", paste(sig_st, collapse = ", "), "\n")

  for (st in sig_st) {
    # Get DEGs for this sub-tissue
    st_degs <- dt[deg_flag == "yes" & predominant_subtissue == st]

    all_genes  <- st_degs$gene_id
    up_genes   <- st_degs[direction == "up", gene_id]
    down_genes <- st_degs[direction == "down", gene_id]

    # Run enrichment for each direction with enough genes
    for (dir_name in c("all", "up", "down")) {
      genes <- switch(dir_name, all = all_genes, up = up_genes, down = down_genes)

      if (length(genes) < MIN_GENES) {
        cat(sprintf("\n  [%s / %s / %s]: %d genes -- skipping (< %d)\n",
                    comp, st, dir_name, length(genes), MIN_GENES))
        next
      }

      label <- paste0(comp, "_", st, "_", dir_name)
      cat(sprintf("\n  [%s / %s / %s]: %d genes\n", comp, st, dir_name, length(genes)))

      res <- run_subtissue_enrichment(genes, universe, label)

      # Record summary
      n_go   <- if (!is.null(res$ego)) nrow(as.data.frame(res$ego)) else 0
      n_kegg <- if (!is.null(res$ekg)) nrow(as.data.frame(res$ekg)) else 0
      n_mod  <- if (!is.null(res$emk)) nrow(as.data.frame(res$emk)) else 0

      enrichment_summary <- rbind(enrichment_summary, data.table(
        comparison = comp, sub_tissue = st, direction = dir_name,
        n_genes = length(genes), n_go_terms = n_go,
        n_kegg_pathways = n_kegg, n_kegg_modules = n_mod
      ))
    }
  }
}

# Save summary
fwrite(enrichment_summary, file.path(OUTPUT_DIR, "enrichment_summary_all.csv"))
cat("\n\n=== Enrichment complete ===\n")
cat("Total runs:", nrow(enrichment_summary), "\n")
cat("Runs with GO results:", sum(enrichment_summary$n_go_terms > 0), "\n")
cat("Runs with KEGG results:", sum(enrichment_summary$n_kegg_pathways > 0), "\n")
cat("Runs with KEGG module results:", sum(enrichment_summary$n_kegg_modules > 0), "\n")

# ---- 5. How to read the output --------------------------------------------
# Each enrichment run produces up to 3 CSV files and 3 dotplot PNGs:
#
# GO_enrichment_<label>.csv       - GO Biological Process terms
# KEGG_enrichment_<label>.csv     - KEGG pathway terms
# KEGG_module_enrichment_<label>.csv - KEGG module terms
#
# Columns in each CSV:
#   Description : pathway/term name
#   GeneRatio   : (your genes in this term) / (total genes in your list)
#   BgRatio     : (genes in this term overall) / (total genes in universe)
#   pvalue      : raw p-value from hypergeometric test
#   p.adjust    : adjusted p-value (Benjamini-Hochberg) -- filter on this < 0.05
#   qvalue      : another FDR estimate -- filter on this < 0.2
#   geneID      : your genes in this term (slash-separated)
#   Count       : number of your genes in this term
#
# For KEGG results, the Description column contains pathway names like
# "DNA replication", "Homologous recombination", etc.
# KEGG module results contain module names like "DNA repair proteins".
