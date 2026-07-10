# Plan: Tissue-Specific Separation of Arabidopsis Radiation DEGs + Educational PDF Manual

## Summary

The input file (`DEG_OSD498_510_radiation_effect.csv`) is a single DESeq2 results table: 23,573 Arabidopsis genes (AGI IDs) with baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, and a yes/no DEG flag (6,942 DEGs). It contains **no tissue column** — the source studies (NASA GeneLab OSD-498/OSD-510) used whole seedlings exposed to ionizing radiation.

To separate DEGs into tissue-specific responses, we will annotate each gene by its tissue-specific expression pattern using the **AtGenExpress developmental expression atlas** (Schmid et al. 2005, *Nature Genetics* [26, 35]), the canonical Arabidopsis tissue expression reference. We compute the **Tau tissue-specificity index** [41] per gene, assign each DEG to its predominant tissue, and organize results hierarchically (broad organs → sub-tissues).

## Deliverables

1. **R script** (`tissue_specific_deg_analysis.R`) — standalone, heavily commented, reproducible
2. **Tissue-specific CSV files** — hierarchical folder structure (broad organ → sub-tissue)
3. **Summary table** — DEG counts per tissue with up/down breakdown
4. **Visualizations** (PNG): bar chart, volcano plots per tissue, heatmap, UpSet diagram
5. **PDF instruction manual** — mixed-level educational document explaining the biology, the code line-by-line, and how to run it

## Approach

### Step 1: Load and validate the DEG dataset
- Read the CSV; confirm 23,573 genes, 8 columns, 6,942 DEGs (flag = "yes")
- Clean column names (first column is unnamed — it's the gene ID)
- Parse AGI IDs (AT{1-5}Gnnnnn format) and validate

### Step 2: Build the tissue expression reference from AtGenExpress
- **Primary approach**: Download the AtGenExpress developmental atlas (GEO: GSE5629 and related series GSE5630/5632/5633) via `GEOquery`, map Affymetrix ATH1 probe IDs to AGI gene IDs using the `ath1121501.db` Bioconductor annotation package, and aggregate expression by tissue
- **Fallback** (if GEO download is unavailable in sandbox): Build a pre-processed tissue expression matrix from the AtGenExpress supplementary data or a curated subset, bundled as a CSV
- Organize samples into a **hierarchical tissue structure**:

| Broad organ | Sub-tissues (from AtGenExpress samples) |
|---|---|
| Root | root tip, elongation zone, mature root, lateral root |
| Leaf/Shoot | cotyledon, hypocotyl, rosette leaf, cauline leaf, stem |
| Flower | sepal, petal, stamen, carpel, pollen, flower bud |
| Seed/Silique | early silique, late silique, dry seed, germinating seed |
| Seedling | young seedling, whole seedling |

### Step 3: Compute tissue-specificity (Tau index)
- For each gene, compute the **Tau index** [41]: τ = Σ(1 − xᵢ/max(x)) / (n − 1), where xᵢ is mean expression in tissue i and n is the number of tissues
- τ ranges from 0 (ubiquitous expression) to 1 (perfectly tissue-specific)
- **Assignment logic**:
  - τ ≥ 0.6 → tissue-specific; assign to the tissue with the highest expression contribution
  - τ < 0.6 → "constitutive/general" (expressed broadly across tissues)
- This produces a hierarchical assignment: each DEG gets a broad organ + sub-tissue label (or "constitutive")

### Step 4: Export tissue-specific CSV files
- Create folder structure: `tissue_specific_degs/{broad_organ}/{sub_tissue}_DEGs.csv`
- Each CSV contains all original DESeq2 columns + tissue, sub_tissue, tau, predominant_tissue_expression
- Also export a `constitutive_DEGs.csv` for broadly-expressed DEGs
- Export a master `all_degs_with_tissue_annotation.csv`

### Step 5: Generate summary table
- `tissue_deg_summary.csv`: tissue, sub_tissue, total_DEGs, upregulated, downregulated, median_log2FC, mean_tau

### Step 6: Visualizations (all PNG)
1. **Bar chart**: DEG counts per broad organ, stacked by up/down regulation
2. **Volcano plots**: one per broad organ (facet or multi-panel), showing log2FC vs −log10(padj) with tissue colors
3. **Heatmap**: top DEGs (by padj) × tissues, showing expression specificity (ComplexHeatmap)
4. **UpSet diagram**: DEG overlaps across broad organs (which DEGs are tissue-specific to multiple organs)

### Step 7: PDF instruction manual
- Use the `pdf-report-generation` skill for a professional Phylo-branded PDF
- **Mixed-level structure**:
  - Main text: intermediate level (assumes basic R, explains biological concepts)
  - Sidebars/callout boxes: beginner-friendly explanations of key terms (DEG, log2FC, padj, Tau index, expression atlas)
  - Advanced notes: parameter choices, alternative approaches, caveats
- **Sections**:
  1. Introduction: The experiment (Arabidopsis + radiation, OSD-498/510 context)
  2. Background: What is differential expression? What is a DESeq2 results table?
  3. The tissue-specificity concept: Why use an expression atlas? What is the Tau index?
  4. Code walkthrough: Section-by-section explanation of every code block
  5. Understanding the outputs: What each CSV and plot tells you
  6. How to run the code: Prerequisites, installation, execution
  7. Interpretation guide: What do tissue-specific radiation responses mean biologically?
  8. Exercises & extensions: Ideas for students to explore further

## Compute/Resource Estimate
- **Input**: 2 MB DEG CSV + ~50-100 MB AtGenExpress atlas download
- **Processing**: 23K genes × ~80 samples — trivial (seconds to minutes)
- **Plots**: standard ggplot2/ComplexHeatmap/UpSetR — seconds each
- **PDF**: ReportLab generation — seconds
- **Execution target**: default machine (worker-0) — no special provisioning needed
- **Estimated total runtime**: 5-15 minutes

## Assumptions
1. AtGenExpress developmental atlas (Schmid et al. 2005) is the appropriate tissue reference for these whole-seedling radiation experiments
2. Tau index threshold of 0.6 distinguishes tissue-specific from constitutive genes (standard in the literature [41])
3. The 6,942 genes flagged "yes" in the DEG column are the DEGs of interest
4. AGI gene IDs in the DEG file map directly to AtGenExpress probe sets via the ATH1 annotation
5. Genes not present in the AtGenExpress atlas will be flagged as "unannotated" rather than dropped
