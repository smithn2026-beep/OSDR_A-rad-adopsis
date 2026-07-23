# Plan: Targeted Analytical Strategy for Arabidopsis Radiation Response Publication

## Summary

Produce a publication-quality PDF strategy document that (1) synthesizes the user's existing GO enrichment results from OSD-498/510, (2) frames a focused literature synthesis on plant radiation biology identifying metadata variables and knowledge gaps, and (3) proposes 6 targeted analytical approaches with concrete R-based implementation guidance, expected outputs, and how each addresses a specific gap for a plant radiation biology publication.

## Data Context (from Phase 1 inspection)

**User's completed work:**
- GO enrichment from OSD-498/510 combined (comparison 5a: "radiation effect")
- Upregulated (226 GO terms, 1,578 genes): defense/immune, hypoxia, programmed cell death, salicylic acid response, hypersensitive response
- Downregulated (167 GO terms, 1,016 genes): photosynthesis, ribosome biogenesis, translation, cell cycle, chloroplast organization, rRNA processing
- Zero gene overlap between up/down — clean directional separation
- Pipeline documented in PDF with planned comparisons: 5b (sog1_1 × WT interaction), 5c (OSD-782 dose-response 10/100 cGy), 5d (OSD-658 GCR 40/80 cGy)

**Study landscape (from literature + GeneLab):**
- OSD-498: Arabidopsis, 6-day-old seedlings, gamma irradiation 100 Gy (Co60, 10 Gy/min), DDR study with sog1_1 vs WT
- OSD-510: Paired study, same design, enables genotype × radiation interaction analysis
- OSD-658: GCR simulation (proton/He/O/Si/Fe), 40 & 80 cGy, 10-day-old seedlings, 3hr post-irradiation (Dixit et al. 2023)
- OSD-782: Low-dose gamma IR, 10 & 100 cGy, 4-week-old rosettes, time course 1/3/24/72hr, multi-omics (RNA-seq + MAPit epigenomics) (Newman et al. 2023)

## Literature Synthesis (focused, ~18 key papers)

Key themes grounded in retrieved literature:
1. **SOG1 as master regulator**: 146 direct target genes spanning DNA repair (HR: RAD51, BRCA1, RAD17), cell cycle (SMR5/7, WEE1, CYCB1;1), AND defense/immune genes (FMO1, SAG101, WRKY50) — explains the user's co-occurring defense + DDR upregulation [14, 12, 15]
2. **Dose-response relationships**: GCR shows dose-dependent DNA repair upregulation; glucosinolate downregulation stronger at 40 cGy than 80 cGy [2]; low-dose gamma shows dose-dependent ethylene signaling and epigenetic changes [21]
3. **Cross-study confounders**: Assay type (microarray vs RNA-seq) is largest confound; radiation treatment is most correlated biological factor; hardware/lighting matter [4]
4. **Spaceflight vs radiation**: Photosynthesis downregulation and defense upregulation are shared across spaceflight and radiation studies [5, 17]
5. **Ethylene signaling**: Key pathway in low-dose IR response, with epigenetic regulation via ERF binding site methylation changes [21]
6. **Glucosinolate-immunity crosstalk**: Downregulation may trade off glutathione for antioxidant defense, with implications for pathogen resistance in space [2]

## PDF Document Structure

### Section 1: Executive Summary
- Current state of analysis, key findings from enrichment, what gaps remain

### Section 2: Literature Synthesis & Metadata Gap Analysis
- Focused review of ~18 papers organized by theme
- Table of metadata variables that influence radiation response (dose, dose rate, radiation type/quality, timepoint, tissue/developmental stage, genotype, growth conditions, assay type)
- Table of additional GeneLab datasets available for integration
- Knowledge gaps the user's data can address

### Section 3: Six Targeted Analytical Approaches

**Approach 1: SOG1 Regulatory Network Overlay**
- Map user's DEGs against the 146 known SOG1 direct targets [14]
- Test whether upregulated defense/immune genes are SOG1-dependent using the 5b comparison (sog1_1 vs WT)
- R implementation: gene set overlap test (Fisher's exact), Venn diagram, SOG1 motif enrichment (CTT(N)7AAG) in promoters of upregulated genes
- Expected output: SOG1 dependency map, motif enrichment plot

**Approach 2: Dose-Response Meta-Analysis Across Studies**
- Integrate OSD-498/510 (100 Gy gamma), OSD-782 (10/100 cGy gamma), OSD-658 (40/80 cGy GCR)
- Model dose-response curves for key DDR genes (RAD51, BRCA1, PARP1/2, TSO2, XRI1)
- R implementation: DESeq2 with dose as continuous covariate, rank-based cross-study normalization (ComBat-seq or RUVseq), dose-response curve fitting
- Expected output: dose-response curves, cross-study DEG overlap heatmap, conserved vs condition-specific gene lists

**Approach 3: Temporal Dynamics Reconstruction**
- Leverage OSD-782 time course (1/3/24/72hr) to map temporal ordering of the user's enriched pathways
- Track when defense, hypoxia, photosynthesis, and translation pathways engage/disengage
- R implementation: maSigPro or ImpulseDE2 for time-series DE, pathway-level temporal profiling, trajectory plots
- Expected output: temporal pathway activation map, early vs late response gene modules

**Approach 4: Cross-Study Consensus Signature**
- Identify a core radiation response signature conserved across all available studies
- Use WGCNA (as in Barker et al. 2023 [4]) to find co-expression modules
- Test signature against spaceflight datasets to distinguish radiation-specific from general spaceflight response
- R implementation: WGCNA across combined normalized datasets, module-trait correlation, consensus module identification
- Expected output: consensus radiation signature gene list, module heatmap, spaceflight vs radiation Venn

**Approach 5: Glucosinolate-Immunity Trade-off Analysis**
- The user's downregulated pathways include defense-related processes; literature shows glucosinolate downregulation under GCR [2]
- Examine whether glucosinolate downregulation trades off with glutathione/antioxidant upregulation
- Test for reciprocal regulation of glucosinolate biosynthesis vs glutathione metabolism
- R implementation: pathway-level expression correlation, gene set enrichment for glutathione metabolism, glucosinolate pathway map overlay
- Expected output: glucosinolate-glutathione trade-off plot, pathway cross-talk network

**Approach 6: Epigenetic Regulatory Landscape Integration**
- Integrate OSD-782 MAPit/bisulfite data (chromatin accessibility + DNA methylation) with user's DEGs
- Test whether upregulated genes show concordant chromatin opening and hypomethylation at ERF binding sites
- Identify transcription factor binding motifs enriched in differentially accessible regions near user's DEGs
- R implementation: overlap DEGs with DAR/DMR-associated genes, motif enrichment (HOMER or monaLisa), integrative visualization
- Expected output: epigenetic-transcriptomic concordance plot, TF motif enrichment in accessible regions

### Section 4: Metadata Variables to Address
- Structured table of variables, why they matter, which analysis addresses them

### Section 5: Recommended Analysis Sequence
- Prioritized order: 1 → 2 → 4 → 3 → 5 → 6
- Dependencies, estimated compute, which can run in parallel

### Section 6: Publication Framework
- Suggested figure list (6-8 figures)
- Suggested table list
- Target journals and framing options

## Implementation

1. Generate the PDF using ReportLab (python) via the pdf-report-generation skill
2. Include 2 summary tables (metadata variables, dataset inventory)
3. Include 1 conceptual figure showing the analytical framework (GenerateImage)
4. Save to /mnt/results/

## Assumptions
- User has access to DEG tables and raw counts for OSD-498/510 (confirmed: "full pipeline access")
- User can download OSD-782 and OSD-658 from GeneLab (confirmed: "can download from GeneLab")
- R is preferred language (user preference)
- PNG is preferred image format (user preference)
- Publication target: plant radiation biology journal (e.g., Plant Cell, Plant Journal, Frontiers in Plant Science, International Journal of Radiation Biology)
