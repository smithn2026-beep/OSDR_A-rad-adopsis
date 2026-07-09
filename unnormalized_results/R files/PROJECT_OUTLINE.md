# Project Outline: Multi-Study Arabidopsis Ionizing Radiation RNA-seq Analysis

## Dataset summary

Six independent NASA GeneLab/OSDR studies, all *Arabidopsis thaliana*, all examining
transcriptional response to ionizing radiation. **158 samples total.**

| Study  | n   | Genotype(s)     | Radiation source       | Dose            | Timepoints                  | Design type |
|--------|-----|-----------------|-------------------------|-----------------|------------------------------|---|
| OSD498 | 16  | WT              | Co-60 gamma             | 100 Gy          | 10, 20, 90, 1440 min          | Time course |
| OSD502 | 8   | WT, myb3r135    | Co-60 gamma             | 100 Gy          | 3h only                       | Genotype x radiation (single tp) |
| OSD508 | 36  | WT, sog1_1      | Co-60 gamma             | 100 Gy          | 0, 10, 20, 45, 90, 180, 360, 720, 1440 min | Time course x genotype |
| OSD510 | 48  | WT, sog1_1      | Co-60 gamma             | 100 Gy          | 20, 90, 180, 360, 720, 1440 min | Time course x genotype ("DREM" follow-up to 508) |
| OSD658 | 14  | WT              | Simulated GCR (mixed particles) | 0, 40, 80 cGy | none (single harvest)        | Dose response |
| OSD782 | 36  | WT              | Cs-137 gamma            | 0, 10, 100 cGy  | 1, 3, 24, 72 h                | Dose x time course |

**Key structural note:** WT is the only genotype present in all 6 studies. sog1_1 appears
only in OSD508/OSD510. myb3r135 appears only in OSD502. Radiation source/dose is NOT
consistent across all studies (Co-60 @ fixed 100Gy in 4 studies; Cs-137 dose-response in
OSD782; simulated GCR dose-response in OSD658) -- this is why subset-based comparisons
(see Phase 4 below) are recommended over one giant pooled model.

---

## Phase 1 — Data gathering ✅ COMPLETE
- [x] Counts table (merged normalized counts, 158 samples x ~23,573 genes)
- [x] ISA-tab metadata (`s_OSD-*.txt`) for all 6 studies, including the dose values
      pulled from `Parameter Value[absorbed radiation dose]` (not in the Factor Value
      column for 4 of the 6 studies — found in the sample table parameters instead)

## Phase 2 — Metadata cleaning ✅ COMPLETE
Inconsistencies identified and resolved:
- Radiation type/dose were entangled inconsistently across studies → split into a
  single standardized `Radiation` factor encoding both type and dose
  (e.g. `gammaCo60_100Gy`, `gammaCs137_10cGy`, `GCR_40cGy`, `none`)
- Genotype label `sog1-1` contained a hyphen → standardized to `sog1_1` to avoid
  parser issues in iDEP and R (hyphens inside factor levels can break automatic
  sample-name parsing)
- Time units were inconsistent (minutes in 4 studies, hours in OSD782, none in
  OSD658) → kept as originally reported (per your preference) but clearly labeled
  with units in every sample name (`10min`, `3h`, `NA`)

## Phase 3 — Human-readable sample naming ✅ COMPLETE
Convention: `{Study}_{Genotype}_{Radiation_Dose}_{Timepoint}_Rep{N}`

Examples:
- `OSD498_WT_gammaCo60_100Gy_10min_Rep1`
- `OSD508_sog1_1_none_0min_Rep2`
- `OSD782_WT_gammaCs137_100cGy_24h_Rep3`
- `OSD658_WT_GCR_80cGy_NA_Rep4`

All 158 names verified unique. Files delivered:
- `renamed_counts.csv` — counts table with original GSM/sample IDs replaced
- `factors_matrix.csv` — standalone factor table (Study, Genotype, Radiation,
  Timepoint, Replicate, OriginalID) for direct use in DESeq2/limma design formulas
- `sample_mapping.csv` — full audit trail: original ID → new name + every factor

## Phase 4 — Analysis (in progress)
**Script delivered:** `multifactor_analysis.R`

Planned structure:
1. **Exploratory PCA** colored by Study — checks whether study/batch effects dominate
   before trusting any pooled comparison
2. **Subset-based DE testing** (recommended primary approach):
   - OSD498 + OSD510 (WT, shared Co-60/100Gy time-course design) → radiation effect
     over time
   - OSD508 + OSD510 (WT vs sog1_1) → genotype x radiation interaction — does sog1_1
     blunt or alter the time-course response?
   - OSD782 alone → Cs-137 dose-response at each timepoint
   - OSD658 alone → simulated GCR dose-response
3. **Full pooled model** (Study + Genotype as covariates) — for exploratory
   clustering/heatmaps only, not for formal cross-study radiation hypothesis tests,
   since radiation type isn't comparable across all 6 studies
4. **iDEP3.0 cross-check** — upload `renamed_counts.csv` + `factors_matrix.csv`,
   confirm grouping/clustering matches the R output, then run pathway/GO enrichment
   on the gene lists from each DESeq2 contrast

## Phase 5 — Pathway / enrichment analysis (not started)
- GO/KEGG enrichment on each DEG list from Phase 4
- Cross-study comparison: which radiation-responsive pathways are shared vs.
  study-specific (e.g. is the SOG1-dependent DNA damage response consistent
  between OSD508 and OSD510?)

---

## Open questions to confirm with your supervisor
1. Should the **full pooled cross-study model** be reported at all, given the
   radiation-type heterogeneity, or used purely for exploratory visualization?
2. For OSD658 (no timepoint factor), should it be excluded from any combined
   time-course analysis, or treated as a single "endpoint" comparable to the
   longest timepoint in other studies?
3. Are normalized counts acceptable for DESeq2, or should we go back and pull the
   **raw/unnormalized counts tables** (`*_RSEM_Unnormalized_Counts.csv` /
   `*_STAR_Unnormalized_Counts.csv`, referenced in the assay files) for a more
   statistically standard DESeq2 run? DESeq2 performs its own normalization
   internally and technically expects raw counts as input.
