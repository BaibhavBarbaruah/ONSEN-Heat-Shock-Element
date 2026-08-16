# ONSEN HSE regulatory-architecture analyses

This repository contains analysis code, metadata and compact source data supporting the revised Biology Open manuscript:

**Heat-responsive ONSEN long terminal repeats integrate heat shock factor motifs, DNA methylation and natural sequence variation in Arabidopsis**  
Manuscript ID: **bio.062799**

Authors: Baibhav R. Barbaruah, Rahmadani P. Airalangga and Hidetaka Ito.

## August 2026 final-revision status

The August 2026 revision uses the following journal-facing numbering:

- main figures: **Fig. 1-Fig. 7**
- supplementary figures: **Fig. S1-Fig. S4**
- supplementary tables: **Table S1-Table S14**

The previous AP2/ERF-only supplementary table is **not part of the final submission**. The editor-requested genome-wide gene/TE differential-expression workbook is final **Table S14**.

### Final supplementary-table order

| Final item | Analysis |
|---|---|
| Table S1 | Native versus designed 49-bp motif effects and family-level non-redundant motif-coordinate placements |
| Table S2 | Arabidopsis HSF PWM models falling below the 0.85 threshold in the designed sequence |
| Table S3 | HSF-focused constrained-mutant and exact-GC robustness control |
| Table S4 | Col-0 terminal-window non-redundant HSF motif-coordinate summary |
| Table S5 | Principal exact-coordinate HSF placements and strict overlap-merged sensitivity analysis |
| Table S6 | Broad strict non-ONSEN TE background comparison |
| Table S7 | Threshold sensitivity and continuous HSF-score analysis |
| Table S8 | Direct non-ONSEN LTR-retrotransposon methylation comparison and sensitivity controls |
| Table S9 | Corrected complete-TAIR10 candidate-window RNA-seq summary |
| Table S10 | Natural-accession candidate abundance, HSF density and structural proxies |
| Table S11 | Natural 49-bp ONSEN-like HSE seed-variant summary |
| Table S12 | Natural seed-variant distinct TF-family PWM-model summary |
| Table S13 | Published accession-evidence comparison |
| Table S14 | Genome-wide gene and transposable-element differential expression after heat stress |

`FINAL_NUMBERING_MAP.tsv` is the authoritative old-to-final map for the August 2026 package.

### Final supplementary-figure order

| Final item | Analysis |
|---|---|
| Fig. S1 | Principal non-redundant HSF motif-coordinate placements versus strict overlap-merged HSF regions |
| Fig. S2 | HSF threshold sensitivity and continuous regional-score analysis |
| Fig. S3 | Additional accession-level HSE architecture and HSF-model compatibility summaries |
| Fig. S4 | Accession candidate abundance, HSF density and structural-proxy summaries |

## Important metric definitions

The **principal motif-count metric** in the final revision is a physical-forward exact-coordinate metric. Reverse-strand predictions are mapped to forward-sequence coordinates and duplicate predictions with identical 1-based start/end coordinates are collapsed across related PWM models and strands. Predictions at different coordinates remain separate.

The **strict overlap-merged** HSF metric transitively merges intervals that overlap by at least 1 nt and is retained only as a stringent sensitivity analysis.

All motif results describe **predicted sequence/PWM compatibility**, not in-vivo binding, occupancy, affinity or functional necessity.

Key corrected values used in the final revision include:

- native 49-bp ONSEN HSE window: **20** principal HSF motif-coordinate placements; designed HSE-disrupted sequence: **0**;
- sixteen Col-0 800-bp terminal candidate windows: **33-49** principal HSF placements per window;
- broad strict TE-background comparison at relative PWM score 0.85: median **49.375 placements/kb** in ONSEN terminal windows versus **6.812 placements/kb** in 1,930 strict non-ONSEN TE regions;
- direct element/copy HSF comparison: eight ONSEN elements versus 794 non-ONSEN LTR-retrotransposons, median **51.2 versus 8.12 placements/kb**;
- direct basal CHH methylation comparison: ONSEN median **47.09%** versus **20.06%** in the annotation-defined non-ONSEN LTR-retrotransposon comparator; the public methylome is from unstressed Col-0 leaves and is not a heat-methylation dataset;
- corrected genome-wide heat-response analysis: 32,833 annotated genes and 31,189 annotated transposable elements; 3,912 genes/4,185 genes were significantly up/down and 692 TEs/397 TEs were significantly up/down under the final significance definition.

## Repository contents

### Workflow

- `00_install_packages.R` — installs/checks required R packages
- `00_run_pipeline.R` — runs the analysis workflow in manuscript order
- `ONSEN_config.R` — portable data/output configuration
- `ONSEN_functions.R` — shared functions
- `01_native_mutated_motif_analysis.R` — native versus designed 49-bp motif analysis
- `02_constrained_mutant_sensitivity.R` — constrained-mutant and exact-GC design-space analysis
- `03_col0_HSF_and_TE_background.R` — Col-0 terminal windows and TE-background analyses
- `03B_threshold_and_continuous_sensitivity.R` — threshold and continuous HSF-score analysis
- `04_nonredundant_HSF_locations.R` — non-redundant/overlap sensitivity analysis
- `05_methylation_analysis.R` — basal methylation analyses
- `06_rnaseq_analysis.R` — RNA-seq analyses
- `07_accession_analysis.R` — accession candidate and natural-variant analyses
- `07B_figure7_identity_and_logo_workflow.R` — Fig. 7 sequence identity and representative motif logos
- `08_write_supplementary_tables.R` — legacy deposited source-sheet rebuild helper; see note below
- `09_write_session_info.R` — captures R/package versions
- `10_validate_manuscript_outputs.R` — repository-level metadata/provenance checks

### Legacy source-sheet numbering

The files under `supplementary_table_source/` were frozen before the final August 2026 deletion/renumbering pass and therefore retain their **legacy source-sheet numbers** for provenance. They must not be interpreted as the final journal table numbers. The final submission Excel workbooks were generated/audited separately from the corrected analysis outputs. Use `FINAL_NUMBERING_MAP.tsv` to translate legacy source-sheet numbers to the final journal numbering.

This distinction intentionally preserves provenance instead of silently renaming historical source files.

## Configuration

Public-facing scripts use configurable data/output roots rather than machine-specific drive paths. Set:

```r
Sys.setenv(
  ONSEN_DATA_ROOT = "path/to/ONSEN_HSE_input_data",
  ONSEN_OUTPUT_ROOT = "path/to/ONSEN_HSE_outputs"
)
```

The complete raw-data workflow requires the external inputs described in `INPUT_PROVENANCE.tsv`.

## Motif-scoring definition

Position frequency matrices were converted to log2 probability-ratio PWMs using:

- pseudocount: `0.8`
- equal A/C/G/T background: `0.25`
- both DNA strands
- relative score: `(raw score - theoretical minimum)/(theoretical maximum - theoretical minimum)`

The primary operational high-confidence threshold is `0.85`. Sensitivity analyses additionally use `0.80`, `0.90` and `0.95`, together with threshold-independent regional score summaries.

Representative Fig. 7D logos are illustrative and are distinct from the JASPAR 2026 HSFC1 model used in the primary HSF-focused scan:

- HSF family: **HSFC1, JASPAR MA1667.2**
- DOF family: **DOF1.8, JASPAR MA0981.2**

## RNA-seq data availability

All six RNA-seq libraries were generated in our laboratory, with three biological replicates per condition.

- Col-0 non-stressed: BioProject **PRJDB39904**; BioSamples **SAMD01789795, SAMD01789796, SAMD01789797**
- Col-0 37°C, 24-h heat stress: BioProject **PRJDB42759**; BioSamples **SAMD01943917, SAMD01943918, SAMD01943919**

Sequencing accessions, public resources and expected filenames are documented in `DATA_AVAILABILITY.md` and `INPUT_PROVENANCE.tsv`.

## Citation and licence

Citation metadata are provided in `CITATION.cff`. Code is released under the MIT License. JASPAR motif data retain their original attribution and licence; see `motifs/README.md`.
