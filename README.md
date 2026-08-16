# ONSEN HSE regulatory-architecture analyses

This repository contains analysis code, metadata, final supplementary workbooks and machine-readable source sheets supporting the revised Biology Open manuscript:

**Heat-responsive ONSEN long terminal repeats integrate heat shock factor motifs, DNA methylation and natural sequence variation in Arabidopsis**  
Manuscript ID: **bio.062799**

Authors: Baibhav R. Barbaruah, Rahmadani P. Airlangga and Hidetaka Ito.

## Final August 2026 submission package

The journal-facing package uses:

- main figures: **Fig. 1-Fig. 7**;
- supplementary figures: **Fig. S1-Fig. S4**;
- supplementary tables: **Table S1-Table S13**.

There is no final Table S14. The AP2/ERF-only workbook and the earlier published-accession-evidence workbook were removed during the final revision. Genome-wide gene and transposable-element differential expression is final **Table S9**.

### Final supplementary-table order

| Final item | Analysis |
|---|---|
| Table S1 | Native versus designed 49-bp motif effects and family-level non-redundant motif-coordinate placements |
| Table S2 | Arabidopsis HSF PWM models falling below the 0.85 threshold in the designed sequence |
| Table S3 | HSF-focused constrained-alternative and exact-GC robustness analyses |
| Table S4 | Col-0 terminal-window non-redundant HSF motif-coordinate summary |
| Table S5 | Principal exact-coordinate HSF placements and strict overlap-merged sensitivity analysis |
| Table S6 | Broad strict non-ONSEN TE background comparison |
| Table S7 | Threshold sensitivity and continuous HSF-score analysis |
| Table S8 | Direct non-ONSEN LTR-retrotransposon methylation comparison and sensitivity controls |
| Table S9 | Genome-wide gene and transposable-element differential expression after heat stress |
| Table S10 | Candidate-window heat-associated RNA-seq signal and the Fig. 5B statistical comparison |
| Table S11 | Natural-accession candidate abundance, HSF density and structural proxies |
| Table S12 | Natural 49-bp ONSEN-like HSE seed-variant accession summary |
| Table S13 | Natural seed-variant distinct TF-family PWM-model summary |

The root `Table_S1.xlsx`-`Table_S13.xlsx` files are exact repository copies of the final supplementary workbooks. Every worksheet is also mirrored as a final-numbered TSV under `supplementary_table_source/`; `SOURCE_SHEET_MANIFEST.tsv` records the workbook-to-sheet mapping.

### Final supplementary-figure order

| Final item | Analysis |
|---|---|
| Fig. S1 | Principal exact-coordinate HSF placements versus strict overlap-merged HSF regions |
| Fig. S2 | HSF threshold sensitivity and continuous regional-score analysis |
| Fig. S3 | Additional natural-sequence-variation and PWM-model compatibility summaries |
| Fig. S4 | Accession candidate abundance, HSF density and structural-proxy summaries |

## Important metric definitions

The principal motif-count metric is a physical-forward exact-coordinate metric. Reverse-strand predictions are mapped to forward-sequence coordinates, and predictions with identical 1-based start/end coordinates are collapsed across related PWM models and strands. Predictions at different coordinates remain separate.

The strict overlap-merged metric transitively merges intervals overlapping by at least 1 nt and is retained only as a stringent sensitivity analysis.

All motif results describe predicted sequence/PWM compatibility, not in-vivo binding, occupancy, affinity or functional necessity.

Key checked values include:

- native 49-bp ONSEN HSE window: **20** principal HSF motif-coordinate placements; designed HSE-disrupted sequence: **0**;
- 16 Col-0 terminal candidate windows: **33-49** principal HSF placements per window;
- broad strict TE-background comparison at 0.85: median **49.375 placements/kb** in ONSEN windows versus **6.812 placements/kb** in 1,930 non-ONSEN TE regions;
- direct element-level HSF comparison: eight ONSEN elements versus 794 non-ONSEN LTR-retrotransposons, median **51.2 versus 8.12 placements/kb**;
- direct basal CHH methylation comparison: ONSEN median **47.09%** versus **20.06%** in 779 non-ONSEN LTR-retrotransposons;
- genome-wide heat-response universe: **32,833 genes** and **31,189 transposable elements**, with 3,912/4,185 genes and 692/397 TEs significantly up/down under the final definition.

## Repository workflow

- `00_run_pipeline.R` — runs analyses in manuscript order;
- `ONSEN_config.R` and `ONSEN_functions.R` — portable configuration and shared functions;
- `01_native_mutated_motif_analysis.R` — native versus designed 49-bp motif analysis;
- `02_constrained_mutant_sensitivity.R` — constrained-alternative and exact-GC design-space analysis;
- `03_col0_HSF_and_TE_background.R` — Col-0 terminal windows and broad TE-background analysis;
- `03B_threshold_and_continuous_sensitivity.R` — four-threshold and continuous-score analyses;
- `03C_direct_LTR_HSF_comparison.R` — direct eight-versus-794 element-level HSF-density comparison;
- `04_nonredundant_HSF_locations.R` — exact-coordinate and strict overlap sensitivity analysis;
- `05A_direct_LTR_methylation_analysis.R` — principal direct element/copy methylation comparison;
- `05_methylation_analysis.R` — broad methylation sensitivity and selected-locus/profile analyses;
- `06_rnaseq_analysis.R` — global gene/TE and candidate-window RNA-seq analyses;
- `07_accession_analysis.R` and `07B_figure7_identity_and_logo_workflow.R` — accession and Fig. 7 analyses;
- `08_write_supplementary_tables.R` — verifies/rebuilds final-numbered supplementary tables from deposited source sheets;
- `10_validate_manuscript_outputs.R` — validates final numbering, tables, metadata and provenance.

Set portable input/output roots as needed:

```r
Sys.setenv(
  ONSEN_DATA_ROOT = "path/to/ONSEN_HSE_input_data",
  ONSEN_OUTPUT_ROOT = "path/to/ONSEN_HSE_outputs"
)
```

The direct Fig. 4A and Fig. 4C scripts recompute element-level statistics when their processed element-level inputs are supplied. Compact final summaries are deposited under `source_data/`; the required processed-input names and provenance are stated in `INPUT_PROVENANCE.tsv`. The final Table S8 workbook is the authoritative submitted statistical record.

## RNA-seq and public methylome

All six RNA-seq libraries were generated in our laboratory, with three biological replicates per condition:

- Col-0 non-stressed: BioProject **PRJDB39904**; BioSamples **SAMD01789795, SAMD01789796, SAMD01789797**;
- Col-0 37°C, 24-h heat stress: BioProject **PRJDB42759**; BioSamples **SAMD01943917, SAMD01943918, SAMD01943919**.

Basal methylation context uses public unstressed Col-0 leaf methylome **GSE43857 / GSM1085222**. It is not a heat-methylation dataset and is not sample-matched to the RNA-seq experiment.

## Citation and licence

Citation metadata are provided in `CITATION.cff`. Code is released under the MIT License. JASPAR motif data retain their original attribution and licence; see `motifs/README.md`.
