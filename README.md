# ONSEN HSE regulatory-architecture analyses

This repository contains the analysis code, compact metadata, journal-numbered supplementary-table source sheets and source data supporting the revised Biology Open manuscript:

**Heat-responsive ONSEN long terminal repeats integrate heat shock factor motifs, DNA methylation and natural sequence variation in Arabidopsis**  
Manuscript ID: **bio.062799**

Authors: Baibhav R. Barbaruah, Rahmadani P. Airalangga and Hidetaka Ito.

## Scope and final journal numbering

The repository is synchronized to the final revision package:

- main figures: **Fig. 1-Fig. 7**
- supplementary figures: **Fig. S1-Fig. S5**
- supplementary tables: **Table S1-Table S14**

The final supplementary order is:

| Final item | Analysis |
|---|---|
| Fig. S1 / Table S4 | constrained-mutant and exact-GC sensitivity |
| Fig. S2 / Table S6 | raw versus non-redundant HSF-compatible locations |
| Fig. S3 / Table S8 | HSF threshold sensitivity and continuous regional scores |
| Fig. S4 | additional accession-level HSE architecture, scaled summaries and HSF-model compatibility |
| Fig. S5 / Table S11 | accession candidate abundance, HSF density and structural proxies |
| Table S1-S3 | native-versus-designed 49-bp motif analyses |
| Table S5 | Col-0 terminal-window HSF summary |
| Table S7 | strict TE-background comparison |
| Table S9 | basal methylation and ordinary-TE control |
| Table S10 | RNA-seq candidate-window analysis |
| Table S12-S13 | natural seed variants and TF-family summaries |
| Table S14 | published accession-evidence synthesis |

`FINAL_NUMBERING_MAP.tsv` records the complete old-to-final supplementary renumbering used during revision.

## Repository contents

### Workflow

- `00_install_packages.R` — installs/checks required R packages
- `00_run_pipeline.R` — runs the analysis workflow in manuscript order
- `ONSEN_config.R` — portable data/output configuration
- `ONSEN_functions.R` — shared functions
- `01_native_mutated_motif_analysis.R` — native versus designed 49-bp motif analysis
- `02_constrained_mutant_sensitivity.R` — 5,000 constrained mutants and complete 5,120-sequence exact-GC design space
- `03_col0_HSF_and_TE_background.R` — Col-0 terminal windows and strict TE background
- `03B_threshold_and_continuous_sensitivity.R` — wrapper for the four-threshold and continuous HSF-score analysis
- `03B_threshold_and_continuous_sensitivity_part1.R` / `part2.R` — implementation of the Reviewer 1.3 analysis
- `04_nonredundant_HSF_locations.R` — merged non-redundant HSF-compatible intervals
- `05_methylation_analysis.R` — basal methylation and ordinary-TE control
- `06_rnaseq_analysis.R` — gene and multimapping-aware candidate-window RNA-seq analysis
- `07_accession_analysis.R` — accession candidate and natural-variant analyses
- `07B_figure7_identity_and_logo_workflow.R` — Fig. 7C sequence identity and Fig. 7D representative motif logos
- `08_write_supplementary_tables.R` — rebuilds Table S1-Table S14 Excel workbooks from deposited TSV source sheets and generated outputs
- `09_write_session_info.R` — captures R/package versions
- `10_validate_manuscript_outputs.R` — checks structure, numbering, checksums and key manuscript values

### Deposited final items

- `supplementary_table_source/` — final journal-numbered source sheets used to rebuild Tables S1-S14
- `source_data/` — exact compact source tables for Fig. 7 and the final threshold/continuous-score analysis
- `motifs/` — exact representative JASPAR PFMs used for the Fig. 7D logo workflow

The journal-ready figure PDFs are submission assets and are not duplicated in this code/data repository. Every final figure and supplementary figure is mapped to its responsible script and source data in `REPRODUCIBILITY_MATRIX.tsv`.

`FILE_CHECKSUMS_SHA256.tsv` records SHA-256 checksums for deposited source sheets, source tables and motif matrices.

## Configuration

The scripts contain no machine-specific drive paths. Set one input root containing the external and processed files listed in `INPUT_PROVENANCE.tsv`:

```r
Sys.setenv(
  ONSEN_DATA_ROOT = "path/to/ONSEN_HSE_input_data",
  ONSEN_OUTPUT_ROOT = "path/to/ONSEN_HSE_outputs"
)
```

Defaults:

- `ONSEN_DATA_ROOT`: repository working directory
- `ONSEN_OUTPUT_ROOT`: `reproduced_outputs` under the repository working directory

Optional execution controls:

```r
Sys.setenv(
  ONSEN_FORCE_RESCAN = "false",
  ONSEN_RUN_LARGE_STEPS = "false",
  ONSEN_MAKE_FIGURES = "true"
)
```

- `ONSEN_FORCE_RESCAN=true` recomputes motif scans instead of preferring processed project outputs.
- `ONSEN_RUN_LARGE_STEPS=true` permits computationally intensive genome-wide scans and alignment/counting steps.
- `ONSEN_MAKE_FIGURES=false` suppresses plotting while retaining tabular analyses.

## Running the workflow

From the repository root:

```r
source("00_install_packages.R")
source("00_run_pipeline.R")
```

Individual order:

```r
source("01_native_mutated_motif_analysis.R")
source("02_constrained_mutant_sensitivity.R")
source("03_col0_HSF_and_TE_background.R")
source("03B_threshold_and_continuous_sensitivity.R")
source("04_nonredundant_HSF_locations.R")
source("05_methylation_analysis.R")
source("06_rnaseq_analysis.R")
source("07_accession_analysis.R")
source("07B_figure7_identity_and_logo_workflow.R")
source("08_write_supplementary_tables.R")
source("09_write_session_info.R")
source("10_validate_manuscript_outputs.R")
```

The complete raw-data workflow requires the external inputs described in `INPUT_PROVENANCE.tsv`. When large inputs are unavailable, the deposited final source sheets and source-data tables preserve the exact numerical revision outputs.

## Validation

Run:

```r
source("10_validate_manuscript_outputs.R")
```

The validator checks:

- all required scripts and metadata;
- final Fig. S1-Fig. S5 and Table S1-Table S14 numbering;
- absence of hard-coded Windows drive paths;
- SHA-256 checksums for deposited source sheets, source tables and motif matrices;
- selected manuscript-level values across Tables S1, S4-S14;
- exact Fig. 7D logo model names and JASPAR IDs.

A successful run ends with:

```text
REPOSITORY VALIDATION PASSED
```

## Motif-scoring definition

Position frequency matrices were converted to log2 probability-ratio PWMs using:

- pseudocount: `0.8`
- equal A/C/G/T background: `0.25`
- both DNA strands
- relative score: `(raw score - theoretical minimum)/(theoretical maximum - theoretical minimum)`

The primary high-confidence threshold is `0.85`. Robustness analyses additionally use `0.80`, `0.90` and `0.95`, together with threshold-independent regional score summaries.

The representative Fig. 7D logos are illustrative and are distinct from the JASPAR 2026 HSFC1 MA1667.3 model used in the primary HSF-focused scan. They use:

- HSF family: **HSFC1, JASPAR MA1667.2**
- DOF family: **DOF1.8, JASPAR MA0981.2**


## Data availability

All six RNA-seq libraries were generated in our laboratory, with three biological replicates per condition.

The RNA-seq data generated in this study have been deposited in DDBJ. The three Col-0 non-stressed libraries are available under BioProject PRJDB39904, with accession numbers SAMD01789795, SAMD01789796 and SAMD01789797. The three Col-0 37°C, 24-h heat-stress libraries are available under BioProject PRJDB42759, with accession numbers SAMD01943917, SAMD01943918 and SAMD01943919.

Sequencing accessions, public resources and expected filenames are documented in `DATA_AVAILABILITY.md` and `INPUT_PROVENANCE.tsv`.

## Citation and license

Citation metadata are provided in `CITATION.cff`. Code is released under the MIT License. JASPAR motif data retain their original attribution and licence; see `motifs/README.md`.
