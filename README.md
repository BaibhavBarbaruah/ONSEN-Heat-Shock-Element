# ONSEN HSE regulatory-architecture analyses

Code, configuration, provenance metadata and final supplementary source data supporting the Biology Open manuscript:

**Heat-responsive ONSEN long terminal repeats integrate heat shock factor motifs, DNA methylation and natural sequence variation in Arabidopsis**

Authors: Baibhav R. Barbaruah, Rahmadani P. Airlangga and Hidetaka Ito.

## Repository scope

The final supplementary set comprises **Fig. S1-Fig. S4** and **Table S1-Table S13**. The formatted `Table_S1.xlsx`-`Table_S13.xlsx` files at the repository root are the journal-facing workbooks. Data-equivalent plain-text mirrors of their worksheets are under `supplementary_table_source/`.

### Supplementary-table mapping

| Table | Analysis |
|---|---|
| S1 | Native versus designed 49-bp motif effects and TF-family non-redundant motif-coordinate placements |
| S2 | Arabidopsis HSF PWM models in the native and designed sequences |
| S3 | HSF-focused constrained-alternative and exact-GC robustness analysis |
| S4 | Col-0 ONSEN terminal-window HSF motif-coordinate summary |
| S5 | Exact-coordinate HSF placements and strict overlap-merged sensitivity analysis |
| S6 | Strict non-ONSEN TE background comparison |
| S7 | Four-threshold and continuous HSF-score sensitivity analyses |
| S8 | Basal DNA-methylation comparisons and sensitivity analyses |
| S9 | Genome-wide gene and transposable-element differential expression after 24-h 37°C heat stress |
| S10 | Candidate-window RNA-seq signal and replicate-summed Fig. 5B statistics |
| S11 | Natural-accession candidate abundance, HSF density and structural proxies |
| S12 | Natural 49-bp ONSEN-like HSE seed variants |
| S13 | Natural seed-variant TF-family PWM-model summary |

## Principal analysis definitions

### HSF motif scanning

JASPAR plant PFMs are converted to log2 probability-ratio PWMs using a pseudocount of 0.8 and equal A/C/G/T background frequencies. Both strands are scanned and predictions are mapped to 1-based forward-sequence coordinates.

The **principal non-redundant HSF metric** is the number of unique `(start, end)` coordinate placements. Predictions from different HSF PWM models or strands with identical start and end coordinates are counted once; predictions at different coordinates remain separate. A stricter sensitivity metric transitively merges HSF intervals whenever they overlap by at least one nucleotide. These are sequence-compatibility metrics, not occupancy measurements.

The primary relative PWM-score threshold is **0.85**. Genome-wide sensitivity analysis additionally evaluates **0.80, 0.90 and 0.95**, together with three threshold-independent regional summaries: maximum HSF relative score, mean of the five highest HSF relative scores and median of the ten per-model maximum scores.

### Col-0 terminal-window comparison

The main genome-wide comparison uses sixteen fixed 800-bp Col-0 ONSEN terminal candidate windows and **1,930** strict non-ONSEN TE regions after annotation and ONSEN/ATCOPIA78 harmonisation. At the 0.85 threshold, the final median densities are **49.375 placements/kb** for ONSEN terminal windows and **6.81238953796413 placements/kb** for the strict non-ONSEN TE background.

### DNA methylation

The public methylome is GEO **GSE43857 / GSM1085222**, generated from unstressed Col-0 leaves. The principal comparison aggregates methylated and total base counts across paired terminal regions at the element/copy level and compares ONSEN with **779 annotation-defined non-ONSEN LTR retrotransposons**. These data describe basal methylation context; they do not measure heat-induced methylation dynamics.

### RNA-seq

All six RNA-seq libraries were generated in our laboratory, with three biological replicates per condition.

- Col-0 non-stressed: BioProject **PRJDB39904**; BioSamples **SAMD01789795, SAMD01789796, SAMD01789797**
- Col-0 37°C, 24-h heat stress: BioProject **PRJDB42759**; BioSamples **SAMD01943917, SAMD01943918, SAMD01943919**

Gene-level differential expression uses DESeq2. The complete gene table retains all **32,833** annotated genes; 25,912 have non-zero signal and 23,922 have non-missing BH-adjusted P values after DESeq2 independent filtering.

Global transposable-element analysis retains all **31,189** Araport11 `transposable_element` loci. Multimapping and overlapping reads are fractionally assigned; expression filtering selects **2,101** loci for limma-voom testing. Final Up/Down calls require BH-adjusted P<0.05 and |log2 fold change|>=1.

Candidate-window measurements use fractional multimapping-aware counts normalized to assigned gene reads. For Fig. 5B, CPM is **summed across the sixteen curated ONSEN terminal windows within each biological replicate** and the three NS sums are compared with the three HS sums by a two-sided Welch's t-test. Table S10 also reports the corresponding Wilcoxon rank-sum result. Candidate-window signal is not interpreted as definitive copy-specific ONSEN expression.

### Natural-accession survey

The accession analysis uses chromosome-level assemblies for An-1, C24, Cvi-0, Eri-1, Kyo, Ler-1 and Sha together with Col-0. Candidate windows are recovered from the Col-0 49-bp HSE seed allowing at most four mismatches. Reported motif-retention and sequence-identity conclusions therefore apply to the detected candidate set rather than every possible ONSEN-related sequence in each accession.

The family-level PWM-model summaries in Fig. 7A,B and Table S13 are analytically distinct from the ten-model Arabidopsis HSF exact-coordinate analyses. They report distinct JASPAR PWM model identities compatible with each 49-bp seed sequence and should not be interpreted as independent motif-coordinate placements or binding sites.

## Workflow

- `00_install_packages.R` — package installation/check
- `00_run_pipeline.R` — master runner
- `ONSEN_config.R` — portable input/output configuration
- `ONSEN_functions.R` — shared sequence, motif and statistics utilities
- `01_native_mutated_motif_analysis.R` — Fig. 1-2; Tables S1-S2
- `02_constrained_mutant_sensitivity.R` — Table S3
- `03_col0_HSF_and_TE_background.R` — Fig. 3; Tables S4 and S6
- `03B_threshold_and_continuous_sensitivity.R` — Fig. S2; Table S7
- `03C_direct_LTR_HSF_comparison.R` — Fig. 4C direct LTR HSF-density comparison
- `04_nonredundant_HSF_locations.R` — Fig. S1; Table S5
- `05A_direct_LTR_methylation_analysis.R` — Fig. 4A; principal Table S8 comparison
- `05_methylation_analysis.R` — Fig. 4B,D and methylation sensitivity analyses
- `06_rnaseq_analysis.R` — Fig. 5; Tables S9-S10
- `07_accession_analysis.R` — Fig. 6-7 and accession supplementary analyses; Tables S11-S13
- `07B_figure7_identity_and_logo_workflow.R` — Fig. 7 sequence-identity/logo workflow
- `08_write_supplementary_tables.R` — data-equivalent workbook rebuilds from final TSV mirrors
- `09_write_session_info.R` — R/package versions
- `10_validate_manuscript_outputs.R` — repository consistency checks

The lightweight default workflow validates and uses deposited processed outputs where a full genome-wide rescan would require large external inputs. Set `ONSEN_RUN_LARGE_STEPS=true` when the inputs listed in `INPUT_PROVENANCE.tsv` are available and a complete re-analysis is intended.

## Configuration

Set external input and output roots with environment variables:

```r
Sys.setenv(
  ONSEN_DATA_ROOT = "path/to/ONSEN_HSE_input_data",
  ONSEN_OUTPUT_ROOT = "path/to/ONSEN_HSE_outputs"
)
```

Large sequencing/reference inputs are not duplicated in the repository. Their expected filenames, public accessions and analysis roles are listed in `INPUT_PROVENANCE.tsv` and `DATA_AVAILABILITY.md`.

## Citation and licence

Citation metadata are provided in `CITATION.cff`. Code is released under the MIT License. JASPAR motif data retain their original attribution and licence; see `motifs/README.md`.
