# Data and code availability

## RNA-seq generated in this study

All six RNA-seq libraries were generated in our laboratory, with three biological replicates per condition.

The three Col-0 non-stressed libraries are available under BioProject **PRJDB39904**, with BioSample accessions **SAMD01789795, SAMD01789796 and SAMD01789797**. The three Col-0 37°C, 24-h heat-stress libraries are available under BioProject **PRJDB42759**, with BioSample accessions **SAMD01943917, SAMD01943918 and SAMD01943919**.

The repository does not duplicate FASTQ/BAM files. The sample-to-accession map and workflow roles are recorded in `RNAseq_sample_metadata_template.csv` and `INPUT_PROVENANCE.tsv`.

The final August 2026 RNA-seq analyses use the corrected complete-TAIR10 alignment/counting universe. Candidate-window measurements are described as fractional multimapping-aware **candidate-window signal**, not copy-specific ONSEN expression. Genome-wide gene differential expression was analysed with DESeq2; the final genome-wide TE workflow uses multimapping-aware fractional counts with expression filtering and limma-voom.

## Public methylome

- GEO series: **GSE43857**
- sample: **GSM1085222**
- material: unstressed *Arabidopsis thaliana* Col-0 leaf
- expected analysis filename: `GSM1085222_mC_calls_Col_0.tsv.gz`

This dataset is used only to describe basal methylation context. It is not a heat-treated methylome and is not interpreted as a measurement of heat-induced methylation change or as a sample-matched methylation-expression dataset.

## Reference resources

- *Arabidopsis thaliana* TAIR10 reference genome
- Araport11 gene and transposable-element annotation
- JASPAR CORE Plants motif collections used by the deposited workflows
- chromosome-level assemblies used for the accession survey

Expected filenames and analysis roles are listed in `INPUT_PROVENANCE.tsv`.

## Natural-accession evidence

The genome scan used a Col-0-derived 49-bp seed allowing up to four mismatches. More divergent ONSEN-related sequences may therefore be absent from the detected candidate set.

Published cross-accession ONSEN evidence is used as independent biological context rather than as functional validation of the present motif predictions. In the final journal numbering, this synthesis is **Table S13**. Historical deposited source filenames may still use the pre-final table number; `FINAL_NUMBERING_MAP.tsv` is authoritative.

## Final August 2026 supplementary package

The journal-facing supplementary package contains **Fig. S1-Fig. S4** and **Table S1-Table S14**. The AP2/ERF-only supplementary table was removed. Final Table S14 contains the editor-requested genome-wide gene and transposable-element differential-expression results, including log2 fold changes and adjusted P values.

Files under `supplementary_table_source/` retain their earlier source-sheet numbering for provenance and are **not** the final journal numbering after the August deletion/renumbering pass. They are preserved to avoid silently rewriting historical source filenames. Use `FINAL_NUMBERING_MAP.tsv` and `REPRODUCIBILITY_MATRIX.tsv` for the final display-item mapping.

## Repository contents

This repository contains:

- portable R scripts and configuration;
- historical deposited supplementary source sheets retained for provenance;
- compact source-data tables for selected figures/analyses;
- representative JASPAR PFMs used by the Fig. 7D workflow;
- script-to-display-item mappings for the final figures and tables;
- input provenance, final-numbering map and metadata validation code.

Large sequencing/reference inputs and journal submission artwork are not duplicated.

## Code repository

https://github.com/BaibhavBarbaruah/ONSEN-Heat-Shock-Element
