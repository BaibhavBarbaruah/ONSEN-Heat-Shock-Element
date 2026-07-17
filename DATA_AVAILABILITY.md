# Data and code availability

## RNA-seq generated in this study

All six RNA-seq libraries were generated in our laboratory, with three biological replicates per condition.

The RNA-seq data generated in this study have been deposited in DDBJ. The three Col-0 non-stressed libraries are available under BioProject PRJDB39904, with accession numbers SAMD01789795, SAMD01789796 and SAMD01789797. The three Col-0 37°C, 24-h heat-stress libraries are available under BioProject PRJDB42759, with accession numbers SAMD01943917, SAMD01943918 and SAMD01943919.

The repository does not duplicate FASTQ/BAM files. The exact sample-to-accession map and workflow roles are recorded in `RNAseq_sample_metadata_template.csv` and `INPUT_PROVENANCE.tsv`.

## Public methylome

- GEO series: **GSE43857**
- sample: **GSM1085222**
- material: unstressed *Arabidopsis thaliana* Col-0 leaf
- expected analysis filename: `GSM1085222_mC_calls_Col_0.tsv.gz`

This dataset is used only to describe basal methylation context. It is not a heat-treated methylome and is not interpreted as a measurement of heat-induced methylation change.

## Reference resources

- *Arabidopsis thaliana* TAIR10 reference genome
- Araport11 gene and transposable-element annotation
- JASPAR 2024 and JASPAR 2026 CORE Plants motif collections
- chromosome-level assemblies used for the accession survey

Exact expected filenames and analysis roles are listed in `INPUT_PROVENANCE.tsv`.

## Natural-accession evidence

The genome scan used a Col-0-derived 49-bp seed allowing up to four mismatches. More divergent ONSEN-related sequences may therefore be absent from the detected candidate set.

Published cross-accession heat-response evidence discussed in the manuscript includes ENA project **PRJEB64476**. It is cited as independent published comparative evidence and is not presented as accession-specific RNA-seq generated in this study. The cited synthesis underlying Table S14 is deposited as `source_data/Table_S14_published_accession_evidence.tsv` and as the journal-numbered TSV source sheet.

## Deposited revision outputs

This repository contains:

- portable R scripts and configuration;
- journal-numbered TSV source sheets used to reconstruct Table S1-Table S14 workbooks;
- compact source-data TSV files for Fig. 7 and Table S8/Fig. S3;
- exact representative JASPAR PFMs used by the Fig. 7D workflow;
- explicit script-to-display-item mappings for all main and supplementary figures;
- input provenance, final-numbering map, checksums and validation code.

Large public/reference inputs and journal submission artwork are not duplicated. The deposited source sheets, compact source data, motif matrices, checksums and reconstruction scripts preserve the exact numerical revision outputs while keeping the repository portable.

## Code repository

https://github.com/BaibhavBarbaruah/ONSEN-Heat-Shock-Element
