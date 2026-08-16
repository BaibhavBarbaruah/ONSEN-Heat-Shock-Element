# Data and code availability

## RNA-seq generated in this study

All six RNA-seq libraries were generated in the authors' laboratory, with three biological replicates per condition.

The three Col-0 non-stressed libraries are available under BioProject **PRJDB39904**, with BioSample accessions **SAMD01789795, SAMD01789796 and SAMD01789797**. The three Col-0 37°C, 24-h heat-stress libraries are available under BioProject **PRJDB42759**, with BioSample accessions **SAMD01943917, SAMD01943918 and SAMD01943919**.

The repository does not duplicate FASTQ or BAM files. The sample-to-accession map and workflow roles are recorded in `RNAseq_sample_metadata_template.csv` and `INPUT_PROVENANCE.tsv`.

The final RNA-seq analyses use the corrected complete-TAIR10 alignment/counting universe. Candidate-window measurements are fractional, multimapping-aware candidate-window signal and are not treated as copy-specific ONSEN expression. Genome-wide gene differential expression was analysed with DESeq2; genome-wide TE analysis used fractional counts, expression filtering and limma-voom.

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

## Natural-accession scope

The genome scan used a Col-0-derived 49-bp seed allowing up to four mismatches. More divergent ONSEN-related sequences may therefore be absent from the detected candidate set. Published cross-accession studies are used as biological context in the manuscript, not as functional validation of the present motif predictions and not as a separate final supplementary table.

## Final supplementary package

The journal-facing package contains **Fig. S1-Fig. S4** and **Table S1-Table S13**. The editor-requested genome-wide gene/transposable-element differential-expression results are Table S9; candidate-window RNA-seq results are Table S10.

This repository deposits:

- portable R scripts and configuration;
- the exact final `Table_S1.xlsx`-`Table_S13.xlsx` workbooks;
- exact worksheet-level TSV exports under `supplementary_table_source/`;
- compact source data for Fig. 7;
- representative JASPAR PFMs used by the Fig. 7D workflow;
- final numbering, provenance and reproducibility metadata.

Large sequencing/reference inputs and journal-ready assembled figure artwork are not duplicated.

## Code repository

https://github.com/BaibhavBarbaruah/ONSEN-Heat-Shock-Element
