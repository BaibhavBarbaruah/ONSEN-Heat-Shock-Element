# Data and code availability

## RNA-seq generated in this study

All six RNA-seq libraries were generated in our laboratory, with three biological replicates per condition.

- Col-0 non-stressed: BioProject **PRJDB39904**; BioSamples **SAMD01789795, SAMD01789796 and SAMD01789797**.
- Col-0 37°C, 24-h heat stress: BioProject **PRJDB42759**; BioSamples **SAMD01943917, SAMD01943918 and SAMD01943919**.

The repository does not duplicate FASTQ or BAM files. The exact sample-to-accession map and expected local filenames are recorded in `RNAseq_sample_metadata_template.csv` and `INPUT_PROVENANCE.tsv`.

The final RNA-seq analyses use the complete TAIR10 counting universe. Candidate-window measurements are fractional, multimapping-aware **candidate-window signal**, not definitive copy-specific ONSEN expression. Genome-wide gene differential expression uses DESeq2; the genome-wide TE workflow uses fractional multimapping-aware counts, expression filtering and limma-voom.

## Public methylome

- GEO series: **GSE43857**
- sample: **GSM1085222**
- material: unstressed *Arabidopsis thaliana* Col-0 leaf
- expected local analysis filename: `GSM1085222_mC_calls_Col_0.tsv.gz`

This dataset is used only to describe basal methylation context. It is not a heat-treated methylome and is not interpreted as heat-induced methylation change or as a sample-matched methylation-expression dataset.

## Reference resources

- *Arabidopsis thaliana* TAIR10 reference genome;
- Araport11 gene and transposable-element annotations;
- JASPAR CORE Plants motif collections used by the deposited workflows;
- chromosome-level assemblies used for the accession survey.

Expected filenames and analysis roles are listed in `INPUT_PROVENANCE.tsv`.

## Final supplementary package

The journal-facing package contains **Fig. S1-Fig. S4** and **Table S1-Table S13**. There is no final Table S14. The earlier AP2/ERF-only and published-accession-evidence workbooks are not part of the final package.

Final Table S9 contains the editor-requested genome-wide gene and transposable-element differential-expression results, including all 32,833 genes and all 31,189 TEs. Final Table S10 contains candidate-window values and an explicit statistics sheet for the Fig. 5B replicate-summed comparison.

The authoritative final Excel workbooks are deposited at the repository root. Plain-text copies of all 36 worksheets are under `supplementary_table_source/` using the same final numbering.

## Code repository

https://github.com/BaibhavBarbaruah/ONSEN-Heat-Shock-Element
