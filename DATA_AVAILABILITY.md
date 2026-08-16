# Data and code availability

## RNA-seq generated in this study

All six RNA-seq libraries were generated in our laboratory, with three biological replicates per condition.

- Col-0 non-stressed: BioProject **PRJDB39904**; BioSamples **SAMD01789795, SAMD01789796 and SAMD01789797**.
- Col-0 37°C, 24-h heat stress: BioProject **PRJDB42759**; BioSamples **SAMD01943917, SAMD01943918 and SAMD01943919**.

FASTQ and BAM files are not duplicated in this repository. The exact sample-to-accession map and expected local filenames are recorded in `RNAseq_sample_metadata_template.csv` and `INPUT_PROVENANCE.tsv`.

Gene-level differential expression uses DESeq2. Global transposable-element differential expression uses fractional multimapping-aware TE counts, expression filtering and limma-voom. Candidate-window measurements are fractional multimapping-aware signal and are not interpreted as definitive copy-specific ONSEN expression.

## Public methylome

- GEO series: **GSE43857**
- sample: **GSM1085222**
- material: unstressed *Arabidopsis thaliana* Col-0 leaf
- expected local analysis filename: `GSM1085222_mC_calls_Col_0.tsv.gz`

The methylome is used to describe basal DNA-methylation context. It is not a heat-treated methylome and is not paired to the whole-seedling RNA-seq experiment.

## Reference resources

The workflows use the TAIR10 reference genome, the gene/TE annotation resources listed in `INPUT_PROVENANCE.tsv`, JASPAR CORE Plants motif data and chromosome-level assemblies for the natural-accession survey. Exact expected filenames and analysis roles are recorded in `INPUT_PROVENANCE.tsv`.

## Supplementary source data

The final supplementary set comprises **Fig. S1-Fig. S4** and **Table S1-Table S13**. Formatted Excel workbooks are at the repository root. Data-equivalent plain-text worksheet mirrors are under `supplementary_table_source/`.

Table S9 contains the complete genome-wide gene and transposable-element differential-expression results. Table S10 contains candidate-window values, replicate-level CPM and the statistical analysis of replicate-summed fractional CPM across the sixteen curated ONSEN terminal candidate windows.

## Code repository

https://github.com/BaibhavBarbaruah/ONSEN-Heat-Shock-Element
