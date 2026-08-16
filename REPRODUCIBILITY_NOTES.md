# Reproducibility and interpretation notes

## Final numbering and authoritative files

The final journal package contains Tables S1-S13 and Figs S1-S4. `FINAL_NUMBERING_MAP.tsv` records the deletion/renumbering history. There is no final Table S14.

The root `Table_S1.xlsx`-`Table_S13.xlsx` workbooks are the authoritative submission tables. `supplementary_table_source/` contains final-numbered TSV mirrors of all 36 worksheets, including the complete global gene and TE results in Table S9. No pre-final AP2/ERF or S14 source sheets are retained.

## Exact-coordinate and overlap-merged HSF metrics

The principal HSF count is the physical-forward exact-coordinate metric: predictions sharing the same 1-based start/end coordinates are collapsed across PWM models and strands. Predictions at different coordinates remain separate.

The strict overlap-merged metric transitively merges any intervals overlapping by at least 1 nt. It is a deliberately stringent sensitivity metric and is reported separately in Fig. S1 and Table S5.

## Threshold and continuous-score analysis

`03B_threshold_and_continuous_sensitivity.R` evaluates relative-score thresholds 0.80, 0.85, 0.90 and 0.95 and three threshold-independent regional metrics. The final summary, region-level and QC outputs are all deposited in Table S7 and its TSV mirrors.

## Direct LTR comparisons

Fig. 4A/Table S8A use ONSEN copy-level units and annotation-defined non-ONSEN LTR-retrotransposon element-level units. Left and right termini are aggregated before testing so termini from the same element are not treated as independent observations. The window-level analysis in Table S8B is a sensitivity analysis.

Fig. 4C compares eight ONSEN elements with 794 annotation-defined non-ONSEN LTR-retrotransposon elements using non-redundant HSF motif-coordinate placement density. This is distinct from the broader 1,930-region strict TE background used in Fig. 3 and Tables S6-S7.

Compact final summaries for both direct comparisons are under `source_data/`. The analysis scripts recompute exact statistics when the processed element-level inputs listed in `INPUT_PROVENANCE.tsv` are available.

## Candidate-window RNA-seq statistics

Table S10C contains fractional CPM for every candidate window and biological replicate. Table S10D sums the 16 curated ONSEN terminal windows within each replicate and reports the two-sided Welch's t-test used in Fig. 5B together with the corresponding Wilcoxon rank-sum result. Candidate-window values are multimapping-aware signal and are not interpreted as definitive copy-specific ONSEN expression.

## Figure 7 logo models

The representative Fig. 7D models are HSFC1 (JASPAR **MA1667.2**) and DOF1.8 (JASPAR **MA0981.2**). Their PFMs and metadata are deposited under `motifs/` and `source_data/`. These logos describe PWM preferences and do not demonstrate occupancy or regulatory function in vivo.

## Methylation and RNA-seq limitations

The methylome is from unstressed Col-0 leaves and is used as basal chromatin context. The RNA-seq experiment used whole seedlings and is not paired to that methylome. No heat-induced methylation change is inferred.

## Accession ascertainment

Accession candidates were recovered with a Col-0-derived 49-bp seed allowing up to four mismatches. Similarity and motif-retention statements apply to the detected candidate set and do not establish conservation across every ONSEN-related sequence in each genome.
