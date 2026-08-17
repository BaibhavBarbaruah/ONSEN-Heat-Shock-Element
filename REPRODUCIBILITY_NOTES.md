# Reproducibility and interpretation notes

## Final supplementary source data

The repository contains the formatted final `Table_S1.xlsx`-`Table_S13.xlsx` workbooks and data-equivalent TSV mirrors under `supplementary_table_source/`. The reproducibility matrix maps each final figure/table to its analysis script and principal inputs.

## Exact-coordinate and overlap-merged HSF metrics

The principal HSF count is the physical-forward exact-coordinate metric: predictions sharing the same 1-based start/end coordinates are collapsed across PWM models and strands. Predictions at different coordinates remain separate.

The strict overlap-merged metric transitively merges intervals overlapping by at least 1 nt. It is a deliberately stringent sensitivity metric and is reported separately in Fig. S1 and Table S5.

## Threshold and continuous-score analysis

`03B_threshold_and_continuous_sensitivity.R` evaluates relative-score thresholds 0.80, 0.85, 0.90 and 0.95 using exact-coordinate HSF placements. It also evaluates three threshold-independent regional metrics: maximum HSF relative PWM score, mean of the five highest HSF relative PWM scores and median of the ten per-model maximum scores. The final comparison contains sixteen ONSEN terminal windows and 1,930 strict non-ONSEN TE regions.

## Direct LTR comparisons

Fig. 4A/Table S8A use ONSEN copy-level units and annotation-defined non-ONSEN LTR-retrotransposon element-level units. Paired terminal regions are aggregated before testing so termini from the same element are not treated as independent observations. The window-level methylation analysis is retained as a sensitivity analysis.

Fig. 4C compares eight ONSEN elements with the annotation-defined non-ONSEN LTR-retrotransposon comparator using non-redundant HSF motif-coordinate placement density. This direct LTR comparison is distinct from the broader 1,930-region strict TE background used in Fig. 3 and Table S6-S7.

## RNA-seq

Gene-level differential expression uses DESeq2. The complete output retains all 32,833 annotated genes. Global TE analysis retains all 31,189 Araport11 `transposable_element` loci; fractional multimapping/overlap-aware TE counts are expression-filtered and 2,101 loci are tested with limma-voom.

For Fig. 5B, fractional CPM is summed across the sixteen curated ONSEN terminal candidate windows within each biological replicate. The principal comparison is a two-sided Welch's t-test on the three NS replicate sums and three HS replicate sums; Table S10 also reports the corresponding Wilcoxon rank-sum result. These measurements are candidate-window signal rather than definitive copy-specific ONSEN expression.

## Figure 7 PWM-model summaries and logo models

Fig. 7A,B and Table S13 summarize distinct JASPAR PWM model identities compatible with each 49-bp seed sequence. These family-level model-identity summaries are analytically distinct from the ten-model Arabidopsis HSF exact-coordinate analyses used for the primary HSF scans; the model-identity totals are not counts of independent motif-coordinate placements or binding sites.

The representative Fig. 7D models are HSFC1 (JASPAR **MA1667.2**) and DOF1.8 (JASPAR **MA0981.2**). Their PFMs and metadata are deposited under `motifs/` and `source_data/`. The logos illustrate PWM sequence preferences and do not demonstrate occupancy or regulatory function in vivo.

## Methylation and RNA-seq limitation

The methylome is from unstressed Col-0 leaves and is used as basal chromatin context. The RNA-seq experiment used whole seedlings and is not paired to that methylome. No heat-induced methylation change is inferred from the methylation dataset used here.

## Accession ascertainment

Accession candidates were recovered with a Col-0-derived 49-bp seed allowing up to four mismatches. Sequence-identity and motif-retention statements therefore apply to the detected candidate set and do not establish conservation across every ONSEN-related sequence in each accession.
