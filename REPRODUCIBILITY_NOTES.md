# Reproducibility and interpretation notes

## Final numbering and deposited tables

The journal submission uses `Table_S1.xlsx`-`Table_S13.xlsx`, Fig. 1-Fig. 7 and Fig. S1-Fig. S4. `FINAL_NUMBERING_MAP.tsv` documents intermediate revision numbering and the final mapping.

The root Excel files are the final audited workbooks. `supplementary_table_source/` contains an exact tab-separated export of every worksheet. `08_write_supplementary_tables.R` can rebuild content-equivalent workbooks from those exports; the tracked Excel files remain authoritative for final formatting and cell types.

## Principal HSF redundancy rule

The principal metric maps all predictions to physical forward-sequence coordinates and collapses predictions with identical 1-based start/end coordinates across HSF PWM models and strands. Predictions at different coordinates remain separate. The strict overlap-merged metric transitively merges HSF intervals that overlap by at least one nucleotide and is reported as a stringent sensitivity analysis.

Both metrics describe computational sequence compatibility, not independent binding events, occupancy or regulatory function in vivo.

## Designed-sequence robustness

Table S3 reports HSF-focused robustness controls for the designed 49-bp sequence. The constrained library contains 5,000 sequences. The exact-GC design space contains the selected design plus 5,119 alternatives; none retained a non-redundant HSF motif-coordinate placement at the primary threshold.

## Final strict TE-background universe

The archived coordinate scan contained 1,942 strict non-ONSEN TE regions. Final annotation and ONSEN/ATCOPIA78 harmonisation removed 12 regions, giving the 1,930-region universe used for reported threshold and continuous-score analyses in Figs. 3/S2 and Tables S6-S7. Table S7H retains the 1,942 count only as archived scan provenance.

## Methylation interpretation

Table S8A is the principal direct comparison between ONSEN terminal candidate regions and annotation-defined non-ONSEN LTR-retrotransposon terminal regions. Table S8B is a window-level sensitivity analysis, and Table S8C is the broader ordinary-TE contextual comparison. The source output did not retain a CG window-level Cliff's delta, so Table S8B reports it as `NA`/`not reported` rather than inferring it.

The methylome is from unstressed Col-0 leaves and is used as basal chromatin context. It is not a heat-treatment methylome and is not sample- or tissue-matched to the whole-seedling RNA-seq experiment.

## RNA-seq interpretation

Table S9 contains the full genome-wide results. DESeq2 produced model estimates for 25,912 genes, raw P values for 25,911 and non-missing BH-adjusted P values for 23,922 after independent filtering. The TE workflow retains 31,189 annotated elements and statistically tests 2,101 after expression filtering.

Table S10 contains candidate-window signal: sixteen curated ONSEN terminal windows, two additional ATCOPIA78 candidate regions and the fixed eight-region HSF-rich non-ONSEN TE subset. Fig. 5C uses `log2(mean HS CPM + 1)` on the x-axis and `max(0, log2[(mean HS CPM + 0.05)/(mean NS CPM + 0.05)])` for point size. These measurements are not interpreted as definitive copy-specific expression.

## Natural-accession ascertainment

Accession candidates were recovered with a Col-0-derived 49-bp seed allowing up to four mismatches. Similarity and motif-retention statements therefore apply to the detected candidate set and do not establish conservation across all ONSEN-related sequences.

Paired HSE/LTR values are structural proxies, not definitive copy-number estimates. In the corrected Fig. S4C/Table S11 summary, Ler has six paired proxies and three unpaired candidate windows.

## Figure 7 logo models

The Fig. 7D representative models are HSFC1 (JASPAR **MA1667.2**) and DOF1.8 (JASPAR **MA0981.2**). Their PFMs and metadata are deposited under `motifs/` and `source_data/Figure7_logo_model_metadata.csv`. The logos illustrate PWM preferences and do not demonstrate binding or regulatory activity in vivo.

## Figure assembly

The analysis scripts generate panel source files and tabular results. Journal-ready assembled figure PDFs are submission assets and are not duplicated in this code/data repository. `REPRODUCIBILITY_MATRIX.tsv` maps each final display item to its responsible workflow and inputs.
