# Changelog

## 1.1.0 — 2026-08-16

- synchronized public-facing repository metadata with the final August 2026 Biology Open revision;
- adopted the physical-forward exact-coordinate HSF motif-placement metric as the principal non-redundant count and retained transitive overlap merging only as a stringent sensitivity analysis;
- recorded the corrected Col-0 ONSEN terminal-window placement values and the revised 1,930-region broad non-ONSEN TE background;
- documented the direct ONSEN versus non-ONSEN LTR-retrotransposon HSF and basal methylation comparisons;
- documented the corrected complete-TAIR10 RNA-seq workflow and genome-wide gene/TE differential-expression analyses;
- removed the AP2/ERF-only supplementary table from the final journal package;
- changed the final supplementary package to Fig. S1-Fig. S4 and Tables S1-S14, with the new global gene/TE differential-expression workbook as Table S14;
- updated `FINAL_NUMBERING_MAP.tsv` and `REPRODUCIBILITY_MATRIX.tsv` to the final journal-facing numbering;
- explicitly marked the files in `supplementary_table_source/` as legacy pre-final-numbering provenance rather than silently renaming historical source files.

## 1.0.1 — 2026-07-17

- synchronized `Table_S1.xlsx` through `Table_S14.xlsx` with the then-current journal workbooks;
- synchronized the deposited TSV source sheets with the corresponding Excel worksheets, including the restored `S9D_selected_outliers` sheet and corrected Table S10 class-level fold-change summary;
- replaced residual mutation terminology in Tables S1-S3 with the manuscript's in-silico-design terminology;
- standardized Table S6 and Table S8 titles, identifiers, symbols and worksheet-source wording;
- improved Table S4 and Table S8 title-row formatting without changing numerical content;
- refreshed all affected SHA-256 records and revalidated the package manifest.

## 1.0.0 — 2026-07-15

- synchronized supplementary references to the July revision numbering;
- deposited journal-numbered TSV source sheets and explicit figure-workflow mappings used at that stage;
- added the four-threshold/continuous HSF-score workflow and source data;
- added non-redundant HSF, constrained-mutant and exact-GC mappings;
- added Fig. 7C sequence-identity source data;
- added exact Fig. 7D representative motif metadata and PFMs (HSFC1 MA1667.2; DOF1.8 MA0981.2);
- updated input provenance, reproducibility matrix, package manifest and validation checks;
- removed reliance on machine-specific drive paths in public-facing scripts.
