# Biology Open/DDBJ final-update change log

## Sequencing provenance

- Confirmed that all six RNA-seq libraries were generated in the authors' laboratory, with three biological replicates per condition.
- Non-stressed Col-0 libraries: BioProject PRJDB39904; BioSamples SAMD01789795-SAMD01789797.
- Col-0 37°C, 24-h heat-stress libraries: BioProject PRJDB42759; BioSamples SAMD01943917-SAMD01943919.
- Removed obsolete DRA/DRR/DRX and pending-accession language.

## Final submission synchronization

- Corrected the co-author spelling to **Rahmadani P. Airlangga**.
- Finalized supplementary numbering as Fig. S1-Fig. S4 and Table S1-Table S13.
- Placed the editor-requested genome-wide gene/TE differential-expression results in Table S9 and candidate-window RNA-seq results in Table S10.
- Reconciled the archived 1,942-region scan with the final 1,930-region strict non-ONSEN TE universe.
- Documented 25,912 gene model estimates, 25,911 raw gene P values and 23,922 non-missing BH-adjusted gene P values after DESeq2 independent filtering.
- Separated sixteen curated ONSEN windows, two additional ATCOPIA78 candidate regions and the fixed eight-region HSF-rich non-ONSEN TE subset in Table S10/Fig. 5C.
- Corrected Fig. S4C/Table S11 for Ler to six paired structural proxies and three unpaired candidate windows.
- Removed the AP2/ERF-only supplementary table and the separate published-accession-evidence table from the final package.

## Repository package

- Synchronized `Table_S1.xlsx`-`Table_S13.xlsx` with the final audited journal workbooks.
- Replaced legacy supplementary-table source sheets with exact final-numbered TSV exports for every worksheet.
- Updated README, data availability, provenance, numbering and reproducibility metadata.
- Updated reconstruction and validation scripts to expect the final thirteen-table package.
- Regenerated the file manifest and SHA-256 inventory.

## Verification

- Confirmed that clean and highlighted manuscripts are text-identical and contain no comments or tracked-change markup.
- Confirmed the 199-word abstract, numerical-order first citations for Tables S1-S13 and complete legends for seven main and four supplementary figures.
- Rendered and inspected all 28 pages of both manuscript versions, all response pages, every figure page and every supplementary-table worksheet.
- Confirmed that the seven main-figure artwork pages are pixel-identical to the supplied source pages and retain their original dimensions.
- Completed 322 automated cross-file content and consistency checks for the submission package.
