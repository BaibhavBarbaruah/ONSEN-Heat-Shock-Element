# Final supplementary-table source sheets

This directory contains one tab-separated export for every worksheet in the final journal-facing `Table_S1.xlsx`-`Table_S13.xlsx` workbooks.

Filenames follow:

`Table_S<table number>__<worksheet name>.tsv`

Each TSV preserves the worksheet's displayed rows, including the title and note rows. The tracked Excel workbooks at the repository root are authoritative for final formatting, merged cells, numeric cell types and column widths. `08_write_supplementary_tables.R` rebuilds content-equivalent workbooks from these source sheets into `ONSEN_OUTPUT_ROOT/rebuilt_final_tables/`.

These are final-numbered sources, not the legacy source numbering used in earlier revision rounds. `../FINAL_NUMBERING_MAP.tsv` documents the historical-to-final mapping.
