# Supplementary-table source sheets

This directory contains tab-separated source sheets deposited during the earlier revision workflow. The filenames preserve their **pre-final August 2026 numbering** for provenance.

They are **not** the final journal numbering after the August 2026 deletion/renumbering pass. In particular, the old AP2/ERF-only Table S3 was removed, subsequent tables shifted, the published accession-evidence synthesis became final Table S13, and a new genome-wide gene/TE differential-expression workbook became final Table S14.

Use `../FINAL_NUMBERING_MAP.tsv` as the authoritative mapping between these historical source numbers and the final journal-facing Tables S1-S14.

`08_write_supplementary_tables.R` now rebuilds these historical source sheets only into a dedicated `legacy_rebuilt_tables` output directory. It must not be used to generate the final submission workbooks.

The historical source sheets remain useful for provenance and for checking the earlier analysis stages. Corrected August 2026 values, analysis roles and final display-item assignments are summarized in the repository `README.md` and `REPRODUCIBILITY_MATRIX.tsv`.
