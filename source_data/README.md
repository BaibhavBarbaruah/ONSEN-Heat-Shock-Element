# Compact source data

This directory contains compact, final-numbered source tables that can be inspected without the large external genome, methylome or sequencing inputs.

- `S8A_threshold_stats.tsv`: four-threshold HSF-density statistics.
- `S8B_continuous_stats.tsv`: threshold-independent regional-score comparisons.
- `S8C_region_metrics.tsv.gz`: generated region-level threshold and continuous-score metrics (written to `ONSEN_OUTPUT_ROOT` by the Reviewer 1.3 workflow; not duplicated here).
- `S8D_QC.tsv`: analysis quality-control and HSF motif inventory.
- `Figure7_variant_metrics.tsv`: mismatch counts and 49-bp sequence identity used for Fig. 7C.
- `Figure7_logo_model_metadata.csv`: exact representative JASPAR models used for Fig. 7D.
- `Table_S13_natural_variant_TF_family.tsv`: compact source table for the natural-variant TF-family comparison.
- `Table_S14_published_accession_evidence.tsv`: cited literature synthesis underlying Table S14.

Large reference inputs and processed files are documented in `../INPUT_PROVENANCE.tsv`.
