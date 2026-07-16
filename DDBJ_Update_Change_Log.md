# Biology Open DDBJ update — exact change log

## Manuscript

- Replaced the RNA-seq deposition text in **Data and resource availability** with the author-supplied wording exactly, including the term **“accession numbers”**, BioProjects PRJDB39904 and PRJDB42759, and BioSamples SAMD01789795–SAMD01789797 and SAMD01943917–SAMD01943919.
- Revised the Introduction's data-source sentence to distinguish public DNA methylation data from the RNA-seq data generated in this study.
- Revised the Discussion's dataset-comparison sentence to identify the study-generated samples accurately as non-stressed and heat-treated whole seedlings rather than a previously published shoot dataset.
- Removed the obsolete pre-update repository commit pin so that the manuscript does not identify an earlier repository state as the audited version.
- Clarified the Acknowledgements so that public genomic and methylation datasets and published expression evidence are distinguished from the six laboratory-generated RNA-seq libraries.
- Changed the Figure 5 legend title from “Reanalysis of heat-associated RNA-seq signal…” to “Heat-associated RNA-seq signal…”; its numerical results and panel descriptions were not changed.
- Corrected the formatting of *Arabidopsis thaliana* in the reference-genome Methods sentence; the wording was not changed.
- Corrected the document Author and Last Modified By metadata to “Baibhav R. Barbaruah.”
- Preserved the statement that the experiment contained three biological replicates per condition and all six libraries were analysed.
- Did not alter the plant-growth, heat-treatment, RNA-extraction, HaploX library-preparation, NovaSeq-sequencing or RNA-seq-analysis Methods.
- Did not alter, reformat or relocate the approved AI-use disclosure. Its final paragraph XML is byte-identical to the newest uploaded manuscript.

## Response to reviewers

- Updated the Handling Editor response on public scripts and sequencing data with laboratory provenance, three biological replicates per condition, and the exact author-supplied DDBJ wording.
- Updated Reviewer 1, Comment 9 with the same laboratory provenance, replicate count and exact DDBJ wording.
- Updated the overview of major revisions to describe the six laboratory-generated libraries and released accession numbers.
- Removed “archived dataset” and unsupported original-study randomisation/blinding/inclusion-reporting language from the overview and journal-compliance response.
- Removed the obsolete pre-update repository commit pin.
- Retained the clarification that the public methylome and study-generated whole-seedling RNA-seq libraries are not sample- or tissue-matched.
- Kept editor/reviewer comment headings with the following response so headings are not stranded or clipped at page boundaries.
- Corrected the document Author and Last Modified By metadata to “Baibhav R. Barbaruah.”

## Official journal checklist

- Added the plant-stock source as “Arabidopsis Col-0: ABRC stock CS60000.”
- Updated the data/resource field to list BioProjects PRJDB39904 and PRJDB42759, their six BioSamples, and the retained GEO methylome accession; no DRA/DRR/DRX or pending identifiers remain.

## Repository package

- `README.md`: added laboratory provenance, three biological replicates per condition and the exact DDBJ wording.
- `DATA_AVAILABILITY.md`: replaced the incorrect public-reanalysis/DRA013053 provenance and DRR/DRX table with laboratory provenance, replicate count and the exact DDBJ wording.
- `RNAseq_sample_metadata_template.csv`: replaced DRA/run/experiment/source-publication fields with BioSample, BioProject and laboratory-provenance fields for all six libraries.
- `INPUT_PROVENANCE.tsv`: replaced public-reanalysis provenance with study-generated metadata, alignments and count-file provenance; retained the distinction between the independent published cross-accession evidence and the six study libraries.
- `06_rnaseq_analysis.R`: replaced validation against DRA013053 DRR/DRX assignments with validation of the six BioSample/BioProject assignments, Col-0 identity, three replicates per condition and laboratory provenance.
- `10_validate_manuscript_outputs.R`: replaced DRA/DRR/DRX checks with BioSample, BioProject, replicate and laboratory-provenance checks, and moved the validation-report/failure gate after those checks so they are actually included in repository validation.
- `FILE_CHECKSUMS_SHA256.tsv`: updated the `INPUT_PROVENANCE.tsv` SHA-256 value and byte count.

## Removed everywhere relevant

- DRA013053 and the associated DRR/DRX identifiers.
- Other DRR, DRX and DRA identifiers for these six libraries.
- `[HS ACCESSIONS PENDING]` and statements that DRA, experiment or run accessions remain pending.
- Statements presenting the six libraries as previously published public RNA-seq or as a reanalysis of Nozawa et al. data.
- The obsolete manuscript/response claim that the pre-update repository commit was the audited repository state.

## Final verification performed

- Rendered the final manuscript (28 pages) and response letter (5 pages) and visually inspected the revised pages; no clipping, displacement or corruption was detected.
- Rendered and visually inspected all three pages of the revised journal checklist.
- Confirmed that the exact DDBJ paragraph occurs once in the manuscript and twice in the response letter.
- Confirmed that the locked AI-use paragraph is XML-identical to the newest uploaded manuscript.
- Confirmed that the experimental RNA-seq Methods are text-identical to the newest uploaded manuscript.
- Confirmed that the final Word files contain no comments, tracked changes, hidden text or highlighting; the manuscript's only field is its balanced, rendered page-number field.
- Scanned the final Word files and repository working tree for obsolete DRA/DRR/DRX identifiers, pending-status phrases, previous-public-RNA-seq provenance and the obsolete repository commit reference; none remain.
- Verified the six-row metadata map, both BioProjects, all six BioSamples, replicate numbering and laboratory provenance programmatically.
- Verified `INPUT_PROVENANCE.tsv` SHA-256 (`62ab75b80b15550e49b533de4c3d89cf7d97ec47671f954f7599475eea247f3b`) and size (7,785 bytes).
- Ran `git diff --check`; no whitespace errors were found.
- Cross-checked the seven main-figure PDFs, the supplementary-figure PDF and all fourteen supplementary workbooks against the manuscript citations, legends and final numbering. Figure 7 has A/B on the top row, C at bottom left and D at bottom right with two vertically stacked motif logos.
- R was not installed in the available workspace, so the R validation script and full analysis pipeline were not executed.

## User action

- Commit/push or merge the updated repository package into the public GitHub repository. The remote repository was not changed from this workspace because the required authenticated GitHub CLI session was unavailable.
- Supply the highlighted manuscript required by the editor and confirm that it is text-identical to the returned clean manuscript apart from revision highlighting. No highlighted manuscript was uploaded in this session, although the response letter says one has been prepared.
- Confirm the truthful reporting of sample-size determination, exclusions, randomisation and blinding before submission. The journal checklist's overall methodology/statistics confirmation is checked, but the manuscript does not currently contain explicit statements for each of these items; no scientific assertion was added without author confirmation.
- Run `10_validate_manuscript_outputs.R` in an R environment with the repository dependencies before publishing the repository update.
