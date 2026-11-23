# Second-author check of extracted records

`second_author_check.csv` (semicolon-separated) lists the 31 extracted records that a second author verified against the source publications and their supplementary data: 26 records sampled at random from 13 benchmarks (seed 20261121), complemented by one random record from each of the remaining five benchmarks.

For each record the file gives the extracted values (`TP`, `FP`, `truth_total`, `recall`, `precision`, `f1_score`), the values found in the source (`source_*` columns, decimal commas as entered) and where they were found (`source_location`). `match_overall` (Y/N), `discrepancy_type`, `action_needed` and `verifier_notes` record the outcome.

One record was discrepant: the Uhrig et al. 2021 spike-in data (record V08) define 10 rather than 20 reference fusions per sample. The seven observations of this dataset were corrected in `data/raw/master_analysis_input.xlsx`.
