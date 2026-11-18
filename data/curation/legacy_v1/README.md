# Curation of the previous extraction workbook (master_analysis_input.xlsx)

These files document the source check of the previous extraction (383 records) carried out on 17 November 2026,
before the re-extraction in `data/raw/benchmark.xlsx`:

- `observation_corrections.csv`: field-level corrections found by an AI-assisted pre-check against the source
  publications (tier `applied` = clear transcription errors with a source value; tier `review` = for author decision).
  The findings were passed to the authors for the re-extraction; they are applied only when the legacy workbook is
  loaded (`scripts/real_simulated_sensitivity.py::apply_observation_corrections`).
- `bio_dataset_map.csv`: biological-dataset identifiers for the legacy dataset IDs, verified against the sources.

The current analyses use `data/processed/analysis_input.xlsx`, built from `data/raw/benchmark.xlsx` by
`scripts/build_analysis_input.py`.
