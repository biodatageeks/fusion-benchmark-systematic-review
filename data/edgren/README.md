# Edgren per-sample, multi-reference re-scoring (R1.2, R1.7, R3.3)

Inputs of `scripts/edgren_reference_sensitivity.py`:

- [`sample_manifest.csv`](sample_manifest.csv): the six tumour cell-line runs of SRP003186 with cell line, layout,
  read length and read pairs (confirmed against `data/raw/edgren_truthsets/manifest.csv` and
  `data/raw/edgren_predictions/README.md`). The normal breast run SRR064437 is not scored.
- [`reference_events.csv`](reference_events.csv): the published 27-, 53- and 99-fusion lists, one row per event and
  list, with cell line, source and evidence; gene aliases are harmonised across lists (`|` = alternatives). The script
  checks every event against the authors' per-run reference sets in `data/raw/edgren_truthsets/`
  (`results/revision/edgren_reference_sensitivity/reference_check.csv`); the only difference is a gene label
  (MED1--UPS32 in the 99-event source table, MED1--USP32 here).
- Caller predictions: `data/raw/edgren_predictions/<workflow version>/<SRA run>/<caller>.csv` for nf-core/rnafusion
  3.0.1 (original run) and 4.1.3 (documented re-run), Arriba, FusionCatcher and STAR-Fusion, unmodified apart from
  TSV to CSV conversion. They are parsed as fusion-report does in nf-core/rnafusion: one gene pair per output row;
  Arriba rows naming several neighbouring genes (intergenic breakpoints) give one gene pair per gene combination.
  Parsed this way, the union of the six 3.0.1 outputs reproduces `data/raw/edgren_rnafusion.xlsx` exactly
  (`raw_vs_workbook_check.csv`).

Run:

```bash
python scripts/edgren_reference_sensitivity.py            # per sample and per cell line, both runs (Figure 7, Table S11)
python scripts/edgren_reference_sensitivity.py --pooled   # pooled 3.0.1 predictions, cell lines ignored (submitted approach); Table S12
python scripts/edgren_reference_sensitivity.py --legacy   # pooled 99-fusion check on data/raw/edgren_rnafusion.xlsx
```

Outputs are written to `results/revision/edgren_reference_sensitivity/` (see the script docstring), the figure to
`results/revision/figures/edgren_rnafusion_reference_lists.pdf` and Supplementary Table S11 to
`results/revision/tables/supp_edgren_runs.tex`.

`sample_manifest_candidates.csv` keeps the evidence used to identify the runs (Uhrig 2021, Haas 2019, FusionCatcher
manual, ENA metadata). The `*_TEMPLATE.csv` files describe the expected columns.
