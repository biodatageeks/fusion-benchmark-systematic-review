# Suggested supplementary/additional files

This manifest lists the repository files that can be referenced as supplementary material or data/code availability items for the manuscript.

## Primary supplementary data

1. `data/raw/master_analysis_input.xlsx`
   - Minimal master workbook used for all manuscript calculations.
   - Contains observation-level performance data, benchmark metadata, dataset metadata, tool-name harmonization, and Edgren truth-set annotations.

2. `data/raw/edgren_rnafusion.xlsx`
   - RNAfusion-derived Arriba, STAR-Fusion, and FusionCatcher calls on the Edgren dataset, together with the 99-fusion truth set.

3. `data/raw/manuscript_publication_selection.xlsx`
   - Consolidated publication selection register for the author-corrected 2015-2025 window.
   - Includes source records and abstracts, search strategies, report-level evidence, staged decisions, and PRISMA reconciliation.
   - AI-assisted proposals remain pending human confirmation; this is not a completed historical screening log or final inclusion list.

## Main derived tables and figures

- `results/real_simulated/`
- `results/read_length/`
- `results/benchmark_statistics/`
- `results/benchmark_reporting/`
- `results/precision_recall_tradeoff/`
- `results/new_tool_performance/`
- `results/edgren_rnafusion_validation/`
- `results/figures/`
- `results/tables/`

## Reproducibility

Run all Python analyses and regenerate derived tables and figures with:

```bash
./run_all.sh
```
