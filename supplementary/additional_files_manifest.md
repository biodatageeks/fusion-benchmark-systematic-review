# Suggested supplementary/additional files

This manifest lists the repository files that can be referenced as supplementary material or data/code availability items for the manuscript.

## Primary supplementary data

1. `data/raw/master_analysis_input.xlsx`
   - Minimal master workbook used for all manuscript calculations.
   - Contains observation-level performance data, benchmark metadata, dataset metadata, tool-name harmonization, and Edgren truth-set annotations.

2. `data/raw/edgren_rnafusion.xlsx`
   - RNAfusion-derived Arriba, STAR-Fusion, and FusionCatcher calls on the Edgren dataset, together with the 99-fusion truth set.

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
