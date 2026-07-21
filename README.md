# RNA-seq fusion detection benchmark meta-analysis

This repository contains curated extraction tables, analysis scripts, derived tables, and figures for a systematic comparative meta-analysis of published benchmarks of RNA-seq-based gene fusion detection tools.

The analysis evaluates how reported fusion detection performance depends on benchmark design, including dataset origin, truth-set definition, benchmark purpose, read length, and reproducibility-related reporting.

## Repository structure

```text
data/
  raw/          Manually curated source spreadsheets used as analysis inputs.
  processed/    CSV exports and enriched tables derived from raw spreadsheets.
scripts/        Python scripts for main analyses and manuscript figures.
scripts/meta_analysis_R/
                Additional exploratory random-effects summaries and forest plots.
results/
  real_simulated/              Mixed-effects real-vs-simulated analyses.
  benchmark_statistics/        Benchmark landscape heatmap and source tables.
  benchmark_reporting/         Reporting/truth-set heatmap and source tables.
  precision_recall_tradeoff/   Tool-level recall-precision summaries.
  new_tool_performance/        New-tool benchmark context summaries.
  figures/                     Exploratory pooled/forest/drift figures.
  tables/                      Exploratory pooled/drift/sensitivity tables.
docs/           Placeholder for manuscript-related notes.
```

## Main input files

- `data/raw/real_simulated.xlsx` - tool-dataset performance table used for real-versus-simulated, leave-one-benchmark-out, read-length, benchmark-purpose, and Edgren sensitivity analyses.
- `data/raw/all_data.xlsx` - broader tool-level performance table used for recall-precision trade-off analyses.
- `data/raw/heatmap.xlsx` - benchmark statistics and dataset-feature table used for the benchmark landscape heatmap.
- `data/raw/reproducibility_and_gold_standard.xlsx` - truth-set, reproducibility, reporting, and benchmark-purpose annotations.
- `data/raw/Benchmarki_Agnieszka.xlsx` - source workbook for the exploratory random-effects workflow.
- `data/raw/edgren_rnafusion.xlsx` - RNAfusion-derived calls from Arriba, STAR-Fusion, and FusionCatcher on the Edgren dataset, together with the 99-fusion truth set.

## Reproducing the Python analyses

Create and activate a Python environment, then install dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Run all main Python analyses:

```bash
./run_all.sh
```

If you want to run individual scripts:

```bash
python scripts/real_simulated_sensitivity.py
python scripts/plot_real_simulated_metrics.py
python scripts/plot_benchmark_statistics_heatmap.py
python scripts/plot_benchmark_reporting_heatmap.py
python scripts/precision_recall_tradeoff.py
python scripts/new_tool_benchmark_performance.py
python scripts/edgren_rnafusion_validation.py
```

## Optional exploratory random-effects workflow

The exploratory random-effects summaries and forest plots require R and the `metafor` package:

```r
install.packages("metafor")
```

Then run:

```bash
python scripts/meta_analysis_R/01_export_from_excel.py
python scripts/meta_analysis_R/02_data_quality_report.py
Rscript scripts/meta_analysis_R/03_meta_analysis.R
python scripts/meta_analysis_R/05_edgren_drift_figure.py
python scripts/meta_analysis_R/04_make_figures.py
```

These analyses are descriptive and are used to quantify between-benchmark heterogeneity. They should not be interpreted as definitive universal rankings of fusion detection tools.

## Edgren dataset annotation

Edgren-derived datasets were annotated as: `2.3`, `5.1`, `6.1`, `6.2`, `7.6`, `7.7`, `8.3`, and `10.2`. These annotations are used for the Edgren truth-set heterogeneity case study and for sensitivity analyses excluding Edgren-derived observations.

## Key outputs used in the manuscript

- `results/benchmark_statistics/benchmark_statistics_heatmap.pdf`
- `results/benchmark_reporting/benchmark_reporting_truthset_heatmap.pdf`
- `results/real_simulated/real_vs_simulated_f1_precision_recall.pdf`
- `results/real_simulated/dataset_type_models.csv`
- `results/real_simulated/dataset_type_models_without_edgren.csv`
- `results/real_simulated/leave_one_benchmark_out_f1.csv`
- `results/precision_recall_tradeoff/tool_precision_recall_summary_real_without_edgren.csv`
- `results/new_tool_performance/new_tool_benchmark_summary.csv`
- `results/edgren_rnafusion_validation/edgren_rnafusion_metrics_summary.csv`
- `results/edgren_rnafusion_validation/edgren_rnafusion_metrics.pdf`
- `results/figures/edgren_drift_scatter.pdf`
- `results/tables/pooled_estimates.csv`

## Data provenance

All values were manually extracted or curated from published benchmark articles and associated supplementary materials. The raw spreadsheets in `data/raw/` are the primary curated data sources. CSV files in `data/processed/` and `results/` are derived outputs.

## Notes for publication

Before public release, add:

- final citation metadata for the associated manuscript;
- a permanent repository archive DOI, e.g. Zenodo;
- a license selected by the authors and institution;
- any final corrections to raw extraction tables made during manuscript revision.
