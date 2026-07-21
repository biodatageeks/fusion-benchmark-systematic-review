#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python3}"

$PYTHON_BIN scripts/real_simulated_sensitivity.py
$PYTHON_BIN scripts/plot_real_simulated_metrics.py
$PYTHON_BIN scripts/plot_edgren_truthset_heterogeneity.py
$PYTHON_BIN scripts/plot_benchmark_statistics_heatmap.py
$PYTHON_BIN scripts/plot_benchmark_reporting_heatmap.py
$PYTHON_BIN scripts/precision_recall_tradeoff.py
$PYTHON_BIN scripts/new_tool_benchmark_performance.py
$PYTHON_BIN scripts/edgren_rnafusion_validation.py

# Optional exploratory random-effects summaries and forest plots.
# Requires R and the R package 'metafor'.
if command -v Rscript >/dev/null 2>&1; then
  $PYTHON_BIN scripts/meta_analysis_R/01_export_from_excel.py
  $PYTHON_BIN scripts/meta_analysis_R/02_data_quality_report.py
  Rscript scripts/meta_analysis_R/03_meta_analysis.R
  $PYTHON_BIN scripts/meta_analysis_R/05_edgren_drift_figure.py
  $PYTHON_BIN scripts/meta_analysis_R/04_make_figures.py
else
  echo "Rscript not found; skipped exploratory random-effects/forest-plot workflow."
fi
