#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python3}"
# fixed hash seed: set iteration order (and so the row order of some Edgren tables) is the same in every run
export PYTHONHASHSEED=0
# fixed timestamp in PDF metadata, so that figures are byte-identical between runs
export SOURCE_DATE_EPOCH=0

$PYTHON_BIN scripts/real_simulated_sensitivity.py
$PYTHON_BIN scripts/crossed_models.py
$PYTHON_BIN scripts/plot_real_simulated_metrics.py
$PYTHON_BIN scripts/read_length_analysis.py
$PYTHON_BIN scripts/plot_edgren_truthset_heterogeneity.py
$PYTHON_BIN scripts/plot_prisma_flow.py
$PYTHON_BIN scripts/plot_benchmark_statistics_heatmap.py
$PYTHON_BIN scripts/plot_benchmark_reporting_heatmap.py
$PYTHON_BIN scripts/precision_recall_tradeoff.py
$PYTHON_BIN scripts/new_tool_benchmark_performance.py
$PYTHON_BIN scripts/edgren_rnafusion_validation.py
$PYTHON_BIN scripts/edgren_per_sample_scoring.py
$PYTHON_BIN scripts/supplementary_tables.py
$PYTHON_BIN scripts/export_source_csv.py
$PYTHON_BIN scripts/exploratory_random_effects.py
