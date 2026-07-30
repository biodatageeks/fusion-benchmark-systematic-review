#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python3}"

$PYTHON_BIN scripts/real_simulated_sensitivity.py
$PYTHON_BIN scripts/plot_real_simulated_metrics.py
$PYTHON_BIN scripts/read_length_analysis.py
$PYTHON_BIN scripts/plot_edgren_truthset_heterogeneity.py
$PYTHON_BIN scripts/plot_benchmark_statistics_heatmap.py
$PYTHON_BIN scripts/plot_benchmark_reporting_heatmap.py
$PYTHON_BIN scripts/precision_recall_tradeoff.py
$PYTHON_BIN scripts/new_tool_benchmark_performance.py
$PYTHON_BIN scripts/edgren_rnafusion_validation.py
$PYTHON_BIN scripts/exploratory_random_effects.py
