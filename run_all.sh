#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python3}"

# Build the analysis input from the curated extraction (data/raw/benchmark.xlsx)
$PYTHON_BIN scripts/build_analysis_input.py

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

# --- Analyses added during revision (BIB-26-2487) ---------------------------------------------
# Curated inputs: data/curation/bio_dataset_map.csv, data/curation/metadata_corrections.csv
$PYTHON_BIN scripts/crossed_dependence_models.py --n-boot "${N_BOOT:-1999}"
$PYTHON_BIN scripts/plot_real_simulated_revised.py
$PYTHON_BIN scripts/descriptive_benchmark_summaries.py
$PYTHON_BIN scripts/tool_tables.py
$PYTHON_BIN scripts/edgren_published_consistency.py
$PYTHON_BIN scripts/edgren_reference_sensitivity.py --legacy
$PYTHON_BIN scripts/edgren_reference_sensitivity.py --pooled
if [[ -f data/edgren/sample_manifest.csv && -f data/edgren/reference_events.csv ]]; then
  $PYTHON_BIN scripts/edgren_reference_sensitivity.py --unit sample
  $PYTHON_BIN scripts/edgren_reference_sensitivity.py --unit cell_line
else
  echo "Skipping per-sample Edgren re-scoring: data/edgren inputs not yet provided (see data/edgren/README.md)."
fi
