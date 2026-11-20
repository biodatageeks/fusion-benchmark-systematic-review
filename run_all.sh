#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python3}"
# Reproducible outputs: fixed hash seed (set ordering in CSV files) and fixed PDF timestamps
export PYTHONHASHSEED=0
export SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-1767225600}"

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
$PYTHON_BIN scripts/edgren_reference_sensitivity.py   # per-sample scoring of the 3.0.1 and 4.1.3 runs (Figure 7, Table S11)

# Checks of scoring rules, input consistency and reported numbers
$PYTHON_BIN -m pytest -q tests
