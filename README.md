# RNA-seq fusion detection benchmark systematic review

This repository contains curated extraction tables, analysis scripts, derived tables, and figures for a systematic comparative review of published benchmarks of RNA-seq-based gene fusion detection tools.

The analysis evaluates how reported fusion detection performance depends on benchmark design, including dataset origin, truth-set definition, benchmark purpose, read length, and reproducibility-related reporting.

## Revision for Briefings in Bioinformatics (BIB-26-2487)

Working documents of the revision (Google Docs; the manuscript text is not kept in this repository, and the documents open only for people the authors have shared them with):

- [Revised manuscript with tracked changes, revised figures and Supplementary Information](https://docs.google.com/document/d/1LUmU_7wDaFPf4BdnOTXALAGnI0l_qqQIpkkNqYDqzy0/edit)
- [Point-by-point response to the Editor and Reviewers](https://docs.google.com/document/d/1MGvnk8UXNP7d1znYfxTnQRBeWuIaCPH-i1hiSPtbC3Q/edit)
- [Revision checklist (meeting 21 November 2026; submission due 27 November 2026)](https://docs.google.com/document/d/1acQmuGHdMkKgqNr77-GbjcNb-acgcAPzkk10qhKlF68/edit)
- [Code and data changes on this branch compared with `main`](https://github.com/biodatageeks/fusion-benchmark-systematic-review/compare/main...revision/bib-26-2487)

The revision analyses are described in [Revision analyses](#revision-analyses-response-to-reviewers-bib-26-2487) below.

## Repository structure

```text
data/
  raw/          Curated input workbooks and archived Edgren workflow execution records.
scripts/        Python scripts for main analyses and manuscript figures.
results/
  real_simulated/              Mixed-effects real-vs-simulated analyses.
  benchmark_statistics/        Benchmark landscape heatmap and source tables.
  benchmark_reporting/         Reporting/truth-set heatmap and source tables.
  precision_recall_tradeoff/   Tool-level recall-precision summaries.
  new_tool_performance/        New-tool benchmark context summaries.
  figures/                     Exploratory heterogeneity and Edgren drift figures.
  tables/                      Exploratory pooled/drift/sensitivity tables.
supplementary/  Suggested supplementary/additional-file manifest.
```

## Main input files

- [`data/raw/benchmark.xlsx`](data/raw/benchmark.xlsx) - curated extraction used from the revision onwards (13 benchmarks; one row per tool x dataset x condition in `Tools_performance`; dataset descriptions and colour-coded shared data sources in `Datasets`).
- [`data/processed/analysis_input.xlsx`](data/processed/analysis_input.xlsx) - analysis input built from `benchmark.xlsx` by [`scripts/build_analysis_input.py`](scripts/build_analysis_input.py) (observations, benchmark metadata, dataset metadata, tool-name map). Per-dataset annotations derived by the build (data origin, real-data subtype, read length, biological-dataset identifier) are written to [`data/processed/dataset_annotations.csv`](data/processed/dataset_annotations.csv) for review.
- [`data/raw/master_analysis_input.xlsx`](data/raw/master_analysis_input.xlsx) - previous extraction (11 benchmarks, 383 records), kept for provenance; its source check is documented in [`data/curation/legacy_v1/`](data/curation/legacy_v1/). It still supplies the benchmark metadata of benchmarks 1-11 and the legacy tool-name map.
- [`data/raw/edgren_rnafusion.xlsx`](data/raw/edgren_rnafusion.xlsx) - curated RNAfusion output workbook (pooled predictions) used for the controlled Edgren/RNAfusion validation analysis.

## Edgren workflow execution records

Logs, reports, traces, timelines, parameters and software versions for the recovered RNAfusion **3.0.1** run and the **4.1.3** rerun are archived in [`data/raw/edgren_execution_logs/`](data/raw/edgren_execution_logs/README.md), with source-file mappings and SHA256 checksums. These are provenance records, separate from the curated prediction workbook and analysis inputs. See the archive README for run status and reproducibility limits.

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
python scripts/read_length_analysis.py
python scripts/plot_benchmark_statistics_heatmap.py
python scripts/plot_benchmark_reporting_heatmap.py
python scripts/precision_recall_tradeoff.py
python scripts/new_tool_benchmark_performance.py
python scripts/edgren_rnafusion_validation.py
python scripts/exploratory_random_effects.py
```

## Revision analyses (response to reviewers, BIB-26-2487)

Added in branch `revision/bib-26-2487`. Raw extraction files are not modified; revision-specific curation is kept in separate, reviewable files.

| Item | Script / file | Output |
|---|---|---|
| Conversion of the curated extraction into the analysis input; derived dataset annotations (origin, subtype, read length, biological-dataset identifier from the shared-source colours and simulated-design variants); rows excluded from models with reasons | [`scripts/build_analysis_input.py`](scripts/build_analysis_input.py), [`data/curation/tool_name_map_additions.csv`](data/curation/tool_name_map_additions.csv), [`data/curation/benchmark_metadata_new_studies.csv`](data/curation/benchmark_metadata_new_studies.csv) | [`data/processed/`](data/processed/) |
| Documented metadata corrections (read length Vicente-Garces 2023, Haas 2019 truth-set type) | [`data/curation/metadata_corrections.csv`](data/curation/metadata_corrections.csv) | Figures 3-4 |
| Source check of the previous extraction (field-level corrections and verified biological-dataset map for the legacy workbook) | [`data/curation/legacy_v1/`](data/curation/legacy_v1/) | provenance only |
| Tool publication year and maintenance status (assessed 2026-11-17) | [`data/curation/tool_release_maintenance.csv`](data/curation/tool_release_maintenance.csv) | Supplementary Table S6 |
| Published Edgren reference lists (TP27, TP53, TP99) with cell lines and evidence; SRA run manifest candidates | [`data/edgren/reference_events.csv`](data/edgren/reference_events.csv), [`data/edgren/sample_manifest_candidates.csv`](data/edgren/sample_manifest_candidates.csv) | Figure 7, Supplementary Tables S11-S12 |
| Cross-classified models (tool x biological dataset x benchmark), benchmark-cluster bootstrap, leave-one-benchmark-out, Edgren exclusion, within-benchmark contrasts | [`scripts/crossed_dependence_models.py`](scripts/crossed_dependence_models.py) | [`results/revision/crossed_models/`](results/revision/crossed_models/) |
| Revised Figure 5 | [`scripts/plot_real_simulated_revised.py`](scripts/plot_real_simulated_revised.py) | [`results/revision/figures/`](results/revision/figures/) |
| Descriptive benchmark-level summaries (purpose, read length, real-data subtype, use contexts) | [`scripts/descriptive_benchmark_summaries.py`](scripts/descriptive_benchmark_summaries.py) | [`results/revision/descriptive/`](results/revision/descriptive/) |
| Tool list and tool profiles with numbers of observations | [`scripts/tool_tables.py`](scripts/tool_tables.py) | [`results/revision/tools/`](results/revision/tools/) |
| Consistency of published TP27 vs TP99 Edgren scores | [`scripts/edgren_published_consistency.py`](scripts/edgren_published_consistency.py) | [`results/revision/edgren_published/`](results/revision/edgren_published/) |
| Per-sample, multi-reference Edgren re-scoring of the same predictions | [`scripts/edgren_reference_sensitivity.py`](scripts/edgren_reference_sensitivity.py) + [`data/edgren/`](data/edgren/) | [`results/revision/edgren_reference_sensitivity/`](results/revision/edgren_reference_sensitivity/) |
| Chronological ordering of Figures 3-4 | `scripts/plot_benchmark_*_heatmap.py` | `results/benchmark_*/` |

LaTeX tables for the Supplement are written to [`results/revision/tables/`](results/revision/tables/). The per-sample Edgren analysis needs the per-sample RNAfusion predictions described in [`data/edgren/README.md`](data/edgren/README.md); the archived execution records in [`data/raw/edgren_execution_logs/`](data/raw/edgren_execution_logs/) do not include caller predictions.

## Edgren dataset annotation

Edgren-derived datasets were annotated as: `2.3`, `5.1`, `6.1`, `6.2`, `7.6`, `7.7`, `8.3`, and `10.2`. These annotations are used for the Edgren truth-set heterogeneity case study and for sensitivity analyses excluding Edgren-derived observations.

## Key outputs used in the manuscript

- [`results/benchmark_statistics/benchmark_statistics_heatmap.pdf`](results/benchmark_statistics/benchmark_statistics_heatmap.pdf)
- [`results/benchmark_reporting/benchmark_reporting_truthset_heatmap.pdf`](results/benchmark_reporting/benchmark_reporting_truthset_heatmap.pdf)
- [`results/real_simulated/real_vs_simulated_f1_precision_recall.pdf`](results/real_simulated/real_vs_simulated_f1_precision_recall.pdf)
- [`results/read_length/read_length_f1.pdf`](results/read_length/read_length_f1.pdf)
- [`results/real_simulated/dataset_type_models.csv`](results/real_simulated/dataset_type_models.csv)
- [`results/real_simulated/dataset_type_models_without_edgren.csv`](results/real_simulated/dataset_type_models_without_edgren.csv)
- [`results/real_simulated/leave_one_benchmark_out_f1.csv`](results/real_simulated/leave_one_benchmark_out_f1.csv)
- [`results/precision_recall_tradeoff/tool_precision_recall_summary_real_without_edgren.csv`](results/precision_recall_tradeoff/tool_precision_recall_summary_real_without_edgren.csv)
- [`results/new_tool_performance/new_tool_benchmark_summary.csv`](results/new_tool_performance/new_tool_benchmark_summary.csv)
- [`results/edgren_rnafusion_validation/edgren_rnafusion_metrics_summary.csv`](results/edgren_rnafusion_validation/edgren_rnafusion_metrics_summary.csv)
- [`results/edgren_rnafusion_validation/edgren_rnafusion_metrics.pdf`](results/edgren_rnafusion_validation/edgren_rnafusion_metrics.pdf)
- [`results/figures/edgren_drift_scatter.pdf`](results/figures/edgren_drift_scatter.pdf)
- [`results/figures/heterogeneity_I2.pdf`](results/figures/heterogeneity_I2.pdf)
- [`results/tables/pooled_estimates.csv`](results/tables/pooled_estimates.csv)

## Data provenance

All values were manually extracted or curated from published benchmark articles and associated supplementary materials. The workbooks in [`data/raw/`](data/raw/) are the primary curated data sources. CSV files and figures in [`results/`](results/) are derived outputs generated by the scripts.

## Notes for publication

Before public release, add:

- final citation metadata for the associated manuscript;
- a permanent repository archive DOI, e.g. Zenodo;
- a license selected by the authors and institution;
- any final corrections to raw extraction tables made during manuscript revision.
