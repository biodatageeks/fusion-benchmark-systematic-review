# Benchmark design shapes reported performance of short-read RNA-seq gene fusion detection tools

Data, code and results for the systematic review and comparative reanalysis *Benchmark design shapes reported performance of short-read RNA-seq gene fusion detection tools: a systematic review and comparative reanalysis* (Ostrowska, Szmurło, Gambin; Briefings in Bioinformatics).

The analysis covers 18 published benchmarks of short-read RNA-seq fusion detection tools. It asks how reported performance depends on benchmark design: dataset origin (real vs simulated), truth-set definition, benchmark purpose, read length and reproducibility-related reporting. A controlled case study re-scores RNAfusion outputs on the Edgren breast cancer cell-line data against four published reference lists.

## Repository structure

```text
data/raw/
  master_analysis_input.xlsx   curated extraction workbook: the input of all analyses
  edgren_rnafusion.xlsx        RNAfusion calls (Arriba, STAR-Fusion, FusionCatcher) on the Edgren data, pooled, with the 99-fusion list
  edgren_predictions/          caller outputs per sequencing run, nf-core/rnafusion 3.0.1 and 4.1.3 (36 CSV, manifest with SHA-256, ENA run metadata)
  edgren_truthsets/            published Edgren reference lists (27, 38, 53 and 99 fusions) per sequencing run
  edgren_execution_logs/       Nextflow logs, execution reports, parameters and software versions of both RNAfusion runs
data/extraction_check/          second-author check of 31 extracted records against the source publications
scripts/                       one script per analysis (run in order by run_all.sh)
tests/                         checks of the input data, the reported counts, Edgren scoring, supplementary tables and models
results/                       outputs of run_all.sh (CSV tables, summaries, figures as PDF and PNG)
```

`master_analysis_input.xlsx` has seven sheets:
- `observations`: 615 tool–dataset–condition observations with TP, FP, FN, precision, recall and F1-score as reported or harmonized;
- `benchmark_metadata`: 18 benchmarks, including purpose, truth-set type and reproducibility annotations;
- `dataset_metadata`: 84 dataset entries, including origin, read length and `biological_dataset_group`, which marks the same biological data used in several benchmarks;
- `tool_name_map`: harmonized tool names;
- the two Edgren sheets and a README sheet.

The same sheets are exported as CSV to `results/source_csv/`.

## Reproducing the analyses

Python 3.10 or later:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
./run_all.sh          # about 5 minutes; rewrites results/
python -m pytest      # about 10 seconds
```

`run_all.sh` regenerates every file in `results/`. CSV tables and PNG figures are identical between runs; PDF files differ only in their creation date.

## Figures and tables of the article

| Article item | File | Script |
|---|---|---|
| Figure 2 (PRISMA flow diagram) | `results/prisma/prisma_flow.pdf` | `plot_prisma_flow.py` |
| Figure 3 (benchmark statistics) | `results/benchmark_statistics/benchmark_statistics_heatmap.pdf` | `plot_benchmark_statistics_heatmap.py` |
| Figure 4 (reporting and truth sets) | `results/benchmark_reporting/benchmark_reporting_truthset_heatmap.pdf` | `plot_benchmark_reporting_heatmap.py` |
| Figure 5 (real vs simulated) | `results/real_simulated/real_vs_simulated_f1_precision_recall.pdf` | `plot_real_simulated_metrics.py`, `crossed_models.py` |
| Figure 6 (truth-set-induced F1 changes) | `results/figures/edgren_drift_scatter.pdf` | `exploratory_random_effects.py` |
| Figure 7 (RNAfusion on the Edgren data) | `results/edgren_per_sample/edgren_metrics_by_list.pdf` | `edgren_per_sample_scoring.py` |
| Supplementary Figure S1 (read length) | `results/read_length/read_length_f1.pdf` | `read_length_analysis.py` |
| Supplementary Figure S2 (recall–precision profiles) | `results/precision_recall_tradeoff/tool_precision_recall_tradeoff.pdf` | `precision_recall_tradeoff.py` |
| Supplementary Figure S3 (I²) | `results/figures/heterogeneity_I2.pdf` | `exploratory_random_effects.py` |
| Supplementary Figure S4 (Edgren truth sets) | `results/real_simulated/edgren_truthset_heterogeneity.pdf` | `plot_edgren_truthset_heterogeneity.py` |
| Supplementary Figure S5 (crossed models) | `results/crossed_models/crossed_models_forest.pdf` | `crossed_models.py` |
| Supplementary Figure S6 (per-run F1) | `results/edgren_per_sample/edgren_per_run_f1.pdf` | `edgren_per_sample_scoring.py` |
| Supplementary Tables S3–S7 | `results/supplementary_tables/S3_sra_runs.csv` … `S7_benchmark_contexts.csv` (LaTeX rows in `latex/`) | `supplementary_tables.py` |

Figure 1 is a schematic. Supplementary Tables S1 and S2 are descriptive and contain no computed values; the PRISMA counts are listed in `plot_prisma_flow.py`.

**Statistical results reported in the text:**
- `results/crossed_models/` (`summary.md`, `real_vs_simulated_models.csv`, `real_subtype_models.csv`, `variance_components.csv`, `sensitivity_lobo_edgren.csv`, `within_benchmark_summary.csv`, `benchmark_purpose_*.csv`);
- `results/real_simulated/dataset_type_models.csv`;
- `results/read_length/`;
- `results/new_tool_performance/`;
- `results/edgren_rnafusion_validation/`;
- `results/edgren_per_sample/`;
- `results/tables/`.

## Main analyses

- **Real vs simulated data:** mixed models of logit-transformed metrics with crossed random intercepts for tool, biological dataset group and benchmark. Confidence intervals come from a cluster bootstrap that resamples whole benchmarks (999 replicates). Sensitivity analyses: beta regression, a model with benchmark as a fixed covariate, leave-one-benchmark-out, exclusion of Edgren-derived observations, within-benchmark contrasts, and real-data subtypes (cell line, patient-derived, mixed).
- **Read length:** Kruskal-Wallis test and a mixed model with read length (reference 100 bp) and the same crossed random intercepts.
- **Benchmark purpose:** descriptive comparison of new-tool and independent benchmarks, overall and for the same tools.
- **Edgren case study:**
  - RNAfusion calls for six SRA runs (SRP003186) from four cell lines;
  - each run scored against the events of its own cell line in the 27-, 38-, 53- and 99-fusion lists (unique undirected gene pairs);
  - gene symbols matched with their previous symbols; MED1--UPS32 in the 99-fusion list is a misspelling of USP32 and is matched as MED1--USP32;
  - outputs of nf-core/rnafusion 3.0.1 and 4.1.3, which used different caller versions and reference resources, scored separately.

## Data provenance

- **Extraction:** all benchmark values were extracted manually from the published articles and their supplementary materials. A second author independently verified a sample of 31 records against the sources (`data/extraction_check/`).
- **Edgren data:** the caller outputs in `data/raw/edgren_predictions/` keep the columns, rows and values of the original outputs (tab-separated files converted to CSV). `data/raw/edgren_execution_logs/README.md` describes what the execution records establish and their limits.

## Citation

Ostrowska I, Szmurło A, Gambin T. Benchmark design shapes reported performance of short-read RNA-seq gene fusion detection tools: a systematic review and comparative reanalysis. *Briefings in Bioinformatics*.

The DOI of the article and of the archived version of this repository (Zenodo) will be added on publication.
