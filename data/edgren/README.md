# Edgren per-sample, multi-reference re-scoring (R1.2, R1.7, R3.3)

Inputs expected by `scripts/edgren_reference_sensitivity.py`. These files must be prepared by the
authors; they cannot be derived from the current repository content.

1. `sample_manifest.csv` (copy `sample_manifest_TEMPLATE.csv`)
   - `sample_id`: sample name used by nf-core/rnafusion
   - `run_accession`: SRA run (SRP003186)
   - `cell_line`: BT-474, KPL-4, MCF-7 or SK-BR-3, spelled as in `reference_events.csv`
   - `prediction_file`: path (relative to this folder) of the per-sample fusion-report CSV with
     columns `Fusion` (GENE1--GENE2) and one column per caller (`arriba`, `fusioncatcher`, `starfusion`);
     a non-empty cell means that the caller reported the fusion in that sample.
2. `reference_events.csv` (copy `reference_events_TEMPLATE.csv`), one row per reference event and list
   - `list_id`: TP27, TP37, TP53, TP80, TP99 (a gene pair present in several lists appears once per list)
   - `gene1`, `gene2`: HGNC symbols; aliases can be separated by `|`
   - `cell_line`: cell line in which the event is reported as true
   - `source`: publication/supplementary table the event was taken from
   - `evidence`: e.g. RT-PCR, Sanger, WGS, literature, reanalysis
3. `nfcore_rnafusion_run/`: provenance of the documented re-run (command, params.json, nextflow.config,
   pipeline_info/ with software versions, execution report, trace and timeline, container digests,
   reference bundle description with checksums).

Run:

```bash
python scripts/edgren_reference_sensitivity.py --unit sample      # default, per sample
python scripts/edgren_reference_sensitivity.py --unit cell_line   # union of runs per cell line
python scripts/edgren_reference_sensitivity.py --legacy           # pooled 99-fusion check on data/raw/edgren_rnafusion.xlsx
```
