# Edgren caller predictions

Predictions for six sequencing runs, three callers and two nf-core/rnafusion workflow versions. Files are organized as `<workflow version>/<SRA run>/<caller>.csv`. Each CSV retains the columns, rows and field values of the corresponding caller output; only the tab-separated format has been converted to comma-separated CSV. No scoring, filtering, pooling or normalization was applied. The manifest records source filenames, SHA-256 hashes of the original files and prediction row counts.

| SRA run | Cell line |
| --- | --- |
| SRR064286 | MCF-7 |
| SRR064287 | KPL-4 |
| SRR064438 | BT-474 |
| SRR064439 | BT-474 |
| SRR064440 | SK-BR-3 |
| SRR064441 | SK-BR-3 |

The workflow versions used different caller versions and reference resources; differences cannot be attributed to workflow version alone. See `../edgren_execution_logs/README.md` for execution provenance and its limitations. These are six sequencing runs from four cell lines. The source transfer was not fully verified; these 36 input tables matched previously recorded SHA-256 hashes and were structurally checked.
