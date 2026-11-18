# Edgren RNAfusion execution records

Recovered execution logs and supporting provenance for the six public Edgren tumour-cell-line RNA-seq runs (SRA project **SRP003186**). Packaged on 18 November 2026. Source files are copied byte for byte; no text was redacted or reconstructed. Private-key, GitHub-token, AWS-access-ID, authorization-header, signed-URL and URL-credential pattern checks found no matches in the selected files. The logs retain the original host/user/path and execution metadata.

## Runs

| Workflow | Recorded execution | Pipeline revision | Nextflow | Successful tasks | Failed tasks |
| --- | --- | --- | --- | ---: | ---: |
| nf-core/rnafusion 3.0.1 | 2 April 2024, recovered `run3` | `98f02be34c` (launched from `master`) | 23.10.0, build 5889 | 105 | 0 |
| nf-core/rnafusion 4.1.3 | 31 October–1 November 2026 | `76ad76e7c3` (tag `4.1.3`; full commit in `metadata/pipeline_commit.txt`) | 25.10.4, build 11173 | 91 | 0 |

Versions are established by the recorded software-version files and main Nextflow logs. The recovered 2024 run is **3.0.1**, not 3.0.2. This archive alone does not establish which run generated the previously curated `../edgren_rnafusion.xlsx` workbook; that mapping requires comparison with its source predictions.

The successful-run logs explicitly record workflow completion. Earlier attempts and reports are retained as separate original files; they should not be confused with the successful executions. Timestamps are preserved as recorded in the source files and may use different time zones. The old params files with different timestamped names have identical SHA256 hashes, as do the two old DAG files; these are preserved rather than counted as different configurations.

## Samples

| SRA run | Cell line |
| --- | --- |
| SRR064286 | MCF-7 |
| SRR064287 | KPL-4 |
| SRR064438 | BT-474 |
| SRR064439 | BT-474 |
| SRR064440 | SK-BR-3 |
| SRR064441 | SK-BR-3 |

These are six sequencing runs from four cell lines, not six independent biological cohorts. The normal breast sample is not part of these runs. Original samplesheets and, for 4.1.3, the ENA/SRA metadata and download-completion record are included.

## Contents

- `rnafusion_3.0.1/nextflow/`: main `.nextflow.log` and the preceding rotated log, stored with visible filenames.
- `rnafusion_3.0.1/pipeline_info/`: execution reports, traces, timelines, DAGs, parameter snapshots, validated samplesheet and software versions, under their original names.
- `rnafusion_4.1.3/nextflow/`: main Nextflow log, stored with a visible filename.
- `rnafusion_4.1.3/logs/`: launcher/resume output, execution report, trace, timeline and retained initial-attempt status files.
- `rnafusion_4.1.3/pipeline_info/`: software versions, parameters and DAG.
- `rnafusion_4.1.3/metadata/`: public input metadata, pipeline commit, download record, pipeline input schema and original script checksums.
- `rnafusion_4.1.3/scripts/`: recorded preparation and launch scripts; their absolute paths describe the original host and require adaptation on a new host.
- `rnafusion_4.1.3/rnafusion.config`, `samplesheet.csv`, `exit_code.txt`, `finished_at.txt`: configuration, input manifest and completion records.
- `manifest.json` and `files.tsv`: source-relative file mapping, file sizes and SHA256 hashes for all 51 payload files; the manifest also records HTML end-tag checks.
- `SHA256SUMS`: hashes of every archived file except the checksum file itself.

## Verification

From this directory on a system with `sha256sum`:

```sh
sha256sum -c SHA256SUMS
```

On macOS, use `shasum -a 256 -c SHA256SUMS` instead. HTML execution reports can be downloaded and opened locally. Their embedded JavaScript may use external resources, so offline rendering can vary.

## Scope and remaining limits

Both recovered 3.0.1 `execution_report_*.html` files appear truncated: each is exactly 2,359,296 bytes and lacks a closing HTML tag. They are retained as recovered evidence, not certified complete reports. The main Nextflow log records successful completion, and the trace, timelines, DAGs and software-version records are available separately. The 4.1.3 report and all archived timeline/DAG HTML files have closing HTML tags; this is a structural check, not full interactive-rendering validation.

This is an execution-provenance archive. FASTQ, genome/annotation bundles, Nextflow work/cache directories and caller predictions are not included. Copying these records does not certify the completeness of the entire local transfer or establish byte-level identity of FASTQ between workflow versions. The 4.1.3 download record documents the checks performed during preparation; this packaging step independently checks the archived records only.

Caller versions and reference resources changed between 3.0.1 and 4.1.3, so their difference cannot be attributed solely to the workflow version. The available metadata indicate GRCh38/Ensembl 102 for the older run and GRCh38/GENCODE v46 resources for the newer run. Reference-bundle checksums and a complete set of immutable container digests have not been independently verified here. These logs therefore improve provenance but do not, by themselves, establish exact end-to-end reproducibility.

For a controlled reference-definition sensitivity analysis, hold predictions and scoring rules fixed within one run/version while changing only the reference list. A comparison of the two workflow versions is a separate experiment.
