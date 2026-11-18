#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
source /data/ccle_fusions/config/paths.env
export TMPDIR="$ROOT/tmp"
export NXF_WORK="$ROOT/work"
mkdir -p "$TMPDIR" "$ROOT/logs" "$ROOT/run" "$NXF_WORK"
cd "$ROOT/run"
finish() {
    code=$?
    printf '%s\n' "$code" > "$ROOT/exit_code.txt"
    date -u +%FT%TZ > "$ROOT/finished_at.txt"
    if [[ $code -eq 0 ]]; then echo completed > "$ROOT/stage.txt"; else echo failed > "$ROOT/stage.txt"; fi
}
trap finish EXIT
echo downloading > "$ROOT/stage.txt"
python3 -u "$ROOT/scripts/prepare_data.py"
echo rnafusion > "$ROOT/stage.txt"
nextflow run nf-core/rnafusion \
    -r 4.1.3 -profile docker -c "$ROOT/rnafusion.config" \
    --input "$ROOT/samplesheet.csv" \
    --tools arriba,starfusion,fusioncatcher,fusionreport \
    --no_cosmic --genomes_base /data/ccle_fusions/references/rnafusion \
    --outdir "$ROOT/results" --publish_dir_mode copy \
    --max_cpus 16 --max_memory '60 GB' --max_time 240.h \
    -with-report "$ROOT/logs/report.html" \
    -with-trace "$ROOT/logs/trace.tsv" \
    -with-timeline "$ROOT/logs/timeline.html" \
    -resume
