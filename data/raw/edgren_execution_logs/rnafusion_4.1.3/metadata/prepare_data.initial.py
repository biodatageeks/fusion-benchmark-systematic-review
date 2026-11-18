#!/usr/bin/env python3
"""Fetch ENA metadata and checksum-verified paired FASTQ for the Edgren runs."""
import concurrent.futures
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess
import time
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
RUNS = {
    "SRR064286": "MCF-7", "SRR064287": "KPL-4",
    "SRR064438": "BT-474", "SRR064439": "BT-474",
    "SRR064440": "SK-BR-3", "SRR064441": "SK-BR-3",
}
FIELDS = (
    "run_accession,sample_accession,experiment_accession,study_accession,"
    "sample_title,experiment_title,instrument_model,library_layout,library_strategy,"
    "library_source,library_selection,read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes"
)


def fetch(url):
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                return response.read()
        except Exception:
            if attempt == 2:
                raise
            time.sleep(5)


def md5(path):
    h = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def download(item):
    url, expected, size = item
    path = ROOT / "fastq" / url.rsplit("/", 1)[1]
    if path.exists():
        if path.stat().st_size == size and md5(path) == expected:
            return path.name
        raise RuntimeError(f"Existing FASTQ fails verification: {path}")
    part = path.with_name(path.name + ".part")
    print(f"Downloading {path.name}: {size:,} bytes", flush=True)
    subprocess.run([
        "curl", "--fail", "--location", "--retry", "5", "--retry-delay", "5",
        "--connect-timeout", "30", "--speed-limit", "1024", "--speed-time", "180",
        "--silent", "--show-error", "--continue-at", "-", "--output", str(part),
        "https://" + url,
    ], check=True)
    if part.stat().st_size != size or md5(part) != expected:
        raise RuntimeError(f"Downloaded FASTQ fails size/MD5 verification: {part}")
    subprocess.run(["gzip", "-t", str(part)], check=True)
    part.rename(path)
    print(f"Verified {path.name}", flush=True)
    return path.name


def main():
    for name in ("metadata", "fastq", "logs"):
        (ROOT / name).mkdir(exist_ok=True)
    rows, downloads = [], []
    for run, cell in RUNS.items():
        url = ("https://www.ebi.ac.uk/ena/portal/api/filereport?"
               f"accession={run}&result=read_run&fields={FIELDS}&format=tsv")
        raw = fetch(url)
        (ROOT / "metadata" / f"{run}.ena.tsv").write_bytes(raw)
        parsed = list(csv.DictReader(io.StringIO(raw.decode()), delimiter="\t"))
        assert len(parsed) == 1, (run, parsed)
        row = parsed[0]
        assert row["run_accession"] == run
        assert row["library_layout"] == "PAIRED"
        assert row["library_source"] == "TRANSCRIPTOMIC"
        assert cell in row["sample_title"], (run, row["sample_title"])
        row["cell_line"] = cell
        experiment = row["experiment_accession"]
        (ROOT / "metadata" / f"{experiment}.xml").write_bytes(
            fetch(f"https://www.ebi.ac.uk/ena/browser/api/xml/{experiment}"))
        paths = row["fastq_ftp"].split(";")
        checksums = row["fastq_md5"].split(";")
        sizes = [int(x) for x in row["fastq_bytes"].split(";")]
        assert len(paths) == len(checksums) == len(sizes) == 2
        assert [p.rsplit("/", 1)[1] for p in paths] == [f"{run}_{i}.fastq.gz" for i in (1, 2)]
        downloads.extend(zip(paths, checksums, sizes))
        rows.append(row)
    with (ROOT / "metadata" / "ena_runs.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=["cell_line"] + FIELDS.split(","), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Total compressed FASTQ: {sum(x[2] for x in downloads):,} bytes", flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(download, downloads))
    with (ROOT / "samplesheet.csv").open("w") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample", "fastq_1", "fastq_2", "strandedness"])
        for run in RUNS:
            # ENA does not specify RNA strand specificity. The schema permits unknown.
            writer.writerow([run, *(str(ROOT / "fastq" / f"{run}_{i}.fastq.gz") for i in (1, 2)), "unknown"])
    (ROOT / "metadata" / "download_complete.json").write_text(json.dumps({
        "runs": list(RUNS), "files": len(downloads),
        "bytes": sum(x[2] for x in downloads), "verification": "size, MD5, gzip -t",
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
