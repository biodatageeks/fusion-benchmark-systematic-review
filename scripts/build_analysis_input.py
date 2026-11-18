"""Build the analysis input workbook from the curated extraction in data/raw/benchmark.xlsx.

data/raw/benchmark.xlsx (sheets `Tools_performance`, `Datasets`) is the primary extraction used from the
revision onwards: 13 benchmarks, one row per tool x dataset x condition, with the shared data sources
colour-coded in the `Datasets` sheet. The analysis scripts expect the schema of the previous
workbook (observations / benchmark_metadata / dataset_metadata / tool_name_map / edgren_truthsets), so this
script converts the extraction into data/processed/analysis_input.xlsx and writes the per-dataset
annotations it derived to data/processed/dataset_annotations.csv for review.

Every classification rule is explicit below. Rows that would double-count other rows (derived aggregates),
long-read data and non-primary truth thresholds are kept in the observations sheet but flagged with
`include_in_models = False` and a reason; the analysis loaders drop them.

Benchmark-level metadata come from the previous workbook (benchmarks 1-11, with the documented corrections in
data/curation/metadata_corrections.csv applied by the plotting scripts) and from
data/curation/benchmark_metadata_new_studies.csv (Liu 2016, Creason 2021); counts of tools and datasets are
recomputed from the new extraction.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
RAW = REPO / "data" / "raw" / "benchmark.xlsx"
LEGACY = REPO / "data" / "raw" / "master_analysis_input.xlsx"
CURATION = REPO / "data" / "curation"
OUT_DIR = REPO / "data" / "processed"
OUT = OUT_DIR / "analysis_input.xlsx"

BENCHMARK_NUMBER = {
    "Haas 2019": 1, "Uhrig 2021": 2, "Kumar 2016": 3, "Balan 2021": 4, "Zhao 2017": 5, "Vu 2018": 6,
    "Singh and Li 2021": 7, "Apostolides 2021": 8, "Zhang 2016": 9, "Davidson 2015": 10,
    "Vicente-Garcés 2023": 11, "Liu 2016": 12, "Creason 2021 (DREAM)": 13,
}
DOI = {12: "10.1093/nar/gkv1234", 13: "10.1016/j.cels.2021.05.021"}

# Shared data sources, colour-coded by the extractor in the `Datasets` sheet (legend in columns F:H).
COLOUR_GROUP = {
    "FFFFD4DF": "EDGREN_BRCA_CELL_LINES", "FFD8CCF0": "HAAS2019_SIMULATIONS", "FFFFE0B2": "BERGER2010_MELANOMA",
    "FFCDE8D5": "BAO2014_GLIOMA_13", "FFFFF0A6": "QIN2015_LNCAP", "FFBFE6F1": "TEMBE2014_COLO829_SPIKEIN",
    "FFBCD3F5": "FUSIONMAP_50_SYNTHETIC", "FFD7E9B4": "ENCODE_MCF7_SRR534293", "FFD8DDE3": "BEERS_NEGATIVE",
}
# Liu 2016 datasets are named differently in Tools_performance and in Datasets.
LIU_NAME = {
    "Edgren (breast cancer; 4 cell lines pooled)": "Edgren",
    "Melanoma (6 samples pooled)": "Berger",
    "Prostate cancer (5 specimens pooled)": "Yu",
    "Type-1A (100X; insert 500 +/- 50 bp)": "Type-1A (Liu 2016)",
    "Type-1B (100X; insert 250 +/- 50 bp)": "Type-1B (Liu 2016)",
    "Type-3B (mean of 5 tissue backgrounds)": "Type-3B (Liu 2016)",
}
# Variants of one simulated design (coverage, insert size or read length) share one biological identifier.
DESIGN_GROUP = [
    (r"^Uhrig 2021: simulated", "UHRIG2021_SIM_H1"),
    (r"^Davidson 2015: MiSeq", "DAVIDSON2015_MISEQ_SIM"),
    (r"^Singh 2021: ChimPipe simulated", "SINGH2021_CHIMSIM"),
    (r"^Type-[13][AB] \(Liu 2016\)", "LIU2016_SIM150"),
    (r"^Balan 2021", "BALAN2021_QIASEQ"),
]
# DatasetType -> analysis origin; anything not listed is excluded with the reason shown.
ORIGIN = {
    "Cell-line": ("Real", "cell line"), "Cell line": ("Real", "cell line"),
    "Patient-derived": ("Real", "patient-derived"),
    "Cell-line + patient-derived cultures": ("Real", "mixed"), "Real": ("Real", "mixed"),
    "Cell line + patient-derived (aggregate)": ("Real", "mixed"),
    "Synthetic": ("Simulated", "simulated"), "Synthetic positive": ("Simulated", "simulated"),
    "Synthetic mixed": ("Simulated", "simulated"), "Synthetic + real background": ("Simulated", "simulated"),
    "Cell-line + spike-in": ("spike-in", "spike-in"), "Synthetic spike-in + cell-line": ("spike-in", "spike-in"),
    "Synthetic spike-in; targeted RNA panel": ("spike-in", "spike-in"),
    "Targeted RNA panel; clinical + spike-in": ("spike-in", "spike-in"),
    "Synthetic negative": ("negative control", "negative control"), "Normal control": ("negative control", "negative control"),
}
# Rows kept in the extraction but not used in models (reason recorded).
EXCLUDE = [
    (lambda r: r.Benchmark == "Haas 2019" and "Consensus" in str(r.Truthset) and not str(r.Truthset).startswith("Consensus >= 7 "),
     "non-primary consensus threshold (primary: >= 7 methods, as in Haas et al.)"),
    (lambda r: r.Dataset == "Zhao 2017: breast cancer lines (three pooled)", "derived aggregate of the three per-line datasets"),
    (lambda r: r.Dataset == "Vicente-Garcés 2023: global aggregate", "aggregate of the cell-line and patient datasets"),
    (lambda r: "PacBio" in str(r.Dataset), "long-read data, outside the review's scope"),
]


def num(value: object) -> float:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return np.nan
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    return float(text) if re.fullmatch(r"-?\d+(\.\d+)?", text) else np.nan


def read_length(value: object) -> tuple[float, str | None]:
    text = str(value) if value is not None else ""
    layout = "SE" if re.search(r"\bSE\b", text) else ("PE" if "PE" in text else None)
    m = re.match(r"^\s*(\d{2,3})\s*(PE|SE)\b", text)
    return (float(m.group(1)) if m else np.nan), layout


def dataset_colours() -> dict[str, str | None]:
    ws = openpyxl.load_workbook(RAW)["Datasets"]
    colours = {}
    for r in range(2, ws.max_row + 1):
        name = ws.cell(r, 1).value
        if name is None:
            continue
        cell = ws.cell(r, 1)
        colours[str(name).strip()] = cell.fill.fgColor.rgb if cell.fill is not None and cell.fill.fill_type else None
    return colours


def tool_map(raw_tools: pd.Series) -> pd.DataFrame:
    legacy = pd.read_excel(LEGACY, sheet_name="tool_name_map")
    legacy_map = {str(a): (b, c, bool(d)) for a, b, c, d in legacy[["tool_reported", "tool_name_clean", "tool_family", "is_variant_or_output"]].itertuples(index=False)}
    legacy_map.update({str(b): (b, c, bool(d)) for _, b, c, d in legacy[["tool_reported", "tool_name_clean", "tool_family", "is_variant_or_output"]].itertuples(index=False)})
    extra = pd.read_csv(CURATION / "tool_name_map_additions.csv")
    legacy_map.update({str(r.tool_reported): (r.tool_name_clean, r.tool_family, bool(r.is_variant_or_output)) for r in extra.itertuples()})
    rows = []
    for name in sorted(raw_tools.astype(str).str.strip().unique()):
        if name not in legacy_map:
            raise KeyError(f"Tool name not mapped: {name!r}; add it to data/curation/tool_name_map_additions.csv")
        clean, family, variant = legacy_map[name]
        rows.append({"tool_reported": name, "tool_name_clean": clean, "tool_family": family, "is_variant_or_output": variant})
    return pd.DataFrame(rows)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    perf = pd.read_excel(RAW, sheet_name="Tools_performance")
    perf["Benchmark"] = perf["Benchmark"].astype(str).str.strip()
    perf["Dataset"] = perf["Dataset"].astype(str).str.strip()
    perf["Tool"] = perf["Tool"].astype(str).str.strip()
    unknown = sorted(set(perf["Benchmark"]) - set(BENCHMARK_NUMBER))
    if unknown:
        raise KeyError(f"Benchmarks without a number: {unknown}")
    colours = dataset_colours()

    # --- datasets ------------------------------------------------------------------------------------
    ds_rows = []
    for (bench, name), g in perf.groupby(["Benchmark", "Dataset"], sort=False):
        b = BENCHMARK_NUMBER[bench]
        sheet_name = LIU_NAME.get(name, name)
        colour = colours.get(sheet_name)
        group = COLOUR_GROUP.get(colour) if colour else None
        if group is None:
            group = next((gid for pat, gid in DESIGN_GROUP if re.search(pat, sheet_name)), None)
        dtype = str(g["DatasetType"].iloc[0]).strip()
        origin, subtype = ORIGIN.get(dtype, ("excluded", "unclassified"))
        rl, layout = read_length(g["Length"].iloc[0])
        ds_rows.append({"benchmark": bench, "benchmark_number": b, "dataset_name": name, "datasets_sheet_name": sheet_name,
                        "dataset_type_reported": dtype, "dataset_origin": origin, "real_subtype": subtype,
                        "read_length_bp": rl, "read_layout": layout, "length_reported": g["Length"].iloc[0],
                        "colour": colour, "shared_source_group": COLOUR_GROUP.get(colour) if colour else None,
                        "bio_dataset_id_raw": group, "is_edgren": group == "EDGREN_BRCA_CELL_LINES"})
    ds = pd.DataFrame(ds_rows)
    ds["dataset_number"] = ds.groupby("benchmark_number").cumcount().add(1).astype(str).radd(ds["benchmark_number"].astype(str) + ".")
    ds["dataset_id"] = "Dataset" + ds["dataset_number"]
    ds["bio_dataset_id"] = ds["bio_dataset_id_raw"].fillna("D_" + ds["dataset_number"])
    ds.drop(columns="bio_dataset_id_raw").to_csv(OUT_DIR / "dataset_annotations.csv", index=False)

    # --- observations --------------------------------------------------------------------------------
    tm = tool_map(perf["Tool"])
    tmap = tm.set_index("tool_reported")
    key = ds.set_index(["benchmark", "dataset_name"])
    legacy_meta = pd.read_excel(LEGACY, sheet_name="benchmark_metadata")
    new_meta = pd.read_csv(CURATION / "benchmark_metadata_new_studies.csv")
    meta = pd.concat([legacy_meta, new_meta[legacy_meta.columns.intersection(new_meta.columns)]], ignore_index=True)
    meta["benchmark_number"] = pd.to_numeric(meta["benchmark_number"]).astype(int)
    legacy_obs = pd.read_excel(LEGACY, sheet_name="observations")
    doi = legacy_obs.groupby("benchmark_number")["doi"].first().to_dict() | DOI
    purpose = meta.set_index("benchmark_number")["benchmark_purpose"].to_dict()
    new_tool = meta.set_index("benchmark_number")["new_tool_annotation"].to_dict()

    obs = []
    for i, r in enumerate(perf.itertuples(index=False), start=1):
        d = key.loc[(r.Benchmark, r.Dataset)]
        p, rc, f1, truth = num(r.Precision), num(r.Recall), num(r.F1), num(r.TrueFusionNumber)
        tp = rc * truth if np.isfinite(rc) and np.isfinite(truth) else np.nan
        fp = tp * (1 - p) / p if np.isfinite(tp) and np.isfinite(p) and p > 0 else np.nan
        reason = next((why for rule, why in EXCLUDE if rule(r)), None)
        if reason is None and d.dataset_origin == "excluded":
            reason = f"dataset type not classified: {d.dataset_type_reported}"
        b = int(d.benchmark_number)
        obs.append({
            "observation_id": f"OBS{i:04d}", "benchmark_number": b, "doi": doi.get(b), "dataset_id": d.dataset_id,
            "dataset_number": d.dataset_number, "dataset_type_reported": d.dataset_type_reported, "dataset_origin": d.dataset_origin,
            "run": r.Notes, "tool_reported": r.Tool, "tool_name_clean": tmap.at[r.Tool, "tool_name_clean"],
            "tool_family": tmap.at[r.Tool, "tool_family"], "tool_version": r.ToolVersion,
            "TP": tp, "FP": fp, "TP_plus_FP": tp + fp if np.isfinite(tp) and np.isfinite(fp) else np.nan, "truth_total": truth,
            "recall": rc, "precision": p, "f1_score": f1, "metric_status": None if np.isfinite(f1) else (str(r.F1) if pd.notna(r.F1) else None),
            "truthset": r.Truthset, "read_length_bp": d.read_length_bp, "read_layout": d.read_layout,
            "has_f1": bool(np.isfinite(f1)), "has_precision_recall": bool(np.isfinite(p) and np.isfinite(rc)),
            "is_edgren": bool(d.is_edgren), "benchmark_purpose": purpose.get(b), "new_tool_annotation": new_tool.get(b),
            "counts_back_calculated": bool(np.isfinite(tp)), "include_in_models": reason is None, "exclusion_reason": reason,
        })
    obs = pd.DataFrame(obs)

    # benchmark-level counts recomputed from the extraction
    counts = obs.groupby("benchmark_number").agg(n_tools=("tool_name_clean", "nunique"), n_datasets=("dataset_id", "nunique"))
    meta["number of tools"] = meta["benchmark_number"].map(counts["n_tools"])
    meta["number of datasets"] = meta["benchmark_number"].map(counts["n_datasets"])

    dsm = ds.rename(columns={"dataset_type_reported": "dataset_type_reported"})[
        ["dataset_id", "dataset_number", "benchmark_number", "dataset_name", "dataset_type_reported", "dataset_origin",
         "real_subtype", "read_length_bp", "read_layout", "is_edgren", "bio_dataset_id", "shared_source_group"]]
    readme = pd.DataFrame({"field": ["source", "built_by", "primary_observation_unit", "excluded_rows"],
                           "value": ["data/raw/benchmark.xlsx", "scripts/build_analysis_input.py", "tool-dataset-condition observation",
                                     "include_in_models = False, reason in exclusion_reason"]})
    with pd.ExcelWriter(OUT) as xw:
        readme.to_excel(xw, sheet_name="README", index=False)
        obs.to_excel(xw, sheet_name="observations", index=False)
        meta.to_excel(xw, sheet_name="benchmark_metadata", index=False)
        dsm.to_excel(xw, sheet_name="dataset_metadata", index=False)
        tm.to_excel(xw, sheet_name="tool_name_map", index=False)
        pd.read_excel(LEGACY, sheet_name="edgren_truthsets").to_excel(xw, sheet_name="edgren_truthsets", index=False)

    print(f"observations: {len(obs)} ({obs.include_in_models.sum()} in models); benchmarks: {obs.benchmark_number.nunique()}; "
          f"datasets: {len(ds)}; biological datasets: {ds.bio_dataset_id.nunique()}; tools: {tm.tool_name_clean.nunique()}")
    print(obs[~obs.include_in_models].exclusion_reason.value_counts().to_string())
    print(ds.groupby("dataset_origin").size().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
