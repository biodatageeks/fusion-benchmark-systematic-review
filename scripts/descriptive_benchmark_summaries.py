"""Benchmark-level descriptive summaries that replace inferential tests (R1.1, R1 minor, R2, R3.2, R3.4).

Outputs (results/revision/descriptive/ and results/revision/tables/):
- per_benchmark_summary.csv / supp_benchmark_summary.tex
    one row per benchmark: year, purpose, observations, datasets, mean/median F1 by data origin.
- benchmark_purpose_descriptive.csv
    new-tool vs independent benchmarks at the benchmark level (no p-values; 2 vs 9 studies).
- read_length_descriptive.csv
    F1 by harmonised read length with the number of benchmarks and biological datasets behind each category.
- real_subtype_descriptive.csv
    real data split into cell lines / patient-derived / unspecified.
- context_strata.csv / supp_context_strata.tex
    evidence organised by measurable use context (R3.2).
- edgren_shared_tools_by_purpose.csv
    the same tools on the same biological dataset (Edgren) in new-tool vs independent benchmarks (R3.4).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from revision_utils import BIO_MAP, REVISION_DIR, fmt, load_observations, write_latex_table
from real_simulated_sensitivity import DEFAULT_INPUT, load_table
from revision_utils import apply_metadata_corrections

OUT = REVISION_DIR / "descriptive"
TABLES = REVISION_DIR / "tables"


def harmonise_read_length(value: float) -> str | None:
    if pd.isna(value):
        return None
    for label, low, high in [("50 bp", 45, 55), ("75 bp", 70, 80), ("100 bp", 95, 105)]:
        if low <= value <= high:
            return label
    return f"{int(value)} bp"


def iqr(series: pd.Series) -> str:
    s = series.dropna()
    return "--" if s.empty else f"{s.quantile(0.25):.2f}--{s.quantile(0.75):.2f}"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    df = load_observations()
    meta = apply_metadata_corrections(load_table(DEFAULT_INPUT, sheet_name="benchmark_metadata"))
    meta["Benchmark"] = pd.to_numeric(meta["benchmark_number"], errors="coerce").astype(int).astype(str)
    f1 = df.dropna(subset=["F1"])

    # --- per benchmark ---------------------------------------------------------------------
    rows = []
    for bench, g in df.groupby("Benchmark"):
        gf = g.dropna(subset=["F1"])
        rows.append({
            "benchmark": int(bench),
            "label": g["Benchmark_label"].iloc[0],
            "year": int(g["Benchmark_year"].iloc[0]),
            "purpose": "new tool" if g["Benchmark_purpose"].iloc[0] == "new_tool" else "independent",
            "n_obs": len(g),
            "n_obs_f1": len(gf),
            "n_tools": g["Tool"].nunique(),
            "n_bio_datasets": g["bio_dataset_id"].nunique(),
            "mean_f1_sim": gf.loc[gf.real == 0, "F1"].mean(),
            "mean_f1_real": gf.loc[gf.real == 1, "F1"].mean(),
            "median_f1_all": gf["F1"].median(),
        })
    per = pd.DataFrame(rows).sort_values(["year", "label"])
    per.to_csv(OUT / "per_benchmark_summary.csv", index=False)
    tex = per.assign(
        **{
            "Benchmark": per.label,
            "Purpose": per.purpose,
            "Obs (F1)": per.apply(lambda r: f"{r.n_obs} ({r.n_obs_f1})", axis=1),
            "Tools": per.n_tools,
            "Datasets": per.n_bio_datasets,
            "Mean F1 sim.": per.mean_f1_sim.map(fmt),
            "Mean F1 real": per.mean_f1_real.map(fmt),
            "Median F1": per.median_f1_all.map(fmt),
        }
    )[["Benchmark", "Purpose", "Obs (F1)", "Tools", "Datasets", "Mean F1 sim.", "Mean F1 real", "Median F1"]]
    tex = tex.replace("NA", "--")
    write_latex_table(tex, TABLES / "supp_benchmark_summary.tex",
                      caption=r"\textbf{Per-benchmark summary.} Observations (with F1 available), distinct tools, "
                              r"biological datasets, and mean F1 by data origin, in chronological order. "
                              r"Benchmark purpose is shown descriptively; no inferential test is reported because only two "
                              r"benchmarks are independent comparisons.",
                      label="tab:supp_benchmark_summary", colspec="llrrrrrr")

    # --- benchmark purpose (benchmark-level) -------------------------------------------------
    purpose = per.groupby("purpose").agg(
        n_benchmarks=("benchmark", "size"),
        n_obs_f1=("n_obs_f1", "sum"),
        mean_of_benchmark_median_f1=("median_f1_all", "mean"),
        min_benchmark_median_f1=("median_f1_all", "min"),
        max_benchmark_median_f1=("median_f1_all", "max"),
    ).reset_index()
    obs_level = f1.groupby("Benchmark_purpose")["F1"].agg(["mean", "median", "size"]).reset_index()
    largest = per.sort_values("n_obs_f1", ascending=False).groupby("purpose").head(1)[["purpose", "label", "n_obs_f1"]]
    purpose = purpose.merge(largest.rename(columns={"label": "largest_benchmark", "n_obs_f1": "largest_benchmark_n_obs_f1"}), on="purpose")
    purpose.to_csv(OUT / "benchmark_purpose_descriptive.csv", index=False)
    obs_level.to_csv(OUT / "benchmark_purpose_observation_level.csv", index=False)

    # --- read length -------------------------------------------------------------------------
    rl = f1.assign(read_length=f1["Read_length"].map(harmonise_read_length)).dropna(subset=["read_length"])
    rlt = rl.groupby("read_length").agg(
        n_obs=("F1", "size"), n_benchmarks=("Benchmark", "nunique"), n_bio_datasets=("bio_dataset_id", "nunique"),
        n_tools=("Tool", "nunique"), median_f1=("F1", "median"), mean_f1=("F1", "mean"),
        share_simulated=("real", lambda s: 1 - s.mean()),
        benchmarks=("Benchmark_label", lambda s: "; ".join(sorted(set(s)))),
    ).reset_index()
    rlt["iqr_f1"] = rl.groupby("read_length")["F1"].apply(iqr).values
    rlt.to_csv(OUT / "read_length_descriptive.csv", index=False)

    # --- real-data subtype -------------------------------------------------------------------
    real = df[df.real == 1]
    sub = real.groupby("real_subtype").agg(
        n_obs=("Tool", "size"), n_obs_f1=("F1", lambda s: s.notna().sum()), n_benchmarks=("Benchmark", "nunique"),
        n_bio_datasets=("bio_dataset_id", "nunique"), median_f1=("F1", "median"), median_precision=("Precision", "median"),
        median_recall=("Recall", "median"),
    ).reset_index()
    sub.to_csv(OUT / "real_subtype_descriptive.csv", index=False)

    # --- use contexts (R3.2) -----------------------------------------------------------------
    targeted = set(meta.loc[meta["targeted"].astype(str).str.lower().str.startswith("yes"), "Benchmark"])

    def context(r: pd.Series) -> str:
        if r.real == 0:
            return "Simulated reads"
        if r.bio_dataset_id == "EDGREN_BRCA_CELL_LINES":
            return "Edgren breast cancer cell lines"
        if r.Benchmark in targeted and r.real_subtype == "patient-derived":
            return "Patient-derived, targeted RNA panel"
        if r.real_subtype == "patient-derived":
            return "Patient-derived, whole transcriptome"
        if r.real_subtype == "cell line":
            return "Other cancer cell lines"
        if r.real_subtype == "mixed":
            return "Real, mixed patient-derived and cell-line material"
        return "Real, material not specified"

    ctx = df.assign(context=df.apply(context, axis=1))
    ct = ctx.groupby("context").agg(
        n_benchmarks=("Benchmark", "nunique"), n_bio_datasets=("bio_dataset_id", "nunique"), n_obs=("Tool", "size"),
        n_tools=("Tool", "nunique"), median_precision=("Precision", "median"), median_recall=("Recall", "median"),
        median_f1=("F1", "median"), benchmarks=("Benchmark_label", lambda s: "; ".join(sorted(set(s)))),
    ).reset_index()
    order = ["Simulated reads", "Edgren breast cancer cell lines", "Other cancer cell lines",
             "Patient-derived, whole transcriptome", "Patient-derived, targeted RNA panel",
             "Real, mixed patient-derived and cell-line material", "Real, material not specified"]
    ct["order"] = ct["context"].map({c: i for i, c in enumerate(order)})
    ct = ct.sort_values("order").drop(columns="order")
    ct.to_csv(OUT / "context_strata.csv", index=False)
    ctex = ct.assign(**{
        "Context": ct.context, "Benchmarks": ct.n_benchmarks, "Datasets": ct.n_bio_datasets, "Obs": ct.n_obs,
        "Tools": ct.n_tools, "Median precision": ct.median_precision.map(fmt), "Median recall": ct.median_recall.map(fmt),
        "Median F1": ct.median_f1.map(fmt),
    })[["Context", "Benchmarks", "Datasets", "Obs", "Tools", "Median precision", "Median recall", "Median F1"]].replace("NA", "--")
    write_latex_table(ctex, TABLES / "supp_context_strata.tex",
                      caption=r"\textbf{Benchmark evidence by measurable use context.} Medians are descriptive and pool "
                              r"different tools, truth sets and scoring rules; contexts with few benchmarks do not support "
                              r"tool recommendations. Real-data subtype follows the curated dataset annotation "
                              r"(\texttt{data/curation/bio\_dataset\_map.csv}).",
                      label="tab:supp_context_strata", colspec="lrrrrrrr", resize=True)

    # --- same tools on Edgren under both purposes (R3.4) -------------------------------------
    ed = f1[f1.bio_dataset_id == "EDGREN_BRCA_CELL_LINES"]
    by = ed.groupby(["Tool", "Benchmark_purpose"])["F1"].median().unstack()
    if {"new_tool", "independent_comparison"}.issubset(by.columns):
        shared = by.dropna(subset=["new_tool", "independent_comparison"]).copy()
        shared["diff_new_minus_independent"] = shared["new_tool"] - shared["independent_comparison"]
        shared.reset_index().to_csv(OUT / "edgren_shared_tools_by_purpose.csv", index=False)
        print(f"Edgren tools evaluated under both purposes: {len(shared)}; "
              f"median F1 higher in new-tool benchmarks for {(shared.diff_new_minus_independent > 0).sum()}")

    # --- biological-dataset map (Supplementary table) -----------------------------------------
    bio = pd.read_csv(BIO_MAP, dtype=str)
    labels = df.drop_duplicates("Benchmark").set_index("Benchmark")["Benchmark_label"]
    bio["Benchmark"] = bio["benchmark_number"].map(labels).fillna(bio["benchmark"])
    bio["order"] = pd.to_numeric(bio["benchmark_number"])
    shared = bio.groupby("bio_dataset_id")["benchmark_number"].transform("nunique").gt(1) | bio.groupby("bio_dataset_id")["dataset_id"].transform("size").gt(1)
    btex = bio[shared].sort_values(["bio_dataset_id", "order", "dataset_id"]).rename(columns={
        "dataset_name": "Extracted dataset", "bio_dataset_id": "Biological dataset", "dataset_origin": "Data type"})[
        ["Biological dataset", "Benchmark", "Extracted dataset", "Data type"]]
    write_latex_table(
        btex, TABLES / "supp_bio_dataset_map.tex",
        caption=rf"\textbf{{Biological-dataset identifiers used to model crossed dependence.}} Extracted datasets that use the "
                rf"same biological samples, the same sequencing runs or the same simulated design share one identifier "
                rf"({bio['dataset_id'].nunique()} extracted datasets correspond to {bio['bio_dataset_id'].nunique()} biological datasets; "
                r"only shared identifiers are listed). Groups follow the shared-source annotation of the extraction workbook, "
                r"plus variants of one simulated design. Full list: \texttt{data/processed/dataset\_annotations.csv}.",
        label="tab:supp_bio_dataset_map", colspec="llp{0.42\\linewidth}l", font=r"\scriptsize", longtable=True)

    # --- user checklist applied to each benchmark (R2.1) ----------------------------------------
    def yes(value: object) -> bool:
        return str(value).strip().lower().startswith("yes")

    def status(value: object) -> str:
        text = str(value).strip().lower()
        if "partial" in text or "limited" in text or "unclear" in text or "minimal" in text:
            return "partial"
        return "yes" if text.startswith("yes") or text == "github" else "no"

    rows = []
    for _, m in meta.sort_values(["year_x", "first_author"]).iterrows():
        origin = [lab for col, lab in [("simulated", "simulated"), ("real (cells)", "cell lines"), ("real (patients)", "patients")] if yes(m[col])]
        protocol = "targeted" if yes(m["targeted"]) and not yes(m["transciptome"]) else ("WTS + targeted" if yes(m["targeted"]) else "WTS")
        layout = "/".join(lab for col, lab in [("PE", "PE"), ("SE", "SE")] if yes(m[col]))
        lengths = [lab for col, lab in [("read lenght 50", "50"), ("read lenght 75", "75"), ("read lenght 100", "100"), ("read lenght 250", "250")] if yes(m[col])]
        introduced = str(m["new_tool_annotation"]).strip()
        rows.append({
            "Benchmark": f"{m['first_author']} {int(m['year_x'])}",
            "1 Data origin": f"{', '.join(origin)}; {protocol}; {layout or '?'}; " + (f"{'/'.join(lengths)} bp" if lengths else "read length n.r."),
            "2 Truth set": str(m["ground_truth_type"]).replace("“", "``").replace("”", "''") + ("; experimentally supported" if yes(m["experimentally_validated"]) else ""),
            "3 Granularity": "--",
            "4 Execution": f"versions {status(m['tool_versions_reported'])}; parameters {status(m['parameters_disclosed'])}; code {status(m['code_available'])}",
            "5 Provenance": "independent" if introduced.lower() in {"no", "nan", ""} else "new tool (" + introduced.replace("yes, ", "") + ")",
        })
    write_latex_table(
        pd.DataFrame(rows), TABLES / "supp_checklist_applied.tex",
        caption=r"\textbf{The user checklist (main text Table~2) applied to each included benchmark.} Entries are derived "
                r"from the extracted benchmark metadata; WTS, whole transcriptome; PE/SE, paired/single-end; n.r., not reported. "
                r"\TODO{add scoring unit and matching rule (dimension 3) per benchmark from the extraction}",
        label="tab:supp_checklist_applied", colspec="p{0.11\\linewidth}p{0.2\\linewidth}p{0.2\\linewidth}p{0.08\\linewidth}p{0.16\\linewidth}p{0.13\\linewidth}",
        font=r"\scriptsize")

    print(per.to_string(index=False))
    print(purpose.to_string(index=False))
    print(rlt.drop(columns="benchmarks").to_string(index=False))
    print(sub.to_string(index=False))
    print(ct.drop(columns="benchmarks").to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
