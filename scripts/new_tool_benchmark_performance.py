from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from real_simulated_sensitivity import load_table


DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "raw"
OUTPUT_DIR = Path(__file__).resolve().parents[1] / "results" / "new_tool_performance"
MASTER_DATA = DATA_DIR / "master_analysis_input.xlsx"
METRICS = ["F1", "Precision", "Recall_Sensitivity"]


def clean_tool(value: object) -> str:
    return re.sub(r"\s+", " ", str(value).strip())


def parse_new_tool(value: object) -> list[str]:
    text = str(value).strip()
    if not text or text.lower() == "no" or text.lower() == "nan":
        return []
    text = re.sub(r"^yes,\s*", "", text, flags=re.IGNORECASE)
    return [part.strip() for part in re.split(r"\s+and\s+|,", text) if part.strip()]


def canonical_family(tool: str) -> str:
    normalized = clean_tool(tool).lower().replace("_", "-")
    normalized_hyphenated = normalized.replace(" ", "-")
    if normalized.startswith("trinityfusion"):
        return "trinityfusion"
    if normalized.startswith("jaffa"):
        return "jaffa"
    if normalized.startswith("integrate"):
        return "integrate"
    if normalized.startswith("metafusion"):
        return "metafusion"
    if normalized == "gfusion":
        return "gfusion"
    if normalized == "star-fusion":
        return "star-fusion"
    if normalized == "arriba" or normalized == "arriba-hc":
        return "arriba"
    if normalized == "seekfusion":
        return "seekfusion"
    if normalized == "fuseq":
        return "fuseq"
    if normalized_hyphenated.startswith("fusion-inpipe"):
        return "fusion-inpipe"
    return normalized


def new_tool_families(tool_names: list[str]) -> set[str]:
    return {canonical_family(tool) for tool in tool_names}


def study_label(row: pd.Series) -> str:
    author = row.get("First_author", row.get("first_author", ""))
    year = row.get("Year", row.get("year", row.get("year_x", "")))
    year_text = str(int(year)) if pd.notna(year) and str(year).strip() else ""
    return f"{author} {year_text}".strip()


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    all_df = load_table(MASTER_DATA, sheet_name="observations").dropna(subset=["benchmark_number", "tool_name_clean"]).copy()
    repro = load_table(MASTER_DATA, sheet_name="benchmark_metadata").dropna(subset=["benchmark_number"]).copy()

    all_df = all_df.rename(
        columns={
            "benchmark_number": "Benchmark",
            "tool_name_clean": "Tool",
            "f1_score": "F1",
            "precision": "Precision",
            "recall": "Recall_Sensitivity",
            "dataset_id": "Dataset",
        }
    )
    all_df["Benchmark"] = all_df["Benchmark"].astype(int)
    all_df["Tool_clean"] = all_df["Tool"].map(clean_tool)
    all_df["Tool_family"] = all_df["Tool_clean"].map(canonical_family)
    for metric in METRICS:
        all_df[metric] = pd.to_numeric(all_df[metric], errors="coerce")

    repro = repro.rename(columns={"benchmark_number": "Benchmark_number", "new_tool_annotation": "New_tool?"})
    repro["Benchmark"] = repro["Benchmark_number"].astype(int)
    repro["New_tool_list"] = repro["New_tool?"].map(parse_new_tool)
    repro["New_tool_families"] = repro["New_tool_list"].map(new_tool_families)
    return all_df, repro


def summarize_benchmark(all_df: pd.DataFrame, repro_row: pd.Series) -> dict[str, object]:
    benchmark = int(repro_row["Benchmark"])
    benchmark_df = all_df[all_df["Benchmark"].eq(benchmark)].copy()
    new_families = set(repro_row["New_tool_families"])
    benchmark_df["Is_new_tool"] = benchmark_df["Tool_family"].isin(new_families)

    tool_summary = (
        benchmark_df.groupby(["Tool_clean", "Tool_family", "Is_new_tool"], observed=False)[METRICS]
        .mean()
        .reset_index()
    )
    tool_summary["F1_rank"] = tool_summary["F1"].rank(ascending=False, method="min")
    tool_summary["Precision_rank"] = tool_summary["Precision"].rank(ascending=False, method="min")
    tool_summary["Recall_rank"] = tool_summary["Recall_Sensitivity"].rank(ascending=False, method="min")

    new_tools = tool_summary[tool_summary["Is_new_tool"]].copy()
    other_tools = tool_summary[~tool_summary["Is_new_tool"]].copy()

    if new_tools.empty:
        return {
            "Benchmark": benchmark,
            "Study": study_label(repro_row),
            "New_tool": "; ".join(repro_row["New_tool_list"]),
            "status": "new tool not found in observations",
        }

    ranked_new_tools = new_tools.dropna(subset=["F1"]).sort_values("F1", ascending=False)
    ranked_all_tools = tool_summary.dropna(subset=["F1"]).sort_values("F1", ascending=False)
    if ranked_new_tools.empty or ranked_all_tools.empty:
        recall_ranked_new = new_tools.dropna(subset=["Recall_Sensitivity"]).sort_values("Recall_Sensitivity", ascending=False)
        recall_ranked_all = tool_summary.dropna(subset=["Recall_Sensitivity"]).sort_values("Recall_Sensitivity", ascending=False)
        if recall_ranked_new.empty or recall_ranked_all.empty:
            return {
                "Benchmark": benchmark,
                "Study": study_label(repro_row),
                "New_tool": "; ".join(repro_row["New_tool_list"]),
                "status": "new tool found but F1 unavailable",
            }
        best_new_recall = recall_ranked_new.iloc[0]
        best_overall_recall = recall_ranked_all.iloc[0]
        return {
            "Benchmark": benchmark,
            "Study": study_label(repro_row),
            "New_tool": "; ".join(repro_row["New_tool_list"]),
            "status": "recall only",
            "n_tools": tool_summary.shape[0],
            "best_new_tool": best_new_recall["Tool_clean"],
            "primary_metric": "Recall_Sensitivity",
            "best_new_primary_value": best_new_recall["Recall_Sensitivity"],
            "best_new_primary_rank": int(best_new_recall["Recall_rank"]),
            "best_overall_tool": best_overall_recall["Tool_clean"],
            "best_overall_primary_value": best_overall_recall["Recall_Sensitivity"],
            "new_tool_is_best_by_primary_metric": bool(best_new_recall["Recall_rank"] == 1),
        }

    best_new = ranked_new_tools.iloc[0]
    best_overall = ranked_all_tools.iloc[0]
    n_tools = tool_summary.shape[0]
    new_mean_f1 = new_tools["F1"].mean()
    other_mean_f1 = other_tools["F1"].mean() if not other_tools.empty else np.nan

    dataset_rows = []
    for dataset, dataset_df in benchmark_df.dropna(subset=["F1"]).groupby("Dataset", observed=False):
        dataset_tool_summary = (
            dataset_df.groupby(["Tool_clean", "Tool_family", "Is_new_tool"], observed=False)["F1"]
            .mean()
            .reset_index()
        )
        dataset_tool_summary["rank"] = dataset_tool_summary["F1"].rank(ascending=False, method="min")
        dataset_new = dataset_tool_summary[dataset_tool_summary["Is_new_tool"]]
        if not dataset_new.empty:
            dataset_rows.append(
                {
                    "Dataset": dataset,
                    "best_new_rank": dataset_new["rank"].min(),
                    "best_new_f1": dataset_new["F1"].max(),
                    "best_overall_f1": dataset_tool_summary["F1"].max(),
                    "n_tools": dataset_tool_summary.shape[0],
                    "new_tool_best": bool(dataset_new["rank"].min() == 1),
                }
            )
    dataset_rank_df = pd.DataFrame(dataset_rows)

    return {
        "Benchmark": benchmark,
        "Study": study_label(repro_row),
        "New_tool": "; ".join(repro_row["New_tool_list"]),
        "status": "ok",
        "n_tools": n_tools,
        "best_new_tool": best_new["Tool_clean"],
        "primary_metric": "F1",
        "best_new_primary_value": best_new["F1"],
        "best_new_primary_rank": int(best_new["F1_rank"]),
        "best_new_f1": best_new["F1"],
        "best_new_f1_rank": int(best_new["F1_rank"]),
        "best_overall_tool": best_overall["Tool_clean"],
        "best_overall_f1": best_overall["F1"],
        "new_tool_is_best_by_mean_f1": bool(best_new["F1_rank"] == 1),
        "new_tool_is_best_by_primary_metric": bool(best_new["F1_rank"] == 1),
        "new_tool_mean_f1": new_mean_f1,
        "other_tools_mean_f1": other_mean_f1,
        "new_minus_other_mean_f1": new_mean_f1 - other_mean_f1 if pd.notna(other_mean_f1) else np.nan,
        "best_new_precision_rank": int(best_new["Precision_rank"]) if pd.notna(best_new["Precision_rank"]) else np.nan,
        "best_new_recall_rank": int(best_new["Recall_rank"]) if pd.notna(best_new["Recall_rank"]) else np.nan,
        "datasets_with_new_tool": int(dataset_rank_df.shape[0]) if not dataset_rank_df.empty else 0,
        "datasets_new_tool_best": int(dataset_rank_df["new_tool_best"].sum()) if not dataset_rank_df.empty else 0,
        "median_dataset_rank": dataset_rank_df["best_new_rank"].median() if not dataset_rank_df.empty else np.nan,
    }


def write_report(summary: pd.DataFrame, path: Path) -> None:
    ok = summary[summary["status"].isin(["ok", "recall only"])].copy()
    lines = ["# New-tool performance within original benchmarks", ""]
    lines.append("Benchmark-level comparison based on mean tool performance within each benchmark.")
    lines.append("")
    if ok.empty:
        lines.append("No new tools were found in all_data.")
    else:
        f1_ok = ok[ok["status"].eq("ok")]
        best_count = int(f1_ok["new_tool_is_best_by_mean_f1"].sum()) if not f1_ok.empty else 0
        lines.append(f"New tool/family ranked first by mean F1 in {best_count} of {f1_ok.shape[0]} new-tool benchmarks with available F1.")
        primary_best = int(ok["new_tool_is_best_by_primary_metric"].sum())
        lines.append(f"Using F1 where available and recall otherwise, the new tool/family ranked first in {primary_best} of {ok.shape[0]} new-tool benchmarks.")
        positive = int((f1_ok["new_minus_other_mean_f1"] > 0).sum()) if not f1_ok.empty else 0
        lines.append(f"New tool/family had higher mean F1 than the average of comparator tools in {positive} of {f1_ok.shape[0]} benchmarks with available F1.")
        lines.append("")
        lines.append("| Benchmark | Study | New tool | best new tool | primary metric/rank | best overall | new mean F1 | others mean F1 | datasets best |")
        lines.append("|---:|---|---|---|---:|---|---:|---:|---:|")
        for _, row in ok.sort_values("Benchmark").iterrows():
            new_mean_f1 = row["new_tool_mean_f1"] if "new_tool_mean_f1" in row and pd.notna(row.get("new_tool_mean_f1")) else np.nan
            other_mean_f1 = row["other_tools_mean_f1"] if "other_tools_mean_f1" in row and pd.notna(row.get("other_tools_mean_f1")) else np.nan
            dataset_best = f"{int(row['datasets_new_tool_best'])}/{int(row['datasets_with_new_tool'])}" if pd.notna(row.get("datasets_new_tool_best")) else "NA"
            lines.append(
                f"| {int(row['Benchmark'])} | {row['Study']} | {row['New_tool']} | {row['best_new_tool']} | "
                f"{row['primary_metric']} {int(row['best_new_primary_rank'])}/{int(row['n_tools'])} | {row['best_overall_tool']} | "
                f"{new_mean_f1:.3f} | {other_mean_f1:.3f} | {dataset_best} |"
            )
    missing = summary[~summary["status"].eq("ok")]
    if not missing.empty:
        lines.append("")
        lines.append("## Notes")
        for _, row in missing.iterrows():
            lines.append(f"- Benchmark {int(row['Benchmark'])}: {row['status']} ({row['New_tool']})")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    all_df, repro = load_inputs()

    repro_new = repro[repro["New_tool_list"].map(bool)].copy()
    summary_rows = [summarize_benchmark(all_df, row) for _, row in repro_new.iterrows()]
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUTPUT_DIR / "new_tool_benchmark_summary.csv", index=False)

    tool_rows = []
    for _, repro_row in repro_new.iterrows():
        benchmark = int(repro_row["Benchmark"])
        benchmark_df = all_df[all_df["Benchmark"].eq(benchmark)].copy()
        families = set(repro_row["New_tool_families"])
        benchmark_df["Is_new_tool"] = benchmark_df["Tool_family"].isin(families)
        tool_summary = (
            benchmark_df.groupby(["Benchmark", "Tool_clean", "Tool_family", "Is_new_tool"], observed=False)[METRICS]
            .mean()
            .reset_index()
        )
        tool_summary["F1_rank"] = tool_summary.groupby("Benchmark")["F1"].rank(ascending=False, method="min")
        tool_rows.append(tool_summary)
    pd.concat(tool_rows, ignore_index=True).to_csv(OUTPUT_DIR / "tool_level_mean_metrics_by_benchmark.csv", index=False)

    write_report(summary, OUTPUT_DIR / "summary.md")
    print((OUTPUT_DIR / "summary.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
