from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(__file__).resolve().parents[1] / ".cache"))

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import ListedColormap

from real_simulated_sensitivity import load_table


INPUT = Path(__file__).resolve().parents[1] / "data" / "raw" / "master_analysis_input.xlsx"
OUTPUT_DIR = Path(__file__).resolve().parents[1] / "results" / "benchmark_reporting"


TRUTH_COLORS = {
    "No defined truth set / wisdom of crowds": "#8C6BB1",
    "Proxy-based": "#D95F02",
    "Published validated fusions": "#1B9E77",
    "RNA + WGS + experimental validation": "#386CB0",
    "Clinically validated/expected fusions": "#E7298A",
    "Unknown": "#BDBDBD",
}

STATUS_COLORS = {
    "no": "#F2F2F2",
    "partial": "#F6C85F",
    "yes": "#4C9F70",
}


def normalize_text(value: object) -> str:
    return str(value).strip() if pd.notna(value) else ""


def truth_type(value: object) -> str:
    text = normalize_text(value)
    text_lower = text.lower()
    if "clinically validated" in text_lower or "clinically expected" in text_lower:
        return "Clinically validated/expected fusions"
    if "wisdom" in text or "No defined" in text:
        return "No defined truth set / wisdom of crowds"
    if "Proxy" in text:
        return "Proxy-based"
    if "Published" in text:
        return "Published validated fusions"
    if "RNA + WGS" in text:
        return "RNA + WGS + experimental validation"
    return text or "Unknown"


def status(value: object) -> str:
    text = normalize_text(value).lower()
    if not text or text == "nan":
        return "no"
    if "partial" in text or "limited" in text or "unclear" in text or "tool only" in text:
        return "partial"
    if "no or minimal" in text:
        return "partial"
    if text in {"yes", "github"}:
        return "yes"
    if text.startswith("yes"):
        return "yes"
    if text.startswith("no"):
        return "no"
    return "partial"


def new_tool_status(value: object) -> str:
    text = normalize_text(value).lower()
    return "no" if text == "no" or not text else "yes"


def reproducibility_status(row: pd.Series) -> str:
    fully_reproducible = status(row["fully_reproducible"])
    workflow = status(row["workflow_cwl_snakemake"])
    container = status(row["docker_singularity"])
    if fully_reproducible == "yes":
        return "yes"
    if "yes" in {workflow, container} or "partial" in {workflow, container, fully_reproducible}:
        return "partial"
    return "no"


def prepare_heatmap_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = load_table(INPUT, sheet_name="benchmark_metadata").dropna(subset=["benchmark_number"]).copy()
    raw["Benchmark_number"] = pd.to_numeric(raw["benchmark_number"], errors="coerce").astype(int)
    raw["Year"] = pd.to_numeric(raw["year_x"], errors="coerce").astype(int)
    raw["First_author"] = raw["first_author"]
    raw["is_new_borderline_benchmark"] = raw["Benchmark_number"].eq(11)
    raw = raw.sort_values(
        ["is_new_borderline_benchmark", "Year", "Benchmark_number"],
        ascending=[True, False, True],
    )
    raw["Benchmark"] = raw["First_author"].astype(str) + " " + raw["Year"].astype(int).astype(str)

    display = pd.DataFrame(
        {
            "Truth-set type": raw["ground_truth_type"].apply(truth_type),
            "Experimental validation": raw["experimentally_validated"].apply(status),
            "Benchmark code": raw["code_available"].apply(status),
            "Tool versions": raw["tool_versions_reported"].apply(status),
            "Parameters": raw["parameters_disclosed"].apply(status),
            "Fully reproducible": raw.apply(reproducibility_status, axis=1),
            "New-tool paper": raw["new_tool_annotation"].apply(new_tool_status),
        },
    )
    display.index = raw["Benchmark"].tolist()

    codes = pd.DataFrame(index=display.index)
    truth_order = list(TRUTH_COLORS)
    for column in display.columns:
        if column == "Truth-set type":
            codes[column] = display[column].fillna("Unknown").map({name: i for i, name in enumerate(truth_order)})
        else:
            codes[column] = display[column].map({"no": 0, "partial": 1, "yes": 2})
    return display, codes


def draw_heatmap(display: pd.DataFrame, codes: pd.DataFrame) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(15.8, 6.2))
    ax.set_xlim(0, len(display.columns))
    ax.set_ylim(0, len(display.index))
    ax.invert_yaxis()

    for row_index, benchmark in enumerate(display.index):
        for col_index, column in enumerate(display.columns):
            value = display.loc[benchmark, column]
            if column == "Truth-set type":
                value = value if pd.notna(value) else "Unknown"
                color = TRUTH_COLORS[value]
                label = ""
            else:
                value = value if pd.notna(value) else "no"
                color = STATUS_COLORS[value]
                label = {"yes": "Y", "partial": "P", "no": "N"}[value]
            rect = mpatches.Rectangle(
                (col_index, row_index),
                1,
                1,
                facecolor=color,
                edgecolor="white",
                linewidth=1.2,
            )
            ax.add_patch(rect)
            if label:
                text_color = "white" if value == "yes" else "#333333"
                ax.text(
                    col_index + 0.5,
                    row_index + 0.5,
                    label,
                    ha="center",
                    va="center",
                    fontsize=10.5,
                    fontweight="bold",
                    color=text_color,
                )

    ax.set_xticks([i + 0.5 for i in range(len(display.columns))])
    ax.set_xticklabels(display.columns, rotation=35, ha="right", fontsize=10.5)
    ax.set_yticks([i + 0.5 for i in range(len(display.index))])
    ax.set_yticklabels(display.index, fontsize=10.5)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_title(
        "Benchmark reporting and truth-set characteristics",
        fontsize=16,
        fontweight="bold",
        pad=18,
    )

    truth_handles = [
        mpatches.Patch(color=color, label=label)
        for label, color in TRUTH_COLORS.items()
    ]
    status_handles = [
        mpatches.Patch(color=STATUS_COLORS["yes"], label="Yes"),
        mpatches.Patch(color=STATUS_COLORS["partial"], label="Partial / limited"),
        mpatches.Patch(color=STATUS_COLORS["no"], label="No"),
    ]
    first_legend = ax.legend(
        handles=truth_handles,
        title="Truth-set type",
        loc="upper right",
        bbox_to_anchor=(0.995, 0.90),
        bbox_transform=fig.transFigure,
        frameon=False,
        fontsize=9.5,
        title_fontsize=10.5,
    )
    ax.add_artist(first_legend)
    ax.legend(
        handles=status_handles,
        title="Reporting status",
        loc="upper right",
        bbox_to_anchor=(0.995, 0.52),
        bbox_transform=fig.transFigure,
        frameon=False,
        fontsize=9.5,
        title_fontsize=10.5,
    )
    fig.subplots_adjust(left=0.11, right=0.73, bottom=0.24, top=0.88)
    return fig


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    display, codes = prepare_heatmap_data()
    display.to_csv(OUTPUT_DIR / "benchmark_reporting_heatmap_display.csv")
    codes.to_csv(OUTPUT_DIR / "benchmark_reporting_heatmap_codes.csv")

    fig = draw_heatmap(display, codes)
    png_path = OUTPUT_DIR / "benchmark_reporting_truthset_heatmap.png"
    pdf_path = OUTPUT_DIR / "benchmark_reporting_truthset_heatmap.pdf"
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")


if __name__ == "__main__":
    main()
