from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(__file__).resolve().parents[1] / ".cache"))

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import colormaps
from matplotlib.colors import Normalize

from real_simulated_sensitivity import load_table


INPUT = Path(__file__).resolve().parents[1] / "data" / "raw" / "heatmap.xlsx"
OUTPUT_DIR = Path(__file__).resolve().parents[1] / "results" / "benchmark_statistics"

METRIC_COLUMNS = ["number of tools", "number of datasets"]
FEATURE_COLUMNS = [
    "tool version",
    "transciptome",
    "targeted",
    "simulated",
    "real (cells)",
    "real (patients)",
    "PE",
    "SE",
    "read lenght 50",
    "read lenght 75",
    "read lenght 100",
    "read lenght 250",
    "time",
    "memory",
]

COLUMN_LABELS = {
    "number of tools": "Tools",
    "number of datasets": "Datasets",
    "tool version": "Tool version",
    "transciptome": "Transcriptome",
    "targeted": "Targeted",
    "simulated": "Simulated",
    "real (cells)": "Real cells",
    "real (patients)": "Real patients",
    "PE": "PE",
    "SE": "SE",
    "read lenght 50": "50 bp",
    "read lenght 75": "75 bp",
    "read lenght 100": "100 bp",
    "read lenght 250": "250 bp",
    "time": "Runtime",
    "memory": "Memory",
}

STATUS_COLORS = {
    "full": "#1F5A99",
    "partial": "#74A9CF",
    "absent": "#D6EAF6",
    "unknown": "#FFFFFF",
}
STATUS_LABELS = {
    "full": "Full/included",
    "partial": "Partial",
    "absent": "Not included",
    "unknown": "Not reported",
}
STATUS_TEXT = {
    "full": "Y",
    "partial": "P",
    "absent": "N",
    "unknown": "?",
}


def normalize_value(value: object) -> str:
    text = str(value).strip().lower() if pd.notna(value) else ""
    if text in {"12", "12.0"}:
        return "full"
    if text in {"8", "8.0"}:
        return "partial"
    if text in {"3", "3.0"}:
        return "absent"
    return "unknown"


def prepare_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = load_table(INPUT).dropna(how="all").copy()
    raw["Benchmark"] = [
        f"B{index + 1} ({int(year)})" if pd.notna(year) else f"B{index + 1}"
        for index, year in enumerate(raw["year"])
    ]
    raw = raw.set_index("Benchmark")

    display = pd.DataFrame(index=raw.index)
    for column in METRIC_COLUMNS:
        display[COLUMN_LABELS[column]] = pd.to_numeric(raw[column], errors="coerce")
    for column in FEATURE_COLUMNS:
        display[COLUMN_LABELS[column]] = raw[column].apply(normalize_value)
    return raw, display


def draw_heatmap(display: pd.DataFrame) -> plt.Figure:
    columns = list(display.columns)
    metric_labels = [COLUMN_LABELS[column] for column in METRIC_COLUMNS]
    metric_values = display[metric_labels].astype(float)
    metric_norm = Normalize(vmin=metric_values.min().min(), vmax=metric_values.max().max())
    metric_cmap = colormaps["Blues"]

    fig, ax = plt.subplots(figsize=(15.8, 6.6))
    ax.set_xlim(0, len(columns))
    ax.set_ylim(0, len(display.index))
    ax.invert_yaxis()

    for row_index, benchmark in enumerate(display.index):
        for col_index, column in enumerate(columns):
            value = display.loc[benchmark, column]
            if column in metric_labels:
                color = metric_cmap(0.25 + 0.68 * metric_norm(float(value)))
                label = str(int(value))
                text_color = "white" if metric_norm(float(value)) > 0.55 else "#1A1A1A"
            else:
                color = STATUS_COLORS[value]
                label = STATUS_TEXT[value]
                text_color = "white" if value == "full" else "#333333"

            rect = mpatches.Rectangle(
                (col_index, row_index),
                1,
                1,
                facecolor=color,
                edgecolor="white",
                linewidth=1.1,
            )
            ax.add_patch(rect)
            ax.text(
                col_index + 0.5,
                row_index + 0.5,
                label,
                ha="center",
                va="center",
                fontsize=10,
                fontweight="bold",
                color=text_color,
            )

    ax.set_xticks([i + 0.5 for i in range(len(columns))])
    ax.set_xticklabels(columns, rotation=35, ha="right", fontsize=10.5)
    ax.set_yticks([i + 0.5 for i in range(len(display.index))])
    ax.set_yticklabels(display.index, fontsize=11)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_title("Benchmark statistics and dataset features", fontsize=16, fontweight="bold", pad=18)

    status_handles = [
        mpatches.Patch(facecolor=STATUS_COLORS[key], edgecolor="#DDDDDD", label=label)
        for key, label in STATUS_LABELS.items()
    ]
    metric_handles = [
        mpatches.Patch(facecolor=metric_cmap(0.30), label="Fewer"),
        mpatches.Patch(facecolor=metric_cmap(0.93), label="More"),
    ]

    first_legend = ax.legend(
        handles=status_handles,
        title="Feature status",
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        frameon=False,
        fontsize=10,
        title_fontsize=11,
    )
    ax.add_artist(first_legend)
    ax.legend(
        handles=metric_handles,
        title="Tools / datasets",
        loc="upper left",
        bbox_to_anchor=(1.01, 0.55),
        frameon=False,
        fontsize=10,
        title_fontsize=11,
    )
    fig.tight_layout()
    return fig


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    raw, display = prepare_data()
    raw.to_csv(OUTPUT_DIR / "benchmark_statistics_heatmap_raw.csv")
    display.to_csv(OUTPUT_DIR / "benchmark_statistics_heatmap_display.csv")

    fig = draw_heatmap(display)
    png_path = OUTPUT_DIR / "benchmark_statistics_heatmap.png"
    pdf_path = OUTPUT_DIR / "benchmark_statistics_heatmap.pdf"
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")


if __name__ == "__main__":
    main()
