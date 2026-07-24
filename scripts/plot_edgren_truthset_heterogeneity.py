from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(__file__).resolve().parents[1] / ".cache"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from real_simulated_sensitivity import DEFAULT_INPUT, DEFAULT_OUTPUT_DIR, load_table, normalize_columns


COLOR_TRUTH = "#5A6FBB"
COLOR_F1 = "#D17C3F"


def jitter_positions(center: float, size: int, width: float = 0.11) -> np.ndarray:
    rng = np.random.default_rng(123 + int(center * 100))
    return center + rng.uniform(-width, width, size=size)


def prepare_edgren_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    df = normalize_columns(load_table(DEFAULT_INPUT))
    if "Is_Edgren_Binary" in df.columns:
        edgren_mask = df["Is_Edgren_Binary"].astype(str).str.lower().str.strip().eq("yes")
    elif "Is_Edgren" in df.columns:
        edgren_mask = df["Is_Edgren"].astype(str).str.lower().str.strip().isin(["yes", "true", "1"])
    else:
        edgren_mask = pd.Series(False, index=df.index)
    edgren = df[edgren_mask].copy()
    edgren = edgren.dropna(subset=["TruthTotal"])
    edgren["Benchmark_label"] = "B" + edgren["Benchmark"].astype(str) + " / " + edgren["Dataset"].astype(str)

    truth_sets = (
        edgren[["Benchmark", "Dataset", "Benchmark_label", "TruthTotal"]]
        .drop_duplicates()
        .sort_values(["Benchmark", "Dataset", "TruthTotal"])
        .reset_index(drop=True)
    )

    f1_data = edgren.dropna(subset=["F1"]).copy()
    tool_counts = f1_data.groupby("Tool").size()
    repeated_tools = tool_counts[tool_counts >= 3].index
    f1_data = f1_data[f1_data["Tool"].isin(repeated_tools)].copy()

    tool_order = (
        f1_data.groupby("Tool")["F1"]
        .median()
        .sort_values(ascending=False)
        .index.tolist()
    )
    f1_data["Tool"] = pd.Categorical(f1_data["Tool"], categories=tool_order, ordered=True)
    return truth_sets, f1_data


def plot_truth_sets(ax: plt.Axes, truth_sets: pd.DataFrame) -> None:
    x = np.arange(len(truth_sets))
    ax.bar(x, truth_sets["TruthTotal"], color=COLOR_TRUTH, alpha=0.82, width=0.68)
    ax.scatter(x, truth_sets["TruthTotal"], color="#27345C", s=28, zorder=3)
    for xpos, value in zip(x, truth_sets["TruthTotal"]):
        ax.text(xpos, value + 2.2, f"{int(value)}", ha="center", va="bottom", fontsize=9)

    ax.set_title("A. Reference-set size", loc="left", fontsize=12, fontweight="bold")
    ax.set_ylabel("True fusions in reference set")
    ax.set_xticks(x)
    ax.set_xticklabels(truth_sets["Benchmark_label"], rotation=35, ha="right")
    ax.set_ylim(0, max(105, truth_sets["TruthTotal"].max() + 15))
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_tool_f1(ax: plt.Axes, f1_data: pd.DataFrame) -> None:
    tools = list(f1_data["Tool"].cat.categories)
    if not tools:
        ax.text(
            0.5,
            0.5,
            "No Edgren-derived tools with at least three F1 observations",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_axis_off()
        return
    positions = np.arange(1, len(tools) + 1)
    values = [f1_data.loc[f1_data["Tool"].eq(tool), "F1"].to_numpy() for tool in tools]

    box = ax.boxplot(
        values,
        positions=positions,
        widths=0.48,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.3},
        whiskerprops={"color": "#444444", "linewidth": 0.9},
        capprops={"color": "#444444", "linewidth": 0.9},
        boxprops={"color": "#444444", "linewidth": 0.9},
    )
    for patch in box["boxes"]:
        patch.set_facecolor(COLOR_F1)
        patch.set_alpha(0.32)

    for xpos, tool in zip(positions, tools):
        tool_values = f1_data.loc[f1_data["Tool"].eq(tool), "F1"].to_numpy()
        ax.scatter(
            jitter_positions(xpos, len(tool_values)),
            tool_values,
            s=28,
            color=COLOR_F1,
            alpha=0.78,
            edgecolor="white",
            linewidth=0.3,
            zorder=3,
        )

    ax.set_title("B. F1-score variation across repeated tools", loc="left", fontsize=12, fontweight="bold")
    ax.set_ylabel("F1-score")
    ax.set_xticks(positions)
    ax.set_xticklabels(tools, rotation=45, ha="right")
    ax.set_ylim(-0.03, 1.03)
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main() -> None:
    output_dir = DEFAULT_OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    truth_sets, f1_data = prepare_edgren_data()
    truth_sets.to_csv(output_dir / "edgren_truthset_figure_truth_sets.csv", index=False)
    f1_data.to_csv(output_dir / "edgren_truthset_figure_f1_data.csv", index=False)

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(12.2, 4.6),
        gridspec_kw={"width_ratios": [1.0, 1.65]},
    )
    plot_truth_sets(axes[0], truth_sets)
    plot_tool_f1(axes[1], f1_data)

    fig.suptitle(
        "Truth-set heterogeneity in benchmarks using the Edgren breast cancer dataset",
        fontsize=13,
        fontweight="bold",
        y=1.02,
    )
    fig.tight_layout()

    png_path = output_dir / "edgren_truthset_heterogeneity.png"
    pdf_path = output_dir / "edgren_truthset_heterogeneity.pdf"
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")


if __name__ == "__main__":
    main()
