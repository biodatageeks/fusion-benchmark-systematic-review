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


METRICS = ["F1", "Precision", "Recall"]
LABELS = {
    "F1": "F1-score",
    "Precision": "Precision",
    "Recall": "Recall",
}
COLORS = {
    "simulated": "#4C78A8",
    "real": "#D65F5F",
}


def jitter_positions(center: float, size: int, width: float = 0.12) -> np.ndarray:
    rng = np.random.default_rng(42 + int(center * 100))
    return center + rng.uniform(-width, width, size=size)


def load_model_results(output_dir: Path) -> pd.DataFrame:
    path = output_dir / "dataset_type_models.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path}. Run real_simulated_sensitivity.py before plotting."
        )
    return pd.read_csv(path).set_index("metric")


def add_metric_panel(ax: plt.Axes, df: pd.DataFrame, models: pd.DataFrame, metric: str) -> None:
    metric_df = df.dropna(subset=[metric, "DatasetType"]).copy()
    groups = [
        metric_df.loc[metric_df["DatasetType"] == "simulated", metric],
        metric_df.loc[metric_df["DatasetType"] == "real", metric],
    ]

    box = ax.boxplot(
        groups,
        positions=[1, 2],
        widths=0.45,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.4},
        whiskerprops={"color": "#444444", "linewidth": 1.0},
        capprops={"color": "#444444", "linewidth": 1.0},
        boxprops={"color": "#444444", "linewidth": 1.0},
    )
    for patch, dataset_type in zip(box["boxes"], ["simulated", "real"]):
        patch.set_facecolor(COLORS[dataset_type])
        patch.set_alpha(0.35)

    for x, dataset_type in zip([1, 2], ["simulated", "real"]):
        values = metric_df.loc[metric_df["DatasetType"] == dataset_type, metric].to_numpy()
        ax.scatter(
            jitter_positions(x, len(values)),
            values,
            s=14,
            color=COLORS[dataset_type],
            alpha=0.55,
            edgecolor="none",
            rasterized=True,
        )

    result = models.loc[metric]
    annotation = (
        rf"$\beta$ = {result['beta_real_vs_simulated']:.3f}"
        "\n"
        rf"95% CI {result['ci_low']:.3f} to {result['ci_high']:.3f}"
        "\n"
        rf"$p$ = {result['p_value']:.2e}"
    )
    ax.text(
        0.04,
        0.05,
        annotation,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.5,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#DDDDDD", "alpha": 0.9},
    )

    ax.set_title(LABELS[metric], fontsize=12, fontweight="bold")
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Simulated", "Real"])
    ax.set_ylim(-0.03, 1.03)
    ax.set_ylabel("Metric value")
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main() -> None:
    output_dir = DEFAULT_OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    df = normalize_columns(load_table(DEFAULT_INPUT))
    models = load_model_results(output_dir)

    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.8), sharey=True)
    for ax, metric in zip(axes, METRICS):
        add_metric_panel(ax, df, models, metric)

    fig.suptitle(
        "Performance estimates are lower for real than simulated RNA-seq datasets",
        fontsize=13,
        fontweight="bold",
        y=1.03,
    )
    fig.tight_layout()

    png_path = output_dir / "real_vs_simulated_f1_precision_recall.png"
    pdf_path = output_dir / "real_vs_simulated_f1_precision_recall.pdf"
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")


if __name__ == "__main__":
    main()
