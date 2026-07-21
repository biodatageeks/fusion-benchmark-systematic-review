"""
Publication-ready figures for the meta-analysis.

Reads processed CSVs from data/processed/ and pooled estimates from
output/tables/, writes PDFs to output/figures/.

Figures:
  1. edgren_drift_slopegraph.pdf — F1 change on same Edgren data when the
     truth set is changed from 27 fusions (Dataset7.6) to 99 (Dataset7.7).
     This is the "killer feature" figure: same predictions, different
     truth sets, dramatic F1 shifts for many tools.
  2. edgren_drift_scatter.pdf — Same data as scatter (F1_27 vs F1_99),
     with y=x reference line.
  3. tool_ranking_f1.pdf — Pooled F1 per tool with 95% CI, sorted
     descending. Provides the "which tool is best" ranking.
  4. real_vs_simulated_f1.pdf — F1 distribution for real vs simulated
     datasets, aggregated across all tools.
  5. f1_by_read_length.pdf — F1 by read length (50/76/100 bp), simulated
     data only, aggregated across tools.
  6. heterogeneity_I2.pdf — I² per tool per metric (bar plot).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
METRICS_CSV = REPO_ROOT / "data" / "processed" / "benchmarks_metrics_enriched.csv"
DATASETS_CSV = REPO_ROOT / "data" / "processed" / "dataset_details.csv"
POOLED_CSV = REPO_ROOT / "results" / "tables" / "pooled_estimates.csv"
DRIFT_CSV = REPO_ROOT / "results" / "tables" / "edgren_truth_set_drift.csv"
FIG_DIR = REPO_ROOT / "results" / "figures"


def savefig(fig, name: str) -> None:
    """Save a figure to PDF for publication."""
    fig.savefig(FIG_DIR / f"{name}.pdf")
    plt.close(fig)
    print(f"  {name}.pdf")

plt.rcParams.update({
    "figure.dpi": 100,
    "savefig.dpi": 200,
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linestyle": "--",
})


def fig_edgren_drift_slopegraph() -> None:
    """Two vertical F1 axes (Dataset7.6 P27, Dataset7.7 P99), lines connect
    tool values. Slope direction shows whether F1 rose or fell."""
    drift = pd.read_csv(DRIFT_CSV).dropna(subset=["f1_27", "f1_99"])
    drift = drift.sort_values("f1_27", ascending=False).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(8, 7))
    x_left, x_right = 0.0, 1.0

    for _, row in drift.iterrows():
        y_left, y_right = row["f1_27"], row["f1_99"]
        color = "tab:red" if y_right < y_left else "tab:blue"
        ax.plot([x_left, x_right], [y_left, y_right],
                color=color, alpha=0.7, linewidth=1.5, marker="o",
                markersize=4)
        ax.text(x_left - 0.05, y_left, row["Tool"],
                ha="right", va="center", fontsize=8)
        ax.text(x_right + 0.05, y_right,
                f"{row['Tool']} ({row['delta_f1']:+.2f})",
                ha="left", va="center", fontsize=8, color=color)

    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(-0.02, 0.65)
    ax.set_xticks([x_left, x_right])
    ax.set_xticklabels(["Truth set = 27 fusions\n(Dataset 7.6)",
                        "Truth set = 99 fusions\n(Dataset 7.7)"])
    ax.set_ylabel("F1 score")
    ax.set_title("Edgren dataset: F1 score varies by truth-set choice alone\n"
                 "(same RNA-seq data, only the reference truth set differs)",
                 fontsize=10)
    ax.grid(axis="y", alpha=0.3)

    # Legend
    ax.plot([], [], color="tab:red", label="F1 decreased with larger truth set")
    ax.plot([], [], color="tab:blue", label="F1 increased with larger truth set")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)

    fig.tight_layout()
    savefig(fig, "edgren_drift_slopegraph")


def fig_edgren_drift_scatter() -> None:
    """Scatter F1_27 vs F1_99 with y=x reference. Distance from diagonal
    equals magnitude of truth-set effect for that tool."""
    drift = pd.read_csv(DRIFT_CSV).dropna(subset=["f1_27", "f1_99"])

    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    lim = max(drift[["f1_27", "f1_99"]].max().max() + 0.05, 0.7)

    ax.plot([0, lim], [0, lim], color="gray", linestyle="--",
            linewidth=1, alpha=0.6, label="y = x (no truth-set effect)")
    ax.scatter(drift["f1_27"], drift["f1_99"],
               s=60, alpha=0.75, edgecolor="black", linewidth=0.5,
               zorder=3)

    for _, row in drift.iterrows():
        ax.annotate(row["Tool"], (row["f1_27"], row["f1_99"]),
                    xytext=(4, 4), textcoords="offset points",
                    fontsize=7.5, alpha=0.85)

    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("F1 with truth set = 27 fusions (Dataset 7.6)")
    ax.set_ylabel("F1 with truth set = 99 fusions (Dataset 7.7)")
    ax.set_title("Truth-set-induced F1 change on identical Edgren data\n"
                 "(16 tools)", fontsize=10)
    ax.set_aspect("equal")
    ax.legend(loc="lower right", fontsize=8)

    fig.tight_layout()
    savefig(fig, "edgren_drift_scatter")


def fig_tool_ranking_f1() -> None:
    """Bar plot: pooled F1 per tool with 95% CI, sorted descending."""
    pooled = pd.read_csv(POOLED_CSV)
    f1 = pooled[pooled["metric"] == "f1"].copy()
    f1 = f1.sort_values("pooled", ascending=True)  # ascending → largest on top

    fig, ax = plt.subplots(figsize=(7, 5))
    y = np.arange(len(f1))
    err_lo = f1["pooled"] - f1["ci_lb"]
    err_hi = f1["ci_ub"] - f1["pooled"]

    ax.errorbar(f1["pooled"], y, xerr=[err_lo, err_hi],
                fmt="o", capsize=3, color="tab:blue",
                markersize=6, elinewidth=1.2)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{t} (k={k})" for t, k in
                        zip(f1["tool"], f1["k"])])
    ax.set_xlabel("Pooled F1 score (95% CI)")
    ax.set_xlim(0, 1)
    ax.axvline(0.5, color="gray", linestyle=":", linewidth=0.8, alpha=0.5)
    ax.set_title("Random-effects pooled F1 across benchmarks\n"
                 "(tools tested in ≥ 3 benchmarks; k = number of data points)",
                 fontsize=10)
    fig.tight_layout()
    savefig(fig, "tool_ranking_f1")


def fig_real_vs_simulated() -> None:
    """Boxplot of F1 for real vs simulated datasets, all tools pooled."""
    df = pd.read_csv(METRICS_CSV)
    df = df.dropna(subset=["F1", "Type of dataset"])

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    groups = df.groupby("Type of dataset")["F1"].apply(list)
    labels = list(groups.index)
    data = [groups[l] for l in labels]

    bp = ax.boxplot(data, tick_labels=labels, patch_artist=True,
                    widths=0.55, showmeans=True,
                    meanprops={"marker": "D", "markersize": 5,
                               "markerfacecolor": "white",
                               "markeredgecolor": "black"})
    colors = ["tab:orange", "tab:green"] * len(labels)
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.5)

    ax.set_ylabel("F1 score")
    ax.set_ylim(0, 1)
    ax.set_title("F1 distribution by dataset type\n"
                 "(all tools, all benchmarks)", fontsize=10)

    for i, arr in enumerate(data, start=1):
        ax.text(i, 1.02, f"n={len(arr)}", ha="center", fontsize=8)

    fig.tight_layout()
    savefig(fig, "real_vs_simulated_f1")


def fig_f1_by_read_length() -> None:
    """F1 by read length (50/76/100 bp) for simulated datasets."""
    metrics = pd.read_csv(METRICS_CSV)
    dsets = pd.read_csv(DATASETS_CSV)
    dsets = dsets.rename(columns={"Dataset number": "Dataset",
                                  "Read lenght": "read_length"})
    merged = metrics.merge(dsets[["Dataset", "read_length"]],
                           on="Dataset", how="left")
    sim = merged[(merged["Type of dataset"] == "simulated")
                 & merged["F1"].notna()
                 & merged["read_length"].notna()].copy()
    sim["read_length"] = sim["read_length"].astype(int)
    lengths = sorted(sim["read_length"].unique())

    fig, ax = plt.subplots(figsize=(6, 4.5))
    data = [sim.loc[sim["read_length"] == L, "F1"].values for L in lengths]
    bp = ax.boxplot(data, tick_labels=[f"{L} bp" for L in lengths],
                    patch_artist=True, widths=0.55, showmeans=True,
                    meanprops={"marker": "D", "markersize": 5,
                               "markerfacecolor": "white",
                               "markeredgecolor": "black"})
    for patch in bp["boxes"]:
        patch.set_facecolor("tab:blue")
        patch.set_alpha(0.4)

    ax.set_ylabel("F1 score")
    ax.set_ylim(0, 1)
    ax.set_xlabel("Read length")
    ax.set_title("F1 by read length (simulated datasets only)", fontsize=10)

    for i, arr in enumerate(data, start=1):
        ax.text(i, 1.02, f"n={len(arr)}", ha="center", fontsize=8)

    fig.tight_layout()
    savefig(fig, "f1_by_read_length")


def fig_heterogeneity() -> None:
    """I² per tool per metric — grouped bar plot."""
    pooled = pd.read_csv(POOLED_CSV)
    tools = pooled[pooled["metric"] == "f1"].sort_values("I2")["tool"].tolist()
    metrics = ["precision", "sensitivity", "f1"]
    colors = {"precision": "tab:orange", "sensitivity": "tab:green",
              "f1": "tab:blue"}

    fig, ax = plt.subplots(figsize=(8, 5))
    n_tools = len(tools)
    bar_width = 0.25
    x = np.arange(n_tools)

    for i, m in enumerate(metrics):
        sub = pooled[pooled["metric"] == m].set_index("tool").reindex(tools)
        ax.bar(x + (i - 1) * bar_width, sub["I2"],
               width=bar_width, label=m, color=colors[m], alpha=0.85)

    ax.axhline(75, color="red", linestyle="--", linewidth=0.8, alpha=0.5,
               label="I² = 75% (Cochrane threshold)")
    ax.set_xticks(x)
    ax.set_xticklabels(tools, rotation=35, ha="right")
    ax.set_ylabel("I² (%)")
    ax.set_ylim(0, 105)
    ax.set_title("Between-benchmark heterogeneity per tool\n"
                 "(I² > 75% ≈ very high; nearly all metrics exceed this)",
                 fontsize=10)
    ax.legend(loc="lower left", fontsize=8)
    fig.tight_layout()
    savefig(fig, "heterogeneity_I2")


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    print("Generating figures:")
    fig_edgren_drift_slopegraph()
    fig_edgren_drift_scatter()
    fig_tool_ranking_f1()
    fig_real_vs_simulated()
    fig_f1_by_read_length()
    fig_heterogeneity()
    print(f"\nAll figures written to {FIG_DIR.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
