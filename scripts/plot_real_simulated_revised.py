"""Revised Figure 5: real versus simulated benchmark settings under crossed dependence.

Panel A: distributions of F1, precision and recall by data origin, annotated with the primary
         cross-classified model (difference on the metric scale, benchmark-bootstrap 95% CI,
         number of benchmarks and biological datasets).
Panel B: one point per benchmark, mean(real) - mean(simulated) within the same benchmark.
Requires crossed_dependence_models.py to have been run.
"""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".mplconfig"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from revision_utils import METRICS, REVISION_DIR, load_observations

IN = REVISION_DIR / "crossed_models"
OUT = REVISION_DIR / "figures"
COLORS = {"simulated": "#4C78A8", "real": "#D65F5F"}
METRIC_COLORS = {"F1": "#333333", "Precision": "#E08E0B", "Recall": "#5B8E7D"}
LABELS = {"F1": "F1-score", "Precision": "Precision", "Recall": "Recall"}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    df = load_observations()
    models = pd.read_csv(IN / "models.csv")
    within = pd.read_csv(IN / "within_benchmark_contrasts.csv")

    fig = plt.figure(figsize=(14, 9.2))
    grid = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.9], hspace=0.42, wspace=0.25)

    for i, metric in enumerate(METRICS):
        ax = fig.add_subplot(grid[0, i])
        d = df.dropna(subset=[metric])
        groups = [d.loc[d.real == 0, metric], d.loc[d.real == 1, metric]]
        box = ax.boxplot(groups, positions=[1, 2], widths=0.45, patch_artist=True, showfliers=False,
                         medianprops={"color": "black", "linewidth": 1.4})
        for patch, key in zip(box["boxes"], ["simulated", "real"]):
            patch.set_facecolor(COLORS[key])
            patch.set_alpha(0.35)
        rng = np.random.default_rng(7 + i)
        for pos, (g, key) in enumerate(zip(groups, ["simulated", "real"]), start=1):
            ax.scatter(pos + rng.uniform(-0.13, 0.13, len(g)), g, s=9, alpha=0.45, color=COLORS[key], linewidths=0)
        r = models[(models.metric == metric) & (models.scale == "metric difference") & models.ci_method.str.startswith("benchmark")].iloc[0]
        ax.text(0.03, 0.03,
                f"Crossed model: {r.estimate:+.2f}\n95% CI {r.ci_low:+.2f} to {r.ci_high:+.2f}\n"
                f"{int(r.n_benchmarks)} benchmarks, {int(r.n_bio_datasets)} datasets",
                transform=ax.transAxes, fontsize=9, va="bottom",
                bbox={"facecolor": "white", "edgecolor": "#CCCCCC", "boxstyle": "round,pad=0.3"})
        ax.set_xticks([1, 2], ["Simulated", "Real"])
        ax.set_ylim(-0.02, 1.02)
        ax.set_title(LABELS[metric], fontsize=12, fontweight="bold")
        if i == 0:
            ax.set_ylabel("Reported value")
        ax.spines[["top", "right"]].set_visible(False)

    ax = fig.add_subplot(grid[1, :])
    order = (within[within.metric == "F1"].sort_values("diff_real_minus_sim")["benchmark_label"].tolist())
    order += [b for b in within["benchmark_label"].unique() if b not in order]
    y = {b: k for k, b in enumerate(order)}
    offsets = {"F1": 0.0, "Precision": -0.22, "Recall": 0.22}
    for metric in METRICS:
        w = within[within.metric == metric]
        ax.scatter(w["diff_real_minus_sim"], [y[b] + offsets[metric] for b in w["benchmark_label"]],
                   color=METRIC_COLORS[metric], s=46, label=LABELS[metric], zorder=3)
    ax.axvline(0, color="#888888", linestyle="--", linewidth=1)
    ax.set_yticks(range(len(order)), order)
    ax.set_xlabel("Within-benchmark difference, mean(real) - mean(simulated)")
    ax.set_title("B. One value per benchmark (benchmarks reporting both real and simulated data)", fontsize=12,
                 fontweight="bold", loc="left")
    ax.legend(frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.text(0.01, 0.935, "A. Distributions by data origin (each point is one tool score; points are not independent)",
             fontsize=12, fontweight="bold")
    fig.savefig(OUT / "real_vs_simulated_revised.pdf", bbox_inches="tight")
    fig.savefig(OUT / "real_vs_simulated_revised.png", dpi=300, bbox_inches="tight")
    print(f"Saved: {OUT / 'real_vs_simulated_revised.pdf'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
