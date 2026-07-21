"""
Generate the Edgren truth-set drift table and figure.

This figure compares Dataset7.6 and Dataset7.7, which represent the same
Edgren RNA-seq data evaluated against two different truth-set definitions
(27 versus 99 reference fusions).
"""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[2] / ".mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(__file__).resolve().parents[2] / ".cache"))

import matplotlib.pyplot as plt
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
INPUT_CSV = REPO_ROOT / "data" / "processed" / "benchmarks_metrics_enriched.csv"
FIG_DIR = REPO_ROOT / "results" / "figures"
TAB_DIR = REPO_ROOT / "results" / "tables"


DISPLAY_NAMES = {
    "pizzly": "Pizzly",
}

LABEL_OFFSETS = {
    "ChimeraScan": (12, -4),
    "Pizzly": (10, 14),
    "STARChip": (12, -18),
    "EricScript": (8, -2),
    "JAFFA": (8, 8),
    "SOAPfuse": (8, 0),
    "ChimPipe": (6, 10),
    "INTEGRATE": (8, -2),
    "FuSeq": (8, 5),
    "InFusion": (8, 6),
    "TopHat-Fusion": (8, -2),
    "MapSplice": (8, 8),
    "FusionCatcher": (8, 8),
    "STAR-Fusion": (8, 8),
    "Arriba": (8, 8),
}


def compute_drift() -> pd.DataFrame:
    df = pd.read_csv(INPUT_CSV)
    cols = ["Tool", "TP", "FP", "FN_recovered", "Precision", "Recall/Sensitivity", "F1"]
    d27 = df.loc[df["Dataset"].eq("Dataset7.6"), cols].copy()
    d99 = df.loc[df["Dataset"].eq("Dataset7.7"), cols].copy()
    d27.columns = ["Tool", "TP_27", "FP_27", "FN_27", "precision_27", "recall_27", "f1_27"]
    d99.columns = ["Tool", "TP_99", "FP_99", "FN_99", "precision_99", "recall_99", "f1_99"]
    drift = d27.merge(d99, on="Tool", how="inner")
    drift["delta_f1"] = drift["f1_99"] - drift["f1_27"]
    drift["delta_precision"] = drift["precision_99"] - drift["precision_27"]
    drift["delta_recall"] = drift["recall_99"] - drift["recall_27"]
    drift["Tool_display"] = drift["Tool"].replace(DISPLAY_NAMES)
    return drift.sort_values("delta_f1")


def plot_drift(drift: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 6.6))
    upper = max(drift[["f1_27", "f1_99"]].max().max() + 0.07, 0.70)
    ax.plot([0, upper], [0, upper], color="#777777", linestyle="--", linewidth=1.0)
    colors = drift["delta_f1"].map(lambda value: "#D55E00" if value < 0 else "#0072B2")
    ax.scatter(
        drift["f1_27"],
        drift["f1_99"],
        s=62,
        color=colors,
        edgecolor="black",
        linewidth=0.55,
        alpha=0.82,
        zorder=3,
    )

    for _, row in drift.iterrows():
        offset = LABEL_OFFSETS.get(row["Tool_display"], (5, 4))
        ax.annotate(
            row["Tool_display"],
            (row["f1_27"], row["f1_99"]),
            xytext=offset,
            textcoords="offset points",
            fontsize=7.2,
            alpha=0.9,
        )

    ax.set_xlim(0, upper)
    ax.set_ylim(0, upper)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("F1 using 27-fusion truth set")
    ax.set_ylabel("F1 using 99-fusion truth set")
    ax.set_title("Truth-set-induced F1 changes on identical Edgren RNA-seq data", fontsize=12)
    ax.grid(color="#E6E6E6", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(
        0.03,
        upper - 0.06,
        "Dashed line: no change",
        fontsize=8,
        color="#555555",
    )
    fig.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / "edgren_drift_scatter.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / "edgren_drift_scatter.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    TAB_DIR.mkdir(parents=True, exist_ok=True)
    drift = compute_drift()
    drift.to_csv(TAB_DIR / "edgren_truth_set_drift.csv", index=False)
    plot_drift(drift)
    print(drift[["Tool", "f1_27", "f1_99", "delta_f1", "precision_27", "precision_99", "recall_27", "recall_99"]])
    print(f"\nSaved: {TAB_DIR / 'edgren_truth_set_drift.csv'}")
    print(f"Saved: {FIG_DIR / 'edgren_drift_scatter.pdf'}")
    print(f"Saved: {FIG_DIR / 'edgren_drift_scatter.png'}")


if __name__ == "__main__":
    main()
