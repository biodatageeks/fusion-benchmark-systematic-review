from __future__ import annotations

import os
import re
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(__file__).resolve().parents[1] / ".cache"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from real_simulated_sensitivity import load_table


DATA_PATH = Path(__file__).resolve().parents[1] / "data" / "raw" / "all_data.xlsx"
OUTPUT_DIR = Path(__file__).resolve().parents[1] / "results" / "precision_recall_tradeoff"
EDGREN_DATASETS = {"2.3", "5.1", "6.1", "6.2", "7.6", "7.7", "8.3", "10.2"}


def clean_tool(value: object) -> str:
    text = re.sub(r"\s+", " ", str(value).strip())
    return text


def canonical_tool(value: object) -> str:
    text = clean_tool(value)
    normalized = text.lower().replace("_", "-")
    if normalized == "soapfuse":
        return "SOAPfuse"
    if normalized == "tophat-fusion":
        return "TopHat-Fusion"
    if normalized == "star-fusion":
        return "STAR-Fusion"
    if normalized == "fusioncatcher":
        return "FusionCatcher"
    if normalized == "chimerascan":
        return "ChimeraScan"
    if normalized == "pizzly":
        return "Pizzly"
    if normalized == "starchip":
        return "STARChip"
    if normalized == "arriba-hc":
        return "Arriba_hc"
    if normalized.startswith("integrate "):
        return text.replace(" ", "")
    return text


def load_clean_data() -> pd.DataFrame:
    df = load_table(DATA_PATH).copy()
    df["Tool"] = df["Tool"].map(canonical_tool)
    df["DatasetType"] = df["Type_of_dataset"].astype(str).str.lower().str.strip()
    df["DatasetID"] = (
        df["Dataset"]
        .astype(str)
        .str.replace("Dataset", "", regex=False)
        .str.strip()
    )
    df["Is_Edgren"] = df["DatasetID"].isin(EDGREN_DATASETS)
    df["Is_real"] = df["DatasetType"].str.startswith("real")
    df["Is_simulated"] = df["DatasetType"].eq("simulated")
    for column in ["Recall_Sensitivity", "Precision", "F1"]:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    return df


def summarize_tools(df: pd.DataFrame, min_observations: int = 3) -> pd.DataFrame:
    paired = df.dropna(subset=["Recall_Sensitivity", "Precision"]).copy()
    summary = (
        paired.groupby("Tool", observed=False)
        .agg(
            n=("Tool", "size"),
            n_benchmarks=("Benchmark", "nunique"),
            mean_recall=("Recall_Sensitivity", "mean"),
            median_recall=("Recall_Sensitivity", "median"),
            mean_precision=("Precision", "mean"),
            median_precision=("Precision", "median"),
            mean_f1=("F1", "mean"),
            median_f1=("F1", "median"),
        )
        .reset_index()
    )
    summary = summary[summary["n"] >= min_observations].copy()
    summary["recall_minus_precision"] = summary["mean_recall"] - summary["mean_precision"]
    summary["recall_precision_ratio"] = summary["mean_recall"] / summary["mean_precision"].replace(0, np.nan)

    recall_q75 = summary["mean_recall"].quantile(0.75)
    precision_q25 = summary["mean_precision"].quantile(0.25)
    recall_q25 = summary["mean_recall"].quantile(0.25)
    precision_q75 = summary["mean_precision"].quantile(0.75)

    def classify(row: pd.Series) -> str:
        if row["mean_recall"] >= recall_q75 and row["mean_precision"] <= precision_q25:
            return "high recall / low precision"
        if row["mean_recall"] <= recall_q25 and row["mean_precision"] >= precision_q75:
            return "low recall / high precision"
        if row["mean_recall"] >= recall_q75 and row["mean_precision"] >= precision_q75:
            return "high recall / high precision"
        if row["mean_recall"] <= recall_q25 and row["mean_precision"] <= precision_q25:
            return "low recall / low precision"
        return "intermediate"

    summary["profile"] = summary.apply(classify, axis=1)
    summary = summary.sort_values(["recall_minus_precision", "mean_recall"], ascending=[False, False])
    return summary


def dataset_level_correlations(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    paired = df.dropna(subset=["Recall_Sensitivity", "Precision"]).copy()
    for label, subset in [
        ("all", paired),
        ("without_edgren", paired[~paired["Is_Edgren"]]),
        ("real", paired[paired["Is_real"]]),
        ("real_without_edgren", paired[paired["Is_real"] & ~paired["Is_Edgren"]]),
        ("simulated", paired[paired["Is_simulated"]]),
    ]:
        if subset.shape[0] < 3:
            continue
        pearson = stats.pearsonr(subset["Recall_Sensitivity"], subset["Precision"])
        spearman = stats.spearmanr(subset["Recall_Sensitivity"], subset["Precision"])
        rows.append(
            {
                "subset": label,
                "n": subset.shape[0],
                "pearson_r": pearson.statistic,
                "pearson_p": pearson.pvalue,
                "spearman_rho": spearman.statistic,
                "spearman_p": spearman.pvalue,
            }
        )
    return pd.DataFrame(rows)


def plot_tradeoff(summary: pd.DataFrame) -> None:
    colors = {
        "high recall / low precision": "#D95F02",
        "low recall / high precision": "#7570B3",
        "high recall / high precision": "#1B9E77",
        "low recall / low precision": "#666666",
        "intermediate": "#A6A6A6",
    }
    fig, ax = plt.subplots(figsize=(8.2, 6.1))
    for profile, subset in summary.groupby("profile", observed=False):
        ax.scatter(
            subset["mean_precision"],
            subset["mean_recall"],
            s=35 + subset["n"] * 7,
            alpha=0.78,
            color=colors[profile],
            edgecolor="white",
            linewidth=0.7,
            label=profile,
        )

    label_profiles = {"high recall / low precision", "low recall / high precision", "high recall / high precision"}
    for _, row in summary.iterrows():
        if row["profile"] in label_profiles or abs(row["recall_minus_precision"]) >= 0.25:
            ax.text(
                row["mean_precision"] + 0.008,
                row["mean_recall"] + 0.008,
                row["Tool"],
                fontsize=7.5,
            )

    ax.plot([0, 1], [0, 1], linestyle="--", color="#444444", linewidth=1.0, alpha=0.7)
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("Mean precision")
    ax.set_ylabel("Mean recall/sensitivity")
    ax.set_title("Tool-level recall--precision profiles across benchmark observations", fontweight="bold")
    ax.grid(color="#E6E6E6", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower left", fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "tool_precision_recall_tradeoff.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "tool_precision_recall_tradeoff.pdf", bbox_inches="tight")


def write_report(summary: pd.DataFrame, correlations: pd.DataFrame) -> None:
    high_recall_low_precision = summary[summary["profile"].eq("high recall / low precision")]
    low_recall_high_precision = summary[summary["profile"].eq("low recall / high precision")]
    top_gap = summary.sort_values("recall_minus_precision", ascending=False).head(10)

    lines = ["# Precision-recall trade-off analysis", ""]
    lines.append("Tool-level means were calculated from observations with both recall/sensitivity and precision available.")
    lines.append("")
    lines.append("## Correlations")
    for _, row in correlations.iterrows():
        lines.append(
            f"- {row['subset']}: n={int(row['n'])}, Pearson r={row['pearson_r']:.3f} (p={row['pearson_p']:.2e}), "
            f"Spearman rho={row['spearman_rho']:.3f} (p={row['spearman_p']:.2e})."
        )
    lines.append("")
    lines.append("## High recall / low precision profile")
    if high_recall_low_precision.empty:
        lines.append("- No tools met the quartile-based high recall / low precision definition.")
    else:
        for _, row in high_recall_low_precision.iterrows():
            lines.append(
                f"- {row['Tool']}: mean recall={row['mean_recall']:.3f}, mean precision={row['mean_precision']:.3f}, "
                f"mean F1={row['mean_f1']:.3f}, n={int(row['n'])}."
            )
    lines.append("")
    lines.append("## Low recall / high precision profile")
    if low_recall_high_precision.empty:
        lines.append("- No tools met the quartile-based low recall / high precision definition.")
    else:
        for _, row in low_recall_high_precision.iterrows():
            lines.append(
                f"- {row['Tool']}: mean recall={row['mean_recall']:.3f}, mean precision={row['mean_precision']:.3f}, "
                f"mean F1={row['mean_f1']:.3f}, n={int(row['n'])}."
            )
    lines.append("")
    lines.append("## Largest recall-minus-precision gaps")
    for _, row in top_gap.iterrows():
        lines.append(
            f"- {row['Tool']}: recall - precision = {row['recall_minus_precision']:.3f}; "
            f"mean recall={row['mean_recall']:.3f}; mean precision={row['mean_precision']:.3f}; n={int(row['n'])}."
        )
    (OUTPUT_DIR / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df = load_clean_data()
    df.to_csv(OUTPUT_DIR / "cleaned_all_data.csv", index=False)
    summary = summarize_tools(df)
    without_edgren_summary = summarize_tools(df[~df["Is_Edgren"]])
    real_summary = summarize_tools(df[df["Is_real"]])
    real_without_edgren_summary = summarize_tools(df[df["Is_real"] & ~df["Is_Edgren"]])
    simulated_summary = summarize_tools(df[df["Is_simulated"]])
    correlations = dataset_level_correlations(df)
    summary.to_csv(OUTPUT_DIR / "tool_precision_recall_summary.csv", index=False)
    without_edgren_summary.to_csv(OUTPUT_DIR / "tool_precision_recall_summary_without_edgren.csv", index=False)
    real_summary.to_csv(OUTPUT_DIR / "tool_precision_recall_summary_real.csv", index=False)
    real_without_edgren_summary.to_csv(OUTPUT_DIR / "tool_precision_recall_summary_real_without_edgren.csv", index=False)
    simulated_summary.to_csv(OUTPUT_DIR / "tool_precision_recall_summary_simulated.csv", index=False)
    correlations.to_csv(OUTPUT_DIR / "precision_recall_correlations.csv", index=False)
    plot_tradeoff(summary)
    write_report(summary, correlations)
    print((OUTPUT_DIR / "summary.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
