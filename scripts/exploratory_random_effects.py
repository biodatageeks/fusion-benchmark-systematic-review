from __future__ import annotations

import math
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(__file__).resolve().parents[1] / ".cache"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from real_simulated_sensitivity import DEFAULT_INPUT, load_table, normalize_columns


REPO_ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = REPO_ROOT / "results" / "tables"
FIG_DIR = REPO_ROOT / "results" / "figures"
MIN_BENCHMARKS_PER_TOOL = 3
BOOT_ITER = 2000


def f1_score(precision: pd.Series, recall: pd.Series) -> pd.Series:
    denominator = precision + recall
    return np.where(denominator.gt(0), 2 * precision * recall / denominator, np.nan)


def bootstrap_f1_variance(tp: float, fp: float, fn: float, seed: int) -> float:
    if pd.isna(tp) or pd.isna(fp) or pd.isna(fn):
        return np.nan
    tp, fp, fn = int(round(tp)), int(round(fp)), int(round(fn))
    n_pred = tp + fp
    n_true = tp + fn
    if n_pred <= 0 or n_true <= 0:
        return np.nan

    precision = tp / n_pred
    recall = tp / n_true
    rng = np.random.default_rng(seed)
    tp_from_predictions = rng.binomial(n_pred, precision, BOOT_ITER)
    tp_from_truth = rng.binomial(n_true, recall, BOOT_ITER)
    tp_boot = np.minimum(tp_from_predictions, tp_from_truth)
    precision_boot = tp_boot / n_pred
    recall_boot = tp_boot / n_true
    denominator = precision_boot + recall_boot
    f1_boot = np.zeros_like(denominator, dtype=float)
    nonzero = denominator > 0
    f1_boot[nonzero] = 2 * precision_boot[nonzero] * recall_boot[nonzero] / denominator[nonzero]
    return float(np.var(f1_boot, ddof=1))


def der_simonian_laird(values: np.ndarray, variances: np.ndarray) -> dict[str, float]:
    mask = np.isfinite(values) & np.isfinite(variances) & (variances > 0)
    values = values[mask]
    variances = variances[mask]
    k = len(values)
    if k < 2:
        return {}

    fixed_weights = 1 / variances
    fixed_mean = np.sum(fixed_weights * values) / np.sum(fixed_weights)
    q_stat = np.sum(fixed_weights * (values - fixed_mean) ** 2)
    df = k - 1
    c_value = np.sum(fixed_weights) - (np.sum(fixed_weights**2) / np.sum(fixed_weights))
    tau2 = max(0.0, (q_stat - df) / c_value) if c_value > 0 else 0.0

    random_weights = 1 / (variances + tau2)
    pooled = np.sum(random_weights * values) / np.sum(random_weights)
    se = math.sqrt(1 / np.sum(random_weights))
    ci_low = pooled - 1.96 * se
    ci_high = pooled + 1.96 * se
    i2 = max(0.0, (q_stat - df) / q_stat) * 100 if q_stat > 0 else 0.0
    q_p = float(stats.chi2.sf(q_stat, df))
    return {
        "k": k,
        "pooled": pooled,
        "ci_lb": ci_low,
        "ci_ub": ci_high,
        "I2": i2,
        "Q": q_stat,
        "Q_p": q_p,
        "tau2": tau2,
    }


def prepare_metrics() -> pd.DataFrame:
    df = normalize_columns(load_table(DEFAULT_INPUT)).copy()
    for column in ["TP", "FP", "FN", "Precision", "Recall", "F1"]:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    df["n_pred"] = df["TP"] + df["FP"]
    df["n_true"] = df["TP"] + df["FN"]
    df["prec_yi"] = np.where(df["n_pred"].gt(0), df["TP"] / df["n_pred"], np.nan)
    df["prec_vi"] = np.where(
        df["n_pred"].gt(0),
        df["prec_yi"] * (1 - df["prec_yi"]) / df["n_pred"],
        np.nan,
    )
    df["sens_yi"] = np.where(df["n_true"].gt(0), df["TP"] / df["n_true"], np.nan)
    df["sens_vi"] = np.where(
        df["n_true"].gt(0),
        df["sens_yi"] * (1 - df["sens_yi"]) / df["n_true"],
        np.nan,
    )
    df["f1_yi"] = f1_score(df["prec_yi"], df["sens_yi"])
    df["f1_vi"] = np.nan
    valid_f1 = df["TP"].notna() & df["FP"].notna() & df["FN"].notna()
    for seed_offset, row_index in enumerate(df.index[valid_f1]):
        row = df.loc[row_index]
        df.loc[row_index, "f1_vi"] = bootstrap_f1_variance(
            row["TP"],
            row["FP"],
            row["FN"],
            seed=42 + seed_offset,
        )
    return df


def select_tools(df: pd.DataFrame) -> list[str]:
    benchmark_counts = df[["Tool", "Benchmark"]].dropna().drop_duplicates().groupby("Tool").size()
    return sorted(benchmark_counts[benchmark_counts >= MIN_BENCHMARKS_PER_TOOL].index.tolist())


def pooled_estimates(df: pd.DataFrame, tools: list[str]) -> pd.DataFrame:
    rows = []
    metric_columns = {
        "precision": ("prec_yi", "prec_vi"),
        "sensitivity": ("sens_yi", "sens_vi"),
        "f1": ("f1_yi", "f1_vi"),
    }
    for tool in tools:
        subset = df[df["Tool"].eq(tool)]
        for metric, (value_col, variance_col) in metric_columns.items():
            result = der_simonian_laird(
                subset[value_col].to_numpy(dtype=float),
                subset[variance_col].to_numpy(dtype=float),
            )
            if result:
                rows.append({"tool": tool, "metric": metric, **result})
    return pd.DataFrame(rows)


def subgroup_real_vs_sim(df: pd.DataFrame, tools: list[str]) -> pd.DataFrame:
    rows = []
    for tool in tools:
        subset = df[
            df["Tool"].eq(tool)
            & df["f1_yi"].notna()
            & df["f1_vi"].notna()
            & df["f1_vi"].gt(0)
            & df["DatasetType"].isin(["real", "simulated"])
        ].copy()
        if len(subset) < 3 or subset["DatasetType"].nunique() < 2:
            continue

        groups = []
        for dataset_type, group in subset.groupby("DatasetType"):
            result = der_simonian_laird(
                group["f1_yi"].to_numpy(dtype=float),
                group["f1_vi"].to_numpy(dtype=float),
            )
            if result:
                groups.append((dataset_type, result))
        if len(groups) < 2:
            continue

        weighted_means = {name: result["pooled"] for name, result in groups}
        delta = weighted_means.get("real", np.nan) - weighted_means.get("simulated", np.nan)
        rows.append(
            {
                "tool": tool,
                "k": len(subset),
                "pooled_f1_real": weighted_means.get("real", np.nan),
                "pooled_f1_simulated": weighted_means.get("simulated", np.nan),
                "delta_real_minus_simulated": delta,
                "direction_consistent": bool(delta < 0) if np.isfinite(delta) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def edgren_sensitivity(df: pd.DataFrame, tools: list[str]) -> pd.DataFrame:
    rows = []
    if "Is_Edgren_Binary" in df.columns:
        is_edgren = df["Is_Edgren_Binary"].astype(str).str.lower().str.strip().eq("yes")
    else:
        is_edgren = df["Is_Edgren"].astype(str).str.lower().str.strip().isin(["true", "yes", "1"])
    for tool in tools:
        full = df[df["Tool"].eq(tool) & df["f1_yi"].notna() & df["f1_vi"].notna() & df["f1_vi"].gt(0)]
        no_edgren = full[~is_edgren.loc[full.index]]
        full_result = der_simonian_laird(full["f1_yi"].to_numpy(float), full["f1_vi"].to_numpy(float))
        no_edgren_result = der_simonian_laird(no_edgren["f1_yi"].to_numpy(float), no_edgren["f1_vi"].to_numpy(float))
        if full_result and no_edgren_result:
            rows.append(
                {
                    "tool": tool,
                    "edgren_rows": int(len(full) - len(no_edgren)),
                    "f1_with": full_result["pooled"],
                    "I2_with": full_result["I2"],
                    "k_with": full_result["k"],
                    "f1_without": no_edgren_result["pooled"],
                    "I2_without": no_edgren_result["I2"],
                    "k_without": no_edgren_result["k"],
                    "delta_f1": no_edgren_result["pooled"] - full_result["pooled"],
                    "delta_I2": no_edgren_result["I2"] - full_result["I2"],
                }
            )
    return pd.DataFrame(rows)


def edgren_truth_set_drift(df: pd.DataFrame) -> pd.DataFrame:
    columns = ["Tool", "TP", "FP", "FN", "prec_yi", "sens_yi", "f1_yi"]
    truth_27 = df[df["Dataset"].eq("Dataset7.6")][columns].copy()
    truth_99 = df[df["Dataset"].eq("Dataset7.7")][columns].copy()
    truth_27.columns = ["Tool", "TP_27", "FP_27", "FN_27", "precision_27", "recall_27", "f1_27"]
    truth_99.columns = ["Tool", "TP_99", "FP_99", "FN_99", "precision_99", "recall_99", "f1_99"]
    drift = truth_27.merge(truth_99, on="Tool", how="inner")
    drift["delta_f1"] = drift["f1_99"] - drift["f1_27"]
    drift["delta_precision"] = drift["precision_99"] - drift["precision_27"]
    drift["delta_recall"] = drift["recall_99"] - drift["recall_27"]
    return drift.sort_values("delta_f1").reset_index(drop=True)


def plot_edgren_drift_scatter(drift: pd.DataFrame) -> None:
    data = drift.dropna(subset=["f1_27", "f1_99"]).copy()
    fig, ax = plt.subplots(figsize=(6.6, 6.4))
    limit = max(0.7, data[["f1_27", "f1_99"]].max().max() + 0.06)
    colors = np.where(data["delta_f1"].lt(0), "#D17C3F", "#5A6FBB")
    ax.plot([0, limit], [0, limit], color="#555555", linestyle="--", linewidth=1, alpha=0.65)
    ax.scatter(data["f1_27"], data["f1_99"], s=58, color=colors, edgecolor="white", linewidth=0.6, zorder=3)
    for _, row in data.iterrows():
        ax.annotate(row["Tool"], (row["f1_27"], row["f1_99"]), xytext=(4, 4), textcoords="offset points", fontsize=7)
    ax.set_xlim(0, limit)
    ax.set_ylim(0, limit)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("F1-score with 27-fusion truth set")
    ax.set_ylabel("F1-score with 99-fusion truth set")
    ax.set_title("Truth-set-induced F1-score changes on identical Edgren data")
    ax.grid(color="#E3E3E3")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "edgren_drift_scatter.png", dpi=600, bbox_inches="tight")
    fig.savefig(FIG_DIR / "edgren_drift_scatter.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_i2(pooled: pd.DataFrame) -> None:
    metric_order = ["precision", "sensitivity", "f1"]
    f1_order = (
        pooled[pooled["metric"].eq("f1")]
        .sort_values("I2", ascending=False)["tool"]
        .tolist()
    )
    pivot = pooled.pivot(index="tool", columns="metric", values="I2").reindex(f1_order)
    fig, ax = plt.subplots(figsize=(8.2, max(4.5, 0.38 * len(pivot) + 1.2)))
    y_pos = np.arange(len(pivot))
    height = 0.24
    colors = {"precision": "#D17C3F", "sensitivity": "#4C9A6A", "f1": "#5A6FBB"}
    for offset, metric in zip([-height, 0, height], metric_order):
        ax.barh(y_pos + offset, pivot[metric], height=height, label=metric, color=colors[metric], alpha=0.86)
    ax.axvline(75, color="#9C2F2F", linestyle="--", linewidth=1, label="I² = 75%")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(pivot.index)
    ax.set_xlim(0, 105)
    ax.set_xlabel("I² (%)")
    ax.set_title("Between-benchmark heterogeneity in exploratory random-effects summaries")
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(axis="x", color="#E3E3E3")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "heterogeneity_I2.png", dpi=600, bbox_inches="tight")
    fig.savefig(FIG_DIR / "heterogeneity_I2.pdf", bbox_inches="tight")
    plt.close(fig)


def write_summary(pooled: pd.DataFrame, subgroup: pd.DataFrame) -> None:
    f1 = pooled[pooled["metric"].eq("f1")]
    i2_over_90 = int((f1["I2"] > 90).sum())
    total = int(len(f1))
    subgroup_consistent = int(subgroup["direction_consistent"].sum()) if not subgroup.empty else 0
    lines = [
        "# Exploratory random-effects summaries",
        "",
        f"Tools with at least {MIN_BENCHMARKS_PER_TOOL} benchmarks: {total}.",
        f"F1 I² exceeded 90% for {i2_over_90} of {total} recurrently evaluated tools.",
        f"Real-versus-simulated subgroup direction was lower for real datasets in {subgroup_consistent} of {len(subgroup)} tools with sufficient data.",
    ]
    (TABLE_DIR / "random_effects_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    df = prepare_metrics()
    df.to_csv(TABLE_DIR / "random_effects_input_effect_sizes.csv", index=False)
    tools = select_tools(df)
    pooled = pooled_estimates(df, tools)
    pooled.to_csv(TABLE_DIR / "pooled_estimates.csv", index=False)
    subgroup = subgroup_real_vs_sim(df, tools)
    subgroup.to_csv(TABLE_DIR / "subgroup_real_vs_sim.csv", index=False)
    edgren = edgren_sensitivity(df, tools)
    edgren.to_csv(TABLE_DIR / "sensitivity_edgren.csv", index=False)
    drift = edgren_truth_set_drift(df)
    drift.to_csv(TABLE_DIR / "edgren_truth_set_drift.csv", index=False)
    plot_edgren_drift_scatter(drift)
    plot_i2(pooled)
    write_summary(pooled, subgroup)
    print((TABLE_DIR / "random_effects_summary.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
