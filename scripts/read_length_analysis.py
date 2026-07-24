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
import statsmodels.formula.api as smf

from real_simulated_sensitivity import DEFAULT_INPUT, load_table, normalize_columns

OUTPUT_DIR = Path(__file__).resolve().parents[1] / "results" / "read_length"
COLOR = "#5A6FBB"


def canonical_read_length(value: float) -> float:
    if pd.isna(value):
        return np.nan
    value = float(value)
    if 49 <= value <= 51:
        return 50
    if 74 <= value <= 76:
        return 75
    if 99 <= value <= 101:
        return 100
    return np.nan


def fit_read_length_model(df: pd.DataFrame) -> pd.DataFrame:
    model_df = df.dropna(subset=["F1", "Read_length", "Benchmark", "Tool"]).copy()
    model_df["Read_length"] = pd.to_numeric(model_df["Read_length"], errors="coerce").map(canonical_read_length)
    model_df = model_df.dropna(subset=["Read_length"]).copy()
    model_df["Read_length"] = model_df["Read_length"].astype(int).astype(str)
    model_df["Read_length"] = pd.Categorical(model_df["Read_length"], categories=["100", "75", "50"])
    model_df["Benchmark"] = model_df["Benchmark"].astype("category")
    model_df["Tool"] = model_df["Tool"].astype("category")

    result = smf.mixedlm(
        "F1 ~ Read_length + C(Benchmark)",
        model_df,
        groups=model_df["Tool"],
    ).fit(reml=True)

    rows = []
    for term in ["Read_length[T.75]", "Read_length[T.50]"]:
        ci = result.conf_int().loc[term]
        rows.append(
            {
                "term": term,
                "beta": result.params[term],
                "ci_low": ci.iloc[0],
                "ci_high": ci.iloc[1],
                "p_value": result.pvalues[term],
                "n": int(model_df.shape[0]),
                "n_tools": int(model_df["Tool"].nunique()),
                "n_benchmarks": int(model_df["Benchmark"].nunique()),
                "converged": bool(result.converged),
            }
        )
    return pd.DataFrame(rows)


def format_p(value: float) -> str:
    if pd.isna(value):
        return "NA"
    exponent = math.floor(math.log10(abs(value))) if value else 0
    return f"{value:.2e}" if exponent <= -4 else f"{value:.4f}"


def plot_read_length_f1(analysis: pd.DataFrame, summary: pd.DataFrame) -> None:
    lengths = [50, 75, 100]
    data = [analysis.loc[analysis["Read_length"].eq(length), "F1"].dropna().to_numpy() for length in lengths]

    fig, ax = plt.subplots(figsize=(5.8, 4.4))
    box = ax.boxplot(
        data,
        tick_labels=[f"{length} bp" for length in lengths],
        patch_artist=True,
        widths=0.55,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.3},
        boxprops={"color": "#444444", "linewidth": 0.9},
        whiskerprops={"color": "#444444", "linewidth": 0.9},
        capprops={"color": "#444444", "linewidth": 0.9},
    )
    for patch in box["boxes"]:
        patch.set_facecolor(COLOR)
        patch.set_alpha(0.33)

    rng = np.random.default_rng(123)
    for xpos, values in enumerate(data, start=1):
        jitter = rng.uniform(-0.08, 0.08, len(values))
        ax.scatter(
            xpos + jitter,
            values,
            s=18,
            color=COLOR,
            alpha=0.55,
            edgecolor="white",
            linewidth=0.25,
            zorder=3,
        )
        ax.text(xpos, 1.03, f"n={len(values)}", ha="center", va="bottom", fontsize=9)

    ax.plot(
        range(1, len(lengths) + 1),
        [summary.loc[summary["Read_length"].eq(length), "mean"].iloc[0] for length in lengths],
        color="#D17C3F",
        marker="D",
        linewidth=1.4,
        markersize=5,
        label="Mean",
    )
    ax.set_ylim(-0.03, 1.10)
    ax.set_ylabel("F1-score")
    ax.set_xlabel("Read length")
    ax.set_title("F1-score distribution by read length")
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()

    png_path = OUTPUT_DIR / "read_length_f1.png"
    pdf_path = OUTPUT_DIR / "read_length_f1.pdf"
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df = normalize_columns(load_table(DEFAULT_INPUT))
    df["Read_length"] = pd.to_numeric(df.get("Read_length"), errors="coerce").map(canonical_read_length)
    analysis = df.dropna(subset=["F1", "Read_length"]).copy()
    analysis["Read_length"] = analysis["Read_length"].astype(int)
    analysis.to_csv(OUTPUT_DIR / "cleaned_read_length_data.csv", index=False)

    summary = (
        analysis.groupby("Read_length", observed=False)["F1"]
        .agg(["count", "mean", "median", "std"])
        .reset_index()
        .sort_values("Read_length")
    )
    summary.to_csv(OUTPUT_DIR / "read_length_f1_summary.csv", index=False)
    plot_read_length_f1(analysis, summary)

    groups = [analysis.loc[analysis["Read_length"].eq(length), "F1"].dropna() for length in [50, 75, 100]]
    h_stat, h_p = stats.kruskal(*groups)
    pd.DataFrame([{"H": h_stat, "p_value": h_p}]).to_csv(OUTPUT_DIR / "read_length_kruskal.csv", index=False)

    model = fit_read_length_model(analysis)
    model.to_csv(OUTPUT_DIR / "read_length_mixed_model.csv", index=False)

    lines = ["# Read-length analysis", ""]
    lines.append(f"Input observations with 50/75/100 bp and F1: {len(analysis)}")
    lines.append("")
    lines.append("## Mean F1 by read length")
    for _, row in summary.iterrows():
        lines.append(
            f"- {int(row['Read_length'])} bp: n={int(row['count'])}, "
            f"mean={row['mean']:.3f}, median={row['median']:.3f}."
        )
    lines.append("")
    lines.append(f"Kruskal-Wallis H={h_stat:.2f}, p={format_p(h_p)}.")
    lines.append("")
    lines.append("## Mixed-effects model, reference 100 bp")
    for _, row in model.iterrows():
        label = "75 bp" if "75" in row["term"] else "50 bp"
        lines.append(
            f"- {label}: beta={row['beta']:.3f}, 95% CI {row['ci_low']:.3f} to {row['ci_high']:.3f}, "
            f"p={format_p(row['p_value'])}."
        )
    (OUTPUT_DIR / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print((OUTPUT_DIR / "summary.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
