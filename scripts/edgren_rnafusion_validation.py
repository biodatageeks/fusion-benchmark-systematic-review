from __future__ import annotations

import os
import re
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(__file__).resolve().parents[1] / ".cache"))

import matplotlib.pyplot as plt
import pandas as pd

from real_simulated_sensitivity import DEFAULT_INPUT, load_table, normalize_columns


INPUT = Path(__file__).resolve().parents[1] / "data" / "raw" / "edgren_rnafusion.xlsx"
OUTPUT_DIR = Path(__file__).resolve().parents[1] / "results" / "edgren_rnafusion_validation"
TOOL_COLUMNS = {
    "Arriba": "arriba",
    "FusionCatcher": "fusioncatcher",
    "STAR-Fusion": "starfusion",
}

ALIASES = {
    "SEPT10": {"SEPT10", "SEPTIN10"},
    "SEPTIN10": {"SEPT10", "SEPTIN10"},
    "GCN1L1": {"GCN1L1", "GCN1"},
    "GCN1": {"GCN1L1", "GCN1"},
    "RP4-791K14.2": {"RP4-791K14.2", "AL035685.1"},
    "AL035685.1": {"RP4-791K14.2", "AL035685.1"},
    "UHRF1BP1": {"UHRF1BP1", "BLTP3A"},
    "BLTP3A": {"UHRF1BP1", "BLTP3A"},
}


def normalize_gene(value: object) -> str:
    text = str(value).strip()
    text = re.sub(r"\s+", "", text)
    return text.upper()


def gene_tokens(value: object) -> frozenset[str]:
    if pd.isna(value):
        return frozenset()
    tokens: set[str] = set()
    for token in str(value).split("|"):
        gene = normalize_gene(token)
        if not gene or gene == "NAN":
            continue
        tokens.add(gene)
        tokens.update(ALIASES.get(gene, set()))
    return frozenset(tokens)


def fusion_to_pair(value: object) -> tuple[frozenset[str], frozenset[str]] | None:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if "--" not in text:
        return None
    left, right = text.split("--", 1)
    left_tokens = gene_tokens(left)
    right_tokens = gene_tokens(right)
    if not left_tokens or not right_tokens:
        return None
    return left_tokens, right_tokens


def pair_key(pair: tuple[frozenset[str], frozenset[str]]) -> tuple[str, str]:
    left, right = pair
    return "|".join(sorted(left)), "|".join(sorted(right))


def undirected_pair_key(pair: tuple[frozenset[str], frozenset[str]]) -> tuple[str, str]:
    left, right = pair_key(pair)
    return tuple(sorted([left, right]))


def pairs_match(
    predicted: tuple[frozenset[str], frozenset[str]],
    truth: tuple[frozenset[str], frozenset[str]],
    directional: bool,
) -> bool:
    pred_left, pred_right = predicted
    truth_left, truth_right = truth
    direct = bool(pred_left & truth_left) and bool(pred_right & truth_right)
    if directional:
        return direct
    reverse = bool(pred_left & truth_right) and bool(pred_right & truth_left)
    return direct or reverse


def best_truth_match(
    predicted: tuple[frozenset[str], frozenset[str]],
    truth_pairs: list[tuple[frozenset[str], frozenset[str]]],
    directional: bool,
) -> int | None:
    for index, truth in enumerate(truth_pairs):
        if pairs_match(predicted, truth, directional=directional):
            return index
    return None


def extract_truth(df: pd.DataFrame) -> pd.DataFrame:
    truth = df[["Gene1", "Gene2"]].dropna().copy()
    truth["truth_pair"] = list(zip(truth["Gene1"].map(gene_tokens), truth["Gene2"].map(gene_tokens)))
    truth["truth_key"] = truth["truth_pair"].map(pair_key)
    truth = truth.drop_duplicates("truth_key").reset_index(drop=True)
    truth["truth_id"] = truth.index + 1
    truth["truth_fusion"] = truth["Gene1"].astype(str) + "--" + truth["Gene2"].astype(str)
    return truth


def extract_predictions(df: pd.DataFrame, tool: str, column: str) -> pd.DataFrame:
    predictions = df.loc[df[column].notna(), ["Fusion", column]].copy()
    predictions["raw_call_id"] = predictions.index + 1
    predictions["pred_pair"] = predictions["Fusion"].map(fusion_to_pair)
    predictions = predictions.dropna(subset=["pred_pair"]).copy()
    predictions["pred_key"] = predictions["pred_pair"].map(pair_key)
    predictions = predictions.drop_duplicates("pred_key").reset_index(drop=True)
    predictions["tool"] = tool
    return predictions[["tool", "Fusion", "pred_pair", "pred_key"]]


def evaluate_tool(
    predictions: pd.DataFrame,
    truth: pd.DataFrame,
    directional: bool,
) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    truth_pairs = truth["truth_pair"].tolist()
    pred = predictions.copy()
    pred["truth_index"] = pred["pred_pair"].map(
        lambda pair: best_truth_match(pair, truth_pairs, directional=directional)
    )
    pred["is_tp"] = pred["truth_index"].notna()
    tp_truth_indices = set(pred.loc[pred["is_tp"], "truth_index"].astype(int))
    fn = truth.loc[~truth.index.isin(tp_truth_indices)].copy()
    tp = pred[pred["is_tp"]].copy()
    fp = pred[~pred["is_tp"]].copy()

    tp_count = len(tp_truth_indices)
    fp_count = len(fp)
    fn_count = len(fn)
    precision = tp_count / (tp_count + fp_count) if tp_count + fp_count else float("nan")
    recall = tp_count / (tp_count + fn_count) if tp_count + fn_count else float("nan")
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else float("nan")

    summary = {
        "tool": predictions["tool"].iloc[0],
        "matching": "directional" if directional else "undirected",
        "unique_predictions": len(pred),
        "TP": tp_count,
        "FP": fp_count,
        "FN": fn_count,
        "precision": precision,
        "recall": recall,
        "F1": f1,
    }
    return summary, tp, fp, fn


def load_paper_reference() -> pd.DataFrame:
    master = normalize_columns(load_table(DEFAULT_INPUT))
    truth_total = pd.to_numeric(master["TruthTotal"], errors="coerce")
    if "Is_Edgren_Binary" in master.columns:
        edgren_mask = master["Is_Edgren_Binary"].astype(str).str.lower().str.strip().eq("yes")
    else:
        edgren_mask = master["Is_Edgren"].astype(str).str.lower().str.strip().isin(["true", "yes", "1"])
    reference = master[
        edgren_mask
        & truth_total.eq(99)
        & master["Tool"].isin(TOOL_COLUMNS)
    ][["Benchmark", "Dataset", "Tool", "Precision", "Recall", "F1"]].copy()
    reference = reference.rename(
        columns={
            "Tool": "tool",
            "Precision": "precision",
            "Recall": "recall",
        }
    )
    reference["source"] = "B" + reference["Benchmark"].astype(str) + ", " + reference["Dataset"].astype(str)
    reference.to_csv(OUTPUT_DIR / "published_edgren_tp99_reference_metrics.csv", index=False)
    return reference


def plot_metrics(summary: pd.DataFrame, paper_reference: pd.DataFrame | None = None) -> None:
    directional = summary[summary["matching"].eq("undirected")].copy()
    tools = directional["tool"].tolist()
    metrics = ["precision", "recall", "F1"]
    colors = {"precision": "#4C78A8", "recall": "#F58518", "F1": "#54A24B"}
    reference_lookup: dict[tuple[str, str], list[float]] = {}
    if paper_reference is not None and not paper_reference.empty:
        for _, row in paper_reference.iterrows():
            for metric in metrics:
                if pd.isna(row[metric]):
                    continue
                reference_lookup.setdefault((row["tool"], metric), []).append(row[metric])

    fig, ax = plt.subplots(figsize=(7.4, 4.9))
    width = 0.23
    x_positions = range(len(tools))
    reference_label_added = False
    for offset, metric in enumerate(metrics):
        values = directional[metric].tolist()
        positions = [x + (offset - 1) * width for x in x_positions]
        ax.bar(positions, values, width=width, label=metric.capitalize(), color=colors[metric])
        for x, value in zip(positions, values):
            ax.text(x, value + 0.015, f"{value:.2f}", ha="center", va="bottom", fontsize=8)
        for tool, x in zip(tools, positions):
            reference_values = reference_lookup.get((tool, metric), [])
            if not reference_values:
                continue
            for reference_value in reference_values:
                ax.hlines(
                    reference_value,
                    x - width * 0.46,
                    x + width * 0.46,
                    colors="black",
                    linewidth=2.2,
                    zorder=5,
                    label="Published Edgren TP99 value(s)" if not reference_label_added else None,
                )
                reference_label_added = True

    ax.set_xticks(list(x_positions))
    ax.set_xticklabels(tools)
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("Metric value")
    ax.legend(frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, 1.16))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "edgren_rnafusion_metrics.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "edgren_rnafusion_metrics.pdf", bbox_inches="tight")
    plt.close(fig)


def write_overlap_summary() -> None:
    all_tp = pd.read_csv(OUTPUT_DIR / "all_true_positives.csv")
    undirected = all_tp[all_tp["matching"].eq("undirected")].copy()
    truth_sets = {
        tool: set(group["truth_index"].astype(int))
        for tool, group in undirected.groupby("tool")
    }
    rows = []
    tools = sorted(truth_sets)
    for tool in tools:
        rows.append({"comparison": tool, "truth_fusions": len(truth_sets[tool])})
    if len(tools) == 3:
        first, second, third = tools
        rows.extend(
            [
                {
                    "comparison": f"{first} ∩ {second}",
                    "truth_fusions": len(truth_sets[first] & truth_sets[second]),
                },
                {
                    "comparison": f"{first} ∩ {third}",
                    "truth_fusions": len(truth_sets[first] & truth_sets[third]),
                },
                {
                    "comparison": f"{second} ∩ {third}",
                    "truth_fusions": len(truth_sets[second] & truth_sets[third]),
                },
                {
                    "comparison": "all three tools",
                    "truth_fusions": len(set.intersection(*truth_sets.values())),
                },
                {
                    "comparison": "any tool",
                    "truth_fusions": len(set.union(*truth_sets.values())),
                },
                {
                    "comparison": "no tool",
                    "truth_fusions": 99 - len(set.union(*truth_sets.values())),
                },
            ]
        )
    pd.DataFrame(rows).to_csv(OUTPUT_DIR / "edgren_rnafusion_overlap_summary.csv", index=False)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df = load_table(INPUT)
    truth = extract_truth(df)
    truth.to_csv(OUTPUT_DIR / "truthset_99_gene_pairs.csv", index=False)

    all_summaries = []
    all_tp = []
    all_fp = []
    all_fn = []
    raw_counts = []

    for tool, column in TOOL_COLUMNS.items():
        raw_counts.append({"tool": tool, "raw_calls": int(df[column].notna().sum())})
        predictions = extract_predictions(df, tool, column)
        predictions.to_csv(OUTPUT_DIR / f"{tool.lower().replace('-', '_')}_unique_predictions.csv", index=False)
        for directional in [True, False]:
            summary, tp, fp, fn = evaluate_tool(predictions, truth, directional=directional)
            all_summaries.append(summary)
            label = summary["matching"]
            tp.assign(tool=tool, matching=label).to_csv(OUTPUT_DIR / f"{tool.lower().replace('-', '_')}_{label}_tp.csv", index=False)
            fp.assign(tool=tool, matching=label).to_csv(OUTPUT_DIR / f"{tool.lower().replace('-', '_')}_{label}_fp.csv", index=False)
            fn.assign(tool=tool, matching=label).to_csv(OUTPUT_DIR / f"{tool.lower().replace('-', '_')}_{label}_fn.csv", index=False)
            all_tp.append(tp.assign(tool=tool, matching=label))
            all_fp.append(fp.assign(tool=tool, matching=label))
            all_fn.append(fn.assign(tool=tool, matching=label))

    summary = pd.DataFrame(all_summaries)
    raw_counts = pd.DataFrame(raw_counts)
    summary = summary.merge(raw_counts, on="tool", how="left")
    summary = summary[
        ["tool", "matching", "raw_calls", "unique_predictions", "TP", "FP", "FN", "precision", "recall", "F1"]
    ]
    summary.to_csv(OUTPUT_DIR / "edgren_rnafusion_metrics_summary.csv", index=False)
    pd.concat(all_tp, ignore_index=True).to_csv(OUTPUT_DIR / "all_true_positives.csv", index=False)
    pd.concat(all_fp, ignore_index=True).to_csv(OUTPUT_DIR / "all_false_positives.csv", index=False)
    pd.concat(all_fn, ignore_index=True).to_csv(OUTPUT_DIR / "all_false_negatives.csv", index=False)
    paper_reference = load_paper_reference()
    plot_metrics(summary, paper_reference)
    write_overlap_summary()

    print(summary.to_string(index=False, float_format=lambda value: f"{value:.3f}"))
    print("\nOverlap summary:")
    print(pd.read_csv(OUTPUT_DIR / "edgren_rnafusion_overlap_summary.csv").to_string(index=False))
    print(f"\nSaved results to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
