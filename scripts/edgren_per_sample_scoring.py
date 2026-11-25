"""Edgren RNAfusion predictions scored per sequencing run against the reference lists 27, 38, 53 and 99.

Each run is scored only against the events of its own cell line. Matching rules are those of
edgren_rnafusion_validation.py (imported, not copied): gene-symbol aliases, "|" alternatives,
directional and undirected matching, one-to-one assignment of predictions to reference events.

Inputs
  data/raw/edgren_predictions/<RNAfusion version>/<SRA run>/<caller>.csv   raw caller outputs (CSV)
  data/raw/edgren_predictions/ena_runs.tsv                                 ENA run metadata
  data/raw/edgren_truthsets/truthset_<N>/<SRA run>_truthset_<N>.csv        reference events per run
Outputs: results/edgren_per_sample/
"""
from __future__ import annotations

import os
import re
from itertools import product
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".mplconfig"))

import matplotlib.pyplot as plt
import pandas as pd

from edgren_rnafusion_validation import evaluate_tool, extract_predictions, extract_truth

ROOT = Path(__file__).resolve().parents[1]
PRED_DIR = ROOT / "data" / "raw" / "edgren_predictions"
TRUTH_DIR = ROOT / "data" / "raw" / "edgren_truthsets"
POOLED_REFERENCE = ROOT / "results" / "edgren_rnafusion_validation" / "edgren_rnafusion_metrics_summary.csv"
OUTPUT_DIR = ROOT / "results" / "edgren_per_sample"

VERSIONS = ["3.0.1", "4.1.3"]
TOOLS = ["Arriba", "FusionCatcher", "STAR-Fusion"]
LISTS = [27, 38, 53, 99]
MATCHING = {"directional": True, "undirected": False}


def arriba_genes(label: object) -> list[str]:
    """Arriba lists the two flanking genes of an intergenic breakpoint as 'A(dist),B(dist)'.
    Each flanking gene becomes a separate candidate, as in the fusion-report table used for
    the pooled analysis (edgren_rnafusion_validation.py) (e.g. BCAS1--B4GALNT2 and BCAS1--B4GALNT2P1)."""
    return [re.sub(r"\(.*\)$", "", part).strip() for part in str(label).split(",") if part.strip()]


def read_calls(version: str, run: str, tool: str) -> pd.DataFrame:
    raw = pd.read_csv(PRED_DIR / version / run / f"{tool}.csv", dtype=str)
    if tool == "Arriba":
        fusions = [f"{a}--{b}" for g1, g2 in zip(raw["#gene1"], raw["gene2"])
                   for a, b in product(arriba_genes(g1), arriba_genes(g2))]
    elif tool == "FusionCatcher":
        fusions = (raw["Gene_1_symbol(5end_fusion_partner)"] + "--" + raw["Gene_2_symbol(3end_fusion_partner)"]).tolist()
    else:
        fusions = raw["#FusionName"].tolist()
    return pd.DataFrame({"Fusion": fusions, "call": True})


def read_truth(size: int, runs: list[str]) -> pd.DataFrame:
    frames = [pd.read_csv(TRUTH_DIR / f"truthset_{size}" / f"{run}_truthset_{size}.csv", dtype=str) for run in runs]
    df = pd.concat(frames, ignore_index=True).rename(columns={"gene1_raw": "Gene1", "gene2_raw": "Gene2"})
    return extract_truth(df)


def score(version: str, tool: str, runs: list[str], truth_runs: list[str], size: int) -> list[dict]:
    calls = pd.concat([read_calls(version, run, tool) for run in runs], ignore_index=True)
    predictions = extract_predictions(calls, tool, "call")
    truth = read_truth(size, truth_runs)
    rows = []
    for label, directional in MATCHING.items():
        summary, *_ = evaluate_tool(predictions, truth, directional=directional)
        summary.update({"rnafusion": version, "truthset": size, "n_reference": len(truth)})
        rows.append(summary)
    return rows


def micro(df: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    out = df.groupby(keys, as_index=False)[["n_reference", "unique_predictions", "TP", "FP", "FN"]].sum()
    out["precision"] = out.TP / (out.TP + out.FP)
    out["recall"] = out.TP / (out.TP + out.FN)
    out["F1"] = 2 * out.precision * out.recall / (out.precision + out.recall)
    return out


def plot(per_run: pd.DataFrame) -> None:
    data = per_run[per_run.matching.eq("undirected")]
    fig, axes = plt.subplots(1, len(LISTS), figsize=(13, 3.8), sharey=True)
    colors = {"3.0.1": "#4C78A8", "4.1.3": "#F58518"}
    for ax, size in zip(axes, LISTS):
        sub = data[data.truthset.eq(size)]
        for i, tool in enumerate(TOOLS):
            for j, version in enumerate(VERSIONS):
                values = sub[sub.tool.eq(tool) & sub.rnafusion.eq(version)].F1
                x = [i + (j - 0.5) * 0.25] * len(values)
                ax.scatter(x, values, s=22, color=colors[version], alpha=0.8,
                           label=f"RNAfusion {version}" if (size == LISTS[0] and i == 0) else None)
        ax.set_xticks(range(len(TOOLS)))
        ax.set_xticklabels(TOOLS, rotation=30, ha="right")
        ax.set_title(f"Reference list {size}")
        ax.set_ylim(0, 1)
        ax.grid(axis="y", color="#E6E6E6")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("F1 per sequencing run (undirected)")
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "edgren_per_run_f1.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "edgren_per_run_f1.pdf", bbox_inches="tight")
    plt.close(fig)


PUBLISHED_TP99 = ROOT / "results" / "edgren_rnafusion_validation" / "published_edgren_tp99_reference_metrics.csv"
LIST_COLORS = {27: "#C6DBEF", 38: "#6BAED6", 53: "#2171B5", 99: "#08306B"}


def plot_by_list(pooled: pd.DataFrame, version: str = "3.0.1") -> None:
    """Main-text figure: precision, recall and F1 per caller for each reference list (each run scored
    against its own cell line's events, summed over the four cell lines), RNAfusion `version`, undirected
    matching. Black markers: published Edgren TP99 values for the same callers."""
    d = pooled[pooled.rnafusion.eq(version) & pooled.matching.eq("undirected") & pooled.pooling.str.startswith("sum over")]
    published = pd.read_csv(PUBLISHED_TP99) if PUBLISHED_TP99.exists() else pd.DataFrame()
    metrics = [("precision", "Precision"), ("recall", "Recall"), ("F1", "F1-score")]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.0), sharey=True)
    width = 0.2
    for ax, (metric, title) in zip(axes, metrics):
        for j, size in enumerate(LISTS):
            vals = [d[(d.tool == t) & (d.truthset == size)][metric].iloc[0] for t in TOOLS]
            xs = [i + (j - 1.5) * width for i in range(len(TOOLS))]
            ax.bar(xs, vals, width=width, color=LIST_COLORS[size], label=f"{size}-fusion list" if metric == "precision" else None)
            for x, v in zip(xs, vals):
                ax.text(x, v + 0.012, f"{v:.2f}", ha="center", va="bottom", fontsize=6.5)
            if size == 99 and not published.empty:
                for x, t in zip(xs, TOOLS):
                    for v in published.loc[published.tool.eq(t), metric].dropna():
                        ax.hlines(v, x - width * 0.46, x + width * 0.46, colors="black", linewidth=2.0, zorder=5)
        ax.set_xticks(range(len(TOOLS)), TOOLS)
        ax.set_title(title, fontweight="bold")
        ax.set_ylim(0, 1.05)
        ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("Metric value")
    handles, labels = axes[0].get_legend_handles_labels()
    if not published.empty:
        handles.append(plt.Line2D([], [], color="black", linewidth=2.0))
        labels.append("Published Edgren TP99 value(s)")
    fig.legend(handles, labels, frameon=False, ncol=5, loc="upper center", bbox_to_anchor=(0.5, 1.06))
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "edgren_metrics_by_list.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "edgren_metrics_by_list.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(TRUTH_DIR / "manifest.csv")
    runs = manifest[["sra_run", "cell_line"]].drop_duplicates().sort_values("sra_run")
    ena = pd.read_csv(PRED_DIR / "ena_runs.tsv", sep="\t")
    sra = runs.merge(ena[["run_accession", "experiment_accession", "sample_accession", "instrument_model",
                          "library_layout", "read_count", "base_count"]],
                     left_on="sra_run", right_on="run_accession", how="left").drop(columns="run_accession")
    sra["read_length_bp"] = sra.base_count / sra.read_count / (2 if (sra.library_layout == "PAIRED").all() else 1)
    sra["events_27_38_53_99"] = [
        "/".join(str(int(manifest[(manifest.sra_run == r) & (manifest.truthset_size == s)].events_for_cell_line.iloc[0]))
                 for s in LISTS) for r in sra.sra_run]
    sra.to_csv(OUTPUT_DIR / "sra_runs.csv", index=False)

    cell_runs = runs.groupby("cell_line").sra_run.apply(list).to_dict()
    per_run, per_cell, pooled = [], [], []
    for version, tool, size in product(VERSIONS, TOOLS, LISTS):
        for run, cell in zip(runs.sra_run, runs.cell_line):
            per_run += [dict(r, sra_run=run, cell_line=cell) for r in score(version, tool, [run], [run], size)]
        for cell, cruns in cell_runs.items():
            per_cell += [dict(r, cell_line=cell, runs="+".join(cruns)) for r in score(version, tool, cruns, cruns[:1], size)]
        # pooled scoring as in edgren_rnafusion_validation.py: all runs' calls against the whole list (no cell-line restriction)
        one_run_per_cell = [r[0] for r in cell_runs.values()]
        pooled += [dict(r, pooling="all runs vs whole list (no cell-line restriction)")
                   for r in score(version, tool, list(runs.sra_run), one_run_per_cell, size)]

    cols = ["rnafusion", "tool", "truthset", "matching"]
    per_run = pd.DataFrame(per_run)[cols + ["sra_run", "cell_line", "n_reference", "unique_predictions", "TP", "FP", "FN", "precision", "recall", "F1"]]
    per_cell = pd.DataFrame(per_cell)[cols + ["cell_line", "runs", "n_reference", "unique_predictions", "TP", "FP", "FN", "precision", "recall", "F1"]]
    restricted = micro(per_cell, cols).assign(pooling="sum over 4 cell lines (each against its own events)")
    pooled = pd.concat([restricted, pd.DataFrame(pooled)], ignore_index=True)[
        cols + ["pooling", "n_reference", "unique_predictions", "TP", "FP", "FN", "precision", "recall", "F1"]]
    per_run.to_csv(OUTPUT_DIR / "per_run_scores.csv", index=False)
    per_cell.to_csv(OUTPUT_DIR / "per_cell_line_scores.csv", index=False)
    pooled.to_csv(OUTPUT_DIR / "pooled_scores.csv", index=False)

    # Sanity check: edgren_rnafusion_validation.py pools the RNAfusion 3.0.1 calls against the 99 list.
    check = pooled[pooled.rnafusion.eq("3.0.1") & pooled.truthset.eq(99) & pooled.pooling.str.startswith("all runs")]
    pooled = pd.read_csv(POOLED_REFERENCE)
    comp = check.merge(pooled, on=["tool", "matching"], suffixes=("_here", "_pooled"))
    comp = comp[["tool", "matching"] + [f"{m}_{s}" for m in ["unique_predictions", "TP", "FP", "FN"] for s in ("here", "pooled")]]
    comp.to_csv(OUTPUT_DIR / "sanity_check_vs_pooled_scoring.csv", index=False)

    plot(per_run)
    plot_by_list(pd.read_csv(OUTPUT_DIR / "pooled_scores.csv"))
    summary = (per_run[per_run.matching.eq("undirected")]
               .groupby(["rnafusion", "tool", "truthset"])
               .agg(F1_min=("F1", "min"), F1_max=("F1", "max"), recall_min=("recall", "min"), recall_max=("recall", "max"))
               .reset_index())
    summary.to_csv(OUTPUT_DIR / "per_run_ranges_undirected.csv", index=False)
    print(sra.to_string(index=False))
    print("\nSanity check vs pooled scoring (3.0.1, list 99):")
    print(comp.to_string(index=False))
    print("\nPer-run ranges (undirected):")
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    print(f"\nSaved results to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
