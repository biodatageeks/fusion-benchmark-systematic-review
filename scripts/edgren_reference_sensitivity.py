"""Same predictions, different Edgren reference lists, scored per sample (R1.2, R3.3).

Only the reference list changes between scorings; predictions, gene-name normalisation,
alias handling, deduplication and the matching rule (undirected gene pair by default) are fixed.
A prediction counts as a true positive only if the fusion is in the reference list for the
same cell line. Duplicates are collapsed within a sample (or within a cell line with
--unit cell_line), never across samples.

Inputs are described in data/edgren/README.md. With --legacy the script scores the pooled
predictions in data/raw/edgren_rnafusion.xlsx against the 99-fusion list stored there; this
reproduces the submitted analysis and checks the implementation.

Outputs (results/revision/edgren_reference_sensitivity/):
  metrics_by_unit.csv     TP/FP/FN/precision/recall/F1 per tool x reference list x sample (or cell line)
  metrics_aggregated.csv  micro (summed counts) and macro (mean of per-unit metrics) aggregates
  rankings.csv            tool rank by F1 under each reference list; kendall_tau.csv between lists
  reference_overlap.csv   pairwise overlap and nesting of the reference lists (event level, per cell line)
"""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from edgren_rnafusion_validation import INPUT as LEGACY_INPUT
from edgren_rnafusion_validation import TOOL_COLUMNS, evaluate_tool, gene_tokens, fusion_to_pair, pair_key
from revision_utils import REPO, REVISION_DIR

DATA = REPO / "data" / "edgren"
OUT = REVISION_DIR / "edgren_reference_sensitivity"


def predictions_from_table(table: pd.DataFrame, tool: str, column: str) -> pd.DataFrame:
    p = table.loc[table[column].notna() & table[column].astype(str).str.strip().ne(""), ["Fusion"]].copy()
    p["pred_pair"] = p["Fusion"].map(fusion_to_pair)
    p = p.dropna(subset=["pred_pair"])
    p["pred_key"] = p["pred_pair"].map(pair_key)
    p = p.drop_duplicates("pred_key").reset_index(drop=True)
    p["tool"] = tool
    return p[["tool", "Fusion", "pred_pair", "pred_key"]]


def reference_table(events: pd.DataFrame) -> pd.DataFrame:
    t = events.dropna(subset=["gene1", "gene2"]).copy()
    t["truth_pair"] = list(zip(t["gene1"].map(gene_tokens), t["gene2"].map(gene_tokens)))
    t["undirected_key"] = t["truth_pair"].map(lambda p: tuple(sorted(pair_key(p))))
    return t.drop_duplicates(["list_id", "cell_line", "undirected_key"]).reset_index(drop=True)


def score(pred: pd.DataFrame, truth: pd.DataFrame) -> dict[str, object]:
    if truth.empty:
        return {"TP": 0, "FP": len(pred), "FN": 0, "precision": np.nan, "recall": np.nan, "F1": np.nan}
    if pred.empty:
        return {"TP": 0, "FP": 0, "FN": len(truth), "precision": np.nan, "recall": 0.0, "F1": np.nan}
    summary, *_ = evaluate_tool(pred, truth.reset_index(drop=True), directional=False)
    return {k: summary[k] for k in ["TP", "FP", "FN", "precision", "recall", "F1"]}


def load_units(unit: str) -> tuple[dict[str, pd.DataFrame], dict[str, str]]:
    manifest = pd.read_csv(DATA / "sample_manifest.csv")
    tables, cell = {}, {}
    for _, row in manifest.iterrows():
        t = pd.read_csv(DATA / row["prediction_file"])
        key = row["sample_id"] if unit == "sample" else row["cell_line"]
        tables[key] = pd.concat([tables[key], t]) if key in tables else t
        cell[key] = row["cell_line"]
    return tables, cell


def overlap_table(ref: pd.DataFrame) -> pd.DataFrame:
    rows = []
    sets = {lid: set(zip(g["cell_line"].fillna("NA"), g["undirected_key"])) for lid, g in ref.groupby("list_id")}
    for a, b in itertools.permutations(sorted(sets), 2):
        rows.append({"list_a": a, "list_b": b, "n_a": len(sets[a]), "n_b": len(sets[b]),
                     "shared": len(sets[a] & sets[b]), "only_in_a": len(sets[a] - sets[b]),
                     "a_nested_in_b": sets[a] <= sets[b]})
    return pd.DataFrame(rows)


def write_supplementary_tables(ref: pd.DataFrame) -> None:
    """Supplementary Table S11 (SRA runs) and S12 (event-level correspondence of the published lists)."""
    from revision_utils import write_latex_table

    tables = REVISION_DIR / "tables"
    manifest = DATA / "sample_manifest.csv"
    if not manifest.exists():
        manifest = DATA / "sample_manifest_candidates.csv"
    if manifest.exists():
        m = pd.read_csv(manifest)
        if "likely_benchmark_run" in m.columns:
            m = m[m["likely_benchmark_run"].astype(str).str.lower().eq("yes")]
        counts = ref.groupby(["cell_line", "list_id"]).size().unstack(fill_value=0)
        m["refs"] = m["cell_line"].map(lambda c: " / ".join(str(int(counts.loc[c, l])) if c in counts.index and l in counts.columns else "0"
                                                              for l in ["TP27", "TP53", "TP99"]))
        m["pairs"] = m.get("read_pairs", pd.Series([np.nan] * len(m))).map(lambda v: "--" if pd.isna(v) else f"{int(v):,}")
        s11 = m[["run_accession", "cell_line", "read_layout", "read_length_bp", "pairs", "refs"]].rename(columns={
            "run_accession": "SRA run", "cell_line": "Cell line", "read_layout": "Layout", "read_length_bp": "Read length (bp)",
            "pairs": "Read pairs", "refs": "Reference fusions (TP27 / TP53 / TP99)"})
        write_latex_table(s11, tables / "supp_edgren_runs.tex",
                          caption=r"\textbf{Edgren RNA-seq runs used in the controlled RNAfusion re-scoring} (SRA project SRP003186). "
                                  r"The six tumour cell-line runs listed in Uhrig et al.; the normal breast run SRR064437 is excluded. "
                                  r"Reference fusions per cell line from the published lists; the 37- and 80-fusion lists were not published as gene lists.",
                          label="tab:supp_edgren_runs", colspec="llllrl")
    wide = (ref.assign(member="yes").pivot_table(index=["cell_line", "undirected_key"], columns="list_id", values="member", aggfunc="first")
            .fillna("--").reset_index())
    first = ref.drop_duplicates(["cell_line", "undirected_key"]).set_index(["cell_line", "undirected_key"])
    wide["Fusion"] = [f"{first.loc[(c, k), 'gene1']}--{first.loc[(c, k), 'gene2']}" for c, k in zip(wide.cell_line, wide.undirected_key)]
    wide["Evidence"] = [str(first.loc[(c, k), "evidence"]).split(";")[0][:45] for c, k in zip(wide.cell_line, wide.undirected_key)]
    lists = [c for c in ["TP27", "TP53", "TP99"] if c in wide.columns]
    s12 = wide.sort_values(["cell_line", "Fusion"])[["cell_line", "Fusion", *lists, "Evidence"]].rename(columns={"cell_line": "Cell line"})
    sets = {l: set(zip(ref.loc[ref.list_id == l, "cell_line"], ref.loc[ref.list_id == l, "undirected_key"])) for l in lists}
    caption = (r"\textbf{Event-level correspondence between published Edgren reference lists.} One row per fusion and cell line "
               rf"({len(wide)} distinct events). TP27 is nested in TP53 and in TP99; TP53 and TP99 share {len(sets['TP53'] & sets['TP99'])} events, "
               rf"{len(sets['TP53'] - sets['TP99'])} are only in TP53 and {len(sets['TP99'] - sets['TP53'])} only in TP99. The 37-fusion "
               r"(Zhao et al.) and 80-fusion (Uhrig et al.) definitions were not published as gene lists. Sources and evidence per event: "
               r"\texttt{data/edgren/reference\_events.csv}.")
    write_latex_table(s12, tables / "supp_edgren_events.tex", caption=caption, label="tab:supp_edgren_events",
                      colspec="lp{0.30\\linewidth}cccp{0.30\\linewidth}", font=r"\scriptsize", longtable=True)


def pooled_multi_reference() -> int:
    """Original pooled predictions (one table for all six runs) scored against each published list."""
    import os
    os.environ.setdefault("MPLBACKEND", "Agg")
    import matplotlib.pyplot as plt

    table = pd.read_excel(LEGACY_INPUT, sheet_name="All")
    ref = reference_table(pd.read_csv(DATA / "reference_events.csv"))
    overlap_table(ref).to_csv(OUT / "reference_overlap.csv", index=False)
    write_supplementary_tables(ref)
    pooled_ref = ref.drop_duplicates(["list_id", "undirected_key"])
    rows = []
    for tool, col in TOOL_COLUMNS.items():
        pred = predictions_from_table(table, tool, col)
        for list_id, lref in pooled_ref.groupby("list_id"):
            rows.append({"tool": tool, "list_id": list_id, "n_reference_events": len(lref), "unique_predictions": len(pred),
                         **score(pred, lref)})
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "pooled_metrics_by_list.csv", index=False)
    ranks = res.pivot(index="tool", columns="list_id", values="F1").rank(ascending=False)
    ranks.to_csv(OUT / "pooled_rankings.csv")
    lists = sorted(res.list_id.unique(), key=lambda s: int(s[2:]))
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.2), sharey=True)
    colors = {"TP27": "#9ecae1", "TP53": "#4292c6", "TP99": "#08519c"}
    tools = list(TOOL_COLUMNS)
    for ax, metric in zip(axes, ["precision", "recall", "F1"]):
        for j, lid in enumerate(lists):
            vals = [res[(res.tool == t) & (res.list_id == lid)][metric].iloc[0] for t in tools]
            xs = [i + (j - 1) * 0.26 for i in range(len(tools))]
            ax.bar(xs, vals, width=0.26, color=colors.get(lid, "#777"), label=f"{lid} ({int(pooled_ref[pooled_ref.list_id == lid].shape[0])} fusions)")
            for x, v in zip(xs, vals):
                ax.text(x, v + 0.015, f"{v:.2f}", ha="center", fontsize=7.5)
        ax.set_xticks(range(len(tools)), tools)
        ax.set_title(metric if metric == "F1" else metric.capitalize(), fontweight="bold")
        ax.set_ylim(0, 1.0)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Value (same predictions)")
    axes[2].legend(frameon=False, fontsize=8, loc="upper right")
    fig.suptitle("Same RNAfusion predictions on the Edgren data, scored against three published reference lists", fontsize=11)
    fig.tight_layout()
    (REVISION_DIR / "figures").mkdir(parents=True, exist_ok=True)
    fig.savefig(REVISION_DIR / "figures" / "edgren_rnafusion_reference_lists.pdf", bbox_inches="tight")
    fig.savefig(REVISION_DIR / "figures" / "edgren_rnafusion_reference_lists.png", dpi=200, bbox_inches="tight")
    print(res.to_string(index=False))
    print(ranks.to_string())
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--unit", choices=["sample", "cell_line"], default="sample")
    parser.add_argument("--legacy", action="store_true", help="pooled 99-fusion check on the submitted data")
    parser.add_argument("--pooled", action="store_true",
                        help="score the pooled original RNAfusion predictions against every published list "
                             "(cell lines ignored), and write Supplementary Tables S11-S12")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    if args.legacy:
        table = pd.read_excel(LEGACY_INPUT, sheet_name="All")
        ref = reference_table(table[["Gene1", "Gene2"]].rename(columns={"Gene1": "gene1", "Gene2": "gene2"})
                              .assign(list_id="TP99", cell_line=np.nan))
        rows = [{"tool": tool, "list_id": "TP99", "unit": "pooled",
                 **score(predictions_from_table(table, tool, col), ref)} for tool, col in TOOL_COLUMNS.items()]
        res = pd.DataFrame(rows)
        res.to_csv(OUT / "legacy_pooled_tp99.csv", index=False)
        print(res.to_string(index=False))
        return 0

    if args.pooled:
        return pooled_multi_reference()

    for needed in ["sample_manifest.csv", "reference_events.csv"]:
        if not (DATA / needed).exists():
            raise SystemExit(f"Missing data/edgren/{needed}; see data/edgren/README.md (templates provided).")

    ref = reference_table(pd.read_csv(DATA / "reference_events.csv"))
    overlap_table(ref).to_csv(OUT / "reference_overlap.csv", index=False)
    tables, cell = load_units(args.unit)
    rows = []
    for unit_id, table in tables.items():
        for tool, col in TOOL_COLUMNS.items():
            if col not in table.columns:
                continue
            pred = predictions_from_table(table, tool, col)
            for list_id, lref in ref.groupby("list_id"):
                truth = lref[lref["cell_line"] == cell[unit_id]]
                rows.append({"unit": unit_id, "cell_line": cell[unit_id], "tool": tool, "list_id": list_id,
                             "n_reference_events": len(truth), **score(pred, truth)})
    by_unit = pd.DataFrame(rows)
    by_unit.to_csv(OUT / f"metrics_by_{args.unit}.csv", index=False)

    micro = by_unit.groupby(["tool", "list_id"])[["TP", "FP", "FN"]].sum().reset_index()
    micro["precision"] = micro.TP / (micro.TP + micro.FP)
    micro["recall"] = micro.TP / (micro.TP + micro.FN)
    micro["F1"] = 2 * micro.precision * micro.recall / (micro.precision + micro.recall)
    macro = by_unit.groupby(["tool", "list_id"])[["precision", "recall", "F1"]].mean().reset_index()
    agg = micro.assign(aggregation="micro").merge(
        macro.rename(columns={c: f"macro_{c}" for c in ["precision", "recall", "F1"]}), on=["tool", "list_id"])
    agg.to_csv(OUT / "metrics_aggregated.csv", index=False)

    ranks = agg.pivot(index="tool", columns="list_id", values="F1").rank(ascending=False)
    ranks.to_csv(OUT / "rankings.csv")
    taus = [{"list_a": a, "list_b": b, "kendall_tau": stats.kendalltau(ranks[a], ranks[b]).statistic}
            for a, b in itertools.combinations(ranks.columns, 2)]
    pd.DataFrame(taus).to_csv(OUT / "kendall_tau.csv", index=False)
    print(agg.to_string(index=False))
    print(ranks.to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
