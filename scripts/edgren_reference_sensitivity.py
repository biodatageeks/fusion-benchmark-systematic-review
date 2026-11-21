"""Same predictions, different Edgren reference lists, scored per sample (R1.2, R3.3).

Only the reference list changes between scorings; predictions, gene-name normalisation,
alias handling, deduplication and the matching rule (undirected gene pair by default) are fixed.
A prediction counts as a true positive only if the fusion is in the reference list for the
same cell line. Duplicates are collapsed within a sample (or within a cell line with
--unit cell_line), never across samples.

Predictions are the per-run caller outputs in data/raw/edgren_predictions/<workflow version>/<SRA run>/
(nf-core/rnafusion 3.0.1 and 4.1.3; Arriba, FusionCatcher, STAR-Fusion). They are parsed as fusion-report
does in nf-core/rnafusion: one gene pair per output row, and Arriba rows that name several neighbouring
genes (intergenic breakpoints) contribute one gene pair per gene combination. Parsed this way, the union of
the six 3.0.1 runs reproduces the pooled workbook data/raw/edgren_rnafusion.xlsx exactly (checked on every
run). Reference events are data/edgren/reference_events.csv and are checked against the per-run reference
sets in data/raw/edgren_truthsets/. Inputs are described in data/edgren/README.md.

With --legacy the script scores the pooled predictions in data/raw/edgren_rnafusion.xlsx against the
99-fusion list stored there (the submitted analysis); --pooled scores them against every published list.

Outputs (results/revision/edgren_reference_sensitivity/):
  metrics_by_sample.csv / metrics_by_cell_line.csv  TP/FP/FN/precision/recall/F1 per version x tool x list x unit
  metrics_aggregated_<unit>.csv  micro (summed counts) and macro (mean of per-unit metrics) aggregates
  rankings_<unit>.csv, kendall_tau_<unit>.csv  tool rank by micro F1 under each list, per workflow version
  arriba_label_sensitivity.csv   Arriba multi-gene labels counted as one prediction with alternative genes
  raw_vs_workbook_check.csv      3.0.1 per-run outputs (pooled) versus the submitted pooled workbook
  run_agreement.csv              gene pairs per run shared by the 3.0.1 and 4.1.3 outputs, per caller
  reference_check.csv            reference_events.csv versus data/raw/edgren_truthsets (event level)
  reference_overlap.csv          pairwise overlap and nesting of the reference lists (event level, per cell line)
"""

from __future__ import annotations

import argparse
import itertools
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from edgren_rnafusion_validation import INPUT as LEGACY_INPUT
from edgren_rnafusion_validation import (TOOL_COLUMNS, evaluate_tool, fusion_to_pair, gene_tokens, pair_key,
                                         pairs_match, undirected_pair_key)
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


RAW = REPO / "data" / "raw"
PRED = RAW / "edgren_predictions"
TRUTHSETS = RAW / "edgren_truthsets"
VERSIONS = ["4.1.3", "3.0.1"]
LISTS = ["TP27", "TP53", "TP99"]


def _arriba_genes(value: object) -> list[str]:
    """Arriba gene field: 'A' or 'A(123),B(456)' for intergenic breakpoints -> ['A', 'B']."""
    genes = [re.sub(r"\(\d+\)$", "", g.strip()) for g in str(value).split(",")]
    return [g for g in genes if g and g != "."]


def raw_fusions(version: str, run: str, tool: str, arriba_labels: str = "expand") -> list[str]:
    """Gene-pair labels (GENE1--GENE2) of one caller output, one per row (duplicates kept)."""
    d = pd.read_csv(PRED / version / run / f"{tool}.csv", dtype=str)
    if tool == "Arriba":
        out = []
        for a, b in zip(d["#gene1"], d["gene2"]):
            left, right = _arriba_genes(a), _arriba_genes(b)
            if arriba_labels == "expand":
                out += [f"{x}--{y}" for x in left for y in right]
            else:  # one prediction, alternative genes on each side
                out.append(f"{'|'.join(left)}--{'|'.join(right)}")
        return out
    if tool == "FusionCatcher":
        return [f"{a}--{b}" for a, b in zip(d.iloc[:, 0], d.iloc[:, 1])]
    return d["#FusionName"].tolist()


def predictions_from_fusions(fusions: list[str], tool: str) -> pd.DataFrame:
    return predictions_from_table(pd.DataFrame({"Fusion": fusions, "call": "x"}), tool, "call")


def sample_manifest() -> pd.DataFrame:
    m = pd.read_csv(DATA / "sample_manifest.csv")
    missing = [(v, r, t) for v in VERSIONS for r in m.run_accession for t in TOOL_COLUMNS
               if not (PRED / v / r / f"{t}.csv").exists()]
    if missing:
        raise SystemExit(f"Missing prediction files, e.g. {missing[:3]}")
    return m


def check_raw_against_workbook(m: pd.DataFrame) -> pd.DataFrame:
    """Pooled 3.0.1 per-run outputs versus the pooled workbook used in the submitted analysis."""
    table = pd.read_excel(LEGACY_INPUT, sheet_name="All")
    rows = []
    for tool, col in TOOL_COLUMNS.items():
        wb = set(predictions_from_table(table, tool, col)["pred_pair"].map(undirected_pair_key))
        raw = set(predictions_from_fusions(sum((raw_fusions("3.0.1", r, tool) for r in m.run_accession), []), tool)
                  ["pred_pair"].map(undirected_pair_key))
        rows.append({"tool": tool, "workbook_pairs": len(wb), "raw_pooled_pairs": len(raw), "shared": len(wb & raw),
                     "only_workbook": len(wb - raw), "only_raw": len(raw - wb), "identical": wb == raw})
    return pd.DataFrame(rows)


def check_reference(ref: pd.DataFrame) -> pd.DataFrame:
    """Event-level comparison of reference_events.csv with the per-run reference sets of the authors."""
    sets = pd.concat([pd.read_csv(f) for f in sorted(TRUTHSETS.glob("truthset_*/*.csv"))])
    sets = sets.drop_duplicates(["truthset_size", "cell_line", "gene1_raw", "gene2_raw"])
    rows = []
    for size, list_id in [(27, "TP27"), (53, "TP53"), (99, "TP99")]:
        theirs = sets[sets.truthset_size == size]
        ours = ref[ref.list_id == list_id]
        for _, r in theirs.iterrows():
            pair = (gene_tokens(r.gene1_raw), gene_tokens(r.gene2_raw))
            hit = ours[(ours.cell_line == r.cell_line) & ours.truth_pair.map(lambda q: pairs_match(pair, q, directional=False))]
            rows.append({"list_id": list_id, "cell_line": r.cell_line, "truthset_event": f"{r.gene1_raw}--{r.gene2_raw}",
                         "reference_event": "" if hit.empty else f"{hit.gene1.iloc[0]}--{hit.gene2.iloc[0]}",
                         "matched": not hit.empty})
        for _, q in ours.iterrows():
            found = theirs[(theirs.cell_line == q.cell_line).to_numpy() & np.array([pairs_match((gene_tokens(a), gene_tokens(b)), q.truth_pair, False)
                                                                for a, b in zip(theirs.gene1_raw, theirs.gene2_raw)], dtype=bool)]
            if found.empty:
                rows.append({"list_id": list_id, "cell_line": q.cell_line, "truthset_event": "",
                             "reference_event": f"{q.gene1}--{q.gene2}", "matched": False})
    return pd.DataFrame(rows)


def run_agreement(m: pd.DataFrame) -> pd.DataFrame:
    """Agreement of the unique gene pairs reported per run by the original (3.0.1) and documented (4.1.3) runs."""
    rows = []
    for tool in TOOL_COLUMNS:
        a_all = b_all = shared = 0
        for run in m.run_accession:
            a = set(predictions_from_fusions(raw_fusions("3.0.1", run, tool), tool)["pred_pair"].map(undirected_pair_key))
            b = set(predictions_from_fusions(raw_fusions("4.1.3", run, tool), tool)["pred_pair"].map(undirected_pair_key))
            a_all, b_all, shared = a_all + len(a), b_all + len(b), shared + len(a & b)
        rows.append({"tool": tool, "pairs_3.0.1": a_all, "pairs_4.1.3": b_all, "shared": shared,
                     "jaccard": shared / (a_all + b_all - shared), "share_of_3.0.1_in_4.1.3": shared / a_all,
                     "share_of_4.1.3_in_3.0.1": shared / b_all})
    return pd.DataFrame(rows)


def per_unit_metrics(m: pd.DataFrame, ref: pd.DataFrame, unit: str, arriba_labels: str = "expand") -> pd.DataFrame:
    rows = []
    for version in VERSIONS:
        groups = m.groupby("run_accession" if unit == "sample" else "cell_line", sort=False)
        for unit_id, g in groups:
            cell_line = g.cell_line.iloc[0]
            for tool in TOOL_COLUMNS:
                fusions = sum((raw_fusions(version, r, tool, arriba_labels) for r in g.run_accession), [])
                pred = predictions_from_fusions(fusions, tool)
                for list_id in LISTS:
                    truth = ref[(ref.list_id == list_id) & (ref.cell_line == cell_line)]
                    rows.append({"workflow_version": version, "unit": unit_id, "cell_line": cell_line, "tool": tool,
                                 "list_id": list_id, "n_reference_events": len(truth), "unique_predictions": len(pred),
                                 **score(pred, truth)})
    return pd.DataFrame(rows)


def aggregate(by_unit: pd.DataFrame) -> pd.DataFrame:
    keys = ["workflow_version", "tool", "list_id"]
    micro = by_unit.groupby(keys, sort=False)[["TP", "FP", "FN"]].sum().reset_index()
    micro["precision"] = micro.TP / (micro.TP + micro.FP)
    micro["recall"] = micro.TP / (micro.TP + micro.FN)
    micro["F1"] = 2 * micro.precision * micro.recall / (micro.precision + micro.recall)
    macro = by_unit.groupby(keys, sort=False)[["precision", "recall", "F1"]].mean().reset_index()
    return micro.merge(macro.rename(columns={c: f"macro_{c}" for c in ["precision", "recall", "F1"]}), on=keys)


def overlap_table(ref: pd.DataFrame) -> pd.DataFrame:
    rows = []
    sets = {lid: set(zip(g["cell_line"].fillna("NA"), g["undirected_key"])) for lid, g in ref.groupby("list_id")}
    for a, b in itertools.permutations(sorted(sets), 2):
        rows.append({"list_a": a, "list_b": b, "n_a": len(sets[a]), "n_b": len(sets[b]),
                     "shared": len(sets[a] & sets[b]), "only_in_a": len(sets[a] - sets[b]),
                     "a_nested_in_b": sets[a] <= sets[b]})
    return pd.DataFrame(rows)


def write_supplementary_tables(ref: pd.DataFrame) -> None:
    """Supplementary Table S12 (event-level correspondence of the published lists); S11 comes from the per-sample analysis."""
    from revision_utils import write_latex_table

    tables = REVISION_DIR / "tables"
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
    ref = ref[ref.list_id.isin(LISTS)].reset_index(drop=True)
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
    fig.savefig(REVISION_DIR / "figures" / "edgren_rnafusion_reference_lists_pooled.pdf", bbox_inches="tight")
    fig.savefig(REVISION_DIR / "figures" / "edgren_rnafusion_reference_lists_pooled.png", dpi=200, bbox_inches="tight")
    print(res.to_string(index=False))
    print(ranks.to_string())
    return 0


def published_edgren_metrics() -> pd.DataFrame:
    """Values published by the included benchmarks for the same callers on the Edgren data, per reference list.

    Only observations with a stated list size (27, 53 or 99) are kept; reported values are used as published
    (Liu et al. 2016 reported F1 only).
    """
    from revision_utils import load_observations

    df = load_observations()
    pub = df[(df.bio_dataset_id == "EDGREN_BRCA_CELL_LINES") & df.TruthTotal.isin([27, 53, 99])
             & df.Tool.isin(TOOL_COLUMNS)].copy()
    pub["list_id"] = "TP" + pub.TruthTotal.astype(int).astype(str)
    pub = pub.rename(columns={"Tool": "tool", "Precision": "precision", "Recall": "recall", "Benchmark_label": "benchmark"})
    pub = pub[["benchmark", "Dataset", "tool", "list_id", "precision", "recall", "F1"]].sort_values(["list_id", "tool", "benchmark"])
    pub.to_csv(OUT / "published_edgren_metrics.csv", index=False)
    return pub.reset_index(drop=True)


def plot_per_sample(agg: pd.DataFrame, published: pd.DataFrame | None = None) -> None:
    """Figure 7: per-sample (micro-averaged) metrics of the same predictions against three reference lists.

    Black markers: values published for the same caller and list (pooled scoring of other runs and versions).
    """
    import os
    os.environ.setdefault("MPLBACKEND", "Agg")
    import matplotlib.pyplot as plt

    tools = list(TOOL_COLUMNS)
    colors = {"TP27": "#9ecae1", "TP53": "#4292c6", "TP99": "#08519c"}
    n_ref = {"TP27": 27, "TP53": 53, "TP99": 99}
    fig, axes = plt.subplots(len(VERSIONS), 3, figsize=(12, 7.6), sharey=True)
    for r, version in enumerate(VERSIONS):
        sub = agg[agg.workflow_version == version]
        for c, metric in enumerate(["precision", "recall", "F1"]):
            ax = axes[r, c]
            marked = False
            for j, lid in enumerate(LISTS):
                vals = [sub[(sub.tool == t) & (sub.list_id == lid)][metric].iloc[0] for t in tools]
                xs = [k + (j - 1) * 0.26 for k in range(len(tools))]
                ax.bar(xs, vals, width=0.26, color=colors[lid], label=f"{n_ref[lid]}-fusion list")
                for x, v in zip(xs, vals):
                    ax.text(x, v + 0.015, f"{v:.2f}", ha="center", fontsize=7)
                if published is None:
                    continue
                for x, t in zip(xs, tools):
                    for v in published[(published.tool == t) & (published.list_id == lid)][metric].dropna():
                        ax.hlines(v, x - 0.12, x + 0.12, colors="black", linewidth=2.2, zorder=5,
                                  label="Published value, same caller and list" if not marked else None)
                        marked = True
            ax.set_xticks(range(len(tools)), tools)
            ax.set_ylim(0, 1.0)
            ax.spines[["top", "right"]].set_visible(False)
            if r == 0:
                ax.set_title(metric if metric == "F1" else metric.capitalize(), fontweight="bold")
        axes[r, 0].set_ylabel(f"nf-core/rnafusion {version}\n(micro-averaged over 6 runs)")
    axes[0, 2].legend(frameon=False, fontsize=8, loc="upper right")
    fig.suptitle("Same RNAfusion predictions on the six Edgren runs, scored per sample against three published reference lists",
                 fontsize=11)
    fig.tight_layout()
    out = REVISION_DIR / "figures"
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "edgren_rnafusion_reference_lists.pdf", bbox_inches="tight")
    fig.savefig(out / "edgren_rnafusion_reference_lists.png", dpi=200, bbox_inches="tight")


def write_table_s11(m: pd.DataFrame, ref: pd.DataFrame, by_sample: pd.DataFrame, agg: pd.DataFrame,
                    pooled: pd.DataFrame | None) -> None:
    """Supplementary Table S11: (a) the six runs with predictions per caller, (b) per-sample scoring."""
    from revision_utils import latex_escape

    counts = ref.groupby(["cell_line", "list_id"]).size().unstack(fill_value=0)
    npred = by_sample.drop_duplicates(["workflow_version", "unit", "tool"]).set_index(["workflow_version", "unit", "tool"]).unique_predictions
    a_rows = []
    for _, r in m.iterrows():
        refs = " / ".join(str(int(counts.loc[r.cell_line, l])) for l in LISTS)
        preds = [f"{npred[('4.1.3', r.run_accession, t)]} / {npred[('3.0.1', r.run_accession, t)]}" for t in TOOL_COLUMNS]
        a_rows.append([r.run_accession, r.cell_line, f"{r.read_layout.lower()}, {r.read_length_bp} bp", f"{int(r.read_pairs):,}", refs, *preds])
    a_head = ["SRA run", "Cell line", "Layout", "Read pairs", "Reference fusions (27 / 53 / 99)",
              *[f"{t} pairs (4.1.3 / 3.0.1)" for t in TOOL_COLUMNS]]
    b_rows = []
    for tool in TOOL_COLUMNS:
        for lid in LISTS:
            row = [tool, lid.replace("TP", "") + " fusions"]
            for version in VERSIONS:
                s = agg[(agg.workflow_version == version) & (agg.tool == tool) & (agg.list_id == lid)].iloc[0]
                row += [f"{int(s.TP)} / {int(s.FP)} / {int(s.FN)}", f"{s.precision:.2f}", f"{s.recall:.2f}", f"{s.F1:.2f}"]
            if pooled is not None:
                q = pooled[(pooled.tool == tool) & (pooled.list_id == lid)].iloc[0]
                row.append(f"{q.F1:.2f}")
            b_rows.append(row)
    b_head = ["Caller", "Reference list"] + [f"{h} ({v})" for v in VERSIONS for h in ["TP / FP / FN", "Precision", "Recall", "F1"]]
    if pooled is not None:
        b_head.append("F1, pooled (3.0.1)")

    def tabular(head, rows, spec):
        body = "\n".join(" & ".join(latex_escape(v) for v in row) + r" \\" for row in rows)
        return (f"\\begin{{tabular}}{{{spec}}}\n\\toprule\n" + " & ".join(latex_escape(h) for h in head) + r" \\" +
                f"\n\\midrule\n{body}\n\\bottomrule\n\\end{{tabular}}")

    caption = (r"\textbf{Edgren RNA-seq runs and per-sample re-scoring of the same RNAfusion predictions} (SRA project SRP003186). "
               r"(a) The six tumour cell-line runs (the normal breast run SRR064437 is excluded), reference fusions for the cell line "
               r"in each published list, and unique gene pairs reported per run by each caller in the documented re-run "
               r"(nf-core/rnafusion 4.1.3) and the recovered original run (3.0.1). (b) Each run scored against the reference "
               r"list of its own cell line; counts summed over the six runs (micro-average). The last column is the submitted scoring of the 3.0.1 predictions "
               r"pooled over runs and cell lines, in which a call repeated in several runs counted once; no call matched only a reference "
               r"fusion of another cell line. "
               r"Undirected gene-pair matching; the 37- and 80-fusion definitions were not published as gene lists.")
    text = ("\\begin{table}[!htbp]\n\\centering\n\\footnotesize\n" + f"\\caption{{{caption}}}\n\\label{{tab:supp_edgren_runs}}\n"
            + "\\textbf{(a) Sequencing runs}\\\\[2pt]\n\\resizebox{\\linewidth}{!}{%\n" + tabular(a_head, a_rows, "lllrlccc") + "}\\\\[8pt]\n"
            + "\\textbf{(b) Per-sample scoring (micro-average over six runs)}\\\\[2pt]\n\\resizebox{\\linewidth}{!}{%\n"
            + tabular(b_head, b_rows, "ll" + "lccc" * len(VERSIONS) + ("c" if pooled is not None else "")) + "}\n\\end{table}\n")
    (REVISION_DIR / "tables" / "supp_edgren_runs.tex").write_text(text)


def per_sample_analysis() -> int:
    m = sample_manifest()
    ref = reference_table(pd.read_csv(DATA / "reference_events.csv"))
    ref = ref[ref.list_id.isin(LISTS)].reset_index(drop=True)
    overlap_table(ref).to_csv(OUT / "reference_overlap.csv", index=False)
    refcheck = check_reference(ref)
    refcheck.to_csv(OUT / "reference_check.csv", index=False)
    rawcheck = check_raw_against_workbook(m)
    rawcheck.to_csv(OUT / "raw_vs_workbook_check.csv", index=False)
    print(rawcheck.to_string(index=False))
    agreement = run_agreement(m)
    agreement.to_csv(OUT / "run_agreement.csv", index=False)
    print(agreement.round(2).to_string(index=False))
    print("reference events not matched:", refcheck[~refcheck.matched].to_dict("records"))

    results = {}
    for unit in ["sample", "cell_line"]:
        by_unit = per_unit_metrics(m, ref, unit)
        by_unit.to_csv(OUT / f"metrics_by_{unit}.csv", index=False)
        agg = aggregate(by_unit)
        agg.to_csv(OUT / f"metrics_aggregated_{unit}.csv", index=False)
        ranks = agg.pivot_table(index=["workflow_version", "tool"], columns="list_id", values="F1").groupby(level=0).rank(ascending=False)
        ranks.to_csv(OUT / f"rankings_{unit}.csv")
        taus = [{"workflow_version": v, "list_a": a, "list_b": b,
                 "kendall_tau": stats.kendalltau(ranks.loc[v, a], ranks.loc[v, b]).statistic}
                for v in VERSIONS for a, b in itertools.combinations(LISTS, 2)]
        pd.DataFrame(taus).to_csv(OUT / f"kendall_tau_{unit}.csv", index=False)
        results[unit] = (by_unit, agg)
        print(f"--- unit = {unit}")
        print(agg.round(3).to_string(index=False))
        print(ranks.to_string())

    alt = aggregate(per_unit_metrics(m, ref, "sample", arriba_labels="alternatives"))
    alt = alt[alt.tool == "Arriba"].merge(results["sample"][1][["workflow_version", "tool", "list_id", "F1", "precision", "FP"]],
                                          on=["workflow_version", "tool", "list_id"], suffixes=("_alternatives", "_expanded"))
    alt.to_csv(OUT / "arriba_label_sensitivity.csv", index=False)
    print(alt[["workflow_version", "list_id", "FP_expanded", "FP_alternatives", "F1_expanded", "F1_alternatives"]].round(3).to_string(index=False))

    pooled_file = OUT / "pooled_metrics_by_list.csv"
    pooled = pd.read_csv(pooled_file) if pooled_file.exists() else None
    by_sample, agg = results["sample"]
    plot_per_sample(agg, published_edgren_metrics())
    write_table_s11(m, ref, by_sample, agg, pooled)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--legacy", action="store_true", help="pooled 99-fusion check on the submitted data")
    parser.add_argument("--pooled", action="store_true",
                        help="score the pooled original RNAfusion predictions against every published list "
                             "(cell lines ignored) and write Supplementary Table S12")
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

    return per_sample_analysis()


if __name__ == "__main__":
    raise SystemExit(main())
