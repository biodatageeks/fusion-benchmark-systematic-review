"""Supplementary tables built from the analysis workbook and the Edgren reference lists.

Outputs (results/supplementary_tables/):
  S3_sra_runs.csv                Edgren sequencing runs used for per-run scoring
  S4_tools.csv                   every tool or tool variant: family, benchmarks, observations
  S5_tool_profiles.csv           mean precision and recall per tool (all data; real data without Edgren)
  S6_edgren_reference_lists.csv  event-level membership of the Edgren reference lists 27, 38, 53 and 99
  S7_benchmark_contexts.csv      benchmark contexts: purpose, data, library, reads, truth set, observations
  summary.md                     counts used in the text
  latex/S3_rows.tex … S7_rows.tex  table rows for the LaTeX supplement (labels and formatting as in the manuscript)
"""
from __future__ import annotations

import glob
from pathlib import Path

import pandas as pd

from edgren_rnafusion_validation import gene_tokens
from plot_benchmark_reporting_heatmap import status, truth_type
from precision_recall_tradeoff import load_clean_data, summarize_tools

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "raw" / "master_analysis_input.xlsx"
TRUTH_DIR = ROOT / "data" / "raw" / "edgren_truthsets"
SRA_RUNS = ROOT / "results" / "edgren_per_sample" / "sra_runs.csv"
OUTPUT_DIR = ROOT / "results" / "supplementary_tables"
LISTS = [27, 38, 53, 99]


def benchmark_labels(bench: pd.DataFrame) -> dict[int, str]:
    return {int(r.benchmark_number): f"{r.first_author} {int(r.year_x)}" for r in bench.itertuples()}


def bench_order(bench: pd.DataFrame) -> list[int]:
    return [int(n) for n in bench.sort_values(["year_x", "benchmark_number"])["benchmark_number"]]


def tools_table(obs: pd.DataFrame, names: pd.DataFrame, bench: pd.DataFrame) -> pd.DataFrame:
    labels, order = benchmark_labels(bench), bench_order(bench)
    rank = {n: i for i, n in enumerate(order)}
    names = names.drop_duplicates("tool_name_clean").set_index("tool_name_clean")
    rows = []
    for tool, g in obs.groupby("tool_name_clean"):
        benches = sorted({int(b) for b in g["benchmark_number"]}, key=rank.get)
        rows.append({
            "tool": tool,
            "family": names["tool_family"].get(tool, g["tool_family"].iloc[0]),
            "variant_or_output": bool(names["is_variant_or_output"].get(tool, False)),
            "n_benchmarks": len(benches),
            "benchmarks": "; ".join(labels[b] for b in benches),
            "n_observations": len(g),
            "n_with_f1": int(pd.to_numeric(g["f1_score"], errors="coerce").notna().sum()),
        })
    return pd.DataFrame(rows).sort_values(["family", "variant_or_output", "tool"], key=lambda s: s if s.dtype == bool else s.astype(str).str.lower())


def profiles_table() -> pd.DataFrame:
    df = load_clean_data()
    cols = ["Tool", "n", "n_benchmarks", "mean_precision", "mean_recall", "profile"]
    everything = summarize_tools(df)[cols]
    real = summarize_tools(df[df["Is_real"] & ~df["Is_Edgren"]])[["Tool", "n", "mean_precision", "mean_recall"]]
    real = real.rename(columns={"n": "n_real_without_edgren", "mean_precision": "mean_precision_real_without_edgren",
                                "mean_recall": "mean_recall_real_without_edgren"})
    out = everything.merge(real, on="Tool", how="left").rename(columns={"Tool": "tool", "n": "n_observations"})
    return out.sort_values("tool", key=lambda s: s.str.lower()).round(3)


def edgren_records() -> list[dict]:
    """One record per (list, cell line, gene pair); runs of the same cell line are collapsed."""
    records: dict[tuple, dict] = {}
    for size in LISTS:
        for path in sorted(glob.glob(str(TRUTH_DIR / f"truthset_{size}" / "*.csv"))):
            for r in pd.read_csv(path).itertuples():
                a, b = gene_tokens(r.gene1_raw), gene_tokens(r.gene2_raw)
                key = (size, r.cell_line, frozenset([a, b]))
                rec = records.setdefault(key, {"list": size, "cell_line": r.cell_line, "a": a, "b": b,
                                               "fusion": f"{r.gene1_raw}::{r.gene2_raw}", "source": str(r.source),
                                               "evidence": set()})
                if pd.notna(r.validation):
                    rec["evidence"].add(str(r.validation))
    return list(records.values())


# Symbols misspelled in a source list, shown corrected in Table S6 (the scoring matches them through ALIASES):
# the 99-fusion list (Davidson 2015) has MED1--UPS32 for MED1--USP32.
MISSPELLED = {"UPS32": "USP32"}


def same_event(x: dict, y: dict) -> bool:
    """Matching rule of the scoring: same cell line, gene aliases, either orientation."""
    return x["cell_line"] == y["cell_line"] and bool(
        (x["a"] & y["a"] and x["b"] & y["b"]) or (x["a"] & y["b"] and x["b"] & y["a"]))


def edgren_lists_table() -> pd.DataFrame:
    recs = edgren_records()
    parent = list(range(len(recs)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(recs)):
        for j in range(i + 1, len(recs)):
            if same_event(recs[i], recs[j]):
                parent[find(i)] = find(j)
    groups: dict[int, list[dict]] = {}
    for i, rec in enumerate(recs):
        groups.setdefault(find(i), []).append(rec)
    rows = []
    for members in groups.values():
        lists = {m["list"] for m in members}
        shown = min(members, key=lambda m: (-m["list"], m["fusion"]))  # name as in the longest list
        fusion = "::".join(MISSPELLED.get(gene, gene) for gene in shown["fusion"].split("::"))
        row = {"cell_line": shown["cell_line"], "fusion": fusion}
        row.update({f"list_{n}": n in lists for n in LISTS})
        row["sources"] = "; ".join(sorted({m["source"] for m in members}))
        row["evidence"] = "; ".join(sorted(set().union(*(m["evidence"] for m in members))))
        row["records_merged"] = len(members)
        rows.append(row)
    out = pd.DataFrame(rows)
    out["first_list"] = out[[f"list_{n}" for n in LISTS]].idxmax(axis=1).str.replace("list_", "").astype(int)
    return out.sort_values(["first_list", "cell_line", "fusion"]).drop(columns="first_list")


def contexts_table(obs: pd.DataFrame, data: pd.DataFrame, bench: pd.DataFrame) -> pd.DataFrame:
    labels = benchmark_labels(bench)
    data = data.drop_duplicates("dataset_id").set_index("dataset_id")
    obs = obs.copy()
    obs["origin"] = obs["dataset_origin"].str.lower().str.strip()
    obs["subtype"] = obs["dataset_id"].map(data["real_subtype"])
    obs["f1"] = pd.to_numeric(obs["f1_score"], errors="coerce")

    def material(r: pd.Series) -> str:
        if r.origin == "real":
            return {"cell line": "cell line", "patient-derived": "patient", "mixed": "mixed"}.get(r.subtype, "real")
        return r.origin

    obs["material"] = obs.apply(material, axis=1)
    rows = []
    for _, b in bench.sort_values(["year_x", "benchmark_number"]).iterrows():
        g = obs[obs["benchmark_number"].astype(int) == int(b.benchmark_number)]
        lengths = [str(n) for n in (50, 75, 100, 250) if str(b[f"read lenght {n}"]).strip().lower() == "yes"]
        layout = "/".join(x for x in ("PE", "SE") if str(b[x]).strip().lower() == "yes")
        library = "/".join(x for x, col in (("whole transcriptome", "transciptome"), ("targeted", "targeted"))
                           if str(b[col]).strip().lower() == "yes")
        counts = g["material"].value_counts()
        rows.append({
            "benchmark": labels[int(b.benchmark_number)],
            "purpose": b.benchmark_purpose,
            "library": library or "not reported",
            "read_length_bp": ", ".join(lengths) or "not reported",
            "layout": layout or "not reported",
            "truth_set_type": truth_type(b.ground_truth_type),
            "parameters_reported": status(b.parameters_disclosed),
            "n_tools": int(g["tool_name_clean"].nunique()),
            "obs_simulated": int(counts.get("simulated", 0)),
            "obs_cell_line": int(counts.get("cell line", 0)),
            "obs_patient": int(counts.get("patient", 0)),
            "obs_mixed": int(counts.get("mixed", 0)),
            "obs_spike_in": int(counts.get("spike-in", 0)),
            "median_f1_simulated": g.loc[g["material"] == "simulated", "f1"].median(),
            "median_f1_real": g.loc[g["material"].isin(["cell line", "patient", "mixed", "real"]), "f1"].median(),
        })
    return pd.DataFrame(rows).round(3)


# ---------------------------------------------------------------------------
# Rows of Supplementary Tables S3-S7 for the LaTeX supplement (tables/S*_rows.tex).
# Only selection, labels and formatting; all values come from the tables above.
# Headers, column widths and captions are defined in supplement.tex.
# ---------------------------------------------------------------------------
LATEX_DIR = OUTPUT_DIR / "latex"
DASH = "--"
PURPOSE_LABELS = {
    "New-tool benchmark": "new tool",
    "Independent comparison": "independent",
    "Simulator development with a two-caller comparison": "simulator (two-caller comparison)",
}
LIBRARY_LABELS = {"whole transcriptome": "WT", "targeted": "targeted", "whole transcriptome/targeted": "WT, targeted"}
TEX_ESCAPE = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_",
              "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}", "|": r"\textbar{}",
              "–": "--", "—": "---"}


def tex(value) -> str:
    return "".join(TEX_ESCAPE.get(ch, ch) for ch in str(value))


def fmt3(value) -> str:
    return DASH if pd.isna(value) else f"{value:.3f}"


def fmt_int(value) -> str:
    return DASH if pd.isna(value) else str(int(value))


def latex_rows(rows: list[list], raw_columns: tuple[int, ...] = ()) -> str:
    """One LaTeX table row per line; cells are escaped unless their index is in raw_columns."""
    lines = []
    for row in rows:
        cells = [str(c) if i in raw_columns else tex(c) for i, c in enumerate(row)]
        lines.append(" & ".join(cells) + r" \\")
    return "\n".join(lines) + "\n"


def supplementary_table_rows(sra: pd.DataFrame, tools: pd.DataFrame, profiles: pd.DataFrame,
                             lists: pd.DataFrame, contexts: pd.DataFrame) -> dict[str, str]:
    s3 = latex_rows([[r.sra_run, r.cell_line, f"{int(r.read_count):,}",
                      " / ".join(str(r.events_27_38_53_99).split("/"))] for r in sra.itertuples()])
    s4 = latex_rows([[r.tool, r.family, r.benchmarks, r.n_observations, r.n_with_f1] for r in tools.itertuples()])
    s5 = latex_rows([[r.tool, r.n_observations, r.n_benchmarks, fmt3(r.mean_precision), fmt3(r.mean_recall), r.profile,
                      fmt_int(r.n_real_without_edgren), fmt3(r.mean_precision_real_without_edgren),
                      fmt3(r.mean_recall_real_without_edgren)] for r in profiles.itertuples()],
                    raw_columns=(3, 4, 6, 7, 8))
    s6 = latex_rows([[r["cell_line"], r["fusion"]] + [r"\checkmark" if r[f"list_{n}"] else "" for n in LISTS]
                     for _, r in lists.iterrows()], raw_columns=(2, 3, 4, 5))
    s7 = latex_rows([[r.benchmark, PURPOSE_LABELS.get(r.purpose, r.purpose), LIBRARY_LABELS.get(r.library, r.library),
                      r.read_length_bp, r.layout.replace("/", ", "), r.truth_set_type, r.n_tools,
                      " / ".join(str(getattr(r, c)) for c in ("obs_simulated", "obs_cell_line", "obs_patient",
                                                              "obs_mixed", "obs_spike_in")),
                      fmt3(r.median_f1_simulated), fmt3(r.median_f1_real)] for r in contexts.itertuples()],
                    raw_columns=(8, 9))
    return {"S3": s3, "S4": s4, "S5": s5, "S6": s6, "S7": s7}


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    book = pd.read_excel(INPUT, sheet_name=None)
    obs, data, bench, names = book["observations"], book["dataset_metadata"], book["benchmark_metadata"], book["tool_name_map"]

    sra = pd.read_csv(SRA_RUNS)
    sra.to_csv(OUTPUT_DIR / "S3_sra_runs.csv", index=False)
    tools = tools_table(obs, names, bench)
    tools.to_csv(OUTPUT_DIR / "S4_tools.csv", index=False)
    profiles = profiles_table()
    profiles.to_csv(OUTPUT_DIR / "S5_tool_profiles.csv", index=False)
    lists = edgren_lists_table()
    lists.to_csv(OUTPUT_DIR / "S6_edgren_reference_lists.csv", index=False)
    contexts = contexts_table(obs, data, bench)
    contexts.to_csv(OUTPUT_DIR / "S7_benchmark_contexts.csv", index=False)

    lines = ["# Supplementary tables: counts", ""]
    lines.append(f"- S4: {len(tools)} tools or variants from {tools['family'].nunique()} families; "
                 f"{int(tools['variant_or_output'].sum())} variants or outputs")
    lines.append(f"- S5: {len(profiles)} tools with at least 3 observations with precision and recall; "
                 f"{int(profiles['n_real_without_edgren'].notna().sum())} also in real data without Edgren")
    in99 = lists["list_99"]
    for n in (27, 38, 53):
        sub = lists[lists[f"list_{n}"]]
        lines.append(f"- S6: list {n}: {len(sub)} events, {int(sub['list_99'].sum())} also in list 99")
    lines.append(f"- S6: {len(lists)} distinct events in total; {int(in99.sum())} in list 99")
    params = contexts["parameters_reported"].value_counts()
    lines.append(f"- S7: parameters reported: yes {params.get('yes', 0)}, partial {params.get('partial', 0)}, "
                 f"no {params.get('no', 0)} (coding as in Figure 4)")
    (OUTPUT_DIR / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))

    LATEX_DIR.mkdir(parents=True, exist_ok=True)
    for name, rows in supplementary_table_rows(sra, tools, profiles, lists, contexts).items():
        (LATEX_DIR / f"{name}_rows.tex").write_text(rows, encoding="utf-8")
    print(f"LaTeX table rows: {LATEX_DIR}")


if __name__ == "__main__":
    main()
