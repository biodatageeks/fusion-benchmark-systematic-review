"""Supplementary tool tables (R1 minor comments m3 and m7).

- tool_list.csv / supp_tool_list.tex        all harmonised tools and variants with reported names and benchmarks.
- tool_profiles.csv / supp_tool_profiles.tex
      mean/median precision, recall and F1 per tool and data origin with the number of
      observations, benchmarks and biological datasets behind each summary.
"""

from __future__ import annotations

import pandas as pd

from real_simulated_sensitivity import DEFAULT_INPUT, load_table
from revision_utils import CURATION, REVISION_DIR, benchmark_labels, fmt, load_observations, write_latex_table

OUT = REVISION_DIR / "tools"
TABLES = REVISION_DIR / "tables"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    df = load_observations()
    name_map = load_table(DEFAULT_INPUT, sheet_name="tool_name_map")
    raw = load_table(DEFAULT_INPUT, sheet_name="observations")

    # tool list: every harmonised name used in the observation table (incl. spike-in rows)
    raw["Benchmark_label"] = pd.to_numeric(raw["benchmark_number"], errors="coerce").astype(int).astype(str)
    labels = benchmark_labels()
    raw["Benchmark_label"] = raw["Benchmark_label"].map(labels).fillna(raw["Benchmark_label"])
    tl = raw.groupby("tool_name_clean").agg(
        reported_names=("tool_reported", lambda s: "; ".join(sorted(set(map(str, s))))),
        tool_family=("tool_family", "first"),
        versions_reported=("tool_version", lambda s: "; ".join(sorted(set(map(str, s.dropna())))) or "--"),
        n_benchmarks=("benchmark_number", "nunique"),
        benchmarks=("Benchmark_label", lambda s: "; ".join(sorted(set(s)))),
    ).reset_index().rename(columns={"tool_name_clean": "tool"})
    variant = name_map.drop_duplicates("tool_name_clean").set_index("tool_name_clean")["is_variant_or_output"]
    tl["variant_or_filtered_output"] = tl["tool"].map(variant).map({True: "yes", False: "no"}).fillna("no")
    # tool age and maintenance (R1 minor comment m6); sources in data/curation/tool_release_maintenance.csv
    maint = pd.read_csv(CURATION / "tool_release_maintenance.csv")[
        ["tool_family", "first_publication_year", "latest_release_date", "last_commit_date", "maintenance_status"]]
    tl = tl.merge(maint, on="tool_family", how="left")
    tl = tl.sort_values(["tool_family", "tool"])
    tl.to_csv(OUT / "tool_list.csv", index=False)
    fam = tl.drop_duplicates("tool_family")
    print("maintenance status by family:", fam["maintenance_status"].value_counts().to_dict())
    show = tl.assign(
        first_publication_year=tl["first_publication_year"].map(lambda v: "--" if pd.isna(v) else str(int(v))),
        maintenance_status=tl["maintenance_status"].fillna("--"),
    )
    write_latex_table(
        show[["tool", "tool_family", "variant_or_filtered_output", "first_publication_year", "maintenance_status",
              "reported_names", "n_benchmarks", "benchmarks"]]
        .rename(columns={"tool": "Tool", "tool_family": "Family", "variant_or_filtered_output": "Variant",
                         "first_publication_year": "Published", "maintenance_status": "Status",
                         "reported_names": "Reported as", "n_benchmarks": "Benchmarks (n)", "benchmarks": "Benchmarks"}),
        TABLES / "supp_tool_list.tex",
        caption=rf"\textbf{{Harmonised fusion detection tools and variants ({tl['tool'].nunique()} entries, "
                rf"{tl['tool_family'].nunique()} tool families).}} Variant marks filtered outputs or alternative modes "
                r"(e.g. Arriba-hc versus Arriba), which are analysed as separate entries. Published: year of the first "
                r"publication of the tool family. Status, assessed on 17 November 2026: active, release or commit since "
                r"November 2023; inactive, none since then. Release and commit dates with sources are in "
                r"\texttt{data/curation/tool\_release\_maintenance.csv}.",
        label="tab:supp_tool_list", colspec="lllllp{0.17\\linewidth}rp{0.24\\linewidth}", font=r"\scriptsize", longtable=True)

    # tool profiles
    rows = []
    for (tool, real), g in df.groupby(["Tool", "real"]):
        pr = g.dropna(subset=["Precision", "Recall"])
        rows.append({
            "tool": tool, "data": "real" if real else "simulated",
            "n_obs": len(g), "n_obs_pr": len(pr), "n_benchmarks": g["Benchmark"].nunique(),
            "n_bio_datasets": g["bio_dataset_id"].nunique(),
            "mean_precision": g["Precision"].mean(), "median_precision": g["Precision"].median(),
            "mean_recall": g["Recall"].mean(), "median_recall": g["Recall"].median(),
            "mean_f1": g["F1"].mean(), "median_f1": g["F1"].median(),
        })
    prof = pd.DataFrame(rows).sort_values(["data", "n_obs", "tool"], ascending=[True, False, True])
    prof.to_csv(OUT / "tool_profiles.csv", index=False)
    show = prof[prof["n_benchmarks"] >= 2].copy()
    show_tex = show.assign(**{
        "Tool": show.tool, "Data": show.data,
        "Obs/bench/datasets": show.apply(lambda r: f"{r.n_obs}/{r.n_benchmarks}/{r.n_bio_datasets}", axis=1),
        "Precision mean (median)": show.apply(lambda r: f"{fmt(r.mean_precision, 2)} ({fmt(r.median_precision, 2)})", axis=1),
        "Recall mean (median)": show.apply(lambda r: f"{fmt(r.mean_recall, 2)} ({fmt(r.median_recall, 2)})", axis=1),
        "F1 mean (median)": show.apply(lambda r: f"{fmt(r.mean_f1, 2)} ({fmt(r.median_f1, 2)})", axis=1),
    })[["Tool", "Data", "Obs/bench/datasets", "Precision mean (median)", "Recall mean (median)", "F1 mean (median)"]]
    show_tex = show_tex.replace({"NA (NA)": "--"}, regex=False)
    write_latex_table(
        show_tex, TABLES / "supp_tool_profiles.tex",
        caption=r"\textbf{Tool-level precision, recall and F1 by data origin.} Tools evaluated in at least two benchmarks "
                r"for the given data origin. Precision on real data is observed precision against each benchmark's "
                r"reference set and is affected by reference incompleteness. Summaries pool different datasets, truth sets "
                r"and scoring rules and are descriptive; the full table is in \texttt{results/revision/tools/tool\_profiles.csv}.",
        label="tab:supp_tool_profiles", colspec="llrlll", font=r"\scriptsize", longtable=True)
    print(f"tools: {tl['tool'].nunique()}, families: {tl['tool_family'].nunique()}, variants: {(tl.variant_or_filtered_output == 'yes').sum()}")
    print(show.head(20).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
