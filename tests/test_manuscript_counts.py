"""Counts stated in the article, recomputed from the input workbook (not from results/).

If a data change moves one of these numbers, the test fails.
"""
from __future__ import annotations

import numpy as np
import pytest


def test_corpus_size(obs, datasets, benchmarks):
    """Abstract and Results 3.1: 18 benchmarks, 53 tools, 80 datasets from 30 biological datasets."""
    assert len(obs) == 615
    assert obs["benchmark_number"].nunique() == 18 == len(benchmarks)
    assert obs["tool_name_clean"].nunique() == 53
    assert obs["tool_family"].nunique() == 40
    used = datasets[datasets["dataset_id"].isin(obs["dataset_id"])]
    assert obs["dataset_id"].nunique() == 80 == len(used)
    assert used["biological_dataset_group"].nunique() == 30
    assert (benchmarks["year_x"].min(), benchmarks["year_x"].max()) == (2012, 2023)


def test_tools_per_benchmark(benchmarks):
    """Results 3.1: 2 to 27 tools per study, mean 8.8, median 7."""
    n = benchmarks["number of tools"]
    assert (n.min(), n.max(), round(n.mean(), 1), n.median()) == (2, 27, 8.8, 7)


def test_most_frequent_tools(obs):
    top = obs.groupby("tool_name_clean")["benchmark_number"].nunique().sort_values(ascending=False)
    assert top.index[0] == "TopHat-Fusion" and set(top.index[1:3]) == {"deFuse", "FusionCatcher"}


def test_real_versus_simulated_observations(obs, normalized):
    """Abstract and Results 3.3: 546 F1; 518 real-vs-simulated (302 simulated, 216 real); 28 spike-in."""
    assert obs["f1_score"].notna().sum() == 546
    f1 = normalized.dropna(subset=["F1"])
    assert len(f1) == 518
    assert f1["DatasetType"].value_counts().to_dict() == {"simulated": 302, "real": 216}
    spike = obs[obs["dataset_origin"].str.lower() == "spike-in"]
    assert spike["f1_score"].notna().sum() == 28
    means = f1.groupby("DatasetType")["F1"].agg(["mean", "median"]).round(3)
    assert means.loc["simulated"].tolist() == [0.685, 0.760]
    assert means.loc["real"].tolist() == [0.374, 0.352]
    for metric, (real, sim) in {"Precision": (0.421, 0.793), "Recall": (0.536, 0.651)}.items():
        m = normalized.dropna(subset=[metric]).groupby("DatasetType")[metric].mean().round(3)
        assert (m["real"], m["simulated"]) == (real, sim)


def test_patient_derived_observations(normalized, datasets):
    """Results 3.3: patient-derived data, 25 F1 observations from three datasets."""
    d = normalized.dropna(subset=["F1"]).merge(datasets[["dataset_id", "real_subtype"]], left_on="Dataset",
                                                right_on="dataset_id")
    patient = d[(d["DatasetType"] == "real") & (d["real_subtype"] == "patient-derived")]
    assert (len(patient), patient["Dataset"].nunique()) == (25, 3)


def test_benchmarks_with_both_origins(normalized):
    """Results 3.3: seven benchmarks reported F1 for both real and simulated data."""
    f1 = normalized.dropna(subset=["F1"])
    assert int(f1.groupby("Benchmark")["DatasetType"].nunique().eq(2).sum()) == 7


def test_read_length_groups(normalized):
    """Results 3.4: 392 observations with 50/75/100 bp; all 107 at 75 bp are simulated."""
    from read_length_analysis import canonical_read_length

    rl = normalized["Read_length"].map(canonical_read_length)
    d = normalized.assign(rl=rl).dropna(subset=["F1", "rl"])
    assert len(d) == 392
    assert d.groupby("rl").size().to_dict() == {50: 160, 75: 107, 100: 125}
    assert set(d.loc[d.rl == 75, "DatasetType"]) == {"simulated"}


@pytest.mark.parametrize("value, expected", [(50, 50), (51, 50), (76, 75), (101, 100), (100, 100), (36, None),
                                             (125, None), (250, None), (None, None)])
def test_canonical_read_length(value, expected):
    from read_length_analysis import canonical_read_length

    out = canonical_read_length(np.nan if value is None else value)
    assert (np.isnan(out) if expected is None else out == expected)


def test_edgren_observations(obs, book):
    """Edgren-derived datasets (Supp. Table S1; Results 3.5)."""
    edgren = obs[obs["is_edgren"].astype(bool)]
    assert edgren["dataset_id"].nunique() == 11 and len(edgren) == 96
    # Ten benchmarks used Edgren data: eight with Edgren observations, Haas 2019 (53-fusion list, no observations) and
    # Uhrig 2021 (Edgren runs scored, but the extracted precision is for the ENCODE MCF-7 sample, Dataset2.3)
    assert edgren["benchmark_number"].nunique() == 8
    assert not obs.loc[obs["dataset_id"] == "Dataset2.3", "is_edgren"].astype(bool).any()
    readme = dict(zip(book["README"]["field"], book["README"]["value"].astype(str)))
    assert set(readme["edgren_dataset_ids"].replace(" ", "").split(",")) == \
        {d.replace("Dataset", "") for d in edgren["dataset_id"].unique()}


def test_benchmark_purpose_counts(benchmarks, normalized):
    """Results 3.6: 13 new-tool and 5 independent benchmarks; four independent ones have real/simulated F1."""
    purpose = benchmarks["benchmark_purpose"].str.lower()
    assert purpose.str.startswith("new-tool").sum() == 13
    assert (~purpose.str.startswith("new-tool")).sum() == 5
    f1 = normalized.dropna(subset=["F1"])
    assert f1.groupby("Benchmark_purpose")["Benchmark"].nunique().to_dict() == \
        {"independent_comparison": 4, "new_tool": 12}


def test_figure3_and_figure4_codings(results):
    """Results 3.1 and 3.2, counted from the committed heatmap source tables (Figures 3 and 4)."""
    stats = results("benchmark_statistics/benchmark_statistics_heatmap_display.csv").set_index("Benchmark")
    count = {c: stats[c].value_counts().to_dict() for c in stats.columns[2:]}
    assert count["Simulated"].get("full") == 13
    assert count["Transcriptome"].get("full") == 17 and count["Targeted"].get("full") == 2
    assert count["SE"].get("full") == 2 and count["PE"].get("full") == 18
    assert [count[f"{n} bp"].get("full") for n in (50, 75, 100, 250)] == [11, 10, 11, 1]
    no_length = (stats[[f"{n} bp" for n in (50, 75, 100, 250)]] == "unknown").all(axis=1).sum()
    assert no_length == 3
    real = (stats["Real cells"] == "full") | (stats["Real patients"] == "full")
    assert list(stats.index[~real]) == ["Tan 2015"]

    rep = results("benchmark_reporting/benchmark_reporting_heatmap_display.csv").set_index("Unnamed: 0")
    assert rep["Truth-set type"].value_counts().to_dict() == {
        "Published validated fusions": 7, "Proxy-based": 6, "RNA + WGS + experimental validation": 3,
        "Clinically validated/expected fusions": 1, "No defined truth set / wisdom of crowds": 1}
    # Results 3.2, Key Points: three benchmarks with materials to rerun the benchmark
    assert rep["Fully reproducible"].eq("yes").sum() == 3
    assert set(rep.index[rep["Fully reproducible"].eq("yes")]) == {"Haas 2019", "Rodriguez-Martin 2017", "Creason 2021"}
    assert rep["Tool versions"].eq("yes").sum() == 11
    assert rep["Parameters"].value_counts().to_dict() == {"yes": 10, "partial": 6, "no": 2}


def test_heatmap_rows_are_chronological(results):
    """Figure 3 and 4 captions: rows ordered by publication year, most recent first."""
    for path, index in [("benchmark_statistics/benchmark_statistics_heatmap_display.csv", "Benchmark"),
                        ("benchmark_reporting/benchmark_reporting_heatmap_display.csv", "Unnamed: 0")]:
        years = results(path)[index].str.extract(r"(\d{4})$")[0].astype(int)
        assert years.is_monotonic_decreasing, path
