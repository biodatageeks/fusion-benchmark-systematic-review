"""Supplementary Tables S3-S7: regenerated from the input, compared with the committed files, and sanity-checked."""
from __future__ import annotations

import pandas as pd
import pytest

import supplementary_tables as st
from conftest import RESULTS

OUT = RESULTS / "supplementary_tables"


@pytest.fixture(scope="module")
def tables(book):
    obs, data, bench, names = book["observations"], book["dataset_metadata"], book["benchmark_metadata"], book["tool_name_map"]
    sra = pd.read_csv(st.SRA_RUNS)
    built = {"S3": sra, "S4": st.tools_table(obs, names, bench), "S5": st.profiles_table(),
             "S6": st.edgren_lists_table(), "S7": st.contexts_table(obs, data, bench)}
    built["latex"] = st.supplementary_table_rows(sra, built["S4"], built["S5"], built["S6"], built["S7"])
    return built


FILES = {"S3": "S3_sra_runs.csv", "S4": "S4_tools.csv", "S5": "S5_tool_profiles.csv",
         "S6": "S6_edgren_reference_lists.csv", "S7": "S7_benchmark_contexts.csv"}


@pytest.mark.parametrize("name", list(FILES))
def test_committed_tables_are_current(tables, name):
    assert tables[name].to_csv(index=False) == (OUT / FILES[name]).read_text(encoding="utf-8")
    assert tables["latex"][name] == (OUT / "latex" / f"{name}_rows.tex").read_text(encoding="utf-8")


@pytest.mark.parametrize("name", list(FILES))
def test_latex_has_one_row_per_table_row(tables, name):
    lines = tables["latex"][name].rstrip("\n").split("\n")
    assert len(lines) == len(tables[name])
    assert all(line.endswith(r" \\") for line in lines)


def test_tex_escaping():
    assert st.tex("A_b & 50% #1 {x} a|b – c") == r"A\_b \& 50\% \#1 \{x\} a\textbar{}b -- c"


def test_s3_sra_runs(tables):
    s3 = tables["S3"]
    assert len(s3) == 6 and s3["sra_run"].is_unique
    assert s3.groupby("cell_line").size().to_dict() == {"BT-474": 2, "KPL-4": 1, "MCF-7": 1, "SK-BR-3": 2}
    assert (s3["read_length_bp"] == 50).all() and (s3["library_layout"] == "PAIRED").all()


def test_s4_tools(tables, obs):
    s4 = tables["S4"]
    assert len(s4) == 53 and s4["tool"].is_unique
    assert s4["family"].nunique() == 40 and int(s4["variant_or_output"].sum()) == 16
    assert s4["n_observations"].sum() == len(obs) == 615
    assert s4["n_with_f1"].sum() == obs["f1_score"].notna().sum() == 546
    assert (s4["benchmarks"].str.count("; ") + 1 == s4["n_benchmarks"]).all()
    per_tool = obs.groupby("tool_name_clean")["benchmark_number"].nunique()
    assert (s4.set_index("tool")["n_benchmarks"] == per_tool.reindex(s4["tool"]).values).all()


def test_s5_profiles(tables, obs):
    s5 = tables["S5"]
    assert len(s5) == 32 and int(s5["n_real_without_edgren"].notna().sum()) == 13
    assert (s5["n_observations"] >= 3).all()
    assert (s5["n_real_without_edgren"].dropna() >= 3).all()
    for c in ["mean_precision", "mean_recall", "mean_precision_real_without_edgren", "mean_recall_real_without_edgren"]:
        assert s5[c].dropna().between(0, 1).all()
    assert set(s5["profile"]) <= {"high recall / low precision", "low recall / high precision",
                                  "high recall / high precision", "low recall / low precision", "intermediate"}
    # recompute one row from the raw observations
    fc = obs[(obs.tool_name_clean == "FusionCatcher")].dropna(subset=["precision", "recall"])
    row = s5.set_index("tool").loc["FusionCatcher"]
    assert row.n_observations == len(fc) and row.n_benchmarks == fc.benchmark_number.nunique()
    assert row.mean_precision == pytest.approx(fc.precision.mean(), abs=5e-4)
    real = fc[(fc.dataset_origin == "Real") & ~fc.is_edgren.astype(bool)]
    assert row.n_real_without_edgren == len(real)
    assert row.mean_recall_real_without_edgren == pytest.approx(real.recall.mean(), abs=5e-4)


def test_s5_tools_with_recall_above_precision_in_real_data(tables):
    """Results 3.4: tools with mean recall at least 0.15 above mean precision (real, without Edgren)."""
    s5 = tables["S5"].dropna(subset=["n_real_without_edgren"])
    gap = s5.set_index("tool").eval("mean_recall_real_without_edgren - mean_precision_real_without_edgren")
    assert set(gap[gap >= 0.15].index) == {"EricScript", "JAFFA", "FuSeq", "SOAPfuse", "FusionMap", "ChimeraScan",
                                           "deFuse", "ChimPipe"}


def test_s6_reference_lists(tables):
    s6 = tables["S6"]
    flags = s6[[f"list_{n}" for n in st.LISTS]]
    assert flags.any(axis=1).all()
    assert flags.sum().to_dict() == {"list_27": 27, "list_38": 38, "list_53": 53, "list_99": 99}
    assert s6["records_merged"].ge(flags.sum(axis=1)).all()
    # every event of the 27-fusion list is in the 99-fusion list (Supp. Table S6 caption)
    assert s6.loc[s6.list_27, "list_99"].all()


def test_s6_no_within_list_merging():
    recs = st.edgren_records()
    for size in st.LISTS:
        same = [r for r in recs if r["list"] == size]
        assert len(same) == size
        for i, x in enumerate(same):
            for y in same[i + 1:]:
                assert not st.same_event(x, y), (size, x["fusion"], y["fusion"])


def test_s7_contexts(tables, obs, benchmarks):
    s7 = tables["S7"]
    assert len(s7) == 18 and s7["benchmark"].is_unique
    years = s7["benchmark"].str.extract(r"(\d{4})$")[0].astype(int)
    assert years.is_monotonic_increasing  # S7 caption: ordered by publication year
    cols = ["obs_simulated", "obs_cell_line", "obs_patient", "obs_mixed", "obs_spike_in"]
    assert s7[cols].sum().sum() == len(obs)
    assert s7["obs_spike_in"].sum() == (obs.dataset_origin.str.lower() == "spike-in").sum()
    assert s7["obs_simulated"].sum() == (obs.dataset_origin == "Simulated").sum()
    assert s7["parameters_reported"].value_counts().to_dict() == {"yes": 10, "partial": 6, "no": 2}
    labels = st.benchmark_labels(benchmarks)
    per_bench = obs.groupby("benchmark_number").size().rename(index=labels)
    assert (s7.set_index("benchmark")[cols].sum(axis=1) == per_bench.reindex(s7["benchmark"]).values).all()
    singh = obs[(obs.benchmark_number == 7) & (obs.dataset_origin == "Real")]
    assert s7.set_index("benchmark").loc["Singh 2021", "median_f1_real"] == pytest.approx(singh.f1_score.median(), abs=5e-4)
