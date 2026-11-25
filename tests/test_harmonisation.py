"""Metric harmonisation (normalize_columns) and observation-level sanity of the input workbook."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from real_simulated_sensitivity import normalize_columns

METRICS = ["precision", "recall", "f1_score"]
# DREAM SMC-RNA (Creason 2021, benchmark 13): the leaderboard reports per-sample-averaged F1, precision and
# recall, so F1 is not the harmonic mean of the reported precision and recall. Documented exception.
F1_NOT_HARMONIC_BENCHMARKS = {13}


def half_unit(value: float) -> float:
    """Half a unit in the last reported decimal (values rounded in the source, e.g. 0.01 in Singh 2021)."""
    text = repr(float(value))
    if "e" in text or "." not in text:
        return 1e-9
    decimals = len(text.split(".")[1].rstrip("0"))
    return 0.5 * 10 ** (-decimals) if decimals <= 4 else 1e-6


def harmonic(p: float, r: float) -> float:
    return 0.0 if p + r == 0 else 2 * p * r / (p + r)


# --------------------------------------------------------------------------- normalize_columns

def toy(**columns) -> pd.DataFrame:
    base = {"benchmark_number": [1, 1, 2], "dataset_id": ["D1", "D2", "D3"], "tool_name_clean": ["A", "B", "A"],
            "dataset_type_reported": ["Simulated", "real (cells)", "spike-in"]}
    base.update(columns)
    return pd.DataFrame(base)


def test_aliases_and_dataset_type_from_origin():
    df = normalize_columns(toy(dataset_origin=["Simulated", "Real", "spike-in"], f1_score=[0.5, 0.4, 0.3],
                               is_edgren=[False, True, False]))
    assert {"Benchmark", "Dataset", "DatasetType", "Tool", "F1"} <= set(df.columns)
    # spike-in is neither real nor simulated and is dropped from the real-vs-simulated analyses
    assert df["DatasetType"].tolist() == ["simulated", "real"]
    assert df["Is_Edgren_Binary"].tolist() == ["no", "yes"]
    assert df["Benchmark"].map(type).eq(str).all()


def test_origin_overrides_reported_type():
    df = normalize_columns(toy(dataset_type_reported=["simulated", "simulated", "simulated"],
                               dataset_origin=["Real", "Simulated", "Real"]))
    assert df["DatasetType"].tolist() == ["real", "simulated", "real"]


def test_missing_metric_columns_are_derived_from_counts():
    df = normalize_columns(toy(dataset_origin=["Simulated", "Real", "Real"], TP=[8, 3, 0], FP=[2, 9, 0],
                               truth_total=[10, 6, 4]))
    assert df["FN"].tolist() == [2, 3, 4]
    assert np.allclose(df["Recall"], [0.8, 0.5, 0.0])
    assert np.allclose(df["Precision"].iloc[:2], [0.8, 0.25]) and np.isnan(df["Precision"].iloc[2])
    assert np.isclose(df["F1"].iloc[0], 0.8) and np.isclose(df["F1"].iloc[1], 2 * 0.25 * 0.5 / 0.75)


def test_reported_metrics_are_kept_not_recomputed():
    """Methods 2.3: published values are retained when reported."""
    df = normalize_columns(toy(dataset_origin=["Simulated", "Real", "Real"], TP=[8, 3, 1], FP=[2, 9, 1],
                               truth_total=[10, 6, 4], precision=[0.81, 0.24, 0.5], recall=[0.79, 0.51, 0.25],
                               f1_score=[0.80, 0.33, np.nan]))
    assert df["Precision"].tolist() == [0.81, 0.24, 0.5]
    assert df["F1"].iloc[:2].tolist() == [0.80, 0.33] and np.isnan(df["F1"].iloc[2])


@pytest.mark.parametrize("raw, expected", [
    ("New-tool benchmark", "new_tool"),
    ("Independent comparison", "independent_comparison"),
    ("new_tool", "new_tool"),
    ("Simulator development with a two-caller comparison", "independent_comparison"),  # Tan 2015
])
def test_benchmark_purpose_normalisation(raw, expected):
    df = normalize_columns(toy(dataset_origin=["Simulated"] * 3, benchmark_purpose=[raw] * 3))
    assert set(df["Benchmark_purpose"]) == {expected}


def test_every_observation_gets_a_two_level_purpose(normalized):
    assert set(normalized["Benchmark_purpose"]) == {"new_tool", "independent_comparison"}


# --------------------------------------------------------------------------- workbook sanity

def test_metrics_in_unit_interval(obs):
    for column in METRICS + ["precision_recalculated", "recall_recalculated", "f1_recalculated"]:
        values = obs[column].dropna()
        assert values.between(0, 1).all(), column


def test_counts_non_negative_and_tp_within_truth(obs):
    for column in ["TP", "FP", "TP_plus_FP", "truth_total", "FN_reported"]:
        assert (obs[column].dropna() >= 0).all(), column
    x = obs.dropna(subset=["TP", "truth_total"])
    assert (x["TP"] <= x["truth_total"] + 1e-9).all()
    x = obs.dropna(subset=["TP", "FP", "TP_plus_FP"])
    assert np.allclose(x["TP"] + x["FP"], x["TP_plus_FP"], atol=0.5)
    x = obs.dropna(subset=["TP", "truth_total", "FN_reported"])
    assert np.allclose(x["truth_total"] - x["TP"], x["FN_reported"], atol=0.5)


def test_precision_and_recall_match_counts(obs):
    x = obs.dropna(subset=["TP", "FP", "precision"])
    x = x[(x.TP + x.FP) > 0]
    assert ((x.TP / (x.TP + x.FP) - x.precision).abs() <= 0.01).all()
    x = obs.dropna(subset=["TP", "truth_total", "recall"])
    assert ((x.TP / x.truth_total - x.recall).abs() <= 0.01).all()


def test_f1_is_harmonic_mean_of_precision_and_recall(obs):
    """Allowing for rounding of the reported values; DREAM is the only documented exception."""
    x = obs.dropna(subset=METRICS)
    bad = []
    for o in x.itertuples():
        dp, dr, df = half_unit(o.precision), half_unit(o.recall), half_unit(o.f1_score)
        lo = harmonic(max(o.precision - dp, 0), max(o.recall - dr, 0))
        hi = harmonic(min(o.precision + dp, 1), min(o.recall + dr, 1))
        if not lo - df - 1e-6 <= o.f1_score <= hi + df + 1e-6:
            bad.append((o.observation_id, o.benchmark_number))
    assert {b for _, b in bad} <= F1_NOT_HARMONIC_BENCHMARKS, bad


def test_observation_keys(obs, datasets, benchmarks):
    assert obs["observation_id"].is_unique
    assert set(obs["dataset_id"]) <= set(datasets["dataset_id"])
    assert set(obs["benchmark_number"]) == set(benchmarks["benchmark_number"])
    prefix = obs["dataset_id"].str.extract(r"^Dataset(\d+)")[0].astype(int)
    assert (prefix == obs["benchmark_number"]).all(), "dataset id must start with its benchmark number"
    assert (obs["dataset_id"].str.replace("Dataset", "") == obs["dataset_number"].astype(str)).all()


def test_observation_metadata_agrees_with_dataset_metadata(obs, datasets):
    m = obs.merge(datasets, on="dataset_id", suffixes=("", "_dm"), validate="many_to_one")
    assert (m["dataset_origin"].str.lower() == m["dataset_origin_dm"].str.lower()).all()
    assert (m["is_edgren"].astype(bool) == m["is_edgren_dm"].astype(bool)).all()
    assert set(datasets.loc[datasets.is_edgren.astype(bool), "dataset_id"]) == \
        set(datasets.loc[datasets.biological_dataset_group == "EDGREN", "dataset_id"])


def test_tool_names_follow_tool_name_map(obs, book):
    names = book["tool_name_map"]
    assert names["tool_reported"].is_unique
    m = obs.merge(names, on="tool_reported", how="left", suffixes=("", "_map"), validate="many_to_one")
    assert m["tool_name_clean_map"].notna().all()
    assert (m["tool_name_clean"] == m["tool_name_clean_map"]).all()
    assert (m["tool_family"] == m["tool_family_map"]).all()


def test_real_simulated_label_is_the_same_in_every_script(obs, normalized):
    """normalize_columns, precision_recall_tradeoff and the S7 context table must classify observations alike."""
    from precision_recall_tradeoff import load_clean_data

    pr = load_clean_data()
    origin = obs["dataset_origin"].str.lower().str.strip()
    assert (pr["Is_real"] == origin.eq("real")).all() and (pr["Is_simulated"] == origin.eq("simulated")).all()
    assert set(normalized.index) == set(origin.index[origin.isin(["real", "simulated"])])
    assert (normalized["DatasetType"] == origin.loc[normalized.index]).all()
    assert (pr["Is_Edgren"] == obs["is_edgren"].astype(bool)).all()
    assert (normalized["Is_Edgren_Binary"].eq("yes") == obs.loc[normalized.index, "is_edgren"].astype(bool)).all()
