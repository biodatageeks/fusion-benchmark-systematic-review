"""Consistency of the analysis input built from data/raw/benchmark.xlsx and regression of the main model counts."""
import numpy as np
import pandas as pd
import pytest

from conftest import REPO

INPUT = REPO / "data" / "processed" / "analysis_input.xlsx"


@pytest.fixture(scope="module")
def obs():
    return pd.read_excel(INPUT, sheet_name="observations")


def test_shape_and_benchmarks(obs):
    assert len(obs) == 1182
    assert sorted(obs.benchmark_number.unique()) == list(range(1, 14))
    assert obs.observation_id.is_unique


def test_metrics_in_unit_interval(obs):
    for c in ["recall", "precision", "f1_score"]:
        v = obs[c].dropna()
        assert ((v >= 0) & (v <= 1)).all(), c


def test_f1_consistent_with_precision_and_recall(obs):
    # Benchmark 13 (SMC-RNA DREAM) reports metrics averaged over samples; the mean F1 is not the harmonic
    # mean of the mean precision and recall, so it is checked separately below.
    d = obs.dropna(subset=["recall", "precision", "f1_score"])
    d = d[d.benchmark_number != 13]
    d = d[(d.recall + d.precision) > 0]
    hm = 2 * d.precision * d.recall / (d.precision + d.recall)
    bad = d[(hm - d.f1_score).abs() > 0.02]
    assert bad.empty, bad[["observation_id", "precision", "recall", "f1_score"]].head(20).to_string()


def test_dream_averaged_f1_close_to_harmonic_mean(obs):
    d = obs[(obs.benchmark_number == 13)].dropna(subset=["recall", "precision", "f1_score"])
    d = d[(d.precision + d.recall) > 0]
    hm = 2 * d.precision * d.recall / (d.precision + d.recall)
    assert len(d) >= 15 and ((hm - d.f1_score) >= -0.005).all() and ((hm - d.f1_score) <= 0.03).all()


def test_counts_consistent_with_rates(obs):
    d = obs.dropna(subset=["TP", "truth_total", "recall"])
    d = d[d.truth_total > 0]
    assert ((d.TP / d.truth_total - d.recall).abs() <= 0.02).all()
    d = obs.dropna(subset=["TP", "FP", "precision"])
    d = d[(d.TP + d.FP) > 0]
    assert ((d.TP / (d.TP + d.FP) - d.precision).abs() <= 0.02).all()


def test_model_rows_have_annotations(obs):
    m = obs[obs.include_in_models.astype(bool)]
    assert m.dataset_origin.notna().all()
    ann = pd.read_csv(REPO / "data" / "processed" / "dataset_annotations.csv", dtype=str)
    assert set(m.dataset_id.astype(str)) <= set(ann.dataset_id)
    assert ann.loc[ann.dataset_id.isin(m.dataset_id.astype(str)), "bio_dataset_id"].notna().all()
    excluded = obs[~obs.include_in_models.astype(bool)]
    assert excluded.exclusion_reason.notna().all()


def test_crossed_model_units_and_point_estimate():
    from revision_utils import load_observations, unit_counts
    from crossed_dependence_models import prepare, crossed_estimates
    d = prepare(load_observations(), "F1")
    assert unit_counts(d) == {"n_obs": 549, "n_tools": 55, "n_bio_datasets": 20, "n_benchmarks": 11}
    est = crossed_estimates(d)
    assert est["m1_diff"] == pytest.approx(-0.4596, abs=1e-3)
    assert est["m2_diff"] == pytest.approx(-0.3516, abs=1e-3)
    stored = pd.read_csv(REPO / "results" / "revision" / "crossed_models" / "models.csv")
    row = stored[(stored.metric == "F1") & stored.model.str.startswith("M1") & (stored.scale == "metric difference")].iloc[0]
    assert row.estimate == pytest.approx(est["m1_diff"], abs=1e-9)


def test_squeeze_maps_into_open_interval():
    from revision_utils import squeeze
    s = squeeze(pd.Series([0.0, 0.5, 1.0]))
    assert (s > 0).all() and (s < 1).all() and s.iloc[1] == pytest.approx(0.5)


def test_latex_escape():
    from revision_utils import latex_escape
    assert latex_escape("a_b & 50% #1 $x^2$ {y} ~") == r"a\_b \& 50\% \#1 \$x\textasciicircum{}2\$ \{y\} \textasciitilde{}"
