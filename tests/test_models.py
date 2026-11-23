"""Model inputs and committed model outputs (no bootstrap refits; each test runs in seconds)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import crossed_models as cm
import exploratory_random_effects as ere


@pytest.fixture(scope="module")
def crossed_data() -> pd.DataFrame:
    return cm.load()


# --------------------------------------------------------------------------- crossed models (M0-M3)

def test_crossed_input(crossed_data):
    d = crossed_data
    assert d["biological_dataset_group"].notna().all()
    assert set(d["DatasetType"]) == {"real", "simulated"}
    assert set(d.loc[d.real == 1, "Subtype"]) == {"cell line", "patient-derived", "mixed"}
    assert set(d.loc[d.real == 0, "Subtype"]) == {"simulated"}


@pytest.mark.parametrize("metric", cm.METRICS)
def test_smithson_verkuilen_transform(crossed_data, metric):
    d = cm.prepare(crossed_data, metric)
    n = len(d)
    assert d["ys"].between(0.5 / n, 1 - 0.5 / n).all()
    assert np.isfinite(d["yl"]).all()
    assert np.allclose(cm.expit(d["yl"]), d["ys"])
    # order-preserving
    assert (np.argsort(d["y"].to_numpy(), kind="stable") == np.argsort(d["ys"].to_numpy(), kind="stable")).all()


def test_m1_point_estimate_reproduces(crossed_data, results):
    """Primary model (M1) refitted once on F1; equals the committed estimate (Abstract, Results 3.3)."""
    d = cm.prepare(crossed_data, "F1")
    est = cm.estimates(d, beta=False)
    committed = results("crossed_models/real_vs_simulated_models.csv")
    row = committed[(committed.metric == "F1") & (committed.model == "M1") & (committed.scale == "metric difference")].iloc[0]
    assert est["m1_diff"] == pytest.approx(row.estimate, abs=1e-6)
    assert round(est["m1_diff"], 3) == -0.407


def test_bootstrap_files_match_reported_intervals(results):
    """CIs in real_vs_simulated_models.csv are the 2.5/97.5 percentiles of the committed bootstrap replicates."""
    models = results("crossed_models/real_vs_simulated_models.csv")
    for metric in cm.METRICS:
        boot = results(f"crossed_models/bootstrap_{metric.lower()}.csv")
        assert len(boot) == 999
        for col, model, scale in [("m1_diff", "M1", "metric difference"), ("m2_diff", "M2", "metric difference"),
                                  ("m3_diff", "M3", "metric difference"), ("m1_logit", "M1", "logit")]:
            r = models[(models.metric == metric) & (models.model == model) & (models.scale == scale)
                       & models.ci_method.str.startswith("benchmark")].iloc[0]
            lo, hi = np.nanpercentile(boot[col], [2.5, 97.5])
            assert (r.ci_low, r.ci_high) == pytest.approx((lo, hi)), (metric, col)
            assert r.ci_method == "benchmark cluster bootstrap (B=999)"


def test_manuscript_model_numbers(results):
    """Numbers quoted in the Abstract and Results 3.3 (estimate, bootstrap CI, Wald p)."""
    m = results("crossed_models/real_vs_simulated_models.csv")

    def get(metric, model, method):
        r = m[(m.metric == metric) & (m.model == model) & (m.scale == "metric difference") & m.ci_method.str.startswith(method)]
        return r.iloc[0]

    f1, f1w = get("F1", "M1", "benchmark"), get("F1", "M1", "Wald")
    assert (round(f1.estimate, 3), round(f1.ci_low, 3), round(f1.ci_high, 3), round(f1w.p, 3)) == (-0.407, -0.566, -0.179, 0.002)
    p, pw = get("Precision", "M1", "benchmark"), get("Precision", "M1", "Wald")
    assert (round(p.estimate, 3), round(p.ci_low, 3), round(p.ci_high, 3)) == (-0.520, -0.774, -0.243) and pw.p < 0.001
    r, rw = get("Recall", "M1", "benchmark"), get("Recall", "M1", "Wald")
    assert (round(r.estimate, 3), round(r.ci_low, 3), round(r.ci_high, 3), round(rw.p, 2)) == (-0.102, -0.298, 0.029, 0.42)
    b = get("F1", "M3", "benchmark")
    assert (round(b.estimate, 3), round(b.ci_low, 3), round(b.ci_high, 3)) == (-0.274, -0.370, -0.183)
    m0 = get("F1", "M0", "Wald")
    assert (round(m0.estimate, 3), round(m0.ci_low, 3), round(m0.ci_high, 3)) == (-0.338, -0.378, -0.297)


def test_manuscript_sensitivity_numbers(results):
    sens = results("crossed_models/sensitivity_lobo_edgren.csv")
    lobo = sens[(sens.metric == "F1") & sens.excluded.str.startswith("benchmark")]
    assert (round(lobo.m1_diff.min(), 3), round(lobo.m1_diff.max(), 3)) == (-0.514, -0.364)
    ed = sens[(sens.metric == "F1") & (sens.excluded == "Edgren-derived observations")].iloc[0]
    assert (round(ed.m1_diff, 3), round(ed.m1_logit_p, 3)) == (-0.405, 0.008)
    sub = results("crossed_models/real_subtype_models.csv")
    pc = sub[(sub.metric == "F1") & (sub.comparison == "patient-derived vs cell line")].iloc[0]
    assert (round(pc.m2_diff, 3), round(pc.m2_ci_low, 3), round(pc.m2_ci_high, 3)) == (-0.077, -0.320, 0.166)
    vc = results("crossed_models/variance_components.csv")
    vc = vc[vc.model == "M1 (logit)"].set_index("metric")
    share = vc[["bench", "bio", "tool"]].div(vc[["bench", "bio", "tool", "residual"]].sum(axis=1), axis=0).mul(100).round()
    assert (share.bio.min(), share.bio.max()) == (22, 34)
    assert (share.bench.min(), share.bench.max()) == (2, 19)
    assert (share.tool.min(), share.tool.max()) == (8, 16)


def test_within_benchmark_contrast():
    d = pd.DataFrame({"Benchmark": ["1"] * 4 + ["2"] * 2, "Tool": ["A", "A", "B", "B", "A", "A"],
                      "real": [0, 1, 0, 1, 0, 0], "y": [0.8, 0.5, 0.6, 0.6, 0.9, 0.7]})
    w = cm.within_benchmark(d, "F1")
    assert w["benchmark"].tolist() == ["1"]  # benchmark 2 has no real data
    r = w.iloc[0]
    assert r.diff_real_minus_sim == pytest.approx(0.55 - 0.7) and r.paired_tool_diff_mean == pytest.approx((-0.3 + 0) / 2)


def test_purpose_descriptive_numbers(results):
    """Results 3.6 (descriptive only)."""
    p = results("crossed_models/benchmark_purpose_descriptive.csv")
    f1 = p[(p.metric == "F1")].set_index(["subset", "purpose"])
    new, ind = f1.loc[("all tools", "new_tool")], f1.loc[("all tools", "independent_comparison")]
    assert (round(new["median"], 3), round(new.q1, 3), round(new.q3, 3)) == (0.722, 0.411, 0.850)
    assert (round(ind["median"], 3), round(ind.q1, 3), round(ind.q3, 3)) == (0.472, 0.209, 0.716)
    assert round(f1.loc[("without the authors' own new tool", "new_tool"), "median"], 3) == 0.658
    same = results("crossed_models/benchmark_purpose_same_tool.csv").set_index(["metric", "subset"])
    a, b = same.loc[("F1", "all tools")], same.loc[("F1", "without the authors' own new tool")]
    assert (a.n_higher_in_new_tool, a.n_tools_in_both, round(a.median_tool_diff_new_minus_independent, 2)) == (18, 21, 0.17)
    assert (b.n_higher_in_new_tool, b.n_tools_in_both) == (17, 20)


def test_new_tool_ranking(results):
    s = results("new_tool_performance/new_tool_benchmark_summary.csv")
    ok = s[s.status == "ok"]
    assert len(ok) == 12 and ok.new_tool_is_best_by_mean_f1.all()
    assert s.loc[s.status == "recall only", "new_tool_is_best_by_primary_metric"].tolist() == [True]


@pytest.mark.parametrize("text, expected", [
    ("yes, STAR-Fusion and TrinityFusion", ["STAR-Fusion", "TrinityFusion"]),
    ("Yes: SimFuse (simulator)", ["SimFuse"]),
    ("Fusion InPipe", ["Fusion InPipe"]),
    ("no", []), ("nan", []),
])
def test_parse_new_tool(text, expected):
    from new_tool_benchmark_performance import parse_new_tool

    assert parse_new_tool(text) == expected


# --------------------------------------------------------------------------- read length and Edgren drift

def test_read_length_numbers(results):
    k = results("read_length/read_length_kruskal.csv").iloc[0]
    assert round(k.H, 2) == 27.86 and k.p_value < 0.001
    # Results 3.4: crossed random intercepts (tool, dataset group, benchmark); neither 50 nor 75 bp differs from 100 bp
    m = results("read_length/read_length_mixed_model.csv").set_index("term")
    r50, r75 = m.loc["Read_length[T.50]"], m.loc["Read_length[T.75]"]
    assert (round(r50["diff"], 3), round(r50.ci_low, 3), round(r50.ci_high, 3), round(r50.p_value, 3)) == (-0.063, -0.181, 0.058, 0.313)
    assert r75.p_value > 0.05 and r75.ci_low < 0 < r75.ci_high
    assert (m["n"] == 392).all() and m.m1_converged.all() and m.m2_converged.all()


def test_edgren_drift_numbers(results):
    """Abstract and Results 3.5: 16 tools, |shift| 0.002-0.214, nine lower and seven higher with the 99 list."""
    d = results("tables/edgren_truth_set_drift.csv")
    assert len(d) == 16 and d.Tool.is_unique
    assert (round(d.delta_f1.abs().min(), 3), round(d.delta_f1.abs().max(), 3)) == (0.002, 0.214)
    assert ((d.delta_f1 < 0).sum(), (d.delta_f1 > 0).sum()) == (9, 7)


# --------------------------------------------------------------------------- DerSimonian-Laird / I2

def test_der_simonian_laird_known_values():
    y = np.array([0.2, 0.5, 0.8])
    v = np.array([0.01, 0.01, 0.01])
    r = ere.der_simonian_laird(y, v)
    q = np.sum((y - y.mean()) ** 2 / 0.01)  # 18
    assert r["Q"] == pytest.approx(q)
    assert r["I2"] == pytest.approx((q - 2) / q * 100)
    tau2 = (q - 2) / (300 - 30000 / 300)
    assert r["tau2"] == pytest.approx(tau2)
    assert r["pooled"] == pytest.approx(0.5)


def test_der_simonian_laird_homogeneous_and_degenerate():
    r = ere.der_simonian_laird(np.array([0.5, 0.5, 0.5]), np.array([0.01, 0.02, 0.03]))
    assert r["I2"] == 0 and r["tau2"] == 0
    assert ere.der_simonian_laird(np.array([0.5]), np.array([0.01])) == {}
    # zero or missing variances are dropped (observations with precision or recall of exactly 0 or 1)
    r = ere.der_simonian_laird(np.array([0.2, 0.5, 1.0, np.nan]), np.array([0.01, 0.01, 0.0, 0.01]))
    assert r["k"] == 2


def test_bootstrap_f1_variance_is_deterministic_and_positive():
    a = ere.bootstrap_f1_variance(40, 10, 20, seed=1)
    assert a == ere.bootstrap_f1_variance(40, 10, 20, seed=1) and a > 0
    assert np.isnan(ere.bootstrap_f1_variance(0, 0, 5, seed=1))


def test_i2_summary_numbers():
    from conftest import RESULTS

    text = (RESULTS / "tables" / "random_effects_summary.md").read_text()
    assert "F1 I² exceeded 90% for 17 of 17" in text
    assert "lower for real datasets in 15 of 15 tools" in text
