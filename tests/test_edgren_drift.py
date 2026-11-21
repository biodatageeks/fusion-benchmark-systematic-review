"""Figure 6: published Singh and Li 2021 scores of the same Edgren predictions against the 27- and 99-fusion lists."""
import numpy as np
import pandas as pd
import pytest

from conftest import REPO

INPUT = REPO / "results" / "tables" / "random_effects_input_effect_sizes.csv"
DRIFT = REPO / "results" / "tables" / "edgren_truth_set_drift.csv"


@pytest.fixture(scope="module")
def drift():
    from exploratory_random_effects import edgren_truth_set_drift

    return edgren_truth_set_drift(pd.read_csv(INPUT)).set_index("Tool")


def test_compares_the_two_edgren_lists(drift):
    # TP + FN is the reference size: the two columns must be the 27- and the 99-fusion lists of Edgren
    assert np.allclose(drift["TP_27"] + drift["FN_27"], 27, atol=0.5)
    assert np.allclose(drift["TP_99"] + drift["FN_99"], 99, atol=0.5)


def test_same_predictions_under_both_lists(drift):
    n27, n99 = drift["TP_27"] + drift["FP_27"], drift["TP_99"] + drift["FP_99"]
    rel = (n99 - n27).abs() / np.maximum(n27, n99)
    # counts are back-calculated from two-decimal precision; for ChimeraScan (precision 0.03 vs 0.04)
    # the rounding alone moves the implied prediction count by ~40%
    assert (rel.drop("ChimeraScan") <= 0.05).all()
    assert (drift["TP_99"] >= drift["TP_27"] - 0.5).all()


def test_values_reported_in_section_3_5(drift):
    assert len(drift) == 16
    # F1 recomputed from the rounded precision and recall gives 0.214; the reported F1 values give 0.22
    assert drift["delta_f1"].abs().max() == pytest.approx(0.214, abs=0.005)
    assert drift["delta_f1"].abs().idxmax() == "MapSplice"
    assert int((drift["delta_f1"] < -1e-9).sum()) == 9
    assert drift.loc["FusionCatcher", ["f1_27", "f1_99"]].tolist() == pytest.approx([0.56, 0.38], abs=0.01)
    assert drift.loc["STAR-Fusion", ["f1_27", "f1_99"]].tolist() == pytest.approx([0.44, 0.47], abs=0.01)


def test_saved_table_matches(drift):
    saved = pd.read_csv(DRIFT).set_index("Tool").sort_index()
    pd.testing.assert_frame_equal(saved, drift.sort_index(), check_exact=False, rtol=1e-9)


def test_tradeoff_edgren_flag_follows_the_data():
    """Supplementary tool profiles without Edgren: rows are excluded by the is_edgren flag, not by dataset number."""
    from precision_recall_tradeoff import load_clean_data
    from real_simulated_sensitivity import DEFAULT_INPUT, load_table

    raw = load_table(DEFAULT_INPUT)
    flagged = raw["is_edgren"].astype(str).str.lower().isin(["true", "yes", "1"])
    df = load_clean_data()
    assert int(df["Is_Edgren"].sum()) == int(flagged.sum()) > 0
    singh = df[df["DatasetID"].isin(["7.5", "7.6", "7.7"])].groupby("DatasetID")["Is_Edgren"].all().to_dict()
    assert singh == {"7.5": True, "7.6": True, "7.7": False}  # 7.7 is Qin 2015 LNCaP, not Edgren
