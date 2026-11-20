"""Scoring rules of the Edgren re-scoring on small synthetic examples (no repository data needed)."""
import pandas as pd
import pytest

import edgren_reference_sensitivity as E
from edgren_rnafusion_validation import fusion_to_pair, gene_tokens


def pred(*fusions):
    return E.predictions_from_fusions(list(fusions), "T")


def truth(*pairs, cell_line="X", list_id="TP99"):
    ev = pd.DataFrame([{"gene1": a, "gene2": b, "cell_line": cell_line, "list_id": list_id} for a, b in pairs])
    return E.reference_table(ev)


def test_undirected_matching_and_counts():
    s = E.score(pred("A--B", "C--D", "E--F"), truth(("B", "A"), ("C", "D"), ("G", "H")))
    assert (s["TP"], s["FP"], s["FN"]) == (2, 1, 1)
    assert s["precision"] == pytest.approx(2 / 3)
    assert s["recall"] == pytest.approx(2 / 3)
    assert s["F1"] == pytest.approx(2 / 3)


def test_duplicates_collapse_and_reverse_orientation_is_one_prediction():
    s = E.score(pred("A--B", "A--B", "B--A"), truth(("A", "B")))
    assert (s["TP"], s["FP"], s["FN"]) == (1, 0, 0)


def test_alias_alternatives_match_either_symbol():
    assert gene_tokens("TTI1|KIAA0406") == frozenset({"TTI1", "KIAA0406"})
    s = E.score(pred("DIDO1--KIAA0406"), truth(("DIDO1", "TTI1|KIAA0406")))
    assert (s["TP"], s["FP"], s["FN"]) == (1, 0, 0)


def test_one_reference_event_is_matched_at_most_once():
    # two different predictions hitting the same event: one TP, the other is reference-unmatched
    s = E.score(pred("SKA2--MYO19", "TRIM37--MYO19"), truth(("SKA2|TRIM37", "MYO19")))
    assert (s["TP"], s["FP"], s["FN"]) == (1, 1, 0)


def test_empty_inputs():
    assert E.score(pred(), truth(("A", "B")))["FN"] == 1
    assert E.score(pred("A--B"), truth().iloc[0:0] if False else E.reference_table(
        pd.DataFrame(columns=["gene1", "gene2", "cell_line", "list_id"])))["FP"] == 1


def test_arriba_gene_field_parsing():
    assert E._arriba_genes("BCAS4") == ["BCAS4"]
    assert E._arriba_genes("ENSG1(36653),ENSG2(198826)") == ["ENSG1", "ENSG2"]
    assert E._arriba_genes(".") == []


@pytest.fixture
def toy_predictions(tmp_path, monkeypatch):
    """Two runs of one cell line X and one run of cell line Y, one workflow version."""
    def write(run, tool, rows):
        d = tmp_path / "v1" / run
        d.mkdir(parents=True, exist_ok=True)
        if tool == "Arriba":
            df = pd.DataFrame(rows, columns=["#gene1", "gene2"])
        elif tool == "FusionCatcher":
            df = pd.DataFrame(rows, columns=["Gene_1_symbol(5end_fusion_partner)", "Gene_2_symbol(3end_fusion_partner)"])
        else:
            df = pd.DataFrame({"#FusionName": [f"{a}--{b}" for a, b in rows]})
        df.to_csv(d / f"{tool}.csv", index=False)

    for run in ["R1", "R2", "R3"]:
        write(run, "FusionCatcher", [])
        write(run, "STAR-Fusion", [])
    write("R1", "Arriba", [("A", "B"), ("C", "D")])                  # X: A--B true in X; C--D true only in Y
    write("R2", "Arriba", [("A", "B"), ("N1(10),N2(20)", "M")])      # X: repeated true call; multi-gene label
    write("R3", "Arriba", [("C", "D")])                              # Y: C--D true in Y
    monkeypatch.setattr(E, "PRED", tmp_path)
    monkeypatch.setattr(E, "VERSIONS", ["v1"])
    manifest = pd.DataFrame({"run_accession": ["R1", "R2", "R3"], "cell_line": ["X", "X", "Y"]})
    ref = pd.concat([truth(("A", "B"), ("E", "F"), cell_line="X"), truth(("C", "D"), cell_line="Y")], ignore_index=True)
    return manifest, ref


def test_per_sample_scoring_uses_own_cell_line_only(toy_predictions):
    manifest, ref = toy_predictions
    by = E.per_unit_metrics(manifest, ref, "sample")
    a = by[(by.tool == "Arriba") & (by.list_id == "TP99")].set_index("unit")
    # R1: A--B TP, C--D is true only in cell line Y -> FP; E--F missed
    assert tuple(a.loc["R1", ["TP", "FP", "FN"]]) == (1, 1, 1)
    # R2: A--B TP again (counted per run); multi-gene label expanded into N1--M and N2--M -> 2 FP
    assert tuple(a.loc["R2", ["TP", "FP", "FN"]]) == (1, 2, 1)
    assert tuple(a.loc["R3", ["TP", "FP", "FN"]]) == (1, 0, 0)
    micro = E.aggregate(by)
    m = micro[(micro.tool == "Arriba") & (micro.list_id == "TP99")].iloc[0]
    assert (m.TP, m.FP, m.FN) == (3, 3, 2)
    assert m.precision == pytest.approx(0.5) and m.recall == pytest.approx(0.6)


def test_cell_line_unit_takes_union_of_runs(toy_predictions):
    manifest, ref = toy_predictions
    by = E.per_unit_metrics(manifest, ref, "cell_line")
    a = by[(by.tool == "Arriba") & (by.list_id == "TP99")].set_index("unit")
    assert tuple(a.loc["X", ["TP", "FP", "FN"]]) == (1, 3, 1)   # A--B once; C--D, N1--M, N2--M


def test_arriba_alternatives_mode_counts_multi_gene_label_once(toy_predictions):
    manifest, ref = toy_predictions
    by = E.per_unit_metrics(manifest, ref, "sample", arriba_labels="alternatives")
    a = by[(by.tool == "Arriba") & (by.list_id == "TP99")].set_index("unit")
    assert tuple(a.loc["R2", ["TP", "FP", "FN"]]) == (1, 1, 1)


def test_fusion_label_parsing():
    assert fusion_to_pair("A--B") == (frozenset({"A"}), frozenset({"B"}))
    assert fusion_to_pair("no separator") is None
