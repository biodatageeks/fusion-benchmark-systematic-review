"""Edgren scoring: gene-symbol matching, directed vs undirected matching, one-to-one assignment, reference lists."""
from __future__ import annotations

import pandas as pd
import pytest

from conftest import ROOT
from edgren_rnafusion_validation import (
    evaluate_tool, extract_predictions, extract_truth, fusion_to_pair, gene_tokens, pair_key, pairs_match,
    undirected_pair_key,
)

TRUTH_DIR = ROOT / "data" / "raw" / "edgren_truthsets"
PRED_DIR = ROOT / "data" / "raw" / "edgren_predictions"
LISTS = [27, 38, 53, 99]
RUNS = {"SRR064286": "MCF-7", "SRR064287": "KPL-4", "SRR064438": "BT-474", "SRR064439": "BT-474",
        "SRR064440": "SK-BR-3", "SRR064441": "SK-BR-3"}


def truth_frame(pairs: list[tuple[str, str]]) -> pd.DataFrame:
    return extract_truth(pd.DataFrame(pairs, columns=["Gene1", "Gene2"]))


def predictions(fusions: list[str]) -> pd.DataFrame:
    return extract_predictions(pd.DataFrame({"Fusion": fusions, "call": True}), "tool", "call")


# --------------------------------------------------------------------------- gene symbols

def test_gene_tokens_aliases_alternatives_and_normalisation():
    assert gene_tokens("SEPT10") == gene_tokens("septin10") == frozenset({"SEPT10", "SEPTIN10"})
    assert gene_tokens(" vmp1 | TMEM49 ") == frozenset({"VMP1", "TMEM49"})
    assert gene_tokens(float("nan")) == frozenset() and gene_tokens("nan") == frozenset()


def test_fusion_to_pair():
    assert fusion_to_pair("BCAS4--BCAS3") == (frozenset({"BCAS4"}), frozenset({"BCAS3"}))
    assert fusion_to_pair("BCAS4-BCAS3") is None and fusion_to_pair("BCAS4--") is None
    assert fusion_to_pair(None) is None


def test_keys_merge_aliases_and_orientation():
    a, b = fusion_to_pair("GCN1L1--MSI1"), fusion_to_pair("MSI1--GCN1")
    assert pair_key(a) != pair_key(b)
    assert undirected_pair_key(a) == undirected_pair_key(b)
    assert pair_key(fusion_to_pair("GCN1--MSI1")) == pair_key(a)


def test_pairs_match_directional_and_undirected():
    truth = fusion_to_pair("RPS6KB1--VMP1|TMEM49")
    assert pairs_match(fusion_to_pair("RPS6KB1--TMEM49"), truth, directional=True)
    assert not pairs_match(fusion_to_pair("TMEM49--RPS6KB1"), truth, directional=True)
    assert pairs_match(fusion_to_pair("TMEM49--RPS6KB1"), truth, directional=False)
    assert not pairs_match(fusion_to_pair("RPS6KB1--SNF8"), truth, directional=False)


# --------------------------------------------------------------------------- evaluate_tool

def test_evaluate_tool_directional_vs_undirected():
    truth = truth_frame([("A", "B"), ("C", "D"), ("G", "H")])
    pred = predictions(["A--B", "B--A", "C|X--D", "E--F"])
    directional, *_ = evaluate_tool(pred, truth, directional=True)
    assert (directional["TP"], directional["FP"], directional["FN"]) == (2, 2, 1)
    undirected, *_ = evaluate_tool(pred, truth, directional=False)
    # B--A collapses with A--B before scoring
    assert undirected["unique_predictions"] == 3
    assert (undirected["TP"], undirected["FP"], undirected["FN"]) == (2, 1, 1)
    assert undirected["precision"] == pytest.approx(2 / 3) and undirected["recall"] == pytest.approx(2 / 3)


def test_evaluate_tool_one_to_one_assignment():
    """Two distinct predictions matching the same reference event give one TP and one FP."""
    truth = truth_frame([("A", "B")])
    summary, tp, fp, fn = evaluate_tool(predictions(["A--B", "A|Z--B"]), truth, directional=True)
    assert (summary["TP"], summary["FP"], summary["FN"]) == (1, 1, 0)
    assert len(tp) == 1 and len(fp) == 1 and fn.empty


def test_evaluate_tool_empty_overlap():
    summary, *_ = evaluate_tool(predictions(["X--Y"]), truth_frame([("A", "B")]), directional=False)
    assert (summary["TP"], summary["FP"], summary["FN"], summary["precision"], summary["recall"]) == (0, 1, 1, 0, 0)


# --------------------------------------------------------------------------- reference lists

@pytest.fixture(scope="module")
def manifest() -> pd.DataFrame:
    return pd.read_csv(TRUTH_DIR / "manifest.csv")


def read_list(size: int, run: str) -> pd.DataFrame:
    return pd.read_csv(TRUTH_DIR / f"truthset_{size}" / f"{run}_truthset_{size}.csv", dtype=str)


def test_manifest_matches_files_and_list_sizes(manifest):
    assert set(manifest["sra_run"]) == set(RUNS)
    assert dict(zip(manifest["sra_run"], manifest["cell_line"])) == RUNS
    for size, g in manifest.groupby("truthset_size"):
        per_cell = g.drop_duplicates("cell_line")
        assert per_cell["events_for_cell_line"].sum() == size
        for r in g.itertuples():
            t = read_list(size, r.sra_run)
            assert len(t) == r.events_for_cell_line
            assert set(t["cell_line"]) == {r.cell_line} and set(t["sra_run"]) == {r.sra_run}
    events = manifest.drop_duplicates(["cell_line", "truthset_size"]).pivot(
        index="cell_line", columns="truthset_size", values="events_for_cell_line")
    assert events.loc["KPL-4", 99] == 3 and events.loc["MCF-7", 99] == 52  # Results 3.5: 3 in KPL-4 to 52 in MCF-7


def test_runs_of_one_cell_line_share_the_reference_list():
    """edgren_per_sample_scoring.py scores a cell line against the list of its first run."""
    for size in LISTS:
        for a, b in [("SRR064438", "SRR064439"), ("SRR064440", "SRR064441")]:
            x, y = (read_list(size, r).drop(columns="sra_run").reset_index(drop=True) for r in (a, b))
            pd.testing.assert_frame_equal(x, y)


def test_no_list_contains_the_same_event_twice():
    """Undirected scoring keeps one prediction per gene pair; a list with A--B and B--A (or alias duplicates)
    in one cell line would make one of the two events unrecoverable."""
    for size in LISTS:
        for run in ["SRR064286", "SRR064287", "SRR064438", "SRR064440"]:
            t = read_list(size, run)
            keys = [undirected_pair_key((gene_tokens(a), gene_tokens(b))) for a, b in zip(t.gene1_raw, t.gene2_raw)]
            assert len(keys) == len(set(keys)), (size, run)


def test_pooled_99_truth_equals_per_run_lists():
    xlsx = extract_truth(pd.read_excel(ROOT / "data" / "raw" / "edgren_rnafusion.xlsx"))
    runs = ["SRR064286", "SRR064287", "SRR064438", "SRR064440"]
    files = pd.concat([read_list(99, r) for r in runs]).rename(columns={"gene1_raw": "Gene1", "gene2_raw": "Gene2"})
    assert len(xlsx) == 99 and set(xlsx["truth_key"]) == set(extract_truth(files)["truth_key"])


# --------------------------------------------------------------------------- caller outputs

def test_arriba_flanking_genes():
    from edgren_per_sample_scoring import arriba_genes

    assert arriba_genes("BCAS1") == ["BCAS1"]
    assert arriba_genes("B4GALNT2(1234),B4GALNT2P1(567)") == ["B4GALNT2", "B4GALNT2P1"]


def test_prediction_files_match_manifest():
    m = pd.read_csv(PRED_DIR / "manifest.csv")
    assert len(m) == 36
    for r in m.itertuples():
        assert len(pd.read_csv(PRED_DIR / r.csv_path, dtype=str)) == r.prediction_rows, r.csv_path


def test_pooled_scoring_reproduces_committed_summary(results):
    """Results 3.5 numbers (pooled RNAfusion 3.0.1 calls vs the 99-fusion list)."""
    from edgren_rnafusion_validation import TOOL_COLUMNS

    df = pd.read_excel(ROOT / "data" / "raw" / "edgren_rnafusion.xlsx")
    truth = extract_truth(df)
    committed = results("edgren_rnafusion_validation/edgren_rnafusion_metrics_summary.csv").set_index(["tool", "matching"])
    for tool, column in TOOL_COLUMNS.items():
        pred = extract_predictions(df, tool, column)
        for directional in (True, False):
            s, *_ = evaluate_tool(pred, truth, directional=directional)
            c = committed.loc[(tool, s["matching"])]
            assert (s["unique_predictions"], s["TP"], s["FP"], s["FN"]) == \
                (c.unique_predictions, c.TP, c.FP, c.FN), (tool, s["matching"])
    und = committed.xs("undirected", level="matching").round(3)
    assert und.loc["Arriba", ["TP", "FP", "precision", "recall", "F1"]].tolist() == [49, 133, 0.269, 0.495, 0.349]
    # MED1--UPS32 of the 99-fusion list matches the USP32--MED1 calls of FusionCatcher and STAR-Fusion (gene alias)
    assert und.loc["FusionCatcher", ["TP", "FP", "precision", "recall", "F1"]].tolist() == [40, 8, 0.833, 0.404, 0.544]
    assert und.loc["STAR-Fusion", ["TP", "FP", "precision", "recall", "F1"]].tolist() == [47, 29, 0.618, 0.475, 0.537]
    overlap = results("edgren_rnafusion_validation/edgren_rnafusion_overlap_summary.csv").set_index("comparison")["truth_fusions"]
    assert (overlap["any tool"], overlap["no tool"]) == (52, 47)


def test_per_run_files_reproduce_pooled_scoring(results):
    """The raw per-run caller outputs pooled over runs give the same counts as the fusion-report workbook."""
    check = results("edgren_per_sample/sanity_check_vs_pooled_scoring.csv")
    for m in ["unique_predictions", "TP", "FP", "FN"]:
        assert (check[f"{m}_here"] == check[f"{m}_pooled"]).all(), m


def test_per_run_scores_reproduce(results):
    from edgren_per_sample_scoring import score

    committed = results("edgren_per_sample/per_run_scores.csv")
    sub = committed[(committed.rnafusion == "3.0.1") & (committed.truthset == 99)]
    for run in ["SRR064286", "SRR064438"]:
        for tool in ["Arriba", "FusionCatcher", "STAR-Fusion"]:
            for r in score("3.0.1", tool, [run], [run], 99):
                c = sub[(sub.sra_run == run) & (sub.tool == tool) & (sub.matching == r["matching"])].iloc[0]
                assert (r["TP"], r["FP"], r["FN"]) == (c.TP, c.FP, c.FN)


# --------------------------------------------------------------------------- Supplementary Table S6 event rule

def test_same_event_rule():
    from supplementary_tables import same_event

    def rec(cell, a, b):
        return {"cell_line": cell, "a": gene_tokens(a), "b": gene_tokens(b)}

    assert same_event(rec("BT-474", "GCN1L1", "MSI1"), rec("BT-474", "MSI1", "GCN1"))
    assert same_event(rec("MCF-7", "RPS6KB1", "VMP1|TMEM49"), rec("MCF-7", "RPS6KB1", "TMEM49"))
    assert not same_event(rec("MCF-7", "A", "B"), rec("BT-474", "A", "B"))
    assert not same_event(rec("MCF-7", "A", "B"), rec("MCF-7", "A", "C"))
