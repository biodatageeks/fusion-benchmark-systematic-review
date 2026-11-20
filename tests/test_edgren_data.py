"""Integrity of the Edgren inputs and regression of the reported Edgren results (repository data)."""
import pandas as pd
import pytest

import edgren_reference_sensitivity as E

LIST_SIZES = {"TP27": 27, "TP53": 53, "TP99": 99}
PER_CELL_LINE = {  # events per cell line in the published lists
    "TP27": {"BT-474": 11, "SK-BR-3": 10, "KPL-4": 3, "MCF-7": 3},
    "TP53": {"BT-474": 23, "SK-BR-3": 12, "KPL-4": 3, "MCF-7": 15},
    "TP99": {"BT-474": 21, "SK-BR-3": 23, "KPL-4": 3, "MCF-7": 52},
}


@pytest.fixture(scope="module")
def ref():
    r = E.reference_table(pd.read_csv(E.DATA / "reference_events.csv"))
    return r[r.list_id.isin(E.LISTS)].reset_index(drop=True)


@pytest.fixture(scope="module")
def manifest():
    return E.sample_manifest()


def test_reference_list_sizes(ref):
    for list_id, n in LIST_SIZES.items():
        assert len(ref[ref.list_id == list_id]) == n
        assert ref[ref.list_id == list_id].groupby("cell_line").size().to_dict() == PER_CELL_LINE[list_id]


def test_reference_nesting(ref):
    ov = E.overlap_table(ref).set_index(["list_a", "list_b"])
    assert ov.loc[("TP27", "TP53"), "a_nested_in_b"] and ov.loc[("TP27", "TP99"), "a_nested_in_b"]
    assert ov.loc[("TP53", "TP99"), "shared"] == 49
    assert ov.loc[("TP53", "TP99"), "only_in_a"] == 4 and ov.loc[("TP99", "TP53"), "only_in_a"] == 50


def test_reference_events_agree_with_authors_per_run_sets(ref):
    chk = E.check_reference(ref)
    unmatched = chk[~chk.matched]
    # the only difference is the gene label MED1--UPS32 (source typo) versus MED1--USP32
    assert sorted(unmatched.truthset_event.tolist() + unmatched.reference_event.tolist()) == ["", "", "MED1--UPS32", "MED1--USP32"]


def test_manifest_matches_authors_cell_lines(manifest):
    iga = pd.read_csv(E.TRUTHSETS / "manifest.csv").drop_duplicates("sra_run").set_index("sra_run").cell_line
    assert len(manifest) == 6
    assert all(iga[r] == c for r, c in zip(manifest.run_accession, manifest.cell_line))


def test_prediction_files_complete_and_row_counts_match_manifest():
    m = pd.read_csv(E.PRED / "manifest.csv")
    assert len(m) == 36
    for _, row in m.iterrows():
        assert len(pd.read_csv(E.PRED / row.csv_path)) == row.prediction_rows, row.csv_path


def test_raw_301_outputs_reproduce_submitted_pooled_workbook(manifest):
    chk = E.check_raw_against_workbook(manifest)
    assert chk.identical.all(), chk


def test_submitted_pooled_tp99_result_is_reproduced():
    table = pd.read_excel(E.LEGACY_INPUT, sheet_name="All")
    ref = E.reference_table(table[["Gene1", "Gene2"]].rename(columns={"Gene1": "gene1", "Gene2": "gene2"})
                            .assign(list_id="TP99", cell_line=None))
    got = {t: E.score(E.predictions_from_table(table, t, c), ref) for t, c in E.TOOL_COLUMNS.items()}
    assert {t: (s["TP"], s["FP"], s["FN"]) for t, s in got.items()} == {
        "Arriba": (49, 133, 50), "FusionCatcher": (39, 9, 60), "STAR-Fusion": (46, 30, 53)}


@pytest.fixture(scope="module")
def per_sample(manifest, ref):
    return E.aggregate(E.per_unit_metrics(manifest, ref, "sample")).set_index(["workflow_version", "tool", "list_id"])


@pytest.mark.parametrize("version,tool,expected", [
    ("4.1.3", "STAR-Fusion", (0.552, 0.678, 0.576)),
    ("4.1.3", "FusionCatcher", (0.472, 0.539, 0.450)),
    ("4.1.3", "Arriba", (0.313, 0.415, 0.394)),
    ("3.0.1", "STAR-Fusion", (0.573, 0.716, 0.588)),
    ("3.0.1", "FusionCatcher", (0.574, 0.619, 0.505)),
    ("3.0.1", "Arriba", (0.300, 0.422, 0.414)),
])
def test_per_sample_f1_reported_in_manuscript(per_sample, version, tool, expected):
    got = tuple(per_sample.loc[(version, tool, l), "F1"] for l in E.LISTS)
    assert got == pytest.approx(expected, abs=5e-4)


def test_precision_rises_and_recall_falls_with_list_size(per_sample):
    for (version, tool), g in per_sample.groupby(level=[0, 1]):
        p = [g.loc[(version, tool, l), "precision"] for l in E.LISTS]
        r = [g.loc[(version, tool, l), "recall"] for l in E.LISTS]
        assert p == sorted(p) and r == sorted(r, reverse=True), (version, tool)


def test_no_call_matches_only_a_reference_of_another_cell_line(manifest, ref):
    from edgren_rnafusion_validation import pairs_match
    for version in E.VERSIONS:
        for tool in E.TOOL_COLUMNS:
            for cl, g in manifest.groupby("cell_line"):
                p = E.predictions_from_fusions(sum((E.raw_fusions(version, r, tool) for r in g.run_accession), []), tool)
                own = ref[ref.cell_line == cl].truth_pair.tolist()
                other = ref[ref.cell_line != cl].truth_pair.tolist()
                for pp in p.pred_pair:
                    if not any(pairs_match(pp, q, False) for q in own):
                        assert not any(pairs_match(pp, q, False) for q in other), (version, tool, cl, pp)
