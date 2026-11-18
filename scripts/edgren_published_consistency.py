"""Were published TP27 and TP99 Edgren scores computed on the same prediction sets? (R3.3)

Two included benchmarks report the same tools on the Edgren data under both the 27-fusion and the
99-fusion reference lists: Vu et al. 2018 (datasets 6.1/6.2) and Singh & Li 2021 (7.6/7.7).
If one prediction set was re-scored, the number of predictions (TP + FP) must be identical under both
lists, and, if the 27 list is nested in the 99 list, TP99 >= TP27 for every tool.
TP and FP of Singh & Li 2021 are back-calculated from rounded precision and recall, so agreement is
assessed with a relative tolerance.
"""

from __future__ import annotations

import pandas as pd

from revision_utils import REVISION_DIR, load_observations

OUT = REVISION_DIR / "edgren_published"
TOLERANCE = 0.05  # relative difference in TP + FP accepted as rounding


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    df = load_observations()
    ed = df[(df["bio_dataset_id"] == "EDGREN_BRCA_CELL_LINES") & df["TruthTotal"].isin([27, 99])]
    pairs = []
    for bench, g in ed.groupby("Benchmark"):
        d27 = sorted(g.loc[g.TruthTotal == 27, "Dataset"].unique()); d99 = sorted(g.loc[g.TruthTotal == 99, "Dataset"].unique())
        if len(d27) == 1 and len(d99) == 1:
            pairs.append((d27[0], d99[0]))
    rows = []
    for d27, d99 in pairs:
        a = df[df.Dataset == d27].drop_duplicates("Tool").set_index("Tool")
        b = df[df.Dataset == d99].drop_duplicates("Tool").set_index("Tool")
        for tool in sorted(set(a.index) & set(b.index)):
            tp27, fp27, tp99, fp99 = a.at[tool, "TP"], a.at[tool, "FP"], b.at[tool, "TP"], b.at[tool, "FP"]
            n27, n99 = tp27 + fp27, tp99 + fp99
            rel = abs(n99 - n27) / max(n27, n99) if max(n27, n99) > 0 else float("nan")
            rows.append({
                "benchmark": a.at[tool, "Benchmark_label"], "tool": tool,
                "TP27": tp27, "FP27": fp27, "predictions27": n27, "F1_27": a.at[tool, "F1"],
                "TP99": tp99, "FP99": fp99, "predictions99": n99, "F1_99": b.at[tool, "F1"],
                "delta_F1": b.at[tool, "F1"] - a.at[tool, "F1"],
                "rel_diff_predictions": rel,
                "same_prediction_count": bool(rel <= TOLERANCE) if pd.notna(rel) else None,
                "TP99_ge_TP27": bool(tp99 >= tp27 - 0.5),
            })
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "tp27_vs_tp99_consistency.csv", index=False)
    lines = ["# Published TP27 vs TP99 scores on Edgren data", ""]
    for bench, g in res.groupby("benchmark"):
        g = g[g["predictions27"] > 0]
        bad = g.loc[g["same_prediction_count"] == False, "tool"].tolist()  # noqa: E712
        lines.append(
            f"- {bench}: {len(g)} tools; prediction count identical within {TOLERANCE:.0%} for "
            f"{int(g['same_prediction_count'].sum())}/{len(g)}"
            + (f" (exceptions: {', '.join(bad)})" if bad else "")
            + f"; TP99 >= TP27 for {int(g['TP99_ge_TP27'].sum())}/{len(g)}; "
            f"|delta F1| range {g['delta_F1'].abs().min():.3f} to {g['delta_F1'].abs().max():.3f}."
        )
    (OUT / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(res.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
