"""
Deeper data-quality inspection of benchmarks_metrics.csv.

Checks whether missing fields can be recovered algebraically:
  - FN from TP and Sensitivity: FN = TP*(1/S - 1)
  - N_true (total true fusions) from TP and Sensitivity: N_true = TP/S
  - Precision from TP/(TP+FP), Sensitivity from TP/(TP+FN)

Also inspects per-benchmark data completeness so we know which benchmarks
can enter the formal meta-analysis (with variance estimation) and which
can only enter as descriptive stats.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
CSV = REPO_ROOT / "data" / "processed" / "benchmarks_metrics.csv"


def main() -> None:
    df = pd.read_csv(CSV)

    print("=" * 70)
    print("PER-BENCHMARK COMPLETENESS OF CRITICAL FIELDS")
    print("=" * 70)

    fields = ["TP", "FP", "FN", "Precision", "Recall/Sensitivity", "F1"]
    per_bm = df.groupby("Benchmark")[fields].apply(
        lambda g: (g.notna().sum() / len(g) * 100).round(1)
    )
    print(per_bm.to_string())

    print("\n" + "=" * 70)
    print("CAN WE RECOVER FN FROM TP AND SENSITIVITY?")
    print("=" * 70)

    tp = df["TP"]
    sens = df["Recall/Sensitivity"]
    mask = tp.notna() & sens.notna() & (sens > 0)
    print(f"Rows where TP and Sensitivity are both available: "
          f"{mask.sum()} / {len(df)}")

    fn_recovered = np.where(mask, tp * (1.0 / sens - 1.0), np.nan)
    df["FN_recovered"] = fn_recovered
    df["N_true_recovered"] = np.where(mask, tp / sens, np.nan)

    print("Distribution of recovered FN (percentiles):")
    fn_vals = pd.Series(fn_recovered).dropna()
    if not fn_vals.empty:
        print(f"  min:    {fn_vals.min():.2f}")
        print(f"  25%:    {fn_vals.quantile(0.25):.2f}")
        print(f"  median: {fn_vals.median():.2f}")
        print(f"  75%:    {fn_vals.quantile(0.75):.2f}")
        print(f"  max:    {fn_vals.max():.2f}")

    # Non-integer FN would be suspicious (FN must be a whole count of fusions).
    non_int = fn_vals[np.abs(fn_vals - fn_vals.round()) > 0.5].shape[0]
    print(f"Recovered FN values that are >0.5 away from integer: {non_int}")
    print("(High count → Sensitivity was reported rounded, so N_true has "
          "some uncertainty. Low count → we can recover FN cleanly.)")

    print("\n" + "=" * 70)
    print("F1 CONSISTENCY WITH TOLERANCE 1e-3 (was 1e-6)")
    print("=" * 70)

    p = df["Precision"]
    r = df["Recall/Sensitivity"]
    f1_stored = df["F1"]
    f1_calc = 2 * p * r / (p + r)
    mask_f1 = p.notna() & r.notna() & f1_stored.notna()
    diff = (f1_stored - f1_calc).abs()
    for tol in [1e-6, 1e-3, 1e-2]:
        agree = ((diff < tol) & mask_f1).sum()
        print(f"  tolerance {tol}: {agree} / {mask_f1.sum()} agree")

    # Rows where F1 disagrees by more than 1e-3 — these need a look.
    bad = df[mask_f1 & (diff >= 1e-3)]
    if not bad.empty:
        print(f"\nRows where F1 differs from 2PR/(P+R) by ≥ 1e-3: {len(bad)}")
        print("First 10:")
        print(bad[["Benchmark", "Tool", "Precision", "Recall/Sensitivity",
                   "F1"]].head(10).to_string())

    print("\n" + "=" * 70)
    print("SUMMARY: WHAT'S USABLE FOR META-ANALYSIS")
    print("=" * 70)

    # A row is fully usable if it has TP, and either FP+FN or Precision+Sensitivity
    has_tp_fp = df["TP"].notna() & df["FP"].notna()
    has_tp_fn = df["TP"].notna() & df["FN"].notna()
    has_tp_sens = df["TP"].notna() & df["Recall/Sensitivity"].notna()
    has_tp_prec = df["TP"].notna() & df["Precision"].notna()

    print(f"Rows with TP and FP:                  {has_tp_fp.sum()} "
          f"({100*has_tp_fp.mean():.1f}%)")
    print(f"Rows with TP and FN (stored):         {has_tp_fn.sum()} "
          f"({100*has_tp_fn.mean():.1f}%)")
    print(f"Rows with TP and Sensitivity:         {has_tp_sens.sum()} "
          f"({100*has_tp_sens.mean():.1f}%)")
    print(f"Rows with TP and Precision:           {has_tp_prec.sum()} "
          f"({100*has_tp_prec.mean():.1f}%)")

    # Rows with enough to reconstruct variance for meta-analysis:
    # precision variance needs TP + FP (or TP + Precision)
    # sensitivity variance needs TP + FN (stored or recoverable from Sens)
    can_do_prec_var = has_tp_fp | has_tp_prec
    can_do_sens_var = has_tp_fn | has_tp_sens
    can_do_both = can_do_prec_var & can_do_sens_var
    print(f"\nRows where precision variance can be estimated: "
          f"{can_do_prec_var.sum()} ({100*can_do_prec_var.mean():.1f}%)")
    print(f"Rows where sensitivity variance can be estimated: "
          f"{can_do_sens_var.sum()} ({100*can_do_sens_var.mean():.1f}%)")
    print(f"Rows where BOTH can be estimated: "
          f"{can_do_both.sum()} ({100*can_do_both.mean():.1f}%)")

    # Save enriched CSV with recovered fields for downstream use.
    out = REPO_ROOT / "data" / "processed" / "benchmarks_metrics_enriched.csv"
    df.to_csv(out, index=False)
    print(f"\nSaved enriched file: {out.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
