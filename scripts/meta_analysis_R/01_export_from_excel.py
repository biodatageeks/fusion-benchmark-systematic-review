"""
Export raw benchmark data from Benchmarki.xlsx into clean, versionable CSVs.

Reads: data/raw/Benchmarki.xlsx
Writes: data/processed/*.csv

Purpose:
  - Provides a stable, git-versionable snapshot of the master Excel spreadsheet
    used for the meta-analysis.
  - Cleans column headers, normalizes types, drops formula artifacts.
  - Reports data quality issues (missing FN, mixed types, etc.) so we know
    what needs to be addressed before running the meta-analysis proper.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_XLSX = REPO_ROOT / "data" / "raw" / "Benchmarki_Agnieszka.xlsx"
OUT_DIR = REPO_ROOT / "data" / "processed"

TOOL_NAME_FIXES = {
    "Chimerascan": "ChimeraScan",
    "pizzly": "Pizzly",
}


def _clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize column names: strip whitespace, fix corrupted first column."""
    new_cols = []
    for c in df.columns:
        c = str(c).strip()
        # The first column in "Sensitivity, Precision" is corrupted:
        # "B+A1:O373enchmark" from a botched paste. Normalize to "Benchmark".
        if "enchmark" in c and c != "Benchmark":
            c = "Benchmark"
        new_cols.append(c)
    df.columns = new_cols
    return df


def _strip_text_values(df: pd.DataFrame) -> pd.DataFrame:
    """Strip leading/trailing whitespace from text columns."""
    for col in df.select_dtypes(include=["object"]).columns:
        df[col] = df[col].map(lambda x: x.strip() if isinstance(x, str) else x)
    return df


def export_metrics() -> pd.DataFrame:
    """
    Master metrics sheet — one row per (benchmark, tool, dataset, run).
    """
    df = pd.read_excel(SRC_XLSX, sheet_name="Sensitivity, Precision")
    df = _clean_columns(df)
    df = _strip_text_values(df)
    if "Tool" in df.columns:
        df["Tool"] = df["Tool"].replace(TOOL_NAME_FIXES)

    # Coerce numeric columns that may be stored as object because of stray text.
    numeric_cols = [
        "TP", "FP", "TP+FP", "PP", "Recall/Sensitivity", "Precision", "F1",
        "Execution time (per run) [h]", "Memory usage (per sample)",
        "Time ranking in the set", "FN",
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Drop fully empty trailing rows.
    df = df.dropna(how="all")

    out = OUT_DIR / "benchmarks_metrics.csv"
    df.to_csv(out, index=False)
    return df


def export_dataset_details() -> pd.DataFrame:
    df = pd.read_excel(SRC_XLSX, sheet_name="Dataset details")
    df = _clean_columns(df)
    df = _strip_text_values(df)
    df = df.dropna(how="all")
    df.to_csv(OUT_DIR / "dataset_details.csv", index=False)
    return df


def export_reproducibility() -> pd.DataFrame:
    df = pd.read_excel(SRC_XLSX, sheet_name="Pre-processing reproducibility")
    df = _clean_columns(df)
    df = _strip_text_values(df)
    df = df.dropna(how="all")
    df.to_csv(OUT_DIR / "reproducibility.csv", index=False)
    return df


def export_gold_standard() -> pd.DataFrame:
    df = pd.read_excel(SRC_XLSX, sheet_name="Gold standard")
    df = _clean_columns(df)
    df = _strip_text_values(df)
    df = df.dropna(how="all")
    df.to_csv(OUT_DIR / "gold_standard.csv", index=False)
    return df


def export_tool_versions() -> pd.DataFrame:
    df = pd.read_excel(SRC_XLSX, sheet_name="Tool versions")
    df = _clean_columns(df)
    df = _strip_text_values(df)
    df = df.dropna(how="all")
    df.to_csv(OUT_DIR / "tool_versions.csv", index=False)
    return df


def export_tool_coverage() -> pd.DataFrame:
    df = pd.read_excel(SRC_XLSX, sheet_name="Tool statistics")
    df = _clean_columns(df)
    df = _strip_text_values(df)
    df = df.dropna(how="all")
    df.to_csv(OUT_DIR / "tool_coverage.csv", index=False)
    return df


def sanity_checks(metrics: pd.DataFrame) -> None:
    print("\n" + "=" * 70)
    print("SANITY CHECKS: benchmarks_metrics.csv")
    print("=" * 70)

    print(f"\nTotal rows: {len(metrics)}")
    print(f"Unique benchmarks: {metrics['Benchmark'].nunique()}")
    print(f"Unique tools: {metrics['Tool'].nunique()}")
    print(f"Unique datasets: {metrics['Dataset'].nunique()}")

    print("\nRows per benchmark:")
    print(metrics.groupby("Benchmark").size().to_string())

    print("\nColumn completeness (non-null count out of total):")
    for col in metrics.columns:
        nn = metrics[col].notna().sum()
        pct = 100 * nn / len(metrics)
        flag = "  MISSING >50%" if pct < 50 else ""
        print(f"  {col:40s} {nn:4d} / {len(metrics)} ({pct:5.1f}%){flag}")

    # Sanity: does F1 in file match 2*P*R/(P+R)?
    p = metrics["Precision"]
    r = metrics["Recall/Sensitivity"]
    f1_calc = 2 * p * r / (p + r)
    diff = (metrics["F1"] - f1_calc).abs()
    mask = metrics["F1"].notna() & f1_calc.notna()
    if mask.any():
        agree = (diff[mask] < 1e-6).sum()
        total = mask.sum()
        print(f"\nF1 formula check: {agree}/{total} rows agree with "
              f"2*P*R/(P+R) within 1e-6")
        mismatches = metrics[mask & (diff >= 1e-6)]
        if not mismatches.empty:
            print(f"Mismatched rows (showing up to 5):")
            print(mismatches[["Benchmark", "Tool", "Precision",
                              "Recall/Sensitivity", "F1"]].head())

    # Can we recompute precision from TP/(TP+FP)?
    tp = metrics["TP"]
    fp = metrics["FP"]
    mask_tpfp = tp.notna() & fp.notna() & ((tp + fp) > 0)
    if mask_tpfp.any():
        p_calc = tp[mask_tpfp] / (tp[mask_tpfp] + fp[mask_tpfp])
        p_stored = metrics.loc[mask_tpfp, "Precision"]
        both = p_calc.notna() & p_stored.notna()
        if both.any():
            diff_p = (p_calc[both] - p_stored[both]).abs()
            agree_p = (diff_p < 1e-3).sum()
            print(f"\nPrecision check: {agree_p}/{both.sum()} rows where "
                  f"TP/(TP+FP) agrees with stored Precision within 1e-3")

    # Rows missing critical fields for meta-analysis
    critical = ["TP", "FP", "FN"]
    missing_any = metrics[critical].isna().any(axis=1).sum()
    print(f"\nRows missing at least one of TP/FP/FN: {missing_any} "
          f"({100 * missing_any / len(metrics):.1f}%)")
    print("This is the ceiling for what can enter formal meta-analysis "
          "with proper variance estimation.")


def main() -> int:
    if not SRC_XLSX.exists():
        print(f"ERROR: source file not found: {SRC_XLSX}", file=sys.stderr)
        return 1

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Exporting sheets...")
    metrics = export_metrics()
    print(f"  benchmarks_metrics.csv  ({len(metrics)} rows)")
    ddet = export_dataset_details()
    print(f"  dataset_details.csv     ({len(ddet)} rows)")
    repro = export_reproducibility()
    print(f"  reproducibility.csv     ({len(repro)} rows)")
    gs = export_gold_standard()
    print(f"  gold_standard.csv       ({len(gs)} rows)")
    tv = export_tool_versions()
    print(f"  tool_versions.csv       ({len(tv)} rows)")
    tc = export_tool_coverage()
    print(f"  tool_coverage.csv       ({len(tc)} rows)")

    sanity_checks(metrics)
    return 0


if __name__ == "__main__":
    sys.exit(main())
