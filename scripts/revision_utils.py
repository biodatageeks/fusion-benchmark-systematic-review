"""Shared helpers for the analyses added in response to the BIB-26-2487 reviews.

Raw extraction files (data/raw/) are left untouched; the analysis input is built by build_analysis_input.py.
Revision-specific curation lives in data/curation/:

- bio_dataset_map.csv       biological-dataset identifier (shared when the same samples or
                            simulated reads are reused by several benchmarks) and real-data subtype.
- metadata_corrections.csv  documented corrections to benchmark-level metadata.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from real_simulated_sensitivity import DEFAULT_INPUT, load_table, normalize_columns

REPO = Path(__file__).resolve().parents[1]
CURATION = REPO / "data" / "curation"
BIO_MAP = REPO / "data" / "processed" / "dataset_annotations.csv"  # built by build_analysis_input.py
CORRECTIONS = CURATION / "metadata_corrections.csv"
REVISION_DIR = REPO / "results" / "revision"
METRICS = ["F1", "Precision", "Recall"]


def benchmark_labels() -> pd.Series:
    """Benchmark number (str) -> 'Author Year'."""
    meta = load_table(DEFAULT_INPUT, sheet_name="benchmark_metadata").dropna(subset=["benchmark_number"])
    numbers = pd.to_numeric(meta["benchmark_number"], errors="coerce").astype(int).astype(str)
    years = pd.to_numeric(meta["year_x"], errors="coerce").astype(int).astype(str)
    return pd.Series((meta["first_author"].astype(str) + " " + years).values, index=numbers.values)


def benchmark_years() -> pd.Series:
    meta = load_table(DEFAULT_INPUT, sheet_name="benchmark_metadata").dropna(subset=["benchmark_number"])
    numbers = pd.to_numeric(meta["benchmark_number"], errors="coerce").astype(int).astype(str)
    return pd.Series(pd.to_numeric(meta["year_x"], errors="coerce").astype(int).values, index=numbers.values)


def apply_metadata_corrections(meta: pd.DataFrame) -> pd.DataFrame:
    """Apply data/curation/metadata_corrections.csv to the benchmark_metadata sheet."""
    if not CORRECTIONS.exists():
        return meta
    meta = meta.copy()
    corrections = pd.read_csv(CORRECTIONS)
    numbers = pd.to_numeric(meta["benchmark_number"], errors="coerce")
    for _, row in corrections.iterrows():
        mask = numbers.eq(int(row["benchmark_number"]))
        if row["column"] not in meta.columns:
            raise KeyError(f"Correction refers to unknown column: {row['column']}")
        meta.loc[mask, row["column"]] = row["new_value"]
    return meta


def load_observations(include_spike_in: bool = False) -> pd.DataFrame:
    """Cleaned observation table (same cleaning as the original analyses) + revision identifiers."""
    raw = load_table(DEFAULT_INPUT)
    df = normalize_columns(raw)
    if include_spike_in:
        raise NotImplementedError("Spike-in observations are excluded by normalize_columns().")
    bio = pd.read_csv(BIO_MAP, dtype=str)
    df = df.merge(bio[["dataset_id", "bio_dataset_id", "real_subtype"]], left_on="Dataset", right_on="dataset_id", how="left")
    missing = df["bio_dataset_id"].isna()
    if missing.any():
        raise ValueError(f"Datasets missing from {BIO_MAP.name}: {sorted(df.loc[missing, 'Dataset'].unique())}")
    df["real"] = df["DatasetType"].eq("real").astype(int)
    df["Benchmark_label"] = df["Benchmark"].map(benchmark_labels())
    df["Benchmark_year"] = df["Benchmark"].map(benchmark_years())
    return df


def unit_counts(df: pd.DataFrame) -> dict[str, int]:
    return {
        "n_obs": int(len(df)),
        "n_tools": int(df["Tool"].nunique()),
        "n_bio_datasets": int(df["bio_dataset_id"].nunique()),
        "n_benchmarks": int(df["Benchmark"].nunique()),
    }


def squeeze(y: pd.Series) -> pd.Series:
    """Smithson-Verkuilen (2006) transformation mapping [0, 1] into (0, 1)."""
    n = y.notna().sum()
    return (y * (n - 1) + 0.5) / n


def logit(p: pd.Series | np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    return np.log(p / (1 - p))


def expit(x: float | np.ndarray) -> float | np.ndarray:
    return 1 / (1 + np.exp(-x))


def fmt(value: float, digits: int = 3) -> str:
    return "NA" if value is None or not np.isfinite(value) else f"{value:.{digits}f}"


def fmt_p(value: float) -> str:
    if value is None or not np.isfinite(value):
        return "NA"
    return f"{value:.2e}" if value < 0.001 else f"{value:.3f}"


def latex_escape(text: object) -> str:
    text = "" if pd.isna(text) else str(text)
    for old, new in [("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"), ("_", r"\_"), ("#", r"\#")]:
        text = text.replace(old, new)
    return text


def write_latex_table(df: pd.DataFrame, path: Path, caption: str, label: str, colspec: str | None = None,
                      font: str = r"\footnotesize", longtable: bool = False, resize: bool = False) -> None:
    """Minimal booktabs table writer (no pandas Styler / jinja2 dependency)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    colspec = colspec or "l" * len(df.columns)
    header = " & ".join(latex_escape(c) for c in df.columns) + r" \\"
    body = "\n".join(" & ".join(latex_escape(v) for v in row) + r" \\" for row in df.itertuples(index=False))
    if longtable:
        text = (
            f"{{{font}\n\\begin{{longtable}}{{{colspec}}}\n\\caption{{{caption}}}\\label{{{label}}}\\\\\n"
            f"\\toprule\n{header}\n\\midrule\n\\endfirsthead\n\\toprule\n{header}\n\\midrule\n\\endhead\n"
            f"{body}\n\\bottomrule\n\\end{{longtable}}\n}}\n"
        )
    else:
        text = (
            f"\\begin{{table}}[!htbp]\n\\centering\n{font}\n\\caption{{{caption}}}\n\\label{{{label}}}\n"
            + (r"\resizebox{\linewidth}{!}{%" + "\n" if resize else "")
            + f"\\begin{{tabular}}{{{colspec}}}\n\\toprule\n{header}\n\\midrule\n{body}\n\\bottomrule\n\\end{{tabular}}"
            + ("}" if resize else "") + "\n\\end{table}\n"
        )
    path.write_text(text)
