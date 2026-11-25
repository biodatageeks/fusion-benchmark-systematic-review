"""Shared fixtures for the analysis tests.

The tests read data/raw/ and the committed results/ only; nothing is written to the repository.
Run from the repository root:  python -m pytest
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
RESULTS = ROOT / "results"
MASTER = ROOT / "data" / "raw" / "master_analysis_input.xlsx"

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


@pytest.fixture(scope="session")
def book() -> dict[str, pd.DataFrame]:
    return pd.read_excel(MASTER, sheet_name=None)


@pytest.fixture(scope="session")
def obs(book) -> pd.DataFrame:
    return book["observations"]


@pytest.fixture(scope="session")
def datasets(book) -> pd.DataFrame:
    return book["dataset_metadata"]


@pytest.fixture(scope="session")
def benchmarks(book) -> pd.DataFrame:
    return book["benchmark_metadata"]


@pytest.fixture(scope="session")
def normalized() -> pd.DataFrame:
    """Observations as the real-vs-simulated analyses see them (real and simulated only)."""
    from real_simulated_sensitivity import load_table, normalize_columns

    return normalize_columns(load_table(MASTER))


@pytest.fixture(scope="session")
def results():
    """Read a committed results file: results('crossed_models/real_vs_simulated_models.csv')."""
    def read(relative: str) -> pd.DataFrame:
        return pd.read_csv(RESULTS / relative)

    return read
