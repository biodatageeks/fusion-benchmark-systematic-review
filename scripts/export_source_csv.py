"""Export every sheet of the analysis workbook to CSV (results/source_csv/), so the source data can be read without Excel."""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "raw" / "master_analysis_input.xlsx"
OUTPUT_DIR = ROOT / "results" / "source_csv"


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for name, sheet in pd.read_excel(INPUT, sheet_name=None).items():
        path = OUTPUT_DIR / f"{name}.csv"
        sheet.to_csv(path, index=False)
        print(f"{path.name}: {sheet.shape[0]} rows x {sheet.shape[1]} columns")


if __name__ == "__main__":
    main()
