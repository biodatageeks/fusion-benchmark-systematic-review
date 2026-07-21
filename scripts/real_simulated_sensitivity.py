from __future__ import annotations

import argparse
import math
import re
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf


DEFAULT_INPUT = Path(__file__).resolve().parents[1] / "data" / "raw" / "real_simulated.xlsx"
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parents[1] / "results" / "real_simulated"


def _column_index(cell_ref: str) -> int:
    letters = re.sub(r"[^A-Z]", "", cell_ref.upper())
    index = 0
    for char in letters:
        index = index * 26 + ord(char) - ord("A") + 1
    return index - 1


def _read_xlsx_without_openpyxl(path: Path, sheet_number: int = 1) -> pd.DataFrame:
    ns = {
        "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
        "rel": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    }
    with zipfile.ZipFile(path) as archive:
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall(".//main:si", ns):
                shared_strings.append("".join(text.text or "" for text in item.findall(".//main:t", ns)))

        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        relationship_targets = {
            rel.attrib["Id"]: rel.attrib["Target"]
            for rel in rels
            if rel.tag.endswith("Relationship")
        }
        sheets = workbook.findall(".//main:sheet", ns)
        if not sheets:
            raise ValueError("No worksheets found in workbook.")
        selected = sheets[sheet_number - 1]
        rel_id = selected.attrib[f"{{{ns['rel']}}}id"]
        target = relationship_targets[rel_id]
        sheet_path = "xl/" + target.lstrip("/")
        if not sheet_path.startswith("xl/worksheets/"):
            sheet_path = "xl/worksheets/" + Path(target).name

        sheet = ET.fromstring(archive.read(sheet_path))
        rows: list[list[object]] = []
        for row in sheet.findall(".//main:sheetData/main:row", ns):
            values: list[object] = []
            for cell in row.findall("main:c", ns):
                ref = cell.attrib.get("r", "")
                col_index = _column_index(ref)
                while len(values) < col_index:
                    values.append(np.nan)

                cell_type = cell.attrib.get("t")
                value_node = cell.find("main:v", ns)
                inline_node = cell.find("main:is/main:t", ns)
                value: object = np.nan
                if cell_type == "s" and value_node is not None:
                    value = shared_strings[int(value_node.text or 0)]
                elif cell_type == "inlineStr" and inline_node is not None:
                    value = inline_node.text or ""
                elif value_node is not None:
                    raw = value_node.text or ""
                    try:
                        value = float(raw)
                        if value.is_integer():
                            value = int(value)
                    except ValueError:
                        value = raw
                values.append(value)
            rows.append(values)

    if not rows:
        return pd.DataFrame()
    width = max(len(row) for row in rows)
    padded = [row + [np.nan] * (width - len(row)) for row in rows]
    header = [str(value).strip() for value in padded[0]]
    return pd.DataFrame(padded[1:], columns=header)


def load_table(path: Path) -> pd.DataFrame:
    try:
        return pd.read_excel(path)
    except ImportError:
        return _read_xlsx_without_openpyxl(path)


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.loc[:, [column for column in df.columns if str(column).strip().lower() != "nan"]].copy()

    aliases = {
        "Type_of_dataset": "DatasetType",
        "type_of_dataset": "DatasetType",
        "dataset_type": "DatasetType",
        "grand_true": "TruthTotal",
        "ground_truth": "TruthTotal",
        "Ground_truth": "TruthTotal",
        "Edgren": "Is_Edgren",
    }
    df = df.rename(columns={column: aliases.get(column, column) for column in df.columns})

    required = ["Benchmark", "Dataset", "DatasetType", "Tool"]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    for column in ["TP", "FP", "FN", "TruthTotal", "Recall", "Precision", "F1"]:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")

    if "FN" not in df.columns and {"TruthTotal", "TP"}.issubset(df.columns):
        df["FN"] = df["TruthTotal"] - df["TP"]

    if "Recall" not in df.columns and {"TP", "FN"}.issubset(df.columns):
        df["Recall"] = df["TP"] / (df["TP"] + df["FN"])

    if "Precision" not in df.columns and {"TP", "FP"}.issubset(df.columns):
        df["Precision"] = df["TP"] / (df["TP"] + df["FP"])

    if "F1" not in df.columns and {"Precision", "Recall"}.issubset(df.columns):
        df["F1"] = 2 * df["Precision"] * df["Recall"] / (df["Precision"] + df["Recall"])

    df["DatasetType"] = df["DatasetType"].astype(str).str.lower().str.strip()
    df = df[df["DatasetType"].isin(["real", "simulated"])].copy()
    df["Benchmark"] = df["Benchmark"].astype(str).str.strip()
    df["Tool"] = df["Tool"].astype(str).str.strip()
    df["Dataset"] = df["Dataset"].astype(str).str.strip()

    if "Is_Edgren" in df.columns:
        edgren = df["Is_Edgren"].astype(str).str.lower().str.strip()
        df["Is_Edgren_Binary"] = np.where(edgren.isin(["yes", "true", "1", "part", "partial"]), "yes", "no")

    if "Benchmark_purpose" in df.columns:
        df["Benchmark_purpose"] = df["Benchmark_purpose"].astype(str).str.lower().str.strip()
    elif "new_tool" in df.columns:
        new_tool = df["new_tool"].astype(str).str.lower().str.strip()
        df["Benchmark_purpose"] = np.where(
            new_tool.eq("no") | new_tool.eq("nan") | new_tool.eq(""),
            "independent_comparison",
            "new_tool",
        )
        df["Introduced_tool"] = (
            df["new_tool"]
            .astype(str)
            .str.replace(r"^yes,\s*", "", regex=True)
            .where(df["Benchmark_purpose"].eq("new_tool"), "")
        )

    return df


def fit_dataset_type_model(df: pd.DataFrame, metric: str) -> dict[str, object]:
    model_df = df.dropna(subset=[metric, "DatasetType", "Tool", "Benchmark"]).copy()
    model_df["DatasetType"] = pd.Categorical(model_df["DatasetType"], categories=["simulated", "real"])
    model_df["Tool"] = model_df["Tool"].astype("category")
    model_df["Benchmark"] = model_df["Benchmark"].astype("category")

    if model_df["DatasetType"].nunique() < 2:
        raise ValueError(f"{metric}: model requires both real and simulated observations.")

    result = smf.mixedlm(
        f"{metric} ~ DatasetType + C(Benchmark)",
        model_df,
        groups=model_df["Tool"],
    ).fit(reml=True)

    term = "DatasetType[T.real]"
    conf_int = result.conf_int().loc[term]
    return {
        "metric": metric,
        "n": int(model_df.shape[0]),
        "n_tools": int(model_df["Tool"].nunique()),
        "n_benchmarks": int(model_df["Benchmark"].nunique()),
        "mean_simulated": model_df.loc[model_df["DatasetType"] == "simulated", metric].mean(),
        "mean_real": model_df.loc[model_df["DatasetType"] == "real", metric].mean(),
        "median_simulated": model_df.loc[model_df["DatasetType"] == "simulated", metric].median(),
        "median_real": model_df.loc[model_df["DatasetType"] == "real", metric].median(),
        "beta_real_vs_simulated": result.params[term],
        "ci_low": conf_int.iloc[0],
        "ci_high": conf_int.iloc[1],
        "p_value": result.pvalues[term],
        "group_var": result.cov_re.iloc[0, 0] if result.cov_re.size else np.nan,
        "converged": bool(result.converged),
    }


def leave_one_benchmark_out(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    rows = []
    for benchmark in sorted(df["Benchmark"].dropna().unique(), key=str):
        subset = df[df["Benchmark"] != benchmark]
        try:
            model = fit_dataset_type_model(subset, metric)
            rows.append({"excluded_benchmark": benchmark, **model})
        except Exception as exc:
            rows.append({"excluded_benchmark": benchmark, "metric": metric, "error": str(exc)})
    return pd.DataFrame(rows)


def benchmark_purpose_analysis(df: pd.DataFrame, metric: str) -> pd.DataFrame | None:
    if "Benchmark_purpose" not in df.columns:
        return None
    model_df = df.dropna(subset=[metric, "Benchmark_purpose"]).copy()
    if model_df["Benchmark_purpose"].nunique() < 2:
        return pd.DataFrame({"note": ["Benchmark_purpose has fewer than two categories."]})
    summary = (
        model_df.groupby("Benchmark_purpose", observed=False)[metric]
        .agg(["count", "mean", "median", "std"])
        .reset_index()
    )
    return summary


def fit_benchmark_purpose_model(df: pd.DataFrame, metric: str) -> dict[str, object] | None:
    if "Benchmark_purpose" not in df.columns:
        return None
    model_df = df.dropna(subset=[metric, "Benchmark_purpose", "DatasetType", "Tool"]).copy()
    model_df = model_df[model_df["Benchmark_purpose"].isin(["independent_comparison", "new_tool"])]
    if model_df["Benchmark_purpose"].nunique() < 2:
        return None
    model_df["Benchmark_purpose"] = pd.Categorical(
        model_df["Benchmark_purpose"],
        categories=["independent_comparison", "new_tool"],
    )
    model_df["DatasetType"] = pd.Categorical(model_df["DatasetType"], categories=["simulated", "real"])
    model_df["Tool"] = model_df["Tool"].astype("category")

    result = smf.mixedlm(
        f"{metric} ~ Benchmark_purpose + DatasetType",
        model_df,
        groups=model_df["Tool"],
    ).fit(reml=True)

    term = "Benchmark_purpose[T.new_tool]"
    conf_int = result.conf_int().loc[term]
    return {
        "metric": metric,
        "n": int(model_df.shape[0]),
        "n_tools": int(model_df["Tool"].nunique()),
        "n_independent_observations": int(model_df["Benchmark_purpose"].eq("independent_comparison").sum()),
        "n_new_tool_observations": int(model_df["Benchmark_purpose"].eq("new_tool").sum()),
        "mean_independent": model_df.loc[
            model_df["Benchmark_purpose"] == "independent_comparison", metric
        ].mean(),
        "mean_new_tool": model_df.loc[model_df["Benchmark_purpose"] == "new_tool", metric].mean(),
        "beta_new_tool_vs_independent": result.params[term],
        "ci_low": conf_int.iloc[0],
        "ci_high": conf_int.iloc[1],
        "p_value": result.pvalues[term],
        "converged": bool(result.converged),
    }


def format_p(value: float) -> str:
    if pd.isna(value):
        return "NA"
    if value == 0:
        return "<1e-300"
    exponent = math.floor(math.log10(abs(value))) if value else 0
    if exponent <= -4:
        return f"{value:.2e}"
    return f"{value:.4f}"


def write_report(
    output_path: Path,
    df: pd.DataFrame,
    model_results: list[dict[str, object]],
    edgren_results: list[dict[str, object]],
    missing: list[str],
) -> None:
    lines = []
    lines.append("# Real vs simulated sensitivity analyses")
    lines.append("")
    lines.append(f"Input observations after cleaning: {len(df)}")
    lines.append(f"Benchmarks: {df['Benchmark'].nunique()}")
    lines.append(f"Tools: {df['Tool'].nunique()}")
    lines.append(f"Datasets: {df['Dataset'].nunique()}")
    lines.append("")
    lines.append("## Dataset type models")
    for result in model_results:
        lines.append("")
        lines.append(f"### {result['metric']}")
        lines.append(
            f"Mean simulated = {result['mean_simulated']:.3f}; "
            f"mean real = {result['mean_real']:.3f}; "
            f"beta real vs simulated = {result['beta_real_vs_simulated']:.3f}; "
            f"95% CI {result['ci_low']:.3f} to {result['ci_high']:.3f}; "
            f"p = {format_p(float(result['p_value']))}."
        )
    if edgren_results:
        lines.append("")
        lines.append("## Excluding Edgren-labelled observations")
        for result in edgren_results:
            lines.append(
                f"{result['metric']}: n = {result['n']}; "
                f"beta real vs simulated = {result['beta_real_vs_simulated']:.3f}; "
                f"95% CI {result['ci_low']:.3f} to {result['ci_high']:.3f}; "
                f"p = {format_p(float(result['p_value']))}."
            )
    if missing:
        lines.append("")
        lines.append("## Analyses not run")
        for item in missing:
            lines.append(f"- {item}")
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Analyze real vs simulated fusion benchmark data.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw = load_table(args.input)
    df = normalize_columns(raw)
    df.to_csv(args.output_dir / "cleaned_real_simulated.csv", index=False)

    available_metrics = [metric for metric in ["F1", "Precision", "Recall"] if metric in df.columns]
    missing_analyses: list[str] = []
    model_results: list[dict[str, object]] = []
    for metric in available_metrics:
        result = fit_dataset_type_model(df, metric)
        model_results.append(result)
    pd.DataFrame(model_results).to_csv(args.output_dir / "dataset_type_models.csv", index=False)

    for metric in available_metrics:
        leave_one_benchmark_out(df, metric).to_csv(
            args.output_dir / f"leave_one_benchmark_out_{metric.lower()}.csv",
            index=False,
        )

    edgren_results: list[dict[str, object]] = []
    if "Is_Edgren_Binary" in df.columns:
        no_edgren = df[df["Is_Edgren_Binary"] != "yes"].copy()
        no_edgren.to_csv(args.output_dir / "cleaned_real_simulated_without_edgren.csv", index=False)
        for metric in available_metrics:
            edgren_results.append(fit_dataset_type_model(no_edgren, metric))
        pd.DataFrame(edgren_results).to_csv(args.output_dir / "dataset_type_models_without_edgren.csv", index=False)
    else:
        missing_analyses.append("Sensitivity without Edgren requires an Edgren/Is_Edgren column.")

    purpose = benchmark_purpose_analysis(df, "F1")
    if purpose is None:
        missing_analyses.append("New-tool vs independent comparison requires a Benchmark_purpose column.")
    else:
        purpose.to_csv(args.output_dir / "benchmark_purpose_f1_summary.csv", index=False)
        purpose_models = [
            fit_benchmark_purpose_model(df, metric)
            for metric in available_metrics
        ]
        purpose_models = [result for result in purpose_models if result is not None]
        if purpose_models:
            pd.DataFrame(purpose_models).to_csv(
                args.output_dir / "benchmark_purpose_models.csv",
                index=False,
            )

    write_report(args.output_dir / "summary.md", df, model_results, edgren_results, missing_analyses)
    print(f"Saved results to: {args.output_dir}")
    print((args.output_dir / "summary.md").read_text(encoding="utf-8"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
