#!/usr/bin/env python3
"""Build Supplementary Table 6 from Comment 14 whole-brain model results."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MODELS = PACKAGE_ROOT / "results/model_results/whole_brain_qmri_models.csv"
DEFAULT_OUTPUT = PACKAGE_ROOT / "results/tables/Supplementary_Table_6.csv"
METRICS = ["T1", "MTR", "FA", "MD"]
SOURCE_OUTCOMES = [
    "EDSS",
    "MSPro",
    "T25FW",
    "9HPT-D",
    "9HPT-ND",
]
TABLE_OUTCOMES = SOURCE_OUTCOMES
EXPECTED_N = {
    "EDSS": {"T1": 80, "MTR": 80, "FA": 79, "MD": 79},
    "MSPro": {"T1": 79, "MTR": 79, "FA": 78, "MD": 78},
    "T25FW": {"T1": 76, "MTR": 76, "FA": 75, "MD": 75},
    "9HPT-D": {"T1": 79, "MTR": 79, "FA": 78, "MD": 78},
    "9HPT-ND": {"T1": 79, "MTR": 79, "FA": 78, "MD": 78},
}
EXPECTED_POSITIVE = {"EDSS": 31, "MSPro": 39}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models", type=Path, default=DEFAULT_MODELS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def inside_package(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    try:
        resolved.relative_to(PACKAGE_ROOT)
    except ValueError as exc:
        raise ValueError(f"{label} must remain inside revisionExtras: {resolved}") from exc
    return resolved


def format_p(value: float) -> str:
    return "< 0.0001" if value < 0.0001 else f"{value:.4f}"


def main() -> None:
    args = parse_args()
    model_path = args.models.expanduser().resolve()
    output_path = inside_package(args.output, "output")
    with model_path.open(newline="", encoding="utf-8") as handle:
        source_rows = list(csv.DictReader(handle))

    if len(source_rows) != 20:
        raise AssertionError("Expected 20 imaging-only whole-brain models")
    if any(row["Adjusted"] != "0" for row in source_rows):
        raise AssertionError("Supplementary Table 6 must use unadjusted models")
    if any(row["Demographics"] != "none" for row in source_rows):
        raise AssertionError("Supplementary Table 6 must contain no covariates")
    if any(row["AnalysisScope"] != "matched_manuscript_case" for row in source_rows):
        raise AssertionError("Whole-brain models must use manuscript-matched cases")
    if {row["Metric"] for row in source_rows} != set(METRICS):
        raise AssertionError("Unexpected metric inventory")
    if {row["Outcome"] for row in source_rows} != set(SOURCE_OUTCOMES):
        raise AssertionError("Unexpected source outcome inventory")
    if any(row["Status"] != "COMPLETED" for row in source_rows):
        raise AssertionError("At least one whole-brain model did not complete")

    selected = [row for row in source_rows if row["Outcome"] in TABLE_OUTCOMES]
    if len(selected) != 20:
        raise AssertionError("Expected 20 imaging-only models for the five table outcomes")
    by_key = {(row["Outcome"], row["Metric"]): row for row in selected}
    if len(by_key) != 20:
        raise AssertionError("Selected outcome/metric keys are not unique")

    table_rows: list[dict[str, str]] = []
    for outcome in TABLE_OUTCOMES:
        for metric in METRICS:
            row = by_key[(outcome, metric)]
            n = int(float(row["N"]))
            if n != EXPECTED_N[outcome][metric]:
                raise AssertionError(f"Cohort N mismatch for {outcome}/{metric}")
            positive = float(row["PositiveN"])
            performance = float(row["Performance"])
            ci_lower = float(row["CI_lower"])
            ci_upper = float(row["CI_upper"])
            if outcome in {"EDSS", "MSPro"}:
                if not all(math.isfinite(value) for value in (positive, ci_lower, ci_upper)):
                    raise AssertionError(f"Missing binary result for {outcome}/{metric}")
                if int(positive) != EXPECTED_POSITIVE[outcome]:
                    raise AssertionError(f"Positive N mismatch for {outcome}/{metric}")
                n_text = f"{n} ({int(positive)})"
                performance_text = (
                    f"AUC {performance:.3f} ({ci_lower:.3f}-{ci_upper:.3f})"
                )
            else:
                n_text = str(n)
                performance_text = f"R2 {performance:.3f}"
            table_rows.append(
                {
                    "Outcome": outcome,
                    "MRI Metric": metric,
                    "Model": "Whole brain",
                    "N": n_text,
                    "Apparent performance": performance_text,
                    "AIC": f"{float(row['AIC']):.2f}",
                    "p (raw)": format_p(float(row["p_raw"])),
                }
            )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(table_rows[0]),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(table_rows)
    print(f"Wrote {len(table_rows)} rows to {output_path}")


if __name__ == "__main__":
    main()
