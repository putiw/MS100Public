#!/usr/bin/env python3
"""Add Figure 5 workbook values to the validated private predictor CSV.

The existing private CSV already holds model predictors, corrected qMRI
figure values, and Figure 4 lesion loads. This adds the 32 combined-hemisphere
values from the local figure-named T1/MTR/FA/MD workbooks. Never writes inside
the public repository or overwrites an existing file.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd


METRICS = ("T1", "MTR", "FA", "MD")
REGIONS = ("Association", "Cerebellar", "Occipitoparietal", "ProjectionBrainstem")
REPO_ROOT = Path(__file__).resolve().parents[3]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original-private", type=Path, required=True)
    parser.add_argument("--figure5-workbook-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.is_relative_to(REPO_ROOT):
        raise ValueError("Per-subject CSV must remain outside the public repository")
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    data = pd.read_csv(args.original_private)
    if len(data) != 133 or data.SubjectID.duplicated().any() or "Subtype" not in data:
        raise AssertionError("Original private predictor cohort changed")
    data = data.set_index("SubjectID", verify_integrity=True)
    expected_ids = None
    for metric in METRICS:
        path = args.figure5_workbook_dir / f"GroupTract{metric}.xlsx"
        figure = pd.read_excel(path).set_index("SubjectID", verify_integrity=True)
        if len(figure) != 132 or not set(figure.index).issubset(data.index):
            raise AssertionError(f"Unexpected Figure 5 cohort: {path}")
        if expected_ids is None:
            expected_ids = set(figure.index)
        elif set(figure.index) != expected_ids:
            raise AssertionError("Figure 5 workbooks have different subject sets")
        for region in REGIONS:
            for tissue in ("Lesion", "NAWM"):
                field = f"{region}_{tissue}"
                column = f"{metric}__Figure5Workbook__{field}"
                if column in data:
                    raise AssertionError(f"New column already exists: {column}")
                data[column] = figure[field].reindex(data.index)
    data = data.reset_index()
    if len(data.columns) != 282 or len(data) != 133:
        raise AssertionError("Expected 133 subjects and 282 fields")
    output.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
        data.to_csv(handle, index=False)
    print(f"Wrote {len(data)} private subject rows and {len(data.columns)} fields to {output}")


if __name__ == "__main__":
    main()
