#!/usr/bin/env python3
"""Build aggregate lesion, exclusion, and exploratory group-test tables.

The one row-per-subject export is deliberately written outside the repository.
All outputs are created exclusively; existing files are never overwritten.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu


METRICS = ("T1", "MTR", "FA", "MD", "NDI", "ODI", "FWF")
FIGURE5_METRICS = METRICS[:4]
TRACT_REGIONS = {
    "Association": "Association",
    "Cerebellar": "Cerebellar",
    "Occipitoparietal": "Occipitoparietal",
    "ProjectionBrainstem": "PB",
}
CLASSICAL_REGIONS = (
    "periventricular", "juxtacortical", "infratentorial", "deepwhitematter"
)
LESION_METRICS = ("LN", "LV", "Lnorm")
MS_GROUPS = ("RRMS", "SPMS", "PPMS")
ALL_GROUPS = MS_GROUPS + ("Control",)
REPO_ROOT = Path(__file__).resolve().parents[3]


def read_table(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    data = pd.read_excel(path) if path.suffix.lower() == ".xlsx" else pd.read_csv(path)
    if "SubjectID" not in data.columns or data.SubjectID.isna().any():
        raise ValueError(f"Missing SubjectID values: {path}")
    if data.SubjectID.duplicated().any():
        raise ValueError(f"Duplicate SubjectID values: {path}")
    return data.set_index("SubjectID", verify_integrity=True)


def subtype(subject: str) -> str:
    for prefix, label in (("sub-RR", "RRMS"), ("sub-SP", "SPMS"),
                          ("sub-PP", "PPMS"), ("sub-C", "Control")):
        if subject.startswith(prefix):
            return label
    raise ValueError(f"Unknown subtype prefix: {subject}")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_csv_exclusive(data: pd.DataFrame, path: Path, *, private: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if private:
        descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            data.to_csv(handle, index=False)
    else:
        with path.open("x", encoding="utf-8", newline="") as handle:
            data.to_csv(handle, index=False)


def holm(values: list[float]) -> list[float]:
    n = len(values)
    order = sorted(range(n), key=lambda i: values[i])
    result = [math.nan] * n
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, min(1.0, (n - rank) * values[index]))
        result[index] = running
    return result


def bh(values: list[float]) -> list[float]:
    n = len(values)
    order = sorted(range(n), key=lambda i: values[i], reverse=True)
    result = [math.nan] * n
    running = 1.0
    for rank, index in enumerate(order):
        running = min(running, values[index] * n / (n - rank))
        result[index] = running
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical", type=Path, required=True)
    parser.add_argument("--metrics-dir", type=Path, required=True)
    parser.add_argument("--motor-dir", type=Path, required=True)
    parser.add_argument("--tract-lesion", type=Path, required=True)
    parser.add_argument("--classical-lesion", type=Path, required=True)
    parser.add_argument("--public-output-dir", type=Path, required=True)
    parser.add_argument("--private-output", type=Path, required=True)
    args = parser.parse_args()

    private_output = args.private_output.resolve()
    if private_output.is_relative_to(REPO_ROOT):
        raise ValueError("Per-subject export must be outside the Git repository")
    outputs = [
        args.public_output_dir / "lesion_load_by_region_subtype.csv",
        args.public_output_dir / "classical_exclusions_by_metric.csv",
        args.public_output_dir / "posthoc_group_tests.csv",
        args.public_output_dir / "analysis_report.json",
        private_output,
    ]
    existing = [str(path) for path in outputs if path.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite: " + ", ".join(existing))

    clinical = read_table(args.clinical)
    tract_lesion = read_table(args.tract_lesion)
    classical_lesion = read_table(args.classical_lesion)
    classical = read_table(args.metrics_dir / "ClassicalRegionAllMetrics.csv")
    if len(clinical) != 132 or len(tract_lesion) != 90 or len(classical_lesion) != 90 or len(classical) != 87:
        raise AssertionError("Frozen source cohort counts changed")
    if set(tract_lesion.index) != set(classical_lesion.index):
        raise AssertionError("Figure 4 lesion tables use different participants")
    ms_ids = [s for s in clinical.index if subtype(s) != "Control"]
    if len(ms_ids) != 89:
        raise AssertionError("Expected 89 clinical MS patients")

    standard = {}
    motor = {}
    for metric in METRICS:
        standard[metric] = read_table(args.metrics_dir / f"GroupTract{metric}_All.csv")
        motor[metric] = read_table(args.motor_dir / f"Motor{metric}_All.csv")
        if len(standard[metric]) != 132 or set(standard[metric].index) != set(motor[metric].index):
            raise AssertionError(f"Standard/Motor source mismatch for {metric}")

    # Figure 4 follows the original radar scripts: all 90 lesion-load rows,
    # true zeros retained, mean and SEM by subtype and anatomical region.
    lesion_rows = []
    figure4_features = []
    for framework, frame, regions in (
        ("Tract", tract_lesion, TRACT_REGIONS),
        ("Classical", classical_lesion, {r: r for r in CLASSICAL_REGIONS}),
    ):
        for region, source_prefix in regions.items():
            for metric in LESION_METRICS:
                column = source_prefix + metric
                if column not in frame.columns:
                    raise KeyError(column)
                values = pd.to_numeric(frame[column], errors="raise")
                figure4_features.append((framework, region, metric, values))
                for group in MS_GROUPS:
                    sample = values[[subtype(s) == group for s in frame.index]].dropna()
                    n = len(sample)
                    lesion_rows.append({
                        "Framework": framework, "Region": region,
                        "Subtype": group, "Metric": metric, "N": n,
                        "Mean": sample.mean(), "SD": sample.std(ddof=1),
                        "SEM": sample.std(ddof=1) / math.sqrt(n) if n else math.nan,
                        "Median": sample.median(), "Q1": sample.quantile(.25),
                        "Q3": sample.quantile(.75), "ZeroCount": int((sample == 0).sum()),
                    })
    lesion_summary = pd.DataFrame(lesion_rows)
    if len(lesion_summary) != 72:
        raise AssertionError("Expected 72 Figure 4 subtype summaries")

    # Structural absence is defined from the four regional lesion counts, not
    # inferred from missing qMRI means (which also have other causes).
    region_counts = classical_lesion.reindex(ms_ids)[[r + "LN" for r in CLASSICAL_REGIONS]]
    if region_counts.isna().any().any():
        raise AssertionError("Classical regional lesion counts are missing")
    lacking_lesion = (region_counts.isna() | (region_counts <= 0)).any(axis=1)
    lacking_lesions_in_all_four = (region_counts <= 0).all(axis=1)
    if (lacking_lesions_in_all_four & ~lacking_lesion).any():
        raise AssertionError("All-four lesion absence must be within any-region absence")
    exclusion_rows = []
    for metric in METRICS:
        predictors = classical.reindex(ms_ids)[[r + metric for r in CLASSICAL_REGIONS]]
        complete = (np.isfinite(predictors) & (predictors != 0)).all(axis=1)
        exclusion_rows.append({
            "Metric": "ISOVF" if metric == "FWF" else metric,
            "SourceMetric": metric,
            "ClinicalMSN": len(ms_ids),
            "ExcludedLackingLesionInAnyClassicalRegionN": int(lacking_lesion.sum()),
            "ExcludedLackingLesionsInAllFourClassicalRegionsN": int(lacking_lesions_in_all_four.sum()),
            "ExcludedOtherClassicalPredictorReasonN": int((~complete & ~lacking_lesion).sum()),
            "TotalUnavailableClassicalPredictorN": int((~complete).sum()),
            "CompleteClassicalPredictorN": int(complete.sum()),
        })
    exclusions = pd.DataFrame(exclusion_rows)
    if exclusions.ClinicalMSN.nunique() != 1 or exclusions.ClinicalMSN.iloc[0] != 89:
        raise AssertionError("Unexpected denominator in exclusion table")

    # Pairwise exploratory post-hoc tests use a two-sided Mann-Whitney U test.
    # Holm adjusts the three or six comparisons for each feature; BH adjusts
    # across all pairwise tests for each complete figure and across both figures.
    test_rows = []
    feature_index = 0

    def add_tests(figure: str, framework: str, tissue: str, metric: str,
                  region: str, values: pd.Series, groups: tuple[str, ...],
                  zero_is_missing: bool) -> None:
        nonlocal feature_index
        feature_index += 1
        samples = {}
        for group in groups:
            member = [subtype(s) == group for s in values.index]
            selected = pd.to_numeric(values[member], errors="raise").dropna()
            if zero_is_missing:
                selected = selected[selected != 0]
            samples[group] = selected.to_numpy(dtype=float)
        if any(len(v) < 2 for v in samples.values()):
            raise ValueError(f"Insufficient group observations: {figure} {framework} {metric} {region}")
        h, p_omnibus = kruskal(*samples.values(), nan_policy="raise")
        local = []
        for group1, group2 in combinations(groups, 2):
            x, y = samples[group1], samples[group2]
            u, p_pair = mannwhitneyu(x, y, alternative="two-sided", method="asymptotic")
            local.append({
                "Figure": figure, "Framework": framework, "Tissue": tissue,
                "Metric": metric, "Region": region,
                "Group1": group1, "Group2": group2,
                "N1": len(x), "N2": len(y),
                "Median1": float(np.median(x)), "Median2": float(np.median(y)),
                "U": float(u), "P_raw": float(p_pair),
                "KW_H": float(h), "KW_P_raw": float(p_omnibus),
                "FeatureIndex": feature_index,
                "ZeroPolicy": "exclude" if zero_is_missing else "retain",
            })
        adjusted = holm([row["P_raw"] for row in local])
        for row, p_adjusted in zip(local, adjusted):
            row["P_HolmWithinFeature"] = p_adjusted
            test_rows.append(row)

    for framework, region, metric, values in figure4_features:
        add_tests("Figure4", framework, "Lesion load", metric, region,
                  values, MS_GROUPS, False)
    for metric in FIGURE5_METRICS:
        frame = standard[metric]
        for region in TRACT_REGIONS:
            for tissue, groups in (("Lesion", MS_GROUPS), ("NAWM", ALL_GROUPS)):
                add_tests("Figure5", "Tract", tissue, metric, region,
                          frame[region + "_" + tissue], groups, True)
        for region in CLASSICAL_REGIONS:
            add_tests("Figure5", "Classical", "Lesion", metric, region,
                      classical[region + metric], MS_GROUPS, True)
    tests = pd.DataFrame(test_rows)
    if len(tests) != 264:
        raise AssertionError(f"Expected 264 pairwise post-hoc rows, got {len(tests)}")
    tests["P_BHAllPairwise"] = bh(tests.P_raw.tolist())
    tests["P_BHWithinFigure"] = np.nan
    tests["KW_P_BHAllFeatures"] = np.nan
    for figure in ("Figure4", "Figure5"):
        rows = tests.index[tests.Figure == figure].tolist()
        tests.loc[rows, "P_BHWithinFigure"] = bh(tests.loc[rows, "P_raw"].tolist())
    first = tests.drop_duplicates("FeatureIndex")
    correction = dict(zip(first.FeatureIndex, bh(first.KW_P_raw.tolist())))
    tests["KW_P_BHAllFeatures"] = tests.FeatureIndex.map(correction)
    tests = tests.drop(columns="FeatureIndex")

    # One wide private CSV holds all raw per-subject model predictors. It also
    # includes the combined lesion/NAWM values plotted in Figure 5 and both
    # lesion-load frameworks needed to reproduce Figure 4 group tests.
    columns = {}
    for metric in METRICS:
        frame = standard[metric]
        for region in TRACT_REGIONS:
            for side in ("L", "R"):
                for tissue in ("Tail", "NAWM"):
                    column = f"{region}{side}_{tissue}"
                    columns[f"{metric}__TractModel__{column}"] = frame[column]
            for tissue in ("Lesion", "NAWM"):
                column = f"{region}_{tissue}"
                columns[f"{metric}__TractFigure__{column}"] = frame[column]
        for column in ("MotorL_Tail", "MotorR_Tail", "MotorL_NAWM", "MotorR_NAWM"):
            columns[f"{metric}__Motor__{column}"] = motor[metric][column]
        for region in CLASSICAL_REGIONS:
            column = region + metric
            columns[f"{metric}__ClassicalMean__{column}"] = classical[column]
    for region, prefix in TRACT_REGIONS.items():
        for metric in LESION_METRICS:
            columns[f"LesionLoad__Tract__{region}__{metric}"] = tract_lesion[prefix + metric]
    for region in CLASSICAL_REGIONS:
        for metric in LESION_METRICS:
            columns[f"LesionLoad__Classical__{region}__{metric}"] = classical_lesion[region + metric]
    private = pd.concat(columns, axis=1, sort=True)
    private.insert(0, "Subtype", [subtype(s) for s in private.index])
    private.index.name = "SubjectID"
    private = private.reset_index().sort_values("SubjectID", kind="stable")
    if private.SubjectID.duplicated().any() or len(private) < 132:
        raise AssertionError("Per-subject export has missing or duplicate rows")

    write_csv_exclusive(lesion_summary, outputs[0])
    write_csv_exclusive(exclusions, outputs[1])
    write_csv_exclusive(tests, outputs[2])
    write_csv_exclusive(private, outputs[4], private=True)
    report = {
        "status": "PASS",
        "figure4_subtype_summary_rows": len(lesion_summary),
        "figure4_source_ms_subjects": len(tract_lesion),
        "clinical_ms_subjects": len(ms_ids),
        "classical_no_lesion_subjects": int(lacking_lesion.sum()),
        "classical_no_lesions_in_all_four_subjects": int(lacking_lesions_in_all_four.sum()),
        "figure4_posthoc_rows": int((tests.Figure == "Figure4").sum()),
        "figure5_posthoc_rows": int((tests.Figure == "Figure5").sum()),
        "private_predictor_rows": len(private),
        "private_predictor_columns": len(private.columns),
        "public_input_sha256": {
            "clinical": sha256(args.clinical),
            "tract_lesion": sha256(args.tract_lesion),
            "classical_lesion": sha256(args.classical_lesion),
            "classical_qmri": sha256(args.metrics_dir / "ClassicalRegionAllMetrics.csv"),
        },
    }
    with outputs[3].open("x", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
