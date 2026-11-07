#!/usr/bin/env python3
"""Recompute figure-source summaries and subtype tests from mounted workbooks.

This script never modifies source workbooks or existing results. Figure 4 has
two explicit cohort policies because the mounted classical lesion-load file
contains one additional RRMS participant relative to the tract file and the
manuscript's n=49 RRMS caption. Figure 4A bars use the 89 tract participants
under either policy; only Figure 4B's classical radar differs. Figure 5
includes all three captioned rows: classical lesions, tract lesions, and
tract NAWM.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import os
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests


METRICS = ("T1", "MTR", "FA", "MD")
LESION_METRICS = ("LN", "LV", "Lnorm")
TRACT_REGIONS = {
    "Association": "Association",
    "Cerebellar": "Cerebellar",
    "Occipitoparietal": "Occipitoparietal",
    "ProjectionBrainstem": "PB",
}
CLASSICAL_REGIONS = (
    "periventricular", "juxtacortical", "infratentorial", "deepwhitematter"
)
MS_GROUPS = ("RRMS", "SPMS", "PPMS")
ALL_GROUPS = MS_GROUPS + ("Control",)
PREFIX = {"RRMS": "sub-RR", "SPMS": "sub-SP", "PPMS": "sub-PP", "Control": "sub-C"}
REPO_ROOT = Path(__file__).resolve().parents[3]


def read_unique(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    frame = pd.read_excel(path) if path.suffix.lower() == ".xlsx" else pd.read_csv(path)
    if "SubjectID" not in frame or frame.SubjectID.isna().any() or frame.SubjectID.duplicated().any():
        raise ValueError(f"Missing or duplicate SubjectID in {path}")
    result = frame.set_index("SubjectID", verify_integrity=True)
    if not result.index.map(lambda value: isinstance(value, str)).all():
        raise ValueError(f"Non-text SubjectID in {path}")
    return result


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sample(values: pd.Series, group: str, *, exclude_zero: bool) -> np.ndarray:
    selected = pd.to_numeric(
        values[values.index.str.startswith(PREFIX[group])], errors="raise"
    ).dropna()
    if exclude_zero:
        selected = selected[selected != 0]
    return selected.to_numpy(dtype=float)


def source_counts(frame: pd.DataFrame) -> tuple[int, int, int, int]:
    return tuple(int(frame.index.str.startswith(PREFIX[group]).sum()) for group in ALL_GROUPS)


def add_feature(
    tests: list[dict],
    summaries: list[dict],
    *,
    figure: str,
    panel: str,
    policy: str,
    framework: str,
    tissue: str,
    metric: str,
    region: str,
    values: pd.Series,
    groups: tuple[str, ...],
    exclude_zero: bool,
    source_hash: str,
) -> None:
    cohorts = {group: sample(values, group, exclude_zero=exclude_zero) for group in groups}
    if any(len(value) < 2 or not np.isfinite(value).all() for value in cohorts.values()):
        raise ValueError(f"Invalid samples for {figure} {panel} {metric} {region}")
    feature = f"{figure}|{panel}|{framework}|{tissue}|{metric}|{region}"
    h, omnibus_p = kruskal(*cohorts.values(), nan_policy="raise")
    for group, observed in cohorts.items():
        summaries.append({
            "Figure": figure,
            "Panel": panel,
            "CohortPolicy": policy,
            "Framework": framework,
            "Tissue": tissue,
            "Metric": metric,
            "Region": region,
            "Subtype": group,
            "N": len(observed),
            "Mean": float(np.mean(observed)),
            "SD": float(np.std(observed, ddof=1)),
            "SEM": float(np.std(observed, ddof=1) / math.sqrt(len(observed))),
            "Median": float(np.median(observed)),
            "Q1": float(np.quantile(observed, .25)),
            "Q3": float(np.quantile(observed, .75)),
            "ZeroCount": int((observed == 0).sum()),
            "ZeroPolicy": "exclude" if exclude_zero else "retain",
            "SourceSHA256": source_hash,
        })
    for first, second in combinations(groups, 2):
        x, y = cohorts[first], cohorts[second]
        u, p = mannwhitneyu(x, y, alternative="two-sided", method="asymptotic")
        tests.append({
            "Figure": figure,
            "Panel": panel,
            "CohortPolicy": policy,
            "Framework": framework,
            "Tissue": tissue,
            "Metric": metric,
            "Region": region,
            "Group1": first,
            "Group2": second,
            "N1": len(x),
            "N2": len(y),
            "Median1": float(np.median(x)),
            "Median2": float(np.median(y)),
            "U": float(u),
            "RankBiserial": float(2 * u / (len(x) * len(y)) - 1),
            "P_raw": float(p),
            "KW_H": float(h),
            "KW_P_raw": float(omnibus_p),
            "FeatureKey": feature,
            "ZeroPolicy": "exclude" if exclude_zero else "retain",
            "SourceSHA256": source_hash,
        })


def adjust_tests(rows: list[dict], policy: str) -> pd.DataFrame:
    result = pd.DataFrame(rows)
    if len(result) != 270 or (result.Figure == "Figure4").sum() != 78 or (result.Figure == "Figure5").sum() != 192:
        raise AssertionError(f"Unexpected test count under {policy}")
    if result.FeatureKey.nunique() != 74:
        raise AssertionError("Expected 74 distinct plotted features")
    result["P_HolmWithinFeature"] = np.nan
    result["P_BHWithinFigure"] = np.nan
    result["P_BHAllFigurePairs"] = multipletests(result.P_raw, method="fdr_bh")[1]
    for _, index in result.groupby("FeatureKey", sort=False).groups.items():
        result.loc[index, "P_HolmWithinFeature"] = multipletests(
            result.loc[index, "P_raw"], method="holm"
        )[1]
    for _, index in result.groupby("Figure", sort=False).groups.items():
        result.loc[index, "P_BHWithinFigure"] = multipletests(
            result.loc[index, "P_raw"], method="fdr_bh"
        )[1]
    omnibus = result.drop_duplicates("FeatureKey")
    q_values = multipletests(omnibus.KW_P_raw, method="fdr_bh")[1]
    result["KW_P_BHAllFeatures"] = result.FeatureKey.map(
        dict(zip(omnibus.FeatureKey, q_values))
    )
    return result.drop(columns="FeatureKey")


def write_csv_exclusive(frame: pd.DataFrame, path: Path, *, private: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if private:
        descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            frame.to_csv(handle, index=False)
    else:
        with path.open("x", encoding="utf-8", newline="") as handle:
            frame.to_csv(handle, index=False)


def build_private(
    prior_path: Path,
    tract_load: pd.DataFrame,
    classical_load: pd.DataFrame,
    tract_qmri: dict[str, pd.DataFrame],
    classical_qmri: pd.DataFrame,
) -> pd.DataFrame:
    previous = pd.read_csv(prior_path).set_index("SubjectID", verify_integrity=True)
    if len(previous) != 133 or len(previous.columns) != 281 or "Subtype" not in previous:
        raise AssertionError("Validated private source has unexpected shape")
    result = previous.copy()
    for region, prefix in TRACT_REGIONS.items():
        for metric in LESION_METRICS:
            column = f"LesionLoad__Tract__{region}__{metric}"
            mounted = tract_load[prefix + metric].reindex(result.index)
            shared = mounted.loc[tract_load.index]
            earlier = previous.loc[tract_load.index, column]
            if not np.allclose(shared, earlier, rtol=0, atol=1e-12, equal_nan=True):
                raise AssertionError(f"Tract lesion-load values changed for shared subjects: {column}")
            result[column] = mounted
    for region in CLASSICAL_REGIONS:
        for metric in LESION_METRICS:
            column = f"LesionLoad__Classical__{region}__{metric}"
            mounted = classical_load[region + metric].reindex(result.index)
            if not np.allclose(mounted, previous[column], rtol=0, atol=1e-12, equal_nan=True):
                raise AssertionError(f"Classical lesion-load values changed: {column}")
            result[column] = mounted
    for metric, frame in tract_qmri.items():
        for region in TRACT_REGIONS:
            for tissue in ("Lesion", "NAWM"):
                column = f"{metric}__Figure5Workbook__{region}_{tissue}"
                mounted = frame[f"{region}_{tissue}"].reindex(result.index)
                if not np.allclose(mounted, previous[column], rtol=0, atol=1e-12, equal_nan=True):
                    raise AssertionError(f"Figure 5 workbook values changed: {column}")
    for metric in METRICS:
        for region in CLASSICAL_REGIONS:
            column = f"{metric}__ClassicalMean__{region}{metric}"
            mounted = classical_qmri[region + metric].reindex(result.index)
            if not np.allclose(mounted, previous[column], rtol=0, atol=1e-12, equal_nan=True):
                raise AssertionError(f"Classical qMRI values changed: {column}")
    result["Figure4TractSourceMember"] = result.index.isin(tract_load.index).astype(int)
    result["Figure4ClassicalSourceMember"] = result.index.isin(classical_load.index).astype(int)
    result = result.reset_index()
    if result.shape != (133, 284) or result.SubjectID.duplicated().any():
        raise AssertionError("Private output shape or uniqueness changed")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stats-dir", type=Path, required=True)
    parser.add_argument("--validated-private", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--private-output", type=Path, required=True)
    args = parser.parse_args()
    if args.private_output.resolve().is_relative_to(REPO_ROOT):
        raise ValueError("Per-subject output must remain outside the public repository")
    outputs = {
        "tests_scripted": args.output_dir / "posthoc_fig4_fig5_mounted_as_scripted.csv",
        "tests_matched": args.output_dir / "posthoc_fig4_fig5_caption_matched_89.csv",
        "bar_scripted": args.output_dir / "figure4a_bar_by_subtype_mounted_as_scripted.csv",
        "bar_matched": args.output_dir / "figure4a_bar_by_subtype_caption_matched_89.csv",
        "load_scripted": args.output_dir / "lesion_load_by_region_subtype_mounted_as_scripted.csv",
        "load_matched": args.output_dir / "lesion_load_by_region_subtype_caption_matched_89.csv",
        "figure5": args.output_dir / "figure5_plotted_values_mounted.csv",
        "private": args.private_output,
    }
    existing = [str(path) for path in outputs.values() if path.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite: " + ", ".join(existing))

    stats = args.stats_dir
    tract_path = stats / "GroupTractLesionLoad.xlsx"
    classical_load_path = stats / "icometrixLesionLoad.xlsx"
    classical_qmri_path = stats / "icometrixLesionStats.xlsx"
    tract_load = read_unique(tract_path)
    classical_load = read_unique(classical_load_path)
    classical_qmri = read_unique(classical_qmri_path)
    tract_qmri = {metric: read_unique(stats / f"GroupTract{metric}.xlsx") for metric in METRICS}
    if len(tract_load) != 89 or source_counts(tract_load) != (49, 23, 17, 0):
        raise AssertionError("Mounted Figure 4 tract cohort changed")
    if len(classical_load) != 90 or source_counts(classical_load) != (50, 23, 17, 0):
        raise AssertionError("Mounted Figure 4 classical cohort changed")
    if not set(tract_load.index).issubset(classical_load.index):
        raise AssertionError("Figure 4 tract cohort is not inside classical cohort")
    if len(classical_qmri) != 87 or source_counts(classical_qmri) != (47, 23, 17, 0):
        raise AssertionError("Mounted Figure 5 classical cohort changed")
    figure5_ids = set(tract_qmri[METRICS[0]].index)
    if any(len(frame) != 132 or set(frame.index) != figure5_ids or source_counts(frame) != (49, 23, 17, 43)
           for frame in tract_qmri.values()):
        raise AssertionError("Mounted Figure 5 tract workbooks have inconsistent cohorts")
    if not set(classical_qmri.index).issubset(figure5_ids):
        raise AssertionError("Classical qMRI cohort is not inside Figure 5 cohort")

    path_hash = {path: sha256(path) for path in (
        tract_path, classical_load_path, classical_qmri_path,
        *(stats / f"GroupTract{metric}.xlsx" for metric in METRICS),
    )}
    tests_by_policy: dict[str, pd.DataFrame] = {}
    bar_by_policy: dict[str, pd.DataFrame] = {}
    load_by_policy: dict[str, pd.DataFrame] = {}
    figure5_summary: pd.DataFrame | None = None
    for policy in ("mounted_as_scripted", "caption_matched_89"):
        rows: list[dict] = []
        summaries: list[dict] = []
        classical_frame = (
            classical_load if policy == "mounted_as_scripted"
            else classical_load.loc[tract_load.index]
        )
        whole_brain_from_regions = classical_load.loc[tract_load.index, [
            region + "LV" for region in CLASSICAL_REGIONS
        ]].sum(axis=1)
        if not np.allclose(
            whole_brain_from_regions.to_numpy(), tract_load.WBLV.to_numpy(),
            rtol=0, atol=1e-8,
        ):
            raise AssertionError("Whole-brain LV does not match classical regional sum")
        bar_summaries: list[dict] = []
        for framework, region, values in (
            ("WholeBrain", "WholeBrain", tract_load.WBLV / 1000),
            ("Tract", "FourTractRegionSum", tract_load[[
                prefix + "LV" for prefix in TRACT_REGIONS.values()
            ]].sum(axis=1) / 1000),
        ):
            add_feature(
                rows, bar_summaries, figure="Figure4", panel="A bar",
                policy=policy, framework=framework, tissue="Lesion volume (mL)",
                metric="LV", region=region, values=values,
                groups=MS_GROUPS, exclude_zero=False, source_hash=path_hash[tract_path],
            )
        if len(bar_summaries) != 6:
            raise AssertionError("Expected six Figure 4A bar summaries")
        bar_by_policy[policy] = pd.DataFrame(bar_summaries)
        for framework, frame, regions, path in (
            ("Tract", tract_load, TRACT_REGIONS, tract_path),
            ("Classical", classical_frame, {region: region for region in CLASSICAL_REGIONS}, classical_load_path),
        ):
            for region, prefix in regions.items():
                for metric in LESION_METRICS:
                    add_feature(
                        rows, summaries, figure="Figure4", panel="B radar",
                        policy=policy, framework=framework, tissue="Lesion load",
                        metric=metric, region=region, values=frame[prefix + metric],
                        groups=MS_GROUPS, exclude_zero=False, source_hash=path_hash[path],
                    )
        figure4_summary = pd.DataFrame(summaries)
        if len(figure4_summary) != 72:
            raise AssertionError("Expected 72 Figure 4 plotted-value summaries")
        load_by_policy[policy] = figure4_summary

        figure5_summaries: list[dict] = []
        for metric in METRICS:
            frame = tract_qmri[metric]
            path = stats / f"GroupTract{metric}.xlsx"
            for region in CLASSICAL_REGIONS:
                add_feature(
                    rows, figure5_summaries, figure="Figure5",
                    panel="top classical lesion", policy=policy,
                    framework="Classical", tissue="Lesion", metric=metric,
                    region=region, values=classical_qmri[region + metric],
                    groups=MS_GROUPS, exclude_zero=True,
                    source_hash=path_hash[classical_qmri_path],
                )
            for region in TRACT_REGIONS:
                for tissue, groups, panel in (
                    ("Lesion", MS_GROUPS, "middle tract lesion"),
                    ("NAWM", ALL_GROUPS, "bottom tract NAWM"),
                ):
                    add_feature(
                        rows, figure5_summaries, figure="Figure5",
                        panel=panel, policy=policy, framework="Tract",
                        tissue=tissue, metric=metric, region=region,
                        values=frame[f"{region}_{tissue}"], groups=groups,
                        exclude_zero=True, source_hash=path_hash[path],
                    )
        if len(figure5_summaries) != 160:
            raise AssertionError("Expected 160 Figure 5 plotted-value summaries")
        current_figure5 = pd.DataFrame(figure5_summaries)
        if figure5_summary is None:
            figure5_summary = current_figure5
        else:
            independent = current_figure5.drop(columns="CohortPolicy")
            first = figure5_summary.drop(columns="CohortPolicy")
            pd.testing.assert_frame_equal(independent, first)
        tests_by_policy[policy] = adjust_tests(rows, policy)

    private = build_private(
        args.validated_private, tract_load, classical_load, tract_qmri, classical_qmri
    )
    write_csv_exclusive(tests_by_policy["mounted_as_scripted"], outputs["tests_scripted"])
    write_csv_exclusive(tests_by_policy["caption_matched_89"], outputs["tests_matched"])
    write_csv_exclusive(bar_by_policy["mounted_as_scripted"], outputs["bar_scripted"])
    write_csv_exclusive(bar_by_policy["caption_matched_89"], outputs["bar_matched"])
    write_csv_exclusive(load_by_policy["mounted_as_scripted"], outputs["load_scripted"])
    write_csv_exclusive(load_by_policy["caption_matched_89"], outputs["load_matched"])
    assert figure5_summary is not None
    write_csv_exclusive(figure5_summary, outputs["figure5"])
    write_csv_exclusive(private, outputs["private"], private=True)
    print("Wrote 270 tests for each Figure 4 cohort policy, six Figure 4A bar")
    print("and 72 Figure 4B regional summaries per policy, 160 Figure 5 summaries,")
    print(f"and {len(private)} private participant rows.")


if __name__ == "__main__":
    main()
