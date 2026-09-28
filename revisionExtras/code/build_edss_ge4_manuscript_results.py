#!/usr/bin/env python3
"""Build an isolated manuscript-format EDSS >=4 results mirror.

The canonical revision/results tree is read-only.  The output tree has the
same inventory and layout; only the six EDSS-threshold-dependent CSV files
are changed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PACKAGE_ROOT.parent
CANONICAL = REPO_ROOT / "revision/results"
QMRI_SENSITIVITY = PACKAGE_ROOT / "work/qmri/edss_threshold_model_results.csv"
LESION_SENSITIVITY = (
    PACKAGE_ROOT / "work/lesion/lesion_edss_threshold_model_results.csv"
)
LESION_BASELINE_QA = (
    PACKAGE_ROOT / "work/qa/lesion/lesion_edss_baseline_reproduction.csv"
)
TARGET = PACKAGE_ROOT / "work/results"
MANIFEST = PACKAGE_ROOT / "work/qa/package_manifest.json"

CANONICAL_FILES = {
    "figures/Figure6.pdf",
    "figures/Supplementary_FigureS1.pdf",
    "figures/Supplementary_FigureS2.pdf",
    "model_results/all_qmri_manuscript_models.csv",
    "model_results/figure6_full_models.csv",
    "model_results/figure6_matched_models.csv",
    "model_results/lesionload_manuscript_models.csv",
    "model_results/supplementary_figure_s1_correlations.json",
    "tables/Main_Table_1.csv",
    "tables/Main_Table_2.csv",
    "tables/Supplementary_Table_3.csv",
    "tables/Supplementary_Table_4.csv",
    "tables/Supplementary_Table_5.csv",
}

CHANGED_FILES = {
    "model_results/all_qmri_manuscript_models.csv",
    "model_results/lesionload_manuscript_models.csv",
    "tables/Main_Table_2.csv",
    "tables/Supplementary_Table_3.csv",
    "tables/Supplementary_Table_4.csv",
    "tables/Supplementary_Table_5.csv",
}

METRICS = ["T1", "MTR", "FA", "MD", "NDI", "ODI", "FWF"]
DISPLAY_METRICS = {
    "T1": "T1",
    "MTR": "MTR",
    "FA": "FA",
    "MD": "MD",
    "NDI": "NDI",
    "ODI": "ODI",
    "FWF": "ISOVF",
}
MODELS = ["Tract-based", "Classical-region"]
LESION_KEYS = [
    ("LN", "Whole Brain"),
    ("LV", "Whole Brain"),
    ("LN", "Tract-based"),
    ("LN", "Classical-region"),
    ("LV", "Tract-based"),
    ("LV", "Classical-region"),
    ("Lnorm", "Tract-based"),
    ("Lnorm", "Classical-region"),
]


def read_dicts(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise AssertionError(f"No CSV header: {path}")
        return list(reader.fieldnames), list(reader)


def write_dicts(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            extrasaction="raise",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_rows(path: Path) -> list[list[str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.reader(handle))


def write_rows(path: Path, rows: list[list[str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tree_hashes(root: Path) -> dict[str, str]:
    return {
        str(path.relative_to(root)): sha256(path)
        for path in sorted(root.rglob("*"))
        if path.is_file() and path.name != ".DS_Store"
    }


def within_package(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    try:
        resolved.relative_to(PACKAGE_ROOT)
    except ValueError as exc:
        raise ValueError(f"{label} must remain inside revisionExtras: {resolved}") from exc
    return resolved


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical-results", type=Path, default=CANONICAL)
    parser.add_argument("--qmri-sensitivity", type=Path, default=QMRI_SENSITIVITY)
    parser.add_argument("--lesion-sensitivity", type=Path, default=LESION_SENSITIVITY)
    parser.add_argument("--lesion-baseline-qa", type=Path, default=LESION_BASELINE_QA)
    parser.add_argument("--output-dir", type=Path, default=TARGET)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    return parser.parse_args()


def compact_number(value: str | float, decimals: int) -> str:
    text = f"{float(value):.{decimals}f}".rstrip("0").rstrip(".")
    return "0" if text == "-0" else text


def p_text(value: str | float) -> str:
    number = float(value)
    return "< 0.0001" if number < 0.0001 else f"{number:.4f}"


def legacy_p_text(value: str | float) -> str:
    text = f"{float(value):.4f}"
    return "< 0.0001" if text == "0.0000" else text


def fmt_auc(row: dict[str, str], separator: str) -> str:
    return (
        f"{float(row['Performance']):.3f} "
        f"({float(row['CI_lower']):.3f}{separator}"
        f"{float(row['CI_upper']):.3f})"
    )


def diff_cell_count(old: list[list[str]], new: list[list[str]]) -> int:
    assert len(old) == len(new)
    assert all(len(a) == len(b) for a, b in zip(old, new))
    return sum(a != b for ra, rb in zip(old, new) for a, b in zip(ra, rb))


def patch_qmri(package: Path) -> dict[str, object]:
    path = package / "model_results/all_qmri_manuscript_models.csv"
    fields, rows = read_dicts(path)
    assert len(rows) == 196

    _, sensitivity_all = read_dicts(QMRI_SENSITIVITY)
    sensitivity = [
        row for row in sensitivity_all if row["ThresholdLabel"] == "sensitivity_ge4"
    ]
    assert len(sensitivity) == 28
    sensitivity_by_key = {
        (row["Metric"], row["Model"], row["Adjusted"]): row
        for row in sensitivity
    }
    assert len(sensitivity_by_key) == 28

    baseline = [
        row for row in sensitivity_all if row["ThresholdLabel"] == "canonical_gt3"
    ]
    baseline_by_key = {
        (row["Metric"], row["Model"], row["Adjusted"]): row for row in baseline
    }
    assert len(baseline_by_key) == 28

    mapped = {
        "N": "N",
        "PositiveN": "PositiveN",
        "Performance": "AUC",
        "CI_lower": "AUC_CI_Lower",
        "CI_upper": "AUC_CI_Upper",
        "Lambda": "Lambda",
        "p_raw": "p_raw",
        "AIC": "AIC",
    }
    edss_count = 0
    for row in rows:
        if row["Outcome"] != "EDSS":
            continue
        key = (row["Metric"], row["Model"], row["Adjusted"])
        assert key in sensitivity_by_key
        source = sensitivity_by_key[key]
        baseline_source = baseline_by_key[key]
        assert row["N"] == baseline_source["N"]
        assert row["PositiveN"] == baseline_source["PositiveN"]
        for destination, origin in mapped.items():
            row[destination] = source[origin]
        if row["Adjusted"] == "1":
            row["p_adj"] = format(min(98.0 * float(row["p_raw"]), 1.0), ".15g")
            assert row["MultiplicityFamilyN"] == "98"
        else:
            row["p_adj"] = "NaN"
            assert row["MultiplicityFamilyN"] == "0"
        edss_count += 1
    assert edss_count == 28
    assert {row["PositiveN"] for row in rows if row["Outcome"] == "EDSS"} == {"31"}
    write_dicts(path, fields, rows)

    adjusted = [row for row in rows if row["Adjusted"] == "1"]
    assert len(adjusted) == 98
    assert sum(float(row["p_adj"]) < 0.05 for row in adjusted) == 44
    edss_adjusted = [row for row in adjusted if row["Outcome"] == "EDSS"]
    assert sum(float(row["p_adj"]) < 0.05 for row in edss_adjusted) == 6
    return {
        "edss_rows_replaced": edss_count,
        "edss_adjusted_significant_m98": 6,
        "all_adjusted_significant_m98": 44,
        "qMRI_N": sorted({int(row["N"]) for row in rows if row["Outcome"] == "EDSS"}),
        "qMRI_positive_N": 31,
    }


def patch_lesion(package: Path) -> dict[str, object]:
    path = package / "model_results/lesionload_manuscript_models.csv"
    fields, rows = read_dicts(path)
    assert len(rows) == 56
    _, all_sensitivity = read_dicts(LESION_SENSITIVITY)
    sensitivity = [
        row for row in all_sensitivity if row["ThresholdLabel"] == "sensitivity_ge4"
    ]
    sensitivity_by_key = {(row["Metric"], row["Model"]): row for row in sensitivity}
    assert len(sensitivity_by_key) == 8

    replaced = 0
    for row in rows:
        if row["Outcome"] != "EDSS":
            continue
        key = (row["Metric"], row["Model"])
        assert key in sensitivity_by_key
        source = sensitivity_by_key[key]
        row["N"] = source["N"]
        row["PositiveN"] = source["PositiveN"]
        row["PerformanceName"] = source["PerformanceName"]
        row["Performance"] = compact_number(source["Performance"], 3)
        row["CI_lower"] = compact_number(source["CI_lower"], 3)
        row["CI_upper"] = compact_number(source["CI_upper"], 3)
        row["p_raw"] = source["p_raw"]
        row["AIC"] = compact_number(source["AIC"], 2)
        replaced += 1
    assert replaced == 8
    assert {row["N"] for row in rows if row["Outcome"] == "EDSS"} == {"89"}
    assert {row["PositiveN"] for row in rows if row["Outcome"] == "EDSS"} == {"32"}
    write_dicts(path, fields, rows)
    return {
        "edss_rows_replaced": replaced,
        "lesion_N": 89,
        "lesion_positive_N": 32,
        "lesion_raw_significant": sum(
            float(row["p_raw"]) < 0.05 for row in rows if row["Outcome"] == "EDSS"
        ),
    }


def qmri_lookup(package: Path) -> dict[tuple[str, str, str], dict[str, str]]:
    _, rows = read_dicts(package / "model_results/all_qmri_manuscript_models.csv")
    return {
        (row["Metric"], row["Model"], row["Adjusted"]): row
        for row in rows
        if row["Outcome"] == "EDSS"
    }


def patch_qmri_auc_table(package: Path, relative: str, adjusted: str, sep: str) -> None:
    path = package / relative
    rows = read_rows(path)
    assert len(rows) == 15 and len(rows[0]) == 5
    lookup = qmri_lookup(package)
    replacement: list[list[str]] = []
    for metric in METRICS:
        tract = lookup[(metric, "Tract-based", adjusted)]
        classical = lookup[(metric, "Classical-region", adjusted)]
        assert tract["N"] == classical["N"]
        assert tract["PositiveN"] == classical["PositiveN"]
        replacement.append(
            [
                f"EDSS-{DISPLAY_METRICS[metric]} (n={tract['N']})",
                fmt_auc(tract, sep),
                fmt_auc(classical, sep),
                f"{float(tract['AIC']):.2f}",
                f"{float(classical['AIC']):.2f}",
            ]
        )
    rows[1:8] = replacement
    write_rows(path, rows)


def patch_supplementary_table_4(package: Path) -> None:
    path = package / "tables/Supplementary_Table_4.csv"
    rows = read_rows(path)
    assert len(rows) == 99 and len(rows[0]) == 6
    lookup = qmri_lookup(package)
    replacement: list[list[str]] = []
    for metric_index, metric in enumerate(METRICS):
        for model_index, model in enumerate(MODELS):
            row = lookup[(metric, model, "1")]
            replacement.append(
                [
                    "EDSS" if metric_index == 0 and model_index == 0 else "",
                    DISPLAY_METRICS[metric] if model_index == 0 else "",
                    "Tract-based" if model == "Tract-based" else "Classic",
                    f"{row['N']} ({row['PositiveN']})",
                    (
                        legacy_p_text(row["p_raw"])
                        if metric in {"T1", "MTR", "FA", "MD"}
                        else p_text(row["p_raw"])
                    ),
                    p_text(row["p_adj"]),
                ]
            )
    assert len(replacement) == 14
    rows[1:15] = replacement
    write_rows(path, rows)


def patch_supplementary_table_5(package: Path) -> None:
    path = package / "tables/Supplementary_Table_5.csv"
    rows = read_rows(path)
    assert len(rows) == 58 and len(rows[0]) == 6
    _, lesion = read_dicts(package / "model_results/lesionload_manuscript_models.csv")
    lookup = {
        (row["Metric"], row["Model"]): row for row in lesion if row["Outcome"] == "EDSS"
    }
    assert len(lookup) == 8
    replacement: list[list[str]] = []
    for index, key in enumerate(LESION_KEYS):
        row = lookup[key]
        replacement.append(
            [
                "EDSS" if index == 0 else "",
                row["Metric"],
                "Classic" if row["Model"] == "Classical-region" else row["Model"],
                f"{row['N']} ({row['PositiveN']})",
                (
                    f"{float(row['Performance']):.3f} "
                    f"[{float(row['CI_lower']):.3f}, {float(row['CI_upper']):.3f}]"
                ),
                f"{float(row['AIC']):.2f}",
            ]
        )
    rows[1:9] = replacement
    write_rows(path, rows)


def validate_package(
    package: Path, canonical_before: dict[str, str]
) -> tuple[dict[str, str], dict[str, object]]:
    canonical_after = tree_hashes(CANONICAL)
    assert canonical_after == canonical_before, "Canonical results changed during build"
    package_hashes = tree_hashes(package)
    assert set(package_hashes) == set(canonical_before), "Package inventory differs"
    actual_changed = {
        relative
        for relative, digest in package_hashes.items()
        if digest != canonical_before[relative]
    }
    assert actual_changed == CHANGED_FILES, actual_changed

    non_edss_checks: dict[str, bool] = {}
    for relative in [
        "model_results/all_qmri_manuscript_models.csv",
        "model_results/lesionload_manuscript_models.csv",
    ]:
        _, old = read_dicts(CANONICAL / relative)
        _, new = read_dicts(package / relative)
        old_non = [row for row in old if row["Outcome"] != "EDSS"]
        new_non = [row for row in new if row["Outcome"] != "EDSS"]
        non_edss_checks[relative] = old_non == new_non
        assert non_edss_checks[relative]

    table_unchanged_slices = {
        "tables/Main_Table_2.csv": 8,
        "tables/Supplementary_Table_3.csv": 8,
        "tables/Supplementary_Table_4.csv": 15,
        "tables/Supplementary_Table_5.csv": 9,
    }
    for relative, start in table_unchanged_slices.items():
        old = read_rows(CANONICAL / relative)
        new = read_rows(package / relative)
        non_edss_checks[relative] = old[start:] == new[start:]
        assert non_edss_checks[relative]

    _, lesion_qa = read_dicts(LESION_BASELINE_QA)
    assert len(lesion_qa) == 8 and all(row["Pass"] == "1" for row in lesion_qa)

    cell_diffs = {
        relative: diff_cell_count(
            read_rows(CANONICAL / relative), read_rows(package / relative)
        )
        for relative in sorted(CHANGED_FILES)
    }
    assert cell_diffs["model_results/all_qmri_manuscript_models.csv"] == 207
    assert cell_diffs["tables/Main_Table_2.csv"] == 28
    assert cell_diffs["tables/Supplementary_Table_3.csv"] == 28
    assert cell_diffs["tables/Supplementary_Table_4.csv"] == 38

    validation = {
        "status": "PASS",
        "canonical_tree_unchanged": True,
        "same_file_inventory": True,
        "file_count": len(package_hashes),
        "changed_file_count": len(actual_changed),
        "unchanged_file_count": len(package_hashes) - len(actual_changed),
        "changed_files": sorted(actual_changed),
        "non_EDSS_rows_unchanged": non_edss_checks,
        "lesion_canonical_baseline_rows_passed": 8,
        "changed_cell_counts": cell_diffs,
    }
    return package_hashes, validation


def main() -> None:
    global CANONICAL, QMRI_SENSITIVITY, LESION_SENSITIVITY
    global LESION_BASELINE_QA, TARGET, MANIFEST

    args = parse_args()
    CANONICAL = args.canonical_results.expanduser().resolve()
    QMRI_SENSITIVITY = args.qmri_sensitivity.expanduser().resolve()
    LESION_SENSITIVITY = args.lesion_sensitivity.expanduser().resolve()
    LESION_BASELINE_QA = args.lesion_baseline_qa.expanduser().resolve()
    TARGET = within_package(args.output_dir, "output-dir")
    MANIFEST = within_package(args.manifest, "manifest")

    assert CANONICAL.is_dir()
    assert QMRI_SENSITIVITY.is_file()
    assert LESION_SENSITIVITY.is_file()
    assert LESION_BASELINE_QA.is_file()
    if TARGET.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {TARGET}")

    canonical_before = tree_hashes(CANONICAL)
    assert set(canonical_before) == CANONICAL_FILES
    TARGET.parent.mkdir(parents=True, exist_ok=True)
    temp_root = Path(tempfile.mkdtemp(prefix="edss_ge4_full_", dir=TARGET.parent))
    stage = temp_root / TARGET.name
    try:
        shutil.copytree(
            CANONICAL,
            stage,
            copy_function=shutil.copy2,
            ignore=shutil.ignore_patterns(".DS_Store"),
        )
        qmri_summary = patch_qmri(stage)
        lesion_summary = patch_lesion(stage)
        patch_qmri_auc_table(stage, "tables/Main_Table_2.csv", "0", "-")
        patch_qmri_auc_table(stage, "tables/Supplementary_Table_3.csv", "1", "–")
        patch_supplementary_table_4(stage)
        patch_supplementary_table_5(stage)
        package_hashes, validation = validate_package(stage, canonical_before)
        os.rename(stage, TARGET)
    finally:
        if temp_root.exists():
            shutil.rmtree(temp_root)

    manifest = {
        "schema_version": 1,
        "status": "PASS",
        "scope": "Complete manuscript-format threshold-only EDSS >=4 sensitivity package",
        "canonical_source": str(CANONICAL),
        "isolated_package": str(TARGET),
        "threshold_change": {
            "canonical": "EDSS > 3.0 (equivalent to EDSS >=3.5 in these data)",
            "sensitivity": "EDSS >= 4.0",
        },
        "cohort_policy": (
            "Hold each canonical analysis cohort fixed and change only the EDSS threshold. "
            "For Supplementary Table 5 this intentionally retains RR047 in the frozen "
            "canonical lesion cohort (N=89); removing RR047 would conflate threshold and "
            "cohort changes."
        ),
        "qMRI": qmri_summary,
        "lesion_load": lesion_summary,
        "multiplicity": (
            "Adjusted qMRI p-values use the manuscript-wide Bonferroni family m=98; "
            "the exploratory m=14 values are not inserted into this package."
        ),
        "validation": validation,
        "canonical_hashes": canonical_before,
        "package_hashes": package_hashes,
    }
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    temporary_manifest = MANIFEST.with_suffix(".json.tmp")
    temporary_manifest.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary_manifest, MANIFEST)
    print(json.dumps(validation, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
