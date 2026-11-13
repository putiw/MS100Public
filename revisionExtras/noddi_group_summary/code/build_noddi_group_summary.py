#!/usr/bin/env python3
"""Build the aggregate NODDI group table (Tail and NAWM) and the Figure S3 plotted values
(lesion mean and NAWM), after the Figure 5 FA equivalence gate."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests


CLUSTERS = ("Association", "Cerebellar", "Occipitoparietal", "ProjectionBrainstem")
GROUPS = ("RRMS", "SPMS", "PPMS", "Control")
LESION_GROUPS = GROUPS[:3]
PREFIX = {"RRMS": "sub-RR", "SPMS": "sub-SP", "PPMS": "sub-PP", "Control": "sub-C"}
METRIC_FILES = {"NDI": "GroupTractNDI_All.csv", "ODI": "GroupTractODI_All.csv", "ISOVF": "GroupTractFWF_All.csv"}
FA_FILE = "GroupTractFA_All.csv"
TABLE_TISSUES = ("Tail", "NAWM")
FIGURE_TISSUES = ("Lesion", "NAWM")
SUMMARY_COLUMNS = [
    "Tissue", "Metric", "Cluster", "Group", "N", "Mean", "SD", "ZeroCount",
    "KW_H", "KW_df", "KW_p", "KW_p_BH",
]
PAIR_COLUMNS = ["Tissue", "Metric", "Cluster", "GroupA", "GroupB", "nA", "nB", "U", "p", "p_BH"]
PLOT_COLUMNS = [
    "Figure", "Panel", "CohortPolicy", "Framework", "Tissue", "Metric", "Region",
    "Subtype", "N", "Mean", "SD", "SEM", "Median", "Q1", "Q3", "ZeroCount",
    "ZeroPolicy", "SourceSHA256",
]


def load_unique(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    if "SubjectID" not in frame or frame.SubjectID.isna().any() or frame.SubjectID.duplicated().any():
        raise ValueError(f"Invalid subject key in {path.name}")
    frame["SubjectID"] = frame.SubjectID.astype(str)
    return frame.set_index("SubjectID", verify_integrity=True)


def groups_for(tissue: str) -> tuple[str, ...]:
    # Controls have no lesions; Tail and NAWM are defined for all four groups.
    return LESION_GROUPS if tissue == "Lesion" else GROUPS


def values_for(frame: pd.DataFrame, group: str, cluster: str, tissue: str) -> np.ndarray:
    column = f"{cluster}_{tissue}"
    if column not in frame:
        raise ValueError(f"Missing aggregate input column {column}")
    selected = pd.to_numeric(frame.loc[frame.index.str.startswith(PREFIX[group]), column], errors="raise")
    result = selected.dropna().to_numpy(dtype=float)
    if not np.isfinite(result).all():
        raise ValueError(f"Non-finite values for {group}/{cluster}/{tissue}")
    return result


def gate_fa(fa: pd.DataFrame, figure5_path: Path) -> tuple[float, list[dict]]:
    frozen = pd.read_csv(figure5_path)
    target = frozen[(frozen.Framework == "Tract") & (frozen.Metric == "FA")]
    observed = []
    for cluster in CLUSTERS:
        for tissue in FIGURE_TISSUES:
            for group in groups_for(tissue):
                matches = target[
                    (target.Region == cluster) & (target.Tissue == tissue) & (target.Subtype == group)
                ]
                if len(matches) != 1:
                    raise ValueError(f"Figure 5 FA gate has {len(matches)} rows for {tissue}/{cluster}/{group}")
                ref = matches.iloc[0]
                x = values_for(fa, group, cluster, tissue)
                # Figure 5's frozen aggregate policy removes exact-zero qMRI observations.
                x = x[x != 0]
                stats = (len(x), float(np.mean(x)), float(np.std(x, ddof=1)), float(np.std(x, ddof=1) / math.sqrt(len(x))))
                expected = (int(ref.N), float(ref.Mean), float(ref.SD), float(ref.SEM))
                delta = max(abs(a - b) for a, b in zip(stats, expected))
                observed.append({"Tissue": tissue, "Cluster": cluster, "Group": group,
                                 "N": stats[0], "MaximumAbsoluteDifference": delta})
    if len(observed) != 28:
        raise AssertionError(f"Expected 28 FA gate cells, found {len(observed)}")
    max_difference = max(row["MaximumAbsoluteDifference"] for row in observed)
    if max_difference > 1e-6:
        raise RuntimeError(f"FA gate failed: maximum absolute difference {max_difference:.12g} > 1e-6")
    return max_difference, observed


def bh(values: list[float]) -> np.ndarray:
    return multipletests(values, method="fdr_bh")[1]


def build_table(sources: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Kruskal-Wallis and pairwise Mann-Whitney for Tail and NAWM; one BH family each."""
    summary_rows: list[dict] = []
    pair_rows: list[dict] = []
    omnibus_keys: list[tuple[str, str, str]] = []
    omnibus_p: list[float] = []
    omnibus_zero_p: list[float] = []
    pair_keys: list[tuple[str, ...]] = []
    pair_p: list[float] = []
    pair_zero_p: list[float] = []
    for metric, frame in sources.items():
        for tissue in TABLE_TISSUES:
            groups = groups_for(tissue)
            for cluster in CLUSTERS:
                cohorts = {group: values_for(frame, group, cluster, tissue) for group in groups}
                if any(len(x) < 3 for x in cohorts.values()):
                    raise ValueError(f"Empty or N < 3 cell: {metric}/{tissue}/{cluster}")
                h, p = kruskal(*(cohorts[g] for g in groups), nan_policy="raise")
                without_zero = {g: x[x != 0] for g, x in cohorts.items()}
                if any(len(x) < 3 for x in without_zero.values()):
                    raise ValueError(f"Zero-excluded sensitivity cell has N < 3: {metric}/{tissue}/{cluster}")
                _, p_zero = kruskal(*(without_zero[g] for g in groups), nan_policy="raise")
                key = (tissue, metric, cluster)
                omnibus_keys.append(key)
                omnibus_p.append(float(p))
                omnibus_zero_p.append(float(p_zero))
                for group in groups:
                    x = cohorts[group]
                    summary_rows.append({
                        "Tissue": tissue, "Metric": metric, "Cluster": cluster, "Group": group,
                        "N": len(x), "Mean": float(np.mean(x)), "SD": float(np.std(x, ddof=1)),
                        "ZeroCount": int(np.count_nonzero(x == 0)), "KW_H": float(h),
                        "KW_df": len(groups) - 1, "KW_p": float(p),
                    })
                for group_a, group_b in combinations(groups, 2):
                    xa, xb = cohorts[group_a], cohorts[group_b]
                    u, p_raw = mannwhitneyu(xa, xb, alternative="two-sided", method="asymptotic")
                    _, p_pair_zero = mannwhitneyu(without_zero[group_a], without_zero[group_b],
                                                  alternative="two-sided", method="asymptotic")
                    pair_keys.append((*key, group_a, group_b))
                    pair_p.append(float(p_raw))
                    pair_zero_p.append(float(p_pair_zero))
                    pair_rows.append({
                        "Tissue": tissue, "Metric": metric, "Cluster": cluster,
                        "GroupA": group_a, "GroupB": group_b, "nA": len(xa), "nB": len(xb),
                        "U": float(u), "p": float(p_raw),
                    })
    if len(omnibus_p) != 24 or len(pair_p) != 144:
        raise AssertionError(f"Unexpected family sizes: {len(omnibus_p)} omnibus, {len(pair_p)} pairwise")
    q_omnibus, q_omnibus_zero = bh(omnibus_p), bh(omnibus_zero_p)
    q_pair, q_pair_zero = bh(pair_p), bh(pair_zero_p)
    q_map = dict(zip(omnibus_keys, q_omnibus))
    for row in summary_rows:
        row["KW_p_BH"] = float(q_map[(row["Tissue"], row["Metric"], row["Cluster"])])
    for row, q in zip(pair_rows, q_pair):
        row["p_BH"] = float(q)
    changed_omnibus = ["|".join(k) for k, a, b in zip(omnibus_keys, q_omnibus < .05, q_omnibus_zero < .05) if a != b]
    changed_pairs = ["|".join(k) for k, a, b in zip(pair_keys, q_pair < .05, q_pair_zero < .05) if a != b]
    sensitivity = {
        "Kruskal_Wallis_BH_significance_changed": changed_omnibus,
        "pairwise_BH_significance_changed": changed_pairs,
        "any_BH_significant_conclusion_changed": bool(changed_omnibus or changed_pairs),
    }
    return (pd.DataFrame(summary_rows, columns=SUMMARY_COLUMNS),
            pd.DataFrame(pair_rows, columns=PAIR_COLUMNS), sensitivity)


def build_plotted_values(sources: dict[str, pd.DataFrame], metrics_dir: Path) -> pd.DataFrame:
    """Means and SEMs drawn in Supplementary Figure S3 (lesion mean and NAWM rows)."""
    rows: list[dict] = []
    for metric, frame in sources.items():
        digest = hashlib.sha256((metrics_dir / METRIC_FILES[metric]).read_bytes()).hexdigest()
        for tissue in FIGURE_TISSUES:
            for cluster in CLUSTERS:
                for group in groups_for(tissue):
                    x = values_for(frame, group, cluster, tissue)
                    if len(x) < 3:
                        raise ValueError(f"Empty or N < 3 cell: {metric}/{tissue}/{cluster}/{group}")
                    rows.append({
                        "Figure": "FigureS3", "Panel": f"{tissue} tract", "CohortPolicy": "Figure5",
                        "Framework": "Tract", "Tissue": tissue, "Metric": metric, "Region": cluster,
                        "Subtype": group, "N": len(x), "Mean": float(np.mean(x)),
                        "SD": float(np.std(x, ddof=1)), "SEM": float(np.std(x, ddof=1) / math.sqrt(len(x))),
                        "Median": float(np.median(x)), "Q1": float(np.quantile(x, .25)),
                        "Q3": float(np.quantile(x, .75)), "ZeroCount": int(np.count_nonzero(x == 0)),
                        "ZeroPolicy": "retain", "SourceSHA256": digest,
                    })
    return pd.DataFrame(rows, columns=PLOT_COLUMNS)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics-dir", type=Path, required=True, help="Directory containing per-subject GroupTract CSVs")
    parser.add_argument("--figure5", type=Path, required=True, help="Frozen Figure 5 plotted-values CSV")
    parser.add_argument("--output-dir", type=Path, required=True, help="Results directory for the aggregate CSVs")
    parser.add_argument("--internal-dir", type=Path, required=True, help="Private report directory outside the repository")
    args = parser.parse_args()
    outputs = {
        "summary": args.output_dir / "Supplementary_Table_NODDI_group_summary.csv",
        "pairwise": args.output_dir / "Supplementary_Table_NODDI_pairwise.csv",
        "plotted": args.output_dir / "figureS3_plotted_values.csv",
    }
    if any(path.exists() for path in outputs.values()):
        raise FileExistsError("Refusing to overwrite an existing aggregate result")

    # The NODDI tables are deliberately not opened until this gate passes.
    fa = load_unique(args.metrics_dir / FA_FILE)
    if len(fa) != 132:
        raise ValueError(f"FA gate input must have 132 Figure 5 subjects, found {len(fa)}")
    fa_difference, gate_cells = gate_fa(fa, args.figure5)

    sources = {metric: load_unique(args.metrics_dir / filename) for metric, filename in METRIC_FILES.items()}
    subject_sets = [set(frame.index) for frame in sources.values()]
    if any(ids != subject_sets[0] for ids in subject_sets[1:]):
        raise ValueError("NDI, ODI and ISOVF inputs do not contain the same subjects")
    if not subject_sets[0].issubset(set(fa.index)):
        raise ValueError("NODDI inputs include subjects outside the Figure 5 cohort")

    summary, pairwise, sensitivity = build_table(sources)
    plotted = build_plotted_values(sources, args.metrics_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(outputs["summary"], index=False)
    pairwise.to_csv(outputs["pairwise"], index=False)
    plotted.to_csv(outputs["plotted"], index=False)

    report = {
        "fa_gate": {"cells": len(gate_cells), "maximum_absolute_difference": fa_difference,
                    "tolerance": 1e-6, "passed": True},
        "noddi_subjects": len(subject_sets[0]),
        "family_sizes": {"Kruskal_Wallis": int(summary.drop_duplicates(["Tissue", "Metric", "Cluster"]).shape[0]),
                         "pairwise_Mann_Whitney": len(pairwise)},
        "zero_excluded_sensitivity": sensitivity,
    }
    args.internal_dir.mkdir(parents=True, exist_ok=True)
    (args.internal_dir / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
