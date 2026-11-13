#!/usr/bin/env python3
"""Build Figure S3 and its aggregate tables after both Figure 5 FA gates."""

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


CLASSICAL_REGIONS = ("periventricular", "juxtacortical", "infratentorial", "deepwhitematter")
TRACT_REGIONS = ("Association", "Cerebellar", "Occipitoparietal", "ProjectionBrainstem")
PANELS = (
    ("Classical", "Lesion", CLASSICAL_REGIONS, "top classical lesion"),
    ("Tract", "Lesion", TRACT_REGIONS, "middle tract lesion"),
    ("Tract", "NAWM", TRACT_REGIONS, "bottom tract NAWM"),
)
GROUPS = ("RRMS", "SPMS", "PPMS", "Control")
PREFIX = {"RRMS": "sub-RR", "SPMS": "sub-SP", "PPMS": "sub-PP", "Control": "sub-C"}
METRIC_FIELDS = {"NDI": "NDI", "ODI": "ODI", "ISOVF": "FWF"}
TOLERANCE = 1e-6
SUMMARY_COLUMNS = [
    "Framework", "Tissue", "Metric", "Region", "Group", "N", "Mean", "SD", "ZeroCount",
    "KW_H", "KW_df", "KW_p", "KW_p_BH",
]
PAIR_COLUMNS = ["Framework", "Tissue", "Metric", "Region", "GroupA", "GroupB", "nA", "nB", "U", "p", "p_BH"]
PLOT_COLUMNS = [
    "Figure", "Panel", "CohortPolicy", "Framework", "Tissue", "Metric", "Region",
    "Subtype", "N", "Mean", "SD", "SEM", "Median", "Q1", "Q3", "ZeroCount",
    "ZeroPolicy", "SourceSHA256",
]


def load_unique(path: Path, columns: list[str] | None = None) -> pd.DataFrame:
    frame = pd.read_csv(path, usecols=columns, dtype={"SubjectID": str})
    if "SubjectID" not in frame or frame.SubjectID.isna().any() or frame.SubjectID.duplicated().any():
        raise ValueError(f"Invalid subject key in {path.name}")
    return frame.set_index("SubjectID", verify_integrity=True)


def groups_for(tissue: str) -> tuple[str, ...]:
    return GROUPS[:3] if tissue == "Lesion" else GROUPS


def values_for(frame: pd.DataFrame, group: str, column: str) -> np.ndarray:
    selected = pd.to_numeric(frame.loc[frame.index.str.startswith(PREFIX[group]), column], errors="raise")
    values = selected.dropna().to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(f"Non-finite values for {group}/{column}")
    return values


def gate_fa(fa: pd.DataFrame, frozen: pd.DataFrame, framework: str, expected_cells: int) -> dict:
    target = frozen[(frozen.Framework == framework) & (frozen.Metric == "FA")]
    checks = []
    for panel_framework, tissue, regions, _ in PANELS:
        if panel_framework != framework:
            continue
        for region in regions:
            column = f"{region}FA" if framework == "Classical" else f"{region}_{tissue}"
            for group in groups_for(tissue):
                matches = target[
                    (target.Region == region) & (target.Tissue == tissue) & (target.Subtype == group)
                ]
                if len(matches) != 1:
                    raise ValueError(f"{framework} FA gate: expected one reference cell for {tissue}/{region}/{group}")
                ref = matches.iloc[0]
                x = values_for(fa, group, column)
                x = x[x != 0]  # Frozen Figure 5 policy for FA.
                if len(x) < 3:
                    raise ValueError(f"{framework} FA gate: N < 3 for {tissue}/{region}/{group}")
                sd = x.std(ddof=1)
                observed = (len(x), x.mean(), sd, sd / math.sqrt(len(x)))
                expected = (ref.N, ref.Mean, ref.SD, ref.SEM)
                differences = np.abs(np.asarray(observed) - np.asarray(expected))
                if not np.isfinite(differences).all():
                    raise ValueError(f"{framework} FA gate: invalid comparison for {tissue}/{region}/{group}")
                checks.append((float(differences.max()), tissue, region, group))
    if len(checks) != expected_cells or len(target) != expected_cells:
        raise ValueError(f"{framework} FA gate: expected {expected_cells} cells")
    maximum, tissue, region, group = max(checks)
    if maximum > TOLERANCE:
        raise ValueError(f"{framework} FA gate failed for {tissue}/{region}/{group}: "
                         f"maximum difference {maximum:.12g} > {TOLERANCE:g}")
    return {"cells": len(checks), "maximum_absolute_difference": maximum,
            "tolerance": TOLERANCE, "passed": True}


def bh(values: list[float]) -> np.ndarray:
    return multipletests(values, method="fdr_bh")[1]


def mw(x: np.ndarray, y: np.ndarray):
    return mannwhitneyu(x, y, alternative="two-sided", method="asymptotic")


def build_results(classical: pd.DataFrame, tracts: dict[str, pd.DataFrame], hashes: dict[str, str]):
    summary_rows = []
    pair_rows = []
    plot_rows = []
    omnibus_keys = []
    omnibus_p = []
    omnibus_zero_p = []
    pair_p = []
    pair_zero_p = []
    for framework, tissue, regions, panel in PANELS:
        groups = groups_for(tissue)
        for metric, field in METRIC_FIELDS.items():
            frame = classical if framework == "Classical" else tracts[metric]
            digest = hashes["Classical" if framework == "Classical" else metric]
            for region in regions:
                column = f"{region}{field}" if framework == "Classical" else f"{region}_{tissue}"
                cohorts = {g: values_for(frame, g, column) for g in groups}
                for group, values in cohorts.items():
                    if len(values) < 3:
                        raise ValueError(f"N < 3 for {framework}/{tissue}/{metric}/{region}/{group}: {len(values)}")
                without_zero = {g: x[x != 0] for g, x in cohorts.items()}
                if any(len(x) < 3 for x in without_zero.values()):
                    raise ValueError(f"Zero-excluded N < 3 for {framework}/{tissue}/{metric}/{region}")
                h, p = kruskal(*cohorts.values(), nan_policy="raise")
                omnibus_keys.append((framework, tissue, metric, region))
                omnibus_p.append(float(p))
                omnibus_zero_p.append(float(kruskal(*without_zero.values(), nan_policy="raise").pvalue))
                feature = dict(Framework=framework, Tissue=tissue, Metric=metric, Region=region)
                for group, x in cohorts.items():
                    n, mean, sd = len(x), float(x.mean()), float(x.std(ddof=1))
                    zero_count = int(np.count_nonzero(x == 0))
                    summary_rows.append(dict(feature, Group=group, N=n, Mean=mean, SD=sd,
                                             ZeroCount=zero_count, KW_H=h, KW_df=len(groups)-1, KW_p=p))
                    plot_rows.append(dict(feature, Figure="FigureS3", Panel=panel, CohortPolicy="Figure5",
                                          Subtype=group, N=n, Mean=mean, SD=sd, SEM=sd/math.sqrt(n),
                                          Median=np.median(x), Q1=np.quantile(x, .25), Q3=np.quantile(x, .75),
                                          ZeroCount=zero_count, ZeroPolicy="retain", SourceSHA256=digest))
                for a, b in combinations(groups, 2):
                    u, p = mw(cohorts[a], cohorts[b])
                    pair_p.append(float(p))
                    pair_zero_p.append(float(mw(without_zero[a], without_zero[b]).pvalue))
                    pair_rows.append(dict(feature, GroupA=a, GroupB=b, nA=len(cohorts[a]), nB=len(cohorts[b]),
                                          U=u, p=p))
    if (len(omnibus_p), len(pair_p)) != (36, 144):
        raise AssertionError(f"Expected 36 omnibus and 144 pairwise tests, got {len(omnibus_p)} and {len(pair_p)}")

    omnibus_q, omnibus_zero_q = bh(omnibus_p), bh(omnibus_zero_p)
    pair_q, pair_zero_q = bh(pair_p), bh(pair_zero_p)
    q_by_feature = dict(zip(omnibus_keys, omnibus_q))
    for row in summary_rows:
        key = (row["Framework"], row["Tissue"], row["Metric"], row["Region"])
        row["KW_p_BH"] = q_by_feature[key]
    for row, q in zip(pair_rows, pair_q):
        row["p_BH"] = q
    sensitivity = {
        "omnibus_maximum_BH_p_difference": float(np.max(np.abs(omnibus_q - omnibus_zero_q))),
        "pairwise_maximum_BH_p_difference": float(np.max(np.abs(pair_q - pair_zero_q))),
        "omnibus_significance_changes_at_0.05": int(np.count_nonzero((omnibus_q < .05) != (omnibus_zero_q < .05))),
        "pairwise_significance_changes_at_0.05": int(np.count_nonzero((pair_q < .05) != (pair_zero_q < .05))),
    }
    return (pd.DataFrame(summary_rows, columns=SUMMARY_COLUMNS), pd.DataFrame(pair_rows, columns=PAIR_COLUMNS),
            pd.DataFrame(plot_rows, columns=PLOT_COLUMNS), sensitivity)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics-dir", type=Path, required=True)
    parser.add_argument("--figure5", type=Path, required=True, help="Frozen Figure 5 plotted-values CSV")
    parser.add_argument("--output-dir", type=Path, required=True, help="Fresh directory for aggregate CSVs")
    parser.add_argument("--internal-dir", type=Path, required=True, help="Private report directory")
    args = parser.parse_args()
    outputs = [args.output_dir / name for name in (
        "Supplementary_Table_NODDI_group_summary.csv", "Supplementary_Table_NODDI_pairwise.csv",
        "figureS3_plotted_values.csv")]
    if any(path.exists() for path in outputs):
        raise FileExistsError("Refusing to overwrite an existing aggregate result")

    frozen = pd.read_csv(args.figure5)
    tract_fa = load_unique(args.metrics_dir / "GroupTractFA_All.csv")
    if len(tract_fa) != 132:
        raise ValueError(f"Tract FA gate input must have 132 subjects, found {len(tract_fa)}")
    tract_gate = gate_fa(tract_fa, frozen, "Tract", 28)
    # Read only FA columns from the combined file before the classical gate.
    classical_fa = load_unique(args.metrics_dir / "ClassicalRegionAllMetrics.csv",
                               ["SubjectID", *(f"{r}FA" for r in CLASSICAL_REGIONS)])
    classical_gate = gate_fa(classical_fa, frozen, "Classical", 12)

    paths = {metric: args.metrics_dir / f"GroupTract{field}_All.csv"
             for metric, field in METRIC_FIELDS.items()}
    paths["Classical"] = args.metrics_dir / "ClassicalRegionNODDI_All.csv"
    tracts = {metric: load_unique(paths[metric]) for metric in METRIC_FIELDS}
    classical = load_unique(paths["Classical"])
    tract_subjects = set(tracts["NDI"].index)
    if any(set(frame.index) != tract_subjects for frame in tracts.values()):
        raise ValueError("Tract NODDI inputs do not contain the same subjects")
    if not tract_subjects.issubset(tract_fa.index) or not set(classical.index).issubset(classical_fa.index):
        raise ValueError("NODDI inputs include subjects outside their Figure 5 cohort")
    hashes = {key: hashlib.sha256(path.read_bytes()).hexdigest() for key, path in paths.items()}
    summary, pairwise, plotted, sensitivity = build_results(classical, tracts, hashes)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for result, path in zip((summary, pairwise, plotted), outputs):
        result.to_csv(path, index=False)

    report = {
        "fa_gates": {"Tract": tract_gate, "Classical": classical_gate},
        "noddi_subjects": {"Tract": len(tract_subjects), "Classical": len(classical)},
        "family_sizes": {"Kruskal_Wallis": len(summary.drop_duplicates(["Framework","Tissue","Metric","Region"])),
                         "pairwise_Mann_Whitney": len(pairwise)},
        "zero_excluded_sensitivity": sensitivity,
    }
    args.internal_dir.mkdir(parents=True, exist_ok=True)
    (args.internal_dir / "fig5_counterpart_analysis.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
