#!/usr/bin/env python3
"""Plot old versus new EDSS-threshold AUCs without modifying source results."""

from __future__ import annotations

import csv
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
CANONICAL = PACKAGE_ROOT.parent / "revision/results/model_results/all_qmri_manuscript_models.csv"
SENSITIVITY = PACKAGE_ROOT / "results/model_results/all_qmri_manuscript_models.csv"
OUTPUT_DIR = PACKAGE_ROOT / "work/figures"

METRICS = ["T1", "MTR", "FA", "MD", "NDI", "ODI", "FWF"]
DISPLAY = ["T1", "MTR", "FA", "MD", "NDI", "ODI", "ISOVF"]
MODELS = ["Tract-based", "Classical-region"]


def load_rows(path: Path) -> dict[tuple[str, str], dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    selected = {
        (row["Metric"], row["Model"]): {
            "auc": float(row["Performance"]),
            "low": float(row["CI_lower"]),
            "high": float(row["CI_upper"]),
        }
        for row in rows
        if row["Outcome"] == "EDSS"
        and row["Adjusted"] == "0"
        and row["Metric"] in METRICS
        and row["Model"] in MODELS
    }
    if len(selected) != 14:
        raise AssertionError(f"Expected 14 EDSS rows in {path}, found {len(selected)}")
    return selected


def series(
    rows: dict[tuple[str, str], dict[str, float]], model: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    values = np.array([rows[(metric, model)]["auc"] for metric in METRICS])
    lower = np.array([rows[(metric, model)]["low"] for metric in METRICS])
    upper = np.array([rows[(metric, model)]["high"] for metric in METRICS])
    return values, values - lower, upper - values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical-models", type=Path, default=CANONICAL)
    parser.add_argument("--sensitivity-models", type=Path, default=SENSITIVITY)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    return parser.parse_args()


def local_output_dir(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    try:
        resolved.relative_to(PACKAGE_ROOT)
    except ValueError as exc:
        raise ValueError(f"Output must remain inside revisionExtras: {resolved}") from exc
    return resolved


def main() -> None:
    args = parse_args()
    canonical = args.canonical_models.expanduser().resolve()
    sensitivity = args.sensitivity_models.expanduser().resolve()
    output_dir = local_output_dir(args.output_dir)
    old = load_rows(canonical)
    new = load_rows(sensitivity)
    x = np.arange(len(METRICS))

    tract_old, tract_old_lo, tract_old_hi = series(old, "Tract-based")
    tract_new, tract_new_lo, tract_new_hi = series(new, "Tract-based")
    classic_old, classic_old_lo, classic_old_hi = series(old, "Classical-region")
    classic_new, classic_new_lo, classic_new_hi = series(new, "Classical-region")

    green = "#168A45"
    blue = "#276FBF"
    neutral = "#6B7280"

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 11,
            "axes.titlesize": 15,
            "axes.labelsize": 12,
            "legend.fontsize": 10,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, ax = plt.subplots(figsize=(10.8, 7.0))
    fig.subplots_adjust(left=0.09, right=0.985, top=0.82, bottom=0.25)

    plot_specs = [
        (tract_old, tract_old_lo, tract_old_hi, green, "-", "o", "Tract-based — old EDSS >3"),
        (tract_new, tract_new_lo, tract_new_hi, green, "--", "s", "Tract-based — new EDSS ≥4"),
        (classic_old, classic_old_lo, classic_old_hi, blue, "-", "o", "Classical-region — old EDSS >3"),
        (classic_new, classic_new_lo, classic_new_hi, blue, "--", "s", "Classical-region — new EDSS ≥4"),
    ]

    for values, lower_error, upper_error, color, linestyle, marker, label in plot_specs:
        ax.errorbar(
            x,
            values,
            yerr=np.vstack([lower_error, upper_error]),
            color=color,
            linestyle=linestyle,
            linewidth=2.25,
            marker=marker,
            markersize=6.2,
            markerfacecolor="white" if linestyle == "--" else color,
            markeredgecolor=color,
            markeredgewidth=1.5,
            elinewidth=1.0,
            capsize=3,
            alpha=0.98,
            label=label,
            zorder=4 if linestyle == "--" else 3,
        )

    ax.axhline(0.5, color=neutral, linestyle=":", linewidth=1.25, alpha=0.85)
    ax.text(
        6.0,
        0.512,
        "Chance AUC = 0.50",
        color=neutral,
        ha="right",
        va="bottom",
        fontsize=9.5,
    )
    ax.set_xticks(x, DISPLAY)
    ax.set_xlim(-0.28, 6.28)
    ax.set_ylim(0.40, 1.025)
    ax.set_yticks(np.arange(0.4, 1.01, 0.1))
    ax.set_xlabel("MRI metric")
    ax.set_ylabel("AUC (95% CI)")
    fig.suptitle(
        "EDSS threshold sensitivity: tract-based and classical-region AUC",
        x=0.5,
        y=0.965,
        fontsize=15,
    )
    fig.text(
        0.5,
        0.915,
        "Old: EDSS >3 (equivalent to ≥3.5 in this cohort); new: EDSS ≥4",
        ha="center",
        va="center",
        fontsize=10.5,
        color="#374151",
    )

    ax.grid(axis="y", color="#D1D5DB", linewidth=0.8, alpha=0.65)
    ax.grid(axis="x", visible=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#9CA3AF")
    ax.spines["bottom"].set_color("#9CA3AF")
    ax.tick_params(colors="#374151")
    fig.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 0.07),
        ncol=2,
        frameon=False,
        handlelength=3.1,
        columnspacing=1.8,
    )

    fig.text(
        0.5,
        0.026,
        "High/low groups: old 39/41 (T1/MTR) or 39/40 (other metrics); "
        "new 31/49 or 31/48.",
        ha="center",
        va="bottom",
        fontsize=9.5,
        color="#4B5563",
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    png = output_dir / "edss_threshold_auc_four_lines.png"
    pdf = output_dir / "edss_threshold_auc_four_lines.pdf"
    fig.savefig(png, dpi=240, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(png)
    print(pdf)


if __name__ == "__main__":
    main()
