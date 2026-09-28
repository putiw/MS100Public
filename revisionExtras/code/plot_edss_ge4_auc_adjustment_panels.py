#!/usr/bin/env python3
"""Plot EDSS-threshold AUCs side by side with and without demographics."""

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


def load_models(path: Path) -> dict[tuple[int, str, str], dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    selected = {
        (int(row["Adjusted"]), row["Metric"], row["Model"]): {
            "auc": float(row["Performance"]),
            "low": float(row["CI_lower"]),
            "high": float(row["CI_upper"]),
        }
        for row in rows
        if row["Outcome"] == "EDSS"
        and row["Metric"] in METRICS
        and row["Model"] in MODELS
    }
    if len(selected) != 28:
        raise AssertionError(f"Expected 28 EDSS rows in {path}, found {len(selected)}")
    return selected


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


def series(
    rows: dict[tuple[int, str, str], dict[str, float]],
    adjusted: int,
    model: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    values = np.array([rows[(adjusted, metric, model)]["auc"] for metric in METRICS])
    lower = np.array([rows[(adjusted, metric, model)]["low"] for metric in METRICS])
    upper = np.array([rows[(adjusted, metric, model)]["high"] for metric in METRICS])
    return values, values - lower, upper - values


def main() -> None:
    args = parse_args()
    canonical = args.canonical_models.expanduser().resolve()
    sensitivity = args.sensitivity_models.expanduser().resolve()
    output_dir = local_output_dir(args.output_dir)
    old = load_models(canonical)
    new = load_models(sensitivity)
    x = np.arange(len(METRICS))

    green = "#168A45"
    blue = "#276FBF"
    neutral = "#6B7280"
    grid = "#D1D5DB"
    text = "#374151"

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 12,
            "legend.fontsize": 10,
            "xtick.labelsize": 10.5,
            "ytick.labelsize": 10.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, axes = plt.subplots(1, 2, figsize=(16.0, 7.1), sharey=True)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.80, bottom=0.25, wspace=0.10)

    panel_titles = [
        "A. Without demographic covariates",
        "B. Adjusted for age, gender, and disease duration",
    ]
    legend_handles = None
    legend_labels = None

    for adjusted, ax, panel_title in zip([0, 1], axes, panel_titles):
        plot_specs = [
            (*series(old, adjusted, "Tract-based"), green, "-", "o", "Tract-based — old EDSS >3"),
            (*series(new, adjusted, "Tract-based"), green, "--", "s", "Tract-based — new EDSS ≥4"),
            (*series(old, adjusted, "Classical-region"), blue, "-", "o", "Classical-region — old EDSS >3"),
            (*series(new, adjusted, "Classical-region"), blue, "--", "s", "Classical-region — new EDSS ≥4"),
        ]

        for values, lower_error, upper_error, color, linestyle, marker, label in plot_specs:
            ax.errorbar(
                x,
                values,
                yerr=np.vstack([lower_error, upper_error]),
                color=color,
                linestyle=linestyle,
                linewidth=2.1,
                marker=marker,
                markersize=5.7,
                markerfacecolor="white" if linestyle == "--" else color,
                markeredgecolor=color,
                markeredgewidth=1.4,
                elinewidth=0.9,
                capsize=2.7,
                alpha=0.98,
                label=label,
                zorder=4 if linestyle == "--" else 3,
            )

        ax.axhline(0.5, color=neutral, linestyle=":", linewidth=1.2, alpha=0.85)
        ax.set_xticks(x, DISPLAY)
        ax.set_xlim(-0.28, 6.28)
        ax.set_ylim(0.40, 1.025)
        ax.set_yticks(np.arange(0.4, 1.01, 0.1))
        ax.set_xlabel("MRI metric")
        ax.set_title(panel_title, pad=13)
        ax.grid(axis="y", color=grid, linewidth=0.8, alpha=0.65)
        ax.grid(axis="x", visible=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#9CA3AF")
        ax.spines["bottom"].set_color("#9CA3AF")
        ax.tick_params(colors=text)

        if legend_handles is None:
            legend_handles, legend_labels = ax.get_legend_handles_labels()

    axes[0].set_ylabel("AUC (95% CI)")
    axes[0].text(
        6.0,
        0.512,
        "Chance AUC = 0.50",
        color=neutral,
        ha="right",
        va="bottom",
        fontsize=9.3,
    )

    fig.suptitle(
        "EDSS threshold sensitivity with and without demographic adjustment",
        x=0.5,
        y=0.965,
        fontsize=16,
    )
    fig.text(
        0.5,
        0.906,
        "Old: EDSS >3 (equivalent to ≥3.5 in this cohort); new: EDSS ≥4",
        ha="center",
        va="center",
        fontsize=10.8,
        color=text,
    )
    fig.legend(
        legend_handles,
        legend_labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.068),
        ncol=2,
        frameon=False,
        handlelength=3.1,
        columnspacing=2.0,
    )
    fig.text(
        0.5,
        0.026,
        "High/low groups: old 39/41 (T1/MTR) or 39/40 (other metrics); "
        "new 31/49 or 31/48.",
        ha="center",
        va="bottom",
        fontsize=9.6,
        color="#4B5563",
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    png = output_dir / "edss_threshold_auc_with_without_demographics.png"
    pdf = output_dir / "edss_threshold_auc_with_without_demographics.pdf"
    fig.savefig(png, dpi=240, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(png)
    print(pdf)


if __name__ == "__main__":
    main()
