#!/usr/bin/env python3
"""Reproduce Comment 14 whole-brain models and Supplementary Table 6."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = PACKAGE_ROOT / "code"
DEFAULT_CONFIG = PACKAGE_ROOT / "helper/config/input_paths.local.json"
INPUT_MANIFEST = PACKAGE_ROOT / "helper/config/comment14_input_manifest.json"
REFERENCE_MODELS = PACKAGE_ROOT / "results/model_results/whole_brain_qmri_models.csv"
REFERENCE_TABLE = PACKAGE_ROOT / "results/tables/Supplementary_Table_6.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--run-dir",
        type=Path,
        help="New output directory inside revisionExtras; defaults to work/run_comment14_<timestamp>.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def inside_package(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    try:
        resolved.relative_to(PACKAGE_ROOT)
    except ValueError as exc:
        raise ValueError(f"{label} must remain inside revisionExtras: {resolved}") from exc
    return resolved


def matlab_quote(path: Path) -> str:
    return "'" + str(path).replace("'", "''") + "'"


def run_logged(command: list[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="", flush=True)
            log.write(line)
        return_code = process.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, command)


def resolve_executable(value: str) -> Path:
    path = Path(value).expanduser()
    if path.is_file():
        return path.resolve()
    discovered = shutil.which(value)
    if discovered is None:
        raise FileNotFoundError(value)
    return Path(discovered).resolve()


def main() -> None:
    args = parse_args()
    config_path = args.config.expanduser().resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    required = {
        "stats_dir",
        "metrics_dir",
        "whole_brain_means_file",
        "matlab_executable",
    }
    missing = sorted(required - set(config))
    if missing:
        raise KeyError(f"Config is missing: {', '.join(missing)}")
    stats_dir = Path(config["stats_dir"]).expanduser().resolve()
    metrics_dir = Path(config["metrics_dir"]).expanduser().resolve()
    clinical_file = stats_dir / "clinicalScore.xlsx"
    classical_file = Path(
        config.get("classical_file", metrics_dir / "ClassicalRegionAllMetrics.csv")
    ).expanduser().resolve()
    whole_brain_file = Path(config["whole_brain_means_file"]).expanduser().resolve()
    matlab = resolve_executable(config["matlab_executable"])

    manifest = json.loads(INPUT_MANIFEST.read_text(encoding="utf-8"))
    input_paths = {
        "clinical_file": clinical_file,
        "whole_brain_means_file": whole_brain_file,
        "classical_file": classical_file,
        "tract_t1_file": metrics_dir / "GroupTractT1_All.xlsx",
        "tract_mtr_file": metrics_dir / "GroupTractMTR_All.xlsx",
        "tract_fa_file": metrics_dir / "GroupTractFA_All.xlsx",
        "tract_md_file": metrics_dir / "GroupTractMD_All.xlsx",
    }
    checked_hashes: dict[str, str] = {}
    for item in manifest["files"]:
        path = input_paths[item["location"]]
        if not path.is_file():
            raise FileNotFoundError(path)
        observed = sha256(path)
        if observed != item["sha256"]:
            raise AssertionError(
                f"Frozen input hash mismatch for {path}: {observed} != {item['sha256']}"
            )
        checked_hashes[str(path)] = observed

    if args.run_dir is None:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        run_dir = PACKAGE_ROOT / "work" / f"run_comment14_{stamp}"
    else:
        run_dir = inside_package(args.run_dir, "run-dir")
    if run_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing run directory: {run_dir}")
    model_dir = run_dir / "results/model_results"
    table_dir = run_dir / "results/tables"
    qa_dir = run_dir / "qa"
    model_dir.mkdir(parents=True)
    table_dir.mkdir(parents=True)

    matlab_expression = (
        f"addpath({matlab_quote(CODE_DIR)}); "
        "run_comment14_whole_brain_models("
        f"'ClinicalFile',{matlab_quote(clinical_file)},"
        f"'WholeBrainFile',{matlab_quote(whole_brain_file)},"
        f"'MetricsDir',{matlab_quote(metrics_dir)},"
        f"'ClassicalFile',{matlab_quote(classical_file)},"
        f"'OutputDir',{matlab_quote(model_dir)},"
        f"'QADir',{matlab_quote(qa_dir)});"
    )
    run_logged(
        [str(matlab), "-batch", matlab_expression],
        run_dir / "logs/matlab.log",
    )

    generated_models = model_dir / "whole_brain_qmri_models.csv"
    generated_table = table_dir / "Supplementary_Table_6.csv"
    run_logged(
        [
            sys.executable,
            str(CODE_DIR / "build_comment14_supplementary_table.py"),
            "--models",
            str(generated_models),
            "--output",
            str(generated_table),
        ],
        run_dir / "logs/build_table.log",
    )

    comparisons = {
        "whole_brain_qmri_models.csv": (generated_models, REFERENCE_MODELS),
        "Supplementary_Table_6.csv": (generated_table, REFERENCE_TABLE),
    }
    output_hashes: dict[str, str] = {}
    for label, (generated, reference) in comparisons.items():
        if not reference.is_file():
            raise FileNotFoundError(f"Packaged reference missing: {reference}")
        generated_hash = sha256(generated)
        reference_hash = sha256(reference)
        if generated_hash != reference_hash:
            raise AssertionError(f"Rebuilt {label} differs from packaged reference")
        output_hashes[label] = generated_hash

    report = {
        "schema_version": 1,
        "status": "PASS",
        "analysis": "Comment 14 imaging-only whole-brain T1/MTR/FA/MD",
        "edss_rule": "EDSS >= 4.0 (revisionExtras sensitivity definition)",
        "cohort_policy": (
            "Joint tract/classical complete-case MS cohorts; no controls."
        ),
        "covariates": "none",
        "run_dir": str(run_dir),
        "input_hashes": checked_hashes,
        "output_hashes": output_hashes,
        "rebuilt_outputs_match_packaged_references": True,
    }
    report_path = qa_dir / "reproduction_report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
