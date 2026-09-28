#!/usr/bin/env python3
"""Reproduce the complete EDSS >=4 sensitivity package in revisionExtras/work."""

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
REPO_ROOT = PACKAGE_ROOT.parent
DEFAULT_CONFIG = PACKAGE_ROOT / "helper/config/input_paths.local.json"
INPUT_MANIFEST = PACKAGE_ROOT / "helper/config/input_manifest.json"
REFERENCE_RESULTS = PACKAGE_ROOT / "results"
CANONICAL_RESULTS = REPO_ROOT / "revision/results"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--run-dir",
        type=Path,
        help="New output directory inside revisionExtras; defaults to a timestamped work run.",
    )
    return parser.parse_args()


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


def inside_package(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    try:
        resolved.relative_to(PACKAGE_ROOT)
    except ValueError as exc:
        raise ValueError(f"{label} must remain inside revisionExtras: {resolved}") from exc
    return resolved


def matlab_quote(value: Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


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


def load_and_verify_inputs(config_path: Path) -> dict[str, Path]:
    config_path = config_path.expanduser().resolve()
    if not config_path.is_file():
        raise FileNotFoundError(
            f"Missing local input config: {config_path}. Copy input_paths.example.json "
            "to input_paths.local.json and set the local paths."
        )
    config = json.loads(config_path.read_text(encoding="utf-8"))
    required = {"metrics_dir", "stats_dir", "classical_lesion_file", "matlab_executable"}
    missing = sorted(required - set(config))
    if missing:
        raise KeyError(f"Config is missing: {', '.join(missing)}")

    paths = {name: Path(config[name]).expanduser().resolve() for name in required}
    paths["plot_python_executable"] = Path(
        config.get("plot_python_executable", sys.executable)
    ).expanduser().resolve()
    if not paths["metrics_dir"].is_dir() or not paths["stats_dir"].is_dir():
        raise FileNotFoundError("metrics_dir and stats_dir must be existing directories")
    if not paths["classical_lesion_file"].is_file():
        raise FileNotFoundError(paths["classical_lesion_file"])
    if not paths["matlab_executable"].is_file():
        discovered = shutil.which(str(paths["matlab_executable"]))
        if discovered is None:
            raise FileNotFoundError(paths["matlab_executable"])
        paths["matlab_executable"] = Path(discovered).resolve()
    if not paths["plot_python_executable"].is_file():
        discovered = shutil.which(str(paths["plot_python_executable"]))
        if discovered is None:
            raise FileNotFoundError(paths["plot_python_executable"])
        paths["plot_python_executable"] = Path(discovered).resolve()
    try:
        subprocess.run(
            [
                str(paths["plot_python_executable"]),
                "-c",
                "import matplotlib, numpy",
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            "plot_python_executable must provide matplotlib and numpy. "
            "Install helper/config/requirements.txt in that environment."
        ) from exc

    manifest = json.loads(INPUT_MANIFEST.read_text(encoding="utf-8"))
    checked: dict[str, str] = {}
    for item in manifest["files"]:
        location = item["location"]
        if location == "classical_lesion_file":
            path = paths[location]
        else:
            path = paths[location] / item["name"]
        if not path.is_file():
            raise FileNotFoundError(path)
        observed = sha256(path)
        if observed != item["sha256"]:
            raise AssertionError(
                f"Frozen input hash mismatch for {path}: {observed} != {item['sha256']}"
            )
        checked[str(path)] = observed
    paths["clinical_file"] = paths["stats_dir"] / "clinicalScore.xlsx"
    paths["tract_lesion_file"] = paths["stats_dir"] / "GroupTractLesionLoad.xlsx"
    paths["config_path"] = config_path
    paths["checked_hashes"] = checked  # type: ignore[assignment]
    return paths


def main() -> None:
    args = parse_args()
    paths = load_and_verify_inputs(args.config)
    if args.run_dir is None:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        run_dir = PACKAGE_ROOT / "work" / f"run_{stamp}"
    else:
        run_dir = inside_package(args.run_dir, "run-dir")
    if run_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing run directory: {run_dir}")
    run_dir.mkdir(parents=True)

    qmri_dir = run_dir / "qmri"
    lesion_dir = run_dir / "lesion"
    qa_qmri = run_dir / "qa/qmri"
    qa_lesion = run_dir / "qa/lesion"
    generated_results = run_dir / "results"
    manifest_path = run_dir / "qa/package_manifest.json"

    canonical_qmri = CANONICAL_RESULTS / "model_results/all_qmri_manuscript_models.csv"
    canonical_lesion = CANONICAL_RESULTS / "model_results/lesionload_manuscript_models.csv"
    matlab_expression = (
        f"addpath({matlab_quote(CODE_DIR)}); "
        "run_edss_ge4_qmri_models("
        f"'StatsDir',{matlab_quote(paths['stats_dir'])},"
        f"'MetricsDir',{matlab_quote(paths['metrics_dir'])},"
        f"'CanonicalModelFile',{matlab_quote(canonical_qmri)},"
        f"'OutputDir',{matlab_quote(qmri_dir)},"
        f"'QADir',{matlab_quote(qa_qmri)}); "
        "run_edss_ge4_lesion_models("
        f"'ClinicalFile',{matlab_quote(paths['clinical_file'])},"
        f"'TractLesionFile',{matlab_quote(paths['tract_lesion_file'])},"
        f"'ClassicalLesionFile',{matlab_quote(paths['classical_lesion_file'])},"
        f"'CanonicalModelFile',{matlab_quote(canonical_lesion)},"
        f"'OutputDir',{matlab_quote(lesion_dir)},"
        f"'QADir',{matlab_quote(qa_lesion)});"
    )
    run_logged(
        [str(paths["matlab_executable"]), "-batch", matlab_expression],
        run_dir / "logs/matlab.log",
    )

    run_logged(
        [
            sys.executable,
            str(CODE_DIR / "build_edss_ge4_manuscript_results.py"),
            "--canonical-results",
            str(CANONICAL_RESULTS),
            "--qmri-sensitivity",
            str(qmri_dir / "edss_threshold_model_results.csv"),
            "--lesion-sensitivity",
            str(lesion_dir / "lesion_edss_threshold_model_results.csv"),
            "--lesion-baseline-qa",
            str(qa_lesion / "lesion_edss_baseline_reproduction.csv"),
            "--output-dir",
            str(generated_results),
            "--manifest",
            str(manifest_path),
        ],
        run_dir / "logs/build.log",
    )

    figures_dir = run_dir / "figures"
    for script in ["plot_edss_ge4_auc.py", "plot_edss_ge4_auc_adjustment_panels.py"]:
        run_logged(
            [
                str(paths["plot_python_executable"]),
                str(CODE_DIR / script),
                "--canonical-models",
                str(canonical_qmri),
                "--sensitivity-models",
                str(generated_results / "model_results/all_qmri_manuscript_models.csv"),
                "--output-dir",
                str(figures_dir),
            ],
            run_dir / "logs" / f"{Path(script).stem}.log",
        )

    canonical_inventory = tree_hashes(CANONICAL_RESULTS)
    expected_hashes = {
        relative: sha256(REFERENCE_RESULTS / relative)
        for relative in canonical_inventory
    }
    observed_hashes = tree_hashes(generated_results)
    if observed_hashes != expected_hashes:
        differing = sorted(
            key
            for key in set(expected_hashes) | set(observed_hashes)
            if expected_hashes.get(key) != observed_hashes.get(key)
        )
        raise AssertionError(f"Rebuilt package differs from committed results: {differing}")

    report = {
        "schema_version": 1,
        "status": "PASS",
        "run_dir": str(run_dir),
        "config": str(paths["config_path"]),
        "input_hashes": paths["checked_hashes"],
        "canonical_results": str(CANONICAL_RESULTS),
        "canonical_results_unchanged": True,
        "rebuilt_results_match_committed_reference": True,
        "result_file_count": len(observed_hashes),
        "generated_figures": sorted(path.name for path in figures_dir.iterdir()),
    }
    report_path = run_dir / "qa/reproduction_report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
