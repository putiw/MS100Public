#!/usr/bin/env python3
"""Fit AMICO-NODDI for every captured ses-01 TractoFlow source subject.

All writes are confined to this temporary analysis tree. Source derivatives and
clinical files are opened read-only. Each subject is isolated in its own process,
with AMICO and BLAS restricted to one thread to preserve the AMICO 2.1.1 RMSE
buffer safeguard. The first subject bootstraps the shared kernel cache; all later
subjects reuse the verified cache.
"""

from __future__ import annotations

import argparse
import atexit
import contextlib
import csv
import hashlib
import json
import os
import pickle
import platform
import shutil
import sys
import time
import traceback
import uuid
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
from typing import Any

import fcntl


ROOT = Path(
    os.environ.get(
        "MS100_NODDI_WORK_ROOT",
        Path(__file__).resolve().parents[1],
    )
).expanduser().resolve()
CONFIG_PATH = Path(
    os.environ.get(
        "MS100_NODDI_CONFIG",
        ROOT / "config" / "full_cohort_config.json",
    )
).expanduser().resolve()
for variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[variable] = "1"
os.environ["DIPY_HOME"] = str(ROOT / "shared" / "dipy_home")

import nibabel as nib
import numpy as np


EXPECTED_NIFTIS = (
    "fit_NDI.nii.gz",
    "fit_ODI.nii.gz",
    "fit_FWF.nii.gz",
    "fit_dir.nii.gz",
    "fit_RMSE.nii.gz",
    "fit_NRMSE.nii.gz",
)
REQUIRED_FINAL = (
    *EXPECTED_NIFTIS,
    "config.pickle",
    "acquisition.scheme",
    "run.log",
    "run_config.json",
    "subject_metrics.json",
    "qc_slice.npz",
)
HASHED_FINAL = tuple(name for name in REQUIRED_FINAL if name != "run_config.json")


@dataclass(frozen=True)
class Preflight:
    subject_id: str
    dwi: str
    mask: str
    bval: str
    bvec: str
    dwi_shape: tuple[int, int, int, int]
    mask_shape: tuple[int, int, int]
    voxel_size: tuple[float, float, float]
    affine: list[list[float]]
    bvec_norm_min: float
    bvec_norm_max: float
    b0_norm_max: float
    mask_voxels: int
    shell_counts: dict[str, int]
    bvals: list[float]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def safe_timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def json_safe(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, tuple):
        return [json_safe(item) for item in value]
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    return value


def load_config() -> dict[str, Any]:
    with CONFIG_PATH.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    if Path(config["output_root"]).resolve() != ROOT.resolve():
        raise RuntimeError("Configured output_root does not match this temporary run tree")
    source_root = Path(config["source_root"]).resolve()
    if source_root == ROOT.resolve() or ROOT.resolve() in source_root.parents:
        raise RuntimeError("Source and output trees must be distinct")
    if ROOT == Path("/") or ROOT == Path.home():
        raise RuntimeError("Unsafe output root")
    return config


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(json_safe(payload), handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, path)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    if not rows:
        temporary.write_text("", encoding="utf-8")
    else:
        fields: list[str] = []
        for row in rows:
            for key in row:
                if key not in fields:
                    fields.append(key)
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
    os.replace(temporary, path)


def acquire_run_lock() -> Any:
    """Hold an advisory single-run lock until the returned handle is closed."""
    path = ROOT / "logs" / "run.lock"
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = path.open("a+", encoding="utf-8")
    try:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as failure:
        handle.close()
        raise RuntimeError(f"Another full-cohort runner holds {path}") from failure
    handle.seek(0)
    handle.truncate()
    handle.write(
        json.dumps(
            {
                "pid": os.getpid(),
                "host": platform.node(),
                "acquired_utc": utc_now(),
            },
            sort_keys=True,
        )
        + "\n"
    )
    handle.flush()
    return handle


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path: Path, include_sha256: bool = False) -> dict[str, Any]:
    stat = path.stat()
    record: dict[str, Any] = {
        "path": str(path),
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }
    if include_sha256:
        record["sha256"] = sha256_file(path)
    return record


def content_record(path: Path) -> dict[str, Any]:
    """Path-independent record for a file that may be atomically relocated."""
    return {
        "size_bytes": int(path.stat().st_size),
        "sha256": sha256_file(path),
    }


def installed_versions() -> dict[str, Any]:
    packages = (
        "dmri-amico",
        "dmri-dicelib",
        "dipy",
        "numpy",
        "scipy",
        "nibabel",
        "threadpoolctl",
    )
    versions: dict[str, str] = {}
    for package in packages:
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = "NOT INSTALLED"
    return {
        "timestamp_utc": utc_now(),
        "python_version": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "packages": versions,
        "thread_environment": {
            key: os.environ.get(key)
            for key in (
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS",
                "NUMEXPR_NUM_THREADS",
                "DIPY_HOME",
            )
        },
    }


def subject_input_paths(source_root: Path, subject_id: str) -> dict[str, Path]:
    subject_root = source_root / subject_id
    return {
        "dwi": subject_root / "Normalize_DWI" / f"{subject_id}__dwi_normalized.nii.gz",
        "mask": subject_root / "Crop_DWI" / f"{subject_id}__b0_mask_cropped.nii.gz",
        "bval": subject_root / "Eddy_Topup" / f"{subject_id}__bval_eddy",
        "bvec": subject_root / "Eddy_Topup" / f"{subject_id}__dwi_eddy_corrected.bvec",
    }


def preflight_subject(config: dict[str, Any], subject_id: str) -> Preflight:
    paths = subject_input_paths(Path(config["source_root"]), subject_id)
    for label, path in paths.items():
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"{subject_id}: missing or empty {label}: {path}")

    acquisition = config["acquisition"]
    expected_volumes = int(acquisition["expected_volumes"])
    dwi_img = nib.load(paths["dwi"])
    mask_img = nib.load(paths["mask"])
    if len(dwi_img.shape) != 4 or dwi_img.shape[3] != expected_volumes:
        raise RuntimeError(f"{subject_id}: expected 4D/132-volume DWI, got {dwi_img.shape}")
    if len(mask_img.shape) != 3 or tuple(dwi_img.shape[:3]) != tuple(mask_img.shape):
        raise RuntimeError(f"{subject_id}: DWI/mask dimensions differ")
    affine_tolerance = float(acquisition["affine_absolute_tolerance"])
    if not np.allclose(dwi_img.affine, mask_img.affine, rtol=0.0, atol=affine_tolerance):
        raise RuntimeError(f"{subject_id}: DWI/mask affines differ")
    dwi_zooms = tuple(float(value) for value in dwi_img.header.get_zooms()[:3])
    mask_zooms = tuple(float(value) for value in mask_img.header.get_zooms()[:3])
    expected_zooms = tuple(float(value) for value in acquisition["voxel_size_mm"])
    if not np.allclose(dwi_zooms, mask_zooms, rtol=0.0, atol=1e-6):
        raise RuntimeError(f"{subject_id}: DWI/mask voxel sizes differ")
    if not np.allclose(dwi_zooms, expected_zooms, rtol=0.0, atol=1e-6):
        raise RuntimeError(f"{subject_id}: expected 2-mm native grid, got {dwi_zooms}")

    bvals = np.loadtxt(paths["bval"], dtype=float).reshape(-1)
    bvecs = np.loadtxt(paths["bvec"], dtype=float)
    if bvals.shape != (expected_volumes,):
        raise RuntimeError(f"{subject_id}: unexpected b-value shape {bvals.shape}")
    if bvecs.shape != (3, expected_volumes):
        raise RuntimeError(f"{subject_id}: unexpected b-vector shape {bvecs.shape}")
    if not np.all(np.isfinite(bvals)) or not np.all(np.isfinite(bvecs)):
        raise RuntimeError(f"{subject_id}: nonfinite b-values or b-vectors")
    unique, counts = np.unique(bvals, return_counts=True)
    observed = {f"{float(shell):g}": int(count) for shell, count in zip(unique, counts, strict=True)}
    expected = {f"{float(shell):g}": int(count) for shell, count in acquisition["expected_shell_counts"].items()}
    if observed != expected:
        raise RuntimeError(f"{subject_id}: shell counts {observed}, expected {expected}")
    norms = np.linalg.norm(bvecs, axis=0)
    diffusion = bvals > 0
    tolerance = float(acquisition["bvec_unit_tolerance"])
    if np.any(np.abs(norms[diffusion] - 1.0) > tolerance):
        raise RuntimeError(f"{subject_id}: diffusion b-vector outside unit tolerance")
    if np.any(norms[~diffusion] > tolerance):
        raise RuntimeError(f"{subject_id}: b0 b-vector outside zero tolerance")
    mask = np.asarray(mask_img.dataobj) > 0
    if not mask.any():
        raise RuntimeError(f"{subject_id}: empty mask")
    mask_voxels = int(mask.sum())
    minimum = int(config["evaluation"]["rmse_buffer_min_voxels_per_thread"])
    if mask_voxels < minimum:
        raise RuntimeError(f"{subject_id}: mask too small for AMICO RMSE buffer safeguard")
    return Preflight(
        subject_id=subject_id,
        dwi=str(paths["dwi"]),
        mask=str(paths["mask"]),
        bval=str(paths["bval"]),
        bvec=str(paths["bvec"]),
        dwi_shape=tuple(int(value) for value in dwi_img.shape),
        mask_shape=tuple(int(value) for value in mask_img.shape),
        voxel_size=dwi_zooms,
        affine=np.asarray(dwi_img.affine, dtype=float).tolist(),
        bvec_norm_min=float(norms[diffusion].min()),
        bvec_norm_max=float(norms[diffusion].max()),
        b0_norm_max=float(norms[~diffusion].max()),
        mask_voxels=mask_voxels,
        shell_counts=observed,
        bvals=[float(value) for value in bvals],
    )


def robust_stats(values: np.ndarray, bounded: bool = False) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    output: dict[str, float] = {"finite_fraction": float(finite.mean())}
    if not finite.any():
        for key in ("min", "p01", "p25", "median", "mean", "sd", "p75", "p99", "max"):
            output[key] = float("nan")
        return output
    finite_values = values[finite]
    output.update(
        {
            "min": float(np.min(finite_values)),
            "p01": float(np.percentile(finite_values, 1)),
            "p25": float(np.percentile(finite_values, 25)),
            "median": float(np.median(finite_values)),
            "mean": float(np.mean(finite_values)),
            "sd": float(np.std(finite_values, ddof=1)),
            "p75": float(np.percentile(finite_values, 75)),
            "p99": float(np.percentile(finite_values, 99)),
            "max": float(np.max(finite_values)),
        }
    )
    if bounded:
        output.update(
            {
                "lower_bound_fraction": float(np.mean(np.isclose(finite_values, 0.0, atol=1e-6))),
                "upper_bound_fraction": float(np.mean(np.isclose(finite_values, 1.0, atol=1e-6))),
                "outside_0_1_fraction": float(
                    np.mean((finite_values < -1e-6) | (finite_values > 1.0 + 1e-6))
                ),
            }
        )
    return output


def configure_evaluation(ae: Any, config: dict[str, Any], atoms_path: Path) -> dict[str, Any]:
    model = config["model"]
    evaluation = config["evaluation"]
    ae.set_model(model["name"])
    ae.set_config("ATOMS_path", str(atoms_path))
    ae.model.set(
        dPar=float(model["dPar_mm2_per_s"]),
        dIso=float(model["dIso_mm2_per_s"]),
        IC_VFs=np.asarray(model["IC_VFs"], dtype=float),
        IC_ODs=np.asarray(model["IC_ODs"], dtype=float),
        isExvivo=bool(model["isExvivo"]),
    )
    ae.set_solver(
        lambda1=float(model["solver_lambda1"]),
        lambda2=float(model["solver_lambda2"]),
    )
    for key in (
        "doNormalizeSignal",
        "doKeepb0Intact",
        "doComputeRMSE",
        "doComputeNRMSE",
        "doSaveModulatedMaps",
        "doMergeB0",
        "doDebiasSignal",
        "doDirectionalAverage",
        "DTI_fit_method",
        "nthreads",
        "BLAS_nthreads",
    ):
        ae.set_config(key, evaluation[key])
    if int(evaluation["nthreads"]) != 1 or int(evaluation["BLAS_nthreads"]) != 1:
        raise RuntimeError("Full-cohort RMSE run is locked to one AMICO and one BLAS thread")
    ae.set_config(
        "rmse_buffer_min_voxels_per_thread",
        int(evaluation["rmse_buffer_min_voxels_per_thread"]),
    )
    ae.set_config("rmse_buffer_safety_note", evaluation["rmse_buffer_safety_note"])
    params = json_safe(ae.model.get_params())
    if not np.isclose(params["dPar"], model["dPar_mm2_per_s"], rtol=0.0, atol=0.0):
        raise RuntimeError("AMICO did not retain dPar")
    if not np.isclose(params["dIso"], model["dIso_mm2_per_s"], rtol=0.0, atol=0.0):
        raise RuntimeError("AMICO did not retain dIso")
    return params


def expected_kernel_names() -> list[str]:
    return [f"A_{index:03d}.npy" for index in range(1, 146)]


def verify_kernel_cache(atoms_path: Path) -> None:
    observed = sorted(path.name for path in atoms_path.glob("A_*.npy"))
    expected = expected_kernel_names()
    if observed != expected:
        raise RuntimeError(f"Kernel cache has {len(observed)} atoms; expected 145")
    if any((atoms_path / name).stat().st_size == 0 for name in expected):
        raise RuntimeError("Kernel cache contains an empty atom")


def kernel_manifest(atoms_path: Path) -> dict[str, dict[str, Any]]:
    verify_kernel_cache(atoms_path)
    return {
        name: {
            "size_bytes": int((atoms_path / name).stat().st_size),
            "sha256": sha256_file(atoms_path / name),
        }
        for name in expected_kernel_names()
    }


def ensure_auxiliary_matrices(config: dict[str, Any], dipy_home: Path) -> None:
    from amico import lut

    lmax = int(config["model"]["lmax"])
    ndirs = int(config["model"]["ndirs"])
    expected = dipy_home / f"AMICO_aux_matrices_lmax={lmax}_ndirs={ndirs}.pickle"
    if not expected.is_file():
        lut.precompute_rotation_matrices(lmax=lmax, ndirs=ndirs)
    if not expected.is_file() or expected.stat().st_size == 0:
        raise RuntimeError(f"AMICO auxiliary matrix missing: {expected}")


def auxiliary_matrix_path(config: dict[str, Any], dipy_home: Path) -> Path:
    return dipy_home / (
        f"AMICO_aux_matrices_lmax={int(config['model']['lmax'])}_"
        f"ndirs={int(config['model']['ndirs'])}.pickle"
    )


def auxiliary_matrix_manifest(config: dict[str, Any], dipy_home: Path) -> dict[str, Any]:
    path = auxiliary_matrix_path(config, dipy_home)
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"AMICO auxiliary matrix missing: {path}")
    return {
        "name": path.name,
        "size_bytes": int(path.stat().st_size),
        "sha256": sha256_file(path),
    }


def extract_metrics(
    output_dir: Path,
    preflight: Preflight,
    mean_b0: np.ndarray,
) -> dict[str, Any]:
    mask_img = nib.load(preflight.mask)
    mask = np.asarray(mask_img.dataobj) > 0
    affine = np.asarray(preflight.affine, dtype=float)
    arrays: dict[str, np.ndarray] = {}
    for name in EXPECTED_NIFTIS:
        path = output_dir / name
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"{preflight.subject_id}: missing output {name}")
        image = nib.load(path)
        expected_shape = (*preflight.mask_shape, 3) if name == "fit_dir.nii.gz" else preflight.mask_shape
        if tuple(image.shape) != tuple(expected_shape):
            raise RuntimeError(f"{preflight.subject_id}: {name} shape {image.shape}, expected {expected_shape}")
        if not np.allclose(image.affine, affine, rtol=0.0, atol=1e-5):
            raise RuntimeError(f"{preflight.subject_id}: {name} differs from native DWI geometry")
        arrays[name] = np.asarray(image.dataobj, dtype=np.float32)

    map_stats = {
        "NDI": robust_stats(arrays["fit_NDI.nii.gz"][mask], bounded=True),
        "ODI": robust_stats(arrays["fit_ODI.nii.gz"][mask], bounded=True),
        "FWF": robust_stats(arrays["fit_FWF.nii.gz"][mask], bounded=True),
        "RMSE": robust_stats(arrays["fit_RMSE.nii.gz"][mask]),
        "NRMSE": robust_stats(arrays["fit_NRMSE.nii.gz"][mask]),
    }
    ndi = np.asarray(arrays["fit_NDI.nii.gz"][mask], dtype=np.float64)
    odi = np.asarray(arrays["fit_ODI.nii.gz"][mask], dtype=np.float64)
    fwf = np.asarray(arrays["fit_FWF.nii.gz"][mask], dtype=np.float64)
    tissue_weights = 1.0 - fwf
    tissue_denominator = float(np.sum(tissue_weights))
    if not np.isfinite(tissue_denominator) or tissue_denominator <= 0:
        raise RuntimeError(f"{preflight.subject_id}: invalid tissue-fraction denominator")
    ndi_weighted = float(np.sum(ndi * tissue_weights) / tissue_denominator)
    odi_weighted = float(np.sum(odi * tissue_weights) / tissue_denominator)

    directions = arrays["fit_dir.nii.gz"][mask, :]
    direction_finite = np.all(np.isfinite(directions), axis=1)
    norms = np.linalg.norm(directions[direction_finite], axis=1)
    norm_stats = robust_stats(norms)
    unit_fraction = float(np.mean(np.isclose(norms, 1.0, atol=1e-3))) if norms.size else 0.0

    failures: list[str] = []
    warnings: list[str] = []
    for name, stats in map_stats.items():
        if stats["finite_fraction"] != 1.0:
            failures.append(f"{name}_nonfinite")
    for name in ("NDI", "ODI", "FWF"):
        if map_stats[name]["outside_0_1_fraction"] > 0:
            failures.append(f"{name}_outside_0_1")
        if map_stats[name]["lower_bound_fraction"] > 0.10:
            warnings.append(f"{name}_lower_bound_gt_10pct")
        if map_stats[name]["upper_bound_fraction"] > 0.10:
            warnings.append(f"{name}_upper_bound_gt_10pct")
    for name in ("RMSE", "NRMSE"):
        if map_stats[name]["min"] < -1e-7:
            failures.append(f"{name}_negative")
    if float(direction_finite.mean()) != 1.0:
        failures.append("direction_nonfinite")
    if unit_fraction < 0.999:
        failures.append("direction_not_unit_length")
    if not 0.0 <= ndi_weighted <= 1.0:
        failures.append("NDI_tissue_weighted_outside_0_1")
    if not 0.0 <= odi_weighted <= 1.0:
        failures.append("ODI_tissue_weighted_outside_0_1")

    row: dict[str, Any] = {
        "subject_id": preflight.subject_id,
        "fit_status": "SUCCESS",
        "technical_qc": "PASS" if not failures else "FAIL",
        "qc_failures": ";".join(failures) if failures else "none",
        "qc_warnings": ";".join(warnings) if warnings else "none",
        "mask_voxels": preflight.mask_voxels,
        "mask_volume_ml": float(preflight.mask_voxels * np.prod(preflight.voxel_size) / 1000.0),
        "NDI_tissue_weighted_mean": ndi_weighted,
        "ODI_tissue_weighted_mean": odi_weighted,
        "FWF_median": map_stats["FWF"]["median"],
        "FWF_mean": map_stats["FWF"]["mean"],
        "mean_tissue_fraction": float(np.mean(tissue_weights)),
        "tissue_weight_sum": tissue_denominator,
        "direction_finite_fraction": float(direction_finite.mean()),
        "direction_norm_p01": norm_stats["p01"],
        "direction_norm_median": norm_stats["median"],
        "direction_norm_p99": norm_stats["p99"],
        "direction_unit_fraction": unit_fraction,
        "NDI_grid_upper_0_99_fraction": float(np.mean(np.isclose(ndi, 0.99, atol=1e-6))),
        "ODI_grid_lower_0_03_fraction": float(np.mean(np.isclose(odi, 0.03, atol=1e-6))),
        "ODI_grid_upper_0_99_fraction": float(np.mean(np.isclose(odi, 0.99, atol=1e-6))),
    }
    for map_name, stats in map_stats.items():
        for statistic, value in stats.items():
            row[f"{map_name}_{statistic}"] = value

    mask_area = mask.sum(axis=(0, 1))
    z_index = int(np.flatnonzero(mask_area == mask_area.max())[0])
    np.savez_compressed(
        output_dir / "qc_slice.npz",
        subject_id=np.asarray(preflight.subject_id),
        native_axial_slice=np.asarray(z_index, dtype=np.int16),
        mask=mask[..., z_index].astype(np.uint8),
        mean_b0=np.asarray(mean_b0[..., z_index], dtype=np.float32),
        NDI=np.asarray(arrays["fit_NDI.nii.gz"][..., z_index], dtype=np.float32),
        ODI=np.asarray(arrays["fit_ODI.nii.gz"][..., z_index], dtype=np.float32),
        FWF=np.asarray(arrays["fit_FWF.nii.gz"][..., z_index], dtype=np.float32),
        NRMSE=np.asarray(arrays["fit_NRMSE.nii.gz"][..., z_index], dtype=np.float32),
    )
    return row


def completed_subject(
    output_dir: Path,
    config: dict[str, Any],
    preflight: Preflight,
) -> bool:
    if not output_dir.is_dir():
        return False
    if any(not (output_dir / name).is_file() or (output_dir / name).stat().st_size == 0 for name in REQUIRED_FINAL):
        return False
    try:
        with (output_dir / "subject_metrics.json").open("r", encoding="utf-8") as handle:
            metrics = json.load(handle)
        with (output_dir / "run_config.json").open("r", encoding="utf-8") as handle:
            run_config = json.load(handle)
        if metrics.get("fit_status") != "SUCCESS":
            return False
        if run_config.get("analysis_config_sha256") != sha256_file(CONFIG_PATH):
            return False
        for label, source in (
            ("dwi", Path(preflight.dwi)),
            ("mask", Path(preflight.mask)),
            ("bval", Path(preflight.bval)),
            ("bvec", Path(preflight.bvec)),
        ):
            recorded = run_config["inputs"][label]
            current = file_record(source, include_sha256=True)
            if recorded != current:
                return False
        for name in HASHED_FINAL:
            output = output_dir / name
            if run_config["outputs"][name] != content_record(output):
                return False
        return bool(config["execution"]["resume_completed_subjects"])
    except (OSError, ValueError, KeyError, TypeError):
        return False


def run_subject_worker(
    config: dict[str, Any],
    preflight_payload: dict[str, Any],
    permit_kernel_generation: bool,
) -> dict[str, Any]:
    os.environ["DIPY_HOME"] = str(ROOT / "shared" / "dipy_home")
    for variable in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[variable] = "1"
    preflight = Preflight(**preflight_payload)
    subject_id = preflight.subject_id
    final_dir = ROOT / "derivatives" / config["session"] / subject_id
    if completed_subject(final_dir, config, preflight):
        with (final_dir / "subject_metrics.json").open("r", encoding="utf-8") as handle:
            metrics = json.load(handle)
        return {"subject_id": subject_id, "status": "SKIPPED_COMPLETE", "metrics": metrics}
    if final_dir.exists():
        preserved = ROOT / "failures" / (
            f"{subject_id}.incomplete_final.{safe_timestamp()}.{uuid.uuid4().hex[:8]}"
        )
        preserved.parent.mkdir(parents=True, exist_ok=True)
        os.replace(final_dir, preserved)

    attempt = ROOT / "work" / f"{subject_id}.{uuid.uuid4().hex}.inprogress"
    output_dir = attempt / "output"
    attempt.mkdir(parents=True, exist_ok=False)
    scheme_path = attempt / f"{subject_id}.scheme"
    log_path = attempt / "run.log"
    atoms_path = ROOT / "shared" / "kernels" / "NODDI"
    started = time.time()
    fit_seconds = float("nan")
    direction_seconds = float("nan")
    model_params: dict[str, Any] = {}
    explicit_config: dict[str, Any] = {}
    try:
        import amico
        from amico.util import fsl2scheme

        if amico.__version__ != "2.1.1":
            raise RuntimeError(f"Expected dmri-amico 2.1.1, got {amico.__version__}")
        if sys.version_info[:2] != (3, 12):
            raise RuntimeError(f"Expected Python 3.12, got {platform.python_version()}")
        amico.set_verbose(2)
        fsl2scheme(
            preflight.bval,
            preflight.bvec,
            str(scheme_path),
            bStep=config["acquisition"]["scheme_b_step"],
        )
        if not scheme_path.is_file() or scheme_path.stat().st_size == 0:
            raise RuntimeError(f"{subject_id}: scheme creation failed")

        with log_path.open("w", encoding="utf-8") as log_handle:
            with contextlib.redirect_stdout(log_handle), contextlib.redirect_stderr(log_handle):
                print(f"MS100 AMICO-NODDI full-cohort subject: {subject_id}")
                print(f"Started UTC: {utc_now()}")
                print(f"Read-only source: {preflight.dwi}")
                print(f"Temporary output: {output_dir}")
                print(json.dumps(installed_versions(), indent=2, sort_keys=True))
                ae = amico.Evaluation(
                    study_path=str(ROOT),
                    subject=subject_id,
                    output_path=str(output_dir),
                )
                model_params = configure_evaluation(ae, config, atoms_path)
                ae.load_data(
                    dwi_filename=preflight.dwi,
                    scheme_filename=str(scheme_path),
                    mask_filename=preflight.mask,
                    b0_thr=float(config["acquisition"]["b0_threshold"]),
                    b0_min_signal=float(config["evaluation"]["b0_min_signal"]),
                    replace_bad_voxels=config["evaluation"]["replace_bad_voxels"],
                )
                mean_b0 = np.asarray(ae.mean_b0s, dtype=np.float32).copy()
                existing = sorted(atoms_path.glob("A_*.npy"))
                if existing:
                    verify_kernel_cache(atoms_path)
                    ae.generate_kernels(
                        regenerate=False,
                        lmax=int(config["model"]["lmax"]),
                        ndirs=int(config["model"]["ndirs"]),
                    )
                elif permit_kernel_generation:
                    ae.generate_kernels(
                        regenerate=True,
                        lmax=int(config["model"]["lmax"]),
                        ndirs=int(config["model"]["ndirs"]),
                    )
                    verify_kernel_cache(atoms_path)
                else:
                    raise RuntimeError("Shared kernel cache is absent outside the bootstrap subject")
                ae.load_kernels()
                ae.fit()
                ae.save_results()
                fit_seconds = float(ae.get_config("fit_time"))
                direction_seconds = float(ae.get_config("dirs_precomputing_time"))
                explicit_config = json_safe(ae.CONFIG)
                del ae
                print(f"Fit completed UTC: {utc_now()}")

        shutil.copy2(scheme_path, output_dir / "acquisition.scheme")
        shutil.copy2(log_path, output_dir / "run.log")
        metrics = extract_metrics(output_dir, preflight, mean_b0)
        metrics["elapsed_seconds"] = float(time.time() - started)
        metrics["fit_seconds"] = fit_seconds
        metrics["direction_estimation_seconds"] = direction_seconds
        write_json(output_dir / "subject_metrics.json", metrics)
        output_records = {
            name: content_record(output_dir / name)
            for name in HASHED_FINAL
        }
        write_json(
            output_dir / "run_config.json",
            {
                "status": "SUCCESS",
                "subject_id": subject_id,
                "analysis_config_sha256": sha256_file(CONFIG_PATH),
                "started_and_completed_in_seconds": float(time.time() - started),
                "fit_seconds": fit_seconds,
                "direction_estimation_seconds": direction_seconds,
                "environment": installed_versions(),
                "inputs": {
                    "dwi": file_record(Path(preflight.dwi), include_sha256=True),
                    "mask": file_record(Path(preflight.mask), include_sha256=True),
                    "bval": file_record(Path(preflight.bval), include_sha256=True),
                    "bvec": file_record(Path(preflight.bvec), include_sha256=True),
                },
                "outputs": output_records,
                "preflight": preflight_payload,
                "model_parameters": model_params,
                "explicit_amico_config": explicit_config,
                "shared_kernel_cache": {
                    "path": str(atoms_path),
                    "atom_count": 145,
                    "generated_if_absent_by_this_subject": permit_kernel_generation,
                },
                "prohibited_operations": {
                    "registration_to_T1": False,
                    "resampling": False,
                    "voxelwise_group_averaging": False,
                },
            },
        )
        for name in REQUIRED_FINAL:
            if not (output_dir / name).is_file() or (output_dir / name).stat().st_size == 0:
                raise RuntimeError(f"{subject_id}: final validation missing {name}")
        final_dir.parent.mkdir(parents=True, exist_ok=True)
        os.replace(output_dir, final_dir)
        shutil.rmtree(attempt)
        return {"subject_id": subject_id, "status": "SUCCESS", "metrics": metrics}
    except BaseException as failure:
        failure_payload = {
            "status": "FAILED",
            "subject_id": subject_id,
            "failed_at_utc": utc_now(),
            "exception_type": type(failure).__name__,
            "exception": str(failure),
            "traceback": traceback.format_exc(),
        }
        write_json(attempt / "failure.json", failure_payload)
        failed_dir = ROOT / "failures" / f"{subject_id}.{safe_timestamp()}.{uuid.uuid4().hex[:8]}"
        failed_dir.parent.mkdir(parents=True, exist_ok=True)
        os.replace(attempt, failed_dir)
        return {
            "subject_id": subject_id,
            "status": "FAILED",
            "failure_dir": str(failed_dir),
            "error": str(failure),
        }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument(
        "--bootstrap-only",
        action="store_true",
        help="Fit/verify only the serialized bootstrap subject, then stop cleanly.",
    )
    parser.add_argument("--max-workers", type=int, default=None)
    args = parser.parse_args()
    config = load_config()
    run_lock = acquire_run_lock()
    atexit.register(run_lock.close)
    workers = args.max_workers or int(config["execution"]["default_subject_processes"])
    if workers < 1 or workers > 4:
        raise RuntimeError("max-workers must be between 1 and 4 for this isolated run")

    source_root = Path(config["source_root"])
    subjects = sorted(path.name for path in source_root.glob("sub-*") if path.is_dir())
    expected_subjects_path = Path(config["expected_source_subjects_file"])
    expected_subjects = sorted(
        line.strip()
        for line in expected_subjects_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    )
    if subjects != expected_subjects:
        missing = sorted(set(expected_subjects) - set(subjects))
        unexpected = sorted(set(subjects) - set(expected_subjects))
        raise RuntimeError(
            f"Captured source cohort changed; missing={missing}, unexpected={unexpected}"
        )
    if config["execution"]["bootstrap_subject"] not in subjects:
        raise RuntimeError("Configured bootstrap subject is absent")
    print(f"Captured {len(subjects)} source subjects; strict preflight started", flush=True)

    preflights: dict[str, Preflight] = {}
    preflight_failures: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=4) as executor:
        future_map = {
            executor.submit(preflight_subject, config, subject_id): subject_id for subject_id in subjects
        }
        for index, future in enumerate(as_completed(future_map), start=1):
            subject_id = future_map[future]
            try:
                preflights[subject_id] = future.result()
            except BaseException as failure:
                preflight_failures.append({"subject_id": subject_id, "error": str(failure)})
            if index % 20 == 0 or index == len(subjects):
                print(
                    f"PREFLIGHT {index}/{len(subjects)} checked; failures={len(preflight_failures)}",
                    flush=True,
                )
    preflight_rows = [asdict(preflights[subject]) for subject in sorted(preflights)]
    write_json(ROOT / "logs" / "preflight.json", {"rows": preflight_rows, "failures": preflight_failures})
    write_csv(
        ROOT / "logs" / "source_manifest.csv",
        [
            {
                "subject_id": subject,
                "source_status": "preflight_pass" if subject in preflights else "preflight_fail",
                "preflight_error": next(
                    (row["error"] for row in preflight_failures if row["subject_id"] == subject),
                    "",
                ),
            }
            for subject in subjects
        ],
    )
    if preflight_failures:
        print(json.dumps(preflight_failures, indent=2), flush=True)
        raise RuntimeError("Strict preflight failed; no subject fitting started")
    reference_bvals = np.asarray(preflights[config["execution"]["bootstrap_subject"]].bvals)
    bval_mismatches = [
        subject
        for subject in subjects
        if not np.array_equal(np.asarray(preflights[subject].bvals), reference_bvals)
    ]
    if bval_mismatches:
        raise RuntimeError(
            "Ordered b-value vectors differ from the bootstrap subject: "
            + ", ".join(bval_mismatches)
        )
    print("Strict preflight passed for all 132 source subjects", flush=True)
    if args.preflight_only:
        return 0

    dipy_home = ROOT / "shared" / "dipy_home"
    dipy_home.mkdir(parents=True, exist_ok=True)
    atoms_path = ROOT / "shared" / "kernels" / "NODDI"
    atoms_path.mkdir(parents=True, exist_ok=True)
    existing_atoms = sorted(atoms_path.glob("A_*.npy"))
    if existing_atoms:
        previous_paths = (
            ROOT / "logs" / "kernel_manifest_final.json",
            ROOT / "logs" / "kernel_manifest_after_bootstrap.json",
        )
        previous_path = next((path for path in previous_paths if path.is_file()), None)
        reusable = False
        if [path.name for path in existing_atoms] == expected_kernel_names() and previous_path:
            try:
                previous = json.loads(previous_path.read_text(encoding="utf-8"))
                reusable = kernel_manifest(atoms_path) == previous
            except (OSError, ValueError, KeyError, TypeError):
                reusable = False
        if not reusable:
            preserved = ROOT / "failures" / (
                f"unverified_kernel_cache.{safe_timestamp()}.{uuid.uuid4().hex[:8]}"
            )
            preserved.parent.mkdir(parents=True, exist_ok=True)
            os.replace(atoms_path, preserved)
            atoms_path.mkdir(parents=True, exist_ok=True)
    import amico

    if amico.__version__ != "2.1.1":
        raise RuntimeError(f"Expected dmri-amico 2.1.1, got {amico.__version__}")
    auxiliary_path = auxiliary_matrix_path(config, dipy_home)
    auxiliary_manifest_path = ROOT / "logs" / "auxiliary_matrix_manifest.json"
    if auxiliary_path.is_file():
        reusable_auxiliary = False
        if auxiliary_manifest_path.is_file():
            try:
                previous_auxiliary = json.loads(
                    auxiliary_manifest_path.read_text(encoding="utf-8")
                )
                reusable_auxiliary = (
                    auxiliary_matrix_manifest(config, dipy_home) == previous_auxiliary
                )
            except (OSError, ValueError, KeyError, TypeError):
                reusable_auxiliary = False
        if not reusable_auxiliary:
            preserved = ROOT / "failures" / (
                f"unverified_auxiliary_matrix.{safe_timestamp()}.{uuid.uuid4().hex[:8]}.pickle"
            )
            preserved.parent.mkdir(parents=True, exist_ok=True)
            os.replace(auxiliary_path, preserved)
    ensure_auxiliary_matrices(config, dipy_home)
    baseline_auxiliary_manifest = auxiliary_matrix_manifest(config, dipy_home)
    write_json(auxiliary_manifest_path, baseline_auxiliary_manifest)
    bootstrap = config["execution"]["bootstrap_subject"]
    with ProcessPoolExecutor(max_workers=1) as bootstrap_executor:
        bootstrap_result = bootstrap_executor.submit(
            run_subject_worker,
            config,
            asdict(preflights[bootstrap]),
            True,
        ).result()
    if bootstrap_result["status"] == "FAILED":
        write_json(ROOT / "logs" / "progress.json", {"results": [bootstrap_result]})
        raise RuntimeError(f"Bootstrap subject failed: {bootstrap_result['error']}")
    baseline_kernel_manifest = kernel_manifest(atoms_path)
    write_json(ROOT / "logs" / "kernel_manifest_after_bootstrap.json", baseline_kernel_manifest)
    results: list[dict[str, Any]] = [bootstrap_result]
    print(f"BOOTSTRAP {bootstrap}: {bootstrap_result['status']}; shared cache verified", flush=True)
    if args.bootstrap_only:
        print("Bootstrap-only validation complete", flush=True)
        return 0

    remaining = [subject for subject in subjects if subject != bootstrap]
    completed = 1
    failed = 0
    with ProcessPoolExecutor(max_workers=workers) as executor:
        future_map = {
            executor.submit(run_subject_worker, config, asdict(preflights[subject]), False): subject
            for subject in remaining
        }
        for future in as_completed(future_map):
            subject_id = future_map[future]
            try:
                result = future.result()
            except BaseException as failure:
                result = {"subject_id": subject_id, "status": "FAILED", "error": str(failure)}
            results.append(result)
            completed += 1
            if result["status"] == "FAILED":
                failed += 1
            write_json(
                ROOT / "logs" / "progress.json",
                {
                    "updated_utc": utc_now(),
                    "subjects_total": len(subjects),
                    "subjects_returned": completed,
                    "failed": failed,
                    "results": sorted(results, key=lambda row: row["subject_id"]),
                },
            )
            print(
                f"PROGRESS {completed}/{len(subjects)}; {subject_id}={result['status']}; failed={failed}",
                flush=True,
            )

    final_kernel_manifest = kernel_manifest(atoms_path)
    if final_kernel_manifest != baseline_kernel_manifest:
        raise RuntimeError("Shared kernel cache changed during read-only reuse")
    write_json(ROOT / "logs" / "kernel_manifest_final.json", final_kernel_manifest)
    if auxiliary_matrix_manifest(config, dipy_home) != baseline_auxiliary_manifest:
        raise RuntimeError("AMICO auxiliary matrix changed during the cohort run")
    status_rows = [
        {
            "subject_id": result["subject_id"],
            "fit_status": result["status"],
            "technical_qc": result.get("metrics", {}).get("technical_qc", ""),
            "error": result.get("error", ""),
            "failure_dir": result.get("failure_dir", ""),
        }
        for result in sorted(results, key=lambda row: row["subject_id"])
    ]
    write_csv(ROOT / "logs" / "run_status.csv", status_rows)
    successful = sum(row["fit_status"] in ("SUCCESS", "SKIPPED_COMPLETE") for row in status_rows)
    technical_failures = sum(row["technical_qc"] == "FAIL" for row in status_rows)
    print(
        "FULL COHORT FIT RETURNED: "
        f"success={successful}, failed={failed}, technical_qc_fail={technical_failures}, "
        f"total={len(subjects)}",
        flush=True,
    )
    run_lock.close()
    return 0 if failed == 0 and technical_failures == 0 and successful == len(subjects) else 1


if __name__ == "__main__":
    raise SystemExit(main())
