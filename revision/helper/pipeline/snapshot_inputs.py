#!/usr/bin/env python3
"""Create a local, read-only-source snapshot for tract/classical NODDI."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


GROUPS = ("Association", "Cerebellar", "Occipitoparietal", "ProjectionBrainstem")
HEMIS = ("L", "R", "")
EXPECTED_GROUP_MASKS_PER_INCLUDED_SUBJECT = len(GROUPS) * len(HEMIS)
ANALYSIS_CONFIG = Path(__file__).resolve().parents[1] / "config" / "analysis_config.json"


@dataclass(frozen=True)
class CopyItem:
    source: Path
    destination: Path
    role: str
    required: bool = True


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def copy_one(item: CopyItem) -> dict[str, object]:
    src = item.source
    dst = item.destination
    result: dict[str, object] = {
        "role": item.role,
        "required": item.required,
        "source": str(src),
        "destination": str(dst),
        "status": "",
        "size_bytes": "",
        "source_mtime_ns": "",
        "source_sha256": "",
        "destination_sha256": "",
        "error": "",
    }
    if not src.is_file():
        if item.required:
            result["status"] = "MISSING_REQUIRED"
            return result

        # A resumed work root must reflect the current source inventory.  In
        # particular, an optional MTR map or grouped mask that has disappeared
        # upstream must not survive locally from an earlier snapshot and be
        # mistaken for a present input.  Symlinks are removed as links (also
        # covering broken symlinks); an unexpected directory fails closed.
        try:
            if dst.is_symlink() or dst.is_file():
                dst.unlink()
                result["status"] = "REMOVED_STALE_OPTIONAL"
            elif dst.exists():
                raise IsADirectoryError(
                    f"optional snapshot destination is not a file: {dst}"
                )
            else:
                result["status"] = "MISSING_OPTIONAL"
        except Exception as exc:
            result["status"] = "COPY_ERROR"
            result["error"] = repr(exc)
        return result

    try:
        src_stat = src.stat()
        source_digest = sha256(src)
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.is_file():
            destination_digest = sha256(dst)
            if destination_digest == source_digest:
                status = "REUSED"
            else:
                temporary = dst.with_name(f".{dst.name}.{uuid.uuid4().hex}.partial")
                shutil.copy2(src, temporary)
                os.replace(temporary, dst)
                status = "UPDATED"
        else:
            temporary = dst.with_name(f".{dst.name}.{uuid.uuid4().hex}.partial")
            shutil.copy2(src, temporary)
            os.replace(temporary, dst)
            status = "COPIED"
        dst_stat = dst.stat()
        if dst_stat.st_size != src_stat.st_size:
            raise IOError(f"size mismatch: source={src_stat.st_size}, destination={dst_stat.st_size}")
        destination_digest = sha256(dst)
        if destination_digest != source_digest:
            raise IOError("SHA-256 mismatch after local snapshot copy")
        result.update(
            status=status,
            size_bytes=dst_stat.st_size,
            source_mtime_ns=src_stat.st_mtime_ns,
            source_sha256=source_digest,
            destination_sha256=destination_digest,
        )
    except Exception as exc:  # fail closed and preserve the manifest entry
        result["status"] = "COPY_ERROR"
        result["error"] = repr(exc)
    return result


def add_file(items: list[CopyItem], bids: Path, snapshot: Path, relative: str, role: str, required: bool = True) -> None:
    rel = Path(relative)
    items.append(CopyItem(bids / rel, snapshot / rel, role, required))


def build_items(
    bids: Path,
    noddi_root: Path,
    snapshot: Path,
    group_mask_overlay: Path | None = None,
) -> tuple[list[str], list[CopyItem]]:
    subjects = sorted(p.name for p in (noddi_root / "derivatives" / "ses-01").glob("sub-*") if p.is_dir())
    if not subjects:
        raise RuntimeError(f"No NODDI subject directories found below {noddi_root}")

    items: list[CopyItem] = []
    for subject in subjects:
        base = f"derivatives/TractoFlow/ses-01/{subject}"
        add_file(items, bids, snapshot, f"{base}/Resample_T1/{subject}__t1_resampled.nii.gz", "t1_reference")
        add_file(items, bids, snapshot, f"{base}/Register_T1/{subject}__output0GenericAffine.mat", "inverse_chain_affine")
        add_file(items, bids, snapshot, f"{base}/Register_T1/{subject}__output1InverseWarp.nii.gz", "inverse_chain_warp")
        add_file(items, bids, snapshot, f"{base}/DTI_Metrics/{subject}__fa.nii.gz", "transform_validation_fa")
        add_file(items, bids, snapshot, f"{base}/Crop_DWI/{subject}__b0_mask_cropped.nii.gz", "noddi_valid_support")
        # Snapshot every qMRI map used by the manuscript so the MATLAB table
        # generator can rebuild the complete seven-metric analysis locally.
        # The source BIDS tree remains read-only throughout the workflow.
        # Two prespecified source participants have no MTR map; their original
        # manuscript rows are blank, so that one optional input is preserved as
        # missing instead of turning a known cohort-level absence into a failure.
        legacy_maps = (
            (
                f"{subject}_ses-01_space-individual_T1map.nii.gz",
                "legacy_t1_map",
                True,
            ),
            (
                f"{subject}_ses-01_space-individual_MTRmap.nii.gz",
                "legacy_mtr_map",
                False,
            ),
            (f"{subject}_space-individual_FA.nii.gz", "legacy_fa_reference", True),
            (f"{subject}_space-individual_MD.nii.gz", "legacy_md_reference", True),
        )
        for filename, role, required in legacy_maps:
            add_file(
                items,
                bids,
                snapshot,
                f"derivatives/maps/{subject}/{filename}",
                role,
                required=required,
            )

        # The classical-region analysis reuses the original manuscript masks.
        # They are required for every MS participant in the 87-row classical
        # cohort; control copies are retained when present but are not analyzed.
        is_ms = not subject.startswith("sub-C")
        add_file(
            items,
            bids,
            snapshot,
            f"derivatives/freesurfer/{subject}/mri/aparc+aseg.mgz",
            "classical_aparc_aseg",
            required=is_ms,
        )
        add_file(
            items,
            bids,
            snapshot,
            f"derivatives/freesurfer/{subject}/mri/ribbon.mgz",
            "classical_ribbon",
            required=is_ms,
        )

        for group in GROUPS:
            for hemi in HEMIS:
                relative = Path(
                    f"derivatives/TractoFlow_post/{subject}/groupTract/{group}{hemi}.nii.gz"
                )
                source = bids / relative
                if group_mask_overlay is not None and not source.is_file():
                    recovered = group_mask_overlay / subject / "groupTract" / relative.name
                    if recovered.is_file():
                        source = recovered
                items.append(
                    CopyItem(
                        source,
                        snapshot / relative,
                        "grouped_tract_mask",
                        required=False,
                    )
                )

        lesion_rel = f"derivatives/lesionMask/{subject}/ses-01/{subject}_ses-01_desc-lesionManual_mask.nii.gz"
        add_file(items, bids, snapshot, lesion_rel, "manual_lesion_mask", required=not subject.startswith("sub-C"))

    stats_files = (
        "clinicalScore.xlsx",
        "GroupTractT1_All.xlsx",
        "GroupTractMTR_All.xlsx",
        "GroupTractFA_All.xlsx",
        "GroupTractMD_All.xlsx",
        "icometrixLesionStats.xlsx",
    )
    for filename in stats_files:
        add_file(
            items,
            bids,
            snapshot,
            f"derivatives/derivativesStats/{filename}",
            "analysis_table",
        )
    return subjects, items


def write_manifest(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "role",
        "required",
        "source",
        "destination",
        "status",
        "size_bytes",
        "source_mtime_ns",
        "source_sha256",
        "destination_sha256",
        "error",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def validate_group_mask_profiles(
    rows: list[dict[str, object]], subjects: list[str]
) -> tuple[int, int, list[str]]:
    study = json.loads(ANALYSIS_CONFIG.read_text(encoding="utf-8"))["study"]
    expected_full_subjects = int(study["expected_imaging_subjects"])
    expected_tract_subjects = int(study["expected_primary_tract_subjects"])
    expected_mask_manifest_hash = str(
        study["expected_group_mask_manifest_sha256"]
    )
    present_statuses = {"COPIED", "UPDATED", "REUSED"}
    counts = {subject: 0 for subject in subjects}
    for row in rows:
        if row["role"] != "grouped_tract_mask":
            continue
        source = Path(str(row["source"]))
        subject = source.parent.parent.name
        if subject not in counts:
            return 0, 0, [f"Unrecognized subject directory in mask snapshot: {subject}"]
        if row["status"] in present_statuses:
            counts[subject] += 1

    complete = sorted(
        subject
        for subject, count in counts.items()
        if count == EXPECTED_GROUP_MASKS_PER_INCLUDED_SUBJECT
    )
    absent = sorted(subject for subject, count in counts.items() if count == 0)
    partial = sorted(
        (subject, count)
        for subject, count in counts.items()
        if count not in {0, EXPECTED_GROUP_MASKS_PER_INCLUDED_SUBJECT}
    )
    errors: list[str] = []
    if len(subjects) != expected_full_subjects:
        errors.append(
            f"Expected {expected_full_subjects} source subjects for grouped-mask profiling, "
            f"found {len(subjects)}"
        )
    if partial:
        errors.append(
            "Partial grouped-mask profiles: "
            + ", ".join(f"{subject}={count}/12" for subject, count in partial)
        )
    if len(complete) != expected_tract_subjects:
        errors.append(
            f"Expected {expected_tract_subjects} complete grouped-mask profiles, "
            f"found {len(complete)}"
        )
    expected_absent = int(study["expected_missing_original_group_mask_subjects"])
    if len(absent) != expected_absent:
        errors.append(
            f"Expected {expected_absent} all-absent grouped-mask profiles, found {len(absent)}"
        )
    manifest_lines: list[str] = []
    for row in rows:
        if row["role"] != "grouped_tract_mask" or row["status"] not in present_statuses:
            continue
        source = Path(str(row["source"]))
        digest = str(row.get("source_sha256") or row.get("destination_sha256") or "")
        if not digest:
            errors.append(f"Missing hash for copied grouped mask: {source.name}")
            continue
        manifest_lines.append(f"{source.parent.parent.name}/{source.name}\t{digest}")
    observed_mask_manifest_hash = hashlib.sha256(
        ("\n".join(sorted(manifest_lines)) + "\n").encode("utf-8")
    ).hexdigest()
    if observed_mask_manifest_hash != expected_mask_manifest_hash:
        errors.append(
            "Grouped-mask files do not match the frozen content-bound manifest fingerprint"
        )
    return len(complete), len(absent), errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bids-root", type=Path, required=True)
    parser.add_argument("--noddi-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument(
        "--group-mask-overlay",
        type=Path,
        help=(
            "Optional read-only directory containing sub-*/groupTract/*.nii.gz; "
            "used only when the corresponding source-tree mask is absent."
        ),
    )
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()

    bids = args.bids_root.resolve()
    noddi_root = args.noddi_root.resolve()
    output_root = args.output_root.resolve()
    group_mask_overlay = (
        args.group_mask_overlay.resolve() if args.group_mask_overlay else None
    )
    snapshot = output_root / "source_snapshot"
    if not bids.is_dir():
        print(f"ERROR: source BIDS root is not mounted: {bids}", file=sys.stderr)
        return 2
    if group_mask_overlay is not None and not group_mask_overlay.is_dir():
        print(
            f"ERROR: grouped-mask overlay is not a directory: {group_mask_overlay}",
            file=sys.stderr,
        )
        return 2
    if str(snapshot).startswith(str(bids) + os.sep):
        print("ERROR: snapshot destination resolves inside the read-only source tree", file=sys.stderr)
        return 2

    subjects, items = build_items(
        bids, noddi_root, snapshot, group_mask_overlay=group_mask_overlay
    )
    print(f"Snapshotting {len(items)} files for {len(subjects)} subjects with {args.workers} workers", flush=True)
    rows: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = {pool.submit(copy_one, item): item for item in items}
        for index, future in enumerate(as_completed(futures), start=1):
            rows.append(future.result())
            if index % 100 == 0 or index == len(items):
                print(f"  completed {index}/{len(items)}", flush=True)

    rows.sort(key=lambda row: (str(row["role"]), str(row["source"])))
    manifest = output_root / "logs" / "source_snapshot_manifest.csv"
    write_manifest(manifest, rows)
    failures = [row for row in rows if row["status"] in {"MISSING_REQUIRED", "COPY_ERROR"}]
    complete_mask_profiles, absent_mask_profiles, mask_profile_errors = (
        validate_group_mask_profiles(rows, subjects)
    )
    summary_path = output_root / "logs" / "source_snapshot_summary.txt"
    summary = [
        f"timestamp_utc={datetime.now(timezone.utc).isoformat()}",
        f"subjects={len(subjects)}",
        f"files_planned={len(rows)}",
        f"files_present={sum(row['status'] in {'COPIED', 'UPDATED', 'REUSED'} for row in rows)}",
        f"required_failures={len(failures)}",
        f"complete_group_mask_profiles={complete_mask_profiles}",
        f"all_absent_group_mask_profiles={absent_mask_profiles}",
        f"group_mask_profile_failures={len(mask_profile_errors)}",
        f"bytes_present={sum(int(row['size_bytes']) for row in rows if row['size_bytes'] != '')}",
        f"manifest={manifest}",
    ]
    summary_path.write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary), flush=True)
    if failures or mask_profile_errors:
        for row in failures[:30]:
            print(f"FAIL {row['status']} {row['source']} {row['error']}", file=sys.stderr)
        for error in mask_profile_errors:
            print(f"FAIL {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
