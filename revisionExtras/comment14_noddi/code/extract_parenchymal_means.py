"""Whole-brain NODDI means from individual-space maps and FreeSurfer labels."""
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import nibabel as nib
import numpy as np
import pandas as pd

METRICS = ('NDI', 'ODI', 'FWF')
PARENCHYMAL_LABELS = (
    2, 7, 8, 10, 11, 12, 13, 16, 17, 18, 26, 28, 41, 46, 47, 49, 50,
    51, 52, 53, 54, 58, 60, 77, 80, 85, 251, 252, 253, 254, 255,
    *range(1001, 1004), *range(1005, 1036),
    *range(2001, 2004), *range(2005, 2036),
)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def extract_subject(job):
    cfg, subject = job
    seg_path = Path(cfg['freesurfer_root']) / subject / 'mri/aparc+aseg.mgz'
    seg_img = nib.load(seg_path)
    seg = np.asanyarray(seg_img.dataobj).astype(np.int32, copy=False)
    parenchyma = np.isin(seg, PARENCHYMAL_LABELS)
    maps = Path(cfg['transformed_maps']) / subject
    support_path = maps / f'{subject}_space-individual_NODDI-support.nii.gz'
    support_img = nib.load(support_path)
    assert support_img.shape == seg_img.shape
    assert np.allclose(support_img.affine, seg_img.affine, atol=1e-4, rtol=0)
    support_data = np.asanyarray(support_img.dataobj)
    assert np.all(np.isfinite(support_data))
    support = support_data > 0
    assert parenchyma.any() and support.any()
    result = {'SubjectID': subject}
    qc, inputs = [], [seg_path, support_path]
    for metric in METRICS:
        path = maps / f'{subject}_space-individual_{metric}.nii.gz'
        img = nib.load(path)
        assert img.shape == seg_img.shape
        assert np.allclose(img.affine, seg_img.affine, atol=1e-4, rtol=0)
        data = np.asanyarray(img.dataobj)
        valid = parenchyma & support & np.isfinite(data) & (data >= 0) & (data <= 1)
        assert valid.any()
        result[metric] = float(np.mean(data[valid], dtype=np.float64))
        qc.append({'SubjectID': subject, 'Metric': metric,
                   'ValidVoxels': int(valid.sum()),
                   'FittedZerosRetained': int((valid & (data == 0)).sum())})
        inputs.append(path)
    provenance = [{'path': str(p), 'sha256': sha256(p)} for p in inputs]
    return result, qc, provenance


def extract(cfg, work_dir):
    """Write participant-level means and provenance only to private work_dir."""
    subjects = pd.read_csv(cfg['cohort_file'])['SubjectID'].astype(str).tolist()
    assert len(subjects) == len(set(subjects)) == 132
    work_dir = Path(work_dir)
    with ProcessPoolExecutor(max_workers=int(cfg.get('workers', 2))) as pool:
        records = []
        for i, record in enumerate(pool.map(extract_subject, [(cfg, s) for s in subjects]), 1):
            records.append(record)
            if i % 10 == 0 or i == len(subjects):
                print(f'Extracted {i}/{len(subjects)} participants', flush=True)
    means = work_dir / 'whole_brain_parenchymal_means_noddi.csv'
    table = pd.DataFrame([r[0] for r in records])
    assert np.isfinite(table[list(METRICS)].to_numpy()).all()
    table.to_csv(means, index=False, float_format='%.17g')
    means.chmod(0o600)
    pd.DataFrame([row for record in records for row in record[1]]).to_csv(
        work_dir / 'extraction_qc_private.csv', index=False)
    manifest = {'extractor_sha256': sha256(__file__),
                'cohort_file_sha256': sha256(cfg['cohort_file']),
                'mean_file_sha256': sha256(means),
                'inputs': [item for record in records for item in record[2]]}
    (work_dir / 'extraction_manifest_private.json').write_text(
        json.dumps(manifest, indent=2) + '\n')
    return means
