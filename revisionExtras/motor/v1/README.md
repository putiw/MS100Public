# Motor analysis: completed-run package

This directory contains the 13-tract Motor analysis requested by the reviewer.
`config/motor.txt` is byte-identical to the supplied tract definition. The
earlier files in the parent `motor/` directory are retained unchanged; this
version writes its own aggregate results under `v1/results/`.

## Reproduce

Create `helper/config/input_paths.local.json` from the example, pointing to
the read-only source BIDS tree, the frozen local map and aggregate snapshots,
MATLAB R2024b, and the local scratch root. Run one command:

The Python environment needs NumPy, SciPy, NiBabel, and openpyxl. MATLAB needs
Statistics and Machine Learning Toolbox.

```bash
python code/run_all_motor.py \
  --config helper/config/input_paths.local.json \
  --stage-dir /path/to/local/scratch/motor_stage \
  --run-dir /path/to/local/scratch/motor_run
```

The runner copies the specified individual masks, hemisphere masks, T1 maps,
and MTR maps from the source to local scratch before processing them. It reuses
the previously validated local FA/MD and NODDI maps, fitting support, and
manual lesion masks. Only the staged local copies are used during mask
construction and extraction. It refuses an existing output directory and
never overwrites staged files of a different size.

The extraction cohort is the 132-subject frozen manuscript imaging cohort in
the standard tract tables. Its overlap with the clinical sheet is 131;
clinical modeling uses that overlap, with endpoint-specific complete cases.
The two imaging subjects with no nonempty Motor component also lack MTR maps;
the report's metric-map missing counts are conditional on a usable Motor mask.

To verify a clean rerun against the packaged outputs, choose a fresh
`--run-dir`, reuse the completed `--stage-dir`, and add
`--reference-dir results`. `--workers 3` may be used to parallelize the
independent per-subject extraction; output order is held fixed.

The new Motor mask is the voxelwise union of available constituents.
Missing tracts are recorded and skipped; a missing tract does not create a
zero-valued metric. The original subject-specific 3×3×3 hemisphere rule splits
the union, with ties assigned left. Geometry, nonzero voxels, partition
identity, and anatomical left/right centroids are checked per subject.

T1 and MD use P90; MTR, FA, and NDI use P10; ODI and FWF use P90.
FWF is displayed as ISOVF. NODDI zeros inside fitted support are valid.
Legacy missing-value behavior follows the manuscript. All extraction is on
the original individual-space grid; no new transformation is applied.

Each metric has a Motor-only model with four predictors and a Standard+Motor
model with the original 16 plus four predictors. The two models share the
same participants within each metric/outcome/adjustment cell. Binary models
use the manuscript joint tract/classical MS cohort; continuous models retain
controls. Each model uses the original MATLAB ridge grid and full-sample
score. EDSS and MSPro AUCs use separate 2,000-resample bootstraps seeded 42.
The result contains 196 rows. AIC and R² follow the existing manuscript
second-stage models. No paired AUC contrast or other new test is run.

Only aggregate CSV and JSON results belong in this package. Subject-level
masks, maps, extracted values, and subject QC stay in local scratch.
