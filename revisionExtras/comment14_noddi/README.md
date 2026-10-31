# Comment 14 whole-brain NODDI

Adds whole-brain NDI, ODI and ISOVF to Supplementary Table 6 for R2 Comment 14.
Means use FreeSurfer aparc+aseg parenchymal labels in individual T1 space and
AMICO NODDI maps transformed by the existing tract pipeline. Lesions are retained.
Finite values in [0,1], including fitted zeros, are averaged within parenchyma;
voxels outside transformed AMICO fitting support are excluded.
Models use the same metric/outcome-specific matched tract/classical MS cohorts
as the existing Supplementary Table 6, with no demographic covariates.

Requires Python with numpy, pandas, nibabel and openpyxl, plus MATLAB R2024b
with Statistics and Machine Learning Toolbox. Copy input_paths.example.json to
a private local config; cohort_file is the 132-row imaging cohort with SubjectID.
Run from this folder:
```bash
python run_comment14_noddi.py --config /path/to/input_paths.local.json --output-dir /path/to/new/results
```

The two combined seven-metric CSVs are in results/. Private intermediates and logs
are written to work_dir outside the public repository. Existing outputs are preserved.
