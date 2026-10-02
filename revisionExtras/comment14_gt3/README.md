# Comment 14 EDSS >3 sensitivity package

This is a parallel, additive package for Comment 14. It preserves the existing
EDSS >=4 files and changes only the EDSS rule to `EDSS > 3`. The models are
whole-brain imaging-only T1, MTR, FA, and MD models with the validated matched
manuscript cohorts, no controls, and no demographic covariates.

Run from this directory with:

```bash
python code/run_all_comment14_gt3.py --config helper/config/input_paths.example.json
```

The runner writes only below `work/`, then compares regenerated CSVs with the
packaged references. No paired comparisons or multiplicity family are added.

