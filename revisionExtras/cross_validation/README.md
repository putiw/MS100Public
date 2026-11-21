# Nested cross-validation of tract-based vs classical models

This folder adds an out-of-sample check to the manuscript models. It does not
change any reported in-sample (apparent) result.

## Why

The manuscript models are evaluated in-sample, and the two frameworks differ in
size (16 tract-based predictors vs 4 classical). In-sample AUC and R2 rise with
the number of predictors even when there is no signal, so the between-framework
difference could in principle reflect model flexibility. This analysis tests that
in two ways: nested cross-validation of both frameworks, and a tract-based model
restricted to four predictors, the same number as the classical model.

## Contents

| Path | What it is |
| --- | --- |
| `code/run_nested_cv.m` | MATLAB implementation (fitclinear / fitrlinear). Main entry point. |
| `code/summarize_nested_cv.m` | Prints headline counts and ranges; compares MATLAB output with the Python reference. |
| `code/python_reference/` | Independent Python implementation (numpy and pandas only) used to produce the reference results. |
| `results/reference_python/cv_results_all.csv` | Aggregate reference results, 30 repeats, all seven metrics. |
| `results/nested_cv_results_matlab.csv` | Written by `run_nested_cv.m` (not present until it is run). |

FWF in file and column names is the isotropic volume fraction (ISOVF in the manuscript).

## Design

- **Cohort:** MS patients only, joint complete cases for the tract-based and
  classical predictors of each metric (the manuscript's matched cohort).
- **Outcomes:** EDSS >= 4.0 and MSPro >= 2 (ridge logistic); log T25FW,
  log 9HPT-D, log 9HPT-ND, SDMT correct and MSFC-SDMT (ridge linear).
- **Models:** `tract16` and `classical4` are the manuscript models.
  `tract4_tail` uses one bilateral tail percentile per cluster (four predictors).
  `tract8_tail` uses the left and right tail percentiles (eight predictors).
- **Outer loop:** 30 repeats of 10-fold cross-validation, stratified for binary
  outcomes, with the same folds for every model of an outcome.
- **Inner loop:** 5-fold on the training part only; lambda chosen from 25 values
  between 1e-4 and 1e4 by held-out deviance (binary) or squared error (continuous).
  Predictors are z-scored using the training fold only.
- **Measures:** cross-validated AUC (pooled out-of-fold scores, averaged over
  repeats) and cross-validated R2 (1 - SSE / SSE of the training-fold mean).
- **Paired difference:** out-of-fold scores averaged over repeats, then a
  participant bootstrap of the difference (class-stratified for binary outcomes).

Two deliberate differences from the manuscript pipeline: lambda is chosen by
deviance, not by the default misclassification rate of `kfoldLoss`, and the
lambda grid is narrower (1e-4 to 1e4 instead of 1e-6 to 1e6).

## Run

MATLAB:

```matlab
addpath('revisionExtras/cross_validation/code');
run_nested_cv('statsDir','/path/to/derivatives/stats','metrics',{'T1'},'reps',2);  % quick test
run_nested_cv('statsDir','/path/to/derivatives/stats');                            % full run
summarize_nested_cv();
```

Python reference:

```bash
cd revisionExtras/cross_validation/code/python_reference
export MS100_STATS_DIR=/path/to/derivatives/stats
export MS100_CV_OUT=/path/to/private/work
for m in T1 MTR FA MD NDI ODI FWF; do python3 run.py $m 30; done
python3 combine.py
```

## Data

Participant-level tables are not included and must not be committed. Only the
aggregate result files in `results/` are public.

## Limits

- Imaging-only models; the covariate-adjusted models are not cross-validated here.
- The whole-brain comparator is not included.
- Bootstrap intervals are computed on out-of-fold predictions and do not capture
  all refitting variability; treat them as approximate.
