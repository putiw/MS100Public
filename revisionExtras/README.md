# Reviewer sensitivity analyses

This directory contains the separate analyses for Reviewer Comments 9 and 14. It is separate from `revision/`; the canonical analysis and its results are read-only inputs and are not replaced.

## Comment 9: EDSS ≥4 sensitivity

The canonical manuscript defines the higher-EDSS group as `EDSS > 3`, which is equivalent to `EDSS ≥3.5` in this cohort. This package changes only that definition to `EDSS ≥4` while retaining each model's canonical subject set, predictors, preprocessing, regularization, bootstrap procedure, and covariate specification.

The aggregate `results/` directory mirrors the layout of `revision/results/`. Exactly six files contain threshold-dependent changes:

- `model_results/all_qmri_manuscript_models.csv`: 28 EDSS rows covering seven qMRI metrics, tract/classical models, and unadjusted/adjusted specifications.
- `model_results/lesionload_manuscript_models.csv`: eight unadjusted EDSS lesion-load rows.
- `tables/Main_Table_2.csv`: seven unadjusted qMRI EDSS rows.
- `tables/Supplementary_Table_3.csv`: seven demographic-adjusted qMRI EDSS rows.
- `tables/Supplementary_Table_4.csv`: 14 demographic-adjusted qMRI EDSS significance rows.
- `tables/Supplementary_Table_5.csv`: eight unadjusted lesion-load EDSS rows.

All other canonical files in `results/` are byte-identical to `revision/results/` and are included to make direct folder-level comparison easy.

## Comment 14: whole-brain qMRI measures

`results/` also contains two new files for the requested whole-brain T1, MTR, FA, and MD analysis:

- `model_results/whole_brain_qmri_models.csv`: the 20 imaging-only whole-brain models used in the table.
- `tables/Supplementary_Table_6.csv`: four whole-brain MRI metrics across the original five clinical outcomes: EDSS, MSPro, T25FW, and dominant and non-dominant 9HPT.

Supplementary Table 6 reports the whole-brain parenchymal mean as the only imaging predictor. No age, sex, disease-duration, or other covariates enter these models. AUC with 95% CI is reported for EDSS and MSPro, and `R2` is reported for continuous outcomes. This matches the imaging-only specification in Supplementary Table 5. The p-values are raw exploratory values and are not added to the manuscript's `m=98` multiplicity family.

This table uses the `revisionExtras` definition `EDSS ≥4`. All controls are excluded. EDSS and MSPro use the exact metric-specific joint tract/classical participant cohorts from Supplementary Table 3. T25FW and both 9HPT outcomes use the analogous matched MS-only cohorts already used for the manuscript's matched tract-versus-classical comparisons. This keeps the whole-brain qMRI results directly comparable to the existing qMRI results rather than allowing extra whole-brain-only participants into the models.

Patient-level source data, whole-brain subject values, and subject-membership QA tables are not included in this public package.

## Table interpretation

- Main Table 2 reports unadjusted qMRI AUC, 95% CI, and AIC for imaging-only tract-based and classical-region models.
- Supplementary Table 3 reports AUC, 95% CI, and AIC for models that combine imaging with age, sex, and disease duration. These combined-model AUCs can exceed Main Table 2 because the demographic variables add predictors; they should not be interpreted as improved imaging-only performance.
- Supplementary Table 4 reports raw and manuscript-wide Bonferroni-adjusted p-values for the adjusted qMRI models. The correction family remains `m=98`.
- Supplementary Table 5 reports unadjusted whole-brain, tract-based, and classical-region lesion-load models.
- Supplementary Table 6 reports imaging-only whole-brain T1, MTR, FA, and MD models for the original five clinical outcomes requested for Comment 14; controls, demographic covariates, SDMT, and MSFC are omitted.

Under EDSS ≥4, qMRI models have 31 positive cases (`N=80` for T1/MTR and `N=79` for the other metrics). Six of the 14 adjusted EDSS qMRI models remain significant after the manuscript-wide correction, compared with eight under the canonical threshold. Across the full adjusted qMRI family, 44 of 98 models are significant, compared with 46 originally.

The lesion-load sensitivity deliberately uses the frozen canonical Supplementary Table 5 cohort (`N=89`, 32 positive cases). RR047 is retained in this one threshold-only comparison because the canonical lesion results included that participant; removing RR047 would combine a threshold change with a cohort change. Seven of eight lesion EDSS models have raw `p<0.05`, compared with eight of eight under the canonical threshold.

## Structure

```text
revisionExtras/
  code/                  Portable MATLAB and Python analysis scripts
  helper/config/         Input-path template, dependency pins, and input hashes
  helper/reference_qa/   Aggregate reference checks and comparison figures
  results/               EDSS ≥4 mirror plus Comment 14 model/table additions
  work/                  Regenerated intermediates and QA; ignored by Git
```

The analysis code is intentionally limited to the two reviewer analyses. It does not duplicate the unrelated NODDI fitting/extraction pipeline in `revision/`.

## Requirements

- MATLAB R2024b with Statistics and Machine Learning Toolbox.
- Python 3.10 or later.
- Python packages in `helper/config/requirements.txt` for the comparison figures.
- The frozen local aggregate inputs listed in `helper/config/input_manifest.json` and `helper/config/comment14_input_manifest.json`.

The raw participant-level inputs cannot be distributed publicly. Each runner verifies its frozen inputs against SHA-256 fingerprints before fitting any model.

## Reproduce the package

1. Copy `helper/config/input_paths.example.json` to `helper/config/input_paths.local.json`.
2. Set the local aggregate-input directories, the frozen whole-brain means file, MATLAB executable, and a Python executable with `matplotlib` and `numpy` in the local config.
3. Install the plotting dependencies if needed:

   ```bash
   python -m pip install -r helper/config/requirements.txt
   ```

4. From `revisionExtras/`, run:

   ```bash
   python code/run_all_edss_ge4.py --config helper/config/input_paths.local.json
   ```

   To reproduce the Comment 14 models and Supplementary Table 6, run:

   ```bash
   python code/run_comment14_whole_brain.py --config helper/config/input_paths.local.json
   ```

The EDSS runner performs the following steps:

1. Verifies every frozen input against `input_manifest.json`.
2. Reproduces all 28 canonical qMRI EDSS rows, then refits them using EDSS ≥4.
3. Reproduces all eight canonical lesion-load EDSS rows, then refits them using EDSS ≥4.
4. Builds a clean manuscript-format mirror under a timestamped `work/run_*/results/` directory.
5. Applies the original manuscript-wide Bonferroni family (`m=98`) to adjusted qMRI models.
6. Regenerates both EDSS AUC comparison figures.
7. Requires the rebuilt 13-file canonical-format result set to match the corresponding files in `results/` exactly.

The Comment 14 runner verifies the frozen clinical, whole-brain, tract, and classical-region inputs; derives the manuscript-matched MS-only cohorts; refits the 20 imaging-only models; rebuilds Supplementary Table 6; and requires both generated CSVs to match the packaged references exactly.

All generated files remain under `revisionExtras/work/`. Both runners refuse to write outside `revisionExtras`; the EDSS package builder also protects the canonical `revision/` directory with before/after content hashes.

## Analysis details

The qMRI and lesion analyses preserve the manuscript settings: ridge logistic regression, 50 lambdas over `10^-6` to `10^6`, 10-fold cross-validation for lambda selection, full-sample refitting, RNG seed 42, and 2,000 bootstrap resamples. The canonical manuscript result fields use the legacy unstratified bootstrap mean and percentile interval. The separate paired tract-versus-classical comparison uses identical class-stratified resamples for the two model scores.

The Comment 14 analysis uses the same binary ridge settings. Continuous whole-brain models use ridge linear regression with 100 lambdas over `10^-6` to `10^6`; T25FW and both 9HPT outcomes are log-transformed. The whole-brain measure is the native-space parenchymal mean defined from the FreeSurfer `aparc+aseg` mask; no maps were resampled.
