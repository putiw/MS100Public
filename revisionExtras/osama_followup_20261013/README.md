# Osama follow-up: current handoff

This is the **single current package** for Osama's 6 October 2026 email.
The already completed Motor and EDSS >3 whole-brain analyses remain in
`../motor/v1/` and `../comment14_gt3/`; he reported that both reproduced.
Earlier follow-up drafts are not in this directory.

## Aggregate results to share

| Request | File |
| --- | --- |
| Paired Motor ΔAUC | [`results/paired_motor_delta_auc.csv`](results/paired_motor_delta_auc.csv) |
| Lesion load by region × subtype (LN, LV, Lnorm) | [`results/lesion_load_by_region_subtype_caption_matched_89.csv`](results/lesion_load_by_region_subtype_caption_matched_89.csv) |
| Figure 4A bar values by subtype (mean, SD, mL) | [`results/figure4a_bar_by_subtype_caption_matched_89.csv`](results/figure4a_bar_by_subtype_caption_matched_89.csv) |
| Classical-model predictor exclusions by metric | [`results/classical_exclusions_by_metric.csv`](results/classical_exclusions_by_metric.csv) |
| Figure 4A/4B and Figure 5 post-hoc group tests | [`results/posthoc_fig4_fig5_caption_matched_89.csv`](results/posthoc_fig4_fig5_caption_matched_89.csv) |
| Figure 4 candidate using the matched 89-person decision, with original styling | [`results/Figure4_matched89_original_style.pdf`](results/Figure4_matched89_original_style.pdf) |

The **one per-subject predictor CSV** is deliberately outside the public Git
repository, with the shared copy at
`/Volumes/projects/MS_osama/MS100/MsBIDS/derivatives/derivativesStats/per_subject_predictors_mounted_figures.csv`.
It has 133 participants, 284 columns, subtype, tract tail/NAWM and classical
regional means for all seven metrics, and Figure 4/5 source values. It is on
the project share and is not included in this public repository.

The Figure 4 PDF is a revision candidate implementing the decision to use
the same 89 imaging participants in both frameworks. It preserves the
original artwork's layout, typography, palette, schematic, and axes while
updating RRMS data marks. It does not replace the original manuscript figure.
The styling archive and rebuilding script remain outside Git at
`/Users/pw1246/Documents/Claude/Projects/MS100/figure4_stylematched_20261014/`.
Numerical bar/polar plots are reproducible with the tracked
[`../code/fig4_plain_matched89.m`](../code/fig4_plain_matched89.m).

## Interpretation and caveats

- The paired file contains 84 comparisons: seven metrics × two binary
  outcomes × adjusted/unadjusted × three requested model pairs. It calls
  `paired_stratified_auc` with 2,000 seeded, class-stratified draws using the
  same sampled participants for both fitted-score AUCs. Models are **not**
  refitted inside bootstrap draws; these are apparent, not held-out AUCs.
- Figure 4 uses the same 89 imaging participants in both frameworks (RRMS 49,
  SPMS 23, PPMS 17). The 72-row region × subtype table describes Figure 4B.
  Figure 4A whole-brain lesion volume equals the sum of four classical-region
  volumes; its tract bar sums four tract-region volumes. Both are in mL, and
  overlapping tract masks may contribute to more than one region. The
  six-row bar table reports mean and SD. The post-hoc table has 270 pairwise
  rows: six for Figure 4A, 72 for Figure 4B, and 192 for Figure 5. Tests are
  exploratory Kruskal-Wallis and two-sided Mann-Whitney comparisons; Holm is
  applied within each plotted feature and BH across each figure and all
  pairwise tests. Figure 5 tests use the mounted plotting-source workbooks.
  Adding Figure 4A to the test family changes cross-feature BH values, not
  the 264 previously computed raw comparisons or their within-feature Holm
  values.
- The classical-exclusion table concerns a separate 89-person **clinical** MS
  cohort. Five lack lesions in **at least one** of four classical regions;
  the dedicated `ExcludedLackingLesionsInAllFourClassicalRegionsN` column is
  **zero for every metric**. Other missing classical qMRI predictors
  account for four more patients for T1/MTR and five for FA/MD/NDI/ODI/FWF.
  These are predictor-availability counts, not outcome-specific model Ns.
- The local TractoFlow derivatives are consistent with tensor-derived maps
  computed on a 1-mm resampled diffusion-coordinate grid and later mapped
  to individual T1 space for sampling. This agrees with TractoFlow's
  [processing steps](https://tractoflow-documentation.readthedocs.io/en/latest/pipeline/steps.html)
  and [1-mm default DWI-resampling option](https://tractoflow-documentation.readthedocs.io/en/latest/pipeline/options.html).
  A run-specific tensor-fit command log was not available, so this is an
  evidence-based inference, not direct command-level proof.

The requested `../motor/results/` placeholder scaffold has been removed.
Completed Motor results remain under `../motor/v1/results/`.

## Reproduction code

All four files in `code/` participate in the current reproduction workflow:

| Code | Role and dependency |
| --- | --- |
| `run_motor_paired_auc_v2.m` | Independent MATLAB entry point for the paired table. Contains `paired_stratified_auc` and its fitting/bootstrap helpers. Requires clinical data, Standard/Classical metric tables, Motor extraction tables, and the three reference result tables supplied as arguments. |
| `build_tables.py` | Python stage 1: produces the final classical-exclusion table and a 133-row, 250-column base private predictor CSV. |
| `augment_private_predictors_v2.py` | Python stage 2: consumes stage 1's private CSV and adds 32 Figure 5 plotting-workbook values, giving 282 columns. |
| `build_mounted_figure_followup.py` | Python stage 3: consumes stage 2's private CSV plus mounted plotting workbooks; produces the final matched Figure 4 summaries, 270 post-hoc comparisons, and the final 284-column private CSV. |

Stage 1 uses the frozen **90-row legacy tract lesion-load workbook**, not
the mounted 89-row tract workbook used by stage 3. Its SHA-256 is
`f2f48e6338b0c6767e49e462a4248dbbb69d7d963f60c594212fa06e81ad8190`.
Its additional `lesion_load_by_region_subtype.csv`, `posthoc_group_tests.csv`,
and report are intermediate historical outputs, not the final Figure 4/5
deliverables listed above. Use the stage 3 files ending in
`caption_matched_89.csv` for this handoff. Python requires pandas, NumPy,
SciPy, statsmodels, and an Excel-reading engine; MATLAB requires Statistics
and Machine Learning Toolbox. Each Python command documents its explicit
path arguments with `--help`.

The complete three-stage Python chain was rerun from source data on
15 October 2026. The four delivered aggregate tables it produces and the
shared private CSV were byte-identical to the current files. The MATLAB
paired table had already passed a clean rerun and canonical reference checks.
These scripts take explicit source/output paths and refuse to overwrite
existing outputs. Intermediate results from earlier stages are not part of
the handoff list above.
