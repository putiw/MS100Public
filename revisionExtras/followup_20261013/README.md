# Follow-up analyses

The `results/` directory contains aggregate tables and a Figure 4 PDF.

| File | Contents |
| --- | --- |
| `paired_motor_delta_auc.csv` | Paired apparent AUC differences for seven metrics, two binary outcomes, adjusted and unadjusted fits, and three model pairs. |
| `lesion_load_by_region_subtype_caption_matched_89.csv` | Figure 4B summaries by region and subtype for lesion number, volume, and normalized volume. |
| `figure4a_bar_by_subtype_caption_matched_89.csv` | Figure 4A whole-brain and tract lesion volume means and SDs in mL. |
| `classical_exclusions_by_metric.csv` | Counts of unavailable classical-region predictors by metric. |
| `posthoc_fig4_fig5_caption_matched_89.csv` | Pairwise group tests for Figures 4A, 4B, and 5. |
| `Figure4_matched89_original_style.pdf` | Figure 4 for the 89 imaging participants used in both frameworks (49 RRMS, 17 PPMS, 23 SPMS). |

## Rebuild the tables

The Python scripts require pandas, NumPy, SciPy, statsmodels, and an Excel
reader such as openpyxl. From this directory, run each script with `--help`
for its input and output arguments. Supply the source workbooks and predictor
tables separately, write aggregate outputs to a fresh directory, and keep
participant-level CSVs outside the repository.

1. Run `python code/build_tables.py --help`. Use the 90-row legacy tract
   lesion-load workbook, the classical lesion-load workbook, the clinical
   workbook, and the metric and Motor tables. This stage produces
   `classical_exclusions_by_metric.csv` and a base predictor CSV.
2. Run `python code/augment_private_predictors_v2.py --help`. Pass the base
   predictor CSV from step 1 and the Figure 5 plotting workbooks to add
   participant-level plotting values for the next stage.
3. Run `python code/build_mounted_figure_followup.py --help`. Pass the predictor
   CSV from step 2 and the plotting workbooks, including the 89-row tract
   lesion-load workbook. The files ending in `caption_matched_89.csv` are
   the matched Figure 4 summaries and the corresponding Figure 4/5 tests.

To rebuild `paired_motor_delta_auc.csv`, use MATLAB with Statistics and
Machine Learning Toolbox. From this directory, add `code/` to the MATLAB path
and call:

```matlab
addpath(fullfile(pwd,'code'));
run_motor_paired_auc_v2(statsDir,metricsDir,motorDir,motorReference, ...
    edssReference,msproReference,outputFile);
```

Supply the clinical workbook in `statsDir`, the Standard and Classical
predictor tables in `metricsDir`, and the Motor tables in `motorDir`.
The three reference arguments are the Motor, EDSS, and MSPro result CSVs.
Set `outputFile` to a new CSV path.

## Open the plain Figure 4 plots

From the repository root in MATLAB, run:

```matlab
addpath(fullfile(pwd,'revision','manuscript_figures','matlab'));
open_figure4();
```

This opens the Figure 4A bars and Figure 4B radar plots from the two matched
aggregate CSVs in `results/`.
