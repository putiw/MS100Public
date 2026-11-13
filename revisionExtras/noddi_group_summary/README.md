# NODDI group summaries

Aggregate group summaries for tract NDI, ODI and ISOVF (the GroupTract FWF
field) and Supplementary Figure S3. No participant-level data are written
here.

## Contents

- `results/Supplementary_Table_NODDI_group_summary.csv`: N, mean and SD per
  group, tract cluster and tissue (Tail and NAWM), with the Kruskal-Wallis
  test for each feature.
- `results/Supplementary_Table_NODDI_pairwise.csv`: two-sided Mann-Whitney
  tests for every pair of groups.
- `results/figureS3_plotted_values.csv`: means and SEMs drawn in Figure S3
  (lesion mean and NAWM), in the Figure 5 plotted-values format.
- `results/Supplementary_FigureS3.pdf`: plain MATLAB export of Figure S3.

Tail is the percentile summary used as the tract predictor in the regression
models: P10 for NDI and P90 for ODI and ISOVF, taken over all voxels of the
tract cluster. It is not restricted to lesion voxels and is defined for
controls as well as MS groups. The figure follows Figure 5 and shows the mean
within lesion voxels (RRMS, SPMS, PPMS) and the mean within NAWM (those
groups plus controls).

Benjamini-Hochberg correction is applied within this analysis only: one
family of 24 Kruskal-Wallis tests and one family of 144 pairwise tests.
Fitted zeros are retained; no cell contains a zero.

## Rebuild

```
python code/build_noddi_group_summary.py \
  --metrics-dir <directory with GroupTractFA_All.csv, GroupTractNDI_All.csv,
                 GroupTractODI_All.csv, GroupTractFWF_All.csv> \
  --figure5 revision/results/model_results/figure5_plotted_values.csv \
  --output-dir <empty output directory> \
  --internal-dir <private directory outside the repository>
```

The script stops before reading NODDI values unless all 28 tract FA cells of
Figure 5 (N, mean, SD, SEM) are reproduced within 1e-6. NODDI values are
available for 130 participants, all within the Figure 5 cohort.

Render the figure in MATLAB:

```matlab
addpath('revisionExtras/noddi_group_summary/code');
render_supplementary_figure_s3( ...
  'revisionExtras/noddi_group_summary/results/figureS3_plotted_values.csv', ...
  fullfile(tempdir, 'Supplementary_FigureS3.pdf'));
```
