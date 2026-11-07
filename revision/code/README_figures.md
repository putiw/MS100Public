# Rebuild the plain MATLAB figures

From the repository root in MATLAB, run:

```matlab
addpath(fullfile(pwd,'revision','code'));
figure5_publication;
figure6_r2_comparison_sdmt_msfc_noddi;
supplementary_figure_s1_correlations_sdmt_msfc;
```

The commands write the following PDFs to `revision/results/figures/`.
The input files are tracked under `revision/results/model_results/`.

| Output | Aggregate input |
| --- | --- |
| `Figure5.pdf` | `figure5_plotted_values.csv` |
| `Figure6.pdf` and `Supplementary_FigureS2.pdf` | `figure6_full_models.csv` and `figure6_matched_models.csv` |
| `Supplementary_FigureS1.pdf` | `supplementary_figure_s1_correlations.json` |

To open the figures in MATLAB without exporting, use the
[`open_figure*` entry points](../manuscript_figures/README.md).
