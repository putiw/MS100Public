# Rebuild the plain Figure 5

From the repository root in MATLAB, run:

```matlab
addpath(fullfile(pwd,'revision','code'));
figure5_publication;
```

The renderer reads the tracked
[`figure5_plotted_values.csv`](../results/model_results/figure5_plotted_values.csv)
through `load_revision_figure_data.m` and writes
`revision/results/figures/Figure5.pdf`. The CSV contains the means and SEMs
plotted for each framework, tissue, metric, region, and subgroup.

To export a separate copy, pass an output directory as the second argument.
An empty first argument selects the tracked aggregate CSV.

```matlab
figure5_publication([],fullfile(tempdir,'figure5-review'));
```
