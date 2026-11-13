# NODDI group summaries and Figure S3

Figure S3 follows the three rows of Figure 5 for NDI, ODI and ISOVF
(the FWF field). The two lesion rows summarize means within lesion voxels
for RRMS, SPMS and PPMS. Each region uses its available observations, with
lesion-free regions recorded as missing. The tract NAWM row also includes
controls. Fitted zeros are retained and counted.

- `results/Supplementary_Table_NODDI_group_summary.csv` gives N, mean, SD,
  zero count and Kruskal-Wallis results by framework, tissue, metric and region.
- `results/Supplementary_Table_NODDI_pairwise.csv` gives two-sided asymptotic
  Mann-Whitney tests for every group pair.
- `results/figureS3_plotted_values.csv` gives the plotted means and SEMs
  in Figure 5's aggregate schema.
- `results/Supplementary_FigureS3.pdf` is the plain MATLAB export.

The figure has two Benjamini-Hochberg families, with 36 Kruskal-Wallis
and 144 pairwise tests. Each panel has its own radial limits, rounded
outward to 0.01 from the minimum and maximum mean ± SEM.

## Rebuild

Run `code/build_noddi_group_summary.py` with `--metrics-dir` pointing to
the `GroupTract{FA,NDI,ODI,FWF}_All.csv`, `ClassicalRegionAllMetrics.csv`
and `ClassicalRegionNODDI_All.csv` inputs. Pass the frozen Figure 5
plotted-values CSV with `--figure5`, a fresh aggregate directory with
`--output-dir` and a private directory with `--internal-dir`.

The script checks N, mean, SD and SEM against all 28 tract and 12 classical
Figure 5 FA cells at a tolerance of 1e-6 before opening the NODDI inputs.
It records the zero-excluded sensitivity results in the private directory.

Add `code/` to the MATLAB path and call
`render_supplementary_figure_s3(plottedCsv, outputPdf)` with the aggregate
plotted-values CSV and a new PDF path.
