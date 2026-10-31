# Rebuild the plain Figure 5

Figure 5 can be rebuilt from this repository with MATLAB R2024b. No original
workbooks, participant-level tables, mounted drive, `tmp/` input or fitted models
are needed. The tracked renderer uses a frozen aggregate table and MATLAB default colours.
Optional manuscript styling lives separately in the ignored `revision/local_style/`
folder; see `README_figures.md` for both workflows.

From the repository root (`MS100Public`):

```matlab
addpath(fullfile(pwd, 'revision', 'code'));
figure5_publication;
```

This writes `revision/results/figures/Figure5.pdf`. To keep the canonical figure
untouched and write a separate review export:

```matlab
figure5_publication([], fullfile(tempdir, 'figure5-review'));
```

The output directory is created if necessary. The default input is resolved from
the script location, so changing the current directory does not change the input.
An explicit alternate aggregate CSV can still be passed as the first argument.
Use the default input to reproduce the committed figure.

## Tracked dependencies

- `revision/code/figure5_publication.m`: twelve-panel compositor.
- `revision/results/model_results/figure5_plotted_values.csv`: 160 frozen
  region/subtype/metric aggregate rows, byte-identical to the audited table used
  for the 31 October 2026 final figure.
- `revision/code/load_revision_figure_data.m`: shared frozen-input loader.
  The older `fig5_T1_radar.m`, `fig5_MTR_radar.m`, `fig5_FA_radar.m` and
  `fig5_MD_radar.m` scripts still require their original workbooks; the
  self-contained aggregate rebuild uses `figure5_publication` instead.

Aggregate CSV SHA-256:

```text
328352ef2c88dccea253090b73755d124a3d5a3951fbc53629af9b824076e218
```

## Data and provenance

Each row describes a group summary, not a participant. It includes the framework,
tissue, metric, region, MS subtype (or Control), sample count, mean, SD, SEM,
median, quartiles, zero count/policy and a source-workbook SHA-256. There are no
subject identifiers, participant-level measurements, machine-specific paths or
contact details. The renderer reads `Framework`, `Tissue`, `Metric`, `Region`,
`Subtype`, `Mean` and `SEM`; the remaining fields preserve the aggregate audit.

The panels are classical lesion, tract lesion and tract NAWM for T1, MTR, FA and
MD. Controls occur only in tract NAWM. Exact-zero exclusions follow the original
Figure 5 convention. T1 means are in seconds, MTR in percent, FA is unitless, and
MD is stored in mm²/s and multiplied by 1,000 for display in 10^-3 mm²/s. Sample
counts and SEMs are frozen per cell; no new cohort selection or model fitting
occurs during rendering.

Before publication, all 160 means/SEMs/counts were checked against recovered
original `GroupTractT1.xlsx`, `GroupTractMTR.xlsx`, `GroupTractFA.xlsx`,
`GroupTractMD.xlsx` and `icometrixLesionStats.xlsx` workbooks. All counts agreed;
maximum mean/SEM differences were floating-point rounding (7.11e-15 / 2.22e-16).
The tract source audit additionally matched all 133 IDs, validity masks and 3,432
valid subject-feature values. The private inputs used for those checks remain
outside this repository and are not required to render the frozen figure. The
classical workbook matched its recorded byte hash; recovered tract workbooks
matched the audited values and membership, rather than the original file bytes.
The suggested later staging copies were not substituted because some T1/FA/MD
summaries and counts differ from the original figure.

## Appearance and validation

The tracked entry point now produces a plain MATLAB version, as requested when
analysis plotting and local manuscript styling were separated. Scientific values,
radial limits, mean/SEM interpretation and units are unchanged. The canonical
tracked PDF is for analysis review; use the ignored local styling entry point for
submission candidates. It uses the original manuscript MS-group hues with grey
controls and is documented in `README_figures.md`.

The earlier reproducibility check confirmed a byte-identical 160-row aggregate
input and a pixel-identical rebuild before the intentional plain/styled split.
That pixel-identity claim describes the earlier version, not the deliberately
restyled outputs. The current split shares one frozen-input loader and is checked
against every source cell without running clinical models.
