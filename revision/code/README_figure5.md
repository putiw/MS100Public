# Rebuild Figure 5

Figure 5 can be rebuilt from this repository with MATLAB R2024b. No original
workbooks, participant-level tables, mounted drive, `tmp/` input or fitted models
are needed. The renderer uses a frozen aggregate table and a Figure-5-only palette.

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
- `helpers/fig5_palette.m`: shared palette used by this compositor and the original
  `fig5_T1_radar.m`, `fig5_MTR_radar.m`, `fig5_FA_radar.m`, `fig5_MD_radar.m` scripts.
  Those older scripts still require their original workbooks; the self-contained
  publication rebuild uses `figure5_publication` instead.

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

The existing twelve-panel layout, radial limits, means, SEM bands and colours are
unchanged by making the input repository-local. The palette is RRMS blue
`#0072B2`, SPMS vermilion `#D55E00`, PPMS orange `#E69F00`, controls grey `#444444`.
It is confined to Figure 5; the global configuration and Figure 1 are unchanged.
The previous final-figure audit checked deuteranopia and protanopia simulations.

The final PDF page is 481 × 595 pt (approximately 169.7 × 209.9 mm), with vector
content and minimum effective text size 7.2 pt, including the MD superscripts.
The repository-local-input rebuild was rendered at 300 dpi and compared with the
previous committed figure; validation details are recorded in the external agent
handoff. The canonical PDF is retained when the rendered output is unchanged.
