# Plain figures and local manuscript styling

The scientific analyses and presentation are separate:

1. Analysis scripts such as `run_all_qmri_manuscript_models.m` produce the frozen
   aggregate results in `revision/results/model_results/`. Those analyses were
   not changed or rerun during the styling separation.
2. `load_revision_figure_data.m` reads those tracked aggregates and preserves their
   model/region/outcome mapping. It does no fitting or styling. Both renderers use
   this same data loader.
3. The tracked plotting entry points below produce plain MATLAB figures. They
   work on a fresh checkout without the local styling folder.
4. Optional manuscript styling is kept in the **ignored**, local-only
   `revision/local_style/` folder. It reads the same frozen values and redraws
   presentation geometry, colours, typography and spacing. It does not calculate
   new statistics or select cohorts.

## Tracked plain outputs

From the repository root in MATLAB R2024b:

```matlab
addpath(fullfile(pwd,'revision','code'));
figure5_publication;
figure6_r2_comparison_sdmt_msfc_noddi;
supplementary_figure_s1_correlations_sdmt_msfc;
```

The outputs in `revision/results/figures/` are now **plain analysis-review
figures, not submission artwork**:

- `Figure5.pdf`
- `Figure6.pdf`
- `Supplementary_FigureS1.pdf`
- `Supplementary_FigureS2.pdf`

The historical function name `figure5_publication` is retained for compatibility;
its tracked default now intentionally uses MATLAB colours and ordinary axes.
The figure functions accept alternate aggregate inputs/output directories as
before. They do not depend on the optional local styling files. The aggregate
inputs are unchanged; see `README_figure5.md` for Figure 5 provenance.

For the plain S1 row labels: A = Association, C = Cerebellar, OP =
Occipitoparietal, BP = Brainstem/Projection, PV = Periventricular, JC =
Juxtacortical, IT = Infratentorial, DW = Deep white matter. Complete region labels
remain available in the frozen data and shared loader.

## Optional local manuscript outputs

On the working laptop, if the ignored styling folder is present:

```matlab
addpath(fullfile(pwd,'revision','local_style'));
make_manuscript_figures;
```

This writes the **styled submission candidates** to
`revision/local_style/output/`, leaving the tracked plain PDFs untouched.
Individual renderers are `style_figure5`, `style_figure6` (Figure 6 and S2) and
`style_s1`; `manuscript_palette` defines the local Figure 5 colours.

The styling folder and its generated outputs are deliberately not tracked by
Git. A new clone will reproduce the plain figures but will not contain the local
styles. Copy that folder separately if another agent or machine needs to edit or
recreate the manuscript artwork. The local handoff records its full paths and a
separate backup archive.

The manuscript styling follows the original figure references: Arial/Helvetica,
bold italic labels where used originally, thin radar lines and translucent SEM
bands, the original side key, blue/red/yellow MS groups and neutral-grey controls,
and the original blue/peach bar comparison. Added outcomes/metrics and the
minimum 7 pt rule prevent an exact copy of the old crowded Figure 6 layout; five
stacked panels are retained. Figure 4 is unchanged under its earlier explicit
instruction not to regenerate it.
