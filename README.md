# MS100Public

MATLAB analysis code for:

**Tract-based Quantitative MRI for Resolving the Clinico-Radiological Paradox in Multiple Sclerosis**

Abdullah O\*, Wen P\*, Abdelrazeq AW\*, Matrosova M\*, Brylev L, Ahmad A, Bryukhov V, Melcher D, Rokers B.

## Overview

This repository contains the public analysis scripts used to generate the tract-based quantitative MRI figures and tables for the manuscript. The code operates on a local BIDS-formatted multiple sclerosis dataset plus derivative spreadsheets and MRI maps. Patient-level clinical and imaging data are not included in this repository.

## Requirements

- MATLAB R2024b or later
- Statistics and Machine Learning Toolbox
- Image Processing Toolbox, for NIfTI I/O
- Parallel Computing Toolbox, for scripts that use `parpool` and `parfor`

## Setup

1. Clone this repository.
2. Edit [configs/config.json](configs/config.json) and set `bidsDir` to the root of your local BIDS dataset.
3. From MATLAB, add the repository to the path:

   ```matlab
   addpath(genpath('/path/to/MS100Public'));
   ```

4. Create the local output directory used by the regression scripts:

   ```matlab
   mkdir(fullfile('/path/to/MS100Public', 'outputs'));
   ```
## Repository Structure

```text
MS100Public/
  configs/
    config.json              Data paths and output filenames
    tract_names.txt          White matter tract names
    tract_groups/            Tract-to-group assignments
  helpers/                   Shared loading, plotting, and statistics helpers
  outputs/                   Local generated outputs; not tracked
  *.m                        Analysis and figure scripts
```


### Main Figures

| Figure | Script | Description |
| --- | --- | --- |
| 1 | `fig1_clinical_distributions.m` | Clinical score distributions and radar plots |
| 2 | Not code-generated | MRI acquisition and preprocessing schematic |
| 3 | Rendered externally | Classical versus tract-based lesion topography |
| 4 | `revisionExtras/code/fig4_plain_matched89.m`; original radar scripts: `fig4_lesionload_tract_radar.m`, `fig4_lesionload_classical_radar.m` | Reproducible 89-person bar and six-radar MATLAB figures |
| 5 | `fig5_T1_radar.m`, `fig5_FA_radar.m`, `fig5_MD_radar.m`, `fig5_MTR_radar.m` | Lesion and normal-appearing white matter qMRI radar plots |
| 6 | `rev_fig6_r2_comparison.m` | R-squared comparison for tract-based versus classical models |

### Tables

| Table | Script | Description |
| --- | --- | --- |
| 1 | `table1_verify.m` | Demographic and clinical summary statistics |
| 2 | `paper_multi_bi_ridge_all.m`, `rev_qMRI.m` | Binary ridge regression summaries for EDSS and MSPro |

### Supplementary Materials

| Item | Script | Description |
| --- | --- | --- |
| Figure S1 | `figS1_correlation_heatmaps.m` | Correlation heatmaps with Bonferroni correction |
| Table S3 | `rev_qMRI.m` | Ridge regression with demographic covariates |
| Table S4 | `paper_bh_p_correction.m`, `rev_p_correct.m` | Bonferroni-Holm corrected p-values |
| Table S5 | `rev_lesionload.m` | Lesion-load ridge regression |

## Configuration Notes

- `configs/tract_names.txt` lists the tract names used by the analysis.
- `configs/tract_groups/*.txt` defines the tract groupings used in group-level analyses.
- `configs/config.json` controls the BIDS root, derivatives folders, output spreadsheet names, and plotting colors.
- Generated regression outputs are written to `outputs/` inside this repository.
- Generated qMRI spreadsheet outputs are written to `<bidsDir>/<statsDir>/` unless a script is run in test mode.

### Figure 4 plain plots on the matched imaging cohort

`revisionExtras/code/fig4_plain_matched89.m` reads `GroupTractLesionLoad.xlsx` and
`icometrixLesionLoad.xlsx` from the configured stats directory. It restricts
the classical file to the 89 participants in the tract file (RRMS 49, PPMS 17,
SPMS 23; `sub-RR047` is the sole classical-only row). It produces two separate
MATLAB figures: a Figure 4A bar chart and a Figure 4B six-polar-panel plot.
The manuscript's Illustrator layout is a separate, local styling step.

From MATLAB, add `revisionExtras/code` to the path and choose a new output
directory for the two vector PDFs:

```matlab
addpath('/path/to/MS100Public/revisionExtras/code');
fig4_plain_matched89('/path/to/new/figure4-plain-output');
```

The function refuses an existing output directory. If the configured volume
is not mounted, pass the tract and classical workbook paths as its second and
third arguments. Figure 4A shows mean ± sample SD in mL, with a 0–50 mL
axis; its tract bars sum the four tract-region lesion volumes and are not a
unique-voxel union. Figure 4B retains real zeros and shows regional mean ±
SEM. Lesion-volume fields are converted from mm³ to mL, and each radial
maximum is 1.1 times its largest plotted mean + SEM, matching the original
MATLAB radar calculation.

## Data Availability

The clinical and imaging data cannot be shared publicly due to patient privacy. De-identified data may be available from the corresponding author on reasonable request, subject to ethics committee approval.

## License

This code is provided for academic research purposes. Please cite the paper if you use this code.
