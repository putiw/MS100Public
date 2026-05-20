# MS100Public

Analysis code for:

**Tract-based Quantitative MRI for Resolving the Clinico-Radiological Paradox in Multiple Sclerosis**

Abdullah O\*, Wen P\*, Abdelrazeq AW\*, Matrosova M\*, Brylev L, Ahmad A, Bryukhov V, Melcher D, Rokers B. *Brain Communications*, 2026.

## Setup

1. Clone this repository.
2. Edit `configs/config.json` and set `bidsDir` to the root of your BIDS-formatted dataset.
3. Add the repository to the MATLAB path:
   ```matlab
   addpath(genpath('/path/to/MS100Public'));
   ```

## Repository Structure

```
MS100Public/
├── configs/                              % Configuration files
│   ├── config.json                       % Data paths (edit bidsDir)
│   ├── tract_names.txt                   % 96 individual white matter tract names
│   └── tract_groups/                     % Tract-to-group assignments (41 tracts → 4 groups)
├── helpers/                              % Shared utility functions
├── outputs/                              % (gitignored) Generated outputs
└── *.m                                   % Analysis and figure scripts
```

## Script-to-Figure/Table Mapping

### Main Figures

| Figure | Script | Description |
|--------|--------|-------------|
| 1 | `fig1_clinical_distributions.m` | Clinical score distributions and radar plots |
| 2 | *(schematic, not code-generated)* | MRI acquisition and preprocessing pipeline |
| 3 | *(rendered in 3D viewer)* | Classical vs tract-based lesion topography |
| 4 | `fig4_lesionload_tract_radar.m`, `fig4_lesionload_classical_radar.m` | Lesion load radar plots |
| 5 | `fig5_T1_radar.m`, `fig5_FA_radar.m`, `fig5_MD_radar.m`, `fig5_MTR_radar.m` | Lesion composition radar plots by qMRI metric |
| 6 | `rev_fig6_r2_comparison.m` | R-squared comparison: tract-based vs classical |

### Tables

| Table | Script | Description |
|-------|--------|-------------|
| 1 | `table1_verify.m` | Demographic and clinical summary statistics |
| 2 | `paper_multi_bi_ridge_all.m`, `rev_qMRI.m` | Binary ridge regression (EDSS, MSPro) |

### Supplementary Materials

| Item | Script | Description |
|------|--------|-------------|
| Figure S1 | `figS1_correlation_heatmaps.m` (calls `rev_corr.m`) | Correlation heatmaps with Bonferroni correction |
| Table S3 | `rev_qMRI.m` | Ridge regression with demographic covariates |
| Table S4 | `paper_bh_p_correction.m`, `rev_p_correct.m` | Bonferroni-Holm corrected p-values |
| Table S5 | `rev_lesionload.m` | Lesion load ridge regression |

### Data Generation (run first)

| Script | Description |
|--------|-------------|
| `paper_gen_excel_stats.m` | Extracts qMRI statistics per tract group |
| `paper_multi_bi_ridge_all.m` | Runs all binary ridge regressions |
| `paper_multi_continuous_ridge_all.m` | Runs all continuous ridge regressions |
| `paper_bh_p_correction.m` | Applies Bonferroni-Holm correction to regression p-values |

## Dependencies

- MATLAB R2024b or later
- Statistics and Machine Learning Toolbox

## Data Availability

The clinical and imaging data cannot be shared publicly due to patient privacy. De-identified data are available from the corresponding author on reasonable request, subject to ethics committee approval.

## License

This code is provided for academic research purposes. Please cite the paper if you use this code.
