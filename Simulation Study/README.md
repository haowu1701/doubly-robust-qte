# Simulation Studies

This folder contains the simulation code for the doubly robust QTE and PTE analyses, including point estimation, inference, figure generation, and computational time comparisons.

Each main subfolder contains a more detailed `README.md` with script-level instructions.

## Folder Structure

``` text
Simulation/
├── QTE/
│   ├── PointEstimation/
│   └── Inference/
├── PTE/
│   ├── PointEstimation/
│   └── PTE inference/
├── trueDGP/
├── Figure S2.1/
├── Figure S3.1/
└── Figure S3.2/
```

## QTE Simulations

`QTE/` contains simulation studies for quantile treatment effects.

- `QTE/PointEstimation/`: point-estimation simulations for QTE. These scripts generate the results reported in the main manuscript Table 1 and Supplementary Tables S3.1-S3.9.
- `QTE/Inference/`: inference simulations for QTE, including empirical variance, sandwich variance, influence-function-based variance, and bootstrap variance. These scripts generate the results reported in the main manuscript Table 2 and Supplementary Table S3.10.

See `QTE/README.md` and the README files inside each QTE subfolder before running the code.

## PTE Simulations

`PTE/` contains simulation studies for probability treatment effects.

- `PTE/PointEstimation/`: point-estimation simulations for PTE. This folder generates simulated datasets and point-estimation summaries.
- `PTE/PTE inference/`: inference simulations for PTE using empirical, bootstrap, sandwich, and influence-function-based methods.

These scripts generate the results reported in Supplementary Table S3.11.

See `PTE/README.md` and the README files inside each PTE subfolder before running the code.

## Figure Scripts

The figure folders contain code for supplementary figures.

- `Figure S2.1/`: code for the quantile interpolation illustration.
- `Figure S3.1/`: code for simulated data visualizations.
- `Figure S3.2/`: code for computational time comparisons.

Each figure folder contains a README describing its scripts and output files.

## True Data-Generating Process

`trueDGP/` contains helper code for generating true-value objects used by the simulation studies.
