# Doubly Robust QTE

This repository contains simulation studies and real-data analysis code for doubly robust estimation of "Doubly Robust Estimators of Quantile Treatment Effects With Semiparametric Cumulative Probability Models".

The code is organized into two main parts:

- `Simulation/`: simulation studies for evaluating point estimation, inference, and computational performance.
- `RealData/`: real-data analysis workflows using synthetic example datasets that mirror the analysis structure.

## Repository Structure

``` text
.
├── RealData/
│   ├── CD4outcome/
│   └── VLoutcome/
└── Simulation/
    ├── QTE/
    ├── PTE/
    ├── trueDGP/
    ├── Figure S2.1/
    ├── Figure S3.1/
    └── Figure S3.2/
```

## Simulation Studies

The `Simulation/` directory contains code for the simulation results reported in the manuscript and supplementary materials.

- `QTE/`: simulation studies for quantile treatment effect estimation and inference.
- `PTE/`: simulation studies for probability treatment effect estimation and inference.
- `trueDGP/`: helper code for generating true values under the simulation data-generating processes.
- `Figure S2.1/`: code for the quantile interpolation illustration.
- `Figure S3.1/`: code for simulated data visualization.
- `Figure S3.2/`: code for computational time comparisons.

Each simulation subfolder includes its own `README.md` with more detailed instructions and file descriptions.

## Real-Data Analyses

The `RealData/` directory contains two analysis examples:

- `CD4outcome/`: analysis using CD4 count as the outcome.
- `VLoutcome/`: analysis using viral load as the outcome.

The original real datasets are not included in this repository because they contain confidential patient-level information. Instead, each real-data analysis folder includes a synthetic dataset and a script that documents how the synthetic data were generated.

Synthetic datasets are intended for reproducing the code workflow and checking that the analysis scripts run correctly. They should not be interpreted as real clinical data.
