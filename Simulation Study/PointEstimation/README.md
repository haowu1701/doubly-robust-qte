# Point Estimation Simulation

This folder contains the simulation code used to produce Table 1 in the main manuscript and Tables S3.1 to S3.9 in the supplementary material.

## Folder Structure

- `nsample500/`: simulations with sample size 500.
- `nsample1000/`: simulations with sample size 1000.
- `QTE0.1/`, `QTE0.25/`, `QTE0.5/`, `QTE0.75/`, `QTE0.9/`: simulations targeting the corresponding quantile treatment effect.

Each `QTE*` folder contains:

- `true_value_generate.R`: generates the true quantiles and true QTE values.
- `Helper_functions.R`: data generation and estimator helper functions.
- `single_sim.R`: runs one simulation replicate.
- `RUN.R`: runs all simulation replicates and summarizes results.
- `haojob.sh`: example SLURM script for running the simulation on a cluster.

## Output Tables

The scripts in this folder generate the point-estimation results reported in:

- Main manuscript: Table 1.
- Supplementary material: Tables S3.1-S3.9.

The exact table layout is generated from the summary files written by `RUN.R`.

## Running the Code

Run each simulation from inside its corresponding folder. Please generate the true-value data first. `RUN.R` expects `true_values.RData` to be available in the same folder.

For example:

```r
setwd("nsample1000/QTE0.5")
source("true_value_generate.R")
source("RUN.R")
```

The required order is:

1. Run `true_value_generate.R` to create `true_values.RData`.
2. Run `RUN.R` to perform the simulation and generate summary files.

