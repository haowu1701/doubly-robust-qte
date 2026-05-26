# Inference Simulation

This folder contains the simulation code used to produce Table 2 in the main manuscript.

## Folder Structure

- `Empirical variance/`: generates the simulated datasets and calculates empirical variances.
- `Sandwich_expectation variance/`: sandwich variance estimator based on taking expectation first.
- `Sandwich_derivation variance/`: sandwich variance estimator based on taking derivation first.
- `IF_based variance/`: influence-function-based variance estimator.
- `Bootstrap variance/`: nonparametric bootstrap with 500 bootstrap replicates.
- `Coverage probability.R`: helper script for summarizing coverage probabilities and variance estimates from generated results.

Each variance-estimator folder contains four simulation scenarios:

- `cc/`: correctly specified CPM and correctly specified propensity score model.
- `cm/`: correctly specified CPM and misspecified propensity score model.
- `mc/`: misspecified CPM and correctly specified propensity score model.
- `mm/`: misspecified CPM and misspecified propensity score model.

Each scenario folder contains some or all of the following scripts:

- `true_value_generate.R`: generates the true quantiles and true QTE values.
- `Helper_functions.R`: data generation and estimator helper functions.
- `variance_helpers.R`: helper functions for variance estimation.
- `single_sim.R`: runs one empirical-variance simulation replicate.
- `RUN.R`: runs empirical-variance simulations and saves simulated datasets.
- `sand_expectation_RUN.R`: runs the expectation-first sandwich variance estimator.
- `sand_derivation_RUN.R`: runs the derivation-first sandwich variance estimator.
- `IF_RUN.R`: runs the influence-function-based variance estimator.
- `Bootstrap_RUN.R`: runs the bootstrap variance estimator.
- `haojob.sh`: example SLURM script for running the corresponding simulation on a cluster.

## Output Tables

The scripts in this folder generate the inference results reported in:

- Main manuscript: Table 2.
- Supplementary material: Table S3.10.

Generated data and results are not included in the repository.

## Running the Code

Run each simulation from inside its corresponding scenario folder.

First, run the empirical-variance simulations. These generate both the empirical-variance summaries and the simulated datasets used by the other variance estimators.

For example:

``` r
setwd("Empirical variance/mc")
source("true_value_generate.R")
source("RUN.R")
```

`RUN.R` creates files such as:

- `true_values.RData`
- `nsim1000_data_n1000_<scenario>.RData`
- `nsim1000_results_n1000_mc.csv`
- `simulation_summary_all_scenarios_long.csv`
- `simulation_summary_all_scenarios_wide.csv`

The empirical-variance simulations use the same generated datasets across `cc`, `cm`, `mc`, and `mm`; the scenario only changes the working model specification. Therefore, before running the other variance estimators, copy or link the generated dataset into the corresponding scenario folder and rename it to the common filename expected by those scripts:

- `nsim1000_data_n1000.RData`
- `true_values.RData`

For example, for scenario `mc`, place those two files in `Sandwich_expectation variance/mc/`, `Sandwich_derivation variance/mc/`, `IF_based variance/mc/`, and `Bootstrap variance/mc/`, then run:

``` r
source("sand_expectation_RUN.R")
source("sand_derivation_RUN.R")
source("IF_RUN.R")
source("Bootstrap_RUN.R")
```

Run the analogous scripts inside `cc/`, `cm/`, `mc/`, and `mm/` to reproduce all Table 2 entries.

## Supplementary Table S3.10

Table S3.10 in the supplementary material uses only the sandwich variance estimators. To reproduce it, run the scripts in:

- `Sandwich_derivation variance/`
- `Sandwich_expectation variance/`

for the additional quantile levels `0.1`, `0.25`, `0.75`, and `0.9`.
