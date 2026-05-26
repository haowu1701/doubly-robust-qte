# Inference Simulation

This folder contains the simulation code used to produce Table 2 in the main manuscript.

## Folder Structure

- `Empirical variance/`: generates the simulated datasets and calculates empirical variances.
- `sandwich_expectation/`: variance estimator based on taking expectation first.
- `sandwich_derivation/`: variance estimator based on taking derivation first.
- `IF_based var/`: influence-function-based variance estimator.
- `BootstrapnBoot500/`: nonparametric bootstrap with 500 bootstrap replicates.
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

Generated data and results are not included in the repository. This includes `.RData`, `.rds`, `.csv`, and cluster `.out` files.

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
- `nsim1000_data_n1000_mc.RData`
- `nsim1000_results_n1000_mc.csv`
- `simulation_summary_all_scenarios_long.csv`
- `simulation_summary_all_scenarios_wide.csv`

Before running the other variance estimators, copy or link the generated dataset and true-value file into the corresponding scenario folder, using the filenames expected by those scripts:

- `nsim1000_data_n1000.RData`
- `true_values.RData`

For example, for scenario `mc`, place those two files in `sandwich_expectation/mc/`, `sandwich_derivation/mc/`, `IF_based var/mc/`, and `BootstrapnBoot500/mc/`, then run:

``` r
source("sand_expectation_RUN.R")
source("sand_derivation_RUN.R")
source("IF_RUN.R")
source("Bootstrap_RUN.R")
```

Run the analogous scripts inside `cc/`, `cm/`, `mc/`, and `mm/` to reproduce all Table 2 entries.
