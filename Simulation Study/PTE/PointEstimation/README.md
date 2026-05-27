# PTE Point Estimation Simulation

This folder contains the point-estimation simulation code for the probability treatment effect (PTE).

## Files

- `true_value_generate.R`: generates the large-sample truth object `true_values.RData`.
- `RUN.R`: main simulation script. It runs all scenarios and writes simulation results and summary tables.
- `single_sim.R`: runs one simulation replicate for a selected scenario.
- `Helper_functions.R`: data-generation, conditional CDF, quantile, and summary helper functions.
- `PTE_helpers.R`: PTE point-estimation helper functions for OR, IPW, and AIPW estimators.
- `haojob.sh`: example SLURM batch script for running `RUN.R` on a cluster.

## Scenarios

The main script runs four scenarios:

- `cc`: correct propensity score model and correct outcome model.
- `cm`: correct outcome model, misspecified propensity score model.
- `mc`: correct propensity score model, misspecified outcome model.
- `mm`: misspecified propensity score model and misspecified outcome model.

## How to Run

Run the scripts from this folder.

First generate the large-sample truth file:

``` bash
Rscript true_value_generate.R
```

This creates:

``` text
true_values.RData
```

Then run the simulation:

``` bash
Rscript RUN.R
```

On a SLURM cluster, submit:

``` bash
sbatch haojob.sh
```

Note: `haojob.sh` assumes `true_values.RData` already exists. Run `true_value_generate.R` first if the truth file has not been generated.
