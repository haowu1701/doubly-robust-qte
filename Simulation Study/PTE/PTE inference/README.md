# PTE Inference

This folder contains simulation code for inference of the probability treatment effect (PTE).

## Inference Methods

`Empirical/` computes empirical simulation summaries across repeated simulated datasets.

`Bootstrap/` computes bootstrap standard errors and confidence intervals for the AIPW PTE estimator. The current scripts use `B = 500` bootstrap resamples.

`Sandwich/` computes sandwich variance estimates for the AIPW PTE estimator.

`Infulence_function/` computes influence-function-based inference for the AIPW PTE estimator.

## Scenarios

Each method is run separately under four scenarios:

| Scenario | Propensity score model | Outcome model       |
|----------|------------------------|---------------------|
| `cc`     | Correctly specified    | Correctly specified |
| `cm`     | Misspecified           | Correctly specified |
| `mc`     | Correctly specified    | Misspecified        |
| `mm`     | Misspecified           | Misspecified        |

## Required Input Files

Before running a scenario folder, place the following files in that folder:

``` text
nsim1000_data_n1000.RData
true_values.RData
```

`nsim1000_data_n1000.RData` should contain the simulated dataset list, usually named `data_list`, from the point-estimation simulation.

`true_values.RData` should contain the true values table used by the inference scripts. The scripts use the row `q1_0.5` to set:

``` r
truth_p1
truth_p0
truth_pte
y_threshold
```

If the point-estimation simulation produces scenario-specific files such as `nsim1000_data_n1000_cc.RData`, copy or rename the relevant file to `nsim1000_data_n1000.RData` inside the matching inference scenario folder before running `RUN.R`.

The scripts assume these inputs already exist. They do not automatically search other folders or generate the inputs during inference.

## Main Scripts

Each scenario folder contains:

| File | Purpose |
|------------------------------------|------------------------------------|
| `RUN.R` | Main script for the inference simulation. |
| `haojob.sh` | Slurm submission script for cluster runs. |
| `Helper_functions.R` | General helper functions used by the simulation. |
| `PTE_helpers.R` | PTE estimation and inference helper functions. |
| `single_sim.R` | Single-dataset simulation helper. |
| `true_value_generate.R` | Optional script for generating `true_values.RData`. |
