# Viral Load Outcome Real-Data Analysis

This directory contains reproducible viral load outcome analysis code using a synthetic dataset. The synthetic data are included so that the public code can be run without distributing private patient-level records.

## Files

- `MQ_PTE_twomodels_VL.Rmd`: marginal viral load outcome analysis. This file estimates marginal quantile treatment effects (QTE) and probability treatment effects (PTE) using AIPW-CPM approaches.
- `Helper_functions.R`: helper functions used by the R Markdown file, including quantile inversion and conditional CDF evaluation from two outcome models.
- `data_model_VL_YA_synthetic.RData`: synthetic analysis dataset used by the R Markdown file. It contains an object named `data_model`.
- `generate_synthetic_vl_data.R`: maintainer script for regenerating the synthetic dataset from the private local data file.

## Data Availability

The original analysis data file, `data_model_VL_YA.RData`, contains real patient-level data and is not included in this repository.

The public R Markdown file loads:

``` r
data_model_VL_YA_synthetic.RData
```

The synthetic file contains only the variables used by the viral load analysis:

``` r
c(
  "Y",
  "A",
  "age",
  "sex",
  "site",
  "route_infection",
  "prior_aids",
  "calendar_year",
  "months_toVLmeasure",
  "cd4_baseline_tran",
  "singledl_rna_baseline_v_tran"
)
```

`generate_synthetic_vl_data.R` requires the private local file `data_model_VL_YA.RData`. Public users do not need to run it unless they have access to the private source data.

## Synthetic Data Generation

The continuous covariates `age`, `months_toVLmeasure`, `cd4_baseline_tran`, and `singledl_rna_baseline_v_tran` are generated jointly from a multivariate normal approximation using the private data's mean vector and covariance matrix, then clipped to the 1st-99th percentile range.

The binary covariates `sex` and `prior_aids` are generated from Bernoulli distributions using empirical event probabilities from the private data.

The categorical covariates `site`, `route_infection`, and `calendar_year` are sampled from empirical category frequencies with a small `+1` smoothing term. `calendar_year` is sampled as a discrete year.

The treatment variable `A` is generated from a fitted logistic propensity model. Predicted treatment probabilities are clipped to `[0.05, 0.95]`, then used for Bernoulli draws.

The viral load outcome `Y` is generated with a two-part model to preserve the strong point mass at the lower detection limit. The script first models whether `Y` equals the observed lower limit, then generates values above the limit from a log-scale outcome regression with residual noise.

## Important Note

The synthetic data do not contain real patient records. Numerical results, fitted models, tables, and figures generated from the synthetic data are for code demonstration only and should not be interpreted as study findings.
