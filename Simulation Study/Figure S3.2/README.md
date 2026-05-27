# Computational Time Figure

This folder contains the code for generating the computational time comparison figure.

## Files

- `RUN.R`: runs the timing simulation and saves `time_summary.csv`.
- `single_sim.R`: runs one simulation replicate and records elapsed time for each estimator.
- `Helper_functions.R`: contains the data-generating function and estimator helper functions.
- `Figures.R`: reads `time_summary.csv` and generates the figure.
- `time_summary.csv`: timing summary used to generate the figure.

## Run

To regenerate the timing summary:

``` bash
Rscript RUN.R
```

To regenerate the figure from the existing timing summary:

``` bash
Rscript Figures.R
```

These files correspond to Figure S3.2 in the Supplementary Material.
