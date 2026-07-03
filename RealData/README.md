# Real-Data Analyses

This directory contains reproducible real-data analysis code for the CD4 and viral load outcomes. The original patient-level analysis data are private and are not included in this repository. Instead, each analysis directory contains a synthetic dataset that allows the public code to run and illustrates the analysis workflow.

## Directory Structure

- `CD4outcome/`: CD4 outcome analyses.
  - `MQ_PTE_CD4.Rmd`: marginal QTE/PTE analysis for the CD4 outcome.
  - `CQ_PTE_CD4.Rmd`: conditional QTE/PTE analysis by sex for the CD4 outcome.
  - `data_model_cd4_YA_synthetic.RData`: synthetic CD4 analysis dataset.
  - `generate_synthetic_cd4_data.R`: maintainer script used to regenerate the synthetic CD4 dataset from the private local source data.
- `VLoutcome/`: viral load outcome analysis.
  - `MQ_PTE_twomodels_VL.Rmd`: marginal QTE/PTE analysis for the viral load outcome.
  - `data_model_VL_YA_synthetic.RData`: synthetic viral load analysis dataset.
  - `generate_synthetic_vl_data.R`: maintainer script used to regenerate the synthetic viral load dataset from the private local source data.

Each outcome-specific folder also contains a `README.md` with more detail about its synthetic data generation and analysis files.

## Data Availability

The private source datasets are not distributed:

``` text
data_model_cd4_YA.RData
data_model_VL_YA.RData
```

The public R Markdown files use synthetic datasets instead:

``` text
CD4outcome/data_model_cd4_YA_synthetic.RData
VLoutcome/data_model_VL_YA_synthetic.RData
```

## Running the Analyses

From within each outcome directory, render the R Markdown files with:

``` r
rmarkdown::render("MQ_PTE_CD4.Rmd")
rmarkdown::render("CQ_PTE_CD4.Rmd")
```

or, for the viral load outcome:

``` r
rmarkdown::render("MQ_PTE_twomodels_VL.Rmd")
```
