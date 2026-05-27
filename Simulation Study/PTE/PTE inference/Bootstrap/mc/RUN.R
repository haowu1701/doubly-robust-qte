library(Matrix)
library(MASS)
library(rms)
library(dplyr)
library(tidyr)

source("Helper_functions.R")
source("PTE_helpers.R")

# Required input files should be generated before running this script.
load("nsim1000_data_n1000.RData")
load("true_values.RData")

truth_p1  = true_values["q1_0.5", "FY1"]
truth_p0  = true_values["q1_0.5", "FY0"]
truth_pte = true_values["q1_0.5", "PTE"]

y_threshold = true_values["q1_0.5","threshold"]
scenario = "mc"
n_run = length(data_list)
#n_run = 4

test_data = generate_data( n = 1000, alpha1 = 0.5, alpha2 = 0.35,
  beta1 = -2, beta2 = 3, delta = 2 )
dd = datadist(test_data); options(datadist = "dd")

## bootstrap estimates
res_bootstrap = run_bootstrap_AIPW_on_datalist_PTE(
  data_list = data_list,
  y_threshold = y_threshold,
  scenario = scenario,
  B = 500,
  n_run = n_run,
  seed = 123
)

## bootstrap SE and percentile CI for each simulation dataset
boot_ci = res_bootstrap |>
  dplyr::filter(success) |>
  dplyr::group_by(sim) |>
  dplyr::summarise(
    boot_se_Fy1 = sd(AIPW_Fy1_hat, na.rm = TRUE),
    boot_se_Fy0 = sd(AIPW_Fy0_hat, na.rm = TRUE),
    boot_se_pte = sd(AIPW_pte_hat, na.rm = TRUE),
    
    boot_lcl_Fy1 = quantile(AIPW_Fy1_hat, 0.025, na.rm = TRUE),
    boot_ucl_Fy1 = quantile(AIPW_Fy1_hat, 0.975, na.rm = TRUE),
    
    boot_lcl_Fy0 = quantile(AIPW_Fy0_hat, 0.025, na.rm = TRUE),
    boot_ucl_Fy0 = quantile(AIPW_Fy0_hat, 0.975, na.rm = TRUE),
    
    boot_lcl_pte = quantile(AIPW_pte_hat, 0.025, na.rm = TRUE),
    boot_ucl_pte = quantile(AIPW_pte_hat, 0.975, na.rm = TRUE),
    
    n_boot_success = dplyr::n(),
    .groups = "drop"
  )


summary_bootstrap = data.frame(
  estimator = "AIPW",
  parameter = c("p1", "p0", "pte"),
  
  var = c(
    mean(boot_ci$boot_se_Fy1^2, na.rm = TRUE),
    mean(boot_ci$boot_se_Fy0^2, na.rm = TRUE),
    mean(boot_ci$boot_se_pte^2, na.rm = TRUE)
  ),
  
  se = c(
    mean(boot_ci$boot_se_Fy1, na.rm = TRUE),
    mean(boot_ci$boot_se_Fy0, na.rm = TRUE),
    mean(boot_ci$boot_se_pte, na.rm = TRUE)
  ),
  
  coverage = c(
    mean(boot_ci$boot_lcl_Fy1 <= truth_p1  & boot_ci$boot_ucl_Fy1 >= truth_p1, na.rm = TRUE),
    mean(boot_ci$boot_lcl_Fy0 <= truth_p0  & boot_ci$boot_ucl_Fy0 >= truth_p0, na.rm = TRUE),
    mean(boot_ci$boot_lcl_pte <= truth_pte & boot_ci$boot_ucl_pte >= truth_pte, na.rm = TRUE)
  )
)

summary_bootstrap

save(res_bootstrap, summary_bootstrap, file = "Bootstrap_variance_results.RData")
