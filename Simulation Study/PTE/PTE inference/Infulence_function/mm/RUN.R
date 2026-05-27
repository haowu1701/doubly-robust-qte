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
scenario = "mm"
n_run = length(data_list)
#n_run = 100

test_data = generate_data( n = 1000, alpha1 = 0.5, alpha2 = 0.35,
  beta1 = -2, beta2 = 3, delta = 2 )
dd = datadist(test_data); options(datadist = "dd")

res_if = run_IF_AIPW_on_datalist_PTE(
  data_list = data_list,
  y_threshold = y_threshold,
  scenario = scenario,
  link = "probit",
  n_run = n_run
)

res_if = res_if |>
  dplyr::mutate(
    cover_Fy1 = ci_Fy1_l <= truth_p1 & ci_Fy1_u >= truth_p1,
    cover_Fy0 = ci_Fy0_l <= truth_p0 & ci_Fy0_u >= truth_p0,
    cover_pte = ci_pte_l <= truth_pte & ci_pte_u >= truth_pte
  )

summary_if = data.frame(
  estimator = "AIPW",
  parameter = c("p1", "p0", "pte"),
  
  true = c(truth_p1, truth_p0, truth_pte),
  
  mean_est = c(
    mean(res_if$AIPW_Fy1_hat, na.rm = TRUE),
    mean(res_if$AIPW_Fy0_hat, na.rm = TRUE),
    mean(res_if$AIPW_pte_hat, na.rm = TRUE)
  ),
  
  bias = c(
    mean(res_if$AIPW_Fy1_hat - truth_p1, na.rm = TRUE),
    mean(res_if$AIPW_Fy0_hat - truth_p0, na.rm = TRUE),
    mean(res_if$AIPW_pte_hat - truth_pte, na.rm = TRUE)
  ),
  
  emp_var = c(
    var(res_if$AIPW_Fy1_hat, na.rm = TRUE),
    var(res_if$AIPW_Fy0_hat, na.rm = TRUE),
    var(res_if$AIPW_pte_hat, na.rm = TRUE)
  ),
  
  emp_sd = c(
    sd(res_if$AIPW_Fy1_hat, na.rm = TRUE),
    sd(res_if$AIPW_Fy0_hat, na.rm = TRUE),
    sd(res_if$AIPW_pte_hat, na.rm = TRUE)
  ),
  
  if_var = c(
    mean(res_if$var_if_Fy1, na.rm = TRUE),
    mean(res_if$var_if_Fy0, na.rm = TRUE),
    mean(res_if$var_if_pte, na.rm = TRUE)
  ),
  
  if_sd = c(
    mean(res_if$se_if_Fy1, na.rm = TRUE),
    mean(res_if$se_if_Fy0, na.rm = TRUE),
    mean(res_if$se_if_pte, na.rm = TRUE)
  ),
  
  coverage = c(
    mean(res_if$cover_Fy1, na.rm = TRUE),
    mean(res_if$cover_Fy0, na.rm = TRUE),
    mean(res_if$cover_pte, na.rm = TRUE)
  )
)

summary_if
summary_if[summary_if$parameter == "pte", ]


save(res_if, summary_if, file = "IF_AIPW_PTE_results.RData")
