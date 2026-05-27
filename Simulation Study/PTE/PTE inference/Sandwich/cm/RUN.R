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
scenario = "cm"
n_run = length(data_list)
#n_run = 4

test_data = generate_data( n = 1000, alpha1 = 0.5, alpha2 = 0.35,
  beta1 = -2, beta2 = 3, delta = 2 )
dd = datadist(test_data); options(datadist = "dd")

res_sandwich = run_sandwich_AIPW_on_datalist_PTE(
  data_list = data_list,
  y_threshold = y_threshold,
  scenario = scenario,
  n_run = n_run
)


res_sandwich = res_sandwich |>
  dplyr::mutate(
    cover_Fy1 = AIPW_Fy1_hat - 1.96 * se_sand_Fy1 <= truth_p1 &
      AIPW_Fy1_hat + 1.96 * se_sand_Fy1 >= truth_p1,
    
    cover_Fy0 = AIPW_Fy0_hat - 1.96 * se_sand_Fy0 <= truth_p0 &
      AIPW_Fy0_hat + 1.96 * se_sand_Fy0 >= truth_p0,
    
    cover_pte = AIPW_pte_hat - 1.96 * se_sand_pte <= truth_pte &
      AIPW_pte_hat + 1.96 * se_sand_pte >= truth_pte
  )

summary_sandwich = data.frame(
  estimator = "AIPW",
  parameter = c("p1", "p0", "pte"),
  
  true = c(truth_p1, truth_p0, truth_pte),
  
  mean_est = c(
    mean(res_sandwich$AIPW_Fy1_hat, na.rm = TRUE),
    mean(res_sandwich$AIPW_Fy0_hat, na.rm = TRUE),
    mean(res_sandwich$AIPW_pte_hat, na.rm = TRUE)
  ),
  
  bias = c(
    mean(res_sandwich$AIPW_Fy1_hat - truth_p1, na.rm = TRUE),
    mean(res_sandwich$AIPW_Fy0_hat - truth_p0, na.rm = TRUE),
    mean(res_sandwich$AIPW_pte_hat - truth_pte, na.rm = TRUE)
  ),
  
  emp_var = c(
    var(res_sandwich$AIPW_Fy1_hat, na.rm = TRUE),
    var(res_sandwich$AIPW_Fy0_hat, na.rm = TRUE),
    var(res_sandwich$AIPW_pte_hat, na.rm = TRUE)
  ),
  
  emp_sd = c(
    sd(res_sandwich$AIPW_Fy1_hat, na.rm = TRUE),
    sd(res_sandwich$AIPW_Fy0_hat, na.rm = TRUE),
    sd(res_sandwich$AIPW_pte_hat, na.rm = TRUE)
  ),
  
  sand_var = c(
    mean(res_sandwich$var_sand_Fy1, na.rm = TRUE),
    mean(res_sandwich$var_sand_Fy0, na.rm = TRUE),
    mean(res_sandwich$var_sand_pte, na.rm = TRUE)
  ),
  
  sand_sd = c(
    mean(res_sandwich$se_sand_Fy1, na.rm = TRUE),
    mean(res_sandwich$se_sand_Fy0, na.rm = TRUE),
    mean(res_sandwich$se_sand_pte, na.rm = TRUE)
  ),
  
  coverage = c(
    mean(res_sandwich$cover_Fy1, na.rm = TRUE),
    mean(res_sandwich$cover_Fy0, na.rm = TRUE),
    mean(res_sandwich$cover_pte, na.rm = TRUE)
  )
)

summary_sandwich
summary_sandwich[summary_sandwich$parameter == "pte", ]


save(res_sandwich, summary_sandwich, file = "Sandwich_variance_results.RData")
