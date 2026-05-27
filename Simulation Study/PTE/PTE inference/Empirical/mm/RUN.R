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
#n_run = 10

test_data = generate_data( n = 1000, alpha1 = 0.5, alpha2 = 0.35,
  beta1 = -2, beta2 = 3, delta = 2 )
dd = datadist(test_data); options(datadist = "dd")

res_empirical = run_empirical_on_datalist_PTE(
  data_list = data_list,
  y_threshold = y_threshold,
  scenario = scenario,
  n_run = n_run
)


## empirical SE based on simulation SD
se_AIPW_Fy1 = sd(res_empirical$AIPW_Fy1_hat, na.rm = TRUE)
se_AIPW_Fy0 = sd(res_empirical$AIPW_Fy0_hat, na.rm = TRUE)
se_AIPW_pte = sd(res_empirical$AIPW_pte_hat, na.rm = TRUE)

se_OR_Fy1 = sd(res_empirical$OR_Fy1_hat, na.rm = TRUE)
se_OR_Fy0 = sd(res_empirical$OR_Fy0_hat, na.rm = TRUE)
se_OR_pte = sd(res_empirical$OR_pte_hat, na.rm = TRUE)

se_IPW_Fy1 = sd(res_empirical$IPW_Fy1_hat, na.rm = TRUE)
se_IPW_Fy0 = sd(res_empirical$IPW_Fy0_hat, na.rm = TRUE)
se_IPW_pte = sd(res_empirical$IPW_pte_hat, na.rm = TRUE)


## empirical CI coverage
res_empirical$AIPW_cover_Fy1 = with(
  res_empirical,
  AIPW_Fy1_hat - 1.96 * se_AIPW_Fy1 <= truth_p1 &
    AIPW_Fy1_hat + 1.96 * se_AIPW_Fy1 >= truth_p1
)

res_empirical$AIPW_cover_Fy0 = with(
  res_empirical,
  AIPW_Fy0_hat - 1.96 * se_AIPW_Fy0 <= truth_p0 &
    AIPW_Fy0_hat + 1.96 * se_AIPW_Fy0 >= truth_p0
)

res_empirical$AIPW_cover_pte = with(
  res_empirical,
  AIPW_pte_hat - 1.96 * se_AIPW_pte <= truth_pte &
    AIPW_pte_hat + 1.96 * se_AIPW_pte >= truth_pte
)


res_empirical$OR_cover_Fy1 = with(
  res_empirical,
  OR_Fy1_hat - 1.96 * se_OR_Fy1 <= truth_p1 &
    OR_Fy1_hat + 1.96 * se_OR_Fy1 >= truth_p1
)

res_empirical$OR_cover_Fy0 = with(
  res_empirical,
  OR_Fy0_hat - 1.96 * se_OR_Fy0 <= truth_p0 &
    OR_Fy0_hat + 1.96 * se_OR_Fy0 >= truth_p0
)

res_empirical$OR_cover_pte = with(
  res_empirical,
  OR_pte_hat - 1.96 * se_OR_pte <= truth_pte &
    OR_pte_hat + 1.96 * se_OR_pte >= truth_pte
)


res_empirical$IPW_cover_Fy1 = with(
  res_empirical,
  IPW_Fy1_hat - 1.96 * se_IPW_Fy1 <= truth_p1 &
    IPW_Fy1_hat + 1.96 * se_IPW_Fy1 >= truth_p1
)

res_empirical$IPW_cover_Fy0 = with(
  res_empirical,
  IPW_Fy0_hat - 1.96 * se_IPW_Fy0 <= truth_p0 &
    IPW_Fy0_hat + 1.96 * se_IPW_Fy0 >= truth_p0
)

res_empirical$IPW_cover_pte = with(
  res_empirical,
  IPW_pte_hat - 1.96 * se_IPW_pte <= truth_pte &
    IPW_pte_hat + 1.96 * se_IPW_pte >= truth_pte
)


summary_empirical = data.frame(
  estimator = rep(c("AIPW", "OR", "IPW"), each = 3),
  parameter = rep(c("p1", "p0", "pte"), times = 3),
  
  true = c(
    truth_p1, truth_p0, truth_pte,
    truth_p1, truth_p0, truth_pte,
    truth_p1, truth_p0, truth_pte
  ),
  
  mean_est = c(
    mean(res_empirical$AIPW_Fy1_hat, na.rm = TRUE),
    mean(res_empirical$AIPW_Fy0_hat, na.rm = TRUE),
    mean(res_empirical$AIPW_pte_hat, na.rm = TRUE),
    
    mean(res_empirical$OR_Fy1_hat, na.rm = TRUE),
    mean(res_empirical$OR_Fy0_hat, na.rm = TRUE),
    mean(res_empirical$OR_pte_hat, na.rm = TRUE),
    
    mean(res_empirical$IPW_Fy1_hat, na.rm = TRUE),
    mean(res_empirical$IPW_Fy0_hat, na.rm = TRUE),
    mean(res_empirical$IPW_pte_hat, na.rm = TRUE)
  ),
  
  bias = c(
    mean(res_empirical$AIPW_Fy1_hat - truth_p1, na.rm = TRUE),
    mean(res_empirical$AIPW_Fy0_hat - truth_p0, na.rm = TRUE),
    mean(res_empirical$AIPW_pte_hat - truth_pte, na.rm = TRUE),
    
    mean(res_empirical$OR_Fy1_hat - truth_p1, na.rm = TRUE),
    mean(res_empirical$OR_Fy0_hat - truth_p0, na.rm = TRUE),
    mean(res_empirical$OR_pte_hat - truth_pte, na.rm = TRUE),
    
    mean(res_empirical$IPW_Fy1_hat - truth_p1, na.rm = TRUE),
    mean(res_empirical$IPW_Fy0_hat - truth_p0, na.rm = TRUE),
    mean(res_empirical$IPW_pte_hat - truth_pte, na.rm = TRUE)
  ),
  
  emp_sd = c(
    sd(res_empirical$AIPW_Fy1_hat, na.rm = TRUE),
    sd(res_empirical$AIPW_Fy0_hat, na.rm = TRUE),
    sd(res_empirical$AIPW_pte_hat, na.rm = TRUE),
    
    sd(res_empirical$OR_Fy1_hat, na.rm = TRUE),
    sd(res_empirical$OR_Fy0_hat, na.rm = TRUE),
    sd(res_empirical$OR_pte_hat, na.rm = TRUE),
    
    sd(res_empirical$IPW_Fy1_hat, na.rm = TRUE),
    sd(res_empirical$IPW_Fy0_hat, na.rm = TRUE),
    sd(res_empirical$IPW_pte_hat, na.rm = TRUE)
  ),
  
  emp_var = c(
    var(res_empirical$AIPW_Fy1_hat, na.rm = TRUE),
    var(res_empirical$AIPW_Fy0_hat, na.rm = TRUE),
    var(res_empirical$AIPW_pte_hat, na.rm = TRUE),
    
    var(res_empirical$OR_Fy1_hat, na.rm = TRUE),
    var(res_empirical$OR_Fy0_hat, na.rm = TRUE),
    var(res_empirical$OR_pte_hat, na.rm = TRUE),
    
    var(res_empirical$IPW_Fy1_hat, na.rm = TRUE),
    var(res_empirical$IPW_Fy0_hat, na.rm = TRUE),
    var(res_empirical$IPW_pte_hat, na.rm = TRUE)
  ),
  
  coverage = c(
    mean(res_empirical$AIPW_cover_Fy1, na.rm = TRUE),
    mean(res_empirical$AIPW_cover_Fy0, na.rm = TRUE),
    mean(res_empirical$AIPW_cover_pte, na.rm = TRUE),
    
    mean(res_empirical$OR_cover_Fy1, na.rm = TRUE),
    mean(res_empirical$OR_cover_Fy0, na.rm = TRUE),
    mean(res_empirical$OR_cover_pte, na.rm = TRUE),
    
    mean(res_empirical$IPW_cover_Fy1, na.rm = TRUE),
    mean(res_empirical$IPW_cover_Fy0, na.rm = TRUE),
    mean(res_empirical$IPW_cover_pte, na.rm = TRUE)
  )
)

summary_empirical
summary_empirical[which(summary_empirical$parameter == "pte"),]


save(res_empirical, summary_empirical, file = "Empirical_variance_results.RData")
