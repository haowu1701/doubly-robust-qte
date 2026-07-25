library(Matrix)
library(MASS)
library(rms)
library(dplyr)
library(tidyr)

source("Helper_functions.R")
source("variance_helpers.R")
load("nsim1000_data_n1000.RData")
load("true_values.RData")
truth_q1  = m_Y1_0.5
truth_q0  = m_Y0_0.5
truth_qte = QTE_0.5_true

p_star= 0.5; scenario = "cm"; n_run  = length(data_list)

test_data = generate_data( n = 1000,alpha1 = 0.5, alpha2 = 0.35, beta1 = -2, beta2 = 3, delta = 2 )
dd = datadist(test_data); options(datadist = "dd")

res = run_IF_variance_on_datalist(
  data_list = data_list,
  scenario  = scenario,
  p_star    = p_star,
  n_run     = n_run
)

res$cover_q1  = with(res, ci_q1_l  <= truth_q1  & ci_q1_u  >= truth_q1)
res$cover_q0  = with(res, ci_q0_l  <= truth_q0  & ci_q0_u  >= truth_q0)
res$cover_qte = with(res, ci_qte_l <= truth_qte & ci_qte_u >= truth_qte)



## overall summary
summary_if = data.frame(
  parameter = c("q1", "q0", "qte"),
  mean_est = c(
    mean(res$q1_hat,  na.rm = TRUE),
    mean(res$q0_hat,  na.rm = TRUE),
    mean(res$qte_hat, na.rm = TRUE)
  ),
  emp_sd = c(
    sd(res$q1_hat,  na.rm = TRUE),
    sd(res$q0_hat,  na.rm = TRUE),
    sd(res$qte_hat, na.rm = TRUE)
  ),
  mean_var_if = c(
    mean(res$var_if_q1,  na.rm = TRUE),
    mean(res$var_if_q0,  na.rm = TRUE),
    mean(res$var_if_qte, na.rm = TRUE)
  ),
  mean_se_if = c(
    mean(res$se_if_q1,  na.rm = TRUE),
    mean(res$se_if_q0,  na.rm = TRUE),
    mean(res$se_if_qte, na.rm = TRUE)
  ),
  coverage = c(
    mean(res$cover_q1,  na.rm = TRUE),
    mean(res$cover_q0,  na.rm = TRUE),
    mean(res$cover_qte, na.rm = TRUE)
  )
)

save(res, summary_if, file = "IF_variance_results_p0.5.RData")
