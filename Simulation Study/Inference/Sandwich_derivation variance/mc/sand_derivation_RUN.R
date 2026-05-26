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

p_star= 0.5; scenario = "mc"; n_run  = length(data_list)
test_data = generate_data( n = 1000,alpha1 = 0.5, alpha2 = 0.35, beta1 = -2, beta2 = 3, delta = 2 )
dd = datadist(test_data); options(datadist = "dd")


res_all = run_sandwich_derivation_on_datalist(
  data_list = data_list,
  scenario  = scenario,   
  p_star = p_star,
  n_run = n_run
)



out = summarize_sandwich_results(
  res_df = res_all,
  truth_q1  = truth_q1,
  truth_q0  = truth_q0,
  truth_qte = truth_qte
)
out$summary
save(res_all, file = paste0("sandwich_results_derivation_", scenario, ".RData"))