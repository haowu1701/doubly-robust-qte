# rm(list = ls())
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(rms)
library(dplyr)
library(tidyr)

nsample = 1000; nsim = 1000; scenarios = "cc"

load("true_values.RData")
source("single_sim.R")
source("Helper_functions.R")

test_data = generate_data( n = nsample,alpha1 = 0.5, alpha2 = 0.35, beta1 = -2, beta2 = 3, delta = 2 )
dd = datadist(test_data); options(datadist = "dd")

all_results = list()
result_names = c(
  "aipw_cpm_q_q1", "aipw_cpm_q_q0", "aipw_cpm_q_QTE",
  "aipw_q1", "aipw_q0", "aipw_QTE"
)

for (sce in scenarios) {
  
  results_sce = matrix(NA, nrow = nsim, ncol = length(result_names))
  colnames(results_sce) = result_names
  
  data_list = vector("list", nsim)
  success_seeds = rep(NA, nsim)
  
  set.seed(37203)
  sim_seeds = sample.int(1e7, nsim * 10)
  
  n_success = 0
  i = 1
  
  while (n_success < nsim) {
    
    set.seed(sim_seeds[i])
    
    out = tryCatch(
      {
        single_sim(
          n = nsample,
          p_star = 0.5,
          alpha1 = 0.5, alpha2 = 0.35,
          beta1 = -2, beta2 = 3,
          delta = 2,
          scenario = sce
        )
      },
      error = function(e) {
        cat("Scenario", sce, "- attempt", i, "ERROR, retry...\n")
        return(NULL)
      }
    )
    
    if (!is.null(out)) {
      n_success = n_success + 1
      
      results_sce[n_success, ] = out$result[result_names]
      data_list[[n_success]] = out$data
      success_seeds[n_success] = sim_seeds[i]
      
      if (n_success %% 2 == 0) {
        cat("Scenario:", sce, "- Completed", n_success, "successful simulations\n")
      }
    }
    
    i = i + 1
  }
  
  write.csv(
    results_sce,
    file = paste0("nsim", nsim, "_results_n", nsample, "_", sce, ".csv"),
    row.names = FALSE
  )
  
  save(
    data_list, success_seeds,
    file = paste0("nsim", nsim, "_data_n", nsample, "_", sce, ".RData")
  )
  
  all_results[[sce]] = results_sce
}

##################################
# Summarize Simulation Statistics 
##################################
all_summaries = lapply(scenarios, function(sce) {
  tmp_res = all_results[[sce]]
  
  res_summary = rbind(
    # AIPW_cpm_q
    y1_aipw_cpm_q = calc_stats(tmp_res[, "aipw_cpm_q_q1"], m_Y1_0.5, nsample = nsample),
    y0_aipw_cpm_q = calc_stats(tmp_res[, "aipw_cpm_q_q0"], m_Y0_0.5, nsample = nsample),
    QTE_aipw_cpm_q = calc_stats(tmp_res[, "aipw_cpm_q_QTE"], QTE_0.5_true, nsample = nsample),
    
    # AIPW_cpm_cdf
    y1_aipw = calc_stats(tmp_res[, "aipw_q1"], m_Y1_0.5, nsample = nsample),
    y0_aipw = calc_stats(tmp_res[, "aipw_q0"], m_Y0_0.5, nsample = nsample),
    QTE_aipw = calc_stats(tmp_res[, "aipw_QTE"], QTE_0.5_true, nsample = nsample)
  )
  
  round(res_summary, 3)
})

names(all_summaries) = scenarios

final_table = do.call(rbind, lapply(names(all_summaries), function(sce) {
  tmp = all_summaries[[sce]]
  data.frame(Scenario = sce, Parameter = rownames(tmp), tmp, row.names = NULL)
}))

final_table[, 3:8] = round(final_table[, 3:8], 3)
write.csv(final_table, "simulation_summary_all_scenarios_long.csv", row.names = FALSE)

############################
#Reshape Data: Long to Wide
############################
final_table2 = final_table %>%
  mutate(
    target = case_when(
      grepl("^y1_", Parameter)  ~ "q1",
      grepl("^y0_", Parameter)  ~ "q0",
      grepl("^QTE_", Parameter) ~ "QTE"
    ),
    method_extracted = sub("^(y1_|y0_|QTE_)", "", Parameter),
    method = case_when(
      method_extracted == "aipw"            ~ "AIPW_CPM_icdf",
      method_extracted == "aipw_cpm_q"      ~ "AIPW_CPM"
    ),
    Truth = case_when(
      target == "q1"  ~ m_Y1_0.5,
      target == "q0"  ~ m_Y0_0.5,
      target == "QTE" ~ QTE_0.5_true
    )
  ) %>%
  mutate(method = factor(method, levels = c(
    "AIPW_CPM", "AIPW_CPM_icdf"
  )))

wide_table = final_table2 %>%
  dplyr::select(
    Scenario, method, target, Truth,
    Estimate, Bias, Variance, MSE, RMSE, MAE, MedAE, Nvar
  ) %>%
  pivot_wider(
    names_from = target,
    values_from = c(Truth, Estimate, Bias, Variance, MSE, RMSE, MAE, MedAE, Nvar),
    names_glue = "{target}_{.value}"
  ) %>%
  dplyr::select(Scenario, method, starts_with("QTE_"), starts_with("q1_"), starts_with("q0_")) %>%
  arrange(Scenario, method) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

write.csv(wide_table, "simulation_summary_all_scenarios_wide.csv", row.names = FALSE)