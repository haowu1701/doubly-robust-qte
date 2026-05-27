# rm(list = ls())
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(rms)
library(dplyr)
library(tidyr)

nsample = 1000; nsim = 1000; scenarios = c("cc", "cm", "mc", "mm")

load("true_values.RData")
source("single_sim.R")
source("Helper_functions.R")
source("PTE_helpers.R")

truth_p1  = true_values["q1_0.5","FY1"]
truth_p0  = true_values["q1_0.5","FY0"]
truth_pte = true_values["q1_0.5","PTE"]

test_data = generate_data( n = nsample,alpha1 = 0.5, alpha2 = 0.35, beta1 = -2, beta2 = 3, delta = 2 )
dd = datadist(test_data); options(datadist = "dd")

all_results = list()
result_names = c( "AIPW_Fy1_hat", "AIPW_Fy0_hat", "AIPW_pte_hat",
                  "OR_Fy1_hat", "OR_Fy0_hat", "OR_pte_hat", 
                  "IPW_Fy1_hat", "IPW_Fy0_hat", "IPW_pte_hat")

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
          y_threshold = true_values["q1_0.5","threshold"],
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
    # AIPW
    p1_AIPW = calc_stats(tmp_res[, "AIPW_Fy1_hat"], truth_p1,  nsample = nsample),
    p0_AIPW = calc_stats(tmp_res[, "AIPW_Fy0_hat"], truth_p0,  nsample = nsample),
    PTE_AIPW = calc_stats(tmp_res[, "AIPW_pte_hat"], truth_pte, nsample = nsample),
    
    # OR
    p1_OR = calc_stats(tmp_res[, "OR_Fy1_hat"], truth_p1,  nsample = nsample),
    p0_OR = calc_stats(tmp_res[, "OR_Fy0_hat"], truth_p0,  nsample = nsample),
    PTE_OR = calc_stats(tmp_res[, "OR_pte_hat"], truth_pte, nsample = nsample),
    
    # IPW
    p1_IPW = calc_stats(tmp_res[, "IPW_Fy1_hat"], truth_p1,  nsample = nsample),
    p0_IPW = calc_stats(tmp_res[, "IPW_Fy0_hat"], truth_p0,  nsample = nsample),
    PTE_IPW = calc_stats(tmp_res[, "IPW_pte_hat"], truth_pte, nsample = nsample)
  )
  
  res_summary
})

names(all_summaries) = scenarios

final_table = do.call(rbind, lapply(names(all_summaries), function(sce) {
  tmp = all_summaries[[sce]]
  data.frame(Scenario = sce, Parameter = rownames(tmp), tmp, row.names = NULL)
}))

final_table[, 3:ncol(final_table)] = final_table[, 3:ncol(final_table)]

write.csv(
  final_table,
  "PTE_simulation_summary_all_scenarios_long.csv",
  row.names = FALSE
)

############################
# Reshape Data: Long to Wide
############################

final_table2 = final_table %>%
  mutate(
    target = case_when(
      grepl("^p1_", Parameter)  ~ "p1",
      grepl("^p0_", Parameter)  ~ "p0",
      grepl("^PTE_", Parameter) ~ "PTE"
    ),
    method = sub("^(p1_|p0_|PTE_)", "", Parameter),
    Truth = case_when(
      target == "p1"  ~ truth_p1,
      target == "p0"  ~ truth_p0,
      target == "PTE" ~ truth_pte
    )
  ) %>%
  mutate(
    method = factor(method, levels = c("AIPW", "OR", "IPW"))
  )

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
  dplyr::select(
    Scenario, method,
    starts_with("PTE_"),
    starts_with("p1_"),
    starts_with("p0_")
  ) %>%
  arrange(Scenario, method) #%>%
 # mutate(across(where(is.numeric), ~ round(., 3)))

write.csv(
  wide_table,
  "PTE_simulation_summary_all_scenarios_wide.csv",
  row.names = FALSE
)
