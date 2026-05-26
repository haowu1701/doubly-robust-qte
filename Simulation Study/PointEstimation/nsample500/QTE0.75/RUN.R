# rm(list = ls())
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(rms)
library(dplyr)
library(tidyr)

# Set parameters and scenarios
nsample = 500 ; nsim = 1000
 scenarios = c("cc", "cm", "mc", "mm")
#scenarios = "cm"

load("true_values.RData")
source("single_sim.R")
source("Helper_functions.R")


test_data = generate_data( n = nsample, alpha1 = 0.5, alpha2 = 0.35, beta1 = -2, beta2 = 3, delta = 2 )
dd = datadist(test_data)
options(datadist = "dd")

#####################
# Run Simulations
####################
all_results = list()

result_names = c(
  "aipw_q_q1", "aipw_q_q0", "aipw_q_QTE",
  "tmle_q1", "tmle_q0", "tmle_QTE",
  "rq_tmle_q1", "rq_tmle_q0", "rq_tmle_QTE",
  "cpm_tmle_q1", "cpm_tmle_q0", "cpm_tmle_QTE",
  "aipw_cpm_q_q1", "aipw_cpm_q_q0", "aipw_cpm_q_QTE",
  "aipw_q1", "aipw_q0", "aipw_QTE",
  "g_q1", "g_q0", "g_QTE",
  "iptw_q1", "iptw_q0", "iptw_QTE",
  "firpo_q1", "firpo_q0", "firpo_QTE",
  "exp_tmle_lm_q1", "exp_tmle_lm_q0", "exp_tmle_lm_QTE",
  "exp_tmle_rq_q1", "exp_tmle_rq_q0", "exp_tmle_rq_QTE",
  "exp_aipw_lm_q1", "exp_aipw_lm_q0", "exp_aipw_lm_QTE",
  "tc_tmle_lm_q1", "tc_tmle_lm_q0", "tc_tmle_lm_QTE",
  "tc_tmle_rq_q1", "tc_tmle_rq_q0", "tc_tmle_rq_QTE",
  "tc_aipw_lm_q1", "tc_aipw_lm_q0", "tc_aipw_lm_QTE",
  "aipw_q1_mis", "aipw_q0_mis", "aipw_QTE_mis",
  "aipw_cpm_q_q1_mis", "aipw_cpm_q_q0_mis", "aipw_cpm_q_QTE_mis",
  "cpm_tmle_q1_mis", "cpm_tmle_q0_mis", "cpm_tmle_QTE_mis"
)

for (sce in scenarios) {
  
  results_sce = matrix(NA, nrow = nsim, ncol = length(result_names))
  colnames(results_sce) = result_names
  
  set.seed(37203)
  sim_seeds = sample.int(1e7, nsim * 10)
  
  n_success = 0
  i = 1
  
  while (n_success < nsim) {
    if (i > length(sim_seeds)) {
      stop("Exceeded maximum retry attempts before reaching nsim successful simulations")
    }
    
    set.seed(sim_seeds[i])
    
    out = tryCatch(
      {
        single_sim(
          n = nsample,
          p_star = 0.75,
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
      results_sce[n_success, ] = out[result_names]
      
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
  
  all_results[[sce]] = results_sce
}

##################################
# Summarize Simulation Statistics 
##################################
all_summaries = lapply(scenarios, function(sce) {
  tmp_res = all_results[[sce]]
  
  res_summary = rbind(
    # AIPW_normality_q
    y1_aipw_q = calc_stats(tmp_res[, "aipw_q_q1"], m_Y1_0.75, nsample = nsample),
    y0_aipw_q = calc_stats(tmp_res[, "aipw_q_q0"], m_Y0_0.75, nsample = nsample),
    QTE_aipw_q = calc_stats(tmp_res[, "aipw_q_QTE"], QTE_0.75_true, nsample = nsample),
    
    # TMLE_normality_q
    y1_tmle = calc_stats(tmp_res[, "tmle_q1"], m_Y1_0.75, nsample = nsample),
    y0_tmle = calc_stats(tmp_res[, "tmle_q0"], m_Y0_0.75, nsample = nsample),
    QTE_tmle = calc_stats(tmp_res[, "tmle_QTE"], QTE_0.75_true, nsample = nsample),
    
    # TMLE_cpm_q
    y1_tmle_cpm = calc_stats(tmp_res[, "cpm_tmle_q1"], m_Y1_0.75, nsample = nsample),
    y0_tmle_cpm = calc_stats(tmp_res[, "cpm_tmle_q0"], m_Y0_0.75, nsample = nsample),
    QTE_tmle_cpm = calc_stats(tmp_res[, "cpm_tmle_QTE"], QTE_0.75_true, nsample = nsample),
    
    # TMLE_rq_q
    y1_tmle_rq = calc_stats(tmp_res[, "rq_tmle_q1"], m_Y1_0.75, nsample = nsample),
    y0_tmle_rq = calc_stats(tmp_res[, "rq_tmle_q0"], m_Y0_0.75, nsample = nsample),
    QTE_tmle_rq = calc_stats(tmp_res[, "rq_tmle_QTE"], QTE_0.75_true, nsample = nsample),
    
    # AIPW_cpm_q
    y1_aipw_cpm_q = calc_stats(tmp_res[, "aipw_cpm_q_q1"], m_Y1_0.75, nsample = nsample),
    y0_aipw_cpm_q = calc_stats(tmp_res[, "aipw_cpm_q_q0"], m_Y0_0.75, nsample = nsample),
    QTE_aipw_cpm_q = calc_stats(tmp_res[, "aipw_cpm_q_QTE"], QTE_0.75_true, nsample = nsample),
    
    # AIPW_cpm_cdf
    y1_aipw = calc_stats(tmp_res[, "aipw_q1"], m_Y1_0.75, nsample = nsample),
    y0_aipw = calc_stats(tmp_res[, "aipw_q0"], m_Y0_0.75, nsample = nsample),
    QTE_aipw = calc_stats(tmp_res[, "aipw_QTE"], QTE_0.75_true, nsample = nsample),
    
    # G-cpm_cdf
    y1_g = calc_stats(tmp_res[, "g_q1"], m_Y1_0.75, nsample = nsample),
    y0_g = calc_stats(tmp_res[, "g_q0"], m_Y0_0.75, nsample = nsample),
    QTE_g = calc_stats(tmp_res[, "g_QTE"], QTE_0.75_true, nsample = nsample),
    
    # IPW_cdf
    y1_iptw = calc_stats(tmp_res[, "iptw_q1"], m_Y1_0.75, nsample = nsample),
    y0_iptw = calc_stats(tmp_res[, "iptw_q0"], m_Y0_0.75, nsample = nsample),
    QTE_iptw = calc_stats(tmp_res[, "iptw_QTE"], QTE_0.75_true, nsample = nsample),
    
    # IPW_Firpo_q
    y1_firpo = calc_stats(tmp_res[, "firpo_q1"], m_Y1_0.75, nsample = nsample),
    y0_firpo = calc_stats(tmp_res[, "firpo_q0"], m_Y0_0.75, nsample = nsample),
    QTE_firpo = calc_stats(tmp_res[, "firpo_QTE"], QTE_0.75_true, nsample = nsample),
    
    # exp_TMLE_normality_q
    y1_tmle_lm_exp = calc_stats(tmp_res[, "exp_tmle_lm_q1"], m_Y1_0.75, nsample = nsample),
    y0_tmle_lm_exp = calc_stats(tmp_res[, "exp_tmle_lm_q0"], m_Y0_0.75, nsample = nsample),
    QTE_tmle_lm_exp = calc_stats(tmp_res[, "exp_tmle_lm_QTE"], QTE_0.75_true, nsample = nsample),
    
    # exp_TMLE_rq_q
    y1_tmle_rq_exp = calc_stats(tmp_res[, "exp_tmle_rq_q1"], m_Y1_0.75, nsample = nsample),
    y0_tmle_rq_exp = calc_stats(tmp_res[, "exp_tmle_rq_q0"], m_Y0_0.75, nsample = nsample),
    QTE_tmle_rq_exp = calc_stats(tmp_res[, "exp_tmle_rq_QTE"], QTE_0.75_true, nsample = nsample),
    
    # exp_AIPW_lm_q
    y1_aipw_lm_exp = calc_stats(tmp_res[, "exp_aipw_lm_q1"], m_Y1_0.75, nsample = nsample),
    y0_aipw_lm_exp = calc_stats(tmp_res[, "exp_aipw_lm_q0"], m_Y0_0.75, nsample = nsample),
    QTE_aipw_lm_exp = calc_stats(tmp_res[, "exp_aipw_lm_QTE"], QTE_0.75_true, nsample = nsample),
    
    # tc_TMLE_normality_q
    y1_tmle_lm_tc = calc_stats(tmp_res[, "tc_tmle_lm_q1"], m_Y1_0.75, nsample = nsample),
    y0_tmle_lm_tc = calc_stats(tmp_res[, "tc_tmle_lm_q0"], m_Y0_0.75, nsample = nsample),
    QTE_tmle_lm_tc = calc_stats(tmp_res[, "tc_tmle_lm_QTE"], QTE_0.75_true, nsample = nsample),
    
    # tc_TMLE_rq_q
    y1_tmle_rq_tc = calc_stats(tmp_res[, "tc_tmle_rq_q1"], m_Y1_0.75, nsample = nsample),
    y0_tmle_rq_tc = calc_stats(tmp_res[, "tc_tmle_rq_q0"], m_Y0_0.75, nsample = nsample),
    QTE_tmle_rq_tc = calc_stats(tmp_res[, "tc_tmle_rq_QTE"], QTE_0.75_true, nsample = nsample),
    
    # tc_AIPW_lm_q
    y1_aipw_lm_tc = calc_stats(tmp_res[, "tc_aipw_lm_q1"], m_Y1_0.75, nsample = nsample),
    y0_aipw_lm_tc = calc_stats(tmp_res[, "tc_aipw_lm_q0"], m_Y0_0.75, nsample = nsample),
    QTE_aipw_lm_tc = calc_stats(tmp_res[, "tc_aipw_lm_QTE"], QTE_0.75_true, nsample = nsample),
    
    # mis-link CPM estimators
    y1_aipw_mis = calc_stats(tmp_res[, "aipw_q1_mis"], m_Y1_0.75, nsample = nsample),
    y0_aipw_mis = calc_stats(tmp_res[, "aipw_q0_mis"], m_Y0_0.75, nsample = nsample),
    QTE_aipw_mis = calc_stats(tmp_res[, "aipw_QTE_mis"], QTE_0.75_true, nsample = nsample),
    
    y1_aipw_cpm_q_mis = calc_stats(tmp_res[, "aipw_cpm_q_q1_mis"], m_Y1_0.75, nsample = nsample),
    y0_aipw_cpm_q_mis = calc_stats(tmp_res[, "aipw_cpm_q_q0_mis"], m_Y0_0.75, nsample = nsample),
    QTE_aipw_cpm_q_mis = calc_stats(tmp_res[, "aipw_cpm_q_QTE_mis"], QTE_0.75_true, nsample = nsample),
    
    y1_tmle_cpm_mis = calc_stats(tmp_res[, "cpm_tmle_q1_mis"], m_Y1_0.75, nsample = nsample),
    y0_tmle_cpm_mis = calc_stats(tmp_res[, "cpm_tmle_q0_mis"], m_Y0_0.75, nsample = nsample),
    QTE_tmle_cpm_mis = calc_stats(tmp_res[, "cpm_tmle_QTE_mis"], QTE_0.75_true, nsample = nsample)
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
      method_extracted == "iptw"            ~ "IPW_icdf",
      method_extracted == "g"               ~ "OR_CPM_cdf",
      method_extracted == "firpo"           ~ "IPW_Firpo_q",
      method_extracted == "aipw_q"          ~ "AIPW",
      method_extracted == "tmle"            ~ "TMLE",
      method_extracted == "aipw_cpm_q"      ~ "AIPW_CPM",
      method_extracted == "tmle_cpm"        ~ "TMLE_CPM",
      method_extracted == "tmle_rq"         ~ "TMLE_cqr",
      method_extracted == "tmle_lm_exp"     ~ "log_TMLE",
      method_extracted == "tmle_rq_exp"     ~ "log_TMLE_cqr",
      method_extracted == "aipw_lm_exp"     ~ "log_AIPW",
      method_extracted == "tmle_lm_tc"      ~ "ct_TMLE",
      method_extracted == "tmle_rq_tc"      ~ "ct_TMLE_cqr",
      method_extracted == "aipw_lm_tc"      ~ "ct_AIPW",
      method_extracted == "aipw_mis"        ~ "AIPW_CPM_icdf_mislink",
      method_extracted == "aipw_cpm_q_mis"  ~ "AIPW_CPM_mislink",
      method_extracted == "tmle_cpm_mis"    ~ "TMLE_CPM_mislink"
    ),
    Truth = case_when(
      target == "q1"  ~ m_Y1_0.75,
      target == "q0"  ~ m_Y0_0.75,
      target == "QTE" ~ QTE_0.75_true
    )
  ) %>%
  mutate(method = factor(method, levels = c(
    "AIPW", "TMLE","TMLE_cqr",
    "log_AIPW", "log_TMLE", "log_TMLE_cqr",
    "ct_AIPW", "ct_TMLE", "ct_TMLE_cqr",
    "AIPW_CPM", "AIPW_CPM_icdf", "TMLE_CPM",
    "AIPW_CPM_mislink", "AIPW_CPM_icdf_mislink",  "TMLE_CPM_mislink",
    "OR_CPM_cdf", "IPW_icdf", "IPW_Firpo_q"
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
