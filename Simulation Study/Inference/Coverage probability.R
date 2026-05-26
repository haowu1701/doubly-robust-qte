############################
# Coverage Probability
############################

## Choose the scenario to summarize: cc, cm, mc, or mm.
scenario = "cc"
p_star = "0.5"

true_values_file = file.path(
  "Empirical variance", scenario, "true_values.RData"
)
empirical_result_file = file.path(
  "Empirical variance", scenario,
  paste0("nsim1000_results_n1000_", scenario, ".csv")
)
sandwich_result_file = file.path(
  "sandwich_expectation", scenario,
  paste0("sandwich_results_expetation_", scenario, ".RData")
)
if_result_file = file.path(
  "IF_based var", scenario, "IF_variance_results_p0.5.RData"
)
bootstrap_result_file = file.path(
  "BootstrapnBoot500", scenario,
  paste0("boot_results_all_p0.5_", scenario, ".rds")
)

if (file.exists(true_values_file)) {
  load(true_values_file)
} else {
  stop("Missing true-value file: ", true_values_file)
}

truth_q1 = get(paste0("m_Y1_", p_star))
truth_q0 = get(paste0("m_Y0_", p_star))
truth_qte = get(paste0("QTE_", p_star, "_true"))

truth = c(
  aipw_cpm_q_q1  = truth_q1,
  aipw_cpm_q_q0  = truth_q0,
  aipw_cpm_q_QTE = truth_qte,
  aipw_q1 = truth_q1,
  aipw_q0 = truth_q0,
  aipw_QTE = truth_qte
)

############################
# Sandwich
############################

summarize_sim_sandwich = function(out) {
  data.frame(
    target = c("q1", "q0", "QTE"),
    true = c(truth_q1, truth_q0, truth_qte),
    mean_est = c(
      mean(out$q1_hat, na.rm = TRUE),
      mean(out$q0_hat, na.rm = TRUE),
      mean(out$qte_hat, na.rm = TRUE)
    ),
    bias = c(
      mean(out$q1_hat - truth_q1, na.rm = TRUE),
      mean(out$q0_hat - truth_q0, na.rm = TRUE),
      mean(out$qte_hat - truth_qte, na.rm = TRUE)
    ),
    emp_sd = c(
      sd(out$q1_hat, na.rm = TRUE),
      sd(out$q0_hat, na.rm = TRUE),
      sd(out$qte_hat, na.rm = TRUE)
    ),
    mean_sand_se = c(
      mean(sqrt(out$var_sand_q1), na.rm = TRUE),
      mean(sqrt(out$var_sand_q0), na.rm = TRUE),
      mean(sqrt(out$var_sand_qte), na.rm = TRUE)
    ),
    coverage = c(
      100 * mean(out$ci_q1_l <= truth_q1 & out$ci_q1_u >= truth_q1, na.rm = TRUE),
      100 * mean(out$ci_q0_l <= truth_q0 & out$ci_q0_u >= truth_q0, na.rm = TRUE),
      100 * mean(out$ci_qte_l <= truth_qte & out$ci_qte_u >= truth_qte, na.rm = TRUE)
    ),
    emp_var = c(
      var(out$q1_hat, na.rm = TRUE),
      var(out$q0_hat, na.rm = TRUE),
      var(out$qte_hat, na.rm = TRUE)
    ),
    mean_sand_var = c(
      mean(out$var_sand_q1, na.rm = TRUE),
      mean(out$var_sand_q0, na.rm = TRUE),
      mean(out$var_sand_qte, na.rm = TRUE)
    )
  )
}

if (file.exists(sandwich_result_file)) {
  load(sandwich_result_file)
  summ_sand = summarize_sim_sandwich(res_all)
  print(summ_sand)
}

############################
# IF Based
############################

if (file.exists(if_result_file)) {
  load(if_result_file)
  
  summary_if = data.frame(
    parameter = c("q1", "q0", "qte"),
    mean_est = c(
      mean(res$q1_hat, na.rm = TRUE),
      mean(res$q0_hat, na.rm = TRUE),
      mean(res$qte_hat, na.rm = TRUE)
    ),
    bias = c(
      mean(res$q1_hat, na.rm = TRUE) - truth_q1,
      mean(res$q0_hat, na.rm = TRUE) - truth_q0,
      mean(res$qte_hat, na.rm = TRUE) - truth_qte
    ),
    emp_sd = c(
      sd(res$q1_hat, na.rm = TRUE),
      sd(res$q0_hat, na.rm = TRUE),
      sd(res$qte_hat, na.rm = TRUE)
    ),
    mean_se_if = c(
      mean(res$se_if_q1, na.rm = TRUE),
      mean(res$se_if_q0, na.rm = TRUE),
      mean(res$se_if_qte, na.rm = TRUE)
    ),
    coverage = c(
      mean(res$cover_q1, na.rm = TRUE),
      mean(res$cover_q0, na.rm = TRUE),
      mean(res$cover_qte, na.rm = TRUE)
    ),
    emp_var = c(
      var(res$q1_hat, na.rm = TRUE),
      var(res$q0_hat, na.rm = TRUE),
      var(res$qte_hat, na.rm = TRUE)
    ),
    mean_var_if = c(
      mean(res$var_if_q1, na.rm = TRUE),
      mean(res$var_if_q0, na.rm = TRUE),
      mean(res$var_if_qte, na.rm = TRUE)
    )
  )
  
  print(summary_if)
}

############################
# Bootstrap
############################

if (file.exists(bootstrap_result_file)) {
  boot_results_all = readRDS(bootstrap_result_file)
  
  cover95_each = lapply(boot_results_all, function(res) {
    ci95 = t(apply(res$boot_reps, 2, function(x) {
      quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
    }))
    
    as.numeric(ci95[, 1] <= truth[colnames(res$boot_reps)] &
                 truth[colnames(res$boot_reps)] <= ci95[, 2])
  })
  cover95_mat = do.call(rbind, cover95_each)
  
  boot_var_all = lapply(boot_results_all, function(res) {
    apply(res$boot_reps, 2, var, na.rm = TRUE)
  })
  boot_var_mat = do.call(rbind, boot_var_all)
  
  print(colMeans(cover95_mat, na.rm = TRUE))
  print(colMeans(boot_var_mat, na.rm = TRUE))
}

############################
# Empirical
############################

summarize_all_estimators = function(data, truth) {
  est_names = names(data)
  
  out = lapply(est_names, function(nm) {
    est = data[[nm]]
    true_val = as.numeric(truth[nm])
    
    mean_est = mean(est, na.rm = TRUE)
    emp_var = var(est, na.rm = TRUE)
    emp_se = sd(est, na.rm = TRUE)
    
    lower = est - 1.96 * emp_se
    upper = est + 1.96 * emp_se
    
    data.frame(
      estimator = nm,
      true = true_val,
      mean_est = mean_est,
      bias = mean_est - true_val,
      emp_var = emp_var,
      emp_se = emp_se,
      coverage = mean(lower <= true_val & upper >= true_val, na.rm = TRUE)
    )
  })
  
  do.call(rbind, out)
}

if (file.exists(empirical_result_file)) {
  data = read.csv(empirical_result_file)
  summ = summarize_all_estimators(data, truth)
  print(summ[c(3, 6), ])
}
