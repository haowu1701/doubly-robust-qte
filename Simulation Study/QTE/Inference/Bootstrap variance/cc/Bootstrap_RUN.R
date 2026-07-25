library(Matrix)
library(MASS)
library(rms)
load("nsim1000_data_n1000.RData")
load("true_values.RData")
source("variance_helpers.R")
source("Helper_functions.R")
truth = c(
  aipw_cpm_q_q1  = m_Y1_0.5,
  aipw_cpm_q_q0  = m_Y0_0.5,
  aipw_cpm_q_QTE = QTE_0.5_true,
  aipw_q1 = m_Y1_0.5,
  aipw_q0 = m_Y0_0.5,
  aipw_QTE = QTE_0.5_true
)


test_data = generate_data( n = 1000,alpha1 = 0.5, alpha2 = 0.35, beta1 = -2, beta2 = 3, delta = 2 )
dd = datadist(test_data); options(datadist = "dd")

#### RUN
boot_results_all = lapply(seq_along(data_list), function(i) {
  cat("Running dataset", i, "\n")
  
  bootstrap_one_dataset(
    data = data_list[[i]],
    p_star = 0.5,
    scenario = "cc",
    B = 500,
    seed = 1000 + i,
    verbose  = TRUE
  )
  
  
})

saveRDS(boot_results_all, file = "boot_results_all_p0.5_cc.rds")


# SUMMARY
## variance
boot_var_all = lapply(boot_results_all, function(res) {
  apply(res$boot_reps, 2, var, na.rm = TRUE)
})
boot_var_mat = do.call(rbind, boot_var_all)



## coverage probability
cover95_each = lapply(boot_results_all, function(res) {
  ci95 = t(apply(res$boot_reps, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
  }))
  as.numeric(ci95[, 1] <= truth[colnames(res$boot_reps)] &
               truth[colnames(res$boot_reps)] <= ci95[, 2])
})
cover95_mat = do.call(rbind, cover95_each)



colMeans(boot_var_mat, na.rm = TRUE)
colMeans(cover95_mat, na.rm = TRUE)
