
set.seed(37203)

BigN = 1e8
alpha1 = 0.5; alpha2 = 0.35
beta1  = -2;  beta2  = 3
delta  = 2

p_vec = c(0.1, 0.25, 0.5, 0.75, 0.9)

get_valid_pvals = function(n_needed, beta1, beta2, delta) {
  X1_samp = rbinom(n_needed, size = 1, prob = 0.5)
  X2_samp = rnorm(n_needed, mean = 0, sd = 1)
  
  Y0_star = beta1 * X1_samp + beta2 * X2_samp + rnorm(n_needed, 0, 1)
  Y1_star = beta1 * X1_samp + beta2 * X2_samp + delta + rnorm(n_needed, 0, 1)
  
  p0 = pnorm(Y0_star)
  p1 = pnorm(Y1_star)
  
  valid_idx = which(
    (p0 > 1e-15 & p0 < 1 - 1e-15) &
      (p1 > 1e-15 & p1 < 1 - 1e-15)
  )
  
  list(p0 = p0[valid_idx], p1 = p1[valid_idx])
}

valid_pvals = get_valid_pvals(BigN, beta1, beta2, delta)

while (length(valid_pvals$p0) < BigN) {
  n_missing = BigN - length(valid_pvals$p0)
  additional = get_valid_pvals(n_missing, beta1, beta2, delta)
  valid_pvals$p0 = c(valid_pvals$p0, additional$p0)
  valid_pvals$p1 = c(valid_pvals$p1, additional$p1)
}

p0_final = valid_pvals$p0[1:BigN]
p1_final = valid_pvals$p1[1:BigN]

Y0_samp = qchisq(p0_final, df = 5)
Y1_samp = qchisq(p1_final, df = 5)

##################################################
# True marginal quantiles
##################################################

q0_true = as.numeric(quantile(Y0_samp, probs = p_vec, names = FALSE))
q1_true = as.numeric(quantile(Y1_samp, probs = p_vec, names = FALSE))

names(q0_true) = paste0("q0_", p_vec)
names(q1_true) = paste0("q1_", p_vec)

##################################################
# Marginal CDF and PTE at q0_p and q1_p
##################################################

true_values = data.frame()

for (j in seq_along(p_vec)) {
  
  p = p_vec[j]
  q0 = q0_true[j]
  q1 = q1_true[j]
  
  # Evaluate marginal CDFs at q0_p
  FY0_q0 = mean(Y0_samp <= q0)
  FY1_q0 = mean(Y1_samp <= q0)
  PTE_q0 = FY1_q0 - FY0_q0
  
  # Evaluate marginal CDFs at q1_p
  FY0_q1 = mean(Y0_samp <= q1)
  FY1_q1 = mean(Y1_samp <= q1)
  PTE_q1 = FY1_q1 - FY0_q1
  
  true_values = rbind(
    true_values,
    data.frame(
      p = p,
      threshold_type = paste0("q0_", p),
      threshold = q0,
      FY0 = FY0_q0,
      FY1 = FY1_q0,
      PTE = PTE_q0
    ),
    data.frame(
      p = p,
      threshold_type = paste0("q1_", p),
      threshold = q1,
      FY0 = FY0_q1,
      FY1 = FY1_q1,
      PTE = PTE_q1
    )
  )
}

print(true_values)

rownames(true_values) = true_values$threshold_type

##################################################
# Save true values
##################################################

save(
  q0_true,
  q1_true,
  true_values,
  file = "true_values.RData"
)
