
library(rms)
set.seed(37203)
BigN = 1e8
alpha1 = 0.5; alpha2 = 0.35  # coefficients for treatment mechanism
beta1 = -2; beta2 = 3       # coefficients for Y generating
delta = 2                  # treatment effect

# Function to generate valid p-values for both Y0_star and Y1_star
get_valid_pvals = function(n_needed, beta1, beta2, delta) {
  # Generate variabels
  X1_samp = rbinom(n_needed, size = 1, prob = 0.5)
  X2_samp = rnorm(n_needed, mean = 0, sd = 1)
  
  # Y0_star = beta1 * X1_samp + beta2 * X2_samp + rlogis(n_needed, location = 0,scale = 3/(pi^2))
  # Y1_star = beta1 * X1_samp + beta2 * X2_samp + delta + rlogis(n_needed,  location = 0,scale = 3/(pi^2))
  
  Y0_star = beta1 * X1_samp + beta2 * X2_samp + rnorm(n_needed, 0,1)
  Y1_star = beta1 * X1_samp + beta2 * X2_samp + delta + rnorm(n_needed, 0,1)
  
  p0 = pnorm(Y0_star); p1 = pnorm(Y1_star)
  
  # Filter to keep only samples where both p-values are strictly within (1e-15, 1-1e-15)
  valid_idx = which((p0 > 1e-15 & p0 < (1 - 1e-15)) & (p1 > 1e-15 & p1 < (1 - 1e-15)))
  
  list(p0 = p0[valid_idx], p1 = p1[valid_idx])
}

###############################################
# Generate valid p-values for both Y0 and Y1 simultaneously
###############################################
valid_pvals = get_valid_pvals(BigN, beta1, beta2, delta)
while (length(valid_pvals$p0) < BigN) {
  n_missing = BigN - length(valid_pvals$p0)
  additional = get_valid_pvals(n_missing, beta1, beta2, delta)
  valid_pvals$p0 = c(valid_pvals$p0, additional$p0)
  valid_pvals$p1 = c(valid_pvals$p1, additional$p1)
}

# Keep only the first BigN valid p-values for both outcomes
p0_final = valid_pvals$p0[1:BigN]
p1_final = valid_pvals$p1[1:BigN]

# Transform the valid p-values to chi-square distributed outcomes (df = 5)
Y0_samp = qchisq(p0_final, df = 5)
Y1_samp = qchisq(p1_final, df = 5)

###########################################
# Results
#########################################

# QTE 0.25
m_Y1_0.25 = as.numeric(quantile(Y1_samp, probs= 0.25, names=FALSE))
m_Y0_0.25 = as.numeric(quantile(Y0_samp, probs= 0.25, names=FALSE))
QTE_0.25_true = m_Y1_0.25 - m_Y0_0.25

cat("q1 when p = 0.25:", m_Y1_0.25, "\n")
cat("q0 when p = 0.25:", m_Y0_0.25, "\n")
cat("QTE when p = 0.25:", QTE_0.25_true, "\n")

log_m_Y1_0.25 = log(m_Y1_0.25)
log_m_Y0_0.25 = log(m_Y0_0.25)
log_QTE_0.25_true = log(QTE_0.25_true)

# QTE 0.5
m_Y1_0.5 = as.numeric(quantile(Y1_samp, probs= 0.5, names=FALSE))
m_Y0_0.5 = as.numeric(quantile(Y0_samp, probs= 0.5, names=FALSE))
QTE_0.5_true = m_Y1_0.5 - m_Y0_0.5

cat("q1 when p = 0.5:", m_Y1_0.5, "\n")
cat("q0 when p = 0.5:", m_Y0_0.5, "\n")
cat("QTE when p = 0.5:", QTE_0.5_true, "\n")

log_m_Y1_0.5 = log(m_Y1_0.5)
log_m_Y0_0.5 = log(m_Y0_0.5)
log_QTE_0.5_true = log(QTE_0.5_true)

# QTE 0.75
m_Y1_0.75 = as.numeric(quantile(Y1_samp, probs= 0.75, names=FALSE))
m_Y0_0.75 = as.numeric(quantile(Y0_samp, probs= 0.75, names=FALSE))
QTE_0.75_true = m_Y1_0.75 - m_Y0_0.75

cat("q1 when p = 0.75:", m_Y1_0.75, "\n")
cat("q0 when p = 0.75:", m_Y0_0.75, "\n")
cat("QTE when p = 0.75:", QTE_0.75_true, "\n")

log_m_Y1_0.75 = log(m_Y1_0.75)
log_m_Y0_0.75 = log(m_Y0_0.75)
log_QTE_0.75_true = log(QTE_0.75_true)

# QTE 0.1
m_Y1_0.1 = as.numeric(quantile(Y1_samp, probs= 0.1, names=FALSE))
m_Y0_0.1 = as.numeric(quantile(Y0_samp, probs= 0.1, names=FALSE))
QTE_0.1_true = m_Y1_0.1 - m_Y0_0.1

cat("q1 when p = 0.1:", m_Y1_0.1, "\n")
cat("q0 when p = 0.1:", m_Y0_0.1, "\n")
cat("QTE when p = 0.1:", QTE_0.1_true, "\n")

log_m_Y1_0.1 = log(m_Y1_0.1)
log_m_Y0_0.1 = log(m_Y0_0.1)
log_QTE_0.1_true = log(QTE_0.1_true)

# QTE 0.9
m_Y1_0.9 = as.numeric(quantile(Y1_samp, probs= 0.9, names=FALSE))
m_Y0_0.9 = as.numeric(quantile(Y0_samp, probs= 0.9, names=FALSE))
QTE_0.9_true = m_Y1_0.9 - m_Y0_0.9

cat("q1 when p = 0.9:", m_Y1_0.9, "\n")
cat("q0 when p = 0.9:", m_Y0_0.9, "\n")
cat("QTE when p = 0.9:", QTE_0.9_true, "\n")

log_m_Y1_0.9 = log(m_Y1_0.9)
log_m_Y0_0.9 = log(m_Y0_0.9)
log_QTE_0.9_true = log(QTE_0.9_true)



# QTE 0.99
m_Y1_0.99 = as.numeric(quantile(Y1_samp, probs= 0.99, names=FALSE))
m_Y0_0.99 = as.numeric(quantile(Y0_samp, probs= 0.99, names=FALSE))
QTE_0.99_true = m_Y1_0.99 - m_Y0_0.99

cat("q1 when p = 0.99:", m_Y1_0.99, "\n")
cat("q0 when p = 0.99:", m_Y0_0.99, "\n")
cat("QTE when p = 0.99:", QTE_0.99_true, "\n")

log_m_Y1_0.99 = log(m_Y1_0.99)
log_m_Y0_0.99 = log(m_Y0_0.99)
log_QTE_0.99_true = log(QTE_0.99_true)

# QTE 0.05
m_Y1_0.05 = as.numeric(quantile(Y1_samp, probs= 0.05, names=FALSE))
m_Y0_0.05 = as.numeric(quantile(Y0_samp, probs= 0.05, names=FALSE))
QTE_0.05_true = m_Y1_0.05 - m_Y0_0.05

cat("q1 when p = 0.05:", m_Y1_0.05, "\n")
cat("q0 when p = 0.05:", m_Y0_0.05, "\n")
cat("QTE when p = 0.05:", QTE_0.05_true, "\n")

log_m_Y1_0.05 = log(m_Y1_0.05)
log_m_Y0_0.05 = log(m_Y0_0.05)
log_QTE_0.05_true = log(QTE_0.05_true)

summary(Y0_samp); summary(Y1_samp)
# 
png("Large_sample_Y0_Y1-original scale.png", width = 800, height = 400)
par(mfrow = c(1, 2))
hist(Y0_samp, main = "Large sample: true Y0-original scale"); hist(Y1_samp, main = "Large sample: true Y1-original scale")
dev.off()

png("Large_sample_Y0_Y1-log scale.png", width = 800, height = 400)
par(mfrow = c(1, 2))
hist(log(Y0_samp), main = "Large sample: true Y0-log scale"); hist(log(Y1_samp), main = "Large sample: true Y1-log scale")
dev.off()


remove(Y1_samp); remove(Y0_samp)

save(m_Y0_0.25, m_Y1_0.25, QTE_0.25_true, 
     m_Y0_0.5, m_Y1_0.5, QTE_0.5_true, 
     m_Y0_0.75, m_Y1_0.75, QTE_0.75_true, 
     m_Y0_0.1, m_Y1_0.1, QTE_0.1_true, 
     m_Y0_0.9, m_Y1_0.9, QTE_0.9_true, 
     m_Y0_0.99, m_Y1_0.99, QTE_0.99_true,
     m_Y0_0.05, m_Y1_0.05, QTE_0.05_true,
     
     log_m_Y0_0.25, log_m_Y1_0.25, log_QTE_0.25_true, 
     log_m_Y0_0.5, log_m_Y1_0.5, log_QTE_0.5_true, 
     log_m_Y0_0.75, log_m_Y1_0.75, log_QTE_0.75_true,
     log_m_Y0_0.1, log_m_Y1_0.1, log_QTE_0.1_true, 
     log_m_Y0_0.9, log_m_Y1_0.9, log_QTE_0.9_true, 
     log_m_Y0_0.99, log_m_Y1_0.99, log_QTE_0.99_true, 
     log_m_Y0_0.05, log_m_Y1_0.05, log_QTE_0.05_true, 
     
     file = "true_values.RData")
