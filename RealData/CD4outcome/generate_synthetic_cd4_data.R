# Generate analysis-compatible synthetic CD4 data for the public GitHub repo.
#
# The private data file is used only to estimate coarse distributions and
# regression relationships. The output rows are newly simulated and are intended
# only to demonstrate that the analysis code runs.


##############################################################################
# 1. File names and variables to keep
##############################################################################

set.seed(20260526)

input_file = "data_model_cd4_YA.RData"
output_file = "data_model_cd4_YA_synthetic.RData"

# These are the variables needed by MQ_PTE_CD4.Rmd.
analysis_vars = c(
  "Y",
  "A",
  "age",
  "sex",
  "site",
  "route_infection",
  "prior_aids",
  "calendar_year",
  "months_toCD4measure",
  "cd4_baseline_tran",
  "singledl_rna_baseline_v_tran"
)


###############################################################################
# 2. Load the private local analysis data
###############################################################################

load(input_file)
real_data = as.data.frame(data_model)
n = nrow(real_data)


###############################################################################
# 3. Small helper functions
###############################################################################

# Restrict generated values to a plausible range.
clip = function(x, lower, upper) {
  pmin(pmax(x, lower), upper)
}

# Use the 1st and 99th percentiles as conservative bounds.
num_bounds = function(x, probs = c(0.01, 0.99)) {
  stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
}

# Avoid zero or missing standard deviations when generating random values.
safe_sd = function(x) {
  sx = stats::sd(x, na.rm = TRUE)
  if (is.na(sx) || sx == 0) 1 else sx
}

# Sample a categorical variable using its empirical distribution.
sample_factor = function(x, n) {
  x = droplevels(as.factor(x))
  lev = levels(x)
  tab = table(x)
  probs = (as.numeric(tab) + 1) / sum(as.numeric(tab) + 1)
  factor(sample(lev, n, replace = TRUE, prob = probs), levels = lev)
}

# Sample a binary variable using its empirical event probability.
sample_binary = function(x, n) {
  p = mean(x == 1, na.rm = TRUE)
  if (is.na(p)) p = 0.5
  stats::rbinom(n, 1, clip(p, 0.05, 0.95))
}


###############################################################################
# 4. Generate continuous covariates jointly
###############################################################################

# These variables are treated as continuous.
continuous_vars = c(
  "age",
  "months_toCD4measure",
  "cd4_baseline_tran",
  "singledl_rna_baseline_v_tran"
)

continuous_data = real_data[, continuous_vars, drop = FALSE]

# Fill missing values before estimating the covariance matrix.
continuous_data = as.data.frame(lapply(continuous_data, function(x) {
  x[is.na(x)] = mean(x, na.rm = TRUE)
  x
}))

mu = colMeans(continuous_data)
sigma = stats::cov(continuous_data)
sigma[is.na(sigma)] = 0
diag(sigma) = pmax(diag(sigma), 1e-6)

# Cholesky decomposition lets us generate correlated normal variables with the
# same approximate mean/covariance structure as the private data.
chol_sigma = tryCatch(
  chol(sigma),
  error = function(e) chol(diag(diag(sigma), nrow(sigma)))
)

synthetic_continuous = matrix(stats::rnorm(n * length(mu)), nrow = n) %*% chol_sigma
synthetic_continuous = sweep(synthetic_continuous, 2, -mu, "-")
synthetic_continuous = as.data.frame(synthetic_continuous)
names(synthetic_continuous) = names(mu)

# Keep each generated continuous variable inside the private data's central
# range to avoid unrealistic public demonstration values.
for (v in names(synthetic_continuous)) {
  bounds = num_bounds(real_data[[v]])
  synthetic_continuous[[v]] = clip(synthetic_continuous[[v]], bounds[1], bounds[2])
}



# Start from one copied row only to preserve column classes and factor levels.
# All analysis variables are overwritten below with synthetic values.
data_model = real_data[rep(1, n), , drop = FALSE]

for (v in names(synthetic_continuous)) {
  data_model[[v]] = synthetic_continuous[[v]]
}

data_model$sex = sample_binary(real_data$sex, n)
data_model$prior_aids = sample_binary(real_data$prior_aids, n)
data_model$site = sample_factor(real_data$site, n)
data_model$route_infection = sample_factor(real_data$route_infection, n)
data_model$calendar_year = as.numeric(as.character(sample_factor(real_data$calendar_year, n)))


##############################################################################
# 5. Generate treatment assignment
##############################################################################

# Fit a simple propensity model on the private data, then use the fitted model
# only to generate treatment probabilities for the synthetic covariates.
treatment_formula = treatment ~ age + sex + site + route_infection + prior_aids +
  calendar_year + cd4_baseline_tran + singledl_rna_baseline_v_tran

treatment_model = tryCatch(
  stats::glm(treatment_formula, data = real_data, family = stats::binomial()),
  error = function(e) NULL
)

if (!is.null(treatment_model)) {
  p_a = stats::predict(treatment_model, newdata = data_model, type = "response")
  p_a = clip(p_a, 0.05, 0.95)
  data_model$A = stats::rbinom(n, 1, p_a)
} else {
  data_model$A = sample_binary(real_data$treatment, n)
}


###############################################################################
# 6. Generate CD4 outcome Y
##############################################################################

# Fit a simple outcome model on log(Y), then add residual noise so the
# synthetic outcome is not a deterministic function of the covariates.
outcome_formula = log(Y) ~ A + age + sex + site + route_infection + prior_aids +
  calendar_year + cd4_baseline_tran + singledl_rna_baseline_v_tran + months_toCD4measure

outcome_model = stats::lm(outcome_formula, data = real_data)
outcome_resid_sd = safe_sd(stats::residuals(outcome_model))

synthetic_log_y = stats::predict(outcome_model, newdata = data_model) +
  stats::rnorm(n, 0, outcome_resid_sd)

y_bounds = num_bounds(real_data$Y)
data_model$Y = round(clip(exp(synthetic_log_y), y_bounds[1], y_bounds[2]))


###############################################################################
# 7. Final checks and save
###############################################################################

missing_analysis_vars = setdiff(analysis_vars, names(data_model))
data_model = data_model[, analysis_vars, drop = FALSE]
data_model$site = droplevels(data_model$site)
data_model$route_infection = droplevels(data_model$route_infection)

save(data_model, file = output_file)
