


##############################################################################
# 1. File names and variables to keep
##############################################################################

set.seed(20260526)

input_file = "data_model_VL_YA.RData"
output_file = "data_model_VL_YA_synthetic.RData"

analysis_vars = c(
  "Y",
  "A",
  "age",
  "sex",
  "site",
  "route_infection",
  "prior_aids",
  "calendar_year",
  "months_toVLmeasure",
  "cd4_baseline_tran",
  "singledl_rna_baseline_v_tran"
)


##############################################################################
# 2. Load the private local analysis data
##############################################################################


load(input_file)

real_data = as.data.frame(data_model)
n = nrow(real_data)


##############################################################################
# 3. Small helper functions
##############################################################################

clip = function(x, lower, upper) {
  pmin(pmax(x, lower), upper)
}

num_bounds = function(x, probs = c(0.01, 0.99)) {
  stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
}

safe_sd = function(x) {
  sx = stats::sd(x, na.rm = TRUE)
  if (is.na(sx) || sx == 0) 1 else sx
}

sample_factor = function(x, n) {
  x = droplevels(as.factor(x))
  lev = levels(x)
  tab = table(x)
  probs = (as.numeric(tab) + 1) / sum(as.numeric(tab) + 1)
  factor(sample(lev, n, replace = TRUE, prob = probs), levels = lev)
}

sample_binary = function(x, n) {
  p = mean(x == 1, na.rm = TRUE)
  if (is.na(p)) p = 0.5
  stats::rbinom(n, 1, clip(p, 0.05, 0.95))
}


##############################################################################
# 4. Generate continuous covariates jointly
##############################################################################

continuous_vars = c(
  "age",
  "months_toVLmeasure",
  "cd4_baseline_tran",
  "singledl_rna_baseline_v_tran"
)

continuous_data = real_data[, continuous_vars, drop = FALSE]
continuous_data = as.data.frame(lapply(continuous_data, function(x) {
  x[is.na(x)] = mean(x, na.rm = TRUE)
  x
}))

mu = colMeans(continuous_data)
sigma = stats::cov(continuous_data)
sigma[is.na(sigma)] = 0
diag(sigma) = pmax(diag(sigma), 1e-6)

chol_sigma = tryCatch(
  chol(sigma),
  error = function(e) chol(diag(diag(sigma), nrow(sigma)))
)

synthetic_continuous = matrix(stats::rnorm(n * length(mu)), nrow = n) %*% chol_sigma
synthetic_continuous = sweep(synthetic_continuous, 2, mu, "+")
synthetic_continuous = as.data.frame(synthetic_continuous)
names(synthetic_continuous) = names(mu)

for (v in names(synthetic_continuous)) {
  bounds = num_bounds(real_data[[v]])
  synthetic_continuous[[v]] = clip(synthetic_continuous[[v]], bounds[1], bounds[2])
}


##############################################################################
# 5. Build the synthetic covariate data frame
##############################################################################

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
# 6. Generate treatment assignment
##############################################################################

treatment_formula = treatment ~ age + sex + site + route_infection + prior_aids +
  calendar_year + cd4_baseline_tran + singledl_rna_baseline_v_tran +
  months_toVLmeasure

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


##############################################################################
# 7. Generate viral load outcome Y
##############################################################################

vl_floor = min(real_data$Y, na.rm = TRUE)

# VL has a large point mass at the lower detection limit. Generate that mass
# separately, then generate positive values above the limit on the log scale.
floor_formula = I(Y == vl_floor) ~ A + age + sex + site + route_infection +
  prior_aids + calendar_year + cd4_baseline_tran +
  singledl_rna_baseline_v_tran + months_toVLmeasure

floor_model = stats::glm(floor_formula, data = real_data, family = stats::binomial())
p_floor = stats::predict(floor_model, newdata = data_model, type = "response")
p_floor = clip(p_floor, 0.05, 0.95)
is_floor = stats::rbinom(n, 1, p_floor) == 1

above_floor = real_data$Y > vl_floor
outcome_formula = log(Y) ~ A + age + sex + site + route_infection + prior_aids +
  calendar_year + cd4_baseline_tran + singledl_rna_baseline_v_tran +
  months_toVLmeasure

outcome_model = stats::lm(outcome_formula, data = real_data[above_floor, , drop = FALSE])
outcome_resid_sd = safe_sd(stats::residuals(outcome_model))

synthetic_log_y = stats::predict(outcome_model, newdata = data_model) +
  stats::rnorm(n, 0, outcome_resid_sd)

y_bounds = num_bounds(real_data$Y[above_floor])
data_model$Y = round(clip(exp(synthetic_log_y), y_bounds[1], y_bounds[2]))
data_model$Y[is_floor] = vl_floor


##############################################################################
# 8. Final checks and save
##############################################################################

missing_analysis_vars = setdiff(analysis_vars, names(data_model))
if (length(missing_analysis_vars) > 0) {
  stop("Missing required analysis variables: ", paste(missing_analysis_vars, collapse = ", "))
}

data_model = data_model[, analysis_vars, drop = FALSE]

data_model$site = droplevels(data_model$site)
data_model$route_infection = droplevels(data_model$route_infection)

save(data_model, file = output_file)

