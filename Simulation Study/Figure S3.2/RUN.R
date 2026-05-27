library(rms)

set.seed(37023)

# Set parameters
nsample = 1000

source("single_sim.R")
source("Helper_functions.R")

# Generate test data for setting up datadist
test_data = generate_data(n = nsample, alpha1 = 0.5, alpha2 = 0.35,
                           beta1 = -2, beta2 = 3, delta = 2)
dd = datadist(test_data)
options(datadist = "dd")

# Run multiple simulations and collect timing results
n_rep = 1000

res_list = vector("list", length = n_rep)

for (i in seq_len(n_rep)) {
  res_list[[i]] = single_sim(
    n = 1000,
    p_star = 0.5,
    alpha1 = 0.5,
    alpha2 = 0.35,
    beta1 = -2,
    beta2 = 3,
    delta = 2,
    scenario = "cc"
  )

  if (i %% 1 == 0) {
    cat("Completed", i, "simulations\n")
  }
}

res_mat = do.call(rbind, res_list)
time_cols = grep("^time_", colnames(res_mat), value = TRUE)
time_mat = res_mat[, time_cols]

time_summary = data.frame(
  method = sub("^time_", "", time_cols),
  mean = colMeans(time_mat),
  sd = apply(time_mat, 2, sd)
)

write.csv(
  time_summary,
  file = "time_summary.csv",
  row.names = FALSE
)
