###############################################################################
# Helper functions for the CD4 real-data analysis
###############################################################################


###############################################################################
# Quantile inversion
###############################################################################

quantile_inversion = function(Y_sorted, CDF_sorted, p_target) {
  ncat = length(Y_sorted)

  denom = CDF_sorted[ncat - 1] - CDF_sorted[1]
  if (abs(denom) < .Machine$double.eps) {
    return(Y_sorted[1])
  }

  weight_yuqi = c(0, (CDF_sorted[-ncat] - CDF_sorted[1]) / denom, 1)
  weight = 1 - weight_yuqi

  weighted_quantile = (1 - weight) * c(Y_sorted[1], Y_sorted) +
    weight * c(Y_sorted, Y_sorted[ncat])

  approx(x = c(0, CDF_sorted), y = weighted_quantile, xout = p_target)$y
}


###############################################################################
# Conditional CDF from two outcome models
###############################################################################

# Compute fitted P(Y <= y_star | X) from one arm-specific rms::orm model.
#
# Args:
#   data: data frame used for prediction; only nrow(data) is used here.
#   y_star: outcome threshold at which to evaluate the fitted CDF.
#   g: arm-specific orm fit, e.g. the model fitted among A = 1 or A = 0.
#   all_g: pooled orm fit whose design matrix contains the full analysis sample.
#   alpha_hat: threshold coefficients from g, usually
#     -g$coefficients[1:(length(g$yunique) - 1)].
#   tolerance: matching tolerance for y_star against g$yunique.
#
# Returns:
#   A list with g_y, the fitted CDF values for all rows in data.
calc_g_twomodels = function(data, y_star, g, all_g, alpha_hat, tolerance) {
  n = nrow(data)

  if (y_star < g$yunique[1]) {
    return(list(g_y = rep(0, n)))
  }

  if (y_star >= g$yunique[length(g$yunique)]) {
    return(list(g_y = rep(1, n)))
  }

  idx_star = max(which(g$yunique <= y_star + tolerance))
  idx_star = min(idx_star, length(alpha_hat))

  myalpha = alpha_hat[idx_star]

  slope_start = length(g$yunique)
  slope_coeff = g$coefficients[slope_start:length(g$coefficients)]

  design_matrix = all_g$x[, colnames(all_g$x) != "A", drop = FALSE]

  coeff_names = names(slope_coeff)
  if (!is.null(coeff_names) && all(coeff_names %in% colnames(design_matrix))) {
    design_matrix = design_matrix[, coeff_names, drop = FALSE]
  } else if (ncol(design_matrix) != length(slope_coeff)) {
    stop("Design matrix columns do not match arm-specific orm coefficients.")
  }

  intercept_pred = as.matrix(design_matrix) %*% slope_coeff
  g_y = plogis(myalpha - intercept_pred)

  list(g_y = as.numeric(g_y))
}
