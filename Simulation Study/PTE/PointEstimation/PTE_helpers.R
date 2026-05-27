#####################################
# PTE point estimation
#####################################

or_cpm_pte = function(data, formY, y_threshold = 7.95) {
  
  # OR model
  g = orm(formY, family = "probit", x = TRUE, y = TRUE, data = data)
  alpha_hat = -g$coefficients[1:(length(g$yunique) - 1)]
  
  # F_Y1(y_threshold)
  g1_vals = cond_cdf_onemodel(
    data = data,
    y_star = y_threshold,
    g = g,
    alpha_hat = alpha_hat,
    A = 1,
    link = "probit"
  )
  FY1 = mean(g1_vals$g_y)
  
  # F_Y0(y_threshold)
  g0_vals = cond_cdf_onemodel(
    data = data,
    y_star = y_threshold,
    g = g,
    alpha_hat = alpha_hat,
    A = 0,
    link = "probit"
  )
  FY0 = mean(g0_vals$g_y)
  
  # PTE on CDF scale
  PTE = FY1 - FY0
  
  out = data.frame(
    threshold = y_threshold,
    FY0 = FY0,
    FY1 = FY1,
    PTE = PTE
  )
  
  return(out)
}



ipw_cpm_pte = function(data, y_threshold = 7.95, formA) {
  
  ps_fit = glm(formA, data = data, family = binomial("logit"))
  pi1 = predict(ps_fit, type = "response"); pi0 = 1 - pi1
  IA1 = as.numeric(data$A == 1); IA0 = as.numeric(data$A == 0)
  
  n = nrow(data); F1_ipw = F0_ipw = numeric(n)
  IY = as.numeric(data$Y <= y_threshold )
  FY1 = mean(IA1 * IY / pi1) 
  FY0 = mean(IA0 * IY / pi0)
  PTE = FY1 - FY0
  
  
  out = data.frame(
    threshold = y_threshold,
    FY0 = FY0,
    FY1 = FY1,
    PTE = PTE
  )
  
  return(out)
}


aipw_cpm_pte = function(data, formA, formY, y_threshold = 7.95) {
  
  # PS model
  ps_fit = glm(formA, data = data, family = binomial("logit"))
  pi1 = predict(ps_fit, type = "response");  pi0 = 1 - pi1
  IA1 = as.numeric(data$A == 1); IA0 = as.numeric(data$A == 0)
  
  # OR model: CPM
  g = orm(formY, family = "probit", x = TRUE, y = TRUE, data = data)
  alpha_hat = -g$coefficients[1:(length(g$yunique) - 1)]
  
  # predicted conditional CDF under A = 1
  g1_vals = cond_cdf_onemodel(
    data = data,
    y_star = y_threshold,
    g = g,
    alpha_hat = alpha_hat,
    A = 1,
    link = "probit"
  )
  F1_hat = g1_vals$g_y
  
  # predicted conditional CDF under A = 0
  g0_vals = cond_cdf_onemodel(
    data = data,
    y_star = y_threshold,
    g = g,
    alpha_hat = alpha_hat,
    A = 0,
    link = "probit"
  )
  F0_hat = g0_vals$g_y
  
  # observed indicator
  IY = as.numeric(data$Y <= y_threshold)
  
  # AIPW estimators of marginal CDFs
  FY1 = mean(F1_hat + IA1 / pi1 * (IY - F1_hat))
  FY0 = mean(F0_hat + IA0 / pi0 * (IY - F0_hat))
  PTE = FY1 - FY0
  
  out = data.frame(
    threshold = y_threshold,
    FY0 = FY0,
    FY1 = FY1,
    PTE = PTE
  )
  
  return(out)
}
###########################
# PS part
#############################
score_logit_glm = function(fit) {
  
  mf = model.frame(fit)
  X = model.matrix(fit)
  y = model.response(mf)
  
  eta = drop(X %*% coef(fit))
  #mu   = plogis(eta) 
  mu = fitted(fit)
  S   = X * ( (y -  mu))     
  attr(S, "total_score") = colSums(S) 
  S
}

dscore_logit_glm =function(fit) {
  X = model.matrix(fit)
  mu = fitted(fit); w = mu * (1 - mu)
  Hobs =-crossprod(X, X * w)     
  Hobs
}
#################################
# CPM part
##################################
expand_a_full = function(mod) {
  a = mod$info.matrix$a
  m = nrow(a)                 # number of cutpoints
  A = matrix(0, m, m)
  diag(A) = a[, 1]            # main diagonal
  if (m > 1) {
    off = a[1:(m - 1), 2]     # off-diagonal
    A[cbind(1:(m - 1), 2:m)] = off
    A[cbind(2:m, 1:(m - 1))] = off
  }
  rownames(A) = mod$info.matrix$iname
  colnames(A) = mod$info.matrix$iname
  A
}

build_info_matrix = function(mod) {
  A = expand_a_full(mod)                 # I_{alpha,alpha}
  B = mod$info.matrix$b                  # I_{beta,beta}
  AB = - mod$info.matrix$ab                # I_{alpha,beta}, alpha_orm = - \alpha
  B = 0.5 * (B + t(B))
  rownames(B) = colnames(B) = mod$info.matrix$xname
  rownames(AB) = mod$info.matrix$iname
  colnames(AB) = mod$info.matrix$xname
  I_full = rbind( cbind(A,AB), cbind(t(AB), B) )
  I_full
}


##################################
# F1, F0 part
###################################

get_EEs_OR_PTE = function(data, formY, y_threshold = 7.95, F1hat, F0hat, 
                          link = c("logistic", "probit")) {
  
  link = match.arg(link)
  
  g = rms::orm(
    formY,
    family = link,
    data = data,
    mscore = TRUE,
    x = TRUE,
    y = TRUE
  )
  
  alpha_hat = -g$coefficients[1:(length(g$yunique) - 1)]
  
  # Conditional CDF under A = 1 at fixed threshold
  g1_vals = cond_cdf_onemodel(
    data = data,
    y_star = y_threshold,
    g = g,
    alpha_hat = alpha_hat,
    A = 1,
    link = link
  )
  
  # Conditional CDF under A = 0 at fixed threshold
  g0_vals = cond_cdf_onemodel(
    data = data,
    y_star = y_threshold,
    g = g,
    alpha_hat = alpha_hat,
    A = 0,
    link = link
  )
  
  cF1_y = g1_vals$g_y
  cF0_y = g0_vals$g_y
  
  
  # Estimating equations for marginal CDFs
  EE1 = cF1_y - F1hat; EE0 = cF0_y - F0hat
  
  # Estimating equation for PTE
  EE_PTE = EE1 - EE0
  
  list(
    EE1 = EE1,
    EE0 = EE0,
    EE_PTE = EE_PTE
  )
  
}



get_EEs_IPW_PTE = function(data, formA, y_threshold = 7.95, F1hat, F0hat) {
  
  ps_fit = glm(formA, data = data, family = binomial("logit"))
  pi1 = predict(ps_fit, type = "response"); pi0 = 1 - pi1
  IA1 = as.numeric(data$A == 1); IA0 = as.numeric(data$A == 0)
  
  IY = as.numeric(data$Y <= y_threshold )
  cF1_y = IA1 * IY / pi1 
  cF0_y = IA0 * IY / pi0
  
  # Estimating equations for marginal CDFs
  EE1 = cF1_y - F1hat; EE0 = cF0_y - F0hat
  
  # Estimating equation for PTE
  EE_PTE = EE1 - EE0
  
  list(
    EE1 = EE1,
    EE0 = EE0,
    EE_PTE = EE_PTE
  )

}

get_EEs_AIPW_PTE = function(data, formA, formY, y_threshold = 7.95, 
                            F1hat, F0hat,
                            link = c("logistic", "probit")) {
  
  link = match.arg(link)
  
  # PS model
  ps_fit = glm(formA, data = data, family = binomial("logit"))
  pi1 = predict(ps_fit, type = "response")
  pi0 = 1 - pi1
  
  IA1 = as.numeric(data$A == 1)
  IA0 = as.numeric(data$A == 0)
  IY  = as.numeric(data$Y <= y_threshold)
  
  # OR model: CPM
  g = rms::orm(
    formY,
    family = link,
    data = data,
    mscore = TRUE,
    x = TRUE,
    y = TRUE
  )
  
  alpha_hat = -g$coefficients[1:(length(g$yunique) - 1)]
  
  # Conditional CDF under A = 1
  g1_vals = cond_cdf_onemodel(
    data = data,
    y_star = y_threshold,
    g = g,
    alpha_hat = alpha_hat,
    A = 1,
    link = link
  )
  
  # Conditional CDF under A = 0
  g0_vals = cond_cdf_onemodel(
    data = data,
    y_star = y_threshold,
    g = g,
    alpha_hat = alpha_hat,
    A = 0,
    link = link
  )
  
  cF1_y = g1_vals$g_y
  cF0_y = g0_vals$g_y
  
  # AIPW estimating functions for marginal CDFs
  AIPW1 = cF1_y + IA1 / pi1 * (IY - cF1_y)
  AIPW0 = cF0_y + IA0 / pi0 * (IY - cF0_y)
  
  EE1 = AIPW1 - F1hat
  EE0 = AIPW0 - F0hat
  
  # Estimating equation for PTE = F1 - F0
  EE_PTE = EE1 - EE0
  
  list(
    EE1 = EE1,
    EE0 = EE0,
    EE_PTE = EE_PTE
  )
}
#######################################
# A1xi
#######################################
Aaxi_cpm_OR_PTE = function(data, formY, y_threshold ,
                       treat_var = "A",
                       link = c("probit", "logistic")) {
  link = match.arg(link)
  cpm_fit = rms::orm( formY, family = link, data = data, mscore = TRUE, x = TRUE, y = TRUE )
  #ps_fit = glm(formA, data = data, family = binomial("logit"))
  
  n = nrow(data)
  y_levels = cpm_fit$yunique; n_alpha = length(y_levels) - 1 
  coef_all = cpm_fit$coefficients
  alpha_hat = -coef_all[1:n_alpha]; beta_hat  = coef_all[-(1:n_alpha)]
  
  X_design = cpm_fit$x
  X1_design = X_design; X0_design = X_design
  X1_design[, treat_var] = 1; X0_design[, treat_var] = 0
  
  # pi1_hat = predict(ps_fit, type = "response"); pi0_hat = 1 - pi1_hat
  # A1 = as.numeric(data[[treat_var]] == 1)
  # A0 = as.numeric(data[[treat_var]] == 0)
  
  compute_one_Axi = function(y_threshold, X_a) {
    # outside support: derivative = 0
    if (y_threshold < y_levels[1] || y_threshold >= y_levels[length(y_levels)]) {
      out = rep(0, length(coef_all))
      names(out) = names(coef_all)
      return(out)
    }
    j = max(which(y_levels <= y_threshold))
    if (j > n_alpha) {
      out = rep(0, length(coef_all))
      names(out) = names(coef_all)
      return(out)
    }
    
    u_hat = alpha_hat[j] - drop(X_a %*% beta_hat)
    f_e = switch( link, probit = dnorm(u_hat), logistic = dlogis(u_hat) )
    
    grad_alpha = matrix(0, nrow = n, ncol = n_alpha)
    grad_alpha[, j] = f_e; grad_beta = -X_a * f_e
    grad_xi = cbind(grad_alpha, grad_beta)
    
    Axi_hat = colMeans( grad_xi)
    names(Axi_hat) = c(names(coef_all)[1:n_alpha], names(beta_hat))
    
    return(Axi_hat)
  }
  
  #w1 = A1 / pi1_hat - 1; w0 = A0 / pi0_hat - 1
  A1xi_hat = -compute_one_Axi(y_threshold, X_a = X1_design)
  A0xi_hat = -compute_one_Axi(y_threshold, X_a = X0_design)
  
  return(list( A1xi = A1xi_hat, A0xi = A0xi_hat ))
}



Aaxi_cpm_AIPW_PTE = function(data, formA, formY, y_threshold,
                             treat_var = "A",
                             link = c("probit", "logistic")) {
  
  link = match.arg(link)
  
  # CPM model
  cpm_fit = rms::orm(
    formY,
    family = link,
    data = data,
    mscore = TRUE,
    x = TRUE,
    y = TRUE
  )
  
  # PS model
  ps_fit = glm(formA, data = data, family = binomial("logit"))
  pi1_hat = as.numeric(predict(ps_fit, type = "response"))
  pi0_hat = 1 - pi1_hat
  
  A1 = as.numeric(data[[treat_var]] == 1)
  A0 = as.numeric(data[[treat_var]] == 0)
  
  n = nrow(data)
  y_levels = cpm_fit$yunique
  n_alpha = length(y_levels) - 1
  
  coef_all = cpm_fit$coefficients
  alpha_hat = -coef_all[1:n_alpha]
  beta_hat  = coef_all[-(1:n_alpha)]
  
  X_design = cpm_fit$x
  
  X1_design = X_design
  X0_design = X_design
  
  X1_design[, treat_var] = 1
  X0_design[, treat_var] = 0
  
  compute_one_Axi = function(y_threshold, X_a, weight_a) {
    
    # outside support: derivative = 0
    if (y_threshold < y_levels[1] || y_threshold >= y_levels[length(y_levels)]) {
      out = rep(0, length(coef_all))
      names(out) = names(coef_all)
      return(out)
    }
    
    j = max(which(y_levels <= y_threshold))
    
    if (j > n_alpha) {
      out = rep(0, length(coef_all))
      names(out) = names(coef_all)
      return(out)
    }
    
    u_hat = alpha_hat[j] - drop(X_a %*% beta_hat)
    
    f_e = switch(
      link,
      probit   = dnorm(u_hat),
      logistic = dlogis(u_hat)
    )
    
    # derivative of F_a(y | X) wrt xi = (alpha, beta)
    grad_alpha = matrix(0, nrow = n, ncol = n_alpha)
    grad_alpha[, j] = f_e
    
    grad_beta = -X_a * f_e
    
    grad_xi = cbind(grad_alpha, grad_beta)
    
    # AIPW A block:
    # A_{a xi} = mean{ (I(A=a)/pi_a - 1) * dF_a/dxi }
    Axi_hat = colMeans(grad_xi * weight_a)
    
    names(Axi_hat) = c(names(coef_all)[1:n_alpha], names(beta_hat))
    
    return(Axi_hat)
  }
  
  # AIPW weights for derivative wrt CPM parameters
  w1 = A1 / pi1_hat - 1
  w0 = A0 / pi0_hat - 1
  
  A1xi_hat = compute_one_Axi(
    y_threshold = y_threshold,
    X_a = X1_design,
    weight_a = w1
  )
  
  A0xi_hat = compute_one_Axi(
    y_threshold = y_threshold,
    X_a = X0_design,
    weight_a = w0
  )
  
  return(list(
    A1xi = A1xi_hat,
    A0xi = A0xi_hat
  ))
}

#################################
# A1psi
####################################
Aapsi_ps_IPW_PTE = function(data, formA, y_threshold)  {
  
  n = nrow(data)
  ps_fit = glm(formA, data = data, family = binomial("logit"))
  X_tilde = model.matrix(ps_fit)
  pi1 = as.numeric(predict(ps_fit, type = "response")); pi0 = 1 - pi1
  IA1 = as.numeric(data$A == 1);IA0 = as.numeric(data$A == 0)
  IY_q1 = as.numeric(data$Y <= y_threshold);  IY_q0 = as.numeric(data$Y <= y_threshold)
  
  w1 = IA1 * ((1 - pi1) / pi1) * IY_q1 
  w0 = -IA0 * ((1 - pi0) / pi0)  * IY_q0 
  
  A1psi = colMeans(X_tilde * w1)
  A0psi = colMeans(X_tilde * w0)
  
  return(list(A1psi = A1psi,A0psi = A0psi))
}



Aapsi_ps_AIPW_PTE = function(data, formA, formY, y_threshold,
                              link = c("logistic", "probit")) {
  
  link = match.arg(link)
  
  # PS model
  ps_fit = glm(formA, data = data, family = binomial("logit"))
  X_tilde = model.matrix(ps_fit)
  
  pi1 = as.numeric(predict(ps_fit, type = "response"))
  pi0 = 1 - pi1
  
  IA1 = as.numeric(data$A == 1)
  IA0 = as.numeric(data$A == 0)
  IY  = as.numeric(data$Y <= y_threshold)
  
  # OR model: CPM
  g = rms::orm(
    formY,
    family = link,
    data = data,
    mscore = TRUE,
    x = TRUE,
    y = TRUE
  )
  
  alpha_hat = -g$coefficients[1:(length(g$yunique) - 1)]
  
  # Conditional CDF under A = 1
  g1_vals = cond_cdf_onemodel(
    data = data,
    y_star = y_threshold,
    g = g,
    alpha_hat = alpha_hat,
    A = 1,
    link = link
  )
  cF1_y = g1_vals$g_y
  
  # Conditional CDF under A = 0
  g0_vals = cond_cdf_onemodel(
    data = data,
    y_star = y_threshold,
    g = g,
    alpha_hat = alpha_hat,
    A = 0,
    link = link
  )
  cF0_y = g0_vals$g_y
  
  # Derivative wrt PS parameters
  # For A = 1:
  # d/dpsi [ IA1 / pi1 * (IY - cF1_y) ]
  # = - IA1 * (1 - pi1) / pi1 * (IY - cF1_y) * X
  #
  # A matrix uses - dU/dtheta, so sign becomes positive.
  w1 = IA1 * ((1 - pi1) / pi1) * (IY - cF1_y)
  
  # For A = 0:
  # d/dpsi [ IA0 / pi0 * (IY - cF0_y) ]
  # = IA0 * pi1 / pi0 * (IY - cF0_y) * X
  #
  # A matrix uses - dU/dtheta, so sign becomes negative.
  w0 = -IA0 * (pi1 / pi0) * (IY - cF0_y)
  
  A1psi = colMeans(X_tilde * w1)
  A0psi = colMeans(X_tilde * w0)
  
  return(list(
    A1psi = A1psi,
    A0psi = A0psi
  ))
}
