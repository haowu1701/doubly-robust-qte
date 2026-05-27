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

########################################
# Emprical inference
########################################

one_run_empirical_from_data = function(data,
                                       y_threshold = 7.95,
                                       scenario = "cc") {
  options(warn = -1)
  
  # scenarios { "cc", "cm", "mc", "mm" }
  if (scenario == "cc") {
    formA  = A ~ X1 + X2
    formY  = Y ~ X1 + X2 + A
    
  } else if (scenario == "cm") {
    formA  = A ~ X1
    formY  = Y ~ X1 + X2 + A
    
  } else if (scenario == "mc") {
    formA  = A ~ X1 + X2
    formY  = Y ~ A + X1
    
  } else if (scenario == "mm") {
    formA  = A ~ X1
    formY  = Y ~ A + X1
    
  } else {
    stop("scenario must be one of: cc, cm, mc, mm")
  }
  
  AIPW = aipw_cpm_pte(data, y_threshold = y_threshold,
                      formA = formA, formY = formY)
  
  OR = or_cpm_pte(data, y_threshold = y_threshold,
                  formY = formY)
  
  IPW = ipw_cpm_pte(data, y_threshold = y_threshold,
                    formA = formA)
  
  result = c(
    AIPW_Fy1_hat = as.numeric(unname(AIPW["FY1"])),
    AIPW_Fy0_hat = as.numeric(unname(AIPW["FY0"])),
    AIPW_pte_hat = as.numeric(unname(AIPW["FY1"])) -
      as.numeric(unname(AIPW["FY0"])),
    
    OR_Fy1_hat = as.numeric(unname(OR["FY1"])),
    OR_Fy0_hat = as.numeric(unname(OR["FY0"])),
    OR_pte_hat = as.numeric(unname(OR["FY1"])) -
      as.numeric(unname(OR["FY0"])),
    
    IPW_Fy1_hat = as.numeric(unname(IPW["FY1"])),
    IPW_Fy0_hat = as.numeric(unname(IPW["FY0"])),
    IPW_pte_hat = as.numeric(unname(IPW["FY1"])) -
      as.numeric(unname(IPW["FY0"]))
  )
  
  return(result)
}


run_empirical_on_datalist_PTE = function(data_list,
                             y_threshold = 7.95,
                             scenario = "cc",
                             n_run = length(data_list)) {
  
  all_results = vector("list", n_run)
  
  for (i in 1:n_run) {
    
    message("Running empirical PTE ", i, " / ", n_run)
    
    all_results[[i]] = tryCatch({
      
      est_i = one_run_empirical_from_data(
        data = data_list[[i]],
        y_threshold = y_threshold,
        scenario = scenario
      )
      
      data.frame(
        sim = i,
        t(est_i),
        success = TRUE,
        error = NA
      )
      
    }, error = function(e) {
      
      data.frame(
        sim = i,
        AIPW_Fy1_hat = NA,
        AIPW_Fy0_hat = NA,
        AIPW_pte_hat = NA,
        OR_Fy1_hat = NA,
        OR_Fy0_hat = NA,
        OR_pte_hat = NA,
        IPW_Fy1_hat = NA,
        IPW_Fy0_hat = NA,
        IPW_pte_hat = NA,
        success = FALSE,
        error = as.character(e$message)
      )
    })
  }
  
  res = bind_rows(all_results)
  attr(res, "n_success") = sum(res$success, na.rm = TRUE)
  
  return(res)
}


########################################
# bootstrap
####################################

one_run_bootstrap_from_data = function(data,
                                        y_threshold = 7.95,
                                        scenario = "cc") {
  options(warn = -1)
  
  # scenarios { "cc", "cm", "mc", "mm" }
  if (scenario == "cc") {
    formA  = A ~ X1 + X2
    formY  = Y ~ X1 + X2 + A
    
  } else if (scenario == "cm") {
    formA  = A ~ X1
    formY  = Y ~ X1 + X2 + A
    
  } else if (scenario == "mc") {
    formA  = A ~ X1 + X2
    formY  = Y ~ A + X1
    
  } else if (scenario == "mm") {
    formA  = A ~ X1
    formY  = Y ~ A + X1
    
  } else {
    stop("scenario must be one of: cc, cm, mc, mm")
  }
  
  AIPW = aipw_cpm_pte(
    data = data,
    y_threshold = y_threshold,
    formA = formA,
    formY = formY
  )
  
  result = c(
    AIPW_Fy1_hat = as.numeric(unname(AIPW["FY1"])),
    AIPW_Fy0_hat = as.numeric(unname(AIPW["FY0"])),
    AIPW_pte_hat = as.numeric(unname(AIPW["FY1"])) -
      as.numeric(unname(AIPW["FY0"]))
  )
  
  return(result)
}


run_bootstrap_AIPW_on_datalist_PTE = function(data_list,
                                              y_threshold = 7.95,
                                              scenario = "cc",
                                              B = 500,
                                              n_run = length(data_list),
                                              seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  all_results = vector("list", n_run)
  
  for (i in 1:n_run) {
    
    message("Running dataset ", i, " / ", n_run)
    
    dat_i = data_list[[i]]
    n_i = nrow(dat_i)
    
    boot_results = vector("list", B)
    
    for (b in 1:B) {
      
      message("  Bootstrap ", b, " / ", B)
      
      boot_results[[b]] = tryCatch({
        
        boot_index = sample(seq_len(n_i), size = n_i, replace = TRUE)
        boot_data = dat_i[boot_index, , drop = FALSE]
        
        est_b = one_run_bootstrap_from_data(
          data = boot_data,
          y_threshold = y_threshold,
          scenario = scenario
        )
        
        data.frame(
          sim = i,
          boot = b,
          t(est_b),
          success = TRUE,
          error = NA
        )
        
      }, error = function(e) {
        
        data.frame(
          sim = i,
          boot = b,
          AIPW_Fy1_hat = NA,
          AIPW_Fy0_hat = NA,
          AIPW_pte_hat = NA,
          success = FALSE,
          error = as.character(e$message)
        )
      })
    }
    
    all_results[[i]] = bind_rows(boot_results)
  }
  
  res = bind_rows(all_results)
  attr(res, "n_success") = sum(res$success, na.rm = TRUE)
  
  return(res)
}

###############################
# sandwich
###############################
one_run_sandwich_from_data = function(data,
                                      y_threshold = 7.95,
                                      scenario = "cc") {
  
  options(warn = -1)
  
  # scenarios { "cc", "cm", "mc", "mm" }
  if (scenario == "cc") {
    formA = A ~ X1 + X2
    formY = Y ~ X1 + X2 + A
    
  } else if (scenario == "cm") {
    formA = A ~ X1
    formY = Y ~ X1 + X2 + A
    
  } else if (scenario == "mc") {
    formA = A ~ X1 + X2
    formY = Y ~ A + X1
    
  } else if (scenario == "mm") {
    formA = A ~ X1
    formY = Y ~ A + X1
    
  } else {
    stop("scenario must be one of: cc, cm, mc, mm")
  }
  
  dd = datadist(data)
  options(datadist = "dd")
  
  n = nrow(data)
  
  ## point estimates
  point = aipw_cpm_pte(
    data = data,
    y_threshold = y_threshold,
    formA = formA,
    formY = formY
  )
  
  Fy1_hat = as.numeric(unname(point["FY1"]))
  Fy0_hat = as.numeric(unname(point["FY0"]))
  pte_hat = Fy1_hat - Fy0_hat
  
  ## PS score
  ps_fit = glm(formA, data = data, family = binomial("logit"))
  Score_PS = score_logit_glm(ps_fit)
  
  ## CPM score
  mod = rms::orm(
    formY,
    family = "probit",
    data = data,
    mscore = TRUE,
    x = TRUE,
    y = TRUE
  )
  
  Score_CPM = mod$mscore
  n_nonslope = mod$non.slopes
  Score_CPM[, 1:n_nonslope] = -Score_CPM[, 1:n_nonslope]
  
  ## estimating equations for F1 and F0
  myEEs = get_EEs_AIPW_PTE(
    data = data,
    formA = formA,
    formY = formY,
    y_threshold = y_threshold,
    F1hat = Fy1_hat,
    F0hat = Fy0_hat
  )
  
  EE1 = myEEs$EE1
  EE0 = myEEs$EE0
  
  ## A matrix blocks
  A_theta1theta1 = -(1 / n) * dscore_logit_glm(ps_fit)
  A_theta2theta2 = -(1 / n) * build_info_matrix(mod)
  
  A_theta3theta3 = 1
  A_theta4theta4 = 1
  
  ## A_{F1, PS} and A_{F0, PS}
  estimate_Aapsi_ps = Aapsi_ps_AIPW_PTE(
    data = data,
    formA = formA,
    formY = formY,
    y_threshold = y_threshold,
    link = "probit"
  )
  
  A_theta3theta1 = estimate_Aapsi_ps$A1psi
  A_theta4theta1 = estimate_Aapsi_ps$A0psi
  
  ## A_{F1, CPM} and A_{F0, CPM}
  estimate_Aaxi_cpm = Aaxi_cpm_AIPW_PTE(
    data = data,
    formA = formA,
    formY = formY,
    y_threshold = y_threshold,
    link = "probit"
  )
  
  A_theta3theta2 = estimate_Aaxi_cpm$A1xi
  A_theta4theta2 = estimate_Aaxi_cpm$A0xi
  
  ## sandwich assembly
  npar_PS = ncol(Score_PS)
  npar_CPM = ncol(Score_CPM)
  Ntheta = npar_PS + npar_CPM + 2
  
  BigPhi = cbind(Score_PS, Score_CPM, EE1, EE0)
  
  A = matrix(0, Ntheta, Ntheta)
  
  idx_PS = 1:npar_PS
  idx_CPM = (npar_PS + 1):(npar_PS + npar_CPM)
  idx_F1 = npar_PS + npar_CPM + 1
  idx_F0 = npar_PS + npar_CPM + 2
  
  ## nuisance blocks
  A[idx_PS, idx_PS] = A_theta1theta1
  A[idx_CPM, idx_CPM] = A_theta2theta2
  
  ## F estimating equation blocks
  A[idx_F1, idx_PS] = A_theta3theta1
  A[idx_F0, idx_PS] = A_theta4theta1
  
  A[idx_F1, idx_CPM] = A_theta3theta2
  A[idx_F0, idx_CPM] = A_theta4theta2
  
  A[idx_F1, idx_F1] = A_theta3theta3
  A[idx_F0, idx_F0] = A_theta4theta4
  
  B = (1 / n) * crossprod(BigPhi)
  
  A_inv = tryCatch(
    solve(A),
    error = function(e) MASS::ginv(A)
  )
  
  var_theta = (1 / n) * A_inv %*% B %*% t(A_inv)
  
  var_sand_p1 = var_theta[idx_F1, idx_F1]
  var_sand_p0 = var_theta[idx_F0, idx_F0]
  cov_sand_p10 = var_theta[idx_F1, idx_F0]
  var_sand_pte = var_sand_p1 + var_sand_p0 - 2 * cov_sand_p10
  
  result = c(
    AIPW_Fy1_hat = Fy1_hat,
    AIPW_Fy0_hat = Fy0_hat,
    AIPW_pte_hat = pte_hat,
    
    var_sand_Fy1 = var_sand_p1,
    var_sand_Fy0 = var_sand_p0,
    cov_sand_Fy10 = cov_sand_p10,
    var_sand_pte = var_sand_pte,
    
    se_sand_Fy1 = sqrt(var_sand_p1),
    se_sand_Fy0 = sqrt(var_sand_p0),
    se_sand_pte = sqrt(var_sand_pte)
  )
  
  return(result)
}


run_sandwich_AIPW_on_datalist_PTE = function(data_list,
                                             y_threshold = 7.95,
                                             scenario = "cc",
                                             n_run = length(data_list),
                                             max_retries = 1L) {
  
  all_results = vector("list", n_run)
  
  for (i in seq_len(n_run)) {
    
    message("Running sandwich AIPW PTE ", i, " / ", n_run)
    
    tries = 0L
    
    repeat {
      
      tries = tries + 1L
      
      res_i = tryCatch({
        
        est_i = one_run_sandwich_from_data(
          data = data_list[[i]],
          y_threshold = y_threshold,
          scenario = scenario
        )
        
        data.frame(
          sim = i,
          t(est_i),
          success = TRUE,
          error = NA
        )
        
      }, error = function(e) {
        
        data.frame(
          sim = i,
          
          AIPW_Fy1_hat = NA,
          AIPW_Fy0_hat = NA,
          AIPW_pte_hat = NA,
          
          var_sand_Fy1 = NA,
          var_sand_Fy0 = NA,
          cov_sand_Fy10 = NA,
          var_sand_pte = NA,
          
          se_sand_Fy1 = NA,
          se_sand_Fy0 = NA,
          se_sand_pte = NA,
          
          success = FALSE,
          error = as.character(e$message)
        )
      })
      
      ok = isTRUE(res_i$success) &&
        all(is.finite(as.numeric(res_i[, c(
          "AIPW_Fy1_hat", "AIPW_Fy0_hat", "AIPW_pte_hat",
          "var_sand_Fy1", "var_sand_Fy0", "var_sand_pte",
          "se_sand_Fy1", "se_sand_Fy0", "se_sand_pte"
        )])))
      
      if (ok || tries >= max_retries) {
        all_results[[i]] = res_i
        break
      }
    }
  }
  
  res = bind_rows(all_results)
  attr(res, "n_success") = sum(res$success, na.rm = TRUE)
  
  return(res)
}

###################################
# IF based variance
#####################################
one_run_IF_PTE_from_data = function(data,
                                    y_threshold = 7.95,
                                    scenario = c("cc", "cm", "mc", "mm"),
                                    link = "probit") {
  
  scenario = match.arg(scenario)
  
  if (scenario == "cc") {
    formA = A ~ X1 + X2
    formY = Y ~ A + X1 + X2
    
  } else if (scenario == "cm") {
    formA = A ~ X1
    formY = Y ~ A + X1 + X2
    
  } else if (scenario == "mc") {
    formA = A ~ X1 + X2
    formY = Y ~ A + X1
    
  } else if (scenario == "mm") {
    formA = A ~ X1
    formY = Y ~ A + X1
  }
  
  dd = datadist(data)
  options(datadist = "dd")
  
  ## point estimates for marginal CDFs and PTE
  point = aipw_cpm_pte(
    data = data,
    y_threshold = y_threshold,
    formA = formA,
    formY = formY
  )
  
  F1_hat = as.numeric(unname(point["FY1"]))
  F0_hat = as.numeric(unname(point["FY0"]))
  pte_hat = F1_hat - F0_hat
  
  ## estimating equations / EIF components for marginal CDFs
  myEEs = get_EEs_AIPW_PTE(
    data = data,
    formA = formA,
    formY = formY,
    y_threshold = y_threshold,
    F1hat = F1_hat,
    F0hat = F0_hat,
    link = link
  )
  
  EE1 = myEEs$EE1
  EE0 = myEEs$EE0
  EE_PTE = myEEs$EE_PTE
  
  ## influence functions for marginal CDFs and PTE
  ## For CDF-scale parameter, no density division is needed
  IF1 = EE1
  IF0 = EE0
  IFPTE = EE_PTE
  
  n = nrow(data)
  
  ## IF-based variance estimates
  var_if_Fy1 = mean((IF1 - mean(IF1, na.rm = TRUE))^2, na.rm = TRUE) / n
  var_if_Fy0 = mean((IF0 - mean(IF0, na.rm = TRUE))^2, na.rm = TRUE) / n
  var_if_pte = mean((IFPTE - mean(IFPTE, na.rm = TRUE))^2, na.rm = TRUE) / n
  
  se_if_Fy1 = sqrt(var_if_Fy1)
  se_if_Fy0 = sqrt(var_if_Fy0)
  se_if_pte = sqrt(var_if_pte)
  
  ci_Fy1_l = F1_hat - 1.96 * se_if_Fy1
  ci_Fy1_u = F1_hat + 1.96 * se_if_Fy1
  
  ci_Fy0_l = F0_hat - 1.96 * se_if_Fy0
  ci_Fy0_u = F0_hat + 1.96 * se_if_Fy0
  
  ci_pte_l = pte_hat - 1.96 * se_if_pte
  ci_pte_u = pte_hat + 1.96 * se_if_pte
  
  c(
    AIPW_Fy1_hat = F1_hat,
    AIPW_Fy0_hat = F0_hat,
    AIPW_pte_hat = pte_hat,
    
    var_if_Fy1 = var_if_Fy1,
    var_if_Fy0 = var_if_Fy0,
    var_if_pte = var_if_pte,
    
    se_if_Fy1 = se_if_Fy1,
    se_if_Fy0 = se_if_Fy0,
    se_if_pte = se_if_pte,
    
    ci_Fy1_l = ci_Fy1_l,
    ci_Fy1_u = ci_Fy1_u,
    ci_Fy0_l = ci_Fy0_l,
    ci_Fy0_u = ci_Fy0_u,
    ci_pte_l = ci_pte_l,
    ci_pte_u = ci_pte_u
  )
}


run_IF_AIPW_on_datalist_PTE = function(data_list,
                                       y_threshold = 7.95,
                                       scenario = c("cc", "cm", "mc", "mm"),
                                       link = "probit",
                                       n_run = length(data_list)) {
  
  scenario = match.arg(scenario)
  n_run = min(n_run, length(data_list))
  
  result_names = c(
    "AIPW_Fy1_hat", "AIPW_Fy0_hat", "AIPW_pte_hat",
    "var_if_Fy1", "var_if_Fy0", "var_if_pte",
    "se_if_Fy1", "se_if_Fy0", "se_if_pte",
    "ci_Fy1_l", "ci_Fy1_u",
    "ci_Fy0_l", "ci_Fy0_u",
    "ci_pte_l", "ci_pte_u"
  )
  
  results = matrix(NA_real_, nrow = n_run, ncol = length(result_names))
  colnames(results) = result_names
  
  for (i in seq_len(n_run)) {
    
    cat("Running IF AIPW PTE dataset", i, "/", n_run, "\n")
    
    out = tryCatch(
      {
        one_run_IF_PTE_from_data(
          data = data_list[[i]],
          y_threshold = y_threshold,
          scenario = scenario,
          link = link
        )
      },
      error = function(e) {
        cat("Dataset", i, "ERROR:", conditionMessage(e), "\n")
        setNames(rep(NA_real_, length(result_names)), result_names)
      }
    )
    
    results[i, ] = out[result_names]
  }
  
  results = as.data.frame(results)
  attr(results, "n_success") = sum(complete.cases(results))
  
  return(results)
}