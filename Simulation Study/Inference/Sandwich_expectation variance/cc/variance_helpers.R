###########################
# PS part
#############################
score_logit_glm = function(fit) {
  
  mf = model.frame(fit)
  X = model.matrix(fit)
  y = model.response(mf)
  
  eta = drop(X %*% coef(fit))
  mu   = plogis(eta)  
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
# q1p, q0p part
###################################
get_EEs = function(data, formA, formY, q1_p, q0_p, p, link = c("logistic", "probit")) {
  link = match.arg(link)
  g = rms::orm( formY, family = link, data = data, mscore = TRUE, x = TRUE, y= TRUE )
  alpha_hat = -g$coefficients[1:(length(g$yunique) - 1)]
  
  g1_vals = cond_cdf_onemodel(data = data,y_star = q1_p, g = g, alpha_hat = alpha_hat,A =1, link = link)
  g0_vals = cond_cdf_onemodel(data = data,y_star = q0_p, g = g, alpha_hat = alpha_hat,A = 0, link = link )
  cF1_q =  g1_vals$g_y; cF0_q =  g0_vals$g_y
  
  ps_mod  = glm(formA, data = data, family = binomial("logit"))
  pi1_hat = predict(ps_mod, type = "response"); pi0_hat = 1 - pi1_hat
  I1 = as.numeric(data$A == 1) ; I0 = as.numeric(data$A == 0)
  
  EE1 = (I1 / pi1_hat) * as.numeric(data$Y <= q1_p) -
    ((I1 - pi1_hat) / pi1_hat) * cF1_q - p
  EE0 = (I0 / pi0_hat) * as.numeric(data$Y <= q0_p) -
    ((I0 - pi0_hat) / pi0_hat) * cF0_q - p
  
  list(EE1 = EE1, EE0 = EE0)
}



counterfactual_dist = function(data, formA, formY, q1_p, q0_p, p, link = c("logistic", "probit")) {
  link = match.arg(link)
  ps_fit  = glm(formA, data = data, family = binomial("logit"))
  pi1 = predict(ps_fit, type = "response");  pi0  = 1 - pi1                            
  IA1 = as.numeric(data$A == 1); IA0 = as.numeric(data$A == 0)
  weight1 = IA1 / pi1 ; weight0 = IA0 / pi0 
  
  g = orm(formY, family = link, x=T, y = T,data = data)
  alpha_hat = -g$coefficients[1:(length(g$yunique) - 1)]

  n = nrow(data); DR1  = numeric(n);DR0  = numeric(n)
  for (i in seq_len(n)) {
    y_star = data$Y[i]
    # A = 1
    IY1 = as.numeric(data$Y <= y_star)
    g1_vals= cond_cdf_onemodel(data = data,y_star = y_star, g = g, alpha_hat = alpha_hat,A =1, link = link )
    DR1[i] = mean( (IA1 / pi1 ) * (IY1 - g1_vals$g_y) + g1_vals$g_y)
    # A = 0
    IY0 = as.numeric(data$Y <= y_star)
    g0_vals= cond_cdf_onemodel(data = data,y_star = y_star, g = g, alpha_hat = alpha_hat,A =0, link = link )
    DR0[i] = mean( (IA0 / pi0) * (IY0 - g0_vals$g_y) + g0_vals$g_y)
  }
  
  ord = order(data$Y); y_levels = Ygrid1 = Ygrid0 = Ygrid = data$Y[ord]
  DR1_sorted = DR1[ord]; DR0_sorted = DR0[ord]
  DR1_sorted = cummax(pmin(pmax(DR1_sorted, 0), 1))
  DR0_sorted = cummax(pmin(pmax(DR0_sorted, 0), 1))
  
  # Compute the marginal CDF by averaging over all units
  mF1 = DR1_sorted; mF0 = DR0_sorted
  bw = bw.SJ(y_levels);gridsize = 5001
  
  lp1_cdf   = KernSmooth::locpoly(y_levels, mF1, drv = 0, bandwidth = bw, gridsize = gridsize)
  lp1_pdf   = KernSmooth::locpoly(y_levels, mF1, drv = 1, bandwidth = bw, gridsize = gridsize)
  lp1_pdf_p = KernSmooth::locpoly(y_levels, mF1, drv = 2, bandwidth = bw, gridsize = gridsize)
  
  lp0_cdf   = KernSmooth::locpoly(y_levels, mF0, drv = 0, bandwidth = bw, gridsize = gridsize)
  lp0_pdf   = KernSmooth::locpoly(y_levels, mF0, drv = 1, bandwidth = bw, gridsize = gridsize)
  lp0_pdf_p = KernSmooth::locpoly(y_levels, mF0, drv = 2, bandwidth = bw, gridsize = gridsize)
  
  F1_q1p    = approx(lp1_cdf$x, lp1_cdf$y, xout = q1_p, rule = 2)$y  
  f1_q      = approx(lp1_pdf$x, lp1_pdf$y, xout = q1_p, rule = 2)$y  
  f1prime_q = approx(lp1_pdf_p$x, lp1_pdf_p$y, xout = q1_p, rule = 2)$y  
  
  F0_q0p    = approx(lp0_cdf$x, lp0_cdf$y, xout = q0_p, rule = 2)$y  
  f0_q      = approx(lp0_pdf$x, lp0_pdf$y, xout = q0_p, rule = 2)$y  
  f0prime_q = approx(lp0_pdf_p$x, lp0_pdf_p$y, xout = q0_p, rule = 2)$y  
  
  eps = .Machine$double.eps^0.5
  y_min = y_levels[1] ; y_max = tail(y_levels, 1)
  zap_outside = function(x, at) {
    x[is.na(x)] = 0
    x[at < y_min | at > y_max] = 0
    x
  }
  
  F1_q1p = zap_outside(F1_q1p,  q1_p); F0_q0p  = zap_outside(F0_q0p,  q0_p)
  f1_q = zap_outside(f1_q, q1_p); f0_q = zap_outside(f0_q, q0_p)
  f1_q = pmax(f1_q, eps); f0_q = pmax(f0_q, eps)
  
  f1prime_q = zap_outside(f1prime_q, q1_p)
  f0prime_q = zap_outside(f0prime_q, q0_p)
  list(
    F1_q1p = F1_q1p, F0_q0p = F0_q0p,#marginal CDF：F1/F0
    f1_q = f1_q, f0_q = f0_q, #first derivative：∂θ3 F = f
    f1prime_q = f1prime_q, f0prime_q = f0prime_q # second derivative：∂θ3 f
  )
}

##################################################
# A1xi and A0xi
##################################################
Aaxi_cpm = function(data, formY, formA, q1_p, q0_p,
                            treat_var = "A",
                            link = c("probit", "logistic")) {
  link = match.arg(link)
  cpm_fit = rms::orm( formY, family = link, data = data, mscore = TRUE, x = TRUE, y = TRUE )
  ps_fit = glm(formA, data = data, family = binomial("logit"))
  
  n = nrow(data)
  y_levels = cpm_fit$yunique; n_alpha = length(y_levels) - 1 
  coef_all = cpm_fit$coefficients
  alpha_hat = -coef_all[1:n_alpha]; beta_hat  = coef_all[-(1:n_alpha)]
  
  X_design = cpm_fit$x
  X1_design = X_design; X0_design = X_design
  X1_design[, treat_var] = 1; X0_design[, treat_var] = 0
  
  pi1_hat = predict(ps_fit, type = "response"); pi0_hat = 1 - pi1_hat
  A1 = as.numeric(data[[treat_var]] == 1)
  A0 = as.numeric(data[[treat_var]] == 0)
  
  compute_one_Axi = function(q_p, X_a, w_a) {
    # outside support: derivative = 0
    if (q_p < y_levels[1] || q_p >= y_levels[length(y_levels)]) {
      out = rep(0, length(coef_all))
      names(out) = names(coef_all)
      return(out)
    }
    j = max(which(y_levels <= q_p))
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
    
    Axi_hat = colMeans(w_a * grad_xi)
    names(Axi_hat) = c(names(coef_all)[1:n_alpha], names(beta_hat))
    
    return(Axi_hat)
  }
  w1 = A1 / pi1_hat - 1; w0 = A0 / pi0_hat - 1
  A1xi_hat = compute_one_Axi(q_p = q1_p, X_a = X1_design, w_a = w1)
  A0xi_hat = compute_one_Axi(q_p = q0_p, X_a = X0_design, w_a = w0)
  
  return(list( A1xi = A1xi_hat, A0xi = A0xi_hat ))
}


##################################################
# A1psi and A0psi
##################################################
Aapsi_ps = function(data, formY, formA, q1_p, q0_p,  link = c("probit", "logistic"))  {
  
  link = match.arg(link); n = nrow(data);  ord = order(data$Y)
  ps_fit = glm(formA, data = data, family = binomial("logit"))
  X_tilde = model.matrix(ps_fit)
  pi1 = as.numeric(predict(ps_fit, type = "response")); pi0 = 1 - pi1
  IA1 = as.numeric(data$A == 1);IA0 = as.numeric(data$A == 0)
  IY_q1 = as.numeric(data$Y <= q1_p);  IY_q0 = as.numeric(data$Y <= q0_p)
  
  g = orm(formY, family = link, x=T, y = T, data = data)
  alpha_hat = -g$coefficients[1:(length(g$yunique) - 1)]
  y_levels = g$yunique
  g1_vals = cond_cdf_onemodel(data = data,y_star = q1_p, g = g, alpha_hat = alpha_hat,A =1, link = link)
  g0_vals = cond_cdf_onemodel(data = data,y_star = q0_p, g = g, alpha_hat = alpha_hat,A = 0, link = link )
  m1_q =  g1_vals$g_y; m0_q =  g0_vals$g_y
  
  w1 = IA1 * ((1 - pi1) / pi1) * (IY_q1 - m1_q)
  w0 = -IA0 * ((1 - pi0) / pi0)  * (IY_q0 - m0_q)
  
  A1psi = colMeans(X_tilde * w1)
  A0psi = colMeans(X_tilde * w0)
  
  return(list(A1psi = A1psi,A0psi = A0psi))
}


#######################################
## Summarize simulation output
#######################################
summarize_sim_sandwich = function(out, n) {
  n_ok = attr(out, "n_success")
  
  truth_q1  = m_Y1_0.5
  truth_q0  = m_Y0_0.5
  truth_qte = QTE_0.5_true
  
  summ = data.frame(
    target        = c("q1", "q0", "QTE"),
    true          = c(truth_q1, truth_q0, truth_qte),
    mean_est      = c(mean(out$q1_hat, na.rm = TRUE),
                      mean(out$q0_hat, na.rm = TRUE),
                      mean(out$qte_hat, na.rm = TRUE)),
    bias          = c(mean(out$q1_hat - truth_q1, na.rm = TRUE),
                      mean(out$q0_hat - truth_q0, na.rm = TRUE),
                      mean(out$qte_hat - truth_qte, na.rm = TRUE)),
    emp_var       = c(var(out$q1_hat, na.rm = TRUE),
                      var(out$q0_hat, na.rm = TRUE),
                      var(out$qte_hat, na.rm = TRUE)),
    mean_sand_var = c(mean(out$var_sand_q1, na.rm = TRUE),
                      mean(out$var_sand_q0, na.rm = TRUE),
                      mean(out$var_sand_qte, na.rm = TRUE)) ,
    n_success     = n_ok
  )
  
  summ
}

#########################################
# take derivation first
###########################################
estimate_cond_density_locpoly2 = function(data, formY, q1_p, q0_p, link = c("probit", "logistic"), gridsize = 5001) {
  link = match.arg(link)
  n = nrow(data); ord = order(data$Y); Ygrid = data$Y[ord]
  
  # OR model
  g = rms::orm(formY, family = link, x = TRUE, y = TRUE, data = data)
  alpha_hat = -g$coefficients[1:(length(g$yunique) - 1)]
  
  OR1 = matrix(NA_real_, n, n); OR0 = matrix(NA_real_, n, n)
  for (j in seq_len(n)) {
    y_star = data$Y[j]
    g1_vals = cond_cdf_onemodel(data = data, y_star = y_star, g = g, 
                                alpha_hat = alpha_hat, A = 1,link = link )
    g0_vals = cond_cdf_onemodel(data = data, y_star = y_star, g = g, 
                                alpha_hat = alpha_hat, A = 0,link = link )
    OR1[, j] = g1_vals$g_y; OR0[, j] = g0_vals$g_y
  }
  
  c_cdf1 = OR1[, ord, drop = FALSE]; c_cdf0 = OR0[, ord, drop = FALSE]
  bw = stats::bw.SJ(Ygrid)
  f_cond1 = numeric(n); f_cond0 = numeric(n)
  
  for (i in seq_len(n)) {
    Fi1 = as.numeric(c_cdf1[i, ]); Fi0 = as.numeric(c_cdf0[i, ])
    lp_pdf1 = KernSmooth::locpoly( x = Ygrid, y = Fi1, drv = 1, bandwidth = bw, gridsize = gridsize )
    lp_pdf0 = KernSmooth::locpoly( x = Ygrid, y = Fi0, drv = 1, bandwidth = bw, gridsize = gridsize )
    
    f_cond1[i] = approx(lp_pdf1$x, lp_pdf1$y, xout = q1_p, rule = 2)$y
    f_cond0[i] = approx(lp_pdf0$x, lp_pdf0$y, xout = q0_p, rule = 2)$y
  }
  
  return(list( f_cond1 = f_cond1, f_cond0 = f_cond0, bw = bw  ))
}


estimate_Aqq_deriv_first = function(data, formA, formY, q1_p, q0_p,
                                    link = c("probit", "logistic"),
                                    h = NULL, gridsize = 5001) {
  
  link = match.arg(link) ; n = nrow(data)
  
  ps_fit = glm(formA, data = data, family = binomial("logit"))
  pi1 = predict(ps_fit, type = "response"); pi0 = 1 - pi1
  A1 = as.numeric(data$A == 1) ; A0 = as.numeric(data$A == 0)
  if (is.null(h)) {
   h = 1.06 * stats::sd(data$Y) * n^(-1/5)
   if (!is.finite(h) || h <= 0) h = 0.5
}

S1 = plogis((q1_p - data$Y) / h); S0 = plogis((q0_p - data$Y) / h)
delta1_approx = (1 / h) * S1 * (1 - S1); delta0_approx = (1 / h) * S0 * (1 - S0)
term1_A11 = -mean(A1 / pi1 * delta1_approx) ; term1_A00 = -mean(A0 / pi0 * delta0_approx)
dens_obj = estimate_cond_density_locpoly2( data = data, formY = formY, q1_p = q1_p, q0_p = q0_p, link = link, gridsize = gridsize )
f_cond1 = dens_obj$f_cond1 ; f_cond0 = dens_obj$f_cond0 ; bw = dens_obj$bw
term2_A11 = mean((pi1 - A1) / pi1 * f_cond1); term2_A00 = mean((pi0 - A0) / pi0 * f_cond0)

A11_hat = term1_A11 + term2_A11; A00_hat = term1_A00 + term2_A00

return(list( A11 = A11_hat, A00 = A00_hat ))
}

##################################
# for bootstrap
###################################
estimate_once_from_data= function(data, p_star = 0.5, scenario = "cc") {
  options(warn = -1)
  # scenarios: cc, cm, mc, mm
  if (scenario == "cc") {
    formA= A ~ X1 + X2
    formY= Y ~ X1 + X2 + A
  } else if (scenario == "cm") {
    formA= A ~ X1
    formY= Y ~ X1 + X2 + A
  } else if (scenario == "mc") {
    formA= A ~ X1 + X2
    formY= Y ~ A + X1
  } else if (scenario == "mm") {
    formA= A ~ X1
    formY= Y ~ A + X1
  } else {
    stop("scenario must be one of: cc, cm, mc, mm")
  }
  
  out= tryCatch({
    result_AIPW = aipw_cpm_cdf(data, p_star, formA = formA, formY = formY, link = "probit")
    result_aipw_cpm_q= aipw_cpm_q(data, p_star, formA = formA, formY = formY, link = "probit")
    c(
      aipw_cpm_q_q1 = as.numeric(result_aipw_cpm_q["q1"]),
      aipw_cpm_q_q0 = as.numeric(result_aipw_cpm_q["q0"]),
      aipw_cpm_q_QTE= as.numeric(result_aipw_cpm_q["QTE"]),
      aipw_q1 = as.numeric(result_AIPW["q1"]),
      aipw_q0 = as.numeric(result_AIPW["q0"]),
      aipw_QTE= as.numeric(result_AIPW["QTE"])
    )
  }, error = function(e) {
    c(
      aipw_cpm_q_q1 = NA_real_, aipw_cpm_q_q0 = NA_real_,
      aipw_cpm_q_QTE= NA_real_, aipw_q1 = NA_real_,
      aipw_q0 = NA_real_, aipw_QTE = NA_real_
    )
  })
  
  return(out)
}



bootstrap_one_dataset = function(data, 
                                 p_star = 0.5, 
                                 scenario = "cc",
                                 B = 300,
                                 seed = NULL,
                                 verbose = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  
  n = nrow(data)
  theta_hat = estimate_once_from_data( data = data,  p_star = p_star,  scenario = scenario)
  est_names = names(theta_hat)
  
  # bootstrap replicates
  boot_mat = matrix(NA_real_, nrow = B, ncol = length(theta_hat))
  colnames(boot_mat) = est_names
  
  for (b in seq_len(B)) {
    idx = sample.int(n, size = n, replace = TRUE)
    boot_data = data[idx, , drop = FALSE]
    
    boot_mat[b, ] = estimate_once_from_data(
      data = boot_data,
      p_star = p_star,
      scenario = scenario
    )
    
    if (verbose && b %% 10 == 0) { cat("completed bootstrap", b, "of", B, "\n") }
  }
  
  return(list(
    point_est = theta_hat,
    boot_reps = boot_mat
  ))
}


###################################
# for IF based variance
#####################################
one_run_IF_from_data = function(data,
                                scenario = c("cc", "cm", "mc", "mm"),
                                p_star = 0.5) {
  
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
  
  dd = datadist(data); options(datadist = "dd")
  point = aipw_cpm_q(data  = data, p = p_star,formA = formA,formY = formY, link  = "probit")
  q1_hat  = unname(point["q1"]); q0_hat  = unname(point["q0"])
  qte_hat = q1_hat - q0_hat
  
  ## estimating equations
  myEEs = get_EEs( data  = data, formA = formA, formY = formY,
                   q1_p  = q1_hat, q0_p = q0_hat,p = p_star,link  = "probit")
  EE1 = myEEs$EE1; EE0 = myEEs$EE0
  
  fa = counterfactual_dist( data = data, formA = formA, formY = formY,
                            q1_p  = q1_hat, q0_p  = q0_hat, p = p_star, link  = "probit")
  f1 = fa$f1_q; f0 = fa$f0_q
  
  ## influence functions
  IF1 = -EE1 / f1; IF0 = -EE0 / f0; IFQ = IF1 - IF0
  n = nrow(data)
  var_if_q1 = mean((IF1 - mean(IF1))^2)/n 
  var_if_q0= mean((IF0 - mean(IF0))^2)/n 
  var_if_qte = mean((IFQ - mean(IFQ))^2)/n
  
  se_if_q1 = sqrt(var_if_q1); se_if_q0 = sqrt(var_if_q0); se_if_qte = sqrt(var_if_qte)
  ci_q1_l  = q1_hat  - 1.96 * se_if_q1
  ci_q1_u  = q1_hat  + 1.96 * se_if_q1
  ci_q0_l  = q0_hat  - 1.96 * se_if_q0
  ci_q0_u  = q0_hat  + 1.96 * se_if_q0
  ci_qte_l = qte_hat - 1.96 * se_if_qte
  ci_qte_u = qte_hat + 1.96 * se_if_qte
  
  c(
    q1_hat = q1_hat,
    q0_hat = q0_hat,
    qte_hat = qte_hat,
    
    var_if_q1 = var_if_q1,
    var_if_q0 = var_if_q0,
    var_if_qte = var_if_qte,
    
    se_if_q1 = se_if_q1,
    se_if_q0 = se_if_q0,
    se_if_qte = se_if_qte,
    
    ci_q1_l = ci_q1_l,
    ci_q1_u = ci_q1_u,
    ci_q0_l = ci_q0_l,
    ci_q0_u = ci_q0_u,
    ci_qte_l = ci_qte_l,
    ci_qte_u = ci_qte_u
  )
}


run_IF_variance_on_datalist = function(data_list,
                                       scenario = c("cc", "cm", "mc", "mm"),
                                       p_star = 0.5,
                                       n_run = length(data_list)) {
  
  scenario = match.arg(scenario)
  n_run = min(n_run, length(data_list))
  
  result_names = c(
    "q1_hat", "q0_hat", "qte_hat",
    "var_if_q1", "var_if_q0", "var_if_qte",
    "se_if_q1", "se_if_q0", "se_if_qte",
    "ci_q1_l", "ci_q1_u",
    "ci_q0_l", "ci_q0_u",
    "ci_qte_l", "ci_qte_u"
  )
  
  results = matrix(NA_real_, nrow = n_run, ncol = length(result_names))
  colnames(results) = result_names
  
  for (i in seq_len(n_run)) {
    cat("Running dataset", i, "\n")
    
    out = tryCatch(
      {
        one_run_IF_from_data(
          data     = data_list[[i]],
          scenario = scenario,
          p_star   = p_star
        )
      },
      error = function(e) {
        cat("Dataset", i, "ERROR:", conditionMessage(e), "\n")
        setNames(rep(NA_real_, length(result_names)), result_names)
      }
    )
    
    results[i, ] = out[result_names]
  }
  
  as.data.frame(results)
}


###################################
# for sandwich variance--take expectation first
#####################################
one_run_sandwich_expectation_onedata = function(data,
                                     scenario = c("cc", "cm", "mc", "mm"),
                                     p_star = 0.5) {
  
  scenario = match.arg(scenario)
  if (scenario== "cc") {
    formA= A ~ X1 + X2             # correct PS
    formY= Y ~ A + X1 + X2          # correct outcome
  } else if (scenario== "cm") {
    formA= A ~ X1           # mis-specified PS
    formY= Y ~ A + X1 + X2          # correct outcome
  } else if (scenario== "mc") {
    formA= A ~ X1 + X2             # correct PS
    formY= Y ~ A + X1           # mis-specified outcome
  } else if (scenario== "mm") {
    formA= A ~  X1              # mis-specified PS
    formY= Y ~ A + X1             # mis-specified outcome
  } else {
    stop("scenario must be one of: cc, cm, mc, mm")
  }
  
  dd   = datadist(data);  options(datadist = "dd")
  
  ##  point estimates
  point   = aipw_cpm_q(data, p_star, formA = formA, formY = formY, link = "probit")
  q1_hat  = unname(point["q1"]); q0_hat  = unname(point["q0"])
  qte_hat = q1_hat - q0_hat
  
  ## score pieces for nuisance parameters
  ## PS model
  ps_fit = glm(formA, data = data, family = binomial("logit"))
  Score_PS = score_logit_glm(ps_fit) # n by the number of parameters in PS model
  
  ## CPM model
  mod = rms::orm( formY, family = "probit", data = data, mscore = TRUE, x = TRUE, y= TRUE )
  Score_CPM = mod$mscore   # n by (non-slope + slope)
  n_nonslope = mod$non.slopes
  Score_CPM[, 1:n_nonslope] = -Score_CPM[, 1:n_nonslope] # Liu Qi, 2017 Sec2.2
  
  ## estimating equations for q1 and q0
  myEEs = get_EEs( data  = data,formA = formA, formY = formY, q1_p  = q1_hat, q0_p  = q0_hat,
                   p  = p_star, link = "probit")
  EE1 = myEEs$EE1; EE0 = myEEs$EE0 # n by 1
  
  ## build Jacobian blocks
  A_theta1theta1 = ( 1/nrow(data) ) * dscore_logit_glm(ps_fit)
  A_theta2theta2 = -( 1/nrow(data) ) * build_info_matrix(mod) #(non-slope + slope) by (non-slope + slope)
  A_3344 = counterfactual_dist(data  = data, formA = formA, formY = formY, 
                               q1_p = q1_hat, q0_p = q0_hat, p = p_star, link = "probit")
  A_theta3theta3 = -A_3344$f1_q ; A_theta4theta4 = -A_3344$f0_q
  ## A_{theta3, theta1} and A_{theta4, theta1}
  estimate_Aapsi_ps = Aapsi_ps( data = data, formY = formY,formA = formA, q1_p = q1_hat,q0_p = q0_hat, link = "probit" )
  A_theta3theta1 = estimate_Aapsi_ps$A1psi; A_theta4theta1 = estimate_Aapsi_ps$A0psi
  ## A_{theta3, theta2} and A_{theta4, theta2}
  estimate_Aaxi_cpm = Aaxi_cpm( data = data, formY = formY,formA = formA, q1_p = q1_hat,q0_p = q0_hat, link = "probit" )
  A_theta3theta2 = estimate_Aaxi_cpm$A1xi; A_theta4theta2 = estimate_Aaxi_cpm$A0xi
  
  ## sandwich assembly
  npar.PS  = ncol(Score_PS); npar.CPM = ncol(Score_CPM)
  N = nrow(data); Ntheta = npar.PS + npar.CPM + 2 # 2 is q1p and q0p
  BigPhi = cbind(Score_PS, Score_CPM, EE1, EE0)
  
  A = matrix(0, Ntheta, Ntheta)
  A[1:npar.PS, 1:npar.PS] = A_theta1theta1
  A[npar.PS + (1:npar.CPM), npar.PS + (1:npar.CPM)] = A_theta2theta2
  A[npar.PS + npar.CPM + 1, 1:npar.PS] = A_theta3theta1
  A[npar.PS + npar.CPM + 1, npar.PS + (1:npar.CPM)] = A_theta3theta2
  A[npar.PS + npar.CPM + 2, 1:npar.PS] = A_theta4theta1
  A[npar.PS + npar.CPM + 2, npar.PS + (1:npar.CPM)] = A_theta4theta2
  A[npar.PS + npar.CPM + 1, npar.PS + npar.CPM + 1] = A_theta3theta3
  A[npar.PS + npar.CPM + 2, npar.PS + npar.CPM + 2] = A_theta4theta4
  
  B = (1 / N) * crossprod(BigPhi)
  A_inv = tryCatch( solve(A),error = function(e) MASS::ginv(A))
  var_theta = (1 / N) * A_inv %*% B %*% t(A_inv)
  
  var_sand_q1  = var_theta[Ntheta - 1, Ntheta - 1]
  var_sand_q0  = var_theta[Ntheta, Ntheta]
  cov_sand_q10 = var_theta[Ntheta, Ntheta - 1]
  var_sand_qte = var_sand_q1 + var_sand_q0 - 2 * cov_sand_q10
  
  c( q1_hat = q1_hat, q0_hat = q0_hat, qte_hat = qte_hat,
    var_sand_q1  = var_sand_q1, var_sand_q0  = var_sand_q0,
    cov_sand_q10 = cov_sand_q10, var_sand_qte = var_sand_qte )
  
}



run_sandwich_expectation_on_datalist = function(data_list,
                                    scenario = c("cc", "cm", "mc", "mm"),
                                    p_star = 0.5,
                                    n_run = NULL) {
  
  scenario = match.arg(scenario)
  
  res_list = vector("list", n_run)
  for (i in seq_len(n_run)) {
    cat("Running dataset", i, "of", n_run, "\n")
    
    out = tryCatch(
      {
        tmp = one_run_sandwich_expectation_onedata(
          data = data_list[[i]],
          scenario = scenario,
          p_star  = p_star
        )
        as.data.frame(as.list(tmp))
      },
      error = function(e) {
        cat("Dataset", i, "ERROR\n")
        data.frame(
          q1_hat       = NA_real_,
          q0_hat       = NA_real_,
          qte_hat      = NA_real_,
          var_sand_q1  = NA_real_,
          var_sand_q0  = NA_real_,
          cov_sand_q10 = NA_real_,
          var_sand_qte = NA_real_,
          stringsAsFactors = FALSE
        )
      }
    )
    
   
    out$se_sand_q1  = sqrt(out$var_sand_q1)
    out$se_sand_q0  = sqrt(out$var_sand_q0)
    out$se_sand_qte = sqrt(out$var_sand_qte)
    
    out$ci_q1_l  = out$q1_hat  - 1.96 * out$se_sand_q1
    out$ci_q1_u  = out$q1_hat  + 1.96 * out$se_sand_q1
    
    out$ci_q0_l  = out$q0_hat  - 1.96 * out$se_sand_q0
    out$ci_q0_u  = out$q0_hat  + 1.96 * out$se_sand_q0
    
    out$ci_qte_l = out$qte_hat - 1.96 * out$se_sand_qte
    out$ci_qte_u = out$qte_hat + 1.96 * out$se_sand_qte
    
    res_list[[i]] = out
  }
  
  bind_rows(res_list) %>%
    dplyr::select(
      q1_hat, q0_hat, qte_hat,
      var_sand_q1, var_sand_q0, cov_sand_q10, var_sand_qte,
      se_sand_q1, se_sand_q0, se_sand_qte,
      ci_q1_l, ci_q1_u, ci_q0_l, ci_q0_u, ci_qte_l, ci_qte_u
    )
}



summarize_sandwich_results = function(res_df,
                                       truth_q1,
                                       truth_q0,
                                       truth_qte) {
  
  res_df2 = res_df %>%
    mutate(
      success = !is.na(q1_hat),
      
      cover_q1  = !is.na(ci_q1_l)  & !is.na(ci_q1_u)  & (ci_q1_l  <= truth_q1  & ci_q1_u  >= truth_q1),
      cover_q0  = !is.na(ci_q0_l)  & !is.na(ci_q0_u)  & (ci_q0_l  <= truth_q0  & ci_q0_u  >= truth_q0),
      cover_qte = !is.na(ci_qte_l) & !is.na(ci_qte_u) & (ci_qte_l <= truth_qte & ci_qte_u >= truth_qte),
      
      bias_q1  = q1_hat  - truth_q1,
      bias_q0  = q0_hat  - truth_q0,
      bias_qte = qte_hat - truth_qte
    )
  
  summary_df = res_df2 %>%
    summarise(
      n_total   = n(),
      n_success = sum(success),
      n_error   = sum(!success),
      
      mean_q1_hat  = mean(q1_hat,  na.rm = TRUE),
      mean_q0_hat  = mean(q0_hat,  na.rm = TRUE),
      mean_qte_hat = mean(qte_hat, na.rm = TRUE),
      
      emp_var_q1  = var(q1_hat,  na.rm = TRUE),
      emp_var_q0  = var(q0_hat,  na.rm = TRUE),
      emp_var_qte = var(qte_hat, na.rm = TRUE),
      
      mean_sand_var_q1  = mean(var_sand_q1,  na.rm = TRUE),
      mean_sand_var_q0  = mean(var_sand_q0,  na.rm = TRUE),
      mean_sand_var_qte = mean(var_sand_qte, na.rm = TRUE),
      
      mean_bias_q1  = mean(bias_q1,  na.rm = TRUE),
      mean_bias_q0  = mean(bias_q0,  na.rm = TRUE),
      mean_bias_qte = mean(bias_qte, na.rm = TRUE),
      
      sd_q1  = sd(q1_hat,  na.rm = TRUE),
      sd_q0  = sd(q0_hat,  na.rm = TRUE),
      sd_qte = sd(qte_hat, na.rm = TRUE),
      
      coverage_q1  = mean(cover_q1,  na.rm = TRUE),
      coverage_q0  = mean(cover_q0,  na.rm = TRUE),
      coverage_qte = mean(cover_qte, na.rm = TRUE)
    )
  
  list(
    replicate_results = res_df2,
    summary = summary_df
  )
}


###################################
# for sandwich variance--take derivation first
#####################################
one_run_sandwich_derivation_onedata = function(data,
                                                scenario = c("cc", "cm", "mc", "mm"),
                                                p_star = 0.5) {
  
  scenario = match.arg(scenario)
  if (scenario== "cc") {
    formA= A ~ X1 + X2             # correct PS
    formY= Y ~ A + X1 + X2          # correct outcome
  } else if (scenario== "cm") {
    formA= A ~ X1           # mis-specified PS
    formY= Y ~ A + X1 + X2          # correct outcome
  } else if (scenario== "mc") {
    formA= A ~ X1 + X2             # correct PS
    formY= Y ~ A + X1           # mis-specified outcome
  } else if (scenario== "mm") {
    formA= A ~  X1              # mis-specified PS
    formY= Y ~ A + X1             # mis-specified outcome
  } else {
    stop("scenario must be one of: cc, cm, mc, mm")
  }
  
  dd   = datadist(data);  options(datadist = "dd")
  
  ##  point estimates
  point   = aipw_cpm_q(data, p_star, formA = formA, formY = formY, link = "probit")
  q1_hat  = unname(point["q1"]); q0_hat  = unname(point["q0"])
  qte_hat = q1_hat - q0_hat
  
  ## score pieces for nuisance parameters
  ## PS model
  ps_fit = glm(formA, data = data, family = binomial("logit"))
  Score_PS = score_logit_glm(ps_fit) # n by the number of parameters in PS model
  
  ## CPM model
  mod = rms::orm( formY, family = "probit", data = data, mscore = TRUE, x = TRUE, y= TRUE )
  Score_CPM = mod$mscore   # n by (non-slope + slope)
  n_nonslope = mod$non.slopes
  Score_CPM[, 1:n_nonslope] = -Score_CPM[, 1:n_nonslope] # Liu Qi, 2017 Sec2.2
  
  ## estimating equations for q1 and q0
  myEEs = get_EEs( data  = data,formA = formA, formY = formY, q1_p  = q1_hat, q0_p  = q0_hat,
                   p  = p_star, link = "probit")
  EE1 = myEEs$EE1; EE0 = myEEs$EE0 # n by 1
  
  ## build Jacobian blocks
  A_theta1theta1 = -( 1/nrow(data) ) * dscore_logit_glm(ps_fit)
  A_theta2theta2 = -( 1/nrow(data) ) * build_info_matrix(mod) #(non-slope + slope) by (non-slope + slope)
  
  A_3344 = estimate_Aqq_deriv_first(data  = data, formA = formA, formY = formY, 
                                    q1_p = q1_hat, q0_p = q0_hat, link = "probit")
  A_theta3theta3 = A_3344$A11 ; A_theta4theta4 = A_3344$A00
  
  
  ## A_{theta3, theta1} and A_{theta4, theta1}
  estimate_Aapsi_ps = Aapsi_ps( data = data, formY = formY,formA = formA, q1_p = q1_hat,q0_p = q0_hat, link = "probit" )
  A_theta3theta1 = estimate_Aapsi_ps$A1psi; A_theta4theta1 = estimate_Aapsi_ps$A0psi
  ## A_{theta3, theta2} and A_{theta4, theta2}
  estimate_Aaxi_cpm = Aaxi_cpm( data = data, formY = formY,formA = formA, q1_p = q1_hat,q0_p = q0_hat, link = "probit" )
  A_theta3theta2 = estimate_Aaxi_cpm$A1xi; A_theta4theta2 = estimate_Aaxi_cpm$A0xi
  
  ## sandwich assembly
  npar.PS  = ncol(Score_PS); npar.CPM = ncol(Score_CPM)
  N = nrow(data); Ntheta = npar.PS + npar.CPM + 2 # 2 is q1p and q0p
  BigPhi = cbind(Score_PS, Score_CPM, EE1, EE0)
  
  A = matrix(0, Ntheta, Ntheta)
  A[1:npar.PS, 1:npar.PS] = A_theta1theta1
  A[npar.PS + (1:npar.CPM), npar.PS + (1:npar.CPM)] = A_theta2theta2
  A[npar.PS + npar.CPM + 1, 1:npar.PS] = A_theta3theta1
  A[npar.PS + npar.CPM + 1, npar.PS + (1:npar.CPM)] = A_theta3theta2
  A[npar.PS + npar.CPM + 2, 1:npar.PS] = A_theta4theta1
  A[npar.PS + npar.CPM + 2, npar.PS + (1:npar.CPM)] = A_theta4theta2
  A[npar.PS + npar.CPM + 1, npar.PS + npar.CPM + 1] = A_theta3theta3
  A[npar.PS + npar.CPM + 2, npar.PS + npar.CPM + 2] = A_theta4theta4
  
  B = (1 / N) * crossprod(BigPhi)
  A_inv <- tryCatch( solve(A),error = function(e) MASS::ginv(A))
  var_theta = (1 / N) * A_inv %*% B %*% t(A_inv)
  
  var_sand_q1  = var_theta[Ntheta - 1, Ntheta - 1]
  var_sand_q0  = var_theta[Ntheta, Ntheta]
  cov_sand_q10 = var_theta[Ntheta, Ntheta - 1]
  var_sand_qte = var_sand_q1 + var_sand_q0 - 2 * cov_sand_q10
  
  c( q1_hat = q1_hat, q0_hat = q0_hat, qte_hat = qte_hat,
     var_sand_q1  = var_sand_q1, var_sand_q0  = var_sand_q0,
     cov_sand_q10 = cov_sand_q10, var_sand_qte = var_sand_qte )
  
}


run_sandwich_derivation_on_datalist = function(data_list,
                                                scenario = c("cc", "cm", "mc", "mm"),
                                                p_star = 0.5,
                                                n_run = NULL) {
  
  scenario = match.arg(scenario)
  
  res_list = vector("list", n_run)
  for (i in seq_len(n_run)) {
    cat("Running dataset", i, "of", n_run, "\n")
    
    out = tryCatch(
      {
        tmp = one_run_sandwich_derivation_onedata(
          data = data_list[[i]],
          scenario = scenario,
          p_star  = p_star
        )
        as.data.frame(as.list(tmp))
      },
      error = function(e) {
        cat("Dataset", i, "ERROR\n")
        data.frame(
          q1_hat       = NA_real_,
          q0_hat       = NA_real_,
          qte_hat      = NA_real_,
          var_sand_q1  = NA_real_,
          var_sand_q0  = NA_real_,
          cov_sand_q10 = NA_real_,
          var_sand_qte = NA_real_,
          stringsAsFactors = FALSE
        )
      }
    )
    
    
    out$se_sand_q1  = sqrt(out$var_sand_q1)
    out$se_sand_q0  = sqrt(out$var_sand_q0)
    out$se_sand_qte = sqrt(out$var_sand_qte)
    
    out$ci_q1_l  = out$q1_hat  - 1.96 * out$se_sand_q1
    out$ci_q1_u  = out$q1_hat  + 1.96 * out$se_sand_q1
    
    out$ci_q0_l  = out$q0_hat  - 1.96 * out$se_sand_q0
    out$ci_q0_u  = out$q0_hat  + 1.96 * out$se_sand_q0
    
    out$ci_qte_l = out$qte_hat - 1.96 * out$se_sand_qte
    out$ci_qte_u = out$qte_hat + 1.96 * out$se_sand_qte
    
    res_list[[i]] = out
  }
  
  bind_rows(res_list) %>%
    dplyr::select(
      q1_hat, q0_hat, qte_hat,
      var_sand_q1, var_sand_q0, cov_sand_q10, var_sand_qte,
      se_sand_q1, se_sand_q0, se_sand_qte,
      ci_q1_l, ci_q1_u, ci_q0_l, ci_q0_u, ci_qte_l, ci_qte_u
    )
}