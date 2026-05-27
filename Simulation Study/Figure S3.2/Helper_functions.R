#####################################
# data generating
######################################

generate_data = function(n, alpha1 = 0.5, alpha2 = 0.35,
                         beta1 = -2, beta2 = 3, delta = 2) {
  Y_all  = c(); A_all  = c(); X1_all = c(); X2_all = c()
  
  while(length(Y_all) < n) {
    
    batch_size = max(n - length(Y_all), 1)
    X1 = rbinom(batch_size, 1, prob = 0.5)
    X2 = rnorm(batch_size, mean = 0, sd = 1)
    ps_true = 1 / (1 + exp(-(alpha1 * X1 + alpha2 * X2)))
    A = rbinom(batch_size, 1, prob = ps_true)
    
    # Compute potential outcomes Y0_star and Y1_star
    Y0_star = beta1 * X1 + beta2 * X2 + rnorm(batch_size, 0, 1)
    Y1_star = beta1 * X1 + beta2 * X2 + delta + rnorm(batch_size, 0, 1)
    
    # Calculate p-values for Y0_star and Y1_star
    p0 = pnorm(Y0_star); p1 = pnorm(Y1_star)
    valid_idx = which((p0 > 1e-15 & p0 < 1 - 1e-15) & (p1 > 1e-15 & p1 < 1 - 1e-15))
    
    if (length(valid_idx) > 0) {
      # Apply chi-square transformation using qchisq for the valid indices
      Y0 = qchisq(p0[valid_idx], df = 5)
      Y1 = qchisq(p1[valid_idx], df = 5)
      Y  = ifelse(A[valid_idx] == 1, Y1, Y0)
      
      
      Y_all  = c(Y_all, Y)
      A_all  = c(A_all, A[valid_idx])
      X1_all = c(X1_all, X1[valid_idx])
      X2_all = c(X2_all, X2[valid_idx])
    }
  }
  
  data.frame(Y = Y_all[1:n], A = A_all[1:n], X1 = X1_all[1:n], X2 = X2_all[1:n])
}



###################################################
#Perform quantile inversion 
#####################################################
quantile_inversion = function(Y_sorted, CDF_sorted, p_target) {
  ncat = length(Y_sorted)
  
  weight_yuqi = c( 0, (CDF_sorted[-ncat] - CDF_sorted[1]) / (CDF_sorted[ncat - 1] - CDF_sorted[1]), 1 )
  
  weight = 1 - weight_yuqi
  #weight = weight_yuqi
  
  weighted_quantile = (1 - weight) * c(Y_sorted[1], Y_sorted) +
    weight * c(Y_sorted, Y_sorted[ncat])
  
  x_vals = c(0, CDF_sorted)
  y_vals = weighted_quantile
  
  
  keep = !duplicated(x_vals, fromLast = TRUE)
  x_vals = x_vals[keep]
  y_vals = y_vals[keep]
  
  res = approx(
    x = x_vals,
    y = y_vals,
    xout = p_target,
    rule = 2
  )
  
  return(res$y)
}

########################################################################
#########################################################################
cond_cdf_onemodel = function(data, y_star, g, alpha_hat, A = 1, tolerance = 1e-8,
                             link = c("logistic", "probit")) {
  link = match.arg(link)
  n = nrow(data); y_unique = g$yunique
  
  if (y_star <= min(y_unique) + tolerance) {
    g_y = rep(0, n)
  } else if (y_star >= max(y_unique) - tolerance) {
    g_y = rep(1, n)
  } else {
    idx_star = max(which(y_unique <= y_star + tolerance))
    myalpha = alpha_hat[idx_star]
    
    design_matrix = g$x ; design_matrix[, "A"] = A
    beta_hat = g$coefficients[colnames(design_matrix)]
    eta_hat = as.matrix(design_matrix) %*% beta_hat
    lp = myalpha - eta_hat
    if (link == "logistic") { g_y = plogis(lp) } else if (link == "probit") { g_y = pnorm(lp)}
  }
  
  return(list(g_y = as.vector(g_y)))
}

#####################################
# IPW_cpm_icdf
######################################
ipw_cpm_cdf = function(data, p_star, formA) {
  
  ps_fit = glm(formA, data = data, family = binomial("logit"))
  pi1 = predict(ps_fit, type = "response"); pi0 = 1 - pi1
  IA1 = as.numeric(data$A == 1); IA0 = as.numeric(data$A == 0)
  
  n = nrow(data); F1_ipw = F0_ipw = numeric(n)
  for (i in 1:n) {
    y_star = data$Y[i]
    IY = as.numeric(data$Y <= y_star)
    F1_ipw[i] = mean(IA1 * IY / pi1) 
    F0_ipw[i] = mean(IA0 * IY / pi0)
  }
  
  ord = order(data$Y);  Ygrid = data$Y[ord]
  F1_ipw_sorted = F1_ipw[ord]; F0_ipw_sorted = F0_ipw[ord]
  F1_ipw_sorted = cummax(pmin(pmax(F1_ipw_sorted, 0), 1))
  F0_ipw_sorted = cummax(pmin(pmax(F0_ipw_sorted, 0), 1))
  
  # all(diff(F1_ipw_sorted) >= 0) && all(F1_ipw_sorted >= 0 & F1_ipw_sorted <= 1)
  # all(diff(F0_ipw_sorted) >= 0) && all(F0_ipw_sorted >= 0 & F0_ipw_sorted <= 1)

  q1 = quantile_inversion(Ygrid, F1_ipw_sorted, p_star)
  q0 = quantile_inversion(Ygrid, F0_ipw_sorted, p_star)
  
  c(q1 = q1, q0 = q0, QTE = q1 - q0)
  
}

###########################################
# G-computation approach
###########################################
gcomp_cpm_cdf = function(data, p_star, formY, link = c("logistic", "probit")) {
  link = match.arg(link)
  #OR model
  g = orm(formY, family = link, x=T, y = T,data = data)
  alpha_hat = -g$coefficients[1:(length(g$yunique) - 1)]
  
  # OR 
  n = nrow(data); OR1  = numeric(n);OR0  = numeric(n)
  
  for (i in seq_len(n)) {
    
    y_star = data$Y[i]
    
    # A = 1
    g1_vals= cond_cdf_onemodel(data = data,y_star = y_star, g = g, alpha_hat = alpha_hat,A =1, link = "probit")
    OR1[i] = mean( g1_vals$g_y, na.rm = T)
    
    # A = 0
    g0_vals= cond_cdf_onemodel(data = data,y_star = y_star, g = g, alpha_hat = alpha_hat,A =0, link = "probit")
    OR0[i] = mean( g0_vals$g_y, na.rm = T)
  }
  
  ord = order(data$Y);  Ygrid = data$Y[ord]
  OR1_sorted = OR1[ord]; OR0_sorted = OR0[ord]
  OR1_sorted = cummax(pmin(pmax(OR1_sorted, 0), 1))
  OR0_sorted = cummax(pmin(pmax(OR0_sorted, 0), 1))
  
  # all(diff(OR1_sorted) >= 0) && all(OR1_sorted >= 0 & OR1_sorted <= 1)
  # all(diff(OR0_sorted) >= 0) && all(OR0_sorted >= 0 & OR0_sorted <= 1)
  
  q1 = quantile_inversion(Ygrid, OR1_sorted, p_star)
  q0 = quantile_inversion(Ygrid, OR0_sorted, p_star)
  
  c(q1 = q1, q0 = q0, QTE = q1 - q0)
}
#################################################
# AIPW_cpm_icdf
#################################################
aipw_cpm_cdf = function(data, p_star, formA, formY, link = c("probit", "logistic")) {
  
  link = match.arg(link)
  # PS model
  ps_fit  = glm(formA, data = data, family = binomial("logit"))
  pi1 = predict(ps_fit, type = "response");  pi0  = 1 - pi1                            
  IA1 = as.numeric(data$A == 1); IA0 = as.numeric(data$A == 0)
  weight1 = IA1 / pi1 ; weight0 = IA0 / pi0 
  
  #OR model
  g = orm(formY, family = link, x=T, y = T,data = data)
  alpha_hat = -g$coefficients[1:(length(g$yunique) - 1)]
  
  # DR 
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
  
  ord = order(data$Y); Ygrid1 = Ygrid0 = Ygrid = data$Y[ord]
  DR1_sorted = DR1[ord]; DR0_sorted = DR0[ord]
  DR1_sorted = cummax(pmin(pmax(DR1_sorted, 0), 1))
  DR0_sorted = cummax(pmin(pmax(DR0_sorted, 0), 1))
  
  # all(diff(DR1_sorted) >= 0) && all(DR1_sorted >= 0 & DR1_sorted <= 1)
  # all(diff(DR0_sorted) >= 0) && all(DR0_sorted >= 0 & DR0_sorted <= 1)
  
  q1 = quantile_inversion(Ygrid, DR1_sorted, p_star)
  q0 = quantile_inversion(Ygrid, DR0_sorted, p_star)
  
  c(q1 = q1, q0 = q0, QTE = q1 - q0)
}


######################################
# From Ivan 2017
#######################################
trim = function(x) pmax(x, 1e-10)

compute.quantile = function(Q, w, q, r) {
  F = function(y)sapply(y, function(x)mean(rowSums((Q <= x) * w)))
  inv = function(qq){
    froot = function(x) (F(x) - qq)
    uniroot(froot, r, extendInt = 'yes')$root
  }
  return(sapply(q, function(qq)inv(qq)))
}


#TMLE
tmle = function(y, t, Q, g, q){
  n = length(y)
  D  = function(y, w, chiq){ 1 / g * ((y <= chiq) - rowSums((Q <= chiq) * w)) }
  w  = matrix(1/dim(Q)[2], ncol = dim(Q)[2], nrow = n)
  h = t
  chiq = compute.quantile(Q, w, q, range(y))
  Do = D(y, w, chiq)
  Dq = D(Q, w, chiq)
  iter     = 1
  crit     = TRUE
  max.iter = 20
  while(crit && iter <= max.iter){
    est.eq = function(eps){
      out = - mean(h * (Do - rowSums(Dq * exp(eps * Dq) * w) /
                           rowSums(exp(eps * Dq) * w)))
      return(out)
    }
    loglik = function(eps){
      out = - mean(h * (eps * Do - log(rowSums(exp(eps * Dq) * w))))
      return(out)
    }
    eps = optim(par = 0, loglik, gr = est.eq, method = 'L-BFGS-B',  lower=-1e2, upper=20)$par
    w = exp(eps * Dq) * w / rowSums(exp(eps * Dq) * w)
    chiq = compute.quantile(Q, w, q, range(y))
    Do = D(y, w, chiq)
    Dq = D(Q, w, chiq)
    iter = iter + 1
    crit = abs(eps)  > 1e-4 / n^0.6
  }
  return(chiq)
}

# Firpo estimator
firpo = function(y, t, Q, g, q){
  library(quantreg)
  h = t / g
  chiq = coef(rq(y ~ 1, weights = h, tau = q))
  names(chiq) = NULL
  return(chiq)
}

#AIPW
aipw_normal_q = function(y, t, Q, g, q){
  n = length(y)
  w = matrix(1/dim(Q)[2], ncol = dim(Q)[2], nrow = n)
  h = t / g
  D  = Vectorize(function(chiq){
    mean(h * ((y <= chiq) - rowSums((Q <= chiq) * w))) +
      mean(rowSums((Q <= chiq) * w) - q)
  })
  chiq = uniroot(D, c(-1000, 1000), extendInt = 'yes')$root
  return(chiq)
}

#####################################################
# CPM based aipw with estimating q directly first
####################################################
aipw_cpm_q = function(data, p_star, formA, formY, link = c("probit", "logistic")) {
  link = match.arg(link)
  n = nrow(data); ord = order(data$Y); Ygrid = data$Y[ord]
  #PS model
  ps_fit  = glm(formA, data = data, family = binomial("logit"))
  pi1 = predict(ps_fit, type = "response"); pi0 = 1 - pi1       
  IA1 = as.numeric(data$A == 1); IA0 = as.numeric(data$A == 0)
  weight1 = IA1 / pi1; weight0 = IA0 / pi0
  
  #OR model
  g = orm(formY, family = link, x=T, y = T,data = data)
  alpha_hat = -g$coefficients[1:(length(g$yunique) - 1)]
  OR1  = matrix(NA, n, n);OR0  = matrix(NA, n, n)
  for (i in 1:n) {
    y_star = data$Y[i]
    g1_vals= cond_cdf_onemodel(data = data,y_star = y_star, g = g, alpha_hat = alpha_hat,A =1, link = link )
    g0_vals= cond_cdf_onemodel(data = data,y_star = y_star, g = g, alpha_hat = alpha_hat,A =0, link = link )
    OR1[,i] = g1_vals$g_y;OR0[,i] = g0_vals$g_y
  }
  c_cdf1 = OR1[,ord] ;  c_cdf0 = OR0[,ord]   
  

  F_DR1 = function(q, A = 1) {
    if (q <= Ygrid[1])        return(0)
    if (q >= Ygrid[n])        return(1)
    IY1 = as.numeric(data$Y <= q)
    j= findInterval(q, Ygrid, rightmost.closed = TRUE)
    DR1 = mean( (IA1 / pi1) * (IY1 - c_cdf1[, j] ) + c_cdf1[, j] )
    DR1
  }
  F_DR0 = function(q, A = 0) {
    if (q <= Ygrid[1])        return(0)
    if (q >= Ygrid[n])        return(1)
    IY0 = as.numeric(data$Y <= q)
    j= findInterval(q, Ygrid, rightmost.closed = TRUE)
    DR0 = mean( (IA0 / pi0 )* (IY0 - c_cdf0[, j] ) + c_cdf0[, j] )
    DR0
  }
  
  lower = min(Ygrid, na.rm = TRUE)
  upper = max(Ygrid, na.rm = TRUE)
  
  f1 = function(q) F_DR1(q, 1) - p_star
  f0 = function(q) F_DR0(q, 0) - p_star
  
  fl1 = f1(lower)
  fu1 = f1(upper)
  fl0 = f0(lower)
  fu0 = f0(upper)
  
  if (!is.finite(fl1) || !is.finite(fu1)) {
    stop("F_DR1 boundary is not finite")
  }
  if (!is.finite(fl0) || !is.finite(fu0)) {
    stop("F_DR0 boundary is not finite")
  }
  if (fl1 * fu1 > 0) {
    stop("F_DR1 has no sign change in interval")
  }
  if (fl0 * fu0 > 0) {
    stop("F_DR0 has no sign change in interval")
  }
  
  q1 = uniroot(f1, lower = lower, upper = upper, tol = 1e-8)$root
  q0 = uniroot(f0, lower = lower, upper = upper, tol = 1e-8)$root
  QTE = q1 - q0
  
  c(q1 = q1, q0 = q0, QTE = QTE)
}

####################################
# Compute Bias, Variance, and MSE
####################################
calc_stats = function(est_vec, true_val, nsample) {
  estimate = mean(est_vec, na.rm =T)
  bias = mean(est_vec - true_val,na.rm =T)
  variance = var(est_vec,na.rm = T)
  Nvar = variance * nsample
  mse = bias^2 + variance
  rmse = sqrt(mse)
  bias_sd = bias/sqrt(variance)
  mae = mean(abs(est_vec - true_val), na.rm = TRUE)
  medae = median(abs(est_vec - true_val), na.rm = TRUE)
  c(Estimate =estimate, Bias = bias, Variance = variance, MSE = mse, RMSE = rmse, MAE = mae,MedAE =  medae, bias_sd = bias_sd,Nvar = Nvar)
}
