single_sim= function(n= 1000, p_star= 0.1,
                       alpha1= 0.5, alpha2= 0.35,
                       beta1= -2, beta2= 3, delta= 2,
                       scenario= "cc") {
  options(warn= -1)
  source("Helper_functions.R")
  
  # scenarios { "cc", "cm", "mc", "mm" }
  if (scenario== "cc") {
    formA  = A ~ X1 + X2             # correct PS
    formY  = Y ~ X1 + X2 + A         # correct outcome
    rq_form = Y ~ X1 + X2
    logformA = A ~ X1 + X2           
    logformY = logY ~ X1 + X2 + A    
    logrq_form = logY ~ X1 + X2
    cformA= A ~ X1 + X2 
    cformY= cY ~ X1 + X2 +A
    crq_form= cY ~ X1 + X2        
  } else if (scenario== "cm") {
    formA  = A ~ X1                  # mis-specified PS
    formY  = Y ~ X1 + X2 + A         # correct outcome
    rq_form = Y ~ X1 + X2
    logformA = A ~ X1                
    logformY = logY ~ X1 + X2 + A    
    logrq_form = logY ~ X1 + X2
    cformA= A ~ X1    
    cformY= cY ~ X1 + X2 +A
    crq_form= cY ~ X1 + X2      
  } else if (scenario== "mc") {
    formA  = A ~ X1 + X2             # correct PS
    formY  = Y ~ A + X1              # mis-specified outcome
    rq_form = Y ~ X1 
    logformA = A ~ X1 + X2          
    logformY = logY ~ A + X1         
    logrq_form = logY ~ X1
    cformA= A ~ X1 + X2 
    cformY= cY ~ A + X1
    crq_form = cY ~ X1            
  } else if (scenario== "mm") {
    formA  = A ~ X1                 # mis-specified PS
    formY  = Y ~ A + X1             # mis-specified outcome
    rq_form = Y ~ X1 
    logformA = A ~ X1             
    logformY = logY ~ A + X1        
    logrq_form = logY ~ X1
    cformA= A ~ X1  
    cformY= cY ~  A + X1
    crq_form= cY ~ X1
  } else {
    stop("scenario must be one of: cc, cm, mc, mm")
  }
  
  # Generate data
  data= generate_data(n, alpha1, alpha2, beta1, beta2, delta)
 
  #########################
  # CPM based estimators, correct link function--probit
  #########################
  # inverse distribution
  result_IPW= ipw_cpm_cdf(data, p_star, formA= formA)
  result_Gcomp= gcomp_cpm_cdf(data, p_star, formY= formY, link = "probit")
  result_AIPW= aipw_cpm_cdf(data,  p_star, formA= formA, formY= formY, link = "probit")
  #directly solving q
  result_aipw_cpm_q = aipw_cpm_q(data,  p_star, formA= formA, formY= formY, link = "probit")
  # TMLE-CPM
  n.quant = 300
  Y = data$Y; A = data$A
  fitT = glm(formA, data = data, family = binomial)
  g1 = trim(predict(fitT, type = 'response'))
  p_seq = seq(1/n.quant, 1 - 1/n.quant, 1/n.quant) ;  n_p = length(p_seq)
  
  ord = order(data$Y);  Y_sorted = data$Y[ord]
  C1 = C0 = matrix(NA, nrow = n, ncol = n)
  cpm_Q1 = cpm_Q0 = matrix(NA, nrow = n, ncol = n_p) 
  #OR model
  g_tmle_cpm = orm(formY, family = "probit", x=T, y = T,data = data)
  alpha_hat = -g_tmle_cpm$coefficients[1:(length( g_tmle_cpm$yunique) - 1)]
  
  for (i in seq_len(n)) {
    y_star = data$Y[i]
    # A = 1
    g1_vals= cond_cdf_onemodel(data = data,y_star = y_star, g = g_tmle_cpm , alpha_hat = alpha_hat,A =1, link = "probit")
    C1[,i] = g1_vals$g_y
    # A = 0
    g0_vals= cond_cdf_onemodel(data = data,y_star = y_star, g = g_tmle_cpm , alpha_hat = alpha_hat,A =0, link = "probit")
    C0[,i] = g0_vals$g_y
  }
  
  for (i in seq_len(n)) {
    cdf1_i = C1[i,]; cdf1_i = cdf1_i[ord]
    cdf0_i = C0[i,]; cdf0_i = cdf0_i[ord]
    cpm_Q1[i, ] = quantile_inversion(Y_sorted = Y_sorted,CDF_sorted = cdf1_i,p_target = p_seq )
    cpm_Q0[i, ] = quantile_inversion(Y_sorted = Y_sorted,CDF_sorted = cdf0_i,p_target = p_seq)
  }
  
  q1_cpm_tmle = tmle(Y, A, cpm_Q1, g1, p_star)
  q0_cpm_tmle =  tmle(Y, 1 - A, cpm_Q0, 1 - g1, p_star)
  QTE_cpm_tmle = q1_cpm_tmle - q0_cpm_tmle
  
  
  #########################
  # CPM based estimators, mislink function--logit
  #########################
  # inverse distribution
  result_AIPW_mis= aipw_cpm_cdf(data,  p_star, formA= formA, formY= formY, link = "logistic")
  #directly solving q
  result_aipw_cpm_q_mis = aipw_cpm_q(data,  p_star, formA= formA, formY= formY, link = "logistic")
  # TMLE-CPM
  C1_mis = C0_mis = matrix(NA, nrow = n, ncol = n)
  cpm_Q1_mis = cpm_Q0_mis = matrix(NA, nrow = n, ncol = n_p) 
  #OR model
  g_tmle_cpm_mis = orm(formY, family = "logistic", x=T, y = T,data = data)
  alpha_hat_mis = - g_tmle_cpm_mis$coefficients[1:(length( g_tmle_cpm_mis$yunique) - 1)]
  
  for (i in seq_len(n)) {
    y_star = data$Y[i]
    # A = 1
    g1_vals= cond_cdf_onemodel(data = data,y_star = y_star, g =  g_tmle_cpm_mis , alpha_hat = alpha_hat_mis,A =1, link = "logistic")
    C1_mis[,i] = g1_vals$g_y
    # A = 0
    g0_vals= cond_cdf_onemodel(data = data,y_star = y_star, g =  g_tmle_cpm_mis, alpha_hat = alpha_hat_mis,A =0, link = "logistic")
    C0_mis[,i] = g0_vals$g_y
  }
  
  for (i in seq_len(n)) {
    cdf1_i_mis = C1_mis[i,]; cdf1_i_mis = cdf1_i_mis[ord]
    cdf0_i_mis = C0_mis[i,]; cdf0_i_mis = cdf0_i_mis[ord]
    
    cpm_Q1_mis[i, ] = quantile_inversion(Y_sorted = Y_sorted,CDF_sorted = cdf1_i_mis,p_target = p_seq)
    cpm_Q0_mis[i, ] = quantile_inversion(Y_sorted = Y_sorted,CDF_sorted = cdf0_i_mis,p_target = p_seq)
  }
  
  q1_cpm_tmle_mis = tmle(Y, A, cpm_Q1_mis, g1, p_star)
  q0_cpm_tmle_mis =  tmle(Y, 1 - A, cpm_Q0_mis, 1 - g1, p_star)
  QTE_cpm_tmle_mis = q1_cpm_tmle_mis - q0_cpm_tmle_mis
  
  ###################################
  #Ivan 2017 estimate q directly 
  ####################################
  n.quant = 300
  Y = data$Y; A = data$A
  
  X1 = data$X1 ;  X2 = data$X2;  covariates = as.data.frame(cbind(X1,X2))
  
  fitT = glm(formA,data = data, family = binomial)
  fitY1 = lm(formY, data = data, subset = (data$A == 1))
  fitY0 = lm(formY, data = data, subset = (data$A == 0))
  
  median1 = predict(fitY1, newdata = data.frame(A = 1, covariates))
  median0 = predict(fitY0, newdata = data.frame(A = 0, covariates))
  Q1 = sapply(seq(1/n.quant, 1 - 1/n.quant, 1/n.quant),
              function(q)qnorm(q, mean = median1, sd = summary(fitY1)$sigma))
  Q0 = sapply(seq(1/n.quant, 1 - 1/n.quant, 1/n.quant),
              function(q)qnorm(q, mean = median0, sd = summary(fitY0)$sigma))
  g1 = trim(predict(fitT, type = 'response'))
  
  # TMLE
  q1_tmle = tmle(Y, A, Q1, g1, p_star)
  q0_tmle =  tmle(Y, 1 - A, Q0, 1 - g1, p_star)
  QTE_tmle = q1_tmle - q0_tmle
  
  #AIPW 
  q1_aipw = aipw_normal_q(Y, A, Q1, g1, p_star)
  q0_aipw = aipw_normal_q(Y, 1 - A, Q0, 1 - g1, p_star)
  QTE_aipw = q1_aipw - q0_aipw
  
  #IPW-Firpo-q
  q1_firpo = firpo(Y, A, Q1, g1, p_star)
  q0_firpo = firpo(Y, 1 - A, Q0, 1 - g1, p_star)
  QTE_firpo = q1_firpo - q0_firpo
  
  ###################################
  # Ivan 2017 estimate q directly -- take log transformation before analysis
  ####################################
  n.quant = 300
  logY = log(data$Y)
  X1 = data$X1 ;  X2 = data$X2;  covariates = as.data.frame(cbind(X1,X2))
  
  logfitT = glm(logformA, data = data, family = binomial)
  logfitY1 = lm(logformY, data = data.frame(logY=logY, covariates), subset = A == 1)
  logfitY0 = lm(logformY, data = data.frame(logY=logY, covariates), subset = A == 0)
  logmedian1 = predict(logfitY1, newdata = data.frame(A = 1, covariates))
  logmedian0 = predict(logfitY0, newdata = data.frame(A = 0, covariates))
  
  logQ1 = sapply(seq(1/n.quant, 1 - 1/n.quant, 1/n.quant),
              function(q)qnorm(q, mean = logmedian1, sd = summary(logfitY1)$sigma))
  logQ0 = sapply(seq(1/n.quant, 1 - 1/n.quant, 1/n.quant),
              function(q)qnorm(q, mean = logmedian0, sd = summary(logfitY0)$sigma))
  logg1 = trim(predict(logfitT, type = 'response'))
  
  # log-TMLE
  log_q1_tmle =  tmle(logY, A, logQ1, logg1, p_star) 
  log_q0_tmle =  tmle(logY, 1 - A, logQ0, 1 - logg1, p_star) 
  log_QTE_tmle =  log_q1_tmle - log_q0_tmle 
  
  exp_q1_tmle = exp(log_q1_tmle)
  exp_q0_tmle = exp(log_q0_tmle)
  exp_QTE_tmle = exp_q1_tmle - exp_q0_tmle

  #log- AIPW
  log_q1_aipw = aipw_normal_q(logY, A, logQ1, logg1, p_star)
  log_q0_aipw = aipw_normal_q(logY, 1 - A, logQ0, 1 - logg1, p_star)
  log_QTE_aipw = log_q1_aipw - log_q0_aipw
  
  exp_q1_aipw = exp(log_q1_aipw)
  exp_q0_aipw = exp(log_q0_aipw)
  exp_QTE_aipw = exp_q1_aipw - exp_q0_aipw
  
  
  ###################################
  # Ivan 2017 estimate q directly --take correct transformation before analysis
  ####################################
  n.quant = 300
  cY = qnorm(pchisq(data$Y, df = 5))
  X1 = data$X1 ;  X2 = data$X2;  covariates = as.data.frame(cbind(X1,X2))
  cfitT = glm(cformA, data = data, family = binomial)
  
  cfitY1 = lm(cformY, data = data.frame(cY=cY, covariates), subset = A == 1)
  cfitY0 = lm(cformY, data = data.frame(cY=cY, covariates), subset = A == 0)
  cmedian1 = predict(cfitY1, newdata = data.frame(A=1, covariates))
  cmedian0 = predict(cfitY0, newdata = data.frame(A=0, covariates))
  
  cQ1 = sapply(seq(1/n.quant, 1 - 1/n.quant, 1/n.quant),
                 function(q)qnorm(q, mean = cmedian1, sd = summary(cfitY1)$sigma))
  cQ0 = sapply(seq(1/n.quant, 1 - 1/n.quant, 1/n.quant),
                 function(q)qnorm(q, mean = cmedian0, sd = summary(cfitY0)$sigma))
  cg1 = trim(predict(cfitT, type = 'response'))
  
  # ct-TMLE
  c_q1_tmle =  tmle(cY, A, cQ1, cg1, p_star) 
  c_q0_tmle =  tmle(cY, 1 - A, cQ0, 1 - cg1, p_star) 
  c_QTE_tmle =  c_q1_tmle - c_q0_tmle 
  tc_q1_tmle = qchisq(pnorm(c_q1_tmle), df = 5)
  tc_q0_tmle = qchisq(pnorm(c_q0_tmle), df = 5)
  tc_QTE_tmle = tc_q1_tmle - tc_q0_tmle

  #ct-AIPW
  c_q1_aipw = aipw_normal_q(cY, A, cQ1, cg1, p_star)
  c_q0_aipw = aipw_normal_q(cY, 1 - A, cQ0, 1 - cg1, p_star)
  c_QTE_aipw = c_q1_aipw - c_q0_aipw
  tc_q1_aipw = qchisq(pnorm(c_q1_aipw), df = 5)
  tc_q0_aipw = qchisq(pnorm(c_q0_aipw), df = 5)
  tc_QTE_aipw = tc_q1_aipw - tc_q0_aipw
  
  ###################################
  # TMLE cqr (conditional quantile regression)
  ####################################
  taus = 1:299/300
  Y = data$Y; A = data$A
  X1 = data$X1 ;  X2 = data$X2;   covariates = as.data.frame(cbind(X1,X2))
  
  fitT = glm(formA, data = data, family = binomial)
  fitY1 = rq(rq_form, data = data, subset = A == 1, tau = taus)
  fitY0 = rq(rq_form, data = data, subset = A == 0,tau = taus)
  
  cqf1 = predict(fitY1,type = "Qhat", newdata = data, stepfun = TRUE)
  cqf1 = lapply(cqf1, rearrange)
  rq_Q1 = t( sapply(cqf1, function(f) f( unique(knots(f)) )) )
  
  cqf0 = predict(fitY0,type = "Qhat", newdata = data, stepfun = TRUE)
  cqf0 = lapply(cqf0, rearrange)
  rq_Q0 = t( sapply(cqf0, function(f) f( unique(knots(f)) )) )
  
  g1 = trim(predict(fitT, type = 'response'))
  
  q1_rq_tmle = tmle(Y, A, rq_Q1, g1, p_star)
  q0_rq_tmle =  tmle(Y, 1 - A, rq_Q0, 1 - g1, p_star)
  QTE_rq_tmle = q1_rq_tmle - q0_rq_tmle
  
  ###################################
  # TMLE cqr--- take log transformation before analysis
  ####################################
  taus = 1:299/300
  
  Y = data$Y; A = data$A
  logY = log(data$Y)
  X1 = data$X1 ;  X2 = data$X2;  covariates = as.data.frame(cbind(X1,X2))
  
  logfitT = glm(logformA, data = data, family = binomial)
  logrqfitY1 = rq(logrq_form, data = data, subset = A == 1, tau = taus)
  logrqfitY0 = rq(logrq_form, data = data, subset = A == 0,tau = taus)
  
  logcqf1 = predict(logrqfitY1,type = "Qhat", newdata = data, stepfun = TRUE)
  logcqf1 = lapply(logcqf1, rearrange)
  logrq_Q1 = t( sapply(logcqf1, function(f) f( unique(knots(f)) )) )
  
  logcqf0 = predict(logrqfitY0,type = "Qhat", newdata = data.frame(Y=Y, A=A, covariates), stepfun = TRUE)
  logcqf0 = lapply(logcqf0, rearrange)
  logrq_Q0 = t( sapply(logcqf0, function(f) f( unique(knots(f)) )) )
  
  logg1 = trim(predict(logfitT, type = 'response'))
  
  logq1_rq_tmle = tmle(logY, A, logrq_Q1, logg1, p_star)
  logq0_rq_tmle =  tmle(logY, 1 - A, logrq_Q0, 1 - logg1, p_star)
  logQTE_rq_tmle = logq1_rq_tmle - logq0_rq_tmle
  
  exp_q1_rq_tmle  = exp(logq1_rq_tmle)
  exp_q0_rq_tmle = exp(logq0_rq_tmle)
  exp_QTE_rq_tmle =  exp_q1_rq_tmle- exp_q0_rq_tmle
  
  
  ###################################
  # TMLE cqr---take correct transformation before analysis
  ####################################
  taus = 1:299/300
  Y = data$Y; A = data$A
  cY = qnorm(pchisq(data$Y, df = 5))
  X1 = data$X1 ;  X2 = data$X2;  covariates = as.data.frame(cbind(X1,X2))
  
  cfitT = glm(cformA, data = data, family = binomial)
  crqfitY1 = rq(crq_form, data = data, subset = A == 1, tau = taus)
  crqfitY0 = rq(crq_form, data = data, subset = A == 0,tau = taus)
  
  ccqf1 = predict(crqfitY1,type = "Qhat", newdata = data, stepfun = TRUE)
  ccqf1 = lapply(ccqf1, rearrange)
  crq_Q1 = t( sapply(ccqf1, function(f) f( unique(knots(f)) )) )
  
  ccqf0 = predict(crqfitY0,type = "Qhat", newdata = data, stepfun = TRUE)
  ccqf0 = lapply(ccqf0, rearrange)
  crq_Q0 = t( sapply(ccqf0, function(f) f( unique(knots(f)) )) )
  
  cg1 = trim(predict(cfitT, type = 'response'))
  
  cq1_rq_tmle = tmle(cY, A, crq_Q1, cg1, p_star)
  cq0_rq_tmle =  tmle(cY, 1 - A, crq_Q0, 1 - cg1, p_star)
  cQTE_rq_tmle = cq1_rq_tmle - cq0_rq_tmle
  
  tc_q1_rq_tmle  = qchisq(pnorm(cq1_rq_tmle), df = 5)
  tc_q0_rq_tmle = qchisq(pnorm(cq0_rq_tmle), df = 5)
  tc_QTE_rq_tmle =  tc_q1_rq_tmle- tc_q0_rq_tmle
  
  # Combine the results into a named vector
  result= c(
    
    aipw_q_q1 = q1_aipw,
    aipw_q_q0 = q0_aipw,
    aipw_q_QTE = QTE_aipw,
    
    tmle_q1 = q1_tmle,
    tmle_q0 = q0_tmle,
    tmle_QTE = QTE_tmle,
    
    rq_tmle_q1 = q1_rq_tmle,
    rq_tmle_q0 = q0_rq_tmle,
    rq_tmle_QTE = QTE_rq_tmle,
    
    cpm_tmle_q1 = q1_cpm_tmle,
    cpm_tmle_q0 = q0_cpm_tmle,
    cpm_tmle_QTE = QTE_cpm_tmle,
    
    aipw_cpm_q_q1 = as.numeric(result_aipw_cpm_q["q1"]),
    aipw_cpm_q_q0 = as.numeric(result_aipw_cpm_q["q0"]),
    aipw_cpm_q_QTE = as.numeric(result_aipw_cpm_q["QTE"]),
    
    aipw_q1= as.numeric(result_AIPW["q1"]),
    aipw_q0= as.numeric(result_AIPW["q0"]),
    aipw_QTE= as.numeric(result_AIPW["QTE"]),
    
    g_q1= as.numeric(result_Gcomp["q1"]),
    g_q0= as.numeric(result_Gcomp["q0"]),
    g_QTE= as.numeric(result_Gcomp["QTE"]),
    
    
    iptw_q1=as.numeric(result_IPW["q1"]),
    iptw_q0=as.numeric(result_IPW["q0"]),
    iptw_QTE= as.numeric(result_IPW["QTE"]),
    
    firpo_q1 = q1_firpo,
    firpo_q0 = q0_firpo,
    firpo_QTE = QTE_firpo,
    
    exp_tmle_lm_q1 = exp_q1_tmle,
    exp_tmle_lm_q0 = exp_q0_tmle,
    exp_tmle_lm_QTE = exp_QTE_tmle,
    
    exp_tmle_rq_q1  = exp_q1_rq_tmle ,
    exp_tmle_rq_q0 = exp_q0_rq_tmle,
    exp_tmle_rq_QTE = exp_QTE_rq_tmle,
    
    exp_aipw_lm_q1  = exp_q1_aipw ,
    exp_aipw_lm_q0 = exp_q0_aipw ,
    exp_aipw_lm_QTE = exp_QTE_aipw ,
    
    tc_tmle_lm_q1 = tc_q1_tmle,
    tc_tmle_lm_q0 = tc_q0_tmle,
    tc_tmle_lm_QTE = tc_QTE_tmle,
    
    tc_tmle_rq_q1  = tc_q1_rq_tmle ,
    tc_tmle_rq_q0 = tc_q0_rq_tmle,
    tc_tmle_rq_QTE = tc_QTE_rq_tmle,
    
    tc_aipw_lm_q1  = tc_q1_aipw ,
    tc_aipw_lm_q0 = tc_q0_aipw ,
    tc_aipw_lm_QTE = tc_QTE_aipw,
    
    aipw_q1_mis = as.numeric(result_AIPW_mis["q1"]),
    aipw_q0_mis = as.numeric(result_AIPW_mis["q0"]),
    aipw_QTE_mis = as.numeric(result_AIPW_mis["QTE"]),
    
    aipw_cpm_q_q1_mis = as.numeric(result_aipw_cpm_q_mis["q1"]),
    aipw_cpm_q_q0_mis = as.numeric(result_aipw_cpm_q_mis["q0"]),
    aipw_cpm_q_QTE_mis = as.numeric(result_aipw_cpm_q_mis["QTE"]),
    
    cpm_tmle_q1_mis = q1_cpm_tmle_mis,
    cpm_tmle_q0_mis = q0_cpm_tmle_mis,
    cpm_tmle_QTE_mis = QTE_cpm_tmle_mis
    
  )
  
  
  return(result)
}
