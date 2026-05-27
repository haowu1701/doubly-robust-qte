single_sim = function(n = 1000, p_star = 0.5,
                      alpha1 = 0.5, alpha2 = 0.35,
                      beta1 = -2, beta2 = 3, delta = 2,
                      scenario = "cc") {
  options(warn = -1)
  
  # Pre-specified models (PS & outcome) based on scenario
  # scenarios: "cc", "cm", "mc", "mm"
  if (scenario == "cc") {
    formA   <- A ~ X1 + X2            # correct PS
    formY   <- Y ~ X1 + X2 + A        # correct outcome
    rq_form <- Y ~ X1 + X2
  } else if (scenario == "cm") {
    formA   <- A ~ X1                 # mis-specified PS
    formY   <- Y ~ X1 + X2 + A        # correct outcome
    rq_form <- Y ~ X1 + X2
  } else if (scenario == "mc") {
    formA   <- A ~ X1 + X2            # correct PS
    formY   <- Y ~ A + X1             # mis-specified outcome
    rq_form <- Y ~ X1
  } else if (scenario == "mm") {
    formA   <- A ~ X1                 # mis-specified PS
    formY   <- Y ~ A + X1             # mis-specified outcome
    rq_form <- Y ~ X1
  } else {
    stop("scenario must be one of: cc, cm, mc, mm")
  }
  
  # Generate data
  data <- generate_data(n, alpha1, alpha2, beta1, beta2, delta)

  # Initialize timing log
  time_log <- list()
  
  #########################
  # 1. CPM-based CDF estimators
  #########################
  time_log$IPW <- system.time({
    result_IPW  <- ipw_cpm_cdf(data, p_star, formA = formA)
  })["elapsed"]
  
  time_log$Gcomp <- system.time({
    result_Gcomp <- gcomp_cpm_cdf(data, p_star, formY = formY)
  })["elapsed"]
  
  time_log$AIPW <- system.time({
    result_AIPW  <- aipw_cpm_cdf(data, p_star, formA = formA, formY = formY)
  })["elapsed"]
  
  
  ###################################
  # 2. Ivan et al. (2017) – estimate Q first
  ###################################
 
    n.quant <- 300
    Y       <- data$Y
    A       <- data$A
    X1      <- data$X1
    X2      <- data$X2
    covariates <- data.frame(X1 = X1, X2 = X2)
    
    fitT   <- glm(formA, data = data, family = binomial)
    fitY1  <- lm(formY, data = data, subset = A == 1)
    fitY0  <- lm(formY, data = data, subset = A == 0)
    
    median1 <- predict(fitY1, newdata = data.frame(A = 1, covariates))
    median0 <- predict(fitY0, newdata = data.frame(A = 0, covariates))
    
    Q1 <- sapply(seq(1/n.quant, 1 - 1/n.quant, 1/n.quant),
                 function(q) qnorm(q, mean = median1, sd = summary(fitY1)$sigma))
    Q0 <- sapply(seq(1/n.quant, 1 - 1/n.quant, 1/n.quant),
                 function(q) qnorm(q, mean = median0, sd = summary(fitY0)$sigma))
    
    g1 <- trim(predict(fitT, type = 'response'))

  
  
  #######
  # 3. TMLE on those Q’s
  #######
  time_log$TMLE <- system.time({
    q1_tmle  <- tmle(Y, A, Q1, g1, p_star)
    q0_tmle  <- tmle(Y, 1 - A, Q0, 1 - g1, p_star)
    QTE_tmle <- q1_tmle - q0_tmle
  })["elapsed"]
  
  
  ####
  # 4. AIPW-normal quantile
  ####
  time_log$AIPW_normal <- system.time({
    q1_aipw  <- aipw_normal_q(Y, A, Q1, g1, p_star)
    q0_aipw  <- aipw_normal_q(Y, 1 - A, Q0, 1 - g1, p_star)
    QTE_aipw <- q1_aipw - q0_aipw
  })["elapsed"]
  
  
  ####
  # 5. Firpo
  ####
  time_log$Firpo <- system.time({
    q1_firpo  <- firpo(Y, A, Q1, g1, p_star)
    q0_firpo  <- firpo(Y, 1 - A, Q0, 1 - g1, p_star)
    QTE_firpo <- q1_firpo - q0_firpo
  })["elapsed"]
  
  
  
  
  ###################################
  # 7. CPM-based TMLE
  ###################################
  
  fitT  <- glm(formA, data = data, family = binomial)
  g1    <- trim(predict(fitT, type = 'response'))
  
  p_seq <- seq(1/300, 1 - 1/300, 1/300)
  n_p <- length(p_seq)
  ord <- order(data$Y)
  Y_sorted <- data$Y[ord]
  C1 <- C0 <- matrix(NA, nrow = n, ncol = n)
  cpm_Q1 <- cpm_Q0 <- matrix(NA, nrow = n, ncol = n_p)

  g_tmle_cpm <- orm(formY, family = "probit", x = TRUE, y = TRUE, data = data)
  alpha_hat <- -g_tmle_cpm$coefficients[1:(length(g_tmle_cpm$yunique) - 1)]

  for (i in seq_len(n)) {
    y_star <- data$Y[i]

    g1_vals <- cond_cdf_onemodel(
      data = data, y_star = y_star, g = g_tmle_cpm,
      alpha_hat = alpha_hat, A = 1, link = "probit"
    )
    C1[, i] <- g1_vals$g_y

    g0_vals <- cond_cdf_onemodel(
      data = data, y_star = y_star, g = g_tmle_cpm,
      alpha_hat = alpha_hat, A = 0, link = "probit"
    )
    C0[, i] <- g0_vals$g_y
  }

  for (i in seq_len(n)) {
    cdf1_i <- C1[i, ]
    cdf1_i <- cdf1_i[ord]
    cdf0_i <- C0[i, ]
    cdf0_i <- cdf0_i[ord]

    cpm_Q1[i, ] <- quantile_inversion(
      Y_sorted = Y_sorted,
      CDF_sorted = cdf1_i,
      p_target = p_seq
    )
    cpm_Q0[i, ] <- quantile_inversion(
      Y_sorted = Y_sorted,
      CDF_sorted = cdf0_i,
      p_target = p_seq
    )
  }
  
  time_log$CPM_TMLE <- system.time({
    q1_cpm_tmle  <- tmle(Y, A, cpm_Q1, g1, p_star)
    q0_cpm_tmle  <- tmle(Y, 1 - A, cpm_Q0, 1 - g1, p_star)
    QTE_cpm_tmle <- q1_cpm_tmle - q0_cpm_tmle
  })["elapsed"]
  
  
  ###################################
  # 8. RQ-based TMLE
  ###################################
  taus <- seq(1/300, 299/300, length.out = 299)
  
  fitT  <- glm(formA, data = data, family = binomial)
  g1    <- trim(predict(fitT, type = 'response'))
  
  fitY1 <- rq(rq_form, data = data, subset = A == 1, tau = taus)
  fitY0 <- rq(rq_form, data = data, subset = A == 0, tau = taus)
  
  cqf1 <- predict(fitY1, type = "Qhat", newdata = data, stepfun = TRUE)
  cqf1 <- lapply(cqf1, rearrange)
  rq_Q1 <- t(sapply(cqf1, function(f) f(unique(knots(f)))))
  
  cqf0 <- predict(fitY0, type = "Qhat", newdata = data.frame(Y = Y, A = A, covariates), stepfun = TRUE)
  cqf0 <- lapply(cqf0, rearrange)
  rq_Q0 <- t(sapply(cqf0, function(f) f(unique(knots(f)))))
  
  time_log$RQ_TMLE <- system.time({
    q1_rq_tmle  <- tmle(Y, A, rq_Q1, g1, p_star)
    q0_rq_tmle  <- tmle(Y, 1 - A, rq_Q0, 1 - g1, p_star)
    QTE_rq_tmle <- q1_rq_tmle - q0_rq_tmle
  })["elapsed"]
  
  
  ####################################################
  # 9. CPM-based AIPW_q
  ####################################################
  time_log$AIPW_CPM_Q <- system.time({
    result_aipw_cpm_q <- aipw_cpm_q(data, p_star, formA = formA, formY = formY)
  })["elapsed"]
  
  
  # Combine all estimates
  result <- c(
    aipw_q_q1       = q1_aipw,
    aipw_q_q0       = q0_aipw,
    aipw_q_QTE      = QTE_aipw,
    
    tmle_q1         = q1_tmle,
    tmle_q0         = q0_tmle,
    tmle_QTE        = QTE_tmle,
    
    firpo_q1        = q1_firpo,
    firpo_q0        = q0_firpo,
    firpo_QTE       = QTE_firpo,
    

    
    rq_tmle_q1      = q1_rq_tmle,
    rq_tmle_q0      = q0_rq_tmle,
    rq_tmle_QTE     = QTE_rq_tmle,
    
    cpm_tmle_q1     = q1_cpm_tmle,
    cpm_tmle_q0     = q0_cpm_tmle,
    cpm_tmle_QTE    = QTE_cpm_tmle,
    
    aipw_cpm_q_q1   = result_aipw_cpm_q["q1"],
    aipw_cpm_q_q0   = result_aipw_cpm_q["q0"],
    aipw_cpm_q_QTE  = result_aipw_cpm_q["QTE"],
    
    AIPW_q1         = result_AIPW["q1"],
    AIPW_q0         = result_AIPW["q0"],
    AIPW_QTE        = result_AIPW["QTE"],
    
    Gcomp_q1        = result_Gcomp["q1"],
    Gcomp_q0        = result_Gcomp["q0"],
    Gcomp_QTE       = result_Gcomp["QTE"],
    
    IPW_q1          = result_IPW["q1"],
    IPW_q0          = result_IPW["q0"],
    IPW_QTE         = result_IPW["QTE"]
  )
  
  # Append timing info
  times <- unlist(time_log)
  names(times) <- paste0("time_", names(times))
  
  return(c(result, times))
}
