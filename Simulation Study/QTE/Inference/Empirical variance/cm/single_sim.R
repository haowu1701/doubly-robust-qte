single_sim= function(n= 1000, p_star= 0.5,
                       alpha1= 0.5, alpha2= 0.35,
                       beta1= -2, beta2= 3, delta= 2,
                       scenario= "cc") {
  options(warn= -1)
  
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
  result_AIPW= aipw_cpm_cdf(data,  p_star, formA= formA, formY= formY, link = "probit")
  #directly solving q
  result_aipw_cpm_q = aipw_cpm_q(data,  p_star, formA= formA, formY= formY, link = "probit")
  
  
  # Combine the results into a named vector
  result= c(
    
    aipw_cpm_q_q1 = as.numeric(result_aipw_cpm_q["q1"]),
    aipw_cpm_q_q0 = as.numeric(result_aipw_cpm_q["q0"]),
    aipw_cpm_q_QTE = as.numeric(result_aipw_cpm_q["QTE"]),
    
    aipw_q1= as.numeric(result_AIPW["q1"]),
    aipw_q0= as.numeric(result_AIPW["q0"]),
    aipw_QTE= as.numeric(result_AIPW["QTE"])
    
  )
  
  
  return(list(
    result = result,
    data   = data
  ))
}
