single_sim= function(n= 1000,  y_threshold = 7.95,
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
 
  AIPW = aipw_cpm_pte(data,  y_threshold = y_threshold, formA = formA, formY = formY)
  OR = or_cpm_pte(data,  y_threshold = y_threshold,  formY = formY)
  IPW = ipw_cpm_pte(data,  y_threshold = y_threshold,  formA = formA)
  
  
  
  # Combine the results into a named vector
  result= c(
    AIPW_Fy1_hat  = as.numeric(unname(AIPW["FY1"])),
    AIPW_Fy0_hat  = as.numeric(unname(AIPW["FY0"])),
    AIPW_pte_hat = as.numeric(unname(AIPW["FY1"])) - as.numeric(unname(AIPW["FY0"])),
    
    OR_Fy1_hat  = as.numeric(unname(OR["FY1"])),
    OR_Fy0_hat  = as.numeric(unname(OR["FY0"])),
    OR_pte_hat = as.numeric(unname(OR["FY1"])) - as.numeric(unname(OR["FY0"])),
    
    IPW_Fy1_hat  = as.numeric(unname(IPW["FY1"])),
    IPW_Fy0_hat  = as.numeric(unname(IPW["FY0"])),
    IPW_pte_hat = as.numeric(unname(IPW["FY1"]))- as.numeric(unname(IPW["FY0"]))
  )
  
  
  return(list(
    result = result,
    data   = data
  ))
}
