#alpha is a number in [1/(R+1),1)

#thetaLow is a vector of length d containing the lower bounds for all d attributes of theta

#thetaHigh is a vector of length d containing the upper bounds for all d attributes of theta

#b is a natural number in [1,d]. We are creating a confidence interval for the bth attribute 
#of theta

#seeds is a vector of length R and is randomly generated

#s is the observed data and is a vector of length d

#T is the exchangeable statistic (what data type???)

#tol is the tolerance, a positive number


confidenceInterval <- function(alpha, thetaLow, thetaHigh, b, seeds, G, s, T = ma_depth, tol) {
  #find beta_init
  #beta_init = argmaxβ{supη{#{Tobs((β,η)) ≥ T(i)((β,η))}}} (pg 20, what does this mean???)
  
  #the interval "B" from the paper would here be [thetaLow[b],thetaHigh[b]]
  
  #Algorithm 2 line 1 says "If there exists beta_init in B" - does this mean that we s
  beta_init = generate_beta_init()
  while((accept(beta_init) == FALSE)) {
    beta_init = generate_beta_init #regenerate beta_init until accept returns true
    #Under what conditions is a valid beta_init guaranteed to exist???
    #Should there be a hard cap on the number of times we generate a beta_init?
  } 
  beta_lower_inf = generate_beta_lower_inf()
  while(accept([inf_B, beta_lower_inf]) == TRUE) {
    beta_lower_inf = generate_beta_lower_inf()
  }
  beta_lower_sup = beta_init
  while(beta_lower_sup - beta_lower_inf > tol) {
    beta_lower_mid = .5*(beta_lower_inf + beta_lower_sup)
    if(accept([beta_lower_inf, beta_lower_mid]) == TRUE) {
      beta_lower_sup = beta_lower_mid
    } else {
      beta_lower_inf = beta_lower_mid
    }
  }
}

















