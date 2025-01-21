
#' p_value
#'
#' This is an internal function that uses `depth.Mahalanobis` from the `ddalpha` package.
#'
#' @importFrom ddalpha depth.Mahalanobis

# the Mahalanobis depth (modified to also take in a parameter theta)
# x: each row of x is a vector whose distance we want to compute, 
# data: an (R+1) by d matrix containing (R+1) vectors
ma_depth <- function(x, data, theta) {
  return(ddalpha::depth.Mahalanobis(x, data))
}


# p_value function 
p_value <- function(lower_bds, upper_bds, seeds, G, s_obs, t_init = NULL, T_stat = ma_depth) {
  # extract the number of seeds R
  R <- dim(seeds)[1]
  d <- length(s_obs)
  
  # a function that generate R simulated values using the seeds and G, store s_obs and s_sim in an R+1 by d matrix
  s <- function(theta) {
    s_values <- rbind(s_obs, G(seeds, theta))
  }
  
  # function that computes and stores the simulated statistics using T_stat
  statistics <- function(theta) {
    s_matrix <- s(theta)
    t_vec <- T_stat(s_matrix, s_matrix, theta)
    
    return(t_vec)
  }
  
  # define the counting function that we feed into optim
  count <- function(theta) {
    t_values <- statistics(theta)
    ct <- sum(t_values[-1] <= t_values[1]) + t_values[1]
    return(-ct)
  }
  
  # pick the midpoint if t_init is not specified
  if (is.null(t_init)) {
    t_init <- (lower_bds + upper_bds)/2
  }
  
  # call the optim function for minimization
  opt <- optim(par = t_init, 
               fn = count,
               method = "L-BFGS-B",
               lower = lower_bds,
               upper = upper_bds)
  print('opt_result')
  print(opt)
  m <- -opt$value
  theta_hat <- opt$par
  
  # compute the p value and return
  p_val <- 1/(R+1) * (min(floor(m), R) + 1)
  
  # compile a list of values to return
  results <- list(p_val = p_val,
                  rank = floor(m)+1,
                  theta_hat = theta_hat)
  
  return(results)
}


