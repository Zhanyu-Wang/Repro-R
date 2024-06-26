

# import this library for Mahalnobis depth
library(ddalpha)



# the Mahalanobis depth (modified to also take in a parameter theta)
# x: each row of x is a vector whose distance we want to compute, 
# data: an (R+1) by d matrix containing (R+1) vectors
ma_depth <- function(x, data, theta) {
  return(depth.Mahalanobis(x, data))
}

# test case for ma_depth
data_test <- matrix(data = rnorm(50),
                 nrow = 10,
                 ncol =5)
x_test <- data_test[1:2,]
ma_depth(x_test, data_test, 0)



# p_value function 
p_value <- function(lower_bds, upper_bds, seeds, G, s_obs, t_init = NULL, T_stat = ma_depth) {
  # extract the number of seeds R
  R <- dim(seeds)[1]
  d <- length(s_obs)
  
  # a function that generate R simulated values using the seeds and G, store s_obs and s_sim in an R+1 by d matrix
  s <- function(theta) {
    s_values <- rbind(s_obs, G(seeds, theta))
    return(s_values)
  }
  
  # function that computes and stores the simulated statistics using T_stat
  statistics <- function(theta) {
    s_matrix <- s(theta)
    t_vec <- T_stat(s_matrix, s_matrix, theta)
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
  
  # finding the largest 
  opt <- optim(par = t_init, 
               fn = count,
               method = "L-BFGS-B",
               lower = lower_bds,
               upper = upper_bds)
  m <- -opt$value
  theta_hat <- opt$par
  
  # compute the p value and return
  p_val <- 1/(R+1) * min(floor(m)+1, R+1)
  
  # compile a list of values to return
  results <- list(p_val = p_val,
                  rank = floor(m)+1,
                  theta_hat = theta_hat)
  
  return(results)
}



# accept function
accept <- function(alpha, lower_bds, upper_bds, seeds, G, s_obs, t_init = NULL, T_stat = ma_depth) {
  # calling the p_value function as a subroutine
  p_val <- p_value(lower_bds, upper_bds, seeds, G, s_obs, t_init, T_stat)$p_val
  
  # return true if we fail to reject
  return(p_val > alpha)
}










