

# import this library for Mahalnobis depth
library(ddalpha)



# the Mahalanobis depth (modified to also take in a parameter theta)
# x: each row of x is a vector whose distance we want to compute, 
# data: an (R+1) by d matrix containing (R+1) vectors
ma_depth <- function(x, data, theta){
  return(depth.Mahalanobis(x, data))
}

# test case for ma_depth
data_test <- matrix(data = rnorm(50),
                 nrow = 10,
                 ncol =5)
x_test <- data_test[1:2,]
ma_depth(x_test, data_test, 0)



# p_value function 
p_value <- function(lower_bd, upper_bd, t_init, seeds, G, s_obs, T_stat = ma_depth){
  # extract the number of seeds R
  dimensions <- dim(seeds)
  R <- dimensions[1]
  d <- dimensions[2]
  
  # a function that generate R simulated values using the seeds and G, store s_obs and s_sim in an R+1 by d matrix
  s <- function(theta){
    s_values <- matrix(s_obs, 
                       nrow = R+1, 
                       ncol = d, 
                       byrow = TRUE)
    for(i in 1:R){
      s_values[i,] <- G(seeds[i-1,], theta)
    }
  }
  
  # function that computes and stores the simulated statistics using T_stat
  statistics <- function(theta){
    t_vec <- T_stat(s(theta), s(theta), theta)
  }
  
  # define the counting function that we feed into optim
  count <- function(theta){
    s_values <- s(theta)
    t_values <- statistics(s_values, s_values, theta)
    ct <- sum(t_values[-1] <= t_values[1]) + t_values[1]
    return(-ct)
  }
  
  # finding the largest 
  max <- -optim(par = t_init, 
                fn = count,
                lower = lower_bd,
                uppper = upper_bd)$value
  
  # compute the p value and return
  p_val <- 1/(R+1) * min(floor(max)+1, R+1)
  return(p_val)
}



# accept function
accept <- function(alpha, lower_bd, upper_bd, t_init, seeds, G, s_obs, T_stat = ma_depth){
  p_val <- p_value(lower_bd, upper_bd, t_init, seeds, G, s_obs, T_stat = ma_depth)
  
  # return true if we fail to reject
  return(p_val > alpha)
}












