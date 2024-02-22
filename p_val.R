

# import this library for Mahalnobis depth
library(ddalpha)



# the Mahalanobis depth 
# x: a vector of dimension d, data: an (R+1) by d matrix containing (R+1) vectors
ma_depth <- function(x, data, theta){
  return(depth.Mahalanobis(x, data))
}



# p_value function 
p_value <- function(range, t_init, alpha, R, seeds, G, s_obs, T_stat = ma_depth){
  s_sim <- function(t, i){
    return(G(seeds[i, ], t))
  }
  
  # function that computes the simulated s values
  data <- function(t){
    data_array <- array(s_obs, c(1, length(s)))
    for(i in seq(R)){
      data_array <- rbind(data_array, s_sim(t, i))
    }
  }
  
  count <- function(t){
    s_values <- data(t)
    t_obs <- T_stat(s_obs, s_values, t)
    max <- 1 + t_obs
    for(i in seq(R)){
      t_i <- T_stat(s_values[1 + i, ], s_values, t)
      if(t_i <= t_obs){
        max <- max + 1
      }
    }
    return(-max)
  }
  
  max <- -optim(par = t_init, 
                fn = count,
                lower = range[,1],
                uppper = range[,2])$value
  
  if(max >= floor(alpha*(R+1)) + 1){
    return(TRUE)
  }else{
    return(FALSE)
  }
}















