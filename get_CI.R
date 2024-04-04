

get_CI <- function(alpha, lower_bds, upper_bds, j, seeds, G, s_obs, tol, T_stat = ma_depth){
  # j indicates that we're computing the confidence interval for the jth parameter
  # tol represents the allowed tolerance on the boundary of our interval
  # T_stat is default to ma_depth as defined in the p_val file
  
  # use the p_value function to identify the best starting point for bisection, if there exits any
  general_search = p_value(lower_bds, upper_bds, seeds, G, s_obs, T_stat)
  
  if (general_search$p_val > alpha) {
    # use the jth coordinate of theta_hat as the starting point for bisection
    beta_init = general_search$theta_hat[j]
    
    # define the sub_search function which returns whether or not a given interval contains a valid point
    sub_search <- function(beta_left, beta_right){
      updated_lower = lower_bds
      updated_lower[j] = beta_left
      updated_upper = upper_bds
      updated_upper[j] = beta_right
      
      # call the accept function to see if (beta_left, beta_right) contains a valid point
      return(accept(alpha, updated_lower, updated_upper, seeds, G, s_obs, T_stat))
    }
    
    # bisection to find the left boundary
    bisect_left <- function(beta_left, beta_right){
      left = beta_left
      right = beta_right
      right_previous = beta_right
      
      # I personally want to change this to tol/2 to make the overall error tol
      while (right-left > tol/2) { 
        mid = (left + right)/2
        if (sub_search(left, right)) {
          right_previous = right
          right = mid
        } else {
          left = right
          right = (right + right_previous)/2
        }
      }
      return(right)
    }
    
    # bisection to find the right boundary
    bisect_right <- function(beta_left, beta_right){
      left = beta_left
      left_previous = beta_left
      right = beta_right
      
      while (right-left > tol/2) { 
        mid = (left + right)/2
        if (sub_search(left, right)) {
          left_previous = left
          left = mid
        } else {
          right = left
          left = (left_previous + left)/2
        }
      }
      return(left)
    }
    
    # use bisect_left and bisect_right to compute the left/right boundaries
    beta_L = bisect_left(lower_bds[j], beta_init)
    beta_R = bisect_right(beta_init, upper_bds[j])
    
    return(c(beta_L, beta_R))
    
  } else {
    return(NULL) # no point in the parameter space is likely, so return the empty set
  }
}









































