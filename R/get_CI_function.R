

#' get_CI
#'
#' Given the observed statistic, this function computes a confidence interval given that the data generating process is known using the supplied random seeds.
#'
#' @param alpha Significance level.
#' @param lower_bds Vector containing the lower bounds for the search space.
#' @param upper_bds Vector containing the upper bounds for the search space
#' @param j Index indicating the parameter of interest.
#' @param seeds Seeds containing all the randomness.
#' @param G Data generating process.
#' @param s_obs Observed statistic.
#' @param tol tolerance of the confidence interval.
#' @param t_init Starting point for the initial search.
#' @param T_stat A distance metric.
#' @return A vector containing the lower and upper boundary of the confidence interval. In the case when no point is accepted in the search space, return NULL.
#' @examples
#' # DP Normal
#' blahblahblah
#' @export


get_CI <- function(alpha, lower_bds, upper_bds, j, seeds, G, s_obs, tol, t_init = NULL, T_stat = ma_depth) {
  # j indicates that we're computing the confidence interval for the jth parameter
  # tol represents the allowed tolerance on the boundary of our interval
  # T_stat is default to ma_depth as defined in the p_val file
  
  # use the p_value function to identify the best starting point for bisection, if there exits any
  general_search <- p_value_function:::p_value(lower_bds, upper_bds, seeds, G, s_obs, t_init, T_stat)
  
  if (general_search$p_val > alpha) {
    # use the jth coordinate of theta_hat as the starting point for bisection
    beta_init <- general_search$theta_hat[j]
    
    # define the sub_search function which returns whether or not a given interval contains a valid point
    sub_search <- function(beta_left, beta_right, beta_init) {
      
      # print search information for tractability
      print('current bisection interval')
      print(c(beta_left, beta_right))
      print('search starting point')
      print(beta_init)
      
      updated_lower <- lower_bds
      updated_lower[j] <- beta_left
      updated_upper <- upper_bds
      updated_upper[j] <- beta_right
      
      # call the p_value function on the interval (beta_left, beta_right) with initial point
      initial_point <- (updated_lower + updated_upper) / 2
      initial_point[j] <- beta_init
      return(p_value_function:::p_value(updated_lower, updated_upper, seeds, G, s_obs, initial_point, T_stat))
    }
    
    # bisection to find the left boundary
    bisect_left <- function(beta_left, beta_right) {
      left <- beta_left
      right <- beta_right
      mid <- (left + right) / 2
      
      while (right-left > tol/2) {
        # use the midpoint as a starting point
        search_result <- sub_search(left, mid, (left + mid) / 2)
        
        if (search_result$p_val > alpha) { # if accepted
          right <- mid
          mid <- (left + right) / 2
        } else { # if no point was found, try again starting from the right boundary
          search_result <- sub_search(left, mid, mid - tol / 5)
          if (search_result$p_val > alpha) {
            right <- mid
            mid <- (left + right) / 2
          } else {
            left <- mid
            mid <- (left + right) / 2
          }
        }
      }
      return(left)
    }
    
    # bisection to find the right boundary
    bisect_right <- function(beta_left, beta_right) {
      left <- beta_left
      right <- beta_right
      mid <- (left + right) / 2
      
      while (right-left > tol/2) {
        # use the midpoint as a starting point
        search_result <- sub_search(mid, right, (mid + right) / 2)
        
        if (search_result$p_val > alpha) {
          left <- mid
          mid <- (left + right) / 2
        } else {
          search_result <- sub_search(mid, right, mid + tol / 5)
          if (search_result$p_val > alpha) {
            left <- mid
            mid <- (left + right) / 2
          } else {
            right <- mid
            mid <- (left + right)/2
          }
        }
      }
      return(right)
    }
    
    # use bisect_left and bisect_right to compute the left/right boundaries
    beta_L <- bisect_left(lower_bds[j], beta_init)
    beta_R <- bisect_right(beta_init, upper_bds[j])
    
    return(c(beta_L, beta_R))
    
  } else {
    return(NULL) # no point in the parameter space is likely, so return the empty set
  }
}


