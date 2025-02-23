
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
#' @param tol Tolerance of the confidence interval.
#' @param t_init Starting point for the initial search.
#' @param T_stat A distance metric (default to ma_depth which measures the Mahalanobis distance).
#' @param print_info Whether or not to print the optim information
#' @param check_input Whether or not to run checks on the function inputs.
#' @return A vector containing the lower and upper boundary of the confidence interval. In the case when no point is accepted in the search space, return NULL.
#' @examples
#' ### Note that the examples may take a few seconds to run.
#' ### Regular Normal
#' set.seed(123)
#' n <- 50 # sample size
#' R <- 200 # Repro sample size
#' alpha <- .05 # significance level
#' tol <- 1e-2 # tolerance for the confidence set
#' s_obs <- c(1.12, 0.67) # the observed sample mean and variance
#' seeds <- matrix(rnorm(R * (n + 2)), nrow = R, ncol = n + 2) # pre-generated seeds
#'
#' # this function computes the repro statistics given the seeds and the parameter
#' s_sample <- function(seeds, theta) {
#' # generate the raw data points
#'   raw_data <- theta[1] + sqrt(theta[2]) * seeds[, 1:n]
#'
#'   # compute the regular statistics
#'   s_mean <- apply(raw_data, 1, mean)
#'   s_var <- apply(raw_data, 1, var)
#'
#'   return(cbind(s_mean, s_var))
#' }
#'
#' lower_bds <- c(-5, 0.01) # lower bounds for the parameter search region
#' upper_bds <- c(5, 5) # upper bounds for the parameter search region
#'
#' # choose j = 1 to get the confidence interval for the mean
#' mean_CI <- get_CI(alpha, lower_bds, upper_bds, 1, seeds, s_sample, s_obs, tol)
#' print(mean_CI) # estimated confidence interval for mean
#' var_CI <- get_CI(alpha, lower_bds, upper_bds, 2, seeds, s_sample, s_obs, tol)
#' print(var_CI) # estimated confidence interval for variance
#'
#' @export


get_CI <- function(alpha, lower_bds, upper_bds, j, seeds, G, s_obs, tol, t_init=NULL, T_stat=ma_depth, print_info=FALSE, check_input=TRUE) {
  # j indicates that we're computing the confidence interval for the jth parameter
  # tol represents the allowed tolerance on the boundary of our interval
  # T_stat is default to ma_depth as defined in the p_val file

  if (isTRUE(check_input)) {
    # input check
    if (!is.numeric(alpha)) {
      stop("Significance level 'alpha' must be a number.")
    } else if (alpha > 1 || alpha < 0) {
      stop("Significance level 'alpha' must be a number between 0 and 1.")
    } else if (!is.numeric(tol)) {
      stop("'tol' must be a positive number.")
    } else if (tol > 1) {
      print("A large 'tol' might lead to inaccuracies in the result.")
    }
  }

  # use the p_value function to identify the best starting point for bisection, if there exits any
  general_search <- p_value(lower_bds, upper_bds, seeds, G, s_obs, t_init, T_stat, print_info, check_input)

  if (general_search$p_val > alpha) {
    # use the jth coordinate of theta_hat as the starting point for bisection
    beta_init <- general_search$theta_hat[j]

    # define the sub_search function which returns whether or not a given interval contains a valid point
    sub_search <- function(beta_left, beta_right, beta_init) {

      # print search information for tractability
      if (isTRUE(print_info)) {
        print('current bisection interval')
        print(c(beta_left, beta_right))
        print('search starting point')
        print(beta_init)
      }

      updated_lower <- lower_bds
      updated_lower[j] <- beta_left
      updated_upper <- upper_bds
      updated_upper[j] <- beta_right

      # call the p_value function on the interval (beta_left, beta_right) with initial point
      initial_point <- (updated_lower + updated_upper) / 2
      initial_point[j] <- beta_init
      return(p_value(updated_lower, updated_upper, seeds, G, s_obs, initial_point, T_stat, print_info, check_input=FALSE))
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


