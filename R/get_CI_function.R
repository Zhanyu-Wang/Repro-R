
#' get_CI
#'
#' Given the observed statistic, this function computes a confidence interval given that the data generating process is known using the supplied random seeds.
#'
#' @param alpha A numeric representing the significance level of the test.
#' @param lower_bds A vector containing the lower bounds for the parameter search space.
#' @param upper_bds A vector containing the upper bounds for the parameter search space.
#' @param parameter_index An integer indicating the parameter of interest.
#' @param seeds A matrix (or array) of seeds for generating artificial statistics.
#' @param generating_fun A function that takes the random seeds above and a parameter in the search space as inputs to generate artificial statistics.
#' @param s_obs A vector representing the observed statistic.
#' @param tol A numeric specifying the tolerance of the confidence interval.
#' @param theta_init A vector specifying the starting point for the initial `optim` search.
#' @param T_stat Default to the Mahalanobis distance. See Vignette for detailed explanation.
#' @param verbose A Boolean variable indicating whether or not to print out the `optim` messages.
#' @param check_input A Boolean variable indicating whether or not to run checks on the function inputs.
#' @return A length-2 vector representing the obtained confidence interval. In the case when no point is accepted in the search space, return `NULL`.
#' @examples
#' ### Note that the examples may take a few seconds to run.
#' ### Regular Normal
#' set.seed(123)
#' n <- 30 # sample size
#' R <- 50 # Repro sample size (should be at least 200 for accuracy in practice)
#' alpha <- .05 # significance level
#' tol <- 0.01 # tolerance for the confidence set (use smaller tolerance in practice)
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
#' # choose parameter_index = 1 to get the confidence interval for the mean
#' mean_CI <- get_CI(alpha, lower_bds, upper_bds, 1, seeds, s_sample, s_obs, tol)
#' print(mean_CI) # estimated confidence interval for mean
#' var_CI <- get_CI(alpha, lower_bds, upper_bds, 2, seeds, s_sample, s_obs, tol)
#' print(var_CI) # estimated confidence interval for variance
#'
#' @export


get_CI <- function(alpha, lower_bds, upper_bds, parameter_index, seeds, generating_fun, s_obs, tol, theta_init=NULL, T_stat=ma_depth, verbose=FALSE, check_input=TRUE) {
  # parameter_index indicates which parameter we're computing the confidence interval for
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
      warning("A large 'tol' might lead to inaccuracies in the result.")
    }
  }

  # use the p_value function to identify the best starting point for bisection, if there exits any
  if (isTRUE(verbose)) {
    message("Initial search")
  }
  general_search <- p_value(lower_bds, upper_bds, seeds, generating_fun, s_obs, theta_init, T_stat, verbose, check_input)

  if (general_search$p_val > alpha) {
    # use the jth coordinate of theta_hat as the starting point for bisection
    beta_init <- general_search$theta_hat[parameter_index]

    # define the sub_search function which returns whether or not a given interval contains a valid point
    sub_search <- function(beta_left, beta_right, beta_init) {

      # print search information for tractability
      if (isTRUE(verbose)) {
        message("Current bisection interval: [", beta_left, ", ", beta_right, "]")
        message("Search starting point: ", beta_init)
      }

      updated_lower <- lower_bds
      updated_lower[parameter_index] <- beta_left
      updated_upper <- upper_bds
      updated_upper[parameter_index] <- beta_right

      # call the p_value function on the interval (beta_left, beta_right) with initial point
      initial_point <- (updated_lower + updated_upper) / 2
      initial_point[parameter_index] <- beta_init
      return(p_value(updated_lower, updated_upper, seeds, generating_fun, s_obs, initial_point, T_stat, verbose, check_input=FALSE))
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
    beta_L <- bisect_left(lower_bds[parameter_index], beta_init)
    beta_R <- bisect_right(beta_init, upper_bds[parameter_index])

    return(c(beta_L, beta_R))

  } else {
    return(NULL) # no point in the parameter space is likely, so return the empty set
  }
}


