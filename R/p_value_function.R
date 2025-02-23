

# the Mahalanobis depth (modified to also take in a parameter theta)
# x: each row of x is a vector whose distance we want to compute,
# data: an (R+1) by d matrix containing (R+1) vectors
ma_depth <- function(x, data, theta) {
  return(ddalpha::depth.Mahalanobis(x, data))
}

#' p_value
#'
#' Given the observed statistic and the given seeds, this function finds the largest p_value and the corresponding parameter
#'
#' @param lower_bds Vector containing the lower bounds for the search space.
#' @param upper_bds Vector containing the upper bounds for the search space
#' @param seeds Seeds containing all the randomness.
#' @param G Data generating function.
#' @param s_obs Observed statistic.
#' @param t_init Starting point for the initial search.
#' @param T_stat A distance metric (default to ma_depth which measures the Mahalanobis distance).
#' @param print_info Whether or not to print the optim information
#' @param check_input Whether or not to run checks on the function inputs.
#' @return This function returns a list containing the minimum p value (p_val) within the search region and the parameter corresponding to it (theta_hat).
#' @examples
#' ### Regular Normal
#' set.seed(123)
#' n <- 50 # sample size
#' R <- 200 # Repro sample size
#' s_obs <- c(1.12, 0.67) # the observed sample mean and variance
#' seeds <- matrix(rnorm(R * (n + 2)), nrow = R, ncol = n + 2) # pre-generated seeds
#'
#' # this function computes the repro statistics given the seeds and the parameter
#' s_sample <- function(seeds, theta) {
#'   # generate the raw data points
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
#' result <- p_value(lower_bds, upper_bds, seeds, s_sample, s_obs)
#' print(result$p_val) # the largest p_value found
#' print(result$theta_hat) # the parameter corresponding to the largest p value
#'
#'
#' @importFrom ddalpha depth.Mahalanobis
#' @export

# p_value function
p_value <- function(lower_bds, upper_bds, seeds, G, s_obs, t_init = NULL, T_stat = ma_depth, print_info = FALSE, check_input = TRUE) {

  seeds_dim = dim(seeds)
  # input tests
  if (isTRUE(check_input)) {
    if (length(lower_bds) != length(upper_bds)) {
      stop("Lengths of inputs 'lower_bds' and 'upper_bds' must match.")
    } else if (any(lower_bds >= upper_bds)) {
      stop("'lower_bds' must be smaller than 'upper_bds' at all entries (elementwise).")
    } else if (length(seeds_dim) != 2) {
      stop("'seeds' must be a 2-dimensional object (either a matrix or an array).")
    } else if (!is.numeric(seeds) || any(is.na(seeds))) {
      stop("'seeds' must be a numeric matrix or array without NA values.")
    } else if (!is.function(G)) {
      stop("'G' must be a function.")
    } else if (length(formals(G)) != 2) {
      stop("'G' must be a function with exactly two inputs. The first one is a matrix or an array, the second one is a vector.")
    } else if (length(s_obs) != length(lower_bds)) {
      stop("'s_obs' must have the same length as 'lower_bds' and 'upper_bds'.")
    }
  }

  # extract the number of seeds R
  R <- seeds_dim[1]
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
  if (isTRUE(print_info)) {
    print("opt result")
    print(opt)
  }
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


