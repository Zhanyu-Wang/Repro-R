

# the Mahalanobis depth (modified to also take in a parameter theta)
# x: each row of x is a vector whose distance we want to compute,
# data: an (R+1) by d matrix containing (R+1) vectors
#' @noRd
ma_depth <- function(x, data, theta) {
  return(ddalpha::depth.Mahalanobis(x, data))
}

#' p_value
#'
#' Given the observed statistic and the given seeds, this function finds the p-value.
#' The method uses simulation-based inference, where having fixed seeds, the parameter is searched which makes the observed statistics most "plausible".
#' In particular, the `T_stat` function measures the "plausibility" of any data point and the procedure maximizes the rank of the observed `T_stat` value relative to the "repro" 'T_stat' values.
#' The p-value is determined from the maximum rank and the corresponding parameter is returned.
#'
#' @param lower_bds A vector containing the lower bounds for the parameter search space.
#' @param upper_bds A vector containing the upper bounds for the parameter search space.
#' @param seeds A matrix (or array) of seeds for generating artificial statistics.
#' @param generating_fun A function that takes the random seeds above and a parameter in the search space as inputs to generate artificial statistics.
#' @param s_obs A vector representing the observed statistic.
#' @param theta_init A vector specifying the starting point for the initial `optim` search.
#' @param T_stat  See Vignette for detailed explanation.
#' @param verbose A Boolean variable indicating whether or not to print out the `optim` messages.
#' @param check_input A Boolean variable indicating whether or not to run checks on the function inputs.
#' @return A list containing the most likely parameter in the search region (`theta_hat`) and its corresponding p-value (`p_val`).
#' @examples
#' ### Regular Normal
#' set.seed(123)
#' n <- 50 # sample size
#' R <- 50 # Repro sample size (should be at least 200 for accuracy in practice)
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
#' lower_bds <- c(-5, 0.01) # lower bounds for null hypothesis region
#' upper_bds <- c(5, 5) # upper bounds for null hypothesis region
#'
#' result <- p_value(lower_bds, upper_bds, seeds, s_sample, s_obs)
#' print(result$p_val) # the largest p_value found
#' print(result$theta_hat) # the parameter corresponding to the largest p value
#'
#'
#' @importFrom ddalpha depth.Mahalanobis
#' @importFrom stats optim
#' @export

# p_value function
p_value <- function(lower_bds, upper_bds, seeds, generating_fun, s_obs, theta_init = NULL, T_stat = ma_depth, verbose = FALSE, check_input = TRUE) {

  seeds_dim = dim(seeds)
  # input tests
  if (isTRUE(check_input)) {
    if (length(lower_bds) != length(upper_bds)) {
      stop("Lengths of inputs 'lower_bds' and 'upper_bds' must match.")
    } else if (any(lower_bds > upper_bds)) {
      stop("'lower_bds' must be smaller than or equal to 'upper_bds' entry-wise.")
    } else if (length(seeds_dim) != 2) {
      stop("'seeds' must be a 2-dimensional object (either a matrix or an array).")
    } else if (!is.numeric(seeds) || any(is.na(seeds))) {
      stop("'seeds' must be a numeric matrix or array without NA values.")
    } else if (!is.function(generating_fun)) {
      stop("'generating_fun' must be a function.")
    } else if (length(formals(generating_fun)) != 2) {
      stop("'generating_fun' must be a function with exactly two inputs. The first one is a matrix or an array, the second one is a vector.")
    }
  }

  # check whether 'lower_bds' equals to 'upper_bds' in any entries
  equal_vec <- lower_bds == upper_bds

  # extract the number of seeds R
  R <- seeds_dim[1]
  d <- length(s_obs)

  # a function that generate R simulated values using the seeds and generating_fun, store s_obs and s_sim in an R+1 by d matrix
  s <- function(theta) {
    s_values <- rbind(s_obs, generating_fun(seeds, theta))
    return(s_values)
  }

  # function that computes and stores the simulated statistics using T_stat
  statistics <- function(theta) {
    s_matrix <- s(theta)
    t_vec <- T_stat(s_matrix, s_matrix, theta)
    return(t_vec)
  }

  # define the counting function that we feed into optim
  count <- function(partial_theta) {
    full_theta <- numeric(length(lower_bds))
    full_theta[!equal_vec] <- partial_theta # for variable positions, assign the partial_theta values
    full_theta <- full_theta + equal_vec * lower_bds # for constant positions, assign the constant value

    t_values <- statistics(full_theta)
    ct <- sum(t_values[-1] <= t_values[1]) + t_values[1]
    return(-ct)
  }

  # pick the midpoint if theta_init is not specified
  if (is.null(theta_init)) {
    theta_init <- (lower_bds + upper_bds)/2
  }

  reduced_t <- theta_init[!equal_vec]
  reduced_lower_bds <- lower_bds[!equal_vec]
  reduced_upper_bds <- upper_bds[!equal_vec]

  # call the optim function for minimization
  opt <- optim(par = reduced_t,
               fn = count,
               method = "L-BFGS-B",
               lower = reduced_lower_bds,
               upper = reduced_upper_bds)
  if (isTRUE(verbose)) {
    message(opt$message)
  }
  m <- -opt$value
  theta_hat <- opt$par

  full_theta_hat <- lower_bds
  full_theta_hat[!equal_vec] <- theta_hat

  # compute the p value and return
  p_val <- 1/(R+1) * (min(floor(m), R) + 1)

  # compile a list of values to return
  results <- list(p_val = p_val,
                  theta_hat = full_theta_hat)

  return(results)
}


