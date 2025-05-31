library(testthat)
library(SimBaRepro)

test_that("p_value: input checks work", {
  # example in the function description
  n <- 50
  R <- 200
  s_obs <- c(1.12, 0.67)
  seeds <- matrix(rnorm(R * (n + 2)), nrow = R, ncol = n + 2)

  s_sample <- function(seeds, theta) {
    raw_data <- theta[1] + sqrt(theta[2]) * seeds[, 1:n]
    s_mean <- apply(raw_data, 1, mean)
    s_var <- apply(raw_data, 1, var)
    return(cbind(s_mean, s_var))
  }

  lower_bds <- c(-5, 0.01)
  upper_bds <- c(5, 5)

  bad_G <- function(x, y, z) {
    return(x)
  }

  # lengths of 'lower_bds' and 'upper_bds' mismatch
  expect_error(p_value(lower_bds = c(0, 1),
                       upper_bds = c(1, 2, 3),
                       seeds = seeds,
                       G = s_sample,
                       s_obs = s_obs),
               "Lengths of inputs 'lower_bds' and 'upper_bds' must match.")

  # 'lower_bds' are not smaller than 'upper_bds'
  expect_error(p_value(lower_bds = c(0, 1),
                       upper_bds =c(2, 0),
                       seeds = seeds,
                       G = s_sample,
                       s_obs = s_obs),
               "'lower_bds' must be smaller than or equal to 'upper_bds' entry-wise.")

  # 'seeds' is not a 2d object
  expect_error(p_value(lower_bds = lower_bds,
                       upper_bds = upper_bds,
                       seeds = array(rnorm(20), c(2, 2, 5)),
                       G = s_sample,
                       s_obs = s_obs),
               "'seeds' must be a 2-dimensional object")

  # 'seeds' contains NA values
  expect_error(p_value(lower_bds = lower_bds,
                       upper_bds = upper_bds,
                       seeds = array(rep(NA, 20), c(4, 5)),
                       G = s_sample,
                       s_obs = s_obs),
               "'seeds' must be a numeric matrix or array without NA values.")

  # 'G' is not a function
  expect_error(p_value(lower_bds = lower_bds,
                       upper_bds = upper_bds,
                       seeds = seeds,
                       G = c(1,2),
                       s_obs = s_obs),
               "'G' must be a function.")

  # 'G' is not a function with two inputs
  expect_error(p_value(lower_bds = lower_bds,
                       upper_bds = upper_bds,
                       seeds = seeds,
                       G = bad_G,
                       s_obs = s_obs),
               "'G' must be a function with exactly two inputs. The first one is a matrix or an array, the second one is a vector.")
})

test_that("p_value example runs without error", {
  set.seed(123)
  n <- 50  # sample size
  R <- 200 # Repro sample size
  s_obs <- c(1.12, 0.67)
  seeds <- matrix(rnorm(R * (n + 2)), nrow = R, ncol = n + 2)

  s_sample <- function(seeds, theta) {
    raw_data <- theta[1] + sqrt(theta[2]) * seeds[, 1:n]
    s_mean <- apply(raw_data, 1, mean)
    s_var <- apply(raw_data, 1, var)
    return(cbind(s_mean, s_var))
  }

  lower_bds <- c(-5, 0.01)
  upper_bds <- c(5, 5)

  result <- p_value(lower_bds, upper_bds, seeds, s_sample, s_obs)

  expect_type(result$p_val, "double")
  expect_true(length(result$theta_hat) == 2)
})
