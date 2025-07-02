library(testthat)
library(SimBaRepro)

test_that("get_CI: input checks", {
  n <- 20
  R <- 100
  alpha <- .05
  tol <- 1e-1
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

  # 'alpha' is not a number
  expect_error(get_CI(alpha = "a",
                      lower_bds = lower_bds,
                      upper_bds = upper_bds,
                      parameter_index = 1,
                      seeds = seeds,
                      generating_fun = s_sample,
                      s_obs = s_obs,
                      tol = tol),
               "Significance level 'alpha' must be a number")

  # 'alpha' is not in between 0 and 1
  expect_error(get_CI(alpha = 1.1,
                      lower_bds = lower_bds,
                      upper_bds = upper_bds,
                      parameter_index = 1,
                      seeds = seeds,
                      generating_fun = s_sample,
                      s_obs = s_obs,
                      tol = tol),
               "Significance level 'alpha' must be a number between 0 and 1.")

  # 'tol' is not a positive number
  expect_error(get_CI(alpha = alpha,
                      lower_bds = lower_bds,
                      upper_bds = upper_bds,
                      parameter_index = 1,
                      seeds = seeds,
                      generating_fun = s_sample,
                      s_obs = s_obs,
                      tol = TRUE),
               "'tol' must be a positive number.")
})

test_that("example in function description runs without issue", {
  set.seed(123)
  # use smaller sample size for faster testing
  n <- 20
  R <- 50
  alpha <- .05
  tol <- 1e-1
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

  mean_CI <- get_CI(alpha, lower_bds, upper_bds, 1, seeds, s_sample, s_obs, tol)
  var_CI <- get_CI(alpha, lower_bds, upper_bds, 2, seeds, s_sample, s_obs, tol)

  expect_type(mean_CI, "double")
  expect_type(var_CI, "double")
  expect_true(mean_CI[1] < mean_CI[2])
  expect_true(var_CI[1] < var_CI[2])
})
