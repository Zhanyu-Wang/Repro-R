library(testthat)
library(SimBaRepro)

test_that("confidence_grid: validates inputs correctly", {
  set.seed(123)
  n <- 50 # sample size
  R <- 200 # Repro sample size
  alpha <- .05 # significance level
  tol <- 1e-2 # tolerance for the confidence set
  s_obs <- c(1.12, 0.67) # the observed sample mean
  seeds <- matrix(rnorm(R * (n + 2)), nrow = R, ncol = n + 2) # pre-generated seeds

  s_sample <- function(seeds, theta) {
    raw_data <- theta[1] + sqrt(theta[2]) * seeds[, 1:n]

    s_mean <- apply(raw_data, 1, mean)
    s_var <- apply(raw_data, 1, var)

    return(cbind(s_mean, s_var))
  }

  lower_bds <- c(0.5, 0.3)
  upper_bds <- c(1.5, 1.3)

  res = 20

  # resolution is not a positive integer
  expect_error(confidence_grid(alpha = alpha,
                               lower_bds = lower_bds,
                               upper_bds = upper_bds,
                               seeds = seeds,
                               generating_fun = s_sample,
                               s_obs = s_obs,
                               tol = tol,
                               resolution = TRUE),
               "'resolution' must be a positive integer")

  # resolution is not a positive integer
  expect_error(confidence_grid(alpha = alpha,
                               lower_bds = lower_bds,
                               upper_bds = upper_bds,
                               seeds = seeds,
                               generating_fun = s_sample,
                               s_obs = s_obs,
                               tol = tol,
                               resolution = 20.5),
               "'resolution' must be a positive integer")
})

test_that("confidence_grid: check example in function description", {
  set.seed(123)
  n <- 50 # sample size
  R <- 100 # Repro sample size
  alpha <- .05 # significance level
  tol <- 1e-1 # tolerance for the confidence set
  s_obs <- c(1.12, 0.67) # the observed sample mean
  seeds <- matrix(rnorm(R * (n + 2)), nrow = R, ncol = n + 2) # pre-generated seeds

  s_sample <- function(seeds, theta) {
    raw_data <- theta[1] + sqrt(theta[2]) * seeds[, 1:n]

    s_mean <- apply(raw_data, 1, mean)
    s_var <- apply(raw_data, 1, var)

    return(cbind(s_mean, s_var))
  }

  lower_bds <- c(0.5, 0.3)
  upper_bds <- c(1.5, 1.3)

  res = 5

  result <- confidence_grid(alpha, lower_bds, upper_bds, seeds, s_sample, s_obs, tol, res)

  expect_type(result, "list")
  expect_true(all(c("ind_array", "search_lower_bds", "search_upper_bds") %in% names(result)))
  expect_true(is.array(result$ind_array))
  expect_type(result$search_lower_bds, "double")
  expect_type(result$search_upper_bds, "double")
})
