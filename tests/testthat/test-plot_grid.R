library(testthat)
library(ggplot2)
library(SimBaRepro)

test_that("plot_grid errors", {
  ind_array <- array(1, dim = c(2, 2))
  lower_bds <- c(0, 0)
  upper_bds <- c(1, 1)

  # 'indicator_array' is not a 2-dimensional object
  expect_error(plot_grid(indicator_array = array(1, dim = c(2, 2, 2)),
                         lower_bds = lower_bds,
                         upper_bds = upper_bds),
               "'indivator_array' needs to be a 2-dimensional object.")

  # 'lower_bds' or 'upper_bds' are not of length 2
  expect_error(plot_grid(indicator_array = ind_array,
                         lower_bds = c(0, 0, 0),
                         upper_bds = c(1, 1, 1)),
               "Lengths of inputs 'lower_bds' and 'upper_bds' must both be 2.")

  # 'lower_bds' needs to be smaller than 'upper_bds'
  expect_error(plot_grid(indicator_array = ind_array,
                         lower_bds = c(0, 1),
                         upper_bds = c(1, 1)),
               "'lower_bds' must be smaller than 'upper_bds' at all entries.")
})

test_that("plot_grid example", {
  ind_array <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
  lower_bds <- c(0, 0)
  upper_bds <- c(1, 1)

  p <- plot_grid(ind_array, lower_bds, upper_bds)
  expect_true(inherits(p, "ggplot"))
})
