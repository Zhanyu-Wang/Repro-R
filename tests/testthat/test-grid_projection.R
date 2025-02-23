library(testthat)
library(SimBaRepro)

test_that("grid_projection returns correct 1D projection for index_set of length 1", {
  ind_arr <- array(c(1, 0, 0, 0, 0, 1, 0, 1), dim = rep(2, 3))
  proj <- grid_projection(ind_arr, 1)

  expect_true(is.vector(proj))
  expect_type(proj, "integer")
})

test_that("grid_projection returns correct 2D projection for index_set of length 2", {
  ind_arr <- array(c(1, 0, 0, 0, 0, 1, 0, 1), dim = rep(2, 3))
  proj <- grid_projection(ind_arr, c(1, 2))

  expect_true(is.matrix(proj) || is.array(proj))
  expect_true(length(dim(proj)) == 2)
  expect_type(proj, "integer")
})

test_that("grid_projection errors on invalid index_set length", {
  indicator_array <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
  expect_error(grid_projection(indicator_array, c(1, 2, 3)),
               "A vector of length 1 or 2 is expected for 'index_set'.")
})
