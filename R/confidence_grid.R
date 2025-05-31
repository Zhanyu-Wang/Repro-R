

#' confidence_grid
#'
#' returns the indicator array
#'
#' @param alpha A numeric representing the significance level of the test.
#' @param lower_bds A vector containing the lower bounds for the parameter search space.
#' @param upper_bds A vector containing the upper bounds for the parameter search space.
#' @param seeds A matrix (or array) of seeds for generating artificial statistics.
#' @param generating_fun A function that takes the random seeds above and a parameter in the search space as inputs to generate artificial statistics.
#' @param s_obs A vector representing the observed statistic.
#' @param tol A numeric specifying the tolerance of the confidence interval.
#' @param resolution An integer specifying the mesh number of the search space.
#' @param theta_init A vector specifying the starting point for the initial `optim` search.
#' @param T_stat Default to the Mahalanobis distance. See Vignette for detailed explanation.
#' @return A list containing an indicator array (`ind_array`) representing the confidence set, the confidence set lower bounds (`updated_lower_bds`), and the confidence set upper bounds (`updated_upper_bds`).
#' @examples
#' ### Note that the examples may take a few seconds to run.
#' ### Regular normal
#' set.seed(123)
#' n <- 50 # sample size
#' R <- 200 # Repro sample size
#' alpha <- .05 # significance level
#' tol <- 1e-2 # tolerance for the confidence set
#' s_obs <- c(1.12, 0.67) # the observed sample mean
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
#' lower_bds <- c(0.5, 0.3) # lower bounds for the parameter search region
#' upper_bds <- c(1.5, 1.3) # upper bounds for the parameter search region
#'
#' resolution = 20  # resolution of the grid
#' result <- confidence_grid(alpha, lower_bds, upper_bds, seeds, s_sample, s_obs, tol, resolution)
#' print(result$ind_array)
#' print(result$search_lower_bds)
#' print(result$search_upper_bds)
#'
#' @export

confidence_grid <- function(alpha, lower_bds, upper_bds, seeds, generating_fun, s_obs, tol, resolution, theta_init = NULL, T_stat = ma_depth) {

  # for each parameter, run the get_CI function to obtain a smaller search bounds
  p <- length(lower_bds)
  search_lower_bds <- rep(0, p)
  search_upper_bds <- rep(0, p)

  # input check
  if (!is.numeric(resolution)) {
    stop("'resolution' must be a positive integer.")
  } else if (resolution %% 1 != 0 || resolution <= 0) {
    stop("'resolution' must be a positive integer.")
  }

  for (j in 1:p) {
    if (j == 1) {
      check_input <- TRUE
    } else {
      check_input <- FALSE
    }
    confidence_interval <- get_CI(alpha, lower_bds, upper_bds, j, seeds, generating_fun, s_obs, tol, theta_init, T_stat, print_info = FALSE, check_input)
    search_lower_bds[j] <- confidence_interval[1]
    search_upper_bds[j] <- confidence_interval[2]
  }

  # resolution gives the number of divisions for each parameter, divide by resolution to get the grid width
  grid_width <- (upper_bds - lower_bds) / resolution

  # an array representing whether or not a grid contains an accepted value
  indicator_array <- array(rep(0, resolution ** p), dim = rep(resolution, p))

  # given search bounds, find the boundary indices in the big indicator array to search over
  lower_indices <- (search_lower_bds - lower_bds) %/% grid_width + 1
  upper_indices <- resolution - (upper_bds - search_upper_bds) %/% grid_width

  # given the indices, this function search the corresponding cube and update the indicator array
  update_array <- function(indices, indicator_array) {
    cube_lower_bds <- lower_bds + (indices - 1) * grid_width
    cube_upper_bds <- cube_lower_bds + grid_width
    mid_point <- (cube_lower_bds + cube_upper_bds) / 2

    if (p_value(cube_lower_bds, cube_upper_bds, seeds, generating_fun, s_obs, mid_point, T_stat, print_info = FALSE, check_input = FALSE)$p_val > alpha) {
      indicator_array <- do.call('[<-', c(list(indicator_array), as.list(indices), 1))
    }

    return(indicator_array)
  }

  # function that recursively loops through all sub-divided cubes within search region,
  loop_thru <- function(curr_pos, indicator_array, former_indices = c()) {
    if (curr_pos <= p) {
      for (i in lower_indices[curr_pos]:upper_indices[curr_pos]) {
        new_indices <- c(former_indices, i)
        indicator_array <- loop_thru(curr_pos + 1, indicator_array, new_indices)
      }
    } else {
      indicator_array <- update_array(former_indices, indicator_array)
    }
    return(indicator_array)
  }

  # call the recursive function defined above to update
  indicator_array <- loop_thru(1, indicator_array)

  return(list(ind_array = indicator_array,
              search_lower_bds = search_lower_bds,
              search_upper_bds = search_upper_bds))
}

#' grid_projection
#'
#' Projects the multidimensional indicator array generated by confidence_grid down to 2d for visualization
#'
#' @param indicator_array An indicator array generated using the `confidence_grid` function.
#' @param index_set A vector containing the indices representing the dimensions to keep.
#' @return A two-dimensional indicator array ready for visualization (`indicator_array` projected onto the subspace specified by `index_set`).
#' @examples
#' ### simple projection
#' ind_arr <- array(c(1, 0, 0, 0, 0, 1, 0, 1), dim = rep(2, 3))
#' print(ind_arr)
#' # project this indicator array onto a 2d subspace by first and second dimension
#' ind_arr_12 <- grid_projection(ind_arr, c(1,2))
#' print(ind_arr_12)
#' ind_arr_13 <- grid_projection(ind_arr, c(1,3))
#' print(ind_arr_13)
#' ind_arr_23 <- grid_projection(ind_arr, c(2,3))
#' print(ind_arr_23)
#' @export

# this function computes the projection of the indicator_array onto 2d or 1d arrays
grid_projection <- function(indicator_array, index_set) {
  if (length(index_set) == 1 | length(index_set) == 2) {
    return(apply(indicator_array, index_set, function(x) as.integer(any(x==1))))
  } else {
    stop("A vector of length 1 or 2 is expected for 'index_set'.")
  }
}


#' plot_grid
#'
#' projects the indicator array generated by confidence_grid down to 2d for visualization
#'
#' @param indicator_array An 2-dimensional indicator array generated by `confidence_grid` or `grid_projection`.
#' @param lower_bds A vector containing the lower bounds for the parameter search space.
#' @param upper_bds A vector containing the upper bounds for the parameter search space.
#' @return A grid plot showing the confidence regions.
#' @examples
#' ### Note that the examples may take a few seconds to run.
#' ### Regular normal
#' set.seed(123)
#' n <- 50 # sample size
#' R <- 200 # Repro sample size
#' alpha <- .05 # significance level
#' tol <- 1e-2 # tolerance for the confidence set
#' s_obs <- c(1.12, 0.67) # the observed sample mean
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
#' lower_bds <- c(0.5, 0.4) # lower bounds for the parameter search region
#' upper_bds <- c(1.5, 1.4) # upper bounds for the parameter search region
#'
#' resolution = 20  # resolution of the grid
#'
#' result <- confidence_grid(alpha, lower_bds, upper_bds, seeds, s_sample, s_obs, tol, resolution)
#' ind_arr <- result$ind_array
#' plot_grid(ind_arr, lower_bds, upper_bds)
#'
#' @import ggplot2
#' @export

# given the grid_projection (1d or 2d) and the bounds, this functions plots the confidence set
plot_grid <- function(indicator_array, lower_bds, upper_bds) {
  dims = dim(indicator_array)

  # check inputs
  if (length(dims) != 2) {
    stop("'indivator_array' needs to be a 2-dimensional object.")
  } else if (length(lower_bds) != 2 || length(upper_bds) != 2) {
    stop("Lengths of inputs 'lower_bds' and 'upper_bds' must both be 2.")
  } else if (any(lower_bds >= upper_bds)) {
    stop("'lower_bds' must be smaller than 'upper_bds' at all entries.")
  }

  # get the dimensions
  n_rows <- dims[1]
  n_cols <- dims[2]

  # Calculate the width and height of each subdivided region
  width <- (upper_bds[1] - lower_bds[1]) / n_rows
  height <- (upper_bds[2] - lower_bds[2]) / n_cols

  # Create a data frame to store the rectangles
  rectangles <- data.frame(
    xmin = numeric(),
    xmax = numeric(),
    ymin = numeric(),
    ymax = numeric(),
    fill = character()
  )

  # Loop through the indicator array and create rectangles
  for (i in 1:n_rows) {
    for (j in 1:n_cols) {
      xmin <- lower_bds[1] + (i - 1) * width
      xmax <- xmin + width
      ymin <- lower_bds[2] + (j - 1) * height
      ymax <- ymin + height
      fill <- ifelse(indicator_array[i, j] == 1, "blue", "white")

      rectangles <- rbind(rectangles, data.frame(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        fill = fill
      ))
    }
  }

  # Plot the rectangles using ggplot2
  ggplot() +
    geom_rect(data = rectangles, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), color = "black") +
    scale_fill_identity() +
    theme_minimal() +
    labs(x = "parameter 1", y = "parameter 2") +
    theme(aspect.ratio = 1)
}
