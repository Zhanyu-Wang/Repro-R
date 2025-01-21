

#' confidence_grid
#'
#' Plots the confidence grid
#'
#' @param alpha Significance level.
#' @param lower_bds Vector containing the lower bounds for the search space.
#' @param upper_bds Vector containing the upper bounds for the search space
#' @param seeds Seeds containing all the randomness.
#' @param G Data generating process.
#' @param s_obs Observed statistic.
#' @param tol tolerance of the confidence interval.
#' @param res resolution of the grid
#' @param t_init Starting point for the initial search.
#' @param T_stat A distance metric, default to Mahalanobis distance.
#' @return An array
#' @examples
#' # example1
#' example1
#' @import ggplot2
#' @export


confidence_grid <- function(alpha, lower_bds, upper_bds, seeds, G, s_obs, tol, res, t_init = NULL, T_stat = ma_depth) {
  
  # for each parameter, run the get_CI function to obtain a confidence interval
  updated_lower_bds <- lower_bds
  updated_upper_bds <- upper_bds
  
  p <- length(lower_bds)
  for (j in 1:p) {
    confidence_interval <- get_CI(alpha, lower_bds, upper_bds, j, seeds, G, s_obs, tol, t_init, T_stat)
    updated_lower_bds[j] <- confidence_interval[1]
    updated_upper_bds[j] <- confidence_interval[2]
  }
  
  # resolution gives the number of divisions for each parameter, divide by res to get the grid width
  grid_width <- (updated_upper_bds - updated_lower_bds) / res
  
  # an array representing whether or not a grid contains an accepted value
  indicator_array <- array(rep(0, res ** p), dim = rep(res, p))
  
  # given the indices, this function search the corresponding cube and update the indicator array
  update_array <- function(indices, indicator_array) {
    cube_lower_bds <- updated_lower_bds + (indices - 1) * grid_width
    cube_upper_bds <- cube_lower_bds + grid_width
    mid_point <- (cube_lower_bds + cube_upper_bds) / 2
    
    if (p_value_function:::p_value(cube_lower_bds, cube_upper_bds, seeds, G, s_obs, mid_point, T_stat)$p_val > alpha) {
      indicator_array <- do.call('[<-', c(list(indicator_array), as.list(indices), 1))
    }
    
    return(indicator_array)
  }
  
  # function that recursively loops through all sub-divided cubes,
  loop_thru <- function(curr_pos, indicator_array, former_indices = c()) {
    if (curr_pos <= p) {
      for (i in 1:res) {
        new_indices <- c(former_indices, i)
        indicator_array <- loop_thru(curr_pos + 1, indicator_array, new_indices)
      }
    } else {
      indicator_array <- update_array(former_indices, indicator_array)
    }
    return(indicator_array)
  }
  
  indicator_array <- loop_thru(1, indicator_array)
  
  return(list(ind_array = indicator_array,
              updated_lower_bds = updated_lower_bds,
              updated_upper_bds = updated_upper_bds))
}


# this function computes the projection of the indicator_array onto 2d or 1d arrays
grid_projection <- function(indicator_array, index_set) {
  if (length(index_set) == 1 | length(index_set) == 2) {
    return(apply(indicator_array, index_set, function(x) any(x==1)) * 1)
  } else {
    print("A vector of length 1 or 2 is expected for the second input.")
  }
}


# given the grid_projection (1d or 2d) and the bounds, this functions plots the confidence set
plot_grid <- function(indicator_array, lower_bds, upper_bds) {
  # get the dimensions
  n_rows <- nrow(indicator_array)
  n_cols <- ncol(indicator_array)
  
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
    labs(x = "X", y = "Y") +
    coord_fixed(ratio = 1)
}

# example
arr <- array(c(0,1,1,1,0,1,0,1,0), dim = c(3,3))
lower_bds <- c(0, 2)
upper_bds <- c(3, 3)

plot_grid(arr, lower_bds, upper_bds)

