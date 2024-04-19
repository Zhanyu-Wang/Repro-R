

source("get_CI.R")


# first identify the smallest hyperbox containing the confidence set
confidence_grid <- function(alpha, lower_bds, upper_bds, seeds, G, s_obs, tol, res, T_stat = ma_depth) {
  
  # resolution gives the number of divisions for each parameter, divide by res to get the grid width
  grid_width <- (upper_bds - lower_bds)/res
  
  # an array representing whether or not a grid contains an accepted value
  indicator_array <- array(0, dim = res)
  
  # j represents which of the parameters we're currectly working with
  sub_search <- function(j, left_index, right_index) {
    updated_lower <- lower_bds
    updated_lower[j] <- lower_bds[j] + (left_index - 1) * grid_width
    updated_upper <- upper_bds
    updated_upper[j] <- lower_bds[j] + right_index * grid_width
    
    # call the accept function to see if the box (left_index, right_index) contains a valid point
    return(accept(alpha, updated_lower, updated_upper, seeds, G, s_obs, NULL, T_stat))
  }
  
  # bisect_left/right subroutines for the jth parameter
  bisect_left <- function(j, left_index, right_index) {
    left <- left_index
    right <- right_index
    right_previous <- right_index
    
    while (left < right + 1) {
      if (sub_search(j, left, right)) {
        right_previous <- right
        right <- floor((left + right)/2)
      } else {
        left <- right + 1
        right <- floor((right + right_previous)/2)
      }
    }
    
    return(left)
  }
  
  bisect_right <- function(j, left_index, right_index) {
    left <- left_index
    right <- right_index
    left_previous <- left_index
    
    while (left < right + 1) {
      if (sub_search(j, left, right)) {
        left_previous <- left
        left <- floor((left + right)/2)
      } else {
        right <- left - 1
        left <- floor((left + left_previous)/2)
      }
    }
    
    return(right)
  }
  
  # get the number of dimensions for the parameter space
  k <- length(lower_bds) 
  
  # create a length k list to store the upper/lower indices
  box <- vector("list", k)
  
  # compute the boundaries for each parameter
  for (j in 1:k) {
    mid_point <- floor(res[j]/2)
    left_boundary <- bisect_left(j, 1, mid_point)
    right_boundary <- bisect_right(j, mid_point + 1, res[j])
    
    # the jth element of the "box" list is a vector containing the left/right index. 
    box[j] <- c(left_boundary, right_boundary)
  }
  
  # now just apply the "accept" function to each grid within the box
  
}







