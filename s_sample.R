s_sample <- function(u, theta) {
  # Extract the mean (location) and std. dev (scale) from theta
  mean <- theta[1]
  sd <- theta[2]
   
  # Determine the number of samples and dimensions from the seed matrix u
  # are determined from the `u` matrix.
  num_samples <- nrow(u)
  num_dimensions <- ncol(u)
   
  # Initialize a matrix to store the samples
  samples <- matrix(0, nrow = num_samples, ncol = num_dimensions)
   
  # Generate samples from a location-scale normal distribution for each dimension 
  # loop iterates over each dimension, generating set of random nums. from a normal distribution 
  # w/  specified mean and std. dev.for that dimension. 
  # These nums are then added to the `samples` matrix.
  for (i in  1 : num_samples) {
    samples[i,] <- qnorm(u[i,], mean = mean, sd = sd)
  }
   
  return(samples) # matrix of simulated samples returned
}




# location-scale normal distribution -> the T statistic could be computed as the sum of sqrd diffs between each simulated value and 
# the observed data, -> summed across all dimensions.

T_calc <- function(x, S, theta) {
  # Extract the mean (location) and standard deviation (scale) from theta
  mean <- theta[1]
  sd <- theta[2]
   
  # Calculate the T statistic for each row in S
  # Here we assume a simple T statistic based on the difference from the mean
  T_values <- apply(S,  1, function(x) {
    diff <- x - mean
    sum(diff * diff)
  })
   
  return(T_values)
}

# assume mahalanobis
