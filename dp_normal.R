

# Testing for differentially-private normal (confidence intervals for both mu and sigma)

source("get_CI.R")

reps <- 1000 # total number of simulations # 1000
upper_clamp <- 3
lower_clamp <- 0
n <- 100 # sample size
repro_size <- 200 # number of synthetic samples
eps <- 1 # privacy guarantee parameter
alpha <- .05 # confidence level
tol <- 10^-2

# search region for mean/variance
lower_bds <- c(-5, 0.1)
upper_bds <- c(5,5)

population_mu <- 1 # population mean
population_sigma <- 1 # population sd

# a function that generates a matrix of seeds
seed_generator <- function() {
  # generate R*(n+2) i.i.d. standard normals
  seeds <- matrix(rnorm(repro_size * (n + 2)), nrow = repro_size, ncol = n + 2)
  
  return(seeds)
}

# a function that computes the private statistics given the seeds and an assumed parameter
s_sample <- function(seeds, theta) {
  # generate the raw data points
  raw_data <- theta[1] + theta[2]*seeds[, 1:n]
  
  # clamp the raw data
  clamped <- pmin(pmax(raw_data, lower_clamp), upper_clamp)
  
  # compute the private statistics
  s_mean <- apply(clamped, 1, mean) + (upper_clamp - lower_clamp) / (n * eps) * seeds[, n+1]
  s_var <- apply(clamped, 1, var) + (upper_clamp - lower_clamp)^2 / (n * eps) * seeds[, n+2]
  
  return(cbind(s_mean, s_var))
}

# generate an observed statistic given the true parameter
s_observed <- function(){
  raw <- population_mu + population_sigma * rnorm(n)
  clamped <- pmin(pmax(raw, lower_clamp), upper_clamp)
  
  # compute the private statistics
  s_mean <- mean(clamped) + (upper_clamp - lower_clamp) / (n * eps) * rnorm(1)
  s_var <- var(clamped) + (upper_clamp - lower_clamp)^2 / (n * eps) * rnorm(1)
  
  return(c(s_mean, s_var))
}



# running a simulation for coverage
simulation <- function() {
  # generate the seeds and the observed statistic
  seeds <- seed_generator()
  s_obs <- s_observed()
  
  # run the get_CI function
  mean_CI <- get_CI2(alpha, lower_bds, upper_bds, 1, seeds, s_sample, s_obs, tol)
  sigma_CI <- get_CI2(alpha, lower_bds, upper_bds, 2, seeds, s_sample, s_obs, tol)
  
  # result (does the CI contain the true parameter?)
  result1 <- FALSE
  result2 <- FALSE
  if (mean_CI[1] <= population_mu && 
      mean_CI[2] >= population_mu) {
    result1 <- TRUE
  }
  if (sigma_CI[1] <= population_sigma && 
      sigma_CI[2] >= population_sigma) {
    result2 <- TRUE
  }
  
  return(cbind(result1, result2))
}

# error instance
start_time <- proc.time()
set.seed(2)
simulation()
end_time <- proc.time()
print(end_time - start_time)






