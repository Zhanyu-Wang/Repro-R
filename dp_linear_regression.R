
# differentially-private simple linear regression

source("get_CI.R")

library(doSNOW)
library(foreach)
library(snow)

reps <- 100 # total number of simulations # 1000
rep_inverse <- 1/reps # one time division
clamp <- 2
n <- 100 # sample size
repro_size <- 200 # number of synthetic samples
eps <- 1 # privacy guarantee parameter
alpha <- .05 # confidence level
tol <- 10^-2

# search region
lower_bds <- c(-0.5, -1, 0.1, 0.5, 0.1)
upper_bds <- c(1.5, 0, 1, 2, 1)

# true parameters
beta1 <- 0 # can change to 0.1, 0.2, ..., 1 to see how the rejection probabilities change
true_beta0 <- -0.5 # true beta 0 
true_x_mean <- 0.5 # true mean for x
true_x_var <- 1 # true var for x
true_eps_var <- 0.25 # true var for noise
mu <- 1 # mu-GDP

# a function that generates a matrix of seeds
seed_generator <- function() {
  # generate R*(n+2) i.i.d. standard normals
  seeds <- matrix(rnorm(repro_size * (2 * n + 5)), nrow = repro_size, ncol = 2 * n + 5)
  
  return(seeds)
}

# a function that computes the private statistics given the seeds and an assumed parameter
s_sample <- function(seeds, theta) {
  # generate the raw data points
  raw_x <- theta[3] + sqrt(theta[4]) * seeds[, 1: n]
  raw_eps <- sqrt(theta[5]) * seeds[, (n + 1):(2 * n)]
  raw_y <- theta[2] + raw_x * theta[1] + raw_eps
  
  # clamped data
  clamped_x <- pmin(pmax(raw_x, -clamp), clamp)
  clamped_x2 <- pmin(pmax(raw_x^2, 0), clamp^2)
  clamped_y <- pmin(pmax(raw_y, -clamp), clamp)
  clamped_y2 <- pmin(pmax(raw_y^2, 0), clamp^2)
  clamped_xy <- pmin(pmax(raw_x * raw_y, -clamp^2), clamp^2)
  
  # compute the private statistics
  s_1 <- apply(clamped_x, 1, mean) + 2 * clamp / (mu * n / sqrt(5)) * seeds[, 2*n + 1]
  s_2 <- apply(clamped_x2, 1, mean) + clamp^2 / (mu * n / sqrt(5)) * seeds[, 2*n + 2]
  s_3 <- apply(clamped_y, 1, mean) + 2 * clamp / (mu * n / sqrt(5)) * seeds[, 2*n + 3]
  s_4 <- apply(clamped_xy, 1, mean) + 2 * clamp^2 / (mu * n / sqrt(5)) * seeds[, 2*n + 4]
  s_5 <- apply(clamped_y2, 1, mean) + clamp^2 / (mu * n / sqrt(5)) * seeds[, 2*n + 5]
  
  return(cbind(s_1, s_2, s_3, s_4, s_5))
}

# generate an observed statistic given the true parameter
s_observed <- function(beta1_value){
  seed_obs <- rnorm(5)
  # generate the raw data points
  raw_x <- true_x_mean + sqrt(true_x_var) * rnorm(n)
  raw_eps <- sqrt(true_eps_var) * rnorm(n)
  raw_y <- true_beta0 + raw_x * beta1_value + raw_eps
  
  # clamped data
  clamped_x <- pmin(pmax(raw_x, -clamp), clamp)
  clamped_x2 <- pmin(pmax(raw_x^2, 0), clamp^2)
  clamped_y <- pmin(pmax(raw_y, -clamp), clamp)
  clamped_y2 <- pmin(pmax(raw_y^2, 0), clamp^2)
  clamped_xy <- pmin(pmax(raw_x * raw_y, -clamp^2), clamp^2)
  
  # compute the private statistics
  s_1 <- mean(clamped_x) + 2 * clamp / (mu * n / sqrt(5)) * seed_obs[1]
  s_2 <- mean(clamped_x2) + clamp^2 / (mu * n / sqrt(5)) * seed_obs[2]
  s_3 <- mean(clamped_y) + 2 * clamp / (mu * n / sqrt(5)) * seed_obs[3]
  s_4 <- mean(clamped_xy) + 2 * clamp^2 / (mu * n / sqrt(5)) * seed_obs[4]
  s_5 <- mean(clamped_y2) + clamp^2 / (mu * n / sqrt(5)) * seed_obs[5]
  
  return(c(s_1, s_2, s_3, s_4, s_5))
}

# testing for beta_1
test <- function(beta1_test) {
  seeds <- seed_generator()
  s_obs <- s_observed(beta1_test)
  
  CI <- get_CI(alpha, lower_bds, upper_bds, 1, seeds, s_sample, s_obs, tol)
  
  reject_or_not <- TRUE
  if (CI[1] < 0 && CI[2] > 0) {
    reject_or_not <- FALSE
  }
  
  return(reject_or_not)
}

# simulation with parallelization
start_time <- proc.time()
numCores <- detectCores() # detect the number of nodes available
cl <- makeCluster(numCores) # create a cluster of all available nodes
registerDoSNOW(cl) # register the cluster

# parallel execution
results <- foreach(i = 1:8,
                   .combine = 'cbind',
                   .packages=c('ddalpha')
                   ) %dopar% {
                     set.seed(i)
                     test(0.0)
                     #c(test(0.0), test(0.1), test(0.2), test(0.4), test(0.6), test(0.8), test(1.0))
                   }

stopCluster(cl) # end the cluster
end_time <- proc.time()

results
print(end_time - start_time)

