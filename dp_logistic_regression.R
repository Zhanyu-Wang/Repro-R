

# differentially-private logistic regression


source("get_CI.R")

library(doSNOW)
library(foreach)
library(snow)



# define the expit function used in logistc regression
expit <- function(x) {
  exp(x) / (1 + exp(x))
}



# sampling (algorithm 5)
sampling <- function(c) {
  U <- 2 * runif(2) - 1
  r <- rgamma(1, shape = 3, rate = c)
  return(r * U)
}



# true parameters
a_true <- 0.5
b_true <- 0.5
beta_0_true <- 0.5
beta_1_true <- 2
repro_size <- 200
alpha <- 0.05

# bounds for the parameters
lower_bds <- c(0.01, 0.01, -5, -5)
upper_bds <- c(0.99, 0.99, 5, 5)
tol <- 10^-2

# epsilon-DP
epsilon <- 1

# algorithm parameters
Delta_inf <- 2 # for objective perturbation
Delta_k <- 1
lambda <- 1/2
q <- 0.85
gamma <- lambda / (exp(epsilon * (1 - q)) - 1) # parameter in objective perturbation
c <- epsilon * q / Delta_inf

# sample size
n <- 200



# a function that determines if a point in R^2 is in predetermined convex hull
is_in_hull <- function(point) {
  is_in <- FALSE
  u1 <- point[1]
  u2 <- point[2]
  
  # First region: u1 <= -1/2
  if (u1 <= -1/2) {
    if ((u1 + 1)^2 - 1 <= u2 && u2 <= -u1^2) {
      is_in <- TRUE
    }
  }
  
  # Second region: -1/2 < u1 <= 1/2
  else if (-1/2 < u1 && u1 <= 1/2) {
    if (u1 - 1/4 <= u2 && u2 <= u1 + 1/4) {
      is_in <- TRUE
    }
  }
  
  # Third region: u1 > 1/2
  else if (u1 > 1/2) {
    if (u1^2 <= u2 && u2 <= 1 - (u1 - 1)^2) {
      is_in <- TRUE
    }
  }
  
  return(is_in)
}



# a function that generates a matrix of seeds
seed_generator <- function(repro_size) {

  # generate R*2n i.i.d. standard uniforms
  unif_seeds <- matrix(runif(repro_size * 2 * n), nrow = repro_size, ncol = 2*n)
  
  # generate K-norm random variables
  K_norm_seeds <- matrix(NA, nrow = repro_size, ncol = 2)
  for (i in 1:repro_size) {
    accepted <- FALSE
    r <- rgamma(1, shape = 3, epsilon / Delta_k)
    
    while (!accepted) {
      U <- (2 * runif(2) - 1) # Delta_inf = 1 for K_norm mechanism
      
      if (is_in_hull(U)) {
        accepted <- TRUE
      }
    }
    
    K_norm_seeds[i,] <- r * U
  }
  
  # The randomly sampled V as in algorithm 5
  obj_pert_seeds <- matrix(NA, nrow = repro_size, ncol = 2)
  for (i in 1:repro_size) {
    obj_pert_seeds[i, ] <- sampling(c)
  }
  
  seeds <- cbind(unif_seeds, K_norm_seeds, obj_pert_seeds)
  
  return(seeds)
}



# loss function in objective perturbation
loss <- function(theta_2half, x_data, y_data) {
  beta0 <- theta_2half[1]
  beta1 <- theta_2half[2]
  
  sum(log(1 + exp(beta0 + beta1 * x_data)) - y_data * (beta0 + beta1 * x_data)) / n
}



# s_sample function
s_sample <- function(seeds, theta) {
  # logistic regression raw data
  z <- qbeta(seeds[, 1:n], theta[1], theta[2])
  x <- 2 * z - 1
  y <- as.integer(expit(theta[3] + theta[4] * x) > seeds[, (n+1):(2*n)])
  
  # the first two privacy statistics
  s_1 <- apply(z, 1, sum) + seeds[, 2*n + 1]
  s_2 <- apply(z^2, 1, sum) + seeds[, 2*n + 2]
  
  # the other two privacy statistics come from objective perturbation
  objective <- function(theta_2half, V) {
    loss(theta_2half, x, y) + gamma / (2 * n) * sum(theta_2half^2) + sum(V * theta_2half) / n
  }
  
  s_34 <- matrix(NA, nrow = repro_size, ncol = 2)
  for (i in 1:repro_size) {
    result <- optim(par = c(1, 1), fn = objective, V = seeds[i, (2*n + 3):(2*n + 4)])
    s_34[i, ] <- result$par
  }
  
  s <- cbind(s_1, s_2, s_34)
  
  return(s)
}


s_observed <- function(integer_seed) {
  set.seed(integer_seed)
  seeds <- seed_generator(1)
  # logistic regression raw data
  z <- qbeta(seeds[1:n], a_true, b_true)
  x <- 2 * z - 1
  y <- as.integer(expit(beta_0_true + beta_1_true * x) > seeds[(n+1):(2*n)])
  
  # the first two privacy statistics
  s_1 <- sum(z) + seeds[2*n + 1]
  s_2 <- sum(z^2) + seeds[2*n + 2]
  
  # the other two privacy statistics come from objective perturbation
  objective <- function(theta_2half, V) {
    loss(theta_2half, x, y) + gamma / (2 * n) * sum(theta_2half^2) + sum(V * theta_2half) / n
  }

  s_34 <- optim(par = c(1, 1), fn = objective, V = seeds[(2*n + 3):(2*n + 4)])$par
  s <- c(s_1, s_2, s_34)
  
  return(s)
}


test_seeds <- seed_generator(repro_size)
s_obs <- s_observed(1)

test_CI1 <- get_CI(alpha, lower_bds, upper_bds, 1, test_seeds, s_sample, s_obs, tol)
test_CI2 <- get_CI(alpha, lower_bds, upper_bds, 2, test_seeds, s_sample, s_obs, tol)
test_CI3 <- get_CI(alpha, lower_bds, upper_bds, 3, test_seeds, s_sample, s_obs, tol)
test_CI4 <- get_CI(alpha, lower_bds, upper_bds, 4, test_seeds, s_sample, s_obs, tol)


test_CI1
test_CI2
test_CI3
test_CI4

# add hessian for faster optimization












