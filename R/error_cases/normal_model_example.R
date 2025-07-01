# Normal Model Example: Bayesian inference for mean and precision
# Based on Section 3 of the proposal document
# Implements both correct and incorrect posterior samplers for validation testing

#' Correct Gibbs Sampler for Normal Model
#' 
#' Model: y_i | mu, omega ~ N(mu, 1/omega) iid
#'        mu ~ N(mu0, tau0^2)  
#'        omega ~ Gamma(nu0/2, nu0*sigma0^2/2)
#' 
#' @param ndraws Number of MCMC draws
#' @param y Observed data vector
#' @param mu0 Prior mean for mu
#' @param tau0_sq Prior variance for mu
#' @param nu0 Prior degrees of freedom for omega
#' @param sigma0_sq Prior scale for omega
#' @param init Initial values c(mu, omega)
#' @return Matrix of posterior draws (ndraws x 2)
correct_gibbs_sampler <- function(ndraws, y, mu0, tau0_sq, nu0, sigma0_sq, 
                                 init = c(mean(y), 1/var(y))) {
  
  # Data summaries
  mean_y <- mean(y)
  var_y <- var(y)
  n <- length(y)
  
  # Preallocate storage and initialize
  THETA <- matrix(0, nrow = ndraws, ncol = 2)
  colnames(THETA) <- c("mu", "omega")
  theta <- init
  
  for (s in 1:ndraws) {
    
    # Draw mean given precision and data
    # Full conditional: mu | omega, y ~ N(mu_n, tau_n^2)
    mu_n <- (mu0/tau0_sq + n * mean_y * theta[2]) / (1/tau0_sq + n * theta[2])
    tau_n_sq <- 1 / (1/tau0_sq + n * theta[2])
    theta[1] <- rnorm(1, mu_n, sqrt(tau_n_sq))
    
    # Draw precision given mean and data  
    # Full conditional: omega | mu, y ~ Gamma(nu_n/2, nu_n*sigma_n^2/2)
    nu_n <- nu0 + n
    sigma_n_sq <- (nu0 * sigma0_sq + (n - 1) * var_y + n * (mean_y - theta[1])^2) / nu_n
    theta[2] <- rgamma(1, nu_n / 2, rate = nu_n * sigma_n_sq / 2)
    
    THETA[s, ] <- theta
  }
  
  return(THETA)
}

#' Incorrect Gibbs Sampler for Normal Model (Missing factor of n)
#' 
#' This version has the error mentioned in the proposal: 
#' forgot to multiply by n when calculating mu_n
#' 
#' @param ndraws Number of MCMC draws
#' @param y Observed data vector  
#' @param mu0 Prior mean for mu
#' @param tau0_sq Prior variance for mu
#' @param nu0 Prior degrees of freedom for omega
#' @param sigma0_sq Prior scale for omega
#' @param init Initial values c(mu, omega)
#' @return Matrix of posterior draws (ndraws x 2)
incorrect_gibbs_sampler <- function(ndraws, y, mu0, tau0_sq, nu0, sigma0_sq,
                                   init = c(mean(y), 1/var(y))) {
  
  # Data summaries
  mean_y <- mean(y)
  var_y <- var(y)
  n <- length(y)
  
  # Preallocate storage and initialize
  THETA <- matrix(0, nrow = ndraws, ncol = 2)
  colnames(THETA) <- c("mu", "omega")
  theta <- init
  
  for (s in 1:ndraws) {
    
    # Draw mean given precision and data
    # ERROR: Missing factor of n in numerator
    mu_n <- (mu0/tau0_sq + mean_y * theta[2]) / (1/tau0_sq + n * theta[2])  # BUG HERE!
    tau_n_sq <- 1 / (1/tau0_sq + n * theta[2])
    theta[1] <- rnorm(1, mu_n, sqrt(tau_n_sq))
    
    # Draw precision given mean and data (this part is correct)
    nu_n <- nu0 + n
    sigma_n_sq <- (nu0 * sigma0_sq + (n - 1) * var_y + n * (mean_y - theta[1])^2) / nu_n
    theta[2] <- rgamma(1, nu_n / 2, rate = nu_n * sigma_n_sq / 2)
    
    THETA[s, ] <- theta
  }
  
  return(THETA)
}

#' Prior sampler for Normal Model
#' 
#' @param mu0 Prior mean for mu
#' @param tau0_sq Prior variance for mu
#' @param nu0 Prior degrees of freedom for omega
#' @param sigma0_sq Prior scale for omega
#' @return Vector c(mu, omega) drawn from prior
normal_prior_sampler <- function(mu0 = 0, tau0_sq = 1, nu0 = 2, sigma0_sq = 1) {
  
  mu <- rnorm(1, mu0, sqrt(tau0_sq))
  omega <- rgamma(1, nu0/2, rate = nu0 * sigma0_sq / 2)
  
  return(c(mu, omega))
}

#' Data simulator for Normal Model
#' 
#' @param theta Parameter vector c(mu, omega)
#' @param n Sample size
#' @return Vector of n observations from N(mu, 1/omega)
normal_data_simulator <- function(theta, n = 10) {
  
  mu <- theta[1]
  omega <- theta[2]
  
  return(rnorm(n, mu, sqrt(1/omega)))
}

#' Create posterior sampler function for Geweke test
#' 
#' @param sampler_type Either "correct" or "incorrect"
#' @param mu0 Prior mean for mu
#' @param tau0_sq Prior variance for mu
#' @param nu0 Prior degrees of freedom for omega
#' @param sigma0_sq Prior scale for omega
#' @param ndraws Number of MCMC draws per call
#' @return Function that takes data y and returns single posterior draw
create_posterior_sampler <- function(sampler_type = "correct", mu0 = 0, tau0_sq = 1, 
                                    nu0 = 2, sigma0_sq = 1, ndraws = 1) {
  
  function(y) {
    if (sampler_type == "correct") {
      draws <- correct_gibbs_sampler(ndraws, y, mu0, tau0_sq, nu0, sigma0_sq)
    } else {
      draws <- incorrect_gibbs_sampler(ndraws, y, mu0, tau0_sq, nu0, sigma0_sq)
    }
    
    # Return last draw (single sample for Geweke test)
    return(draws[ndraws, ])
  }
}

#' Test functions for normal model validation
#' 
#' @return List of test functions for Geweke validation
normal_test_functions <- function() {
  
  list(
    # Test mu alone
    "mu" = function(theta, y) theta[1],
    
    # Test omega alone  
    "omega" = function(theta, y) theta[2],
    
    # Test product mu * omega (tests dependency)
    "mu_times_omega" = function(theta, y) theta[1] * theta[2],
    
    # Test sum of squared residuals
    "sum_sq_resid" = function(theta, y) sum((y - theta[1])^2),
    
    # Test data mean
    "data_mean" = function(theta, y) mean(y),
    
    # Test data variance
    "data_var" = function(theta, y) var(y),
    
    # Test interaction between parameter and data
    "mu_times_ybar" = function(theta, y) theta[1] * mean(y)
  )
} 