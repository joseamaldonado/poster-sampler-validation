# SBC-compatible Normal Model Implementation
# Adapts the normal model example for Simulation-Based Calibration testing

source("R/error_cases/normal_model_example.R")

#' Create SBC-compatible posterior sampler for normal model
#' 
#' @param sampler_type Either "correct" or "incorrect"
#' @param mu0 Prior mean for mu
#' @param tau0_sq Prior variance for mu
#' @param nu0 Prior degrees of freedom for omega
#' @param sigma0_sq Prior scale for omega
#' @return Function that takes data and number of draws, returns posterior samples
create_sbc_posterior_sampler <- function(sampler_type = "correct", mu0 = 0, tau0_sq = 1,
                                        nu0 = 2, sigma0_sq = 1) {
  
  function(y, n_draws = 100) {
    if (sampler_type == "correct") {
      return(correct_gibbs_sampler(n_draws, y, mu0, tau0_sq, nu0, sigma0_sq))
    } else {
      return(incorrect_gibbs_sampler(n_draws, y, mu0, tau0_sq, nu0, sigma0_sq))
    }
  }
}

#' Create SBC-compatible prior sampler for normal model
#' 
#' @param mu0 Prior mean for mu
#' @param tau0_sq Prior variance for mu
#' @param nu0 Prior degrees of freedom for omega
#' @param sigma0_sq Prior scale for omega
#' @return Function that generates prior samples
create_sbc_prior_sampler <- function(mu0 = 0, tau0_sq = 1, nu0 = 2, sigma0_sq = 1) {
  function() {
    normal_prior_sampler(mu0, tau0_sq, nu0, sigma0_sq)
  }
}

#' Create SBC-compatible data simulator for normal model
#' 
#' @param n Sample size
#' @return Function that simulates data given parameters
create_sbc_data_simulator <- function(n = 10) {
  function(theta) {
    mu <- theta[1]
    omega <- theta[2]
    sigma <- sqrt(1/omega)
    rnorm(n, mean = mu, sd = sigma)
  }
}

#' Run SBC validation for normal model samplers
#' 
#' @param sampler_type Either "correct" or "incorrect"
#' @param n_simulations Number of SBC rounds
#' @param n_posterior_draws Number of posterior draws per round
#' @param n_data Sample size for synthetic data
#' @return SBC test results
run_normal_sbc <- function(sampler_type = "correct", n_simulations = 100, 
                          n_posterior_draws = 100, n_data = 10) {
  
  cat("Running SBC validation for", sampler_type, "sampler...\\n")
  
  # Create samplers and simulator
  posterior_sampler <- create_sbc_posterior_sampler(sampler_type)
  prior_sampler <- create_sbc_prior_sampler()
  data_simulator <- create_sbc_data_simulator(n_data)
  
  # Load SBC functions
  source("R/validation_methods/sbc_test.R")
  
  # Run SBC test
  sbc_results <- sbc_test(
    posterior_sampler = posterior_sampler,
    prior_sampler = prior_sampler,
    data_simulator = data_simulator,
    n_simulations = n_simulations,
    n_posterior_draws = n_posterior_draws,
    parameter_names = c("mu", "omega")
  )
  
  return(sbc_results)
} 