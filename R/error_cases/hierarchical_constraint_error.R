# Hierarchical Model Constraint Error Case
# Documented pattern: Incorrect handling of positive constraints in hierarchical models
# Based on common Stan forum posts about non-positive variance parameter errors

#' Simulate Hierarchical Normal Data
#' 
#' Standard hierarchical model with group-level effects
#' y_ij ~ N(theta_j, sigma^2)
#' theta_j ~ N(mu, tau^2)
#' 
#' @param J Number of groups
#' @param n_per_group Observations per group
#' @param mu Population mean
#' @param tau Between-group standard deviation
#' @param sigma Within-group standard deviation
#' @return List with data and true parameters
simulate_hierarchical_data <- function(J = 8, n_per_group = 10, mu = 0, tau = 1, sigma = 1) {
  
  # Group means
  theta <- rnorm(J, mean = mu, sd = tau)
  
  # Observations
  y <- list()
  group <- c()
  for (j in 1:J) {
    y_j <- rnorm(n_per_group, mean = theta[j], sd = sigma)
    y[[j]] <- y_j
    group <- c(group, rep(j, n_per_group))
  }
  
  y_combined <- unlist(y)
  
  list(
    y = y_combined,
    group = group,
    y_by_group = y,
    J = J,
    n_per_group = n_per_group,
    theta_true = theta,
    mu_true = mu,
    tau_true = tau,
    sigma_true = sigma
  )
}

#' Incorrect Hierarchical Sampler (Constraint Handling Error)
#' 
#' Common error: Allowing variance parameters to become negative or zero
#' This leads to sampling failures and incorrect posteriors
#' 
#' @param ndraws Number of MCMC draws
#' @param data Output from simulate_hierarchical_data()
#' @param init Initial values
#' @param proposal_sd Proposal standard deviations
#' @return Matrix of posterior draws
incorrect_hierarchical_sampler <- function(ndraws, data, 
                                         init = list(mu = 0, log_tau = 0, log_sigma = 0),
                                         proposal_sd = list(mu = 0.2, log_tau = 0.1, log_sigma = 0.1)) {
  
  y <- data$y
  group <- data$group
  J <- data$J
  n <- length(y)
  
  draws <- matrix(NA, ndraws, J + 3)  # theta_1, ..., theta_J, mu, tau, sigma
  colnames(draws) <- c(paste0("theta_", 1:J), "mu", "tau", "sigma")
  
  # Initialize
  current_mu <- init$mu
  current_log_tau <- init$log_tau
  current_log_sigma <- init$log_sigma
  current_theta <- aggregate(y, by = list(group), FUN = mean)$x  # Start with group means
  
  n_accepted <- list(mu = 0, log_tau = 0, log_sigma = 0, theta = rep(0, J))
  
  for (i in 1:ndraws) {
    
    # Update mu
    proposal_mu <- current_mu + rnorm(1, 0, proposal_sd$mu)
    # Accept unconditionally for simplicity (flat prior)
    current_mu <- proposal_mu
    n_accepted$mu <- n_accepted$mu + 1
    
    # Update log_tau (INCORRECT: Can still sample negative tau due to transformation error)
    proposal_log_tau <- current_log_tau + rnorm(1, 0, proposal_sd$log_tau)
    
    # BUG: Incorrect handling - sometimes use log_tau directly as tau
    if (runif(1) < 0.1) {  # 10% of the time, make the constraint error
      proposed_tau <- proposal_log_tau  # ERROR: Using log value as tau (can be negative!)
    } else {
      proposed_tau <- exp(proposal_log_tau)  # Correct transformation
    }
    
    if (proposed_tau > 0) {  # Only check positivity, but bug can still create negative values
      current_log_tau <- proposal_log_tau
      n_accepted$log_tau <- n_accepted$log_tau + 1
    }
    
    # Update log_sigma (similar constraint error)
    proposal_log_sigma <- current_log_sigma + rnorm(1, 0, proposal_sd$log_sigma)
    
    # BUG: Same transformation error as tau
    if (runif(1) < 0.1) {
      proposed_sigma <- proposal_log_sigma  # ERROR: Using log value as sigma
    } else {
      proposed_sigma <- exp(proposal_log_sigma)
    }
    
    if (proposed_sigma > 0) {
      current_log_sigma <- proposal_log_sigma
      n_accepted$log_sigma <- n_accepted$log_sigma + 1
    }
    
    # Update theta_j (can fail when tau or sigma are incorrectly negative)
    current_tau <- exp(current_log_tau)
    current_sigma <- exp(current_log_sigma)
    
    # This can fail with negative variance parameters
    for (j in 1:J) {
      y_j <- y[group == j]
      n_j <- length(y_j)
      
      # Posterior for theta_j | data, mu, tau, sigma
      precision_prior <- 1 / (current_tau^2 + 1e-8)  # Add small constant to prevent division by zero
      precision_data <- n_j / (current_sigma^2 + 1e-8)
      
      posterior_precision <- precision_prior + precision_data
      posterior_mean <- (precision_prior * current_mu + precision_data * mean(y_j)) / posterior_precision
      posterior_sd <- 1 / sqrt(posterior_precision)
      
      # This can sample extreme values when variances are wrong
      current_theta[j] <- rnorm(1, posterior_mean, posterior_sd)
      n_accepted$theta[j] <- n_accepted$theta[j] + 1
    }
    
    draws[i, ] <- c(current_theta, current_mu, current_tau, current_sigma)
  }
  
  attr(draws, "acceptance_rates") <- list(
    mu = n_accepted$mu / ndraws,
    tau = n_accepted$log_tau / ndraws,
    sigma = n_accepted$log_sigma / ndraws,
    theta = n_accepted$theta / ndraws
  )
  
  return(draws)
}

#' Correct Hierarchical Sampler (Proper Constraint Handling)
#' 
#' Implements proper transformations and constraint handling for variance parameters
#' 
#' @param ndraws Number of MCMC draws
#' @param data Output from simulate_hierarchical_data()
#' @param init Initial values
#' @param proposal_sd Proposal standard deviations
#' @return Matrix of posterior draws
correct_hierarchical_sampler <- function(ndraws, data,
                                       init = list(mu = 0, log_tau = 0, log_sigma = 0),
                                       proposal_sd = list(mu = 0.2, log_tau = 0.1, log_sigma = 0.1)) {
  
  y <- data$y
  group <- data$group
  J <- data$J
  n <- length(y)
  
  draws <- matrix(NA, ndraws, J + 3)
  colnames(draws) <- c(paste0("theta_", 1:J), "mu", "tau", "sigma")
  
  # Initialize
  current_mu <- init$mu
  current_log_tau <- init$log_tau
  current_log_sigma <- init$log_sigma
  current_theta <- aggregate(y, by = list(group), FUN = mean)$x
  
  n_accepted <- list(mu = 0, log_tau = 0, log_sigma = 0, theta = rep(0, J))
  
  for (i in 1:ndraws) {
    
    # Update mu
    proposal_mu <- current_mu + rnorm(1, 0, proposal_sd$mu)
    current_mu <- proposal_mu
    n_accepted$mu <- n_accepted$mu + 1
    
    # Update log_tau (CORRECT: Always use proper transformations)
    proposal_log_tau <- current_log_tau + rnorm(1, 0, proposal_sd$log_tau)
    
    # CORRECT: Always transform properly and add bounds checking
    proposed_tau <- exp(proposal_log_tau)
    
    # Add reasonable bounds to prevent extreme values
    if (proposed_tau > 1e-6 && proposed_tau < 100) {
      current_log_tau <- proposal_log_tau
      n_accepted$log_tau <- n_accepted$log_tau + 1
    }
    
    # Update log_sigma (CORRECT version)
    proposal_log_sigma <- current_log_sigma + rnorm(1, 0, proposal_sd$log_sigma)
    proposed_sigma <- exp(proposal_log_sigma)
    
    if (proposed_sigma > 1e-6 && proposed_sigma < 100) {
      current_log_sigma <- proposal_log_sigma
      n_accepted$log_sigma <- n_accepted$log_sigma + 1
    }
    
    # Update theta_j (robust to parameter values)
    current_tau <- exp(current_log_tau)
    current_sigma <- exp(current_log_sigma)
    
    for (j in 1:J) {
      y_j <- y[group == j]
      n_j <- length(y_j)
      
      # Robust posterior computation
      precision_prior <- 1 / current_tau^2
      precision_data <- n_j / current_sigma^2
      
      posterior_precision <- precision_prior + precision_data
      posterior_mean <- (precision_prior * current_mu + precision_data * mean(y_j)) / posterior_precision
      posterior_sd <- 1 / sqrt(posterior_precision)
      
      # Bounds checking for extreme posterior values
      posterior_sd <- min(posterior_sd, 10)  # Prevent extreme values
      
      current_theta[j] <- rnorm(1, posterior_mean, posterior_sd)
      n_accepted$theta[j] <- n_accepted$theta[j] + 1
    }
    
    draws[i, ] <- c(current_theta, current_mu, current_tau, current_sigma)
  }
  
  attr(draws, "acceptance_rates") <- list(
    mu = n_accepted$mu / ndraws,
    tau = n_accepted$log_tau / ndraws,
    sigma = n_accepted$log_sigma / ndraws,
    theta = n_accepted$theta / ndraws
  )
  
  return(draws)
}

#' Create validation samplers for hierarchical constraint error
#' 
#' @return List with correct and incorrect sampler functions
create_hierarchical_constraint_samplers <- function() {
  
  correct_sampler <- function(data, n_draws = 1000) {
    return(correct_hierarchical_sampler(n_draws, data))
  }
  
  incorrect_sampler <- function(data, n_draws = 1000) {
    return(incorrect_hierarchical_sampler(n_draws, data))
  }
  
  list(
    correct = correct_sampler,
    incorrect = incorrect_sampler,
    error_description = "Hierarchical constraint error: Improper handling of positive variance constraints"
  )
}

#' Test hierarchical constraint handling errors
#' 
#' @param J Number of groups
#' @param n_per_group Observations per group
#' @param n_posterior Number of posterior draws
#' @return Comparison of correct vs incorrect constraint handling
test_hierarchical_constraint_error <- function(J = 8, n_per_group = 10, n_posterior = 1000) {
  
  cat("=== Testing Hierarchical Constraint Handling Error ===\n")
  cat("Common MCMC error: Improper handling of positive variance constraints\n\n")
  
  # Generate test data
  data <- simulate_hierarchical_data(J, n_per_group, mu = 2, tau = 1.5, sigma = 0.8)
  
  cat("True parameters:\n")
  cat("  mu:", data$mu_true, "\n")
  cat("  tau:", data$tau_true, "\n")
  cat("  sigma:", data$sigma_true, "\n")
  cat("  Groups:", J, "with", n_per_group, "obs each\n\n")
  
  # Create samplers
  samplers <- create_hierarchical_constraint_samplers()
  
  # Test both samplers
  cat("Testing correct sampler (proper constraint handling)...\n")
  correct_result <- samplers$correct(data, n_posterior)
  
  cat("Testing incorrect sampler (constraint bugs)...\n")
  incorrect_result <- tryCatch({
    suppressWarnings(samplers$incorrect(data, n_posterior))
  }, error = function(e) {
    cat("ERROR in incorrect sampler:", e$message, "\n")
    return(NULL)
  })
  
  # Analyze results
  correct_estimates <- list(
    mu = mean(correct_result[, "mu"]),
    tau = mean(correct_result[, "tau"]),
    sigma = mean(correct_result[, "sigma"])
  )
  
  cat("\nCorrect sampler estimates:\n")
  cat("  mu:", correct_estimates$mu, "\n")
  cat("  tau:", correct_estimates$tau, "\n")
  cat("  sigma:", correct_estimates$sigma, "\n")
  
  if (!is.null(incorrect_result)) {
    incorrect_estimates <- list(
      mu = mean(incorrect_result[, "mu"]),
      tau = mean(incorrect_result[, "tau"]),
      sigma = mean(incorrect_result[, "sigma"])
    )
    
    cat("\nIncorrect sampler estimates:\n")
    cat("  mu:", incorrect_estimates$mu, "\n")
    cat("  tau:", incorrect_estimates$tau, "\n")
    cat("  sigma:", incorrect_estimates$sigma, "\n")
    
    # Check for constraint violations
    negative_tau <- sum(incorrect_result[, "tau"] <= 0)
    negative_sigma <- sum(incorrect_result[, "sigma"] <= 0)
    
    cat("\nConstraint violations in incorrect sampler:\n")
    cat("  Negative/zero tau:", negative_tau, "out of", n_posterior, "\n")
    cat("  Negative/zero sigma:", negative_sigma, "out of", n_posterior, "\n")
    
  } else {
    cat("\nIncorrect sampler FAILED due to constraint handling errors\n")
  }
  
  cat("\nComparison to true values:\n")
  cat("  Parameter    True      Correct   Error\n")
  cat(sprintf("  mu         %6.3f    %6.3f    %6.3f\n", 
              data$mu_true, correct_estimates$mu, abs(correct_estimates$mu - data$mu_true)))
  cat(sprintf("  tau        %6.3f    %6.3f    %6.3f\n", 
              data$tau_true, correct_estimates$tau, abs(correct_estimates$tau - data$tau_true)))
  cat(sprintf("  sigma      %6.3f    %6.3f    %6.3f\n", 
              data$sigma_true, correct_estimates$sigma, abs(correct_estimates$sigma - data$sigma_true)))
  
  return(invisible(list(
    data = data,
    correct_result = correct_result,
    incorrect_result = incorrect_result
  )))
} 