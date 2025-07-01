# Bessel Function Overflow Error Case  
# Documented pattern of computational errors in Stan/Boost math library
# Based on common Stan forum posts about Bessel function convergence failures

#' Skellam Distribution Implementation (Source of Bessel Overflow)
#' 
#' The Skellam distribution uses modified Bessel functions and is prone to
#' overflow errors when lambda parameters become large, causing Stan sampling failures
#' 
#' @param k Integer difference (data)
#' @param lambda1 First Poisson rate parameter  
#' @param lambda2 Second Poisson rate parameter
#' @param safe_mode If TRUE, use numerically stable implementation
#' @return Log probability density
skellam_lpmf <- function(k, lambda1, lambda2, safe_mode = FALSE) {
  
  if (!safe_mode) {
    # Unstable implementation prone to overflow
    # This is the version that causes Bessel function failures
    norm_term <- -(lambda1 + lambda2)
    power_term <- (k/2) * log(lambda1/lambda2)
    bessel_arg <- 2 * sqrt(lambda1 * lambda2)
    
    # This can overflow when lambda1, lambda2 become large
    bessel_val <- besselI(bessel_arg, nu = k, expon.scaled = FALSE)
    log_bessel <- log(bessel_val)
    
    return(norm_term + power_term + log_bessel)
    
  } else {
    # Numerically stable implementation
    # Uses log-space calculations and scaled Bessel functions
    norm_term <- -(lambda1 + lambda2)
    power_term <- (k/2) * log(lambda1/lambda2)
    bessel_arg <- 2 * sqrt(lambda1 * lambda2)
    
    # Use expon.scaled=TRUE to prevent overflow
    log_bessel_scaled <- log(besselI(bessel_arg, nu = k, expon.scaled = TRUE))
    bessel_correction <- -bessel_arg  # Correction for scaling
    
    return(norm_term + power_term + log_bessel_scaled + bessel_correction)
  }
}

#' Simulate Skellam Data (Difference of Poissons)
#' 
#' @param n Sample size
#' @param lambda1 First Poisson rate
#' @param lambda2 Second Poisson rate
#' @return Vector of integer differences
simulate_skellam_data <- function(n, lambda1, lambda2) {
  x1 <- rpois(n, lambda1)
  x2 <- rpois(n, lambda2) 
  return(x1 - x2)
}

#' Unstable Posterior Sampler for Skellam Model
#' 
#' This sampler is prone to Bessel function overflow errors when
#' parameter proposals become large during MCMC exploration
#' 
#' @param ndraws Number of MCMC draws
#' @param y Observed Skellam data
#' @param init_log_lambda Initial log-lambda values
#' @param proposal_sd Standard deviation for random walk proposals
#' @return Matrix of posterior draws (ndraws x 2)
unstable_skellam_sampler <- function(ndraws, y, init_log_lambda = c(0, 0), 
                                   proposal_sd = 0.5) {
  
  n <- length(y)
  draws <- matrix(NA, ndraws, 2)
  colnames(draws) <- c("log_lambda1", "log_lambda2")
  
  # Initialize
  current_log_lambda <- init_log_lambda
  current_lambda <- exp(current_log_lambda)
  
  # Compute initial log-likelihood (prone to overflow)
  current_loglik <- tryCatch({
    sum(sapply(y, function(k) skellam_lpmf(k, current_lambda[1], current_lambda[2], safe_mode = FALSE)))
  }, error = function(e) -Inf)
  
  n_accepted <- 0
  
  for (i in 1:ndraws) {
    # Propose new values
    proposal_log_lambda <- current_log_lambda + rnorm(2, 0, proposal_sd)
    proposal_lambda <- exp(proposal_log_lambda)
    
    # Compute proposed log-likelihood (THIS CAN FAIL WITH BESSEL OVERFLOW)
    proposal_loglik <- tryCatch({
      sum(sapply(y, function(k) skellam_lpmf(k, proposal_lambda[1], proposal_lambda[2], safe_mode = FALSE)))
    }, error = function(e) {
      # Simulate the documented Stan error behavior
      if (runif(1) < 0.1) {  # 10% chance of reporting the error
        warning("Bessel function overflow error: Series evaluation exceeded 1000000 iterations")
      }
      return(-Inf)
    })
    
    # Metropolis acceptance
    if (is.finite(proposal_loglik)) {
      log_ratio <- proposal_loglik - current_loglik
      if (log(runif(1)) < log_ratio) {
        current_log_lambda <- proposal_log_lambda
        current_lambda <- proposal_lambda
        current_loglik <- proposal_loglik
        n_accepted <- n_accepted + 1
      }
    }
    
    draws[i, ] <- current_log_lambda
  }
  
  attr(draws, "acceptance_rate") <- n_accepted / ndraws
  return(draws)
}

#' Stable Posterior Sampler for Skellam Model  
#' 
#' This sampler uses numerically stable Bessel function calculations
#' and parameter constraints to avoid overflow errors
#' 
#' @param ndraws Number of MCMC draws
#' @param y Observed Skellam data
#' @param init_log_lambda Initial log-lambda values
#' @param proposal_sd Standard deviation for random walk proposals
#' @param max_log_lambda Maximum allowed log-lambda (constraint)
#' @return Matrix of posterior draws (ndraws x 2)
stable_skellam_sampler <- function(ndraws, y, init_log_lambda = c(0, 0),
                                 proposal_sd = 0.3, max_log_lambda = 5) {
  
  n <- length(y)
  draws <- matrix(NA, ndraws, 2)
  colnames(draws) <- c("log_lambda1", "log_lambda2")
  
  # Initialize
  current_log_lambda <- init_log_lambda
  current_lambda <- exp(current_log_lambda)
  
  # Compute initial log-likelihood (stable version)
  current_loglik <- sum(sapply(y, function(k) {
    skellam_lpmf(k, current_lambda[1], current_lambda[2], safe_mode = TRUE)
  }))
  
  n_accepted <- 0
  
  for (i in 1:ndraws) {
    # Propose new values with constraint
    proposal_log_lambda <- current_log_lambda + rnorm(2, 0, proposal_sd)
    
    # Apply constraint to prevent extreme values
    proposal_log_lambda <- pmin(proposal_log_lambda, max_log_lambda)
    proposal_lambda <- exp(proposal_log_lambda)
    
    # Compute proposed log-likelihood (stable version)
    proposal_loglik <- sum(sapply(y, function(k) {
      skellam_lpmf(k, proposal_lambda[1], proposal_lambda[2], safe_mode = TRUE)
    }))
    
    # Metropolis acceptance
    log_ratio <- proposal_loglik - current_loglik
    if (log(runif(1)) < log_ratio) {
      current_log_lambda <- proposal_log_lambda
      current_lambda <- proposal_lambda
      current_loglik <- proposal_loglik
      n_accepted <- n_accepted + 1
    }
    
    draws[i, ] <- current_log_lambda
  }
  
  attr(draws, "acceptance_rate") <- n_accepted / ndraws
  return(draws)
}

#' Create validation samplers for Bessel overflow error case
#' 
#' @return List with stable and unstable sampler functions
create_bessel_overflow_samplers <- function() {
  
  unstable_sampler <- function(y, n_draws = 1000) {
    return(unstable_skellam_sampler(n_draws, y))
  }
  
  stable_sampler <- function(y, n_draws = 1000) {
    return(stable_skellam_sampler(n_draws, y))
  }
  
  list(
    unstable = unstable_sampler,
    stable = stable_sampler,
    error_description = "Bessel function overflow: Common Stan/Boost math computational error"
  )
}

#' Test Bessel function overflow error patterns
#' 
#' @param n_data Sample size for test data
#' @param n_posterior Number of posterior draws
#' @param lambda1_true True first lambda parameter
#' @param lambda2_true True second lambda parameter
#' @return Comparison of stable vs unstable samplers
test_bessel_overflow_error <- function(n_data = 100, n_posterior = 1000,
                                     lambda1_true = 2.0, lambda2_true = 1.5) {
  
  cat("=== Testing Bessel Function Overflow Error ===\n")
  cat("Documented computational error in Stan/Boost math library\n\n")
  
  # Generate test data
  y <- simulate_skellam_data(n_data, lambda1_true, lambda2_true)
  
  cat("True parameters:\n")
  cat("  lambda1:", lambda1_true, "(log:", log(lambda1_true), ")\n")
  cat("  lambda2:", lambda2_true, "(log:", log(lambda2_true), ")\n")
  cat("  Data range:", range(y), "\n\n")
  
  # Create samplers
  samplers <- create_bessel_overflow_samplers()
  
  # Test unstable sampler (prone to overflow)
  cat("Testing unstable sampler (prone to Bessel overflow)...\n")
  unstable_result <- tryCatch({
    suppressWarnings(samplers$unstable(y, n_posterior))
  }, error = function(e) {
    cat("ERROR in unstable sampler:", e$message, "\n")
    return(NULL)
  })
  
  # Test stable sampler
  cat("Testing stable sampler (numerically robust)...\n")
  stable_result <- samplers$stable(y, n_posterior)
  
  # Compare results
  if (!is.null(unstable_result)) {
    unstable_lambda <- exp(unstable_result)
    unstable_mean <- colMeans(unstable_lambda)
    unstable_acceptance <- attr(unstable_result, "acceptance_rate")
    
    cat("\nUnstable sampler results:\n")
    cat("  lambda1 estimate:", unstable_mean[1], "\n")
    cat("  lambda2 estimate:", unstable_mean[2], "\n") 
    cat("  Acceptance rate:", unstable_acceptance, "\n")
  } else {
    cat("\nUnstable sampler FAILED due to computational errors\n")
  }
  
  stable_lambda <- exp(stable_result)
  stable_mean <- colMeans(stable_lambda)
  stable_acceptance <- attr(stable_result, "acceptance_rate")
  
  cat("\nStable sampler results:\n")
  cat("  lambda1 estimate:", stable_mean[1], "\n")
  cat("  lambda2 estimate:", stable_mean[2], "\n")
  cat("  Acceptance rate:", stable_acceptance, "\n")
  
  cat("\nComparison to true values:\n")
  cat("  Parameter    True      Stable    Error\n")
  cat(sprintf("  lambda1    %6.3f    %6.3f    %6.3f\n", 
              lambda1_true, stable_mean[1], abs(stable_mean[1] - lambda1_true)))
  cat(sprintf("  lambda2    %6.3f    %6.3f    %6.3f\n", 
              lambda2_true, stable_mean[2], abs(stable_mean[2] - lambda2_true)))
  
  return(invisible(list(
    data = y,
    true_params = c(lambda1_true, lambda2_true),
    unstable_result = unstable_result,
    stable_result = stable_result
  )))
} 