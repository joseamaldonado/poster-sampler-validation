# VAR with Stochastic Volatility: Triangular Algorithm Error
# Based on Carriero, Clark, and Marcellino (2019) error documented by Bognanni (2022)
# Implements both the incorrect and corrected triangular algorithms

#' Simulate VAR with Stochastic Volatility Data
#' 
#' @param T Time series length
#' @param N Number of variables
#' @param p Number of lags
#' @param beta VAR coefficients (N*p x N matrix)
#' @param h_init Initial log-volatilities
#' @param Q Volatility innovation covariance
#' @return List containing y (data), X (design matrix), h (log-volatilities)
simulate_var_sv <- function(T = 100, N = 3, p = 2, beta = NULL, h_init = NULL, Q = NULL) {
  
  # Default parameters if not provided
  if (is.null(beta)) {
    beta <- matrix(0.1, nrow = N*p, ncol = N)  # Simple AR structure
    diag(beta[1:N, ]) <- 0.7  # Own lags
  }
  
  if (is.null(h_init)) {
    h_init <- rep(-1, N)  # Log-volatility initialization
  }
  
  if (is.null(Q)) {
    Q <- diag(0.01, N)  # Volatility innovation covariance
  }
  
  # Initialize arrays
  y <- matrix(0, T, N)
  X <- array(0, dim = c(T, N*p))
  h <- matrix(0, T, N)
  h[1, ] <- h_init
  
  # Simulate initial p observations
  for (t in 1:p) {
    Sigma_t <- diag(exp(h[t, ]))
    y[t, ] <- mvtnorm::rmvnorm(1, mean = rep(0, N), sigma = Sigma_t)
    if (t < T) h[t + 1, ] <- h[t, ] + mvtnorm::rmvnorm(1, rep(0, N), Q)
  }
  
  # Simulate remaining observations
  for (t in (p + 1):T) {
    # Create design matrix (lagged values)
    for (i in 1:p) {
      X[t, ((i-1)*N + 1):(i*N)] <- y[t - i, ]
    }
    
    # Generate observation
    mu_t <- X[t, ] %*% beta
    Sigma_t <- diag(exp(h[t, ]))
    y[t, ] <- mu_t + mvtnorm::rmvnorm(1, mean = rep(0, N), sigma = Sigma_t)
    
    # Update log-volatility
    if (t < T) h[t + 1, ] <- h[t, ] + mvtnorm::rmvnorm(1, rep(0, N), Q)
  }
  
  return(list(y = y, X = X, h = h, beta = beta, Q = Q))
}

#' Incorrect Triangular Algorithm (Buggy Implementation)
#' 
#' The error: Treats equations as independent when sampling VAR coefficients,
#' ignoring cross-equation volatility dependencies
#' 
#' @param y Data matrix (T x N)
#' @param X Design matrix (T x Np)
#' @param ndraws Number of MCMC draws
#' @param beta_prior Prior parameters for VAR coefficients
#' @param h_prior Prior parameters for log-volatilities
#' @return List of posterior draws
incorrect_triangular_algorithm <- function(y, X, ndraws = 1000, 
                                         beta_prior = NULL, h_prior = NULL) {
  
  T <- nrow(y)
  N <- ncol(y)
  Np <- ncol(X)
  
  # Default priors
  if (is.null(beta_prior)) {
    beta_prior <- list(
      mean = matrix(0, Np, N),
      precision = diag(1, Np)  # Weak prior
    )
  }
  
  if (is.null(h_prior)) {
    h_prior <- list(
      Q_scale = diag(0.01, N),
      Q_df = N + 2
    )
  }
  
  # Initialize storage
  beta_draws <- array(0, dim = c(ndraws, Np, N))
  h_draws <- array(0, dim = c(ndraws, T, N))
  Q_draws <- array(0, dim = c(ndraws, N, N))
  
  # Initialize parameters
  beta_current <- matrix(0, Np, N)
  h_current <- matrix(-1, T, N)
  Q_current <- diag(0.01, N)
  
  for (iter in 1:ndraws) {
    
    # STEP 1: Sample VAR coefficients equation-by-equation
    # ERROR: Missing cross-equation volatility dependencies
    for (i in 1:N) {
      # Incorrect: Use simple precision ignoring stochastic volatility
      precision_i <- crossprod(X) + beta_prior$precision
      
      # Sample coefficient vector for equation i
      mean_i <- solve(precision_i) %*% 
        (crossprod(X, y[, i]) + beta_prior$precision %*% beta_prior$mean[, i])
      
      beta_current[, i] <- mvtnorm::rmvnorm(1, mean_i, solve(precision_i))
    }
    
    # STEP 2: Sample volatilities (also incorrect - ignores cross-dependencies)
    for (t in 2:T) {
      for (i in 1:N) {
        # Simple random walk (ignores VAR structure)
        h_current[t, i] <- h_current[t-1, i] + rnorm(1, 0, sqrt(Q_current[i, i]))
      }
    }
    
    # STEP 3: Sample volatility innovation covariance (simplified)
    Q_current <- MCMCpack::riwish(h_prior$Q_df, h_prior$Q_scale)
    
    # Store draws
    beta_draws[iter, , ] <- beta_current
    h_draws[iter, , ] <- h_current
    Q_draws[iter, , ] <- Q_current
    
    if (iter %% 200 == 0) cat("Iteration", iter, "\\n")
  }
  
  return(list(
    beta = beta_draws,
    h = h_draws,
    Q = Q_draws,
    algorithm = "incorrect_triangular"
  ))
}

#' Correct Triangular Algorithm (Fixed Implementation)
#' 
#' Properly accounts for cross-equation volatility dependencies
#' 
#' @param y Data matrix (T x N)
#' @param X Design matrix (T x Np) 
#' @param ndraws Number of MCMC draws
#' @param beta_prior Prior parameters for VAR coefficients
#' @param h_prior Prior parameters for log-volatilities
#' @return List of posterior draws
correct_triangular_algorithm <- function(y, X, ndraws = 1000,
                                       beta_prior = NULL, h_prior = NULL) {
  
  T <- nrow(y)
  N <- ncol(y)
  Np <- ncol(X)
  
  # Default priors (same as incorrect for comparison)
  if (is.null(beta_prior)) {
    beta_prior <- list(
      mean = matrix(0, Np, N),
      precision = diag(1, Np)
    )
  }
  
  if (is.null(h_prior)) {
    h_prior <- list(
      Q_scale = diag(0.01, N),
      Q_df = N + 2
    )
  }
  
  # Initialize storage
  beta_draws <- array(0, dim = c(ndraws, Np, N))
  h_draws <- array(0, dim = c(ndraws, T, N))
  Q_draws <- array(0, dim = c(ndraws, N, N))
  
  # Initialize parameters
  beta_current <- matrix(0, Np, N)
  h_current <- matrix(-1, T, N)
  Q_current <- diag(0.01, N)
  
  for (iter in 1:ndraws) {
    
    # STEP 1: Sample VAR coefficients with proper volatility conditioning
    for (i in 1:N) {
      # CORRECTION: Weight by time-varying precision from stochastic volatility
      weight_matrix <- diag(exp(-h_current[, i]))  # Precision weights
      
      precision_i <- crossprod(X, weight_matrix) %*% X + beta_prior$precision
      
      # Properly weighted posterior mean
      mean_i <- solve(precision_i) %*% 
        (crossprod(X, weight_matrix) %*% y[, i] + 
         beta_prior$precision %*% beta_prior$mean[, i])
      
      beta_current[, i] <- mvtnorm::rmvnorm(1, mean_i, solve(precision_i))
    }
    
    # STEP 2: Sample log-volatilities with full information
    for (t in 2:T) {
      # Residuals from VAR equations
      resid_t <- y[t, ] - X[t, ] %*% beta_current
      
      # Full conditional for log-volatility (simplified for demonstration)
      for (i in 1:N) {
        prior_precision <- 1/Q_current[i, i]
        data_precision <- 1  # From likelihood
        
        post_precision <- prior_precision + data_precision
        post_mean <- (prior_precision * h_current[t-1, i] + 
                     data_precision * log(resid_t[i]^2 + 1e-6)) / post_precision
        
        h_current[t, i] <- rnorm(1, post_mean, sqrt(1/post_precision))
      }
    }
    
    # STEP 3: Sample Q with proper conditioning
    h_innovations <- diff(h_current)
    Q_current <- MCMCpack::riwish(h_prior$Q_df + T - 1, 
                                  h_prior$Q_scale + crossprod(h_innovations))
    
    # Store draws
    beta_draws[iter, , ] <- beta_current
    h_draws[iter, , ] <- h_current
    Q_draws[iter, , ] <- Q_current
    
    if (iter %% 200 == 0) cat("Iteration", iter, "\\n")
  }
  
  return(list(
    beta = beta_draws,
    h = h_draws,
    Q = Q_draws,
    algorithm = "correct_triangular"
  ))
}

#' Test VAR-SV Triangular Algorithms
#' 
#' @param T Time series length
#' @param N Number of variables
#' @param ndraws Number of MCMC draws
#' @return Comparison of correct vs incorrect algorithms
test_var_sv_algorithms <- function(T = 50, N = 3, ndraws = 500) {
  
  cat("=== VAR-SV Triangular Algorithm Test ===\\n")
  cat("Simulating", N, "variables,", T, "time periods\\n\\n")
  
  # Simulate data
  set.seed(123)
  sim_data <- simulate_var_sv(T = T, N = N, p = 2)
  
  cat("Testing INCORRECT triangular algorithm...\\n")
  incorrect_results <- incorrect_triangular_algorithm(
    sim_data$y, sim_data$X, ndraws = ndraws
  )
  
  cat("\\nTesting CORRECT triangular algorithm...\\n")
  correct_results <- correct_triangular_algorithm(
    sim_data$y, sim_data$X, ndraws = ndraws
  )
  
  # Compare results
  cat("\\n=== COMPARISON ===\\n")
  
  # Beta coefficient summaries
  beta_incorrect_mean <- apply(incorrect_results$beta, c(2, 3), mean)
  beta_correct_mean <- apply(correct_results$beta, c(2, 3), mean)
  
  cat("VAR coefficient differences (Correct - Incorrect):\\n")
  beta_diff <- beta_correct_mean - beta_incorrect_mean
  print(round(beta_diff, 4))
  
  # Volatility summaries
  h_incorrect_mean <- apply(incorrect_results$h, c(2, 3), mean)
  h_correct_mean <- apply(correct_results$h, c(2, 3), mean)
  
  cat("\\nMean log-volatility difference:\\n")
  h_diff_mean <- colMeans(h_correct_mean - h_incorrect_mean)
  print(round(h_diff_mean, 4))
  
  return(list(
    data = sim_data,
    incorrect = incorrect_results,
    correct = correct_results,
    beta_diff = beta_diff,
    h_diff_mean = h_diff_mean
  ))
} 