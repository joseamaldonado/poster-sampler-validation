# Simulation-Based Calibration (SBC) Implementation
# Based on Talts et al. (2020) and related Stan community work

#' Simulation-Based Calibration Test
#' 
#' Implements SBC validation by checking that posterior quantiles are uniformly distributed
#' when parameters are drawn from the prior and data is simulated accordingly.
#' 
#' @param posterior_sampler Function that returns posterior samples given data
#' @param prior_sampler Function that generates samples from the prior
#' @param data_simulator Function that simulates data given parameters
#' @param n_simulations Number of SBC simulation rounds
#' @param n_posterior_draws Number of posterior draws per simulation
#' @param parameter_names Names of parameters for labeling
#' @param seed Random seed for reproducibility
#' @return List containing SBC test results and diagnostics
sbc_test <- function(posterior_sampler, prior_sampler, data_simulator,
                    n_simulations = 1000, n_posterior_draws = 100,
                    parameter_names = NULL, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  cat("Running SBC test with", n_simulations, "simulations...\n")
  cat("Using", n_posterior_draws, "posterior draws per simulation\n\n")
  
  # Store results
  sbc_results <- list()
  
  for (sim in 1:n_simulations) {
    if (sim %% 100 == 0) cat("Simulation", sim, "of", n_simulations, "\n")
    
    # Step 1: Sample "true" parameters from prior
    theta_true <- prior_sampler()
    if (is.vector(theta_true)) {
      theta_true <- matrix(theta_true, nrow = 1)
    }
    
    # Step 2: Simulate data given true parameters
    y_sim <- data_simulator(theta_true)
    
    # Step 3: Sample from posterior given simulated data
    posterior_draws <- posterior_sampler(y_sim, n_posterior_draws)
    if (is.vector(posterior_draws)) {
      posterior_draws <- matrix(posterior_draws, nrow = 1)
    }
    
    # Step 4: Calculate rank statistics
    # For each parameter, count how many posterior draws are less than true value
    n_params <- ncol(posterior_draws)
    ranks <- numeric(n_params)
    
    for (p in 1:n_params) {
      ranks[p] <- sum(posterior_draws[, p] < theta_true[p])
    }
    
    sbc_results[[sim]] <- list(
      theta_true = theta_true,
      ranks = ranks,
      simulation = sim
    )
  }
  
  # Process results
  cat("Processing SBC results...\n")
  
  # Extract rank statistics for all parameters
  n_params <- length(sbc_results[[1]]$ranks)
  if (is.null(parameter_names)) {
    parameter_names <- paste0("param_", 1:n_params)
  }
  
  rank_matrix <- matrix(0, nrow = n_simulations, ncol = n_params)
  colnames(rank_matrix) <- parameter_names
  
  for (sim in 1:n_simulations) {
    rank_matrix[sim, ] <- sbc_results[[sim]]$ranks
  }
  
  # Analyze uniformity of ranks
  analysis_results <- list()
  
  for (p in 1:n_params) {
    param_name <- parameter_names[p]
    ranks_p <- rank_matrix[, p]
    
    # Expected: uniform distribution on {0, 1, ..., n_posterior_draws}
    # Bin the ranks and test for uniformity
    breaks <- seq(-0.5, n_posterior_draws + 0.5, by = 1)
    hist_result <- hist(ranks_p, breaks = breaks, plot = FALSE)
    observed_counts <- hist_result$counts
    
    # Chi-square test for uniformity
    expected_count <- n_simulations / (n_posterior_draws + 1)
    chi_sq_stat <- sum((observed_counts - expected_count)^2 / expected_count)
    chi_sq_p_value <- 1 - pchisq(chi_sq_stat, df = n_posterior_draws)
    
    # Kolmogorov-Smirnov test against uniform
    # Transform ranks to [0,1] scale
    ranks_scaled <- ranks_p / n_posterior_draws
    ks_test <- ks.test(ranks_scaled, "punif", 0, 1)
    
    # Summary statistics
    rank_mean <- mean(ranks_p)
    rank_sd <- sd(ranks_p)
    expected_mean <- n_posterior_draws / 2
    expected_sd <- sqrt(n_posterior_draws * (n_posterior_draws + 1) / 12)
    
    analysis_results[[param_name]] <- list(
      ranks = ranks_p,
      histogram = list(breaks = breaks, counts = observed_counts),
      chi_sq_test = list(statistic = chi_sq_stat, p_value = chi_sq_p_value),
      ks_test = ks_test,
      summary_stats = list(
        mean = rank_mean, sd = rank_sd,
        expected_mean = expected_mean, expected_sd = expected_sd,
        z_score_mean = (rank_mean - expected_mean) / (expected_sd / sqrt(n_simulations))
      )
    )
  }
  
  return(list(
    rank_matrix = rank_matrix,
    analysis = analysis_results,
    settings = list(
      n_simulations = n_simulations,
      n_posterior_draws = n_posterior_draws,
      parameter_names = parameter_names
    ),
    raw_results = sbc_results
  ))
}

#' Plot SBC Results
#' 
#' @param sbc_results Output from sbc_test
#' @param parameter_name Name of parameter to plot (if NULL, plots all)
plot_sbc_results <- function(sbc_results, parameter_name = NULL) {
  
  settings <- sbc_results$settings
  n_params <- length(settings$parameter_names)
  
  if (is.null(parameter_name)) {
    # Plot all parameters
    n_cols <- min(3, n_params)
    n_rows <- ceiling(n_params / n_cols)
    par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))
    
    for (param in settings$parameter_names) {
      plot_single_sbc(sbc_results, param)
    }
    
    par(mfrow = c(1, 1))
  } else {
    # Plot single parameter
    plot_single_sbc(sbc_results, parameter_name)
  }
}

#' Plot single parameter SBC results
plot_single_sbc <- function(sbc_results, parameter_name) {
  
  if (!(parameter_name %in% names(sbc_results$analysis))) {
    stop("Parameter not found in SBC results")
  }
  
  analysis <- sbc_results$analysis[[parameter_name]]
  settings <- sbc_results$settings
  
  # Histogram of ranks
  hist_data <- analysis$histogram
  expected_count <- settings$n_simulations / (settings$n_posterior_draws + 1)
  
  # Create histogram
  barplot(hist_data$counts, 
          names.arg = 0:settings$n_posterior_draws,
          main = paste("SBC Histogram:", parameter_name),
          xlab = "Rank", ylab = "Count",
          col = "lightblue", border = "black")
  
  # Add expected count line
  abline(h = expected_count, col = "red", lwd = 2, lty = 2)
  
  # Add test results as text
  chi_sq_p <- round(analysis$chi_sq_test$p_value, 4)
  ks_p <- round(analysis$ks_test$p.value, 4)
  
  legend("topright", 
         legend = c(paste("Chi-sq p =", chi_sq_p),
                   paste("KS p =", ks_p),
                   paste("Expected =", round(expected_count, 1))),
         col = c("black", "black", "red"),
         lty = c(NA, NA, 2),
         pch = c(NA, NA, NA),
         lwd = c(NA, NA, 2),
         bty = "n")
}

#' Summarize SBC Test Results
#' 
#' @param sbc_results Output from sbc_test
#' @return Data frame with summary statistics
summarize_sbc_results <- function(sbc_results) {
  
  parameter_names <- sbc_results$settings$parameter_names
  n_params <- length(parameter_names)
  
  summary_df <- data.frame(
    parameter = character(n_params),
    chi_sq_statistic = numeric(n_params),
    chi_sq_p_value = numeric(n_params),
    ks_statistic = numeric(n_params),
    ks_p_value = numeric(n_params),
    mean_rank = numeric(n_params),
    sd_rank = numeric(n_params),
    z_score_mean = numeric(n_params),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n_params) {
    param_name <- parameter_names[i]
    analysis <- sbc_results$analysis[[param_name]]
    
    summary_df[i, ] <- list(
      parameter = param_name,
      chi_sq_statistic = analysis$chi_sq_test$statistic,
      chi_sq_p_value = analysis$chi_sq_test$p_value,
      ks_statistic = analysis$ks_test$statistic,
      ks_p_value = analysis$ks_test$p.value,
      mean_rank = analysis$summary_stats$mean,
      sd_rank = analysis$summary_stats$sd,
      z_score_mean = analysis$summary_stats$z_score_mean
    )
  }
  
  return(summary_df)
} 