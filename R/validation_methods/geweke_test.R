# Geweke (2004) Joint Distribution Test Implementation
# Tests posterior samplers by comparing two simulators for the joint distribution

#' Joint Distribution Test for Posterior Samplers
#' 
#' Implements Geweke (2004) validation method comparing marginal-conditional 
#' vs successive-conditional simulators for the joint distribution p(y, theta)
#' 
#' @param posterior_sampler Function that generates posterior samples given data
#' @param prior_sampler Function that generates samples from the prior
#' @param data_simulator Function that simulates data given parameters
#' @param test_functions List of functions to apply to joint samples
#' @param M Number of simulation draws
#' @param seed Random seed for reproducibility
#' @return List containing test results and diagnostics
geweke_joint_test <- function(posterior_sampler, prior_sampler, data_simulator, 
                             test_functions, M = 10000, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  cat("Running Geweke joint distribution test with", M, "draws...\n")
  
  # Method 1: Marginal-conditional simulator
  # Draw parameters from prior, then simulate data conditional on parameters
  cat("Method 1: Marginal-conditional simulation...\n")
  marginal_conditional <- vector("list", M)
  
  for (m in 1:M) {
    if (m %% 1000 == 0) cat("  Draw", m, "of", M, "\n")
    
    # Sample from prior
    theta <- prior_sampler()
    
    # Simulate data conditional on parameters  
    y <- data_simulator(theta)
    
    marginal_conditional[[m]] <- list(theta = theta, y = y)
  }
  
  # Method 2: Successive-conditional simulator (Gibbs-style)
  # Alternate between posterior sampling and data simulation
  cat("Method 2: Successive-conditional simulation...\n")
  successive_conditional <- vector("list", M)
  
  # Initialize with draw from prior
  theta <- prior_sampler()
  y <- data_simulator(theta)
  
  for (m in 1:M) {
    if (m %% 1000 == 0) cat("  Draw", m, "of", M, "\n")
    
    # Sample from posterior given current data
    theta <- posterior_sampler(y)
    
    # Simulate new data given new parameters
    y <- data_simulator(theta)
    
    successive_conditional[[m]] <- list(theta = theta, y = y)
  }
  
  # Apply test functions and compare distributions
  cat("Applying test functions and comparing distributions...\n")
  test_results <- list()
  
  for (i in seq_along(test_functions)) {
    func_name <- names(test_functions)[i]
    test_func <- test_functions[[i]]
    
    cat("  Testing function:", func_name, "\n")
    
    # Apply test function to both simulation methods
    marginal_values <- sapply(marginal_conditional, function(x) test_func(x$theta, x$y))
    successive_values <- sapply(successive_conditional, function(x) test_func(x$theta, x$y))
    
    # Kolmogorov-Smirnov test for distributional equality
    ks_test <- ks.test(marginal_values, successive_values)
    
    # Quantile-quantile comparison
    qq_comparison <- list(
      marginal_quantiles = quantile(marginal_values, probs = seq(0, 1, 0.01)),
      successive_quantiles = quantile(successive_values, probs = seq(0, 1, 0.01))
    )
    
    # Moment comparisons
    moments <- list(
      marginal = list(
        mean = mean(marginal_values),
        var = var(marginal_values),
        skew = moments::skewness(marginal_values),
        kurt = moments::kurtosis(marginal_values)
      ),
      successive = list(
        mean = mean(successive_values),
        var = var(successive_values), 
        skew = moments::skewness(successive_values),
        kurt = moments::kurtosis(successive_values)
      )
    )
    
    test_results[[func_name]] <- list(
      ks_test = ks_test,
      qq_comparison = qq_comparison,
      moments = moments,
      values = list(marginal = marginal_values, successive = successive_values)
    )
  }
  
  return(list(
    test_results = test_results,
    marginal_conditional = marginal_conditional,
    successive_conditional = successive_conditional,
    summary = summarize_test_results(test_results)
  ))
}

#' Summarize Geweke Test Results
#' 
#' @param test_results Output from geweke_joint_test
#' @return Summary data frame
summarize_test_results <- function(test_results) {
  
  summary_df <- data.frame(
    test_function = character(),
    ks_statistic = numeric(),
    ks_p_value = numeric(),
    mean_diff = numeric(),
    var_diff = numeric(),
    max_qq_diff = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (func_name in names(test_results)) {
    result <- test_results[[func_name]]
    
    # Calculate differences
    mean_diff <- result$moments$marginal$mean - result$moments$successive$mean
    var_diff <- result$moments$marginal$var - result$moments$successive$var
    max_qq_diff <- max(abs(result$qq_comparison$marginal_quantiles - 
                          result$qq_comparison$successive_quantiles))
    
    summary_df <- rbind(summary_df, data.frame(
      test_function = func_name,
      ks_statistic = result$ks_test$statistic,
      ks_p_value = result$ks_test$p.value,
      mean_diff = mean_diff,
      var_diff = var_diff,
      max_qq_diff = max_qq_diff,
      stringsAsFactors = FALSE
    ))
  }
  
  return(summary_df)
}

#' Plot Geweke Test Results
#' 
#' @param geweke_results Output from geweke_joint_test
#' @param test_function_name Name of test function to plot
plot_geweke_results <- function(geweke_results, test_function_name) {
  
  if (!(test_function_name %in% names(geweke_results$test_results))) {
    stop("Test function not found in results")
  }
  
  result <- geweke_results$test_results[[test_function_name]]
  
  # Create Q-Q plot
  par(mfrow = c(1, 2))
  
  # Q-Q plot
  qqplot(result$values$marginal, result$values$successive,
         main = paste("Q-Q Plot:", test_function_name),
         xlab = "Marginal-conditional",
         ylab = "Successive-conditional",
         pch = 19, cex = 0.5)
  abline(0, 1, col = "red", lwd = 2)
  
  # Overlay histograms
  hist(result$values$marginal, breaks = 50, alpha = 0.5, 
       col = "blue", freq = FALSE,
       main = paste("Distribution Comparison:", test_function_name),
       xlab = "Value")
  hist(result$values$successive, breaks = 50, alpha = 0.5,
       col = "red", freq = FALSE, add = TRUE)
  legend("topright", c("Marginal-conditional", "Successive-conditional"),
         fill = c("blue", "red"), alpha = 0.5)
  
  par(mfrow = c(1, 1))
} 