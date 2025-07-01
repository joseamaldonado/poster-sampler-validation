# Example: Geweke Test for Normal Model Posterior Samplers
# Demonstrates validation of correct vs incorrect Gibbs samplers

# Load required functions
source("R/validation_methods/geweke_test.R")
source("R/error_cases/normal_model_example.R")

#' Run Geweke validation example
#' 
#' @param M Number of simulation draws (default 5000 for demo)
#' @param n Sample size for synthetic data
#' @return List with results for correct and incorrect samplers
run_geweke_normal_example <- function(M = 5000, n = 10) {
  
  cat("=== Geweke Validation Example: Normal Model ===\n\n")
  
  # Model hyperparameters
  mu0 <- 0
  tau0_sq <- 1  
  nu0 <- 2
  sigma0_sq <- 1
  
  # Create samplers and simulators
  prior_sampler <- function() normal_prior_sampler(mu0, tau0_sq, nu0, sigma0_sq)
  data_simulator <- function(theta) normal_data_simulator(theta, n)
  
  correct_sampler <- create_posterior_sampler("correct", mu0, tau0_sq, nu0, sigma0_sq)
  incorrect_sampler <- create_posterior_sampler("incorrect", mu0, tau0_sq, nu0, sigma0_sq)
  
  # Test functions
  test_functions <- normal_test_functions()
  
  # Run Geweke tests
  cat("Testing CORRECT sampler:\n")
  cat("=======================\n")
  correct_results <- geweke_joint_test(
    posterior_sampler = correct_sampler,
    prior_sampler = prior_sampler, 
    data_simulator = data_simulator,
    test_functions = test_functions,
    M = M,
    seed = 12345
  )
  
  cat("\nTesting INCORRECT sampler:\n")
  cat("=========================\n")
  incorrect_results <- geweke_joint_test(
    posterior_sampler = incorrect_sampler,
    prior_sampler = prior_sampler,
    data_simulator = data_simulator, 
    test_functions = test_functions,
    M = M,
    seed = 12345
  )
  
  # Display results
  cat("\n=== RESULTS SUMMARY ===\n")
  cat("\nCorrect sampler results:\n")
  print(correct_results$summary)
  
  cat("\nIncorrect sampler results:\n")
  print(incorrect_results$summary)
  
  # Interpretation
  cat("\n=== INTERPRETATION ===\n")
  cat("For a CORRECT sampler, we expect:\n")
  cat("- High p-values (> 0.05) for Kolmogorov-Smirnov tests\n")
  cat("- Small differences in means and variances\n") 
  cat("- Small maximum Q-Q plot deviations\n\n")
  
  cat("For an INCORRECT sampler, we expect:\n")
  cat("- Low p-values (< 0.05) indicating distributional differences\n")
  cat("- Larger differences in moments\n")
  cat("- Larger Q-Q plot deviations\n\n")
  
  # Flag potential issues
  flag_issues <- function(results, sampler_name) {
    summary_df <- results$summary
    issues <- summary_df$ks_p_value < 0.05
    
    if (any(issues)) {
      cat("WARNING: ", sampler_name, " sampler shows issues in test functions:\n")
      problem_tests <- summary_df$test_function[issues]
      cat("  ", paste(problem_tests, collapse = ", "), "\n")
    } else {
      cat("GOOD: ", sampler_name, " sampler passes all tests\n")
    }
  }
  
  flag_issues(correct_results, "Correct")
  flag_issues(incorrect_results, "Incorrect")
  
  return(list(
    correct = correct_results,
    incorrect = incorrect_results,
    hyperparameters = list(mu0 = mu0, tau0_sq = tau0_sq, nu0 = nu0, sigma0_sq = sigma0_sq),
    sample_size = n
  ))
}

#' Plot comparison of correct vs incorrect sampler results
#' 
#' @param results Output from run_geweke_normal_example
plot_geweke_comparison <- function(results) {
  
  # Plot results for key test functions
  test_funcs_to_plot <- c("mu", "omega", "mu_times_omega")
  
  par(mfrow = c(length(test_funcs_to_plot), 2), mar = c(4, 4, 3, 1))
  
  for (func_name in test_funcs_to_plot) {
    
    # Correct sampler
    correct_result <- results$correct$test_results[[func_name]]
    qqplot(correct_result$values$marginal, correct_result$values$successive,
           main = paste("Correct:", func_name),
           xlab = "Marginal-conditional", ylab = "Successive-conditional",
           pch = 19, cex = 0.3)
    abline(0, 1, col = "red", lwd = 2)
    
    # Incorrect sampler  
    incorrect_result <- results$incorrect$test_results[[func_name]]
    qqplot(incorrect_result$values$marginal, incorrect_result$values$successive,
           main = paste("Incorrect:", func_name), 
           xlab = "Marginal-conditional", ylab = "Successive-conditional",
           pch = 19, cex = 0.3)
    abline(0, 1, col = "red", lwd = 2)
  }
  
  par(mfrow = c(1, 1))
}

# Run the example if this script is executed directly
if (interactive() || !exists("skip_example")) {
  
  cat("Loading packages...\n")
  if (!require(moments, quietly = TRUE)) {
    install.packages("moments")
    library(moments)
  }
  
  cat("Running Geweke validation example...\n")
  example_results <- run_geweke_normal_example(M = 2000)  # Smaller M for demo
  
  cat("\nCreating plots...\n")
  plot_geweke_comparison(example_results)
  
  cat("\nExample complete! Check the Q-Q plots to see the difference.\n")
  cat("Red line shows perfect agreement - deviations indicate sampler errors.\n")
} 