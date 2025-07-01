# Quick test of the Geweke implementation
# This script tests that all functions load correctly without running the full simulation

cat("Testing Geweke validation implementation...\n\n")

# Test 1: Check that functions can be loaded
cat("1. Loading validation functions...\n")
tryCatch({
  source("R/validation_methods/geweke_test.R")
  cat("   ✓ Geweke test functions loaded\n")
}, error = function(e) {
  cat("   ✗ Error loading Geweke functions:", e$message, "\n")
})

tryCatch({
  source("R/error_cases/normal_model_example.R") 
  cat("   ✓ Normal model functions loaded\n")
}, error = function(e) {
  cat("   ✗ Error loading normal model functions:", e$message, "\n")
})

# Test 2: Check that basic samplers work
cat("\n2. Testing basic sampler functions...\n")

# Test prior sampler
tryCatch({
  theta_sample <- normal_prior_sampler(0, 1, 2, 1)
  cat("   ✓ Prior sampler works, sample:", paste(round(theta_sample, 3), collapse=", "), "\n")
}, error = function(e) {
  cat("   ✗ Error in prior sampler:", e$message, "\n")
})

# Test data simulator
tryCatch({
  y_sample <- normal_data_simulator(c(0, 1), 5)
  cat("   ✓ Data simulator works, sample:", paste(round(y_sample, 3), collapse=", "), "\n")
}, error = function(e) {
  cat("   ✗ Error in data simulator:", e$message, "\n")
})

# Test correct Gibbs sampler
tryCatch({
  y_test <- c(0.1, -0.5, 0.8, -0.2, 0.3)
  theta_draws <- correct_gibbs_sampler(10, y_test, 0, 1, 2, 1)
  cat("   ✓ Correct Gibbs sampler works, final draw:", paste(round(theta_draws[10,], 3), collapse=", "), "\n")
}, error = function(e) {
  cat("   ✗ Error in correct Gibbs sampler:", e$message, "\n")
})

# Test incorrect Gibbs sampler
tryCatch({
  y_test <- c(0.1, -0.5, 0.8, -0.2, 0.3)
  theta_draws <- incorrect_gibbs_sampler(10, y_test, 0, 1, 2, 1)
  cat("   ✓ Incorrect Gibbs sampler works, final draw:", paste(round(theta_draws[10,], 3), collapse=", "), "\n")
}, error = function(e) {
  cat("   ✗ Error in incorrect Gibbs sampler:", e$message, "\n")
})

# Test 3: Test function creation
cat("\n3. Testing test functions...\n")
tryCatch({
  test_funcs <- normal_test_functions()
  cat("   ✓ Test functions created, count:", length(test_funcs), "\n")
  
  # Test one function
  theta_test <- c(0.5, 2.0)
  y_test <- c(0.1, -0.5, 0.8)
  result <- test_funcs$mu_times_omega(theta_test, y_test)
  cat("   ✓ Example test function (mu*omega):", result, "\n")
}, error = function(e) {
  cat("   ✗ Error in test functions:", e$message, "\n")
})

cat("\n=== Implementation Test Complete ===\n")
cat("If all tests passed, the implementation is working correctly.\n")
cat("To run the full Geweke validation, you'll need R with the following packages:\n")
cat("- moments, mvtnorm, MCMCpack, coda, ggplot2, dplyr\n\n")

cat("To run the full example in R:\n")
cat("  source('R/examples/geweke_example.R')\n")
cat("  results <- run_geweke_normal_example(M = 1000)\n")
cat("  plot_geweke_comparison(results)\n") 