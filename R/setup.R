# Setup script for posterior sampler validation project
# Installs required packages and sets up environment

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Core packages for Bayesian computation
required_packages <- c(
  "mvtnorm",        # Multivariate normal distributions  
  "MCMCpack",       # MCMC algorithms
  "coda",           # MCMC diagnostics
  "bayesplot",      # Bayesian plotting
  "ggplot2",        # Plotting
  "dplyr",          # Data manipulation
  "MASS",           # Additional statistical functions
  "forecast",       # Time series (for VAR examples)
  "vars",           # Vector autoregressions
  "matrixStats",    # Matrix operations
  "parallel"        # Parallel computing
)

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all packages
lapply(required_packages, library, character.only = TRUE)

cat("Setup complete! All required packages loaded.\n")
cat("Project directory structure created.\n")
cat("Ready to run validation examples.\n") 