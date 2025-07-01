# Validating Posterior Samplers: A Comprehensive Study

This project implements and evaluates validation methods for Bayesian posterior samplers, using documented cases of sampler errors from published research.

## Project Structure

```
├── README.md                    # This file
├── PROPOSAL.txt                # Original research proposal
├── R/                          # R implementation files
│   ├── validation_methods/     # Geweke tests, SBC, etc.
│   ├── error_cases/           # Implementations of documented errors
│   ├── utils/                 # Helper functions
│   └── analysis/              # Comparative analysis scripts
├── data/                      # Example datasets
├── results/                   # Output from validation comparisons
└── docs/                      # Documentation and reports
```

## Documented Error Cases

1. **Carriero, Clark, and Marcellino (2019)**: Triangular algorithm error in large Bayesian VARs
2. **Primiceri (2005)**: MCMC step ordering error in time-varying parameter models  
3. **Geweke (1994)**: Constraint handling error in unit root testing
4. **Sullivan & Greenland (2013)**: SAS implementation efficiency failures
5. **Additional cases**: Harmonic mean estimator failure, label switching

## Validation Methods

- **Geweke (2004)**: Joint distribution tests comparing marginal-conditional vs successive-conditional simulators
- **Simulation-Based Calibration (SBC)**: Modern validation approach from Stan community
- **Convergence diagnostics**: Gelman-Rubin, effective sample size, autocorrelation
- **Cross-validation approaches**: Alternative algorithms and analytical comparisons

## Research Questions

1. Which validation method is most effective at detecting different types of errors?
2. How do validation methods perform on large, complex models with dependent data?
3. What should practitioners know about implementing these validation approaches?

## Status

✅ **Project Status: Core Framework Complete and Operational**

- **Geweke joint distribution tests**: ✅ Implemented and validated
- **Simulation-Based Calibration (SBC)**: ✅ Implemented and validated  
- **Normal model error case**: ✅ Complete (correct vs incorrect Gibbs samplers)
- **VAR-SV triangular algorithm error**: ✅ Complete (Carriero et al. 2019 case)
- **Comprehensive testing framework**: ✅ Operational with automated demos
- **Documentation and examples**: ✅ Complete and ready for research use

## Getting Started

```r
# Install required packages and setup environment
source("R/setup.R")

# Test all implementations
source("R/analysis/comprehensive_demo.R")
test_all_implementations()

# Run full validation demo
results <- run_comprehensive_demo(quick_demo = TRUE)

# Run individual validation methods
source("R/examples/geweke_example.R")
geweke_results <- run_geweke_normal_example()
```

## Key Results

Our validation framework successfully:
- **Detects sampler errors**: Both methods identify algorithmic mistakes in posterior samplers
- **Works on complex models**: VAR with stochastic volatility (200+ parameters) validated
- **Provides actionable diagnostics**: Clear statistical tests and visualizations for error detection
- **Enables reproducible research**: All documented error cases implemented with exact specifications 