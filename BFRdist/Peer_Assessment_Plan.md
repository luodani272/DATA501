# Peer Assessment Plan: BFRdist Package

**Assessees:** Jennifer Hanna and Daniel Luo  
**Assessor:** Taruna  

## 1. Repository Link and Installation Instructions

The `BFRdist` package is hosted on GitHub and contains a frozen release ready for peer review. To install and load the package, please run the following commands in your R console using the `remotes` package:   

```r
# Install the remotes package if you do not already have it
install.packages("remotes")

# Install the BFRdist package directly from the GitHub repository
remotes::install_github("luodani272/DATA501")

# Load the package
library(BFRdist)
```

## 2. Overview of Implemented Functions

The `BFRdist` package provides computational and visual tools for the BFR distribution (a model for double-bounded data on the unit interval). At this stage, all major functionality, C++ integration, and S3 object-oriented methods are fully implemented:   

* **Core Distribution Suite:** `dbfr()` (density), `pbfr()` (cumulative probability), `qbfr()` (quantiles), and `rbfr()` (random generation).
* **Summary Statistics:** `galton_skewness()` and `moors_kurtosis()` calculate robust, order statistic-based measures of asymmetry and tail heaviness.
* **Parameter Estimation:** `fit_bfr()` estimates shape parameters from empirical data. It utilizes a fast C++ backend (`estimateBFR_cpp`) for initial estimations and supports both Maximum Likelihood Estimation (MLE) and the Method of Moments (MME).
* **S3 Object-Oriented Ecosystem:** `fit_bfr()` returns a custom `bfr` object. This object is supported by `summary.bfr()` (for standard errors and goodness-of-fit metrics), `predict.bfr()`, and `plot.bfr()`.   
* **Diagnostic Plotting:** The `plot.bfr()` method simulates parametric bootstrap replications to generate a diagnostic Pearson-like $SK_G \times KU_M$ diagram with bivariate normal confidence ellipses.

## 3. Package Testing Plan (How to Check the Functions)

Please follow these steps to verify the package structure, test the functions, and attempt to break the code.

### Step A: Verify Documentation & Vignettes
1. Run `?fit_bfr` and `?plot.bfr` to ensure the Roxygen2 help files load correctly in the RStudio help pane and that all arguments are properly documented.
2. Run `browseVignettes("BFRdist")` to read the compiled introductory vignette detailing the mathematical formulation and package workflows.   

### Step B: Verify Core Functions & Output
Run the following script to test the standard probability functions and confirm they return the expected numeric outputs:   

```r
# Set parameters
alpha_val <- 1.5
beta_val <- 2.0

# 1. Test Density and Cumulative Probability
d_val <- dbfr(0.5, alpha = alpha_val, beta = beta_val)
p_val <- pbfr(0.5, alpha = alpha_val, beta = beta_val)

# 2. Test Random Generation
set.seed(123)
synthetic_data <- rbfr(500, alpha = alpha_val, beta = beta_val)

# 3. Test Summary Statistics
skg <- galton_skewness(synthetic_data)
kum <- moors_kurtosis(synthetic_data)
```

### Step C: Verify S3 Methods and OOP Implementation
Test the primary modeling function and ensure the resulting S3 object triggers the correct generic methods:   

```r
# 1. Fit the synthetic data using MLE
model_fit <- fit_bfr(synthetic_data, method = "mle")

# 2. Test S3 Methods
summary(model_fit)
predict(model_fit, n.ahead = 5)

# 3. Test Diagnostic Plot (Should generate a SK_G x KU_M diagram)
plot(model_fit)
```

### Step D: Injecting Errors (Breaking the Functions)
To verify the robustness of the package, please attempt to inject the following errors and report whether the package catches them with appropriate warning/error messages:

1. **Out-of-Bounds Data:** The BFR distribution is strictly defined on the $(0,1)$ interval. Pass a vector containing negative numbers, numbers greater than 1, or `NA` values into `fit_bfr()` or `galton_skewness()`.
2. **Invalid Shape Parameters:** Pass $\alpha = -1$ or $\beta = 0$ into the `dbfr()` or `rbfr()` functions.
3. **Invalid Method Argument:** Run `fit_bfr(synthetic_data, method = "bayesian")` to test if `match.arg` correctly rejects unsupported estimation methods.