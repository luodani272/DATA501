# BFRdist: The BFR Distribution for Bounded Measurements

**DATA501 Software Package**  
**Authors:** Daniel Luo and Jennifer Hanna  
**Supervisor:** Alejandro Frery  

The `BFRdist` package provides a robust computational environment for the BFR distribution, an  analytically tractable model designed for double-bounded data confined to the unit interval $(0, 1)$. It overcomes the limitations of traditional bounded models (like the Beta and Kumaraswamy distributions) and is highly adaptable for analyzing proportions, entropy, and remote sensing metrics like the Generalized Radar Vegetation Index (GRVI).

## Key Features

*   **Core Distribution Suite:** Standard `d/p/q/r` functions (`dbfr`, `pbfr`, `qbfr`, `rbfr`) for probability density, cumulative distribution, quantiles, and random generation.
*   **High-Performance Estimation:** Fast parameter estimation utilizing a compiled C++ backend (`Rcpp`), accessible via `fit_bfr()` using Maximum Likelihood Estimation (MLE) or Method of Moments (MME). By utilizing compiled C++ code for these intensive mathematical iterations, the package avoids R's native memory overhead and significantly reduces computation time when working with large datasets.
*   **Statistical Utilities:** Robust, order statistic-based measures of asymmetry and tail heaviness (`galton_skewness()`, `moors_kurtosis()`).
*   **Object-Oriented Diagnostics:** Custom S3 methods (`print`, `summary`, `predict`, `plot`) for goodness-of-fit metrics and Pearson-like $SK_G \times KU_M$ diagnostic diagrams.

---

## Installation

### For Standard Users (Recommended)
You can install the development version of `BFRdist` directly from our GitHub repository using the `pak` package. Because the package sits in a subdirectory of the main repository, use the following command:

```R
# Install pak if you haven't already
# install.packages("pak")

# Install BFRdist from the repository subdirectory
pak::pak("luodani272/DATA501/BFRdist")
```

### For Package Developers
If you want to modify the source code or build the package locally, clone the repository and use `devtools`:

```bash
# 1. Clone the repository in your terminal
git clone [https://github.com/luodani272/DATA501.git](https://github.com/luodani272/DATA501.git)
```

```R
# 2. Open the BFRdist project in RStudio, then run:
install.packages("devtools")
Rcpp::compileAttributes() # Re-compile C++ code if modified
devtools::load_all()      # Load the package into your session
devtools::document()      # Update roxygen2 documentation
```

---

## Quick Start Example

```R
library(BFRdist)

# Load the built-in multi-temporal agricultural dataset
data(grvi_data)

# Extract July 24th GRVI measurements for Canola
canola_july <- grvi_data$GRVI_24_july[grvi_data$Crop == "Canola"]

# The BFR distribution strictly requires data on the open interval (0, 1).
# Filter out any exact 0.0 or 1.0 boundary values before fitting.
canola_clean <- canola_july[canola_july > 0 & canola_july < 1]

# Fit the BFR distribution using Maximum Likelihood Estimation (MLE)
bfr_model <- fit_bfr(canola_clean, method = "mle")

# View a clean summary of the parameter estimates
print(bfr_model)

# View standard errors and goodness-of-fit metrics (Log-Likelihood, AIC, BIC)
summary(bfr_model)

# Generate a diagnostic plot with confidence ellipses
plot(bfr_model)
```

## Input and Output Specifications
*   **Input (`x`):** A numeric vector of observations. All values must be strictly bounded within the open interval $(0, 1)$. The functions will reject `NA`, `NaN`, `Inf`, and values exactly equal to $0$ or $1$.
*   **Output:** The `fit_bfr()` function returns a custom `bfr` S3 object containing:
    *   Estimated `alpha` and `beta` parameters.
    *   The raw data vector.
    *   The estimation method used (`mle` or `mme`).
    *   Goodness-of-fit metrics (Log-Likelihood, AIC, BIC).

## Limitations
*   The package requires data strictly strictly within $(0, 1)$. For datasets containing exact $0$s or $1$s (e.g., zero-inflated proportions), a transformation or a zero-inflated mixture model is required before using `BFRdist`.
*   The Method of Moments (MME) estimator relies on sample quartiles, which may be sensitive to extreme sparsity in very small sample sizes. 

## Documentation

For a comprehensive theoretical overview and mathematical formulation, please read the package vignette:

```R
vignette("Introduction_to_BFRdist", package = "BFRdist")
```

## Licence
This package is released under the MIT License.