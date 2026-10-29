**Status:** Active / Submitted for DATA 501  
**Funding:** This project was developed as academic coursework for DATA 501 at Victoria University of Wellington. No external funding was received.

# BFRdist: The BFR Distribution for Bounded Measurements

**DATA501 Software Package**  
**Authors:** Daniel Luo and Jennifer Hanna  
**Supervisor:** Alejandro Frery  

The `BFRdist` package provides a robust computational environment for the novel BFR distribution, an  analytically tractable model designed for double-bounded data confined to the unit interval $(0, 1)$. The BFR distribution, named after the authors of its foundational paper, Avik Bhattacharya, Alejandro C. Frery, and Andrea A. Rey, narrows the gap between the limitations of trading flexibility for analytic convenience in traditional bounded models (Beta and Kumaraswamy distributions) and is highly adaptable for analyzing proportions, entropy, and remote sensing metrics like the Generalized Radar Vegetation Index (GRVI).
Reports are in the BFRDist/inst/text/reports directory.

## Key Features

*   **Core Distribution Suite:** Standard `d/p/q/r` functions (`dbfr`, `pbfr`, `qbfr`, `rbfr`) for probability density, cumulative distribution, quantiles, and random generation.
*   **High-Performance Estimation:** Fast parameter estimation utilizing a compiled C++ backend (`Rcpp`), accessible via `fit_bfr()` using Quantile Estimation (QE) or Maximum Likelihood Estimation (MLE). By utilizing compiled C++ code for these intensive mathematical iterations, the package avoids R's native memory overhead and significantly reduces computation time when working with large datasets.
*   **Statistical Utilities:** Robust, order statistic-based measures of asymmetry and tail heaviness (`galton_skewness()`, `moors_kurtosis()`).
*   **Object-Oriented Diagnostics:** Custom S3 methods (`print`, `summary`, `predict`, `plot`) for goodness-of-fit metrics and Pearson-like $SK_G \times KU_M$ diagnostic diagrams.

---

## Installation

### For Standard Users (Recommended)
You can install the development version of `BFRdist` directly from our GitHub repository using the `pak` package. Because the package sits in a subdirectory of the main repository, use the following command:

**Note for Windows/macOS users:** Because `BFRdist` compiles C++ code during installation, Windows users must have Rtools installed. macOS users may need to install Xcode Command Line Tools by running `xcode-select --install` in their terminal.

```R
# Install pak if you haven't already
# install.packages("pak")

# Install BFRdist from the repository subdirectory
options(pkg.build_vignettes = TRUE)
pak::pak("luodani272/DATA501/BFRdist")
```

### For Package Developers
If you want to modify the source code or build the package locally, clone the repository and use `devtools`:

```bash
# Clone the repository in your terminal
git clone [https://github.com/luodani272/DATA501.git](https://github.com/luodani272/DATA501.git)
```

```R
# Open the BFRdist project in RStudio, then run:
install.packages("devtools")
Rcpp::compileAttributes() # Re-compile C++ code if modified
devtools::load_all()      # Load the package into your session
devtools::document()      # Update roxygen2 documentation
```

---

## Quick Start Example

```R
library(BFRdist)

# Load the built-in multi-temporal agricultural dataset (Bhattacharya et al., 2026)
data(grvi_data)

# Extract early-season June 6th GRVI for Soybeans 
soy_june <- grvi_data$GRVI_06Jun[grvi_data$Crop == "Soybeans"]

# Filter strictly to (0, 1) to satisfy BFR constraints
soy_clean <- soy_june[soy_june > 0 & soy_june < 1]

# Fit the BFR distribution using Quantile Estimation (QE) or Maximum Likelihood Estimation (MLE)
bfr_model_qe <- fit_bfr(soy_clean, method = "qe")
bfr_model_mle <- fit_bfr(soy_clean, method = "mle")

# View a clean summary of the parameter estimates
print(bfr_model_qe)
print(bfr_model_mle)

# View standard errors and goodness-of-fit metrics (Log-Likelihood, AIC, BIC)
summary(bfr_model_qe)
summary(bfr_model_mle)

# Generate a diagnostic plot with confidence ellipses
plot(bfr_model_qe)

# Generate the fitted density curve over the empirical histogram
plot(bfr_model_qe, type = "density")
```

## Input and Output Specifications
*   **Input (`x`):** A numeric vector of observations. All values must be strictly bounded within the open interval $(0, 1)$. The functions will reject `NA`, `NaN`, `Inf`, and values exactly equal to $0$ or $1$.
*   **Output:** The `fit_bfr()` function returns a custom `bfr` S3 object containing:
    *   Estimated `alpha` and `beta` parameters.
    *   The raw data vector.
    *   The estimation method used (`qe` or `mme`).
    *   Goodness-of-fit metrics (Log-Likelihood, AIC, BIC).
    
## Computational Platform & Reproducibility 

To ensure full reproducibility and replicability of the models and stochastic simulations (such as the parametric bootstrap confidence ellipses), this package was developed and tested on the following computational platform:

*   **Operating System:** Windows 11 x64 (build 26200)
*   **Architecture:** x86_64-w64-mingw32 / x64
*   **R Version:** 4.6.1 (2026-06-24 ucrt)
*   **IDE:** RStudio 2026.10.0.233
*   **C++ Backend:** `Rcpp` version 1.1.2 
*   **Random Number Generation (RNG):** `Mersenne-Twister` algorithm (Inversion and Rejection methods). 

**Note on Stochastic Simulation:** Any functions relying on random generation (such as `rbfr()` or the bootstrap replications in `plot.bfr()`) require a user-defined seed (e.g., `set.seed(123)`) prior to execution to guarantee exactly reproducible outputs.

**Note on Open Data:** To comply with open data standards, the built-in `grvi_data` dataset originates from standard, machine-readable CSV files, avoiding proprietary formats like Excel.

**Note on File Structure** 
Following the reproducible research framework by Frery et al. (2020), this package centralizes raw project artifacts, datasets, and documentation within the `inst/` directory to facilitate comprehensive review and evaluation. Specifically:

- **`inst/extdata/csv/`**: Houses the raw, machine-readable CSV datasets.

- **`inst/figures/`**: Contains raw design artifacts and vector schematics (such as Draw.io source files).

- **`inst/images/`**: Stores supplementary visual assets.

- **`inst/text/`**: Contains auxiliary project reports and documentation.

## Unit Testing

This package includes a comprehensive suite of automated unit tests using the `testthat` framework to verify mathematical boundaries (e.g., ensuring CDF bounds are exactly 0 to 1) and confirm C++ estimator accuracy. 

To run the test suite locally, clone the repository and execute:

```R
devtools::test()
```

## Documentation

For a comprehensive theoretical overview and mathematical formulation, please read the package vignette:

```R
vignette("intro-to-bfr", package = "BFRdist")
```

## Known Limitations

These are expected behaviors based on the mathematical formulation of the BFR distribution:

* **Strict Unit Interval:** The BFR distribution is defined strictly on the open interval $(0, 1)$. Field data (like GRVI or proportions) containing exact `0` or `1` boundaries must be filtered or adjusted prior to fitting, otherwise `fit_bfr()` will throw an out-of-bounds error.
* **MLE Computation Time:** While the C++ quantile estimator is near-instantaneous, Maximum Likelihood Estimation (MLE) utilizes `stats::optim()` to iteratively traverse the log-likelihood surface. For extremely large multitemporal datasets, MLE will take noticeably longer to converge. This method is intentionally included to allow for direct baseline comparisons and for scenarios where traditional statistical precision is prioritized over execution speed.

## Bug Reports 
If you encounter any issues or have feature requests, please open an issue on our [GitHub Issues page](https://github.com/luodani272/DATA501/issues).

## Citing the Package
To cite the BFRdist package in publications, please use:

```R
citation("BFRdist")
```

## Licence
This package is released under the MIT License.

## References

* Bhattacharya, A., Frery, A. C., & Rey, A. A. (2026). A Statistical Framework for Bounded Measurements with Applications to SAR Data. *IEEE Transactions on Geoscience and Remote Sensing*.
* Bhattacharya, A., Frery, A. C., & Rey, A. A. (2026). GRVI data. Personal communication (October).
* Frery, A. C., Gomez, L., & Medeiros, A. C. (2020). A Badging System for Reproducibility and Replicability in Remote Sensing Research. *IEEE Journal of Selected Topics on Applied Earth Observations and Remote Sensing*, 13, 4988-4995. https://doi.org/10.1109/JSTARS.2020.3019418
