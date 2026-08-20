#include <Rcpp.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// Linear interpolation corresponding to R quantile type = 7.
// For probability p and sample size n:
// h = 1 + (n - 1) * p
// j = floor(h)
// gamma = h - j
// Q(p) = (1-gamma) * x[j] + gamma * x[j+1]
//
// [[Rcpp::export]]
double quantile_type7(NumericVector x, double p) {

    int n = x.size();

    if (n == 0) {
        stop("Cannot calculate a quantile from an empty vector.");
    }

    if (p < 0.0 || p > 1.0) {
        stop("Probability must be between 0 and 1.");
    }

    if (n == 1) {
        return x[0];
    }

    // h = 1 + (n - 1)p
    double h = 1.0 + (n - 1.0) * p;

    // Convert from R's 1-based indexing to C++'s 0-based indexing
    double index = h - 1.0;

    int lower = static_cast<int>(std::floor(index));
    int upper = static_cast<int>(std::ceil(index));

    // Boundary cases
    if (lower < 0) {
        lower = 0;
    }

    if (upper >= n) {
        upper = n - 1;
    }

    // Interpolation fraction
    double gamma = index - lower;

    return (1.0 - gamma) * x[lower] +
           gamma * x[upper];
}

// The BFR quantile relationship can be written as:
//
// log((1-x)/x)
//     = log(alpha)
//       - (1/beta) * log(tan(pi*p/2))

// [[Rcpp::export]]
NumericVector estimateBFR_cpp(NumericVector x) {

    // Input validation

    if (x.size() < 2) {
        stop("At least two observations are required.");
    }

    // Check all observations
    for (int i = 0; i < x.size(); ++i) {

        if (NumericVector::is_na(x[i])) {
            stop("Missing values are not allowed.");
        }

        if (!R_finite(x[i])) {
            stop("All observations must be finite.");
        }

        if (x[i] <= 0.0 || x[i] >= 1.0) {
            stop("All observations must be strictly between 0 and 1.");
        }
    }

    // Copy and sort the data

    NumericVector sorted_x = clone(x);

    std::sort(sorted_x.begin(), sorted_x.end());

    // Calculate sample quartiles

    double q25 = quantile_type7(sorted_x, 0.25);
    double q75 = quantile_type7(sorted_x, 0.75);

    // Transform the sample quantiles

    //
    // y = log((1-x)/x)
    //
    // z = log(tan(pi*p/2))
    //
    // The theoretical relationship is:
    //
    // y = log(alpha) - (1/beta) * z
    //

    double y25 = std::log((1.0 - q25) / q25);
    double y75 = std::log((1.0 - q75) / q75);

    double z25 =
        std::log(std::tan(M_PI * 0.25 / 2.0));

    double z75 =
        std::log(std::tan(M_PI * 0.75 / 2.0));

    // Estimate slope

    double denominator = z75 - z25;

    if (std::abs(denominator) < 1e-14) {
        stop("Unable to estimate parameters: denominator is too small.");
    }

    double slope =
        (y75 - y25) / denominator;

    // From:
    //
    // slope = -1 / beta
    //
    // therefore:
    //
    // beta = -1 / slope

    if (std::abs(slope) < 1e-14) {
        stop("Unable to estimate beta: estimated slope is too close to zero.");
    }

    double beta_hat = -1.0 / slope;

    // Estimate alpha

    // From:
    //
    // y = log(alpha) - z/beta
    //
    // therefore:
    //
    // log(alpha) = y + z/beta

    double log_alpha =
        y25 + z25 / beta_hat;

    double alpha_hat =
        std::exp(log_alpha);

    // Parameter validation

    if (!R_finite(alpha_hat) || alpha_hat <= 0.0) {
        stop("Estimated alpha is not a valid positive value.");
    }

    if (!R_finite(beta_hat) || beta_hat <= 0.0) {
        stop("Estimated beta is not a valid positive value.");
    }

    // Return estimates

    return NumericVector::create(
        _["alpha"] = alpha_hat,
        _["beta"] = beta_hat,
        _["q25"] = q25,
        _["q75"] = q75
    );
}