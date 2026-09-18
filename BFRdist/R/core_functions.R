#' The BFR Distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the BFR distribution with parameters `alpha` and `beta`.
#'
#' @details
#' The BFR distribution with shape parameters \eqn{\alpha > 0} and \eqn{\beta > 0} 
#' has the probability density function (PDF):
#' \deqn{f(x) = \frac{2\beta}{\pi} \frac{\alpha^\beta x^{\beta - 1} (1 - x)^{\beta - 1}}{\alpha^{2\beta} x^{2\beta} + (1 - x)^{2\beta}}}
#' for \eqn{0 < x < 1}.
#' 
#' The cumulative distribution function (CDF) is given by:
#' \deqn{F(q) = \frac{2}{\pi} \arctan\left( \frac{\alpha q}{1 - q} \right)^\beta}
#' 
#' The quantile function is the inverse of the CDF:
#' \deqn{Q(p) = \left( 1 + \alpha \left( \tan\left(\frac{\pi p}{2}\right) \right)^{-1/\beta} \right)^{-1}}
#' 
#' @useDynLib BFRdist, .registration = TRUE
#'
#' @param x,q vector of quantiles (values must be between 0 and 1).
#' @param p vector of probabilities (values must be between 0 and 1).
#' @param n number of observations. If `length(n) > 1`, the length is taken to be the number required.
#' @param alpha positive shape parameter.
#' @param beta positive shape parameter.
#' 
#' @return 
#' `dbfr` gives the density, `pbfr` gives the distribution function (CDF), 
#' `qbfr` gives the quantile function, and `rbfr` generates random deviates.
#' 
#' The length of the result is determined by `n` for `rbfr`, and is the maximum 
#' of the lengths of the numerical arguments for the other functions.
#' 
#' @import Rcpp
#' @name BFRdist
#' 
#' @examples
#' # Set shape parameters
#' alpha_val <- 1.5
#' beta_val <- 2.0
#' 
#' # Calculate the density at x = 0.5
#' dbfr(0.5, alpha = alpha_val, beta = beta_val)
#' 
#' # Calculate the cumulative probability up to q = 0.5
#' pbfr(0.5, alpha = alpha_val, beta = beta_val)
#' 
#' # Find the 50th percentile (median)
#' qbfr(0.5, alpha = alpha_val, beta = beta_val)
#' 
#' # Generate 10 random observations
#' set.seed(123)
#' r_data <- rbfr(10, alpha = alpha_val, beta = beta_val)
#' print(r_data)
#' 
#' # Plot the density curve across the unit interval
#' x_seq <- seq(0.01, 0.99, length.out = 100)
#' y_dens <- dbfr(x_seq, alpha = alpha_val, beta = beta_val)
#' plot(x_seq, y_dens, type = "l", col = "blue", lwd = 2,
#'      main = "BFR Density (alpha=1.5, beta=2)",
#'      xlab = "x", ylab = "Density")
NULL

#' @rdname BFRdist
#' @export
dbfr <- function(x, alpha, beta) {
  if (any(alpha <= 0) || any(beta <= 0)) stop("alpha and beta must be positive.")
  
  # The PDF is only defined for 0 < x < 1
  ifelse(x > 0 & x < 1,
         (2 * beta / pi) * ((alpha^beta * x^(beta - 1) * (1 - x)^(beta - 1)) / 
                              (alpha^(2 * beta) * x^(2 * beta) + (1 - x)^(2 * beta))),
         0)
}

#' @rdname BFRdist
#' @export
pbfr <- function(q, alpha, beta) {
  if (any(alpha <= 0) || any(beta <= 0)) stop("alpha and beta must be positive.")
  
  # The CDF is 0 for x <= 0 and 1 for x >= 1
  ifelse(q <= 0, 0,
         ifelse(q >= 1, 1,
                (2 / pi) * atan(((alpha * q) / (1 - q))^beta)))
}

#' @rdname BFRdist
#' @export
qbfr <- function(p, alpha, beta) {
  if (any(alpha <= 0) || any(beta <= 0)) stop("alpha and beta must be positive.")
  if (any(p < 0 | p > 1, na.rm = TRUE)) stop("Probabilities p must be between 0 and 1.")
  
  # Quantile function
  (1 + alpha * (tan(pi * p / 2))^(-1 / beta))^(-1)
}

#' @rdname BFRdist
#' @export
rbfr <- function(n, alpha, beta) {
  if (any(alpha <= 0) || any(beta <= 0)) stop("alpha and beta must be positive.")
  
  # Random generation using the inverse transform method
  u <- stats::runif(n)
  qbfr(u, alpha, beta)
}