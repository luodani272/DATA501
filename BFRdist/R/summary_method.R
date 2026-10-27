#' Summary Method for BFR Model Fits
#'
#' Calculates and outputs standard errors and goodness-of-fit metrics 
#' for a fitted BFR model.
#'
#' @details
#' This function is an S3 method for objects of class \code{bfr}. The output 
#' depends on the estimation method used during fitting:
#' 
#' \itemize{
#'   \item \strong{Maximum Likelihood Estimation (MLE):} Provides a complete 
#'     summary including parameter estimates, approximated standard errors 
#'     (calculated by numerically computing and inverting the Hessian matrix), 
#'     and full goodness-of-fit metrics (Log-Likelihood, AIC, and BIC).
#'   \item \strong{Quantile Estimator (QE):} Provides the parameter estimates 
#'     and the Log-Likelihood of those estimates. Because QE is a closed-form 
#'     algebraic estimator rather than a likelihood-based optimization, standard 
#'     errors cannot be derived via the Hessian matrix and are returned as \code{NA}. 
#'     Information criteria (AIC/BIC) are also typically omitted in the print output 
#'     for QE fits.
#' }
#'
#' @param object An object of S3 class \code{bfr}, typically the result of \code{fit_bfr}.
#' @param ... Additional arguments affecting the summary produced. (Currently unused).
#' 
#' @return An S3 object of class \code{summary.bfr} containing:
#' \itemize{
#'   \item \code{estimate}: The fitted shape parameters (\code{alpha} and \code{beta}).
#'   \item \code{se}: The standard errors of the estimates (Numeric for MLE, \code{NA} for QE).
#'   \item \code{method}: The estimation method used (\code{"mle"} or \code{"qe"}).
#'   \item \code{loglik}: The Log-Likelihood of the fitted model.
#'   \item \code{aic}: The Akaike Information Criterion.
#'   \item \code{bic}: The Bayesian Information Criterion.
#'   \item \code{n}: The number of observations in the data.
#' }
#' 
#' @seealso \code{\link{fit_bfr}} for model fitting, and \code{\link{print.summary.bfr}} for the print method.
#' 
#' @examples
#' \donttest{
#' # Generate synthetic data
#' set.seed(123)
#' raw_data <- rbfr(n = 100, alpha = 1.5, beta = 2.5)
#' 
#' # Fit using the primary Quantile Estimator (QE) for speed
#' # Note: Standard errors will be NA.
#' bfr_qe <- fit_bfr(raw_data, method = "qe")
#' qe_summary <- summary(bfr_qe)
#' qe_summary
#' 
#' # Fit using MLE when standard errors and complete goodness-of-fit metrics are required
#' bfr_mle <- fit_bfr(raw_data, method = "mle")
#' mle_summary <- summary(bfr_mle)
#' mle_summary
#' }
#' 
#' @method summary bfr
#' @export
#' 
summary.bfr <- function(object, ...) {
  x <- object$data
  n <- length(x)
  a <- unname(object$estimate["alpha"])
  b <- unname(object$estimate["beta"])
  
  # Calculate Log-Likelihood
  pdf_vals <- dbfr(x, alpha = a, beta = b)
  # Filter out zeros to prevent log(0) issues
  ll <- sum(log(pdf_vals[pdf_vals > 0])) 
  
  # Goodness-of-Fit Metrics
  k <- 2 # number of parameters
  aic <- 2 * k - 2 * ll
  bic <- k * log(n) - 2 * ll
  
  # Approximate Standard Errors using the Hessian matrix (for MLE only)
  se <- c(alpha = NA_real_, beta = NA_real_)
  if (object$method == "mle") {
    nll <- function(par) {
      -sum(log(dbfr(x, par[1], par[2])))
    }
    # Calculate Hessian numerically
    hess <- stats::optimHess(c(a, b), nll)
    
    # Invert the Hessian to get the covariance matrix
    cov_mat <- tryCatch(solve(hess), error = function(e) matrix(NA, 2, 2))
    se <- sqrt(diag(cov_mat))
    names(se) <- c("alpha", "beta")
  }
  
  # Bundle into a summary object
  res <- list(
    estimate = object$estimate,
    se = se,
    method = object$method,
    loglik = ll,
    aic = aic,
    bic = bic,
    n = n
  )
  class(res) <- "summary.bfr"
  return(res)
}

