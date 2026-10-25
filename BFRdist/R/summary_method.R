#' Summary Method for BFR Model Fits
#'
#' Calculates and outputs standard errors and goodness-of-fit metrics 
#' for a fitted BFR model.
#'
#' @details
#' This function is an S3 method for objects of class \code{bfr}. It calculates 
#' goodness-of-fit metrics including the Log-Likelihood, Akaike Information Criterion (AIC), 
#' and Bayesian Information Criterion (BIC). For models fitted using 
#' Maximum Likelihood Estimation (MLE), it also approximates the standard errors 
#' of the parameter estimates by numerically computing and inverting the Hessian matrix.
#'
#' @param object An object of S3 class \code{bfr}, typically the result of \code{fit_bfr}.
#' @param ... Additional arguments affecting the summary produced. (Currently unused).
#' 
#' @return An S3 object of class \code{summary.bfr} containing:
#' \itemize{
#'   \item \code{estimate}: The fitted shape parameters (\code{alpha} and \code{beta}).
#'   \item \code{se}: The standard errors of the estimates (for MLE fits).
#'   \item \code{method}: The estimation method used.
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
#' # 1. Generate synthetic data
#' set.seed(123)
#' raw_data <- rbfr(n = 100, alpha = 1.5, beta = 2.5)
#' 
#' # 2. Fit the BFR model to the data using MLE
#' bfr_model <- fit_bfr(raw_data, method = "mle")
#' 
#' # 3. Generate the summary statistics
#' model_summary <- summary(bfr_model)
#' 
#' # The print method is called automatically when typing the object name:
#' model_summary
#' }
#' 
#' @method summary bfr
#' @export
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

