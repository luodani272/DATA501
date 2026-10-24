#' Fit the BFR distribution
#'
#' Estimates the shape parameters (\code{alpha} and \code{beta}) of the BFR distribution 
#' from empirical data (such as Generalized Radar Vegetation Index datasets) using 
#' Maximum Likelihood Estimation (MLE) or the Method of Moments (MME).
#'
#' @details
#' Because the BFR distribution is bounded on the unit interval, all input data 
#' must be strictly between 0 and 1. 
#' 
#' The estimation process utilizes a two-step approach:
#' 1. Initial parameter estimates are computed quickly using a C++ implementation of 
#'    the quantile-based estimator.
#' 2. These initial values are passed into \code{stats::optim()} using the "L-BFGS-B" 
#'    algorithm to refine the estimates via either Maximum Likelihood (MLE) or 
#'    the Method of Moments (MME).
#'
#' @param x A numeric vector of observations. Values must be strictly greater than 0 and less than 1. Missing (\code{NA}) and infinite values are not permitted.
#' @param method A character string specifying the estimation method. Options are \code{"mle"} (Maximum Likelihood) or \code{"mme"} (Method of Moments).
#'
#' @return An S3 object of class \code{bfr} containing the following components:
#' \itemize{
#'   \item \code{estimate}: A named numeric vector containing the fitted \code{alpha} and \code{beta} parameters.
#'   \item \code{method}: The estimation method used (\code{"mle"} or \code{"mme"}).
#'   \item \code{data}: The original numeric vector of data used for fitting.
#'   \item \code{optim_results}: The complete output list from the internal \code{stats::optim()} call.
#' }
#'
#' @seealso \code{\link{BFRdist}} for the core distribution functions.
#' 
#' @examples
#' # 1. Generate synthetic BFR data (e.g., simulating GRVI observations)
#' set.seed(42)
#' simulated_data <- rbfr(n = 200, alpha = 2, beta = 3)
#' 
#' # 2. Fit the data using Maximum Likelihood Estimation (MLE)
#' fit_ml <- fit_bfr(simulated_data, method = "mle")
#' print(fit_ml$estimate)
#' 
#' # 3. Fit the data using the Method of Moments (MME)
#' fit_moments <- fit_bfr(simulated_data, method = "mme")
#' print(fit_moments$estimate)
#' 
#' @export
fit_bfr <- function(x, method = c("mle", "mme")) {
  method <- match.arg(method)
  
  # Data validation checks
  if (any(is.na(x)) || any(is.infinite(x))) stop("NA and Inf values are not permitted.")
  if (any(x <= 0) || any(x >= 1)) stop("All input data must be strictly greater than 0 and less than 1.")
  
  # Use C++ quantile estimator for initial starting values
  init_est <- estimateBFR(x) 
  init_alpha <- init_est["alpha"]
  init_beta <- init_est["beta"]
  
  # Maximum Likelihood Estimation
  if (method == "mle") {
    nll <- function(par) {
      a <- par[1]; b <- par[2]
      if (a <= 0 || b <= 0) return(Inf)
      pdf_vals <- dbfr(x, alpha = a, beta = b)
      if (any(pdf_vals <= 0, na.rm = TRUE)) return(Inf)
      -sum(log(pdf_vals))
    }
    opt <- stats::optim(par = c(init_alpha, init_beta), fn = nll, 
                 method = "L-BFGS-B", lower = c(1e-5, 1e-5))
    
    # Method of Moments Estimation
  } else if (method == "mme") {
    m1 <- mean(x)
    m2 <- mean(x^2)
    
    moment_obj <- function(par) {
      a <- par[1]; b <- par[2]
      if (a <= 0 || b <= 0) return(Inf)
      
      e_x <- tryCatch(stats::integrate(function(u) u * dbfr(u, a, b), 0, 1)$value, error = function(e) NA)
      e_x2 <- tryCatch(stats::integrate(function(u) u^2 * dbfr(u, a, b), 0, 1)$value, error = function(e) NA)
      
      if (is.na(e_x) || is.na(e_x2)) return(Inf)
      (e_x - m1)^2 + (e_x2 - m2)^2
    }
    opt <- stats::optim(par = c(init_alpha, init_beta), fn = moment_obj, 
                 method = "L-BFGS-B", lower = c(1e-5, 1e-5))
  }
  
  # Construct and return the S3 object
  out <- list(
    estimate = c(alpha = unname(opt$par[1]), beta = unname(opt$par[2])),
    method = method,
    data = x,
    optim_results = opt
  )
  class(out) <- "bfr"
  return(out)
}