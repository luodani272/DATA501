#' Fit the BFR distribution
#'
#' Estimates the shape parameters (\code{alpha} and \code{beta}) of the BFR distribution 
#' from empirical data (such as Generalized Radar Vegetation Index datasets) using 
#' Maximum Likelihood Estimation (MLE) or the Quantile Estimator (QE).
#'
#' @details
#' Because the BFR distribution is bounded on the unit interval, all input data 
#' must be strictly between 0 and 1. 
#' 
#' The estimation process utilizes a highly efficient C++ implementation of 
#' the quantile-based estimator (QE) to calculate exact parameter values based on 
#' sample quartiles. If \code{method = "qe"} is selected, these values are returned instantly. 
#' If \code{method = "mle"} is selected, these QE values are passed into \code{stats::optim()} 
#' as highly accurate starting guesses to rapidly refine the estimates via Maximum Likelihood.
#'
#' @param x A numeric vector of observations. Values must be strictly greater than 0 and less than 1. Missing (\code{NA}) and infinite values are not permitted.
#' @param method A character string specifying the estimation method. Options are \code{"mle"} (Maximum Likelihood) or \code{"qe"} (Quantile Estimator).
#'
#' @return An S3 object of class \code{bfr} containing the following components:
#' \itemize{
#'   \item \code{estimate}: A named numeric vector containing the fitted \code{alpha} and \code{beta} parameters.
#'   \item \code{method}: The estimation method used (\code{"mle"} or \code{"qe"}).
#'   \item \code{data}: The original numeric vector of data used for fitting.
#'   \item \code{optim_results}: The complete output list from the internal \code{stats::optim()} call (will be \code{NULL} if \code{method = "qe"}).
#' }
#'
#' @seealso \code{\link{BFRdist}} for the core distribution functions.
#' 
#' @examples
#' # Generate synthetic BFR data 
#' set.seed(123)
#' simulated_data <- rbfr(n = 200, alpha = 2, beta = 3)
#' 
#' # Fit the data using the fast Quantile Estimator (QE)
#' fit_qe <- fit_bfr(simulated_data, method = "qe")
#' print(fit_qe$estimate)
#' 
#' # Fit the data using Maximum Likelihood Estimation (MLE)
#' fit_ml <- fit_bfr(simulated_data, method = "mle")
#' print(fit_ml$estimate)
#' 
#' @references
#' Bhattacharya, A., Frery, A. C., & Rey, A. A. (2026). A Statistical Framework 
#' for Bounded Measurements with Applications to SAR Data. IEEE Transactions 
#' on Geoscience and Remote Sensing.
#' 
#' @export
fit_bfr <- function(x, method = c("mle", "qe")) {
  method <- match.arg(method)
  
  # Data validation checks
  if (any(is.na(x)) || any(is.infinite(x))) stop("NA and Inf values are not permitted.")
  if (any(x <= 0) || any(x >= 1)) stop("All input data must be strictly greater than 0 and less than 1.")
  
  # Use C++ quantile estimator for initial starting values (or final QE values)
  init_est <- estimateBFR(x) 
  init_alpha <- init_est["alpha"]
  init_beta <- init_est["beta"]
  
  # Quantile Estimator (QE) - Return instantly without optimization
  if (method == "qe") {
    out <- list(
      estimate = c(alpha = unname(init_alpha), beta = unname(init_beta)),
      method = method,
      data = x,
      optim_results = NULL
    )
    class(out) <- "bfr"
    return(out)
  }
  
  # Maximum Likelihood Estimation (MLE)
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
    
    # Construct and return the MLE S3 object
    out <- list(
      estimate = c(alpha = unname(opt$par[1]), beta = unname(opt$par[2])),
      method = method,
      data = x,
      optim_results = opt
    )
    class(out) <- "bfr"
    return(out)
  }
}