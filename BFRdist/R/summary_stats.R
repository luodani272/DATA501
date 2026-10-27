#' Galton Skewness
#'
#' Calculates Galton's measure of skewness for a numeric vector based on sample quartiles.
#'
#' @details
#' Galton's skewness is an order statistic-based measure of asymmetry. 
#' Unlike traditional skewness metrics that rely on the second and third central moments, 
#' Galton's skewness depends solely on the first, second (median), and third quartiles. 
#' This makes it significantly more robust to outliers and applicable to distributions 
#' where theoretical moments may not exist. 
#' 
#' By default, this function computes sample quantiles using R's default Type 7 algorithm.
#'
#' @param x A numeric vector of observations.
#' @param na.rm Logical. Should missing values be removed before computation? Default is \code{FALSE}.
#'
#' @return A numeric value representing the Galton skewness. Returns \code{NA_real_} 
#'   with a warning if the interquartile range (the denominator) is zero.
#'   
#' @seealso \code{\link{moors_kurtosis}} for the complementary kurtosis measure, 
#'   and \code{\link{plot.bfr}} for their application in distributional diagrams.
#' 
#' @examples
#' set.seed(123)
#' data <- rnorm(100, mean = 0, sd = 1)
#' galton_skewness(data)
#' 
#' @export
galton_skewness <- function(x, na.rm = FALSE) {
  if (na.rm) {
    x <- stats::na.omit(x)
  }
  
  # Calculate required quantiles (Type 7 is R's default)
  q <- stats::quantile(x, probs = c(0.25, 0.5, 0.75), type = 7)
  
  denominator <- q[[3]] - q[[1]]
  
  # Avoid division by zero if Q3 == Q1
  if (denominator == 0) {
    warning("Interquartile range is zero. Returning NA.")
    return(NA_real_)
  }
  
  (q[[3]] + q[[1]] - 2 * q[[2]]) / denominator
}

#' Moors Kurtosis
#'
#' Calculates Moors' measure of kurtosis for a numeric vector based on sample octiles.
#'
#' @details
#' Moors' kurtosis is an order statistic-based measure of peakedness and tail heaviness. 
#' Unlike traditional kurtosis metrics that rely on the second and fourth central moments, 
#' Moors' kurtosis depends exclusively on octiles (the 12.5th, 25th, 37.5th, 62.5th, 75th, 
#' and 87.5th percentiles). This provides a robust alternative to moment-based 
#' kurtosis that is less sensitive to extreme outliers.
#' 
#' By default, this function computes sample quantiles using R's default Type 7 algorithm.
#'
#' @param x A numeric vector of observations.
#' @param na.rm Logical. Should missing values be removed before computation? Default is \code{FALSE}.
#'
#' @return A numeric value representing the Moors kurtosis. Returns \code{NA_real_} 
#'   with a warning if the interquartile range (the denominator) is zero.
#'   
#' @seealso \code{\link{galton_skewness}} for the complementary skewness measure, 
#'   and \code{\link{plot.bfr}} for their application in distributional diagrams.
#' 
#' @examples
#' set.seed(123)
#' data <- rnorm(100, mean = 0, sd = 1)
#' moors_kurtosis(data)
#' 
#' @export
moors_kurtosis <- function(x, na.rm = FALSE) {
  if (na.rm) {
    x <- stats::na.omit(x)
  }
  
  # Calculate required quantiles (Type 7 is R's default)
  q <- stats::quantile(x, probs = c(0.125, 0.25, 0.375, 0.625, 0.75, 0.875), type = 7)
  
  denominator <- q[[5]] - q[[2]]
  
  # Avoid division by zero if Q3 == Q1
  if (denominator == 0) {
    warning("Interquartile range is zero. Returning NA.")
    return(NA_real_)
  }
  
  (q[[6]] - q[[4]] + q[[3]] - q[[1]]) / denominator
}