#' Galton Skewness
#'
#' Calculates the Galton skewness (Bowley's skewness) for a numeric vector.
#'
#' @param x A numeric vector of observations.
#' @param na.rm Logical. Should missing values be removed? Default is FALSE.
#'
#' @return The Galton skewness value as a numeric. Returns NA if the 
#'   interquartile range is zero.
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
#' Calculates the Moors kurtosis for a numeric vector.
#'
#' @param x A numeric vector of observations.
#' @param na.rm Logical. Should missing values be removed? Default is FALSE.
#'
#' @return The Moors kurtosis value as a numeric. Returns NA if the 
#'   interquartile range is zero.
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