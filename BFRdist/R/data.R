#' Multi-Temporal GRVI and DpRVI Measurements for Agricultural Crops
#'
#' A real-world dataset containing Generalized Radar Vegetation Index (GRVI) 
#' and Dual-polarization Radar Vegetation Index (DpRVI) measurements across 
#' multiple dates. The dataset includes 2013 total observations distributed 
#' across three crop types: Spring Wheat (750), Soybeans (637), and Canola (626).
#'
#' @format A data frame with 2013 rows and 12 variables:
#' \describe{
#'   \item{S.NO}{Sample identification number}
#'   \item{GRVI_06Jun}{GRVI measurement on June 6}
#'   \item{GRVI_23Jun}{GRVI measurement on June 23}
#'   \item{GRVI_24_july}{GRVI measurement on July 24}
#'   \item{GRVI_17Aug}{GRVI measurement on August 17}
#'   \item{GRVI_27Aug}{GRVI measurement on August 27}
#'   \item{Dprvi_06_June}{DpRVI measurement on June 6}
#'   \item{Dprvi_23_June}{DpRVI measurement on June 23}
#'   \item{Dprvi_24_July}{DpRVI measurement on July 24}
#'   \item{DprVI_17_Aug}{DpRVI measurement on August 17}
#'   \item{Dprvi_27_Aug}{DpRVI measurement on August 27}
#'   \item{Crop}{Type of crop (Canola, Soybeans, Spring Wheat)}
#' }
#' 
#' @details 
#' GRVI and DpRVI are bounded remote sensing metrics used to monitor crop canopy 
#' growth and structure using synthetic aperture radar (SAR). 
#' 
#' **Important:** Some late-season measurements (e.g., Spring Wheat on August 27) 
#' contain values of exactly 0.0. Because the BFR distribution is strictly defined 
#' on the open interval (0, 1), users must filter or transform these boundary 
#' values before applying `fit_bfr()`.
#'
#' @source
#' The original multi-temporal dataset was provided via personal communication.
#' 
#' Bhattacharya, A., Frery, A. C., & Rey, A. A. (2026). GRVI data [Data set]. Personal communication.
"grvi_data"