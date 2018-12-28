#' @import ggraph
#' @import ggplot2
#' @importFrom stats predict
#' @importFrom stats quantile
#' @importFrom stats sd
#' @importFrom stats smooth.spline
#' @importFrom utils data
NULL

if (getRversion() >= "2.15.1") utils::globalVariables(
    c(".", "abundance", "BA", "best_E", "censusdate", "coeff_name", "const"))