#' @title linear rescaling function
#' @description rescale the input vector to mean = 0, variance = 1
#'
#' @param x A numeric vector.
#' @param na.rm A logical (default = TRUE) that affects how NAs are dealt with
#'
#' @return A rescaled vector the same length as `x`
#'
#' @export
norm_rescale <- function(x, na.rm = TRUE) {(x - mean(x, na.rm = na.rm)) / sd(x, na.rm = na.rm)}
