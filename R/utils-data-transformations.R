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

#' @title make_surrogate_annual_spline
#' @description Generate surrogate time series by computing a mean seasonal
#'   trend for each year, and shuffling the residuals. This differs from
#'   \code{\link[rEDM]{make_surrogate_seasonal}} in that the data are not
#'   sampled uniformly at the same time each year. Thus, we also ask for the
#'   `day_of_year` as input, in order to compute the mean seasonal trend.
#' @param day_of_year a vector of the numerical day of year, e.g. January 1 = 1,
#'   January 2 = 2, December 31 = 365 (leap years are a bit funny, but I don't
#'   think it should have a large effect)
#' @inheritParams rEDM::make_surrogate_data
#'
#' @return A matrix where each column is a separate surrogate with the same
#'   length as `ts`.
#'
#' @importFrom stats smooth.spline predict
#' @export
make_surrogate_annual_spline <- function(ts, num_surr = 100, day_of_year = NULL)
{
    if (is.data.frame(ts))
    {
        ts <- ts[[1]]
    }
    if (is.null(day_of_year))
    {
        day_of_year <- seq_along(ts)
    }
    
    # filter out NA first
    n <- length(day_of_year)
    idx <- which(is.finite(day_of_year) & is.finite(ts))
    day_of_year <- day_of_year[idx]
    ts <- ts[idx]
    
    xx <- c(day_of_year - 365, day_of_year, day_of_year + 365)
    yy <- c(ts, ts, ts)
    seasonal_F <- smooth.spline(xx, yy)
    seasonal_cyc <- predict(seasonal_F, day_of_year)$y
    seasonal_resid <- ts - seasonal_cyc
    
    vals <- matrix(unlist(
        lapply(seq(num_surr), function(i) {
            seasonal_cyc + sample(seasonal_resid, n)
        })
    ), ncol = num_surr)
    
    out <- matrix(NA, nrow = n, ncol = num_surr)
    out[idx, ] <- vals
}
