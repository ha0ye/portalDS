#' @title Get a data.frame of simulated abundances (3-species food chain)
#' @description Simulate time series from the 3-species food chain model of 
#'   Hastings & Powell (1991) <\url{https://doi.org/10.2307/1940591}>. The 
#'   chosen initial state and parameters should give chaotic dynamics that 
#'   produce a "teacup" attractor. 
#' @param params model parameters 
#' @param initial_state initial conditions
#' @param sample_times the time values at which to make observations
#' @param ... remaining arguments to be passed to \code{\link[deSolve]{ode}}
#'
#' @return a matrix (and `deSolve`) object with the observations (and times)
#' 
#' @export
simulate_3sp_food_chain <- function(params = c(a_1 = 2.5, b_1 = 3.2, 
                                               a_2 = 0.1, b_2 = 2.0, 
                                               d_1 = 0.2, d_2 = 0.015), 
                                    initial_state = c(x = 0.8, y = 0.2, z = 8), 
                                    sample_times = seq(0, 5000, by = 5), 
                                    ...)
{
    
    model_equations <- function(t, state, parameters)
    {
        with(as.list(c(state, parameters)), {
            f1x_y <- (a_1*x) / (1 + b_1*x) * y
            f2y_z <- (a_2*y) / (1 + b_2*y) * z
            list(c(x*(1-x) - f1x_y, 
                   f1x_y - f2y_z - d_1*y, 
                   f2y_z - d_2*z))
        })
    }
    
    deSolve::ode(y = initial_state, 
                 times = sample_times, 
                 func = model_equations, 
                 parms = params, 
                 ...)
}
