#' @title Get a data.frame of simulated abundances (3-stage flour beetles)
#' @description Simulate time series from the LPA (larvae-pupae-adults) flour
#'   beetle model of Dennis et al. (2001) <\url{https://doi.org/10.1890/0012-9615(2001)071[0277:ECACDI]2.0.CO;2}>
#'   . The initial conditions follow the experimental setup of Dennis et al.
#'   (2001) and the parameters are the maximum-likelihood estimates from Table 1
#'   of the same, with u_a and c_pa adjusted according to the description so as
#'   to give chaotic behavior.
#' @param params model parameters
#' @param initial_state initial conditions
#' @param sample_times the time values at which to make observations
#' @param ... remaining arguments to be passed to \code{\link[deSolve]{ode}}
#'
#' @return a matrix (and `deSolve`) object with the observations (and times)
#'
#' @export
simulate_LPA_flour_beetles <- function(params = c(
                                         b = 10.67,
                                         u_l = 0.1955,
                                         u_a = 0.96,
                                         c_el = 0.01647,
                                         c_ea = 0.01313,
                                         c_pa = 0.35
                                       ),
                                       initial_state = c(L = 250, P = 5, A = 100),
                                       sample_times = seq(0, 1000),
                                       ...) {
  model_equations <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      list(c(
        b * A * exp(-c_el * L - c_ea * A),
        L * (1 - u_l),
        P * exp(-c_pa * A) + A * (1 - u_a)
      ))
    })
  }

  deSolve::ode(
    y = initial_state,
    times = sample_times,
    func = model_equations,
    parms = params,
    method = "iteration",
    ...
  )
}
