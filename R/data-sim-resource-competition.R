#' @title Get a data.frame of simulated abundances (resource competition)
#' @aliases simulate_resource_competition_12sp
#' @description Simulate time series from the resource competition model of
#'   Huisman & Weissing (1999) <\url{https://doi.org/10.1038/46540}>. The
#'   initial state, parameters, and perturbations are for the 5-species system
#'   competing on 3 resources as described in Huisman & Weissing (2001)
#'   <\url{https://doi.org/10.1086/319929}>.
#'
#' `simulate_resource_competition_12sp()` is a simulation with the same
#'   structure, but is the 12-species system competing on 5 resources, with
#'   introductions of species 6-12 at later times.
#' @param params model parameters
#' @param initial_state initial conditions
#' @param sample_times the time values at which to make observations
#' @param ... remaining arguments to be passed to \code{\link[deSolve]{ode}}
#'
#' @return a matrix (and `deSolve`) object with the observations (and times)
#'
#' @export
simulate_resource_competition <- function(params = list(
                                            r = rep.int(1, 5),
                                            m = rep.int(0.25, 5),
                                            D = 0.25,
                                            K = matrix(c(
                                              0.20, 0.05, 0.50, 0.05, 0.50,
                                              0.15, 0.06, 0.05, 0.50, 0.30,
                                              0.15, 0.50, 0.30, 0.06, 0.05
                                            ), nrow = 3, byrow = TRUE),
                                            C = matrix(c(
                                              0.20, 0.10, 0.10, 0.10, 0.10,
                                              0.10, 0.20, 0.10, 0.10, 0.20,
                                              0.10, 0.10, 0.20, 0.20, 0.10
                                            ), nrow = 3, byrow = TRUE)
                                          ),
                                          initial_state = c(
                                            R = c(10, 10, 10),
                                            N = c(0.1, 0.1, 0.1, 0.1, 0.1)
                                          ),
                                          sample_times = seq(0, 3000, by = 1),
                                          ...) {
  num_species <- length(grep("N[0-9]+", names(initial_state)))
  num_resources <- length(grep("R[0-9]+", names(initial_state)))
  S <- utils::head(initial_state, num_resources)

  model_equations <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      R <- utils::head(state, num_resources)
      N <- utils::tail(state, num_species)
      mu <- r * apply(R / (K + R), 2, min)
      list(c(
        D * (S - R) - C %*% (mu * N),
        N * (mu - m)
      ))
    })
  }

  deSolve::ode(
    y = initial_state,
    times = sample_times,
    func = model_equations,
    parms = params,
    ...
  )
}

#' @rdname simulate_resource_competition
#'
#' @param events the data.frame describing the introductions of later species
#' @export
sim_resource_competition_12sp <- function(params = list(
                                            r = rep.int(1, 12),
                                            m = rep.int(0.25, 12),
                                            D = 0.25,
                                            K = matrix(c(
                                              0.39, 0.34, 0.30, 0.24, 0.23, 0.41, 0.20, 0.45, 0.14, 0.15, 0.38, 0.28,
                                              0.22, 0.39, 0.34, 0.30, 0.27, 0.16, 0.15, 0.05, 0.38, 0.29, 0.37, 0.31,
                                              0.27, 0.22, 0.39, 0.34, 0.30, 0.07, 0.11, 0.05, 0.38, 0.41, 0.24, 0.25,
                                              0.30, 0.24, 0.22, 0.39, 0.34, 0.28, 0.12, 0.13, 0.27, 0.33, 0.04, 0.41,
                                              0.34, 0.30, 0.22, 0.20, 0.39, 0.40, 0.50, 0.26, 0.12, 0.29, 0.09, 0.16
                                            ), nrow = 5, byrow = TRUE),
                                            C = matrix(c(
                                              0.04, 0.04, 0.07, 0.04, 0.04, 0.22, 0.10, 0.08, 0.02, 0.17, 0.25, 0.03,
                                              0.08, 0.08, 0.08, 0.10, 0.08, 0.14, 0.22, 0.04, 0.18, 0.06, 0.20, 0.04,
                                              0.10, 0.10, 0.10, 0.10, 0.14, 0.22, 0.24, 0.12, 0.03, 0.24, 0.17, 0.01,
                                              0.05, 0.03, 0.03, 0.03, 0.03, 0.09, 0.07, 0.06, 0.03, 0.03, 0.11, 0.05,
                                              0.07, 0.09, 0.07, 0.07, 0.07, 0.05, 0.24, 0.05, 0.08, 0.10, 0.02, 0.04
                                            ), nrow = 5, byrow = TRUE)
                                          ),
                                          initial_state = c(
                                            R = c(6, 10, 14, 4, 9),
                                            N = c(0.11, 0.12, 0.13, 0.14, 0.15, 0, 0, 0, 0, 0, 0, 0)
                                          ),
                                          sample_times = seq(0, 20000, by = 1),
                                          events = list(data = data.frame(
                                            var = c(
                                              "N6", "N7", "N8",
                                              "N9", "N10",
                                              "N11", "N12"
                                            ),
                                            time = c(
                                              1000, 1000, 1000,
                                              3000, 3000,
                                              5000, 5000
                                            ),
                                            value = c(
                                              0.1, 0.1, 0.1,
                                              0.1, 0.1,
                                              0.1, 0.1
                                            ),
                                            method = c(
                                              "rep", "rep", "rep",
                                              "rep", "rep",
                                              "rep", "rep"
                                            )
                                          )),
                                          ...) {
  sim_resource_competition_12sp(
    params = params,
    initial_state = initial_state,
    sample_times = sample_times,
    events = events,
    ...
  )
}
