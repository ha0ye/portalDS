#' @import ggraph
#' @import ggplot2
#' @importFrom stats predict
#' @importFrom stats quantile
#' @importFrom stats sd
#' @importFrom stats smooth.spline
#' @importFrom utils data
NULL

if (getRversion() >= "2.15.1") utils::globalVariables(
    c(".", "1", "2", "abundance", "BA", "best_E", "censusdate", "coeff_name", 
      "const", "data_type", "delta_rho", "E", "lambda", 
      "lib_column", "lib_size", "mae", "name", "newmoonnumber", "nplots", 
      "ntraps", "num_pred", "period", "predictor", "quantile_q", "rho", 
      "rho_minus_upper_q_null", "rmse", "simplex_out", "SO", "species", 
      "surr_CI", "target", "target_column", "upper_q", "Var1", "Var2", 
      "variable", "xmax", "xmin", "y", "ymax", "ymin"))