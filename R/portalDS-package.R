#' @import ggraph
#' @import ggplot2
#' @importFrom rlang .data
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    # c(".", "1", "2", "abundance", "BA", "best_E", "censusdate", "coeff_name",
    #   "const", "data_type", "delta_rho", "E", "lambda",
    #   "lib_column", "lib_size", "mae", "newmoonnumber", "nplots",
    #   "ntraps", "num_pred", "period", "predictor", "quantile_q", "rho",
    #   "rho_minus_upper_q_null", "rmse", "simplex_out", "SO", "species",
    #   "surr_CI", "target", "target_column", "upper_q", "Var1", "Var2",
    #   "block", "ccm_func", "ccm_links", "ccm_params", "ccm_results",
    #   "eigen_decomp", "from_idx", "get_smap_coefficients", "simplex_results",
    #   "smap_coeffs", "smap_matrices", "to_idx")
  )
}

#' @title Time series for the Maizuru Bay fish community
#' @author Reiji Masuda
#' @description Time series of twice-monthly visual census data of the Maizuru
#'   fish community (only the 15 species that had total observation count >
#'   1000). As used in "Fluctuating interaction network and time-varying
#'   stability of a natural fish community" (Ushio et al. 2018)
#' Data are available from https://zenodo.org/record/1181937
#' @format
#' \describe{
#'   \item{\code{censusdate}}{date of sampling}
#'   \item{\code{Aurelia.sp}}{survey data on Jellyfish}
#'   \item{\code{Engraulis.japonicus}}{survey data on \emph{Engraulis japonicus}}
#'   \item{\code{Plotosus.lineatus}}{survey data on \emph{Plotosus lineatus}}
#'   \item{\code{Sebastes.inermis}}{survey data on \emph{Sebastes inermis}}
#'   \item{\code{Trachurus.japonicus}}{survey data on \emph{Trachurus japonicus}}
#'   \item{\code{Girella.punctata}}{survey data on \emph{Girella punctata}}
#'   \item{\code{Pseudolabrus.sieboldi}}{survey data on \emph{Pseudolabrus sieboldi}}
#'   \item{\code{Halichoeres.poecilopterus}}{survey data on \emph{Halichoeres poecilopterus}}
#'   \item{\code{Halichoeres.tenuispinnis}}{survey data on \emph{Halichoeres tenuispinnis}}
#'   \item{\code{Chaenogobius.gulosus}}{survey data on \emph{Chaenogobius gulosus}}
#'   \item{\code{Pterogobius.zonoleucus}}{survey data on \emph{Pterogobius zonoleucus}}
#'   \item{\code{Tridentiger.trigonocephalus}}{survey data on \emph{Tridentiger trigonocephalus}}
#'   \item{\code{Siganus.fuscescens}}{survey data on \emph{Siganus fuscescens}}
#'   \item{\code{Sphyraena.pinguis}}{survey data on \emph{Sphyraena pinguis}}
#'   \item{\code{Rudarius.ercodes}}{survey data on \emph{Rudarius ercodes}}
#' }
"maizuru_block"
