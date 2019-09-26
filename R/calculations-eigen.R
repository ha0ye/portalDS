#' @title Compute the eigen-decomposition of the smap matrices

#' @param smap_matrices A list with the matrix of smap-coefficients at each
#'   time point \code{\link{compute_smap_matrices}}
#' @return A list with two elements:
#'   \describe{
#'     \item{`values`}{a list of the eigenvalues (a vector) for each time point}
#'     \item{`vectors`}{a list of the eigenvectors (a matrix, each column is an
#'       eigenvector) for each time point}
#'   }
#'
#' @export
compute_eigen_decomp <- function(smap_matrices) {
  eigen_decomp <- purrr::map(smap_matrices, function(J) {
    if (any(is.na(J))) {
      return(c("values" = NA, "vectors" = NA))
    }
    out <- eigen(J)
    rownames(out$vectors) <- rownames(J)
    return(out)
  })
  return(purrr::transpose(eigen_decomp, .names = c("values", "vectors")))
}
