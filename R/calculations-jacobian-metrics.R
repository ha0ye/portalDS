#' @title Compute eigen-decompositions of the Jacobian matrices
#' 
#' @param smap_matrices A list with the Jacobian matrix (of smap-coefficients) 
#'   at each time point, resulting from [compute_smap_matrices()]
#' @return A list with two elements:
#'   \describe{
#'     \item{`values`}{a list of the eigenvalues (a vector) for each time point}
#'     \item{`vectors`}{a list of the eigenvectors (a matrix, each column is an
#'       eigenvector) for each time point}
#'   }
#'
#' @export
compute_eigen_decomp <- function(smap_matrices)
{
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

#' @title Compute singular value decompositions of the Jacobian matrics
#' 
#' @details The full Jacobian resulting from [compute_smap_matrices()] is of 
#'   the form, J = 
#'   \tabular{rrrrr}{
#'     C^0 \tab C^1 \tab ... \tab C^(d-1) \tab C^d\cr
#'       I \tab   0 \tab ... \tab       0 \tab   0\cr
#'       0 \tab   I \tab ... \tab       0 \tab   0\cr
#'     ... \tab ... \tab ... \tab     ... \tab ...\cr
#'       0 \tab   0 \tab ... \tab       I \tab   0
#'    }
#'   Note that this maps from the column vector \[N(t) N(t-1) ... N(t-d)\]^T to 
#'   the column vector \[N(t+1) N(t) ... N(t-(d-1))\]^T. However, the only 
#'   relevant components for our purposes are the rows which map the column 
#'   vector \[N(t) N(t-1) ... N(t-d)\]^T to \[N(t+1)\]^T, let this be J_s.
#'   
#'   Thus, we extract this portion of the Jacobian for applying SVD.
#' 
#' @inheritParams compute_eigen_decomp
#' @param s the number of species in the system (optional parameter to restrict 
#'   the analysis just to the portions of the Jacobian that are relevant for the 
#'   forecasts)
#' 
#' @return A list with three elements:
#'   \describe{
#'     \item{`d`}{a list of the singular values (a vector) for each time point}
#'     \item{`u`}{a list of the left singular vectors (a matrix, each column is 
#'       an axis in the output space) for each time point}
#'     \item{`v`}{a list of the right singular vectors (a matrix, each column 
#'       is an axis in the input space) for each time point}
#'   }
#'   
#' @export
compute_svd_decomp <- function(smap_matrices, s = NULL)
{
  if (is.null(s))
  {
    # find first full matrix
    matrix_vars <- lapply(smap_matrices, rownames)
    idx <- min(which(!vapply(matrix_vars, is.null, FALSE)))
    
    # compute s from rownames
    regex_splits <- stringr::str_match(rownames(smap_matrices[[idx]]), "(.+)_(\\d+)")
    s <- max(which(regex_splits[, 3] == "0"))
    var_names <- regex_splits[seq(s), 2]
  }
  
  svd_out <- purrr::map(smap_matrices, function(J) {
    if (any(is.na(J))) {
      return(c("d" = NA, "u" = NA, "v" = NA))
    }
    
    out <- svd(J[seq_len(s), ])
    rownames(out$u) <- var_names
    return(out)
  })
  return(purrr::transpose(svd_out, .names = c("d", "u", "v")))
}

#' @title Compute volume contraction of the Jacobian matrics
#' 
#' @description We define the volume contraction as the rate at which an 
#'   infinitesimal volume scales through the Jacobian. This can be calculated 
#'   as the determinant of the Jacobian (i.e. the product of the eigenvalues or 
#'   singular values, which represent the scaling along various directions in 
#'   the state space).
#' 
#' @details See [compute_svd_decomp()] for details on extracting the portion of 
#'   the Jacobian used for calculating the determinant.
#'   
#'   We do this in order to account for the low-rank of the full Jacobian, which 
#'   otherwise results in a determinant of 0.
#'   
#'   We then compute the pseudo-determinant as |det(J_s %*% t(J_s))|.
#' @inheritParams compute_svd_decomp
#' 
#' @return a numeric vector of the volume contraction values
#' 
#' @export
compute_volume_contraction <- function(smap_matrices, s = NULL)
{
  if (is.null(s))
  {
    # find first full matrix
    matrix_vars <- lapply(smap_matrices, rownames)
    idx <- min(which(!vapply(matrix_vars, is.null, FALSE)))
    
    # compute s from rownames
    regex_splits <- stringr::str_match(rownames(smap_matrices[[idx]]), "(.+)_(\\d+)")
    s <- max(which(regex_splits[, 3] == "0"))
  }
  
  vc_out <- purrr::map_dbl(smap_matrices, function(J) {
    if (any(is.na(J))) {
      return(NA_real_)
    }
    
    J_lags_trimmed <- J[seq_len(s), ]
    out <- sqrt(det(J_lags_trimmed %*% t(J_lags_trimmed)))
    return(out)
  })
  return(vc_out)
}

#' @title Compute total variance of the Jacobian matrics
#' 
#' @description We define the total variance as the variance that results when 
#'   applying the Jacobian to a point with zero mean and unit variance. This can 
#'   be calculated as the trace of the Jacobian multiplied by its transpose. The 
#'   diagonal elements then represent the variance of each dimension, with off-
#'   diagonal elements being the covariance.
#' 
#' @details See [compute_svd_decomp()] for details on extracting the portion of 
#'   the Jacobian used for calculating the determinant.
#'   
#'   We do this in order to account for the low-rank of the full Jacobian, which 
#'   otherwise results in a determinant of 0.
#'   
#'   We then compute the pseudo-determinant as Tr(J_s %*% t(J_s)).
#' @inheritParams compute_svd_decomp
#' 
#' @return a numeric vector of the total variance values
#' 
#' @export
compute_total_variance <- function(smap_matrices, s = NULL)
{
  if (is.null(s))
  {
    # find first full matrix
    matrix_vars <- lapply(smap_matrices, rownames)
    idx <- min(which(!vapply(matrix_vars, is.null, FALSE)))
    
    # compute s from rownames
    regex_splits <- stringr::str_match(rownames(smap_matrices[[idx]]), "(.+)_(\\d+)")
    s <- max(which(regex_splits[, 3] == "0"))
  }
  
  tv_out <- purrr::map_dbl(smap_matrices, function(J) {
    if (any(is.na(J))) {
      return(NA_real_)
    }
    
    J_lags_trimmed <- J[seq_len(s), ]
    out <- sum(diag(J_lags_trimmed %*% t(J_lags_trimmed)))
    return(out)
  })
  return(tv_out)
}