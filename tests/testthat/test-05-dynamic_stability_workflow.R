context("Check full dynamic stability workflow")

test_that("Compute dynamic stability using block_3sp example data", {
    set.seed(616382247)

    # setup data
    data("block_3sp", package = "rEDM")
    var_names <- c("x", "y", "z")
    block <- setNames(block_3sp[, c("time", "x_t", "y_t", "z_t")], 
                      c("censusdate", var_names))
    
    # simplex results
    expect_error(simplex_results <- compute_simplex(block, 
                                                    E_list = 3:5, 
                                                    num_surr = 4, 
                                                    surrogate_method = "random_shuffle"), NA)
    # round simplex outputs for hash
    simplex_results$simplex_out <- lapply(simplex_results$simplex_out, round, 4)
    expect_known_hash(simplex_results, "ae7a4100e9")

    # ccm results
    expect_error(ccm_results <- compute_ccm(simplex_results, 
                                            num_samples = 10), NA)
    expect_equal(dim(ccm_results), c(450, 9))
    expect_setequal(as.character(ccm_results$lib_column), var_names)
    expect_setequal(as.character(ccm_results$target_column), var_names)
    expect_setequal(as.character(ccm_results$data_type), c("actual", "surrogate"))
    expect_setequal(as.character(ccm_results$data_type), c("actual", "surrogate"))
    expect_setequal(ccm_results$lib_size, seq(10, 100, 10))
    expect_known_hash(dplyr::mutate_at(ccm_results, c("rho", "mae", "rmse"), round, 4), 
                      "40a8eb1668")
    
    # ccm links
    expect_error(ccm_links <- compute_ccm_links(ccm_results), NA)
    expect_known_hash(ccm_links, "4351514392")

    # s-map coefficients
    expect_error(smap_coeffs <- compute_smap_coeffs(block, ccm_links), NA)
    expect_known_hash(lapply(smap_coeffs, round, 4), "3f88e8ddf9")
    
    expect_error(smap_matrices <- compute_smap_matrices(smap_coeffs, ccm_links), NA)
    expect_known_hash(lapply(smap_matrices, round, 4), "5393993b55")
    
    expect_error(eigen_decomp <- compute_eigen_decomp(smap_matrices), NA)
    expect_error(eigenvalues <- eigen_decomp$values, NA)
    expect_known_hash(lapply(eigenvalues, round, 4), "c5a2f2aaca")
    
    expect_error(eigenvectors <- eigen_decomp$vectors, NA)
    expect_known_hash(which(is.na(eigenvectors)), "ac6cf5b073")
    expect_error(ev <- eigenvectors[!is.na(eigenvectors)], NA)
    expect_known_hash(lapply(ev, function(x) {
        s <- ifelse(Re(x[1, ]) == 0, 1, sign(Re(x[1, ])))
        s <- matrix(s, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
        round(s * x, 4)}), 
        "68244b7ef9")
})

test_that("Verify compute_dynamic_stability function", {
    set.seed(616382247)
    
    # setup data
    data("block_3sp", package = "rEDM")
    var_names <- c("x", "y", "z")
    block <- setNames(block_3sp[, c("time", "x_t", "y_t", "z_t")], 
                      c("censusdate", var_names))
    
    # simplex results
    expect_error(results <- compute_dynamic_stability(block, 
                                                      E_list = 3:5, 
                                                      surrogate_method = "random_shuffle", 
                                                      num_surr = 4, 
                                                      num_samples = 10), NA)
    expect_is(results, "list")
    expect_true(all(c("block", "simplex_results", "ccm_results", 
                      "ccm_links", "smap_matrices", "eigenvalues", 
                      "eigenvectors") %in% names(results)))
    
    simplex_results <- results$simplex_results
    ccm_results <- results$ccm_results
    ccm_links <- results$ccm_links
    smap_coeffs <- results$smap_coeffs
    smap_matrices <- results$smap_matrices
    eigenvalues <- results$eigenvalues
    eigenvectors <- results$eigenvectors
    
    simplex_results$simplex_out <- lapply(simplex_results$simplex_out, round, 4)
    expect_known_hash(simplex_results, "ae7a4100e9")
    expect_known_hash(dplyr::mutate_at(ccm_results, c("rho", "mae", "rmse"), round, 4), 
                      "40a8eb1668")
    expect_known_hash(ccm_links, "4351514392")
    expect_known_hash(lapply(smap_coeffs, round, 4), "3f88e8ddf9")
    expect_known_hash(lapply(smap_matrices, round, 4), "5393993b55")
    expect_known_hash(lapply(eigenvalues, round, 4), "c5a2f2aaca")
    expect_known_hash(which(is.na(eigenvectors)), "ac6cf5b073")
    expect_error(ev <- eigenvectors[!is.na(eigenvectors)], NA)
    expect_known_hash(lapply(ev, function(x) {
        s <- ifelse(Re(x[1, ]) == 0, 1, sign(Re(x[1, ])))
        s <- matrix(s, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
        round(s * x, 4)}), 
        "68244b7ef9")
})
