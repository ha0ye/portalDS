context("Check functions to build plans for dynamic stability workflow")

test_that("build_ccm_plan produces a correct plan", {
    skip_on_covr()
    expect_error(ccm_plan <- build_ccm_plan(), NA)
    expect_is(ccm_plan, c("drake_plan", "tbl_df"))
    
    targets <- c("ccm_func", "ccm_params", "ccm_results")
    expect_equal(vapply(targets, function(patt) {grep(patt, ccm_plan$target)}, 0), 
                 seq_along(targets), check.names = FALSE)
    expect_known_hash(as.character(ccm_plan$command), "731a7491d4")
})

test_that("check full dynamic stability plan using block_3sp example data", {
    skip_on_covr()
    # setup
    data("block_3sp", package = "rEDM")
    block <- setNames(block_3sp[, c("time", "x_t", "y_t", "z_t")], 
                      c("censusdate", "x", "y", "z"))
    targets <- c("simplex_results", 
                 "ccm_func", "ccm_params", "ccm_results", 
                 "ccm_links", "smap_coeffs", "smap_matrices", 
                 "eigen_decomp", "eigenvalues", "eigenvectors")
    
    # check plan
    expect_error(my_plan <- build_dynamic_stability_plan(E_list = 3:5, 
                                                         surrogate_method = "random_shuffle", 
                                                         num_surr = 4, 
                                                         lib_sizes = seq(10, 100, 10), 
                                                         num_samples = 10), NA)
    expect_is(my_plan, c("drake_plan", "tbl_df"))
    expect_equal(vapply(targets, function(patt) {grep(patt, my_plan$target)}, 0), 
                 seq_along(targets), check.names = FALSE)
    expect_known_hash(as.character(my_plan$command), "1307b57d7a")
})

test_that("running full dynamic stability plan using block_3sp example data", {
    data("block_3sp", package = "rEDM")
    block <- setNames(block_3sp[, c("time", "x_t", "y_t", "z_t")], 
                      c("censusdate", "x", "y", "z"))
    targets <- c("simplex_results", 
                 "ccm_func", "ccm_params", "ccm_results", 
                 "ccm_links", "smap_coeffs", "smap_matrices", 
                 "eigen_decomp", "eigenvalues", "eigenvectors")
    
    # check plan
    my_plan <- build_dynamic_stability_plan(E_list = 3:5, 
                                            surrogate_method = "random_shuffle", 
                                            num_surr = 4, 
                                            lib_sizes = seq(10, 100, 10), 
                                            num_samples = 10)

    # run plan
    future::plan(future.callr::callr)
    drake::clean(destroy = TRUE)
    expect_error(drake::make(my_plan, seed = 42), NA)
    
    # inspect targets and check seed
    expect_error(drake::loadd(), NA)
    expect_true(all(targets %in% ls()))
    expect_equal(diagnose(simplex_results)$seed, 616382247)
    
    # check targets
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
