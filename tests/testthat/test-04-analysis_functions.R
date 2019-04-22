context("Check dynamic stability analysis functions work")

data_path <- system.file("extdata", "portal_block.RDS",
                         package = "portalDS", mustWork = TRUE)
portal_block <- readRDS(data_path)

test_that("compute_simplex works as expected", {
    var_names <- setdiff(names(portal_block), "censusdate")
    E_list <- 3:5
    expect_error(simplex_out <- compute_simplex(portal_block, 
                                                E_list = E_list, 
                                                num_surr = 4), NA)
    # check columns
    expect_setequal(c("species", "data", "simplex_out", "best_E", "surrogate_data"), names(simplex_out))
    
    # check species values
    expect_setequal(simplex_out$species, var_names)
    
    # check random data value
    idx <- sample(NROW(simplex_out), 1)
    simplex_out_data <- simplex_out[[idx, "data"]][[1]]
    original_portal_data <- portal_block[[simplex_out$species[idx]]]
    expect_equal(simplex_out_data, original_portal_data)
    
    # check simplex_out
    simplex_result <- simplex_out[[idx, "simplex_out"]]
    expect_equal(simplex_result$E, E_list)
    simplex_results <- simplex_out$simplex_out
    expect_known_hash(lapply(simplex_results, round, 4), "cbd599464a")
    
    # check best_E
    best_E <- simplex_out[["best_E"]]
    expect_true(all(best_E >= min(E_list)))
    expect_true(all(best_E <= max(E_list)))
    
    # check surrogate_data
    expect_error(surr_data <- simplex_out[[idx, "surrogate_data"]], NA)
    expect_equal(dim(surr_data), c(NROW(portal_block), 4))
    expect_equal(vapply(simplex_out$surrogate_data, dim, c(0, 0)), 
                 matrix(c(NROW(portal_block), 4), ncol = length(var_names), nrow = 2))
})

test_that("identify_ccm_links works as expected", {
    data(maizuru_block)
    var_names <- setdiff(names(maizuru_block), "censusdate")
    
    data_path <- system.file("extdata", "maizuru_ccm_results.RDS",
                             package = "portalDS", mustWork = TRUE)
    maizuru_ccm_results <- readRDS(data_path)
    
    expect_error(ccm_links <- identify_ccm_links(maizuru_ccm_results), NA)
    expect_equal(dim(ccm_links), c(32, 5))
    expect_setequal(var_names, as.character(ccm_links$lib_column))
    expect_setequal(var_names,  as.character(ccm_links$target_column))
    expect_known_hash(ccm_links, "dd99bafb30")
})

test_that("compute_smap_coeffs and compute_smap_matrices work as expected", {
    data(maizuru_block)
    var_names <- setdiff(names(maizuru_block), "censusdate")
    
    data_path <- system.file("extdata", "maizuru_ccm_links.RDS",
                             package = "portalDS", mustWork = TRUE)
    maizuru_ccm_links <- readRDS(data_path)
    expect_error(smap_coeffs <- compute_smap_coeffs(maizuru_block, maizuru_ccm_links), NA)
    smap_coeffs <- lapply(smap_coeffs, round, 4)
    expect_known_hash(smap_coeffs, "1c7a16172f")

    expect_error(smap_matrices <- compute_smap_matrices(smap_coeffs, maizuru_ccm_links), NA)
    expect_is(smap_matrices, "list")
    expect_equal(length(smap_matrices), NROW(maizuru_block))
    expect_error(matrix_sizes <- lapply(smap_matrices, dim), NA)
    expect_error(idx <- vapply(matrix_sizes, is.null, TRUE), NA)
    expect_equal(do.call(rbind, matrix_sizes[!idx]), 
                 matrix(312, nrow = sum(!idx), ncol = 2))
    expect_known_hash(smap_matrices, "ca491bb58d")
})


## check remaining dynamic stability workflow functions
# data_path <- system.file("extdata", "maizuru_smap_matrices.RDS",
#                          package = "portalDS", mustWork = TRUE)
# maizuru_smap_matrices <- readRDS(data_path)
