context("Check dynamic stability analysis functions work")

data(maizuru_block)
var_names <- setdiff(names(maizuru_block), "censusdate")

test_that("compute_simplex works as expected", {
    E_list <- 3:5
    num_surr <- 4
    expect_error(simplex_out <- compute_simplex(maizuru_block, 
                                                E_list = E_list, 
                                                num_surr = num_surr), NA)
    # check columns
    expect_setequal(c("species", "data", "simplex_out", "best_E", "surrogate_data"), names(simplex_out))
    
    # check species values
    expect_setequal(simplex_out$species, var_names)
    
    # check random data value
    idx <- sample(NROW(simplex_out), 1)
    simplex_out_data <- simplex_out[[idx, "data"]][[1]]
    original_data <- maizuru_block[[simplex_out$species[idx]]]
    expect_equal(simplex_out_data, original_data)
    
    # check simplex_out
    simplex_result <- simplex_out[[idx, "simplex_out"]]
    expect_equal(simplex_result$E, E_list)
    simplex_results <- simplex_out$simplex_out
    expect_known_hash(lapply(simplex_results, round, 4), "5073eb853f")
    
    # check best_E
    best_E <- simplex_out[["best_E"]]
    expect_true(all(best_E >= min(E_list)))
    expect_true(all(best_E <= max(E_list)))
    
    # check surrogate_data
    expect_error(surr_data <- simplex_out[[idx, "surrogate_data"]], NA)
    expect_equal(dim(surr_data), c(NROW(maizuru_block), num_surr))
    expect_equal(vapply(simplex_out$surrogate_data, dim, c(0, 0)), 
                 matrix(c(NROW(maizuru_block), num_surr), ncol = length(var_names), nrow = 2))
})

test_that("identify_ccm_links works as expected", {
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
    data_path <- system.file("extdata", "maizuru_ccm_links.RDS",
                             package = "portalDS", mustWork = TRUE)
    maizuru_ccm_links <- readRDS(data_path)
    
    expect_error(smap_coeffs <- compute_smap_coeffs(maizuru_block, maizuru_ccm_links), NA)
    expect_equal(length(smap_coeffs), length(var_names))
    best_E <- maizuru_ccm_links %>% 
        dplyr::select(lib_column, E) %>% 
        dplyr::distinct() %>% 
        dplyr::pull(E)
    expect_equal(vapply(smap_coeffs, dim, c(0, 0)), 
                 matrix(c(rep.int(NROW(maizuru_block), length(var_names)), 
                          best_E + 1), 
                        byrow = TRUE, nrow = 2, 
                        dimnames = list(NULL, sort(var_names))))
    smap_coeffs <- lapply(smap_coeffs, round, 4)
    expect_known_hash(smap_coeffs, "4f48e52412")
    
    expect_error(smap_matrices <- compute_smap_matrices(smap_coeffs, maizuru_ccm_links), NA)
    expect_is(smap_matrices, "list")
    expect_equal(length(smap_matrices), NROW(maizuru_block))
    expect_error(matrix_sizes <- lapply(smap_matrices, dim), NA)
    expect_error(idx <- vapply(matrix_sizes, is.null, TRUE), NA)
    expect_equal(do.call(rbind, matrix_sizes[!idx]), 
                 matrix(312, nrow = sum(!idx), ncol = 2))
    expect_known_hash(smap_matrices, "51e878ba99")
})

test_that("compute_smap_coeffs and compute_smap_matrices work as expected", {
    data_path <- system.file("extdata", "maizuru_smap_matrices.RDS",
                             package = "portalDS", mustWork = TRUE)
    maizuru_smap_matrices <- readRDS(data_path)
    
    expect_error(eigen_decomp <- compute_eigen_decomp(maizuru_smap_matrices), NA)
    expect_is(eigen_decomp, "list")
    expect_equal(names(eigen_decomp), c("values", "vectors"))
    
    expect_error(maizuru_eigenvalues <- eigen_decomp$values, NA)
    expect_is(maizuru_eigenvalues, "list")
    expect_equal(length(maizuru_eigenvalues), NROW(maizuru_block))
    expect_error(idx <- vapply(maizuru_eigenvalues, anyNA, TRUE), NA)
    expect_error(eigenvalue_matrix <- do.call(rbind, maizuru_eigenvalues[!idx]), NA)
    expect_equal(dim(eigenvalue_matrix), c(sum(!idx), 312))
    expect_known_hash(eigenvalue_matrix, "cbf145998b")
    
    expect_error(maizuru_eigenvectors <- eigen_decomp$vectors, NA)
    expect_is(maizuru_eigenvectors, "list")
    expect_equal(length(maizuru_eigenvectors), NROW(maizuru_block))
    expect_error(idx <- vapply(maizuru_eigenvectors, is.null, TRUE), NA)
    expect_equal(vapply(maizuru_eigenvectors[!idx], dim, c(0, 0)), 
                 matrix(312, nrow = 2, ncol = sum(!idx), 
                        dimnames = list(NULL, names(maizuru_eigenvectors[!idx]))))
    expect_known_hash(maizuru_eigenvectors[!idx], "8d49244739")
    
    expect_equal(names(maizuru_eigenvalues), 
                 names(maizuru_eigenvectors))
})