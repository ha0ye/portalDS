library(portalDS)

block <- make_portal_block(filter_q = 0.5)
saveRDS(block, here::here("data/portal_block_50.RDS"))

results_file <- here::here("output/portal_ds_results_50.RDS")
compute_dynamic_stability(block, results_file)

results <- readRDS(results_file)

portal_smap_coeffs <- plot_smap_coeffs(results$smap_matrices)
portal_eigenvalues <- plot_eigenvalues(results$eigenvalues)
portal_eigenvectors <- plot_eigenvectors(results$eigenvectors)

plot_eigenvalues(results$eigenvalues)
plot_eigenvectors(results$eigenvectors)

results$svd_decomp <- compute_svd_decomp(results$smap_matrices)

plot_svd_values(results$svd_decomp$d)
plot_svd_vectors(results$svd_decomp$u)
