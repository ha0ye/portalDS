context("Plotting function")

# vdiffr::manage_cases(filter = "plotting")

data(maizuru_block)
var_names <- setdiff(names(maizuru_block), "censusdate")

test_that("check plotting of time series", {
  expect_error(ts_plot <- plot_time_series(maizuru_block), NA)

  lower <- seq.Date(as.Date("2002-12-01"), as.Date("2013-12-01"), "1 year")
  upper <- seq.Date(as.Date("2003-03-01"), as.Date("2014-03-01"), "1 year")
  expect_error(ts_plot_winters <- add_regime_shift_highlight(ts_plot, lower, upper), NA)

  skip_on_travis()
  vdiffr::expect_doppelganger("Maizuru Time Series", ts_plot)
  vdiffr::expect_doppelganger("Maizuru Time Series (Winters)", ts_plot_winters)
})

test_that("check plotting of the interaction network", {
  data_path <- system.file("extdata", "maizuru_ccm_links.RDS",
    package = "portalDS", mustWork = TRUE
  )
  maizuru_ccm_links <- readRDS(data_path)

  expect_error(maizuru_network <- plot_network(maizuru_ccm_links), NA)

  expect_error(palette <- maizuru_network$palette, NA)
  expect_setequal(names(palette), var_names)

  expect_error(graph <- maizuru_network$graph, NA)
  expect_s3_class(graph, c("layour_tbl_graph", "layout_ggraph", "data.frame"))
  expect_equal(dim(graph), c(length(var_names), 6))

  skip_on_travis()
  vdiffr::expect_doppelganger("Maizuru Interaction Network", maizuru_network$plot)
})

data_path <- system.file("extdata", "maizuru_smap_matrices.RDS",
  package = "portalDS", mustWork = TRUE
)
maizuru_smap_matrices <- readRDS(data_path)

test_that("check plotting of the smap coefficients", {
  expect_error(maizuru_smap_coeffs <- plot_smap_coeffs(maizuru_smap_matrices), NA)

  skip_on_travis()
  vdiffr::expect_doppelganger("Maizuru S-map Coefficients", maizuru_smap_coeffs)
})

test_that("check plotting of eigenvalues and eigenvectors", {
  eigen_decomp <- compute_eigen_decomp(maizuru_smap_matrices)

  maizuru_eigenvalues <- eigen_decomp$values
  expect_error(eigenvalue_plot <- plot_eigenvalues(maizuru_eigenvalues), NA)

  expect_error(eigenvalue_plot_complex <- plot_eigenvalues(maizuru_eigenvalues,
    highlight_complex = TRUE,
    num_values = 2
  ), NA)

  maizuru_eigenvectors <- eigen_decomp$vectors
  expect_error(eigenvector_plot <- plot_eigenvectors(maizuru_eigenvectors,
    num_values = 2
  ), NA)

  expect_error(eigenvector_plot_IPR <- plot_eigenvectors(maizuru_eigenvectors,
    num_values = 2,
    add_IPR = TRUE
  ), NA)

  skip_on_travis()
  vdiffr::expect_doppelganger("Maizuru Eigenvalues", eigenvalue_plot)
  vdiffr::expect_doppelganger("Maizuru Eigenvalues (complex)", eigenvalue_plot_complex)
  vdiffr::expect_doppelganger("Maizuru Eigenvectors", eigenvector_plot)
  vdiffr::expect_doppelganger("Maizuru Eigenvectors (IPR)", eigenvector_plot_IPR)
})
