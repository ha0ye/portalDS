context("Plotting function")

# vdiffr::manage_cases(filter = "plotting")

data(maizuru_block)
var_names <- setdiff(names(maizuru_block), "censusdate")

test_that("check plotting of the interaction network", {
    data_path <- system.file("extdata", "maizuru_ccm_links.RDS",
                             package = "portalDS", mustWork = TRUE)
    maizuru_ccm_links <- readRDS(data_path)
    
    expect_error(maizuru_network <- plot_network(maizuru_ccm_links), NA)
    
    expect_error(palette <- maizuru_network$palette, NA)
    expect_setequal(names(palette), var_names)
    
    expect_error(graph <- maizuru_network$graph, NA)
    expect_s3_class(graph, c("layour_tbl_graph", "layout_ggraph", "data.frame"))
    expect_equal(dim(graph), c(length(var_names), 6))
    
    vdiffr::expect_doppelganger("Maizuru Interaction Network", maizuru_network$plot)
})

test_that("check plotting of the smap coefficients", {
    data_path <- system.file("extdata", "maizuru_smap_matrices.RDS",
                             package = "portalDS", mustWork = TRUE)
    maizuru_smap_matrices <- readRDS(data_path)
    
    expect_error(maizuru_smap_coeffs <- plot_smap_coeffs(maizuru_smap_matrices), NA)
    vdiffr::expect_doppelganger("Maizuru S-map Coefficients", maizuru_smap_coeffs)
})
