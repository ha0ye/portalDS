library(portalDS)
library(drake)

expose_imports("portalDS")

data("block_3sp", package = "rEDM")

simplex_plan <- drake_plan(
  block = setNames(
    block_3sp[, c("time", "x_t", "y_t", "z_t")],
    c("censusdate", "x", "y", "z")
  ),
  simplex_results = compute_simplex(
    block = block,
    E_list = 3:5,
    surrogate_method = "random_shuffle",
    num_surr = 20
  )
)

my_plan <- dplyr::bind_rows(
  simplex_plan,
  drake_plan(ccm_results = compute_ccm(simplex_results,
    lib_sizes = seq(10, 100, by = 10),
    random_libs = TRUE, num_samples = 10,
    replace = TRUE, RNGseed = 42,
    silent = TRUE
  ))
)

my_plan_2 <- dplyr::bind_rows(
  simplex_plan,
  build_ccm_plan(
    lib_sizes = seq(10, 100, by = 10),
    random_libs = TRUE, num_samples = 10,
    replace = TRUE, RNGseed = 42,
    silent = TRUE
  )
)

future::plan(NULL)
tictoc::tic()
clean()
make(my_plan)
tictoc::toc()

future::plan(future::multiprocess)
tictoc::tic()
clean()
make(my_plan)
tictoc::toc()

future::plan(future.callr::callr)
tictoc::tic()
clean()
make(my_plan)
tictoc::toc()

future::plan(future.callr::callr)
tictoc::tic()
clean()
make(my_plan_2)
tictoc::toc()
