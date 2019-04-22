library(portalDS)

expose_imports("portalDS")

data("block_3sp", package = "rEDM")
block <- block_3sp[, c("time", "x_t", "y_t", "z_t")]
var_names <- c("x", "y", "z")
names(block) <- c("censusdate", var_names)

my_plan <- build_dynamic_stability_plan(block, 
                                        E_list = 3:5, 
                                        surrogate_method = "random_shuffle", 
                                        num_surr = 4, 
                                        lib_sizes = seq(10, 100, by = 10), 
                                        num_samples = 10)

drake::make(my_plan)
