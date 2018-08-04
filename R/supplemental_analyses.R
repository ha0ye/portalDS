#### run with species that are present at least 1/3 of the time ----
if (FALSE) #  DM, DO, DS, NA, OL, OT, PB, PE, PF, PP, RM
{
    block <- make_portal_block(filter_q = 1/3)
    saveRDS(block, here::here("data/portal_block_33.RDS"))
    
    results_file <- here::here("output/portal_ds_results_33.RDS")
    compute_dynamic_stability(block, results_file)
    
    plot_smap_coeffs(results_file, 
                     here::here("figures/portal_smap_values_33.pdf"))
    plot_eigenvalues(results_file, num_values = 2, 
                     here::here("figures/portal_eigenvalues_33.pdf"))
}

#### run on LDA topics ----
if (FALSE)
{
    block_file <- here::here("data/portal_block_lda.RDS")
    if (!file.exists(block_file))
    {
        raw_rodent_data <- portalr::abundance(time = "all", 
                                              plots = c(2, 4, 8, 11, 12, 14, 17, 22), 
                                              effort = TRUE, 
                                              na_drop = TRUE)
        
        # summarize by each newmmonnumber, and for only the control plots we want
        block <- raw_rodent_data %>% 
            filter(censusdate < "2015-04-18") %>% 
            select(-period)
        
        # check that effort is equal across samples
        stopifnot(length(unique(block$ntraps)) == 1)
        
        # generate LDA
        LDA_models <- LDATS::parLDA(data = select(block, BA:SO), ntopics = 2:5, nseeds = 200, ncores = 2)
        best_LDA <- LDATS::LDA_select(LDA_models, correction = TRUE)
        topic_block <- data.frame(best_LDA@gamma)
        names(topic_block) <- paste0("Topic", seq(NCOL(topic_block)))
        topic_block <- cbind(select(block, -(BA:SO)), topic_block)
        
        # add in NAs for unsampled newmoonnumbers and interpolate
        block <- topic_block %>%
            complete(newmoonnumber = full_seq(newmoonnumber, 1), fill = list(NA)) %>%
            mutate_at(vars(-newmoonnumber, -ntraps), forecast::na.interp) %>%
            mutate(censusdate = as.Date(as.numeric(censusdate), origin = "1970-01-01")) %>% 
            select(-newmoonnumber, -ntraps)
        
        saveRDS(block, file = block_file)
    }
    
    block <- readRDS(block_file)
    results_file <- here::here("output/portal_ds_results_lda.RDS")
    compute_dynamic_stability(block, results_file)
    
    plot_smap_coeffs(results_file, 
                     here::here("figures/portal_smap_values_lda.pdf"))
    plot_eigenvalues(results_file, num_values = 2, 
                     here::here("figures/portal_eigenvalues_lda.pdf"))
}


#### misc. graphs ----
if (FALSE) 
{
    make_combined_network(plot_file = here::here("figure/portal_networks_all.pdf"))
}

#### run with 50% species, but fully connected graph ----
if (FALSE)
{
    block <- readRDS("data/portal_block_50.RDS")
    
    results <- readRDS(here::here("output/portal_ds_results_50.RDS"))
    results$ccm_links <- tidyr::complete(results$ccm_links, lib_column, target_column, 
                                         fill = list(E = 1, delta_rho = NA, 
                                                     rho_minus_upper_q_null = NA)) %>%
        group_by(lib_column) %>%
        mutate(E = max(E)) %>%
        ungroup()
    
    results$smap_matrices <- NULL
    results$eigenvalues <- NULL
    
    results_file <- here::here("output/portal_ds_results_50_full_network.RDS")    
    saveRDS(results, results_file)
    compute_dynamic_stability(block, results_file)
    
    plot_smap_coeffs(results_file, 
                     here::here("figures/portal_smap_values_50_full_network.pdf"))
    plot_eigenvalues(results_file, num_values = 2, 
                     here::here("figures/portal_eigenvalues_50_full_network.pdf"))
    plot_eigenvectors(results_file, num_values = 2, 
                      here::here("figures/portal_eigenvectors_50_full_network.pdf"))
}

#### run with 50% species, but with biomass ----
if (FALSE)
{
    block <- make_portal_block(filter_q = 0.5, output = "biomass")
    saveRDS(block, here::here("data/portal_block_50_biomass.RDS"))
    
    results_file <- here::here("output/portal_ds_results_50_biomass.RDS")
    compute_dynamic_stability(block, results_file, rescale = FALSE)
    
    plot_smap_coeffs(results_file, 
                     here::here("figures/portal_smap_values_50_biomass.pdf"))
    plot_eigenvalues(results_file, num_values = 2, 
                     here::here("figures/portal_eigenvalues_50_biomass.pdf"))
    plot_eigenvectors(results_file, num_values = 2, 
                      here::here("figures/portal_eigenvectors_50_biomass.pdf"))
}

#### run with 50% species, but unscaled ----
if (FALSE)
{
    orig_results_file <- here::here("output/portal_ds_results_50.RDS")
    orig_results <- readRDS(orig_results_file)
    
    results_file <- here::here("output/portal_ds_results_50_unscaled.RDS")
    results <- orig_results[c("block", "simplex_results", "ccm_results")]
    saveRDS(results, file = results_file)
    
    block <- results$block
    compute_dynamic_stability(block, results_file, rescale = FALSE)
    
    plot_smap_coeffs(results_file, 
                     here::here("figures/portal_smap_values_50_unscaled.pdf"))
    plot_eigenvalues(results_file, num_values = 2, 
                     here::here("figures/portal_eigenvalues_50_unscaled.pdf"))
    plot_eigenvectors(results_file, num_values = 2, 
                      here::here("figures/portal_eigenvectors_50_unscaled.pdf"))
}

#### run with 50% species, but rolling s-map models forward----
if (FALSE)
{
    orig_results_file <- here::here("output/portal_ds_results_50.RDS")
    orig_results <- readRDS(orig_results_file)
    
    results_file <- here::here("output/portal_ds_results_50_rolling.RDS")
    results <- orig_results[c("block", "simplex_results", "ccm_results")]
    saveRDS(results, file = results_file)
    
    block <- results$block
    compute_dynamic_stability(block, results_file, rolling_forecast = TRUE)
    
    plot_smap_coeffs(results_file, 
                     here::here("figures/portal_smap_values_50_rolling.pdf"))
    plot_eigenvalues(results_file, num_values = 2, 
                     here::here("figures/portal_eigenvalues_50_rolling.pdf"))
    plot_eigenvectors(results_file, num_values = 2, 
                      here::here("figures/portal_eigenvectors_50_rolling.pdf"))
}

#### run with species that are present at least 1/2 of the time and PB ----
if (FALSE) #  DM, DO, DS, NA, OL, OT, PP, PB
{
    block <- make_portal_block(filter_q = 0.45) %>% 
        select(-PE)
    saveRDS(block, here::here("data/portal_block_50_pb.RDS"))
    
    results_file <- here::here("output/portal_ds_results_50_pb.RDS")
    compute_dynamic_stability(block, results_file)
    
    plot_smap_coeffs(results_file, 
                     here::here("figures/portal_smap_values_50_pb.pdf"))
    plot_eigenvalues(results_file, num_values = 2, 
                     here::here("figures/portal_eigenvalues_50_pb.pdf"))
    plot_eigenvectors(results_file, num_values = 2, 
                      here::here("figures/portal_eigenvectors_50_pb.pdf"))
}