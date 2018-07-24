#### setup ----

library(tidyverse)
library(lubridate)
library(rEDM)
library(parallel)

# devtools::install_github("weecology/portalr")
# devtools::install_github("weecology/LDATS")

source(here::here("R/dynamic_stability_functions.R"))
source(here::here("R/plotting_functions.R"))

#### base run ----
if (FALSE) #  using all species in the dataset and default options
{
    
    block <- make_portal_block()
    saveRDS(block, here::here("data/portal_block.RDS"))
    
    results_file <- here::here("output/portal_ds_results.RDS")
    compute_dynamic_stability(block, results_file)
    
    results <- readRDS(results_file)
    plot_smap_coeffs(results$smap_matrices, here::here("figures/portal_smap_values.pdf"))
    plot_eigenvalues(results$eigenvalues, here::here("figures/portal_eigenvalues.pdf"))
}

#### revised run (50%) ----
if (FALSE) #  only species which appear at least 50% of the time)
{          #  DM, DO, DS, NA, OL, OT, PP
    
    block <- make_portal_block(filter_q = 0.5)
    saveRDS(block, here::here("data/portal_block_50.RDS"))
    
    results_file <- here::here("output/portal_ds_results_50.RDS")
    compute_dynamic_stability(block, results_file)
    
    results <- readRDS(results_file)
    plot_smap_coeffs(results$smap_matrices, here::here("figures/portal_smap_values_50.pdf"))
    plot_eigenvalues(results$eigenvalues, here::here("figures/portal_eigenvalues_50.pdf"))
    plot_eigenvalues(results$eigenvalues, here::here("figures/portal_eigenvalues_50_highlight.pdf"), highlight_shifts = TRUE)
    
    portal_network_50 <- results$ccm_results %>% 
        make_network_from_ccm_results()
    ggsave(here::here("figures/portal_network_50.pdf"),
           portal_network_50$plot, width = 8, height = 6)
    
    plot_eigenvalues(results$eigenvalues, num_values = 4)
}

#### revised run (33%) ----
if (FALSE) #  only species which appear at least 33% of the time)
{          #  DM, DO, DS, NA, OL, OT, PB, PE, PF, PP, RM
    
    block <- make_portal_block(filter_q = 1/3)
    saveRDS(block, here::here("data/portal_block_33.RDS"))
    
    results_file <- here::here("output/portal_ds_results_33.RDS")
    compute_dynamic_stability(block, results_file)
    
    results <- readRDS(results_file)
    plot_smap_coeffs(results$smap_matrices, here::here("figures/portal_smap_values_33.pdf"))
    plot_eigenvalues(results$eigenvalues, here::here("figures/portal_eigenvalues_33.pdf"))
    plot_eigenvalues(results$eigenvalues, here::here("figures/portal_eigenvalues_33_highlight.pdf"), highlight_shifts = TRUE)
}

#### run on LDA topics ----
if (FALSE) # do our own block generation using LDA topics
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
    
    results <- readRDS(results_file)
    plot_smap_coeffs(results$smap_matrices, here::here("figures/portal_smap_values_lda.pdf"))
    plot_eigenvalues(results$eigenvalues, here::here("figures/portal_eigenvalues_lda.pdf"))
    plot_eigenvalues(results$eigenvalues, here::here("figures/portal_eigenvalues_lda_highlight.pdf"), highlight_shifts = TRUE)
}

#### misc. graphs ----

make_combined_network(plot_file = here::here("figure/portal_networks_all.pdf"))

#### modified network runs ----
if (FALSE) # run with 50% species, but fully connected graph
{          #  DM, DO, DS, NA, OL, OT, PP
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

    results <- readRDS(results_file)
    plot_smap_coeffs(results$smap_matrices, here::here("figures/portal_smap_values_50_full_network.pdf"))
    plot_eigenvalues(results$eigenvalues, here::here("figures/portal_eigenvalues_50_full_network.pdf"))
    plot_eigenvalues(results$eigenvalues, here::here("figures/portal_eigenvalues_50_full_network_highlight.pdf"), highlight_shifts = TRUE)
}

#### revised run (50%, biomass) ----

if (FALSE)
{
    block <- make_portal_block(filter_q = 0.5, output = "biomass")
    saveRDS(block, here::here("data/portal_block_50_biomass.RDS"))
    
    results_file <- here::here("output/portal_ds_results_50_biomass.RDS")
    compute_dynamic_stability(block, results_file, rescale = FALSE)
    
    results <- readRDS(results_file)
    plot_smap_coeffs(results$smap_matrices, here::here("figures/portal_smap_values_50_biomass.pdf"))
    plot_eigenvalues(results$eigenvalues, here::here("figures/portal_eigenvalues_50_biomass.pdf"))
    plot_eigenvalues(results$eigenvalues, here::here("figures/portal_eigenvalues_50_biomass_highlight.pdf"), highlight_shifts = TRUE)
}