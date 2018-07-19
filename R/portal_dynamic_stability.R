#### setup ----

library(tidyverse)
library(lubridate)
library(rEDM)
library(parallel)
# devtools::install_github("weecology/portalr")
# devtools::install_github("weecology/LDATS")

source(here::here("R", "dynamic_stability_functions.R"))
source(here::here("R", "plotting_functions.R"))

#### function for full dynamic stability analysis ----

compute_portal_dynamic_stability <- function(proc_data_file = "data/portal_block.RDS", 
                                             simplex_results_file = "output/portal_simplex_results.RDS", 
                                             ccm_results_file = "output/portal_ccm_results.RDS", 
                                             smap_matrices_file = "output/portal_smap_matrices.RDS", 
                                             eigenvalues_file = "output/portal_eigenvalues.RDS", 
                                             filter_q = NULL, 
                                             file_suffix = NULL, 
                                             num_cores = 2)
{

    if (is.null(file_suffix))
    {
        file_suffix <- ""
        if (!is.null(filter_q))
            file_suffix <- paste0("_", as.character(floor(filter_q * 100)))
    }
    proc_data_file <- sub(".RDS", paste0(file_suffix, ".RDS"), proc_data_file)
    simplex_results_file <- sub(".RDS", paste0(file_suffix, ".RDS"), simplex_results_file)
    ccm_results_file <- sub(".RDS", paste0(file_suffix, ".RDS"), ccm_results_file)
    smap_matrices_file <- sub(".RDS", paste0(file_suffix, ".RDS"), smap_matrices_file)
    eigenvalues_file <- sub(".RDS", paste0(file_suffix, ".RDS"), eigenvalues_file)
    
    # pre-processing
    if (!file.exists(proc_data_file))
    {
        # options here are:
        #   time = "all" (allow us to do correct interpolation and accouting for seasonality)
        #   level = "plot" (allow us to pull out abundances on the plots we want)
        #   effort = TRUE (so we can check effort)
        #   na_drop = TRUE (ignore periods where sampling did not occur)
        raw_rodent_data <- portalr::abundance(time = "all", 
                                              level = "plot", 
                                              effort = TRUE, 
                                              na_drop = TRUE)
        
        # summarize by each newmmonnumber, and for only the control plots we want
        block <- raw_rodent_data %>% 
            filter(plot %in% c(2, 4, 8, 11, 12, 14, 17, 22), 
                   censusdate < "2015-04-18") %>% 
            select(-treatment, -plot, -period) %>%
            group_by(newmoonnumber, censusdate) %>%
            summarize_all(sum) %>%
            ungroup()
        
        # check that effort is equal across samples
        stopifnot(length(unique(block$ntraps)) == 1)
        
        if (!is.null(filter_q))
        {
            species_list <- block %>%
                gather(species, abundance, BA:SO) %>%
                group_by(species) %>%
                summarize(quantile_q = quantile(abundance, 1 - filter_q)) %>%
                filter(quantile_q > 0) %>%
                pull(species)
            
            block <- block %>% 
                select(c("newmoonnumber", "censusdate", "ntraps", species_list))
        }
        
        # add in NAs for unsampled newmoonnumbers and interpolate
        block <- block %>%
            complete(newmoonnumber = full_seq(newmoonnumber, 1), fill = list(NA)) %>%
            mutate_at(vars(-newmoonnumber, -ntraps), forecast::na.interp) %>%
            mutate(censusdate = as.Date(as.numeric(censusdate), origin = "1970-01-01")) %>% 
            select(-newmoonnumber, -ntraps)
        
        saveRDS(block, proc_data_file)
    }
    block <- readRDS(proc_data_file)
    block_dates <- block$censusdate
    block_yearday <- yday(block$censusdate)
    block <- block %>% select(-censusdate) # drop the time column
    
    # identify best embedding dimension and generate surrogate data
    if (!file.exists(simplex_results_file))
    {
        simplex_results <- compute_simplex(block, E = 1:16)
        
        simplex_results$surrogate_data <- 
            map(seq(NROW(simplex_results)), function(i) {
                make_surrogate_annual_spline(block_yearday, 
                                             simplex_results[[i, "data"]]$abundance, 
                                             num_surr = 200)
            })
        saveRDS(simplex_results, simplex_results_file)
    }
    simplex_results <- readRDS(simplex_results_file)
    
    if (!file.exists(ccm_results_file))
    {
        max_ts_length <- max(map_int(simplex_results$data, NROW))
        lib_vec <- c(6, 12, 24, 40, 80, 140, 220, 320, max_ts_length)
        ccm_results <- compute_CCM(simplex_results, lib_sizes = lib_vec, 
                                   num_samples = 200, num_cores = num_cores)
        saveRDS(ccm_results, ccm_results_file)
    }
    ccm_results <- readRDS(ccm_results_file)
    
    if (!file.exists(smap_matrices_file))
    {
        ccm_links <- identify_ccm_links(ccm_results)
        smap_coeffs <- compute_smap_coeffs(block, ccm_links)
        smap_matrices <- compute_smap_matrices(smap_coeffs, ccm_links)
        stopifnot(length(smap_matrices) == length(block_dates))
        names(smap_matrices) <- block_dates
        saveRDS(smap_matrices, smap_matrices_file)
    }
    smap_matrices <- readRDS(smap_matrices_file)
    
    if (!file.exists(eigenvalues_file))
    {
        eigenvalues <- compute_eigenvalues(smap_matrices)
        saveRDS(eigenvalues, eigenvalues_file)
    }
}

#### base run ----
#  using all species in the dataset and default options
proc_data_file <- here::here("data", "portal_block.RDS")
simplex_results_file <- here::here("output", "portal_simplex_results.RDS")
ccm_results_file <- here::here("output", "portal_ccm_results.RDS")
smap_matrices_file <- here::here("output", "portal_smap_matrices.RDS")
smap_plot_file <- here::here("output", "portal_smap_values.pdf")
eigenvalues_file <- here::here("output", "portal_eigenvalues.RDS")
eigenvalues_plot_file <- here::here("output", "portal_dynamic_stability.pdf")

compute_portal_dynamic_stability(proc_data_file, 
                                 simplex_results_file, 
                                 ccm_results_file, 
                                 smap_matrices_file, 
                                 eigenvalues_file)

plot_eigenvalues(eigenvalues_file, 
                 eigenvalues_plot_file)

plot_smap_coeffs(smap_matrices_file, 
                 smap_plot_file)

#### revised run (50%) ----
#  only species which appear at least 50% of the time)
#  DM, DO, DS, NA, OL, OT, PP

compute_portal_dynamic_stability(filter_q = 0.5, 
                                 num_cores = 1)

proc_data_file <- here::here("data", "portal_block_50.RDS")
smap_matrices_file <- here::here("output", "portal_smap_matrices_50.RDS")
smap_plot_file <- here::here("output", "portal_smap_values_50.pdf")
eigenvalues_file <- here::here("output", "portal_eigenvalues_50.RDS")
eigenvalues_plot_file <- here::here("output", "portal_dynamic_stability_50.pdf")

plot_eigenvalues(eigenvalues_file, 
                 eigenvalues_plot_file)

plot_eigenvalues(eigenvalues_file, 
                 plot_file = here::here("output", "portal_dynamic_stability_50_highlight.pdf"), 
                 highlight_shifts = TRUE)

plot_smap_coeffs(smap_matrices_file, 
                 smap_plot_file)

#### revised run (33%) ----
#  only species which appear at least 33% of the time)
#  DM, DO, DS, NA, OL, OT, PB, PE, PF, PP, RM

compute_portal_dynamic_stability(filter_q = 1/3,
                                 num_cores = 1)

proc_data_file <- here::here("data", "portal_block_33.RDS")
smap_matrices_file <- here::here("output", "portal_smap_matrices_33.RDS")
smap_plot_file <- here::here("output", "portal_smap_values_33.pdf")
eigenvalues_file <- here::here("output", "portal_eigenvalues_33.RDS")
eigenvalues_plot_file <- here::here("output", "portal_dynamic_stability_33.pdf")

plot_eigenvalues(eigenvalues_file, 
                 eigenvalues_plot_file)

plot_smap_coeffs(smap_matrices_file, 
                 smap_plot_file)

#### run on LDA topics ----

proc_data_file <- here::here("data", "portal_block_lda.RDS")
smap_matrices_file <- here::here("output", "portal_smap_matrices_lda.RDS")
smap_plot_file <- here::here("output", "portal_smap_values_lda.pdf")
eigenvalues_file <- here::here("output", "portal_eigenvalues_lda.RDS")
plot_file <- here::here("output", "portal_dynamic_stability_lda.pdf")

# do our own block generation using LDA topics
if (!file.exists(proc_data_file))
{
    raw_rodent_data <- portalr::abundance(time = "all", 
                                          level = "plot", 
                                          effort = TRUE, 
                                          na_drop = TRUE)
    
    # summarize by each newmmonnumber, and for only the control plots we want
    block <- raw_rodent_data %>% 
        filter(plot %in% c(2, 4, 8, 11, 12, 14, 17, 22), 
               censusdate < "2015-04-18") %>% 
        select(-treatment, -plot, -period) %>%
        group_by(newmoonnumber, censusdate) %>%
        summarize_all(sum) %>%
        ungroup()
    
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
    
    saveRDS(block, proc_data_file)
}

compute_portal_dynamic_stability(file_suffix = "_lda")

plot_eigenvalues(eigenvalues_file, 
                 eigenvalues_plot_file)

plot_smap_coeffs(smap_matrices_file, 
                 smap_plot_file)

