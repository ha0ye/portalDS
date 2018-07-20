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
#  using all species in the dataset and default options

results_file <- here::here("output/portal_ds_results.RDS")

block <- make_portal_block()
compute_dynamic_stability(block, results_file)

results <- readRDS(results_file)
plot_smap_coeffs(results$smap_matrices)
plot_eigenvalues(results$eigenvalues)

#### revised run (50%) ----
#  only species which appear at least 50% of the time)
#  DM, DO, DS, NA, OL, OT, PP

results_file <- here::here("output/portal_ds_results_50.RDS")

block <- make_portal_block(filter_q = 0.5)
compute_dynamic_stability(block, results_file)

results <- readRDS(results_file)
plot_smap_coeffs(results$smap_matrices)
plot_eigenvalues(results$eigenvalues)
plot_eigenvalues(results$eigenvalues, highlight_shifts = TRUE)

#### revised run (33%) ----
#  only species which appear at least 33% of the time)
#  DM, DO, DS, NA, OL, OT, PB, PE, PF, PP, RM

results_file <- here::here("output/portal_ds_results_33.RDS")

block <- make_portal_block(filter_q = 1/3)
compute_dynamic_stability(block, results_file)

results <- readRDS(results_file)
plot_smap_coeffs(results$smap_matrices)
plot_eigenvalues(results$eigenvalues)
plot_eigenvalues(results$eigenvalues, highlight_shifts = TRUE)

#### run on LDA topics ----

results_file <- here::here("output/portal_ds_results_lda.RDS")

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
}
compute_dynamic_stability(block, results_file)

results <- readRDS(results_file)
plot_smap_coeffs(results$smap_matrices)
plot_eigenvalues(results$eigenvalues)
plot_eigenvalues(results$eigenvalues, highlight_shifts = TRUE)

#### make network figure ----
portal_network <- readRDS(here::here("output", "portal_ccm_results.RDS")) %>% 
    make_network_from_ccm_results()

portal_network_50 <- readRDS(here::here("output", "portal_ccm_results_50.RDS")) %>%
    make_network_from_ccm_results(palette = portal_network$palette, 
                                  existing_graph = portal_network$graph)

portal_network_33 <- readRDS(here::here("output", "portal_ccm_results_33.RDS")) %>%
    make_network_from_ccm_results(palette = portal_network$palette, 
                                  existing_graph = portal_network$graph)

combined_network_plot <- plot_grid(portal_network$plot, 
                                            portal_network_50$plot, 
                                            portal_network_33$plot, 
                                            labels = c("all species", 
                                                       "species present >= 50%", 
                                                       "species present >= 33%"), 
                                            align = "v", axis = "l", 
                                            ncol = 1)
# print(combined_network_plot)

ggsave(here::here("output", "portal_interaction_networks.pdf"),
       combined_network_plot, width = 6, height = 12)

portal_network_50 <- readRDS(here::here("output", "portal_ccm_results_50.RDS")) %>%
    make_network_from_ccm_results()
ggsave(here::here("output", "portal_interaction_network_50.pdf"),
       portal_network_50$plot, width = 8, height = 6)
