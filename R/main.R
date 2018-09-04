#### setup ----

library(tidyverse)
library(lubridate)
library(rEDM)
library(parallel)

# devtools::install_github("weecology/portalr")
# devtools::install_github("weecology/LDATS")

source(here::here("R/dynamic_stability_functions.R"))
source(here::here("R/plotting_functions.R"))

#### default run using all species in the dataset ----
if (FALSE)
{
    block <- make_portal_block()
    saveRDS(block, here::here("data/portal_block.RDS"))
    
    results_file <- here::here("output/portal_ds_results.RDS")
    compute_dynamic_stability(block, results_file)
}

#### run with species that are present at least 1/2 of the time ----
if (FALSE) #  DM, DO, DS, NA, OL, OT, PP
{
    block <- make_portal_block(filter_q = 0.5)
    saveRDS(block, here::here("data/portal_block_50.RDS"))
    
    results_file <- here::here("output/portal_ds_results_50.RDS")
    compute_dynamic_stability(block, results_file)
    
    portal_network_50 <- make_network_from_ccm_results(results_file)
    ggsave(here::here("figures/portal_network_50.pdf"),
           portal_network_50$plot, width = 8, height = 6)
}

#### run with species that are present at least 1/2 of the time (K-rat exclosures) ----
if (FALSE) #  DM, NA, OT, PE, PF, PP, RM
{
    block <- make_portal_block(filter_q = 0.5, 
                               plots = c(3, 6, 13, 15, 18, 19, 20, 21))
    saveRDS(block, here::here("data/portal_block_50_krr.RDS"))
    
    results_file <- here::here("output/portal_ds_results_50_krr.RDS")
    compute_dynamic_stability(block, results_file)
}