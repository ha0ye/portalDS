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
    
    plot_smap_coeffs(results_file, 
                     here::here("figures/portal_smap_values.pdf"))
    plot_eigenvalues(results_file, 
                     here::here("figures/portal_eigenvalues.pdf"))
    plot_eigenvectors(results_file, 
                     here::here("figures/portal_eigenvectors.pdf"))
}

#### run with species that are present at least 1/2 of the time ----
if (FALSE) #  DM, DO, DS, NA, OL, OT, PP
{
    block <- make_portal_block(filter_q = 0.5)
    saveRDS(block, here::here("data/portal_block_50.RDS"))
    
    results_file <- here::here("output/portal_ds_results_50.RDS")
    compute_dynamic_stability(block, results_file)
    
    plot_smap_coeffs(results_file, 
                     here::here("figures/portal_smap_values_50.pdf"))
    plot_eigenvalues(results_file, num_values = 2, 
                     here::here("figures/portal_eigenvalues_50.pdf"))
    plot_eigenvectors(results_file, num_values = 2, 
                      here::here("figures/portal_eigenvectors_50.pdf"))
    
    portal_network_50 <- make_network_from_ccm_results(results_file)
    ggsave(here::here("figures/portal_network_50.pdf"),
           portal_network_50$plot, width = 8, height = 6)
}

block <- readRDS(here::here("data/portal_block_50.RDS"))

## plot time series ----
to_plot <- block %>%
    gather(species, abundance, -censusdate)
n <- length(unique(to_plot$species))
portal_time_series <- to_plot %>%
    ggplot(aes(x = censusdate, y = abundance, color = species)) + 
    geom_line(size = 1) + 
    scale_x_date(limits = as.Date(c("2005-03-01", "2011-01-01")), 
                 date_breaks = "1 year", date_labels = "%Y", expand = c(0.01, 0)) + 
    scale_color_manual(values = viridis(n)) + 
    labs(x = NULL, y = "relative abundance") + 
#    guides(color = FALSE) + 
    theme_bw(base_size = 20, base_family = "Helvetica", 
             base_line_size = 1) + 
    theme(panel.grid.minor = element_line(size = 1))
ggsave(here::here("figures/portal_time_series.pdf"),
       portal_time_series, width = 9, height = 4.5)