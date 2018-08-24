#### setup ----

library(tidyverse)
library(lubridate)
library(rEDM)
library(parallel)

# devtools::install_github("weecology/portalr")
# devtools::install_github("weecology/LDATS")

source(here::here("R/dynamic_stability_functions.R"))
source(here::here("R/plotting_functions.R"))

#### run using Maizuru Bay fish community ----
if (FALSE)
{
    block <- readRDS(here::here("data/maizuru_block.RDS"))
    results_file <- here::here("output/maizuru_ds_results.RDS")
    compute_dynamic_stability(block, results_file)
    results <- readRDS(results_file)
    
    ## plot time series ----
    maizuru_time_series <- plot_time_series(block) + 
        scale_x_date(limits = as.Date(c("2005-03-01", "2011-01-01")), 
                     date_breaks = "1 year", date_labels = "%Y", expand = c(0.01, 0)) + 
        scale_color_viridis(discrete = TRUE) + 
        guides(color = FALSE)
    ggsave(here::here("figures/maizuru_time_series.pdf"),
           maizuru_time_series, width = 9, height = 4.5)
    
    ## plot network ----
    my_graph <- results$ccm_links %>%
        filter(lib_column != target_column) %>% 
        arrange(target_column) %>%
        select(target_column, lib_column) %>%
        graph_from_data_frame(vertices = levels(.$target_column))
    
    palette <- viridis(length(V(my_graph)))
    names(palette) <- as_ids(V(my_graph))
    
    my_graph <- create_layout(my_graph, layout = layout)
    
    maizuru_network <- ggraph(my_graph) + 
        geom_edge_fan(edge_width = 1, start_cap = circle(0.3, "inches"), 
                      end_cap = circle(0.3, "inches"),
                      arrow = arrow(angle = 20, type = "closed", 
                                    length = unit(0.15, "inches"))) + 
        geom_node_circle(aes(r = 0.08, fill = name)) + 
        theme_graph(foreground = "black", fg_text_colour = "white", 
                    background = "transparent") + 
        coord_fixed() + 
        theme(text = element_text(family = "Helvetica"),
              panel.border = element_rect(color = NA, fill = NA)) + 
        scale_fill_manual(values = palette) + 
        guides(fill = guide_legend(title = "Species"))
    ggsave(here::here("figures/maizuru_network.pdf"),
           maizuru_network, width = 8, height = 6)
    
    ## plot subset of time series ----
    sp_list <- c("Rudarius.ercodes", "Trachurus.japonicus")
    palette <- viridis(NCOL(block)-1)[match(sp_list, sort(colnames(select(block, -censusdate))))]
    maizuru_time_series_subset <- 
        select(block, c("censusdate", sp_list)) %>%
        plot_time_series() + 
        scale_x_date(limits = as.Date(c("2005-03-01", "2011-01-01")), 
                     date_breaks = "1 year", date_labels = "%Y", expand = c(0.01, 0)) + 
        scale_color_manual(values = palette) + 
        guides(color = FALSE) + 
        facet_wrap(~species, ncol = 1)
    
    ggsave(here::here("figures/maizuru_time_series_subset.pdf"),
           maizuru_time_series_subset, width = 9, height = 5)
    
    ## plot subset of network ----
    existing_graph <- my_graph # save old graph layout
    
    my_graph <- results$ccm_links %>%
        filter(lib_column != target_column) %>% 
        arrange(target_column) %>%
        select(target_column, lib_column) %>%
        filter(target_column %in% sp_list, lib_column %in% sp_list) %>%
        graph_from_data_frame()
    
    my_graph <- create_layout(my_graph, layout = layout)
    idx <- match(my_graph$name, existing_graph$name)
    my_graph$x <- existing_graph$x[idx]
    my_graph$y <- existing_graph$y[idx]
    
    maizuru_network_subset <- ggraph(my_graph) + 
        geom_edge_fan(edge_width = 1, start_cap = circle(0.3, "inches"), 
                      end_cap = circle(0.3, "inches"),
                      arrow = arrow(angle = 20, type = "closed", 
                                    length = unit(0.15, "inches"))) + 
        geom_node_circle(aes(r = 0.08, fill = name)) + 
        theme_graph(foreground = "black", fg_text_colour = "white", 
                    background = "transparent") + 
        coord_fixed() + 
        theme(text = element_text(family = "Helvetica"),
              panel.border = element_rect(color = NA, fill = NA)) + 
        scale_fill_manual(values = palette) + 
        guides(fill = guide_legend(title = "Species"))
    ggsave(here::here("figures/maizuru_network_subset.pdf"),
           maizuru_network_subset, width = 8, height = 6)
    
    ## plot eigenvalues ----
    ds_plot <- results$eigenvalues %>%
        plot_eigenvalues() + 
        scale_color_manual(values = NA) + 
        scale_x_date(limits = as.Date(c("2005-03-01", "2011-01-01")), 
                     date_breaks = "1 year", date_labels = "%Y", expand = c(0.01, 0))
        
    ggsave(here::here("figures/maizuru_ds_blank.pdf"),
           ds_plot, width = 10, height = 5)
    ggsave(here::here("figures/maizuru_ds_threshold.pdf"),
           ds_plot + geom_hline(yintercept = 1.0, size = 1, linetype = 2), 
           width = 10, height = 5)
    ggsave(here::here("figures/maizuru_ds.pdf"),
           ds_plot + scale_color_manual(values = "black") + 
               geom_hline(yintercept = 1.0, size = 1, linetype = 2), 
           width = 10, height = 5)
}
