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
    to_plot <- block %>%
        mutate_at(vars(-censusdate), 
                  function(x) {(x - min(x, na.rm = TRUE)) / 
                          (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))}) %>%
        gather(species, abundance, -censusdate)
    n <- length(unique(to_plot$species))
    maizuru_time_series <- to_plot %>%
        ggplot(aes(x = censusdate, y = abundance, color = species)) + 
        geom_line(size = 1) + 
        scale_x_date(limits = as.Date(c("2005-03-01", "2011-01-01")), 
                     date_breaks = "1 year", date_labels = "%Y", expand = c(0.01, 0)) + 
        scale_color_manual(values = viridis(n)) + 
        labs(x = NULL, y = "relative abundance") + 
        guides(color = FALSE) + 
        theme_bw(base_size = 20, base_family = "Helvetica", 
                 base_line_size = 1) + 
        theme(panel.grid.minor = element_line(size = 1))
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
    maizuru_time_series_subset <- to_plot %>% 
        filter(species %in% sp_list) %>%
        ggplot(aes(x = censusdate, y = abundance, color = species)) + 
        facet_wrap(~species, ncol = 1) + 
        geom_line(size = 1) + 
        scale_x_date(limits = as.Date(c("2005-03-01", "2011-01-01")), 
                     date_breaks = "1 year", date_labels = "%Y", expand = c(0.01, 0)) + 
        scale_color_manual(values = viridis(n)[match(sp_list, sort(unique(to_plot$species)))]) +  
        labs(x = NULL, y = "relative abundance") + 
        guides(color = FALSE) + 
        theme_bw(base_size = 20, base_family = "Helvetica", 
                 base_line_size = 1) + 
        theme(panel.grid.minor = element_line(size = 1))
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
    eigenvalues <- results$eigenvalues
    num_values <- 1
    eigenvalue_dist <- map_dfr(seq(eigenvalues), function(i) {
        lambda <- eigenvalues[[i]]
        if (any(is.na(lambda)))
            return(data.frame())
        lambda <- sort(abs(lambda), decreasing = TRUE)
        data.frame(lambda = lambda, censusdate = names(eigenvalues)[i], rank = seq(lambda))
    })
    eigenvalue_dist$censusdate <- as.Date(eigenvalue_dist$censusdate)
    
    # construct plot
    ds_plot <- eigenvalue_dist %>%
        filter(rank <= num_values) %>% 
        ggplot(aes(x = censusdate, y = lambda, color = as.factor(rank))) + 
        scale_color_manual(values = viridis(7, option = "A")[c(1, 3, 4, 5, 6)]) +
        scale_x_date(limits = as.Date(c("2005-03-01", "2011-01-01")), 
                     date_breaks = "1 year", date_labels = "%Y", expand = c(0.01, 0)) + 
        labs(x = NULL, y = "dynamic stability \n(higher is more unstable)") +
        theme_bw(base_size = 20, base_family = "Helvetica", 
                 base_line_size = 1) + 
        theme(panel.grid.minor = element_line(size = 1)) + 
        guides(color = FALSE)
    ggsave(here::here("figures/maizuru_ds_blank.pdf"),
           ds_plot, width = 10, height = 5)
    ggsave(here::here("figures/maizuru_ds_threshold.pdf"),
           ds_plot + geom_hline(yintercept = 1.0, size = 1, linetype = 2) + 
               geom_blank(), 
           width = 10, height = 5)
    ggsave(here::here("figures/maizuru_ds.pdf"),
           ds_plot + geom_line(size = 1) + 
               geom_hline(yintercept = 1.0, size = 1, linetype = 2), 
           width = 10, height = 5)
}