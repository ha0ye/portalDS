library(viridis)
library(cowplot)

source(here::here("R", "dynamic_stability_functions.R"))
source(here::here("R", "plotting_functions.R"))

make_network_from_ccm_results <- function(ccm_results, 
                                          palette = NULL, 
                                          palette_option = "plasma", 
                                          existing_graph = NULL)
{
    my_graph <- ccm_results %>% 
        identify_ccm_links() %>%
        arrange(target_column) %>%
        select(target_column, lib_column) %>%
        graph_from_data_frame()
    
    if (is.null(palette))
    {
        palette <- viridis(length(V(my_graph)), option = palette_option)
        names(palette) <- as_ids(V(my_graph))
    }
    
    my_graph <- create_layout(my_graph, "circle")
    
    if (!is.null(existing_graph))
    {
        idx <- match(my_graph$name, existing_graph$name)
        my_graph$x <- existing_graph$x[idx]
        my_graph$y <- existing_graph$y[idx]
    }
    
    my_plot <- ggraph(my_graph) + 
        geom_edge_fan(edge_width = 0.5) + 
        geom_node_circle(aes(r = 0.05, fill = name)) + 
        theme_graph(foreground = "black", fg_text_colour = "white", 
                    background = "transparent") + 
        coord_fixed() + 
        theme(text = element_text(family = "Helvetica"),
              panel.border = element_rect(color = NA, fill = NA)) + 
        scale_fill_manual(values = palette)
    
    return(list(plot = my_plot, 
                palette = palette, 
                graph = my_graph))
}

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
